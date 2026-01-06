/*==============================================================================
 * PARMETIS GRAPH BUILDER FOR ESMF LOCSTREAM
 * 
 * This replaces the GDAL-based shapefile reader with the ParMETIS graph
 * builder. Add this to ESMCI_LocStream_F.C
 *============================================================================*/

// Add to includes at top of ESMCI_LocStream_F.C:
// #include "shapefile_parmetis_fortran.c"  // or link library

/*------------------------------------------------------------------------------
 * FTN_X(c_esmc_shapefile_to_graph) - Read shapefile and create graph
 *
 * This function replaces multiple GDAL calls:
 *   - c_esmc_gdal_getnfeatures
 *   - c_esmc_gdal_shpinquire  
 *   - c_esmc_gdal_shpgetcoords
 *   - ExtractPolylineConnectivity
 *
 * All in one MPI-distributed call.
 *----------------------------------------------------------------------------*/
void FTN_X(c_esmc_shapefile_to_graph)(
    char *filename,
    int *local_pet,
    int *pet_count,
    double *tolerance,
    int *totalpoints,      // OUT: Total nodes across all ranks
    int *localcount,       // OUT: Local node count on this rank
    int *totaldims,        // OUT: Number of dimensions (always 2)
    double *coordX,        // OUT: Pre-allocated array for X coords [localcount]
    double *coordY,        // OUT: Pre-allocated array for Y coords [localcount]
    int *nodedist,         // OUT: Node distribution [pet_count+1]
    int *nedges,           // OUT: Number of local edges
    int *xadj,             // OUT: CSR row pointers [localcount+1]
    int *adjncy,           // OUT: CSR adjacency list [nedges]
    int *rc,
    ESMCI_FortranStrLenArg filename_l) {

#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_shapefile_to_graph()"

  // Initialize return code
  if (rc) *rc = ESMC_RC_NOT_IMPL;
  
#ifdef ESMF_GDAL
  
  /*--------------------------------------------------------------------------
   * Call ParMETIS graph builder
   *------------------------------------------------------------------------*/
  
  // Convert MPI communicator
  MPI_Comm comm = MPI_COMM_WORLD;  // ESMF uses MPI_COMM_WORLD
  
  // Create graph structure
  ParmetisGraph graph;
  
  // Call the graph builder
  int result = shapefile_to_parmetis_graph(filename, comm, &graph, *tolerance);
  
  if (result != 0) {
    if (*local_pet == 0) {
      printf("ERROR: shapefile_to_parmetis_graph failed for %s\n", filename);
    }
    if (rc) *rc = ESMC_RC_FILE_READ;
    return;
  }
  
  /*--------------------------------------------------------------------------
   * Copy results to Fortran arrays
   *------------------------------------------------------------------------*/
  
  // Set output scalars
  *totalpoints = (int)graph.nodedist[*pet_count];  // Total nodes
  *localcount = (int)graph.nnodes;                 // Local nodes
  *totaldims = 2;                                   // Always 2D
  *nedges = (int)graph.nedges;                      // Local edges
  
  // Copy node distribution
  for (int i = 0; i <= *pet_count; i++) {
    nodedist[i] = (int)graph.nodedist[i];
  }
  
  // Copy coordinates (if arrays are pre-allocated)
  // Note: Fortran will allocate based on localcount from first call
  if (coordX != NULL && coordY != NULL && *localcount > 0) {
    for (int i = 0; i < *localcount; i++) {
      coordX[i] = graph.node_x[i];
      coordY[i] = graph.node_y[i];
    }
  }
  
  // Copy CSR structure (if arrays are pre-allocated)
  if (xadj != NULL && *localcount > 0) {
    for (int i = 0; i <= *localcount; i++) {
      xadj[i] = (int)graph.xadj[i];
    }
  }
  
  if (adjncy != NULL && *nedges > 0) {
    for (int i = 0; i < *nedges; i++) {
      adjncy[i] = (int)graph.adjncy[i];
    }
  }
  
  /*--------------------------------------------------------------------------
   * Cleanup
   *------------------------------------------------------------------------*/
  free_parmetis_graph(&graph);
  
  if (*local_pet == 0) {
    printf("Successfully created graph from %s\n", filename);
    printf("  Total nodes: %d\n", *totalpoints);
    printf("  Total edges: %d (local on rank 0)\n", *nedges);
  }
  
#endif  // ESMF_GDAL
  
  // Return success
  if (rc) *rc = ESMF_SUCCESS;
  return;
}


/*------------------------------------------------------------------------------
 * FTN_X(c_esmc_shapefile_inquire) - Two-stage call: first get sizes
 *
 * This follows ESMF pattern of calling once to get sizes, then allocating,
 * then calling again with allocated arrays.
 *----------------------------------------------------------------------------*/
void FTN_X(c_esmc_shapefile_inquire)(
    char *filename,
    int *local_pet,
    int *pet_count,
    double *tolerance,
    int *totalpoints,      // OUT: Total nodes across all ranks
    int *localcount,       // OUT: Local node count on this rank
    int *nedges,           // OUT: Number of local edges
    int *totaldims,        // OUT: Number of dimensions (always 2)
    int *nodedist,         // OUT: Node distribution [pet_count+1]
    int *rc,
    ESMCI_FortranStrLenArg filename_l) {

#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_shapefile_inquire()"

  // Initialize return code
  if (rc) *rc = ESMC_RC_NOT_IMPL;
  
#ifdef ESMF_GDAL
  
  // Convert MPI communicator
  MPI_Comm comm = MPI_COMM_WORLD;
  
  // Create graph structure
  ParmetisGraph graph;
  
  // Call the graph builder
  int result = shapefile_to_parmetis_graph(filename, comm, &graph, *tolerance);
  
  if (result != 0) {
    if (*local_pet == 0) {
      printf("ERROR: shapefile_to_parmetis_graph failed for %s\n", filename);
    }
    if (rc) *rc = ESMC_RC_FILE_READ;
    return;
  }
  
  // Return sizes only
  *totalpoints = (int)graph.nodedist[*pet_count];
  *localcount = (int)graph.nnodes;
  *nedges = (int)graph.nedges;
  *totaldims = 2;
  
  // Copy node distribution
  for (int i = 0; i <= *pet_count; i++) {
    nodedist[i] = (int)graph.nodedist[i];
  }
  
  // Store graph data in static storage for subsequent retrieval
  // (Alternative: could store in a map keyed by filename)
  // For now, we'll just free and let the next call rebuild
  // This is inefficient but matches ESMF's pattern
  
  free_parmetis_graph(&graph);
  
#endif  // ESMF_GDAL
  
  if (rc) *rc = ESMF_SUCCESS;
  return;
}


/*------------------------------------------------------------------------------
 * FTN_X(c_esmc_shapefile_getdata) - Second stage: get actual data
 *
 * Called after arrays have been allocated based on sizes from inquire.
 *----------------------------------------------------------------------------*/
void FTN_X(c_esmc_shapefile_getdata)(
    char *filename,
    int *local_pet,
    int *pet_count,
    double *tolerance,
    int *localcount,       // IN: Size of allocated arrays
    int *nedges,           // IN: Number of edges (for validation)
    double *coordX,        // OUT: X coordinates [localcount]
    double *coordY,        // OUT: Y coordinates [localcount]
    int *xadj,             // OUT: CSR row pointers [localcount+1]
    int *adjncy,           // OUT: CSR adjacency list [nedges]
    int *rc,
    ESMCI_FortranStrLenArg filename_l) {

#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_shapefile_getdata()"

  // Initialize return code
  if (rc) *rc = ESMC_RC_NOT_IMPL;
  
#ifdef ESMF_GDAL
  
  // Convert MPI communicator
  MPI_Comm comm = MPI_COMM_WORLD;
  
  // Create graph structure
  ParmetisGraph graph;
  
  // Call the graph builder (yes, we have to re-read - inefficient but safe)
  int result = shapefile_to_parmetis_graph(filename, comm, &graph, *tolerance);
  
  if (result != 0) {
    if (*local_pet == 0) {
      printf("ERROR: shapefile_to_parmetis_graph failed for %s\n", filename);
    }
    if (rc) *rc = ESMC_RC_FILE_READ;
    return;
  }
  
  // Validate sizes match
  if ((int)graph.nnodes != *localcount) {
    if (*local_pet == 0) {
      printf("ERROR: Node count mismatch! Expected %d, got %lld\n", 
             *localcount, (long long)graph.nnodes);
    }
    free_parmetis_graph(&graph);
    if (rc) *rc = ESMC_RC_ARG_SIZE;
    return;
  }
  
  if ((int)graph.nedges != *nedges) {
    if (*local_pet == 0) {
      printf("ERROR: Edge count mismatch! Expected %d, got %lld\n",
             *nedges, (long long)graph.nedges);
    }
    free_parmetis_graph(&graph);
    if (rc) *rc = ESMC_RC_ARG_SIZE;
    return;
  }
  
  // Copy data to Fortran arrays
  for (int i = 0; i < *localcount; i++) {
    coordX[i] = graph.node_x[i];
    coordY[i] = graph.node_y[i];
  }
  
  for (int i = 0; i <= *localcount; i++) {
    xadj[i] = (int)graph.xadj[i];
  }
  
  for (int i = 0; i < *nedges; i++) {
    adjncy[i] = (int)graph.adjncy[i];
  }
  
  // Cleanup
  free_parmetis_graph(&graph);
  
#endif  // ESMF_GDAL
  
  if (rc) *rc = ESMF_SUCCESS;
  return;
}
