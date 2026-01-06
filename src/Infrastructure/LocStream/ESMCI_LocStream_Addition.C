/*==============================================================================
 * ADD THIS TO ESMCI_LocStream_F.C
 * 
 * Insert before the final closing brace (before line 1151)
 * This provides the Fortran interface to the ParMETIS graph builder
 *============================================================================*/

#ifdef ESMF_GDAL

/*------------------------------------------------------------------------------
 * ParMETIS graph builder integration
 * 
 * This replaces the GDAL-based shapefile reader with the ParMETIS graph
 * builder which provides:
 *   - Network topology (CSR format)
 *   - Node coordinates
 *   - Distributed across MPI ranks
 *   - Single function call
 *----------------------------------------------------------------------------*/

// Include the ParMETIS graph builder
// Either include the source directly or link the library
#include "shapefile_parmetis_fortran.c"
// OR if compiled separately:
// extern int shapefile_to_parmetis_graph(const char*, MPI_Comm, ParmetisGraph*, double);
// extern void free_parmetis_graph(ParmetisGraph*);

void FTN_X(c_esmc_shapefile_to_graph_f)(
    char *filename,
    int *comm_int,
    void **nodedist_ptr,     // OUT: Node distribution array pointer
    void **xadj_ptr,          // OUT: CSR row pointers pointer
    void **adjncy_ptr,        // OUT: CSR adjacency list pointer  
    void **node_x_ptr,        // OUT: Node X coordinates pointer
    void **node_y_ptr,        // OUT: Node Y coordinates pointer
    idx_t *nnodes,            // OUT: Number of local nodes
    idx_t *nedges,            // OUT: Number of local edges
    double *tolerance,        // IN: Coordinate tolerance
    int *ierr,                // OUT: Error code
    ESMCI_FortranStrLenArg filename_l) {

#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_shapefile_to_graph_f()"

  // Initialize return code
  *ierr = ESMC_RC_NOT_IMPL;
  
  // Convert Fortran MPI communicator to C
  MPI_Comm comm = MPI_Comm_f2c(*comm_int);
  
  // Get rank for debug output
  int rank;
  MPI_Comm_rank(comm, &rank);
  
  // Allocate graph structure
  ParmetisGraph *graph = (ParmetisGraph*)malloc(sizeof(ParmetisGraph));
  if (graph == NULL) {
    if (rank == 0) {
      printf("ERROR: Failed to allocate ParmetisGraph structure\n");
    }
    *ierr = ESMC_RC_MEM;
    return;
  }
  
  // Call the ParMETIS graph builder
  int result = shapefile_to_parmetis_graph(filename, comm, graph, *tolerance);
  
  if (result != 0) {
    if (rank == 0) {
      printf("ERROR: shapefile_to_parmetis_graph failed for %s\n", filename);
    }
    free(graph);
    *ierr = ESMC_RC_FILE_READ;
    return;
  }
  
  // Return pointers to Fortran
  *nodedist_ptr = (void*)graph->nodedist;
  *xadj_ptr = (void*)graph->xadj;
  *adjncy_ptr = (void*)graph->adjncy;
  *node_x_ptr = (void*)graph->node_x;
  *node_y_ptr = (void*)graph->node_y;
  
  // Return sizes
  *nnodes = graph->nnodes;
  *nedges = graph->nedges;
  
  // Free the structure but NOT the arrays (Fortran will use them)
  free(graph);
  
  if (rank == 0) {
    printf("Successfully created graph from %s\n", filename);
    printf("  Local nodes (rank 0): %lld\n", (long long)*nnodes);
    printf("  Local edges (rank 0): %lld\n", (long long)*nedges);
  }
  
  *ierr = ESMF_SUCCESS;
  return;
}


void FTN_X(c_esmc_free_graph_arrays_f)(
    void **nodedist_ptr,
    void **xadj_ptr,
    void **adjncy_ptr,
    void **node_x_ptr,
    void **node_y_ptr,
    int *rc) {

#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_free_graph_arrays_f()"

  // Initialize return code
  if (rc) *rc = ESMC_RC_NOT_IMPL;
  
  // Free all arrays
  if (*nodedist_ptr) free(*nodedist_ptr);
  if (*xadj_ptr) free(*xadj_ptr);
  if (*adjncy_ptr) free(*adjncy_ptr);
  if (*node_x_ptr) free(*node_x_ptr);
  if (*node_y_ptr) free(*node_y_ptr);
  
  // Null out pointers
  *nodedist_ptr = NULL;
  *xadj_ptr = NULL;
  *adjncy_ptr = NULL;
  *node_x_ptr = NULL;
  *node_y_ptr = NULL;
  
  if (rc) *rc = ESMF_SUCCESS;
  return;
}

#endif  // ESMF_GDAL
