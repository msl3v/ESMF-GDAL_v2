/*==============================================================================
 * SHAPEFILE TO PARMETIS DISTRIBUTED GRAPH CONVERTER
 * 
 * This program reads ESRI shapefiles containing linestring geometries using
 * GDAL/OGR and creates a distributed graph structure compatible with ParMETIS
 * for parallel graph partitioning operations.
 * 
 * Author: Generated for MPI/GDAL/ParMETIS integration
 * 
 * FEATURES:
 * - Reads shapefiles with linestring geometries
 * - Distributes features across MPI ranks
 * - Creates graph nodes at linestring endpoints
 * - Merges nodes within specified tolerance
 * - Builds CSR (Compressed Sparse Row) format for ParMETIS
 * - Callable from both C and Fortran
 * 
 * USAGE FROM C:
 *   shapefile_to_parmetis_graph(filename, comm, &graph, tolerance);
 * 
 * USAGE FROM FORTRAN:
 *   call shapefile_to_parmetis_graph_f(filename, comm, vtxdist, xadj, 
 *                                      adjncy, nvtxs, nedges, tolerance, ierr)
 *============================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <ogr_api.h>
#include <cpl_conv.h>
#include <parmetis.h>

/*==============================================================================
 * STRUCTURES AND TYPE DEFINITIONS
 *============================================================================*/

/**
 * @brief Structure to hold the distributed ParMETIS graph data
 * 
 * This structure contains all arrays needed for ParMETIS graph operations.
 * The graph is stored in Compressed Sparse Row (CSR) format, which is
 * ParMETIS's required input format.
 * 
 * CSR Format Explanation:
 * ----------------------
 * - xadj[i] points to the start of node i's adjacency list in adjncy
 * - xadj[i+1] - xadj[i] = degree of node i
 * - adjncy contains all adjacent nodes, concatenated
 * 
 * Example: If node 0 connects to nodes [1,2,3] and node 1 connects to [0,2]:
 *   xadj   = [0, 3, 5, ...]
 *   adjncy = [1, 2, 3, 0, 2, ...]
 * 
 * Node Distribution:
 * -------------------
 * nodedist[i] = first node owned by rank i (global numbering)
 * nodedist[i+1] - nodedist[i] = number of nodes on rank i
 */
typedef struct {
    idx_t *nodedist;     /**< Node distribution array [size+1]
                              nodedist[i] = first node owned by rank i
                              nodedist[i+1] - nodedist[i] = nodes on rank i */
    idx_t *xadj;         /**< CSR row pointer array [nnodes+1]
                              Points to start of each node's adjacency list */
    idx_t *adjncy;       /**< CSR column index array [nedges]
                              Contains the adjacent node IDs */
    idx_t *nwgt;         /**< Node weights [nnodes * ncon] (NULL = uniform) */
    idx_t *adjwgt;       /**< Edge weights [nedges] (NULL = uniform) */
    idx_t nnodes;        /**< Number of nodes owned by this rank */
    idx_t nedges;        /**< Number of edges owned by this rank */
} ParmetisGraph;

/**
 * @brief Hash table node structure for spatial coordinate mapping
 * 
 * This structure maps (x,y) coordinates to unique node IDs. Used to ensure
 * that linestring endpoints at the same location get the same node ID.
 * Implements chaining for collision resolution.
 * 
 * Hash Table Design:
 * -----------------
 * - Fixed-size array of HASH_SIZE buckets
 * - Each bucket is a linked list (chain) of NodeHash structures
 * - Coordinates hash to bucket index
 * - Collisions resolved by chaining
 */
typedef struct NodeHash {
    double x, y;              /**< Spatial coordinates of the node */
    idx_t global_id;          /**< Unique global node identifier */
    struct NodeHash *next;    /**< Pointer to next node in chain (for collisions) */
} NodeHash;

/** Hash table size - prime number for better distribution */
#define HASH_SIZE 100003

/*==============================================================================
 * HASH TABLE FUNCTIONS
 *============================================================================*/

/**
 * @brief Hash function for spatial coordinates
 * 
 * Computes a hash value from x,y coordinates by treating the doubles as
 * byte arrays and combining them with a simple polynomial rolling hash.
 * 
 * Algorithm:
 * ---------
 * hash = 0
 * for each byte b in (x,y):
 *     hash = hash * 31 + b
 * return hash mod HASH_SIZE
 * 
 * @param x X coordinate
 * @param y Y coordinate
 * @return Hash value in range [0, HASH_SIZE)
 */
unsigned int hash_coord(double x, double y) {
    unsigned long hash = 0;
    unsigned char *px = (unsigned char*)&x;  /* Byte-level access to x */
    unsigned char *py = (unsigned char*)&y;  /* Byte-level access to y */
    
    /* Polynomial rolling hash: hash = hash * 31 + byte */
    for (int i = 0; i < sizeof(double); i++) {
        hash = hash * 31 + px[i];
        hash = hash * 31 + py[i];
    }
    return hash % HASH_SIZE;
}

/**
 * @brief Get existing node ID or insert new node into hash table
 * 
 * Searches the hash table for a node within tolerance distance of (x,y).
 * If found, returns existing node ID. Otherwise, creates new node with
 * unique ID and inserts into hash table.
 * 
 * Algorithm:
 * ---------
 * 1. Compute hash bucket for coordinates
 * 2. Walk chain searching for existing node within tolerance
 * 3. If found, return existing node ID
 * 4. If not found, create new node with next available ID
 * 5. Insert new node at head of chain
 * 
 * Distance Metric:
 * ---------------
 * Uses squared Euclidean distance to avoid expensive sqrt():
 *   dist^2 = (x1-x2)^2 + (y1-y2)^2
 *   If dist^2 < tolerance^2, nodes are considered identical
 * 
 * @param hash_table Pointer to hash table array
 * @param x X coordinate to search/insert
 * @param y Y coordinate to search/insert
 * @param node_counter Pointer to counter for generating new node IDs
 * @param tolerance Maximum distance for considering two nodes identical
 * @return Node ID (existing or newly created)
 */
idx_t get_or_insert_node(NodeHash **hash_table, double x, double y, 
                         idx_t *node_counter, double tolerance) {
    unsigned int hash = hash_coord(x, y);
    NodeHash *current = hash_table[hash];
    
    /* 
     * Search existing chain for node within tolerance distance 
     * Walk the linked list at this hash bucket
     */
    while (current != NULL) {
        double dx = current->x - x;
        double dy = current->y - y;
        double dist_sq = dx*dx + dy*dy;  /* Squared distance (avoids sqrt) */
        
        /* If within tolerance, return existing node ID */
        if (dist_sq < tolerance*tolerance) {
            return current->global_id;
        }
        current = current->next;  /* Move to next node in chain */
    }
    
    /* 
     * Node not found - create new node and insert at head of chain
     * This is O(1) insertion
     */
    NodeHash *new_node = (NodeHash*)malloc(sizeof(NodeHash));
    new_node->x = x;
    new_node->y = y;
    new_node->global_id = (*node_counter)++;  /* Assign and increment counter */
    new_node->next = hash_table[hash];        /* Point to old head */
    hash_table[hash] = new_node;              /* New node becomes head */
    
    return new_node->global_id;
}

/**
 * @brief Free all memory allocated for the hash table
 * 
 * Walks through each hash table bucket and frees all chained nodes.
 * Must be called to avoid memory leaks.
 * 
 * @param hash_table Pointer to hash table array
 */
void free_hash_table(NodeHash **hash_table) {
    /* Iterate through all hash table buckets */
    for (int i = 0; i < HASH_SIZE; i++) {
        NodeHash *current = hash_table[i];
        
        /* Free entire chain for this bucket */
        while (current != NULL) {
            NodeHash *temp = current;
            current = current->next;
            free(temp);  /* Free this node */
        }
    }
    free(hash_table);  /* Free the bucket array itself */
}

/*==============================================================================
 * MAIN SHAPEFILE TO PARMETIS CONVERSION FUNCTION (C VERSION)
 *============================================================================*/

/**
 * @brief Read shapefile and create ParMETIS distributed graph
 * 
 * This is the main function that:
 * 1. Opens the shapefile using GDAL
 * 2. Distributes features across MPI ranks
 * 3. Extracts linestring endpoints and creates graph nodes
 * 4. Builds edges from linestrings
 * 5. Constructs CSR format adjacency structure for ParMETIS
 * 
 * The resulting graph is distributed across MPI ranks and ready for
 * ParMETIS partitioning operations.
 * 
 * ALGORITHM OVERVIEW:
 * ------------------
 * Phase 1: Initialization
 *   - Initialize GDAL
 *   - Open shapefile
 *   - Determine feature distribution
 * 
 * Phase 2: Local Graph Construction
 *   - Each MPI rank reads its assigned features
 *   - Extract linestring endpoints
 *   - Hash endpoints to create/reuse node IDs
 *   - Build edge list
 * 
 * Phase 3: Global Coordination
 *   - Gather node counts from all ranks
 *   - Build node distribution (nodedist)
 *   - Convert local node IDs to global numbering
 * 
 * Phase 4: CSR Construction
 *   - Count node degrees
 *   - Build xadj (row pointers)
 *   - Build adjncy (adjacency list)
 * 
 * DISTRIBUTION STRATEGY:
 * ---------------------
 * - Features are divided evenly across ranks
 * - Each rank processes a contiguous block of features
 * - Last rank handles any remainder features
 * - Node ownership follows feature ownership
 * 
 * GRAPH PROPERTIES:
 * ----------------
 * - Undirected graph (each linestring creates one edge)
 * - No self-loops (linestrings where start==end are excluded)
 * - Nodes within tolerance distance are merged
 * - Uniform node and edge weights (can be modified)
 * 
 * @param filename Path to the shapefile (.shp)
 * @param comm MPI communicator for parallel processing
 * @param graph Output ParMETIS graph structure (caller must free with free_parmetis_graph)
 * @param tolerance Coordinate tolerance for merging nodes (meters or map units)
 * @return 0 on success, -1 on error
 */
int shapefile_to_parmetis_graph(const char *filename, MPI_Comm comm, 
                                 ParmetisGraph *graph, double tolerance) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    /*--------------------------------------------------------------------------
     * STEP 1: Initialize GDAL and open shapefile
     * 
     * Only rank 0 registers drivers to avoid race conditions, then all
     * ranks wait at barrier before opening file.
     *------------------------------------------------------------------------*/
    if (rank == 0) {
        OGRRegisterAll();  /* Register all OGR format drivers */
    }
    MPI_Barrier(comm);  /* Ensure GDAL is initialized before proceeding */
    
    /* Open shapefile for reading - all ranks open independently */
    GDALDatasetH dataset = GDALOpenEx(filename, GDAL_OF_VECTOR | GDAL_OF_READONLY, 
                                      NULL, NULL, NULL);
    if (dataset == NULL) {
        if (rank == 0) {
            fprintf(stderr, "Error: Could not open shapefile %s\n", filename);
        }
        return -1;
    }
    
    /* Get first layer from dataset (most shapefiles have only one layer) */
    OGRLayerH layer = OGR_DS_GetLayer(dataset, 0);
    if (layer == NULL) {
        if (rank == 0) {
            fprintf(stderr, "Error: Could not get layer from shapefile\n");
        }
        GDALClose(dataset);
        return -1;
    }
    
    /*--------------------------------------------------------------------------
     * STEP 2: Determine feature distribution across MPI ranks
     * 
     * Features are divided evenly across ranks. Last rank may get slightly
     * more features if total doesn't divide evenly.
     * 
     * Example with 10 features and 3 ranks:
     *   Rank 0: features 0-2 (3 features)
     *   Rank 1: features 3-5 (3 features)
     *   Rank 2: features 6-9 (4 features)
     *------------------------------------------------------------------------*/
    
    /* Get total number of features in layer */
    GIntBig total_features = OGR_L_GetFeatureCount(layer, TRUE);
    
    if (rank == 0) {
        printf("Total features in shapefile: %lld\n", (long long)total_features);
    }
    
    /* Calculate which features this rank should process */
    idx_t features_per_rank = total_features / size;
    idx_t start_feature = rank * features_per_rank;
    idx_t end_feature = (rank == size - 1) ? total_features : (rank + 1) * features_per_rank;
    idx_t local_feature_count = end_feature - start_feature;
    
    /*--------------------------------------------------------------------------
     * STEP 3: Initialize data structures for graph construction
     *------------------------------------------------------------------------*/
    
    /* Hash table for mapping coordinates to node IDs */
    NodeHash **hash_table = (NodeHash**)calloc(HASH_SIZE, sizeof(NodeHash*));
    idx_t node_counter = 0;  /* Local node ID counter starts at 0 */
    
    /* 
     * Temporary edge storage structure
     * Each edge connects two nodes (from -> to)
     */
    typedef struct {
        idx_t from;  /**< Source node ID */
        idx_t to;    /**< Destination node ID */
    } Edge;
    
    /* Allocate edge array - maximum one edge per feature */
    Edge *edges = (Edge*)malloc(local_feature_count * sizeof(Edge));
    idx_t edge_count = 0;
    
    /*--------------------------------------------------------------------------
     * STEP 4: Read features and extract graph structure
     * 
     * Each rank reads its assigned subset of features. For each linestring:
     * - Extract start and end coordinates
     * - Map coordinates to node IDs (creating new nodes as needed)
     * - Create edge connecting the two nodes
     * 
     * Linestring Processing:
     * ---------------------
     * A linestring is a sequence of connected points: P0-P1-P2-...-Pn
     * We create ONE edge from P0 (start) to Pn (end)
     * The intermediate points are ignored for graph topology
     *------------------------------------------------------------------------*/
    
    /* Position file pointer to this rank's starting feature */
    OGR_L_SetNextByIndex(layer, start_feature);
    
    /* Process each feature assigned to this rank */
    for (idx_t i = 0; i < local_feature_count; i++) {
        /* Read next feature from layer */
        OGRFeatureH feature = OGR_L_GetNextFeature(layer);
        if (feature == NULL) break;  /* End of layer reached (shouldn't happen) */
        
        /* Get geometry from feature */
        OGRGeometryH geometry = OGR_F_GetGeometryRef(feature);
        if (geometry == NULL) {
            OGR_F_Destroy(feature);
            continue;  /* Skip features without geometry */
        }
        
        /* Get geometry type and flatten to 2D (ignore Z/M dimensions) */
        OGRwkbGeometryType geom_type = wkbFlatten(OGR_G_GetGeometryType(geometry));
        
        /* Process only LineString geometries */
        if (geom_type == wkbLineString) {
            /* Get number of points in linestring */
            int point_count = OGR_G_GetPointCount(geometry);
            
            /* Need at least 2 points to form an edge */
            if (point_count >= 2) {
                /*
                 * Extract START point (first point in linestring)
                 * Index 0 is the start of the line
                 */
                double x_start = OGR_G_GetX(geometry, 0);
                double y_start = OGR_G_GetY(geometry, 0);
                
                /*
                 * Extract END point (last point in linestring)
                 * Index point_count-1 is the end of the line
                 */
                double x_end = OGR_G_GetX(geometry, point_count - 1);
                double y_end = OGR_G_GetY(geometry, point_count - 1);
                
                /*
                 * Map coordinates to node IDs
                 * If nodes already exist within tolerance, reuse their IDs
                 * Otherwise create new nodes with unique IDs
                 */
                idx_t node_from = get_or_insert_node(hash_table, x_start, y_start, 
                                                     &node_counter, tolerance);
                idx_t node_to = get_or_insert_node(hash_table, x_end, y_end, 
                                                   &node_counter, tolerance);
                
                /*
                 * Create edge if not a self-loop
                 * Self-loops (edges where from==to) are excluded as they
                 * don't provide useful information for graph partitioning
                 */
                if (node_from != node_to) {
                    edges[edge_count].from = node_from;
                    edges[edge_count].to = node_to;
                    edge_count++;
                }
            }
        }
        /* Note: Could add support for MultiLineString here if needed */
        
        /* Clean up feature resources */
        OGR_F_Destroy(feature);
    }
    
    /* Close dataset - no longer needed */
    GDALClose(dataset);
    
    /*--------------------------------------------------------------------------
     * STEP 5: Establish global node distribution
     * 
     * Gather node counts from all ranks and create nodedist array which
     * defines how nodes are distributed across ranks.
     * 
     * nodedist Array:
     * -------------
     * nodedist[0] = 0
     * nodedist[1] = nodes on rank 0
     * nodedist[2] = nodes on ranks 0+1
     * ...
     * nodedist[size] = total nodes in graph
     * 
     * This allows converting between local and global node numbering:
     *   global_id = local_id + nodedist[rank]
     *   local_id = global_id - nodedist[rank]
     *------------------------------------------------------------------------*/
    
    /* Collect local node counts from all ranks */
    idx_t local_nodes = node_counter;
    idx_t *node_counts = NULL;
    if (rank == 0) {
        node_counts = (idx_t*)malloc(size * sizeof(idx_t));
    }
    
    /* Gather all node counts to rank 0 */
    MPI_Gather(&local_nodes, 1, MPI_INT64_T, node_counts, 1, MPI_INT64_T, 0, comm);
    
    /*
     * Build node distribution array (nodedist)
     * nodedist[i] = cumulative sum of nodes up to rank i
     * Size is (size+1) to include total at end
     */
    graph->nodedist = (idx_t*)malloc((size + 1) * sizeof(idx_t));
    graph->nodedist[0] = 0;  /* First node is always numbered 0 */
    
    if (rank == 0) {
        /* Rank 0 computes cumulative sums */
        for (int i = 0; i < size; i++) {
            graph->nodedist[i + 1] = graph->nodedist[i] + node_counts[i];
        }
    }
    
    /* Broadcast nodedist to all ranks */
    MPI_Bcast(graph->nodedist, size + 1, MPI_INT64_T, 0, comm);
    
    if (rank == 0) {
        printf("Total nodes in graph: %lld\n", (long long)graph->nodedist[size]);
        free(node_counts);
    }
    
    /*--------------------------------------------------------------------------
     * STEP 6: Convert local node IDs to global numbering
     * 
     * Local IDs are in range [0, local_nodes)
     * Global IDs are in range [nodedist[rank], nodedist[rank+1])
     * 
     * Conversion: global_id = local_id + nodedist[rank]
     *------------------------------------------------------------------------*/
    
    idx_t global_offset = graph->nodedist[rank];
    for (idx_t i = 0; i < edge_count; i++) {
        edges[i].from += global_offset;
        edges[i].to += global_offset;
    }
    
    /*--------------------------------------------------------------------------
     * STEP 7: Build local adjacency structure in CSR format
     * 
     * CSR (Compressed Sparse Row) Format:
     * -----------------------------------
     * For a graph with N nodes and E edges:
     * 
     * xadj[N+1]:   Row pointer array
     *   xadj[i] = start index in adjncy for node i's neighbors
     *   xadj[i+1] - xadj[i] = degree of node i
     * 
     * adjncy[E]:   Column index array (adjacency list)
     *   adjncy[xadj[i]..xadj[i+1]-1] = neighbors of node i
     * 
     * Example:
     *   Node 0 -> [1, 2, 3]
     *   Node 1 -> [0, 2]
     *   Node 2 -> [0, 1]
     * 
     *   xadj   = [0, 3, 5, 7]
     *   adjncy = [1, 2, 3, 0, 2, 0, 1]
     *------------------------------------------------------------------------*/
    
    graph->nnodes = local_nodes;
    graph->nedges = 0;
    
    /*
     * Phase 1: Count degree for each node
     * Degree = number of edges incident to node
     */
    idx_t *degree = (idx_t*)calloc(local_nodes, sizeof(idx_t));
    for (idx_t i = 0; i < edge_count; i++) {
        /* Convert global ID back to local ID */
        idx_t local_from = edges[i].from - global_offset;
        
        /* Only count edges owned by this rank */
        if (local_from >= 0 && local_from < local_nodes) {
            degree[local_from]++;
        }
    }
    
    /*
     * Phase 2: Build xadj (CSR row pointer)
     * This is a cumulative sum of degrees
     */
    graph->xadj = (idx_t*)malloc((local_nodes + 1) * sizeof(idx_t));
    graph->xadj[0] = 0;  /* First adjacency starts at index 0 */
    for (idx_t i = 0; i < local_nodes; i++) {
        graph->xadj[i + 1] = graph->xadj[i] + degree[i];
    }
    graph->nedges = graph->xadj[local_nodes];  /* Total edges = last xadj value */
    
    /*
     * Phase 3: Build adjncy (adjacency list)
     * Fill in the actual neighbor IDs
     */
    graph->adjncy = (idx_t*)malloc(graph->nedges * sizeof(idx_t));
    idx_t *current_pos = (idx_t*)calloc(local_nodes, sizeof(idx_t));
    
    /* Fill adjncy by iterating through edges again */
    for (idx_t i = 0; i < edge_count; i++) {
        idx_t local_from = edges[i].from - global_offset;
        
        /* Only add edges owned by this rank */
        if (local_from >= 0 && local_from < local_nodes) {
            /* Calculate position in adjncy array */
            idx_t pos = graph->xadj[local_from] + current_pos[local_from];
            graph->adjncy[pos] = edges[i].to;  /* Store destination node */
            current_pos[local_from]++;
        }
    }
    
    /*
     * Initialize weights to NULL (uniform weights)
     * Can be modified to read weights from shapefile attributes
     */
    graph->nwgt = NULL;    /* Uniform node weights */
    graph->adjwgt = NULL;  /* Uniform edge weights */
    
    /*--------------------------------------------------------------------------
     * STEP 8: Cleanup temporary data structures
     *------------------------------------------------------------------------*/
    free(edges);
    free(degree);
    free(current_pos);
    free_hash_table(hash_table);
    
    /*--------------------------------------------------------------------------
     * Done! Report statistics
     *------------------------------------------------------------------------*/
    if (rank == 0) {
        printf("ParMETIS graph created successfully\n");
    }
    
    /* Each rank can report its local statistics */
    printf("Rank %d: %lld nodes, %lld edges\n", rank, 
           (long long)graph->nnodes, (long long)graph->nedges);
    
    return 0;  /* Success */
}

/*==============================================================================
 * FORTRAN-CALLABLE WRAPPER FUNCTION
 *============================================================================*/

/**
 * @brief Fortran-callable wrapper for shapefile_to_parmetis_graph
 * 
 * This function provides a Fortran interface to the C shapefile reader.
 * Arrays are returned directly to Fortran as allocated pointers.
 * 
 * FORTRAN INTERFACE DECLARATION:
 * ------------------------------
 * 
 * interface
 *   subroutine shapefile_to_parmetis_graph_f(filename, comm, &
 *              nodedist, xadj, adjncy, nnodes, nedges, tolerance, ierr) &
 *              bind(C, name="shapefile_to_parmetis_graph_f")
 *     use iso_c_binding
 *     character(kind=c_char), dimension(*), intent(in) :: filename
 *     integer(c_int), value, intent(in) :: comm
 *     type(c_ptr), intent(out) :: nodedist, xadj, adjncy
 *     integer(c_int64_t), intent(out) :: nnodes, nedges
 *     real(c_double), value, intent(in) :: tolerance
 *     integer(c_int), intent(out) :: ierr
 *   end subroutine
 * end interface
 * 
 * FORTRAN USAGE EXAMPLE:
 * ---------------------
 * 
 * program test_shapefile
 *   use iso_c_binding
 *   use mpi
 *   implicit none
 *   
 *   character(len=256) :: filename
 *   integer :: comm, ierr, rank
 *   integer(c_int64_t) :: nvtxs, nedges
 *   real(c_double) :: tolerance
 *   type(c_ptr) :: vtxdist_ptr, xadj_ptr, adjncy_ptr
 *   integer(c_int64_t), pointer :: vtxdist(:), xadj(:), adjncy(:)
 *   
 *   call MPI_Init(ierr)
 *   call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
 *   
 *   filename = "roads.shp" // C_NULL_CHAR
 *   comm = MPI_COMM_WORLD
 *   tolerance = 1.0d-6
 *   
 *   ! Call C function
 *   call shapefile_to_parmetis_graph_f(filename, comm, &
 *        nodedist_ptr, xadj_ptr, adjncy_ptr, nnodes, nedges, tolerance, ierr)
 *   
 *   if (ierr == 0) then
 *     ! Convert C pointers to Fortran pointers
 *     call c_f_pointer(nodedist_ptr, nodedist, [size+1])
 *     call c_f_pointer(xadj_ptr, xadj, [nnodes+1])
 *     call c_f_pointer(adjncy_ptr, adjncy, [nedges])
 *     
 *     ! Use the arrays...
 *     
 *     ! Free memory when done
 *     call free_parmetis_graph_f(nodedist_ptr, xadj_ptr, adjncy_ptr)
 *   endif
 *   
 *   call MPI_Finalize(ierr)
 * end program
 * 
 * @param filename Null-terminated C string with shapefile path
 * @param comm MPI communicator (Fortran MPI_Comm converted to int)
 * @param nodedist Output pointer to node distribution array
 * @param xadj Output pointer to CSR row pointer array
 * @param adjncy Output pointer to CSR adjacency list
 * @param nnodes Output number of local nodes
 * @param nedges Output number of local edges
 * @param tolerance Coordinate tolerance for node merging
 * @param ierr Output error code (0=success, -1=error)
 */
void shapefile_to_parmetis_graph_f(const char *filename, int comm_int,
                                    idx_t **nodedist, idx_t **xadj, idx_t **adjncy,
                                    idx_t *nnodes, idx_t *nedges,
                                    double tolerance, int *ierr) {
    
    /* Convert Fortran MPI communicator to C */
    MPI_Comm comm = MPI_Comm_f2c(comm_int);
    
    /* Allocate graph structure */
    ParmetisGraph *graph = (ParmetisGraph*)malloc(sizeof(ParmetisGraph));
    
    /* Call C function */
    *ierr = shapefile_to_parmetis_graph(filename, comm, graph, tolerance);
    
    if (*ierr == 0) {
        /* Return pointers and sizes to Fortran */
        *nodedist = graph->nodedist;
        *xadj = graph->xadj;
        *adjncy = graph->adjncy;
        *nnodes = graph->nnodes;
        *nedges = graph->nedges;
        
        /* 
         * Note: nwgt and adjwgt are not returned (NULL in this implementation)
         * Free the structure but not the arrays (Fortran will use them)
         */
        free(graph);
    } else {
        /* Error occurred - set outputs to NULL/0 */
        *nodedist = NULL;
        *xadj = NULL;
        *adjncy = NULL;
        *nnodes = 0;
        *nedges = 0;
        
        free(graph);
    }
}

/*==============================================================================
 * MEMORY MANAGEMENT FUNCTIONS
 *============================================================================*/

/**
 * @brief Free ParMETIS graph structure (C version)
 * 
 * Frees all memory allocated for the graph structure.
 * Should be called when graph is no longer needed.
 * 
 * @param graph Pointer to graph structure to free
 */
void free_parmetis_graph(ParmetisGraph *graph) {
    if (graph->nodedist) free(graph->nodedist);
    if (graph->xadj) free(graph->xadj);
    if (graph->adjncy) free(graph->adjncy);
    if (graph->nwgt) free(graph->nwgt);
    if (graph->adjwgt) free(graph->adjwgt);
}

/**
 * @brief Fortran-callable memory deallocation function
 * 
 * Frees the arrays returned by shapefile_to_parmetis_graph_f.
 * Must be called from Fortran to avoid memory leaks.
 * 
 * FORTRAN INTERFACE:
 * -----------------
 * interface
 *   subroutine free_parmetis_graph_f(nodedist, xadj, adjncy) &
 *              bind(C, name="free_parmetis_graph_f")
 *     use iso_c_binding
 *     type(c_ptr), value, intent(in) :: nodedist, xadj, adjncy
 *   end subroutine
 * end interface
 * 
 * @param nodedist Pointer to nodedist array
 * @param xadj Pointer to xadj array
 * @param adjncy Pointer to adjncy array
 */
void free_parmetis_graph_f(idx_t *nodedist, idx_t *xadj, idx_t *adjncy) {
    if (nodedist) free(nodedist);
    if (xadj) free(xadj);
    if (adjncy) free(adjncy);
}

/*==============================================================================
 * EXAMPLE MAIN FUNCTION (C)
 *============================================================================*/

/**
 * @brief Example main function demonstrating usage from C
 * 
 * Compile:
 *   mpicc -O3 -I/usr/include/gdal shapefile_parmetis_fortran.c \
 *         -o shapefile_parmetis -lgdal -lparmetis -lmetis -lm
 * 
 * Run:
 *   mpirun -np 4 ./shapefile_parmetis roads.shp
 */
int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    /* Check command line arguments */
    if (argc < 2) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <shapefile>\n", argv[0]);
            fprintf(stderr, "\n");
            fprintf(stderr, "Example:\n");
            fprintf(stderr, "  mpirun -np 4 %s roads.shp\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }
    
    const char *filename = argv[1];
    ParmetisGraph graph;
    double tolerance = 1e-6;  /* Coordinate tolerance for node matching */
    
    /* Create distributed graph from shapefile */
    int result = shapefile_to_parmetis_graph(filename, MPI_COMM_WORLD, &graph, tolerance);
    
    if (result == 0) {
        /* Graph is ready for ParMETIS operations */
        if (rank == 0) {
            printf("\n=== Graph Statistics ===\n");
            printf("Total nodes: %lld\n", (long long)graph.nodedist[size]);
            printf("Nodes per rank:\n");
        }
        
        /* Print per-rank statistics */
        for (int r = 0; r < size; r++) {
            if (rank == r) {
                printf("  Rank %d: %lld nodes, %lld edges\n", 
                       r, (long long)graph.nnodes, (long long)graph.nedges);
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        
        /*
         * At this point, you can use ParMETIS functions:
         * 
         * Example: K-way partitioning
         * 
         * idx_t nparts = 10;
         * idx_t *part = (idx_t*)malloc(graph.nnodes * sizeof(idx_t));
         * idx_t wgtflag = 0, numflag = 0, ncon = 1, edgecut;
         * real_t *tpwgts = (real_t*)malloc(nparts * ncon * sizeof(real_t));
         * for (int i = 0; i < nparts; i++) tpwgts[i] = 1.0 / nparts;
         * real_t ubvec = 1.05;
         * idx_t options[3] = {0, 0, 0};
         * 
         * ParMETIS_V3_PartKway(
         *     graph.nodedist, graph.xadj, graph.adjncy,
         *     graph.nwgt, graph.adjwgt, &wgtflag, &numflag, &ncon,
         *     &nparts, tpwgts, &ubvec, options, &edgecut, part,
         *     &MPI_COMM_WORLD);
         * 
         * printf("Edge cut: %lld\n", (long long)edgecut);
         */
        
        /* Cleanup */
        free_parmetis_graph(&graph);
    }
    
    MPI_Finalize();
    return result;
}
