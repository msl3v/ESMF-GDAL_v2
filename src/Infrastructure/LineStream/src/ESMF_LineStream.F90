! $Id$
!
! Earth System Modeling Framework
! Copyright (c) 2002-2025, University Corporation for Atmospheric Research,
! Massachusetts Institute of Technology, Geophysical Fluid Dynamics
! Laboratory, University of Michigan, National Centers for Environmental
! Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
! NASA Goddard Space Flight Center.
! Licensed under the University of Illinois-NCSA License.
!
!==============================================================================
#define ESMF_FILENAME "ESMF_LineStream.F90"
!==============================================================================
!
!     ESMF LineStream module
module ESMF_LineStreamMod
  !
  !==============================================================================
  !
  ! This file contains the LineStream class definition and all LineStream
  ! class methods.
  !
  ! A LineStream represents a collection of line segments distributed across
  ! PETs. Unlike a LocStream (which holds zero-dimensional points), a
  ! LineStream holds one-dimensional segments defined by pairs of node
  ! indices into a shared node coordinate table. Nodes that are shared
  ! between segments (e.g. at river junctions) are stored exactly once.
  !
  ! The primary DistGrid (segDistgrid) distributes *segments* across PETs
  ! so that Fields built on a LineStream are indexed per-segment. A
  ! secondary DistGrid (nodeDistgrid) distributes nodes. Each PET owns
  ! the segments assigned to it and holds coordinate data for all nodes
  ! referenced by those segments (including ghost copies of remote nodes).
  !
  !------------------------------------------------------------------------------
  ! INCLUDES
#include "ESMF.h"
  !------------------------------------------------------------------------------
  !
  ! !USES:
  use ESMF_UtilTypesMod
  use ESMF_UtilMod
  use ESMF_BaseMod
  use ESMF_VMMod
  use ESMF_LogErrMod
  use ESMF_IOUtilMod

  use ESMF_ArraySpecMod
  use ESMF_LocalArrayMod
  use ESMF_DELayoutMod
  use ESMF_DistGridMod
  use ESMF_RHandleMod
  use ESMF_ArrayMod
  use ESMF_ArrayCreateMod
  use ESMF_ArrayGetMod
  use ESMF_InitMacrosMod

#ifdef ESMF_GDAL
  use iso_c_binding
#endif
  implicit none

  !------------------------------------------------------------------------------
  ! !PRIVATE TYPES:
  private

  !------------------------------------------------------------------------------
  ! ! ESMF_LineStreamType
  ! ! Internal data structure for the LineStream class.
  ! !
  ! ! Design rationale:
  ! !   - segDistgrid distributes segments (the primary entity for Fields)
  ! !   - nodeDistgrid distributes the unique nodes
  ! !   - Node coordinates are stored as ESMF_Arrays (one per dimension)
  ! !     on the nodeDistgrid, following the same pattern as LocStream keys
  ! !   - Segment connectivity is stored as an ESMF_Array of shape
  ! !     (2, localNumSegments) on the segDistgrid, holding global node
  ! !     indices for (startNode, endNode) of each segment
  ! !   - Ghost node coordinates are stored as plain Fortran arrays for
  ! !     remote nodes referenced by local segments but owned by other PETs

  type ESMF_LineStreamType
#ifndef ESMF_NO_SEQUENCE
     sequence
#endif

     type(ESMF_Base)             :: base

     ! --- Distribution ---
     type(ESMF_DistGrid)         :: segDistgrid      ! distributes segments across PETs
     type(ESMF_DistGrid)         :: nodeDistgrid     ! distributes owned nodes across PETs
     logical                     :: destroySegDistgrid
     logical                     :: destroyNodeDistgrid
     type(ESMF_Index_Flag)       :: indexflag
     type(ESMF_CoordSys_Flag)    :: coordSys

     ! --- Counts ---
     integer                     :: coordDim          ! number of coordinate dimensions (2 or 3)
     integer                     :: localNodeCount    ! owned nodes on this PET
     integer                     :: localSegCount     ! segments on this PET
     integer                     :: ghostNodeCount    ! ghost nodes on this PET (remote endpoints)
     integer                     :: totalNodeCount    ! global total unique nodes
     integer                     :: totalSegCount     ! global total segments

     ! --- Node coordinates (ESMF_Arrays on nodeDistgrid) ---
     ! One Array per coordinate dimension, each of shape (localNodeCount).
     ! For spherical: nodeCoordsArrays(1) = Lon, nodeCoordsArrays(2) = Lat
     ! For cartesian: nodeCoordsArrays(1) = X,   nodeCoordsArrays(2) = Y
     ! Optional 3rd dimension for 3D coordinates.
     type(ESMF_Array), pointer   :: nodeCoordsArrays(:)  ! size = coordDim
     logical, pointer            :: destroyNodeCoords(:)  ! size = coordDim

     ! --- Ghost node coordinates (plain Fortran arrays) ---
     ! Coordinates of remote nodes needed by local segments but owned
     ! by other PETs. Indexed 1:ghostNodeCount.
     real(ESMF_KIND_R8), pointer :: ghostNodeCoordsX(:)
     real(ESMF_KIND_R8), pointer :: ghostNodeCoordsY(:)
     real(ESMF_KIND_R8), pointer :: ghostNodeCoordsZ(:)   ! only if coordDim==3
     integer, pointer            :: ghostNodeGlobalIds(:)  ! global IDs of ghost nodes

     ! --- Segment connectivity (ESMF_Arrays on segDistgrid) ---
     ! Two 1D Arrays, each of shape (localSegCount), holding 1-based
     ! global node IDs for start and end nodes of each segment.
     ! This follows the same one-Array-per-attribute pattern as LocStream keys.
     type(ESMF_Array)            :: segStartNodeArray   ! start node global IDs (I4)
     type(ESMF_Array)            :: segEndNodeArray     ! end node global IDs (I4)
     logical                     :: destroySegStart
     logical                     :: destroySegEnd

     ! --- Per-segment attributes (ESMF_Arrays on segDistgrid) ---
     type(ESMF_Array)            :: segLengthArray    ! segment lengths (R8)
     logical                     :: destroySegLength
     type(ESMF_Array)            :: segMaskArray       ! integer mask per segment
     logical                     :: destroySegMask
     logical                     :: hasMask

     ! --- Node-to-segment mapping (plain Fortran, for adjacency queries) ---
     ! CSR format: for owned node i, segments referencing it are at
     !   nodeSegIdx( nodeSegPtr(i) : nodeSegPtr(i+1)-1 )
     integer, pointer            :: nodeSegPtr(:)      ! size = localNodeCount+1
     integer, pointer            :: nodeSegIdx(:)      ! size = sum of valences

     ESMF_INIT_DECLARE
  end type ESMF_LineStreamType

  !------------------------------------------------------------------------------
  ! ! ESMF_LineStream
  ! ! The LineStream data structure that is passed between implementation and
  ! ! calling languages.

  type ESMF_LineStream
#ifndef ESMF_NO_SEQUENCE
     sequence
#endif
     type(ESMF_LineStreamType), pointer :: lstypep
     ESMF_INIT_DECLARE
  end type ESMF_LineStream

  !------------------------------------------------------------------------------
  ! !PUBLIC TYPES:
  public ESMF_LineStream
  public ESMF_LineStreamType   ! For internal use only

  public operator(==)
  public operator(/=)

  public ESMF_LineStreamIsCreated
  public ESMF_LineStreamCreate
  public ESMF_LineStreamGet
  public ESMF_LineStreamGetBounds
  public ESMF_LineStreamDestroy
  public ESMF_LineStreamDestruct        ! for ESMF garbage collection
  public ESMF_LineStreamGetNodeCoords
  public ESMF_LineStreamGetSegConn
  public ESMF_LineStreamGetSegLength
  public ESMF_LineStreamPrint

  !EOPI

  !------------------------------------------------------------------------------
  character(*), parameter, private :: version = &
       '$Id$'

  !==============================================================================
  !
  ! INTERFACE BLOCKS
  !
  !==============================================================================

  interface ESMF_LineStreamCreate

     module procedure ESMF_LineStreamCreateFromFile

     !EOPI
  end interface ESMF_LineStreamCreate

  interface operator(==)
     module procedure ESMF_LineStreamEQ
  end interface operator(==)

  interface operator(/=)
     module procedure ESMF_LineStreamNE
  end interface operator(/=)


  ! C interface for shapefile reading (reuse existing GDAL/ParMETIS infrastructure)
#ifdef ESMF_GDAL
  interface
    subroutine shapefile_to_parmetis_graph_f(filename, comm, &
         nodedist, xadj, adjncy, node_x, node_y, nnodes, nedges, tolerance, ierr) &
         bind(C, name="shapefile_to_parmetis_graph_f")
      use iso_c_binding
      character(kind=c_char), dimension(*), intent(in) :: filename
      integer(c_int), value, intent(in) :: comm
      type(c_ptr), intent(out) :: nodedist, xadj, adjncy, node_x, node_y
      integer(c_int), intent(out) :: nnodes, nedges
      real(c_double), value, intent(in) :: tolerance
      integer(c_int), intent(out) :: ierr
    end subroutine

    subroutine free_parmetis_graph_f(nodedist, xadj, adjncy) &
         bind(C, name="free_parmetis_graph_f")
      use iso_c_binding
      type(c_ptr), value, intent(in) :: nodedist, xadj, adjncy
    end subroutine

    subroutine c_free(ptr) bind(C, name="free")
      use iso_c_binding
      type(c_ptr), value, intent(in) :: ptr
    end subroutine
  end interface
#endif

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamEQ"
  function ESMF_LineStreamEQ(linestream1, linestream2)
    logical :: ESMF_LineStreamEQ
    type(ESMF_LineStream), intent(in) :: linestream1
    type(ESMF_LineStream), intent(in) :: linestream2

    ESMF_LineStreamEQ = associated(linestream1%lstypep, linestream2%lstypep)
  end function ESMF_LineStreamEQ


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamNE"
  function ESMF_LineStreamNE(linestream1, linestream2)
    logical :: ESMF_LineStreamNE
    type(ESMF_LineStream), intent(in) :: linestream1
    type(ESMF_LineStream), intent(in) :: linestream2

    ESMF_LineStreamNE = .not. associated(linestream1%lstypep, linestream2%lstypep)
  end function ESMF_LineStreamNE


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamIsCreated"
!BOP
! !IROUTINE: ESMF_LineStreamIsCreated - Check whether a LineStream object has been created
!
! !INTERFACE:
  function ESMF_LineStreamIsCreated(linestream, keywordEnforcer, rc)
!
! !RETURN VALUE:
    logical :: ESMF_LineStreamIsCreated
!
! !ARGUMENTS:
    type(ESMF_LineStream), intent(in)            :: linestream
    type(ESMF_KeywordEnforcer), optional:: keywordEnforcer
    integer,               intent(out), optional :: rc
!
! !DESCRIPTION:
!   Returns .TRUE. if the {\tt linestream} has been created.
!EOP

    if (present(rc)) rc = ESMF_SUCCESS
    ESMF_LineStreamIsCreated = .false.
    if (ESMF_LineStreamGetInit(linestream)==ESMF_INIT_CREATED) &
         ESMF_LineStreamIsCreated = .true.

  end function ESMF_LineStreamIsCreated


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamCreateFromFile"
!BOP
! !IROUTINE: ESMF_LineStreamCreate - Create a LineStream from a shapefile
!
! !INTERFACE:
  function ESMF_LineStreamCreateFromFile(filename, keywordEnforcer, &
       indexflag, coordSys, tolerance, name, rc)
!
! !RETURN VALUE:
    type(ESMF_LineStream) :: ESMF_LineStreamCreateFromFile
!
! !ARGUMENTS:
    character(len=*),           intent(in)            :: filename
    type(ESMF_KeywordEnforcer), optional:: keywordEnforcer ! must use keywords below
    type(ESMF_Index_Flag),      intent(in),  optional :: indexflag
    type(ESMF_CoordSys_Flag),   intent(in),  optional :: coordSys
    real(ESMF_KIND_R8),         intent(in),  optional :: tolerance
    character(len=*),           intent(in),  optional :: name
    integer,                    intent(out), optional :: rc
!
! !DESCRIPTION:
!     Create a new {\tt ESMF\_LineStream} object from a shapefile containing
!     polyline features. The shapefile is read using GDAL, nodes are
!     deduplicated using a coordinate tolerance, and the resulting node/segment
!     graph is distributed across PETs.
!
!     The shapefile must contain line or polyline geometry. Each polyline
!     feature is decomposed into individual two-point segments. Shared
!     endpoints (e.g. at river junctions) are stored once in the node table.
!
!     The arguments are:
!     \begin{description}
!     \item[filename]
!          Name of the shapefile (.shp) to read.
!     \item[{[indexflag]}]
!          Flag that indicates how the DE-local indices are to be defined.
!          Defaults to {\tt ESMF\_INDEX\_DELOCAL}.
!     \item[{[coordSys]}]
!          The coordinate system. Defaults to {\tt ESMF\_COORDSYS\_SPH\_DEG}.
!     \item[{[tolerance]}]
!          Coordinate tolerance for node deduplication. Two points within
!          this distance are considered the same node. Defaults to 1.0e-6.
!     \item[{[name]}]
!          Name of the LineStream object.
!     \item[{[rc]}]
!          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
!     \end{description}
!
!EOP

    ! Local variables
    type(ESMF_LineStream)              :: linestream
    type(ESMF_LineStreamType), pointer :: lstypep
    type(ESMF_VM)                      :: vm
    type(ESMF_Index_Flag)              :: indexflagLocal
    type(ESMF_CoordSys_Flag)           :: coordSysLocal
    real(ESMF_KIND_R8)                 :: tolLocal
    integer                            :: localrc
    integer                            :: PetNo, PetCnt
    integer                            :: mpi_comm

    ! DistGrid construction
    integer, allocatable :: segCountsPerPet(:)
    integer, allocatable :: nodeCountsPerPet(:)
    integer, pointer     :: segDeBlockList(:,:,:)
    integer, pointer     :: nodeDeBlockList(:,:,:)
    integer              :: segMinIndex(1), segMaxIndex(1)
    integer              :: nodeMinIndex(1), nodeMaxIndex(1)
    integer              :: currMin, i, j, k

    ! Node coordinate and segment connectivity working arrays
    real(ESMF_KIND_R8), allocatable :: localNodeX(:), localNodeY(:)
    integer, allocatable           :: localSegSrc(:), localSegDst(:)
    real(ESMF_KIND_R8), allocatable :: localSegLen(:)

    ! Ghost node tracking
    integer, allocatable :: ghostIds(:)
    integer              :: ghostCount, maxGhost
    logical              :: found
    integer              :: gIdx

    ! For segment length calculation
    real(ESMF_KIND_R8) :: x1, y1, x2, y2, dx, dy
    real(ESMF_KIND_R8), parameter :: DEG2RAD = 3.14159265358979323846d0 / 180.0d0
    real(ESMF_KIND_R8), parameter :: EARTH_RADIUS = 6371000.0d0  ! meters

    ! ParMETIS data structures
#ifdef ESMF_GDAL
    type(c_ptr) :: nodedist_ptr, xadj_ptr, adjncy_ptr, node_x_ptr, node_y_ptr
    integer(c_int32_t), pointer :: nodedist(:), xadj(:), adjncy(:)
    real(c_double), pointer :: node_x(:), node_y(:)
    integer(c_int) :: nnodes, nedges, c_ierr
    character(kind=c_char, len=256) :: c_filename
    real(c_double) :: tol_val

    ! Edge construction from CSR
    integer :: edge_local, my_node_start, my_node_end
    integer :: src_global, dst_global
    integer :: localSegs, totalSegs, totalNodes
#endif

    ! -----------------------------------------------------------------
    ! Initialize return code
    ! -----------------------------------------------------------------
    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    nullify(ESMF_LineStreamCreateFromFile%lstypep)

#ifdef ESMF_GDAL

    ! -----------------------------------------------------------------
    ! Set defaults for optional arguments
    ! -----------------------------------------------------------------
    if (present(indexflag)) then
       indexflagLocal = indexflag
    else
       indexflagLocal = ESMF_INDEX_DELOCAL
    endif

    if (present(coordSys)) then
       coordSysLocal = coordSys
    else
       coordSysLocal = ESMF_COORDSYS_SPH_DEG
    endif

    if (present(tolerance)) then
       tolLocal = tolerance
    else
       tolLocal = 1.0d-6
    endif


    ! -----------------------------------------------------------------
    ! Get VM information
    ! -----------------------------------------------------------------
    call ESMF_VMGetCurrent(vm, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    call ESMF_VMGet(vm, localPet=PetNo, petCount=PetCnt, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    call ESMF_VMGet(vm, mpiCommunicator=mpi_comm, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return


    ! -----------------------------------------------------------------
    ! Read shapefile via C/GDAL and build ParMETIS distributed graph
    ! This reuses the same C function as LocStreamCreateFromParmetis.
    ! Returns:
    !   nodedist(PetCnt+1) - node distribution across PETs (0-based CSR)
    !   xadj(nnodes+1)     - CSR row pointers for adjacency (0-based)
    !   adjncy(nedges)     - CSR column indices (0-based global node IDs)
    !   node_x(nnodes)     - longitude of local nodes
    !   node_y(nnodes)     - latitude of local nodes
    !   nnodes             - number of local nodes on this PET
    !   nedges             - number of local adjacency entries (edges from local nodes)
    ! -----------------------------------------------------------------
    tol_val = tolLocal
    c_filename = trim(filename) // C_NULL_CHAR

    call shapefile_to_parmetis_graph_f(c_filename, mpi_comm, &
         nodedist_ptr, xadj_ptr, adjncy_ptr, node_x_ptr, node_y_ptr, &
         nnodes, nedges, tol_val, c_ierr)

    if (c_ierr /= 0) then
       call ESMF_LogSetError(rcToCheck=ESMF_FAILURE, &
            msg="Failed to read shapefile and build distributed graph", &
            ESMF_CONTEXT, rcToReturn=rc)
       return
    endif

    ! Convert C pointers to Fortran pointers
    call c_f_pointer(nodedist_ptr, nodedist, [PetCnt+1])
    call c_f_pointer(node_x_ptr, node_x, [nnodes])
    call c_f_pointer(node_y_ptr, node_y, [nnodes])
    call c_f_pointer(xadj_ptr, xadj, [nnodes+1])
    call c_f_pointer(adjncy_ptr, adjncy, [nedges])

    ! Compute global range of nodes owned by this PET (0-based)
    my_node_start = int(nodedist(PetNo+1))
    my_node_end   = int(nodedist(PetNo+2)) - 1

    totalNodes = int(nodedist(PetCnt+1))

    ! -----------------------------------------------------------------
    ! Convert CSR adjacency to segment list.
    !
    ! The ParMETIS graph stores edges bidirectionally (if node A connects
    ! to node B, both A->B and B->A appear). To avoid duplicate segments,
    ! we only create a segment when src_global < dst_global.
    ! Each segment is owned by the PET that owns the source node
    ! (the node with the smaller global ID).
    ! -----------------------------------------------------------------

    ! First pass: count local segments
    localSegs = 0
    do i = 0, int(nnodes)-1
       src_global = my_node_start + i   ! 0-based global node ID
       do j = int(xadj(i+1)), int(xadj(i+2))-1
          dst_global = int(adjncy(j+1))  ! 0-based global node ID
          if (src_global < dst_global) then
             localSegs = localSegs + 1
          endif
       enddo
    enddo

    ! Allocate working arrays for local segments
    allocate(localSegSrc(localSegs), localSegDst(localSegs))
    allocate(localSegLen(localSegs))

    ! Second pass: fill segment arrays (using 1-based global node IDs)
    edge_local = 0
    do i = 0, int(nnodes)-1
       src_global = my_node_start + i   ! 0-based
       do j = int(xadj(i+1)), int(xadj(i+2))-1
          dst_global = int(adjncy(j+1))  ! 0-based
          if (src_global < dst_global) then
             edge_local = edge_local + 1
             localSegSrc(edge_local) = src_global + 1   ! convert to 1-based
             localSegDst(edge_local) = dst_global + 1   ! convert to 1-based
          endif
       enddo
    enddo

    ! Copy local node coordinates
    allocate(localNodeX(nnodes), localNodeY(nnodes))
    localNodeX(1:nnodes) = node_x(1:nnodes)
    localNodeY(1:nnodes) = node_y(1:nnodes)

    ! Free C-allocated memory (we've copied everything to Fortran arrays)
    ! Note: free_parmetis_graph_f only frees nodedist, xadj, adjncy.
    ! node_x and node_y are separate C allocations that we need to free
    ! via the C library. Since there is no dedicated free function for
    ! these, we use the iso_c_binding interface to call C free() directly.
    call free_parmetis_graph_f(nodedist_ptr, xadj_ptr, adjncy_ptr)
    call c_free(node_x_ptr)
    call c_free(node_y_ptr)


    ! -----------------------------------------------------------------
    ! Identify ghost nodes: destination nodes of local segments that
    ! are owned by other PETs (their global ID is outside my_node_start..my_node_end).
    ! We need their coordinates for segment length calculation and
    ! for PointList creation during regridding.
    ! -----------------------------------------------------------------
    maxGhost = localSegs  ! upper bound on ghost count
    allocate(ghostIds(maxGhost))
    ghostCount = 0

    do i = 1, localSegs
       ! Check if the destination node is remote
       ! (src is always local because we only keep src < dst with src owned locally)
       dst_global = localSegDst(i)   ! 1-based
       if (dst_global < my_node_start+1 .or. dst_global > my_node_end+1) then
          ! Check if we already have this ghost
          found = .false.
          do k = 1, ghostCount
             if (ghostIds(k) == dst_global) then
                found = .true.
                exit
             endif
          enddo
          if (.not. found) then
             ghostCount = ghostCount + 1
             ghostIds(ghostCount) = dst_global
          endif
       endif
    enddo


    ! -----------------------------------------------------------------
    ! Gather segment counts across PETs and compute total
    ! -----------------------------------------------------------------
    allocate(segCountsPerPet(PetCnt))
    call ESMF_VMAllGather(vm, sendData=(/localSegs/), &
         recvData=segCountsPerPet, count=1, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    totalSegs = sum(segCountsPerPet)


    ! -----------------------------------------------------------------
    ! Create segment DistGrid (primary distribution for Fields)
    ! -----------------------------------------------------------------
    segMinIndex(1) = 1
    segMaxIndex(1) = totalSegs

    allocate(segDeBlockList(1,2,PetCnt))
    currMin = 1
    do i = 1, PetCnt
       segDeBlockList(1,1,i) = currMin
       segDeBlockList(1,2,i) = currMin + segCountsPerPet(i) - 1
       currMin = segDeBlockList(1,2,i) + 1
    enddo


    ! -----------------------------------------------------------------
    ! Create node DistGrid (for node coordinate Arrays)
    ! -----------------------------------------------------------------
    allocate(nodeCountsPerPet(PetCnt))
    call ESMF_VMAllGather(vm, sendData=(/int(nnodes)/), &
         recvData=nodeCountsPerPet, count=1, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    nodeMinIndex(1) = 1
    nodeMaxIndex(1) = totalNodes

    allocate(nodeDeBlockList(1,2,PetCnt))
    currMin = 1
    do i = 1, PetCnt
       nodeDeBlockList(1,1,i) = currMin
       nodeDeBlockList(1,2,i) = currMin + nodeCountsPerPet(i) - 1
       currMin = nodeDeBlockList(1,2,i) + 1
    enddo


    ! -----------------------------------------------------------------
    ! Allocate and initialize the LineStreamType
    ! -----------------------------------------------------------------
    nullify(lstypep)
    allocate(lstypep, stat=localrc)
    if (ESMF_LogFoundAllocError(localrc, &
         msg="Allocating LineStreamType object", &
         ESMF_CONTEXT, rcToReturn=rc)) return

    ! Initialize pointer members to null
    nullify(lstypep%nodeCoordsArrays)
    nullify(lstypep%destroyNodeCoords)
    nullify(lstypep%ghostNodeCoordsX)
    nullify(lstypep%ghostNodeCoordsY)
    nullify(lstypep%ghostNodeCoordsZ)
    nullify(lstypep%ghostNodeGlobalIds)
    nullify(lstypep%nodeSegPtr)
    nullify(lstypep%nodeSegIdx)

    lstypep%indexflag = indexflagLocal
    lstypep%coordSys  = coordSysLocal
    lstypep%coordDim  = 2   ! shapefile gives us 2D (lon/lat)

    lstypep%localNodeCount  = int(nnodes)
    lstypep%localSegCount   = localSegs
    lstypep%ghostNodeCount  = ghostCount
    lstypep%totalNodeCount  = totalNodes
    lstypep%totalSegCount   = totalSegs

    lstypep%hasMask         = .false.


    ! -----------------------------------------------------------------
    ! Create the DistGrids
    ! -----------------------------------------------------------------
    lstypep%segDistgrid = ESMF_DistGridCreate( &
         minIndex=segMinIndex, maxIndex=segMaxIndex, &
         deBlockList=segDeBlockList, &
         indexflag=indexflagLocal, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return
    lstypep%destroySegDistgrid = .true.

    lstypep%nodeDistgrid = ESMF_DistGridCreate( &
         minIndex=nodeMinIndex, maxIndex=nodeMaxIndex, &
         deBlockList=nodeDeBlockList, &
         indexflag=indexflagLocal, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return
    lstypep%destroyNodeDistgrid = .true.

    deallocate(segDeBlockList, nodeDeBlockList)
    deallocate(segCountsPerPet, nodeCountsPerPet)


    ! -----------------------------------------------------------------
    ! Create node coordinate Arrays on nodeDistgrid
    ! Following the LocStream pattern: one ESMF_Array per coord dimension.
    ! -----------------------------------------------------------------
    allocate(lstypep%nodeCoordsArrays(lstypep%coordDim))
    allocate(lstypep%destroyNodeCoords(lstypep%coordDim))

    ! Longitude / X
    lstypep%nodeCoordsArrays(1) = ESMF_ArrayCreate( &
         lstypep%nodeDistgrid, localNodeX, &
         datacopyflag=ESMF_DATACOPY_VALUE, &
         indexflag=indexflagLocal, &
         name="ESMF:NodeCoordX", rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return
    lstypep%destroyNodeCoords(1) = .true.

    ! Latitude / Y
    lstypep%nodeCoordsArrays(2) = ESMF_ArrayCreate( &
         lstypep%nodeDistgrid, localNodeY, &
         datacopyflag=ESMF_DATACOPY_VALUE, &
         indexflag=indexflagLocal, &
         name="ESMF:NodeCoordY", rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return
    lstypep%destroyNodeCoords(2) = .true.


    ! -----------------------------------------------------------------
    ! Store ghost node information
    ! Note: Ghost coordinates must be obtained via communication from
    ! the owning PETs. For this initial implementation, we allocate
    ! the arrays and populate them via an AllToAllV exchange.
    ! -----------------------------------------------------------------
    if (ghostCount > 0) then
       allocate(lstypep%ghostNodeGlobalIds(ghostCount))
       allocate(lstypep%ghostNodeCoordsX(ghostCount))
       allocate(lstypep%ghostNodeCoordsY(ghostCount))
       lstypep%ghostNodeGlobalIds(1:ghostCount) = ghostIds(1:ghostCount)

       ! TODO: Implement ghost node coordinate exchange.
       ! This requires an AllToAllV or point-to-point communication
       ! where each PET sends the coordinates of nodes requested by
       ! other PETs. For now, initialize to zero; the exchange should
       ! be implemented before segment length computation is valid for
       ! segments with remote endpoints.
       lstypep%ghostNodeCoordsX(:) = 0.0d0
       lstypep%ghostNodeCoordsY(:) = 0.0d0
    else
       nullify(lstypep%ghostNodeGlobalIds)
       nullify(lstypep%ghostNodeCoordsX)
       nullify(lstypep%ghostNodeCoordsY)
    endif
    nullify(lstypep%ghostNodeCoordsZ)  ! no 3D for shapefiles

    deallocate(ghostIds)


    ! -----------------------------------------------------------------
    ! Create segment connectivity Arrays on segDistgrid
    ! Two 1D integer Arrays, one for start nodes and one for end nodes,
    ! following the same one-Array-per-attribute pattern as LocStream keys.
    ! -----------------------------------------------------------------
    lstypep%segStartNodeArray = ESMF_ArrayCreate( &
         lstypep%segDistgrid, localSegSrc, &
         datacopyflag=ESMF_DATACOPY_VALUE, &
         indexflag=indexflagLocal, &
         name="SegStartNode", rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return
    lstypep%destroySegStart = .true.

    lstypep%segEndNodeArray = ESMF_ArrayCreate( &
         lstypep%segDistgrid, localSegDst, &
         datacopyflag=ESMF_DATACOPY_VALUE, &
         indexflag=indexflagLocal, &
         name="SegEndNode", rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return
    lstypep%destroySegEnd = .true.


    ! -----------------------------------------------------------------
    ! Compute segment lengths
    ! For segments where both endpoints are local, compute directly.
    ! For segments with a ghost endpoint, length requires ghost coords
    ! (marked TODO above). For now, compute for local-only segments.
    ! -----------------------------------------------------------------
    do i = 1, localSegs
       src_global = localSegSrc(i)  ! 1-based
       dst_global = localSegDst(i)  ! 1-based

       ! Get source coords (always local since src < dst and src is owned)
       x1 = localNodeX(src_global - my_node_start)  ! local index (1-based)
       y1 = localNodeY(src_global - my_node_start)

       ! Get destination coords
       if (dst_global >= my_node_start+1 .and. dst_global <= my_node_end+1) then
          ! Local destination
          x2 = localNodeX(dst_global - my_node_start)
          y2 = localNodeY(dst_global - my_node_start)
       else
          ! Ghost destination - use ghost coordinates (TODO: after exchange)
          ! Find ghost index
          x2 = 0.0d0
          y2 = 0.0d0
          if (associated(lstypep%ghostNodeGlobalIds)) then
             do k = 1, lstypep%ghostNodeCount
                if (lstypep%ghostNodeGlobalIds(k) == dst_global) then
                   x2 = lstypep%ghostNodeCoordsX(k)
                   y2 = lstypep%ghostNodeCoordsY(k)
                   exit
                endif
             enddo
          endif
       endif

       ! Compute great-circle distance for spherical coords
       if (coordSysLocal == ESMF_COORDSYS_SPH_DEG) then
          call haversine_distance(y1, x1, y2, x2, localSegLen(i))
       else
          ! Cartesian distance
          dx = x2 - x1
          dy = y2 - y1
          localSegLen(i) = sqrt(dx*dx + dy*dy)
       endif
    enddo

    ! Create segment length Array
    lstypep%segLengthArray = ESMF_ArrayCreate( &
         lstypep%segDistgrid, localSegLen, &
         datacopyflag=ESMF_DATACOPY_VALUE, &
         indexflag=indexflagLocal, &
         name="SegLength", rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return
    lstypep%destroySegLength = .true.

    lstypep%destroySegMask = .false.   ! no mask created yet


    ! -----------------------------------------------------------------
    ! Build node-to-segment adjacency (CSR format)
    ! For each owned node, list which local segments reference it.
    ! -----------------------------------------------------------------
    block
       integer, allocatable :: valence(:)
       integer :: nodeLocalIdx, pos

       allocate(lstypep%nodeSegPtr(int(nnodes)+1))
       allocate(valence(int(nnodes)))
       valence(:) = 0

       ! Count valence: how many local segments reference each local node
       do i = 1, localSegs
          src_global = localSegSrc(i)  ! 1-based
          if (src_global >= my_node_start+1 .and. src_global <= my_node_end+1) then
             nodeLocalIdx = src_global - my_node_start  ! 1-based local
             valence(nodeLocalIdx) = valence(nodeLocalIdx) + 1
          endif
          dst_global = localSegDst(i)  ! 1-based
          if (dst_global >= my_node_start+1 .and. dst_global <= my_node_end+1) then
             nodeLocalIdx = dst_global - my_node_start  ! 1-based local
             valence(nodeLocalIdx) = valence(nodeLocalIdx) + 1
          endif
       enddo

       ! Build pointer array (exclusive prefix sum)
       lstypep%nodeSegPtr(1) = 1
       do i = 1, int(nnodes)
          lstypep%nodeSegPtr(i+1) = lstypep%nodeSegPtr(i) + valence(i)
       enddo

       ! Allocate index array
       allocate(lstypep%nodeSegIdx(lstypep%nodeSegPtr(int(nnodes)+1)-1))

       ! Reset valence as position counter
       valence(:) = 0

       ! Fill index array
       do i = 1, localSegs
          src_global = localSegSrc(i)
          if (src_global >= my_node_start+1 .and. src_global <= my_node_end+1) then
             nodeLocalIdx = src_global - my_node_start
             pos = lstypep%nodeSegPtr(nodeLocalIdx) + valence(nodeLocalIdx)
             lstypep%nodeSegIdx(pos) = i
             valence(nodeLocalIdx) = valence(nodeLocalIdx) + 1
          endif
          dst_global = localSegDst(i)
          if (dst_global >= my_node_start+1 .and. dst_global <= my_node_end+1) then
             nodeLocalIdx = dst_global - my_node_start
             pos = lstypep%nodeSegPtr(nodeLocalIdx) + valence(nodeLocalIdx)
             lstypep%nodeSegIdx(pos) = i
             valence(nodeLocalIdx) = valence(nodeLocalIdx) + 1
          endif
       enddo

       deallocate(valence)
    end block


    ! -----------------------------------------------------------------
    ! Clean up working arrays
    ! -----------------------------------------------------------------
    deallocate(localNodeX, localNodeY)
    deallocate(localSegSrc, localSegDst, localSegLen)


    ! -----------------------------------------------------------------
    ! Initialize ESMF_Base and finalize the LineStream object
    ! -----------------------------------------------------------------
    call ESMF_BaseCreate(lstypep%base, "LineStream", name, 0, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    linestream%lstypep => lstypep

    ! Register with ESMF garbage collection
    ! Note: ESMF_ID_LINESTREAM would need to be registered in ESMF_InitMacrosMod.
    ! For now, use the base registration mechanism.
    ! call c_ESMC_VMAddFObject(linestream, ESMF_ID_LINESTREAM%objectID)

    ESMF_LineStreamCreateFromFile = linestream

    ESMF_INIT_SET_CREATED(ESMF_LineStreamCreateFromFile)

    if (present(rc)) rc = ESMF_SUCCESS

#else
    ! GDAL not available
    call ESMF_LogSetError(rcToCheck=ESMF_RC_LIB_NOT_PRESENT, &
         msg="ESMF not compiled with GDAL support - shapefile reading unavailable", &
         ESMF_CONTEXT, rcToReturn=rc)
    if (present(rc)) rc = ESMF_RC_LIB_NOT_PRESENT
#endif

  end function ESMF_LineStreamCreateFromFile


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamGet"
!BOP
! !IROUTINE: ESMF_LineStreamGet - Get information from a LineStream
!
! !INTERFACE:
  subroutine ESMF_LineStreamGet(linestream, keywordEnforcer, &
       segDistgrid, nodeDistgrid, localSegCount, localNodeCount, &
       ghostNodeCount, totalSegCount, totalNodeCount, &
       coordDim, coordSys, indexflag, localDECount, name, rc)
!
! !ARGUMENTS:
    type(ESMF_LineStream),      intent(in)            :: linestream
    type(ESMF_KeywordEnforcer), optional:: keywordEnforcer ! must use keywords below
    type(ESMF_DistGrid),        intent(out), optional :: segDistgrid
    type(ESMF_DistGrid),        intent(out), optional :: nodeDistgrid
    integer,                    intent(out), optional :: localSegCount
    integer,                    intent(out), optional :: localNodeCount
    integer,                    intent(out), optional :: ghostNodeCount
    integer,                    intent(out), optional :: totalSegCount
    integer,                    intent(out), optional :: totalNodeCount
    integer,                    intent(out), optional :: coordDim
    type(ESMF_CoordSys_Flag),   intent(out), optional :: coordSys
    type(ESMF_Index_Flag),      intent(out), optional :: indexflag
    integer,                    intent(out), optional :: localDECount
    character(len=*),           intent(out), optional :: name
    integer,                    intent(out), optional :: rc
!
! !DESCRIPTION:
!   Query an {\tt ESMF\_LineStream} for various information.
!EOP

    type(ESMF_LineStreamType), pointer :: lstypep
    type(ESMF_DELayout) :: delayout
    integer :: localrc

    localrc = ESMF_RC_NOT_IMPL
    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    ESMF_INIT_CHECK_DEEP(ESMF_LineStreamGetInit, linestream, rc)

    lstypep => linestream%lstypep

    if (present(segDistgrid))    segDistgrid   = lstypep%segDistgrid
    if (present(nodeDistgrid))   nodeDistgrid  = lstypep%nodeDistgrid
    if (present(localSegCount))  localSegCount = lstypep%localSegCount
    if (present(localNodeCount)) localNodeCount= lstypep%localNodeCount
    if (present(ghostNodeCount)) ghostNodeCount= lstypep%ghostNodeCount
    if (present(totalSegCount))  totalSegCount = lstypep%totalSegCount
    if (present(totalNodeCount)) totalNodeCount= lstypep%totalNodeCount
    if (present(coordDim))       coordDim      = lstypep%coordDim
    if (present(coordSys))       coordSys      = lstypep%coordSys
    if (present(indexflag))      indexflag     = lstypep%indexflag

    if (present(localDECount)) then
       call ESMF_DistGridGet(lstypep%segDistgrid, delayout=delayout, rc=localrc)
       if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
       call ESMF_DELayoutGet(delayout, localDeCount=localDECount, rc=localrc)
       if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
    endif

    if (present(name)) then
       call ESMF_GetName(lstypep%base, name, localrc)
       if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
    endif

    if (present(rc)) rc = ESMF_SUCCESS

  end subroutine ESMF_LineStreamGet


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamGetBounds"
!BOP
! !IROUTINE: ESMF_LineStreamGetBounds - Get DE-local segment index bounds
!
! !INTERFACE:
  subroutine ESMF_LineStreamGetBounds(linestream, keywordEnforcer, &
       localDE, exclusiveLBound, exclusiveUBound, rc)
!
! !ARGUMENTS:
    type(ESMF_LineStream), intent(in)            :: linestream
    type(ESMF_KeywordEnforcer), optional:: keywordEnforcer
    integer,               intent(in),  optional :: localDE
    integer,               intent(out), optional :: exclusiveLBound
    integer,               intent(out), optional :: exclusiveUBound
    integer,               intent(out), optional :: rc
!
! !DESCRIPTION:
!   Get the DE-local exclusive bounds of the segment index space.
!EOP

    type(ESMF_LineStreamType), pointer :: lstypep
    integer :: localrc, lDE

    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    ESMF_INIT_CHECK_DEEP(ESMF_LineStreamGetInit, linestream, rc)

    lstypep => linestream%lstypep

    lDE = 0
    if (present(localDE)) lDE = localDE

    if (present(exclusiveLBound)) then
       call c_ESMC_locstreamgetelbnd(lstypep%segDistgrid, lDE, lstypep%indexflag, &
            exclusiveLBound, localrc)
       if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
    endif

    if (present(exclusiveUBound)) then
       call c_ESMC_locstreamgeteubnd(lstypep%segDistgrid, lDE, lstypep%indexflag, &
            exclusiveUBound, localrc)
       if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
    endif

    if (present(rc)) rc = ESMF_SUCCESS

  end subroutine ESMF_LineStreamGetBounds


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamGetNodeCoords"
!BOP
! !IROUTINE: ESMF_LineStreamGetNodeCoords - Get node coordinate Arrays
!
! !INTERFACE:
  subroutine ESMF_LineStreamGetNodeCoords(linestream, coordDimIdx, &
       coordArray, rc)
!
! !ARGUMENTS:
    type(ESMF_LineStream), intent(in)  :: linestream
    integer,               intent(in)  :: coordDimIdx   ! 1=X/Lon, 2=Y/Lat, 3=Z/Radius
    type(ESMF_Array),      intent(out) :: coordArray
    integer,               intent(out), optional :: rc
!
! !DESCRIPTION:
!   Get a reference to the node coordinate Array for the specified
!   coordinate dimension. The returned Array is built on the node
!   DistGrid and has one entry per owned node.
!EOP

    type(ESMF_LineStreamType), pointer :: lstypep
    integer :: localrc

    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    ESMF_INIT_CHECK_DEEP(ESMF_LineStreamGetInit, linestream, rc)

    lstypep => linestream%lstypep

    if (coordDimIdx < 1 .or. coordDimIdx > lstypep%coordDim) then
       call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_OUTOFRANGE, &
            msg="coordDimIdx out of range", &
            ESMF_CONTEXT, rcToReturn=rc)
       return
    endif

    coordArray = lstypep%nodeCoordsArrays(coordDimIdx)

    if (present(rc)) rc = ESMF_SUCCESS

  end subroutine ESMF_LineStreamGetNodeCoords


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamGetSegConn"
!BOP
! !IROUTINE: ESMF_LineStreamGetSegConn - Get segment connectivity Array
!
! !INTERFACE:
  subroutine ESMF_LineStreamGetSegConn(linestream, startNodeArray, &
       endNodeArray, rc)
!
! !ARGUMENTS:
    type(ESMF_LineStream), intent(in)  :: linestream
    type(ESMF_Array),      intent(out), optional :: startNodeArray
    type(ESMF_Array),      intent(out), optional :: endNodeArray
    integer,               intent(out), optional :: rc
!
! !DESCRIPTION:
!   Get references to the segment connectivity Arrays. Each is a 1D
!   integer Array of shape (localSegCount) on the segment DistGrid,
!   holding 1-based global node IDs.
!EOP

    type(ESMF_LineStreamType), pointer :: lstypep

    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    ESMF_INIT_CHECK_DEEP(ESMF_LineStreamGetInit, linestream, rc)

    lstypep => linestream%lstypep
    if (present(startNodeArray)) startNodeArray = lstypep%segStartNodeArray
    if (present(endNodeArray))   endNodeArray   = lstypep%segEndNodeArray

    if (present(rc)) rc = ESMF_SUCCESS

  end subroutine ESMF_LineStreamGetSegConn


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamGetSegLength"
!BOP
! !IROUTINE: ESMF_LineStreamGetSegLength - Get segment length Array
!
! !INTERFACE:
  subroutine ESMF_LineStreamGetSegLength(linestream, lengthArray, rc)
!
! !ARGUMENTS:
    type(ESMF_LineStream), intent(in)  :: linestream
    type(ESMF_Array),      intent(out) :: lengthArray
    integer,               intent(out), optional :: rc
!
! !DESCRIPTION:
!   Get a reference to the segment length Array. One R8 value per segment.
!EOP

    type(ESMF_LineStreamType), pointer :: lstypep

    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    ESMF_INIT_CHECK_DEEP(ESMF_LineStreamGetInit, linestream, rc)

    lstypep => linestream%lstypep
    lengthArray = lstypep%segLengthArray

    if (present(rc)) rc = ESMF_SUCCESS

  end subroutine ESMF_LineStreamGetSegLength


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamDestroy"
!BOP
! !IROUTINE: ESMF_LineStreamDestroy - Release resources held by a LineStream
!
! !INTERFACE:
  subroutine ESMF_LineStreamDestroy(linestream, keywordEnforcer, &
       noGarbage, rc)
!
! !ARGUMENTS:
    type(ESMF_LineStream),      intent(inout)         :: linestream
    type(ESMF_KeywordEnforcer), optional:: keywordEnforcer
    logical,                    intent(in),  optional :: noGarbage
    integer,                    intent(out), optional :: rc
!
! !DESCRIPTION:
!   Destroys an {\tt ESMF\_LineStream} object and releases all associated
!   resources.
!EOP

    integer :: localrc

    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    ESMF_INIT_CHECK_DEEP(ESMF_LineStreamGetInit, linestream, rc)

    ! Destruct internals
    call ESMF_LineStreamDestruct(linestream%lstypep, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    ! Deallocate the type pointer
    deallocate(linestream%lstypep)
    nullify(linestream%lstypep)

    ESMF_INIT_SET_DELETED(linestream)

    if (present(rc)) rc = ESMF_SUCCESS

  end subroutine ESMF_LineStreamDestroy


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamDestruct"
  subroutine ESMF_LineStreamDestruct(lstypep, rc)
    type(ESMF_LineStreamType), pointer :: lstypep
    integer, intent(out), optional     :: rc

    integer :: localrc, i
    type(ESMF_Status) :: status

    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    call ESMF_BaseGetStatus(lstypep%base, status, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    if (status .eq. ESMF_STATUS_READY) then

       ! Destroy node coordinate Arrays
       if (associated(lstypep%nodeCoordsArrays)) then
          do i = 1, lstypep%coordDim
             if (lstypep%destroyNodeCoords(i)) then
                call ESMF_ArrayDestroy(lstypep%nodeCoordsArrays(i), rc=localrc)
                if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
                     ESMF_CONTEXT, rcToReturn=rc)) return
             endif
          enddo
          deallocate(lstypep%nodeCoordsArrays)
          deallocate(lstypep%destroyNodeCoords)
       endif

       ! Destroy segment connectivity Arrays
       if (lstypep%destroySegStart) then
          call ESMF_ArrayDestroy(lstypep%segStartNodeArray, rc=localrc)
          if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
               ESMF_CONTEXT, rcToReturn=rc)) return
       endif
       if (lstypep%destroySegEnd) then
          call ESMF_ArrayDestroy(lstypep%segEndNodeArray, rc=localrc)
          if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
               ESMF_CONTEXT, rcToReturn=rc)) return
       endif

       ! Destroy segment length Array
       if (lstypep%destroySegLength) then
          call ESMF_ArrayDestroy(lstypep%segLengthArray, rc=localrc)
          if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
               ESMF_CONTEXT, rcToReturn=rc)) return
       endif

       ! Destroy segment mask Array
       if (lstypep%destroySegMask) then
          call ESMF_ArrayDestroy(lstypep%segMaskArray, rc=localrc)
          if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
               ESMF_CONTEXT, rcToReturn=rc)) return
       endif

       ! Destroy DistGrids
       if (lstypep%destroySegDistgrid) then
          call ESMF_DistGridDestroy(lstypep%segDistgrid, rc=localrc)
          if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
               ESMF_CONTEXT, rcToReturn=rc)) return
       endif

       if (lstypep%destroyNodeDistgrid) then
          call ESMF_DistGridDestroy(lstypep%nodeDistgrid, rc=localrc)
          if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
               ESMF_CONTEXT, rcToReturn=rc)) return
       endif

       ! Deallocate ghost node arrays
       if (associated(lstypep%ghostNodeCoordsX))  deallocate(lstypep%ghostNodeCoordsX)
       if (associated(lstypep%ghostNodeCoordsY))  deallocate(lstypep%ghostNodeCoordsY)
       if (associated(lstypep%ghostNodeCoordsZ))  deallocate(lstypep%ghostNodeCoordsZ)
       if (associated(lstypep%ghostNodeGlobalIds)) deallocate(lstypep%ghostNodeGlobalIds)

       ! Deallocate node-to-segment adjacency
       if (associated(lstypep%nodeSegPtr)) deallocate(lstypep%nodeSegPtr)
       if (associated(lstypep%nodeSegIdx)) deallocate(lstypep%nodeSegIdx)

    endif

    ! Mark object invalid
    call ESMF_BaseSetStatus(lstypep%base, ESMF_STATUS_INVALID, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    if (present(rc)) rc = ESMF_SUCCESS

  end subroutine ESMF_LineStreamDestruct


!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamPrint"
!BOP
! !IROUTINE: ESMF_LineStreamPrint - Print LineStream information
!
! !INTERFACE:
  subroutine ESMF_LineStreamPrint(linestream, keywordEnforcer, rc)
!
! !ARGUMENTS:
    type(ESMF_LineStream),      intent(in)            :: linestream
    type(ESMF_KeywordEnforcer), optional:: keywordEnforcer
    integer,                    intent(out), optional :: rc
!
! !DESCRIPTION:
!   Print diagnostic information about a LineStream.
!EOP

    type(ESMF_LineStreamType), pointer :: lstypep
    character(len=ESMF_MAXSTR) :: nameStr

    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    ESMF_INIT_CHECK_DEEP(ESMF_LineStreamGetInit, linestream, rc)

    lstypep => linestream%lstypep

    call ESMF_GetName(lstypep%base, nameStr, rc=rc)

    write(*, '(A)')      '=============================='
    write(*, '(A,A)')    '  LineStream: ', trim(nameStr)
    write(*, '(A,I10)')  '  Total nodes:    ', lstypep%totalNodeCount
    write(*, '(A,I10)')  '  Total segments: ', lstypep%totalSegCount
    write(*, '(A,I10)')  '  Local nodes:    ', lstypep%localNodeCount
    write(*, '(A,I10)')  '  Local segments: ', lstypep%localSegCount
    write(*, '(A,I10)')  '  Ghost nodes:    ', lstypep%ghostNodeCount
    write(*, '(A,I6)')   '  Coord dim:      ', lstypep%coordDim
    write(*, '(A)')      '=============================='

    if (present(rc)) rc = ESMF_SUCCESS

  end subroutine ESMF_LineStreamPrint


!------------------------------------------------------------------------------
! Private helper: Haversine distance between two points in degrees
!------------------------------------------------------------------------------
  subroutine haversine_distance(lat1_deg, lon1_deg, lat2_deg, lon2_deg, dist_m)
    real(ESMF_KIND_R8), intent(in)  :: lat1_deg, lon1_deg, lat2_deg, lon2_deg
    real(ESMF_KIND_R8), intent(out) :: dist_m

    real(ESMF_KIND_R8) :: lat1, lon1, lat2, lon2
    real(ESMF_KIND_R8) :: dlat, dlon, a, c
    real(ESMF_KIND_R8), parameter :: DEG2RAD = 3.14159265358979323846d0 / 180.0d0
    real(ESMF_KIND_R8), parameter :: R = 6371000.0d0  ! Earth radius in meters

    lat1 = lat1_deg * DEG2RAD
    lon1 = lon1_deg * DEG2RAD
    lat2 = lat2_deg * DEG2RAD
    lon2 = lon2_deg * DEG2RAD

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    a = sin(dlat/2.0d0)**2 + cos(lat1) * cos(lat2) * sin(dlon/2.0d0)**2
    c = 2.0d0 * atan2(sqrt(a), sqrt(1.0d0 - a))

    dist_m = R * c

  end subroutine haversine_distance


!------------------------------------------------------------------------------
! Standardized initialization support
!------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamGetInit"
  function ESMF_LineStreamGetInit(d)
    type(ESMF_LineStream), intent(in), optional :: d
    ESMF_INIT_TYPE :: ESMF_LineStreamGetInit

    if (present(d)) then
      ESMF_LineStreamGetInit = ESMF_INIT_GET(d)
    else
      ESMF_LineStreamGetInit = ESMF_INIT_CREATED
    endif

  end function ESMF_LineStreamGetInit

end module ESMF_LineStreamMod
