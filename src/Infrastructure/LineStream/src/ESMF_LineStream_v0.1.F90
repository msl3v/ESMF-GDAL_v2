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
  ! class method.
  !
  !------------------------------------------------------------------------------
  ! INCLUDES
#include "ESMF.h"
  !------------------------------------------------------------------------------
  !
  !BOPI
  ! !MODULE: ESMF_LineStreamMod - Combine physical field metadata, data and grid
  !
  ! !DESCRIPTION:
  ! The code in this file implements the {\tt ESMF\_LineStream} class, which 
  ! represents a single scalar or vector field.  {\tt ESMF\_LineStream}s associate
  ! a metadata description expressed as a set of {\tt ESMF\_Attributes} with
  ! a data {\tt ESMF\_Array} and an {\tt ESMF\_Grid}.
  ! 
  ! This type is implemented in Fortran 90.
  !
  !------------------------------------------------------------------------------
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
  use ESMF_GridMod
  use ESMF_GridUtilMod
  use ESMF_RHandleMod
  use ESMF_StaggerLocMod
  use ESMF_ArrayMod
  use ESMF_ArrayBundleMod
  use ESMF_ArrayCreateMod
  use ESMF_ArrayGetMod
  use ESMF_InitMacrosMod
  use ESMF_MeshMod
  use ESMF_IOScripMod
  use ESMF_IOUGridMod

#ifdef ESMF_GDAL
  use iso_c_binding
#endif
  implicit none

  !------------------------------------------------------------------------------
  ! !PRIVATE TYPES:
  private

  !------------------------------------------------------------------------------
  ! ! ESMF_LineStreamType
  ! ! Definition of the LineStream class.


  type ESMF_LineStreamType
#ifndef ESMF_NO_SEQUENCE
     sequence
#endif

     !private
     type(ESMF_Base)                      :: base             ! base class object
     logical                              :: destroyDistgrid 
     type(ESMF_DistGrid)                  :: distgrid         ! description of index space of Arrays
     type(ESMF_Index_Flag)                :: indexflag
     type(ESMF_CoordSys_Flag)             :: coordSys
     integer                              :: numNodes         ! Number of unique local nodes
     integer                              :: numSegments      ! Number of local line segments
     type(ESMF_Array), pointer            :: nodeCoords(:,:)  ! (nDims,numNodes)
     type(ESMF_Array), pointer            :: segments(:,:)    ! (2,numSegments); point to indices in nodeCoords
!<>     integer                              :: keyCount         ! Number of keys
!<>     character(len=ESMF_MAXSTR), pointer  :: keyNames(:)      ! Names
!<>     character(len=ESMF_MAXSTR), pointer  :: keyUnits(:)      ! Units
!<>     character(len=ESMF_MAXSTR), pointer  :: keyLongNames(:)  ! Long names
!<>     logical, pointer                     :: destroyKeys(:)   ! if we're responsible for destroying key array
!<>     type (ESMF_Array), pointer           :: keys(:)          ! Contents
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
     !private       
     type (ESMF_LineStreamType), pointer :: lstypep
     ESMF_INIT_DECLARE
  end type ESMF_LineStream

  !------------------------------------------------------------------------------
  ! !PUBLIC TYPES:
  public ESMF_LineStream
  public ESMF_LineStreamType ! For internal use only

  public operator(==)
  public operator(/=)

!<>  public ESMF_LineStreamIsCreated
!<>  public ESMF_LineStreamValidate           ! Check internal consistency
!<>  public ESMF_LineStreamCreate
!<>  public ESMF_LineStreamGet
!<>  public ESMF_LineStreamGetBounds
!<>  public ESMF_LineStreamDeserialize
!<>  public ESMF_LineStreamSerialize
!<>  public ESMF_LineStreamDestroy
!<>  public ESMF_LineStreamDestruct           ! for ESMF garbage collection
!<>  public ESMF_LineStreamPrint              ! Print contents of a LineStream
!<>  public ESMF_LineStreamGetKey
!<>  public ESMF_LineStreamAddKey
!<>  public ESMF_LineStreamMatch
!<>  ! ====== MSL: Network connectivity methods ======
!<>  public ESMF_LineStreamGetConnectivity
!<>  public ESMF_LineStreamGetEdge
!<>  public ESMF_LineStreamGetNeighbors
!<>  public ESMF_LineStreamGetBoundaryFlags
!<>  public ESMF_LineStreamAddGhostEdges
!<>  public ESMF_LineStreamGetParmetisData
  ! ===============================================

  ! - ESMF-internal methods:
!<>  public ESMF_LineStreamTypeGetInit        ! For Standardized Initialization
!<>  public ESMF_LineStreamTypeInit           ! For Standardized Initialization
!<>  public ESMF_LineStreamTypeValidate
!<>  public ESMF_LineStreamGetInit            ! For Standardized Initialization

  !EOPI

  ! !PRIVATE MEMBER FUNCTIONS:

  !------------------------------------------------------------------------------
  ! The following line turns the CVS identifier string into a printable variable.
  character(*), parameter, private :: version = &
       '$Id$'

  !==============================================================================
  !
  ! INTERFACE BLOCKS
  !
  !==============================================================================

  ! -------------------------- ESMF-public method -------------------------------
  !BOPI
  ! !IROUTINE: ESMF_LineStreamCreate -- Generic interface

  ! !INTERFACE:
  interface ESMF_LineStreamCreate

     ! !PRIVATE MEMBER FUNCTIONS:
     !
     module procedure ESMF_LineStreamCreateFromDG
     module procedure ESMF_LineStreamCreateFromLocal
!<>     module procedure ESMF_LineStreamCreateFromNewDG
!<>     module procedure ESMF_LineStreamCreateReg
!<>     module procedure ESMF_LineStreamCreateIrreg
!<>     module procedure ESMF_LineStreamCreateByBkgMesh
!<>     module procedure ESMF_LineStreamCreateByBkgGrid
     module procedure ESMF_LineStreamCreateFromFile
!<>     module procedure ESMF_LineStreamCreateFromParmetis

     ! !DESCRIPTION: 
     ! This interface provides a single entry point for the various 
     !  types of {\tt ESMF\_LineStreamCreate} functions.   
     !EOPI 
  end interface ESMF_LineStreamCreate


  ! -------------------------- ESMF-public method -------------------------------
  !BOPI
  ! !IROUTINE: ESMF_LineStreamGetKey -- Generic interface

  ! !INTERFACE:
  interface ESMF_LineStreamGetKey

     ! !PRIVATE MEMBER FUNCTIONS:
     !
!<>     module procedure ESMF_LineStreamGetKeyI4
!<>     module procedure ESMF_LineStreamGetKeyR4
!<>     module procedure ESMF_LineStreamGetKeyR8
!<>     module procedure ESMF_LineStreamGetKeyArray  
!<>     module procedure ESMF_LineStreamGetKeyInfo

     ! !DESCRIPTION: 
     ! This interface provides a single entry point for the various 
     !  types of {\tt ESMF\_LineStreamGetKey} functions.   
     !EOPI 
  end interface ESMF_LineStreamGetKey

  ! -------------------------- ESMF-public method -------------------------------
  !BOPI
  ! !IROUTINE: ESMF_LineStreamAddKey -- Generic interface

  ! !INTERFACE:
  interface ESMF_LineStreamAddKey

     ! !PRIVATE MEMBER FUNCTIONS:
     !
!<>     module procedure ESMF_LineStreamAddKeyAlloc
!<>     module procedure ESMF_LineStreamAddKeyArray
!<>     module procedure ESMF_LineStreamAddKeyI4
!<>     module procedure ESMF_LineStreamAddKeyR4
!<>     module procedure ESMF_LineStreamAddKeyR8


     ! !DESCRIPTION: 
     ! This interface provides a single entry point for the various 
     !  types of {\tt ESMF\_LineStreamAddKey} functions.   
     !EOPI 
  end interface ESMF_LineStreamAddKey

  !===============================================================================
  ! LineStreamOperator() interfaces
  !===============================================================================

  ! -------------------------- ESMF-public method -------------------------------
  !BOP
  ! !IROUTINE: ESMF_LineStreamAssignment(=) - LineStream assignment
  !
  ! !INTERFACE:
  !   interface assignment(=)
  !   locstream1 = locstream2
  !
  ! !ARGUMENTS:
  !   type(ESMF_LineStream) :: locstream1
  !   type(ESMF_LineStream) :: locstream2
  !
  !
  ! !STATUS:
  ! \begin{itemize}
  ! \item\apiStatusCompatibleVersion{5.2.0r}
  ! \end{itemize}
  !
  ! !DESCRIPTION:
  !   Assign locstream1 as an alias to the same ESMF LineStream object in memory
  !   as locstream2. If locstream2 is invalid, then locstream1 will be equally invalid after
  !   the assignment.
  !
  !   The arguments are:
  !   \begin{description}
  !   \item[locstream1]
  !     The {\tt ESMF\_LineStream} object on the left hand side of the assignment.
  !   \item[locstream2]
  !     The {\tt ESMF\_LineStream} object on the right hand side of the assignment.
  !   \end{description}
  !
  !EOP
  !------------------------------------------------------------------------------


  ! -------------------------- ESMF-public method -------------------------------
  !BOP
  ! !IROUTINE: ESMF_LineStreamOperator(==) - LineStream equality operator
  !
  ! !INTERFACE:
  interface operator(==)
     !   if (locstream1 == locstream2) then ... endif
     !             OR
     !   result = (locstream1 == locstream2)
     ! !RETURN VALUE:
     !   logical :: result
     !
     ! !ARGUMENTS:
     !   type(ESMF_LineStream), intent(in) :: locstream1
     !   type(ESMF_LineStream), intent(in) :: locstream2
     !
     !
     ! !STATUS:
     ! \begin{itemize}
     ! \item\apiStatusCompatibleVersion{5.2.0r}
     ! \end{itemize}
     !
     ! !DESCRIPTION:
     !   Test whether locstream1 and locstream2 are valid aliases to the same ESMF
     !   LineStream object in memory. For a more general comparison of two ESMF LineStreams,
     !   going beyond the simple alias test, the ESMF\_LineStreamMatch() function (not yet
     !   implemented) must be used.
     !
     !   The arguments are:
     !   \begin{description}
     !   \item[locstream1]
     !     The {\tt ESMF\_LineStream} object on the left hand side of the equality
     !     operation.
     !   \item[locstream2]
     !     The {\tt ESMF\_LineStream} object on the right hand side of the equality
     !     operation.
     !   \end{description}
     !
     !EOP
!<>     module procedure ESMF_LineStreamEQ

  end interface operator(==)
  !------------------------------------------------------------------------------


  ! -------------------------- ESMF-public method -------------------------------
  !BOP
  ! !IROUTINE: ESMF_LineStreamOperator(/=) - LineStream not equal operator
  !
  ! !INTERFACE:
  interface operator(/=)
     !   if (locstream1 /= locstream2) then ... endif
     !             OR
     !   result = (locstream1 /= locstream2)
     ! !RETURN VALUE:
     !   logical :: result
     !
     ! !ARGUMENTS:
     !   type(ESMF_LineStream), intent(in) :: locstream1
     !   type(ESMF_LineStream), intent(in) :: locstream2
     !
     !
     ! !STATUS:
     ! \begin{itemize}
     ! \item\apiStatusCompatibleVersion{5.2.0r}
     ! \end{itemize}
     !
     ! !DESCRIPTION:
     !   Test whether locstream1 and locstream2 are {\it not} valid aliases to the
     !   same ESMF LineStream object in memory. For a more general comparison of two ESMF
     !   LineStreams, going beyond the simple alias test, the ESMF\_LineStreamMatch() function
     !   (not yet implemented) must be used.
     !
     !   The arguments are:
     !   \begin{description}
     !   \item[locstream1]
     !     The {\tt ESMF\_LineStream} object on the left hand side of the non-equality
     !     operation.
     !   \item[locstream2]
     !     The {\tt ESMF\_LineStream} object on the right hand side of the non-equality
     !     operation.
     !   \end{description}
     !
     !EOP
!<>     module procedure ESMF_LineStreamNE

  end interface operator(/=)
  !------------------------------------------------------------------------------

! C interface
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
    
    subroutine free_parmetis_graph_f(nodedist, xadj, adjncy, node_x, node_y) &
         bind(C, name="free_parmetis_graph_f")
      use iso_c_binding
      type(c_ptr), value, intent(in) :: nodedist, xadj, adjncy, node_x, node_y
    end subroutine

!    subroutine c_esmc_gdal_shpreadcoords(filename, local_pet, pet_count, &
!         localpoints, totaldims, cx_ptr, cy_ptr, rc) &
!         bind(C, name="c_esmc_gdal_shpreadcoords")
!      use iso_c_binding
!      character(kind=c_char), dimension(*), intent(in) :: filename
!      integer(c_int), intent(in)  :: local_pet, pet_count
!      integer(c_int), intent(out) :: localpoints, totaldims
!      type(c_ptr), intent(out)    :: cx_ptr, cy_ptr
!      integer(c_int), intent(out) :: rc
!    end subroutine
!
    subroutine c_esmc_gdal_shpfreecoords(cx_ptr, cy_ptr) &
         bind(C, name="c_esmc_gdal_shpfreecoords_c")
      use iso_c_binding
      type(c_ptr), intent(inout) :: cx_ptr, cy_ptr
    end subroutine
  end interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamCreateFromFile"
  !BOP
  ! !IROUTINE: ESMF_LineStreamCreate - Create a new LineStream from a grid file
  !\label{locstream:createfromfile}
  ! !INTERFACE:
  ! Private name: call using ESMF_LineStreamCreate()
  function ESMF_LineStreamCreateFromFile(filename, keywordEnforcer, &
       fileformat, varname, indexflag, centerflag, name, rc)
    !
    ! !RETURN VALUE:
    type(ESMF_LineStream) :: ESMF_LineStreamCreateFromFile

    !
    ! !ARGUMENTS:
    character (len=*),          intent(in)           :: filename
    type(ESMF_KeywordEnforcer), optional:: keywordEnforcer ! must use keywords below
    type(ESMF_FileFormat_Flag), intent(in), optional :: fileformat
    character(len=*),           intent(in), optional :: varname
    type(ESMF_Index_Flag),      intent(in), optional :: indexflag
    logical,                    intent(in), optional :: centerflag
    character (len=*),          intent(in), optional :: name
    integer,                    intent(out),optional :: rc

    ! For GDAL/METIS
#ifdef ESMF_GDAL
    type(c_ptr) :: nodedist_ptr, xadj_ptr, adjncy_ptr, node_x_ptr, node_y_ptr
    integer(c_int32_t), pointer :: nodedist(:), xadj(:), adjncy(:)
    real(c_double), pointer :: node_x(:), node_y(:)
    integer(c_int) :: nnodes, nedges
    integer :: c_ierr
    character(kind=c_char, len=256) :: c_filename
    integer :: mpi_comm
    real(c_double) :: tol_val
#endif

    ! !DESCRIPTION:
    !     Create a new {\tt ESMF\_LineStream} object and add the coordinate keys and mask key
    !     to the LineStream using the coordinates defined in a grid file.  Currently, it 
    !     supports the SCRIP format, the ESMF unstructured grid format and the UGRID format.
    !     For a 2D or 3D grid in ESMF or UGRID format, it can construct the LineStream using either 
    !     the center coordinates or the corner coordinates.  For a SCRIP format grid file, the
    !     LineStream can only be constructed using the center coordinates.  In
    !     addition, it supports 1D network topology in UGRID format.  When
    !     construction a LineStream using a 1D UGRID, it always uses node
    !     coordinates (i.e., corner coordinates). 
    !
    !     The arguments are:
    !     \begin{description}
    !     \item[filename]
    !          Name of grid file to be used to create the location stream.  
    !     \item[{[fileformat]}]
    !     The file format.  The valid options are {\tt ESMF\_FILEFORMAT\_SCRIP},
    !     {\tt ESMF\_FILEFORMAT\_ESMFMESH}, and {\tt ESMF\_FILEFORMAT\_UGRID}.
    !      Please see section~\ref{const:fileformatflag} for a detailed description of the options.
    !     If not specified, the default is {\tt ESMF\_FILEFORMAT\_SCRIP}.
    !     \item[{[varname]}]
    !         An optional variable name stored in the UGRID file to be used to
    !         generate the mask using the missing value of the data value of
    !         this variable.  The first two dimensions of the variable has to be the
    !         the longitude and the latitude dimension and the mask is derived from the
    !         first 2D values of this variable even if this data is 3D, or 4D array. If not 
    !         specified, no mask is used for a UGRID file.
    !     \item[{[indexflag]}]
    !          Flag that indicates how the DE-local indices are to be defined.
    !          Defaults to {\tt ESMF\_INDEX\_DELOCAL}, which indicates
    !          that the index range on each DE starts at 1. See Section~\ref{const:indexflag}
    !          for the full range of options. 
    !     \item[{[centerflag]}]
    !          Flag that indicates whether to use the center coordinates to construct the location stream.
    !          If true, use center coordinates, otherwise, use the corner coordinates.  If not specified,
    !          use center coordinates as default.  For SCRIP files, only center coordinate 
    !          is supported.
    !     \item[{[name]}]
    !          Name of the location stream
    !     \item[{[rc]}]
    !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
    !   \end{description}
    !
    !EOP

#ifdef ESMF_NETCDF
    integer :: totalpoints,totaldims
    type(ESMF_VM) :: vm
    integer :: numDim, buf(1), msgbuf(3)
    integer :: localrc
    integer :: PetNo, PetCnt
    type(ESMF_Index_Flag) :: indexflagLocal
    real(ESMF_KIND_R8), pointer :: coordX(:), coordY(:), coordZ(:)
    real(ESMF_KIND_R8), pointer :: coord2D(:,:), varbuffer(:)
    integer(ESMF_KIND_I4), pointer :: imask(:)
    integer :: starti, count, localcount, index
    integer :: remain, i, j
    integer :: meshid
    real(ESMF_KIND_R8) :: missingvalue
    type(ESMF_CoordSys_Flag) :: coordSys
    type(ESMF_LineStream) :: locStream
    type(ESMF_FileFormat_Flag) :: localfileformat
    logical :: localcenterflag, haveface
    character(len=16) :: units, location
#ifdef ESMF_GDAL
    integer :: maxid, minid
    integer(c_int), pointer :: gFIDs(:)
    integer, allocatable :: FIDs(:)
    ! C pointer returns from c_esmc_gdal_shpreadcoords
    type(c_ptr) :: cx_ptr, cy_ptr
    real(ESMF_KIND_R8), pointer :: c_coordX(:), c_coordY(:)
    ! Variables for shapefile connectivity
    type(ESMF_LineStreamType), pointer :: lstypep
#endif

    if (present(indexflag)) then
       indexflagLocal=indexflag
    else
       indexflagLocal=ESMF_INDEX_DELOCAL
    endif

    if (present(fileformat)) then
       localfileformat = fileformat
    else
       localfileformat = ESMF_FILEFORMAT_SCRIP
    endif

    if (present(centerflag)) then
       localcenterflag = centerflag
    else
       localcenterflag = .TRUE.
    endif

    if (localcenterflag) then
       location = 'face'
    else
       location = 'node'
    endif
    if (localfileformat == ESMF_FILEFORMAT_GRIDSPEC) then
       call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, &
            msg="Create LineStream from a GRIDSPEC file is not supported.", &
            ESMF_CONTEXT, rcToReturn=rc)
       return
    endif

    if ( localfileformat == ESMF_FILEFORMAT_SCRIP .and. &
         .NOT. localcenterflag) then
       call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, &
            msg="Only allow center coordinates if the file is in SCRIP format", &
            ESMF_CONTEXT, rcToReturn=rc)
       return
    endif

#if 0
    if (localfileformat /= ESMF_FILEFORMAT_UGRID .and. &
         (present(meshname) .or. present(varname))) then
       call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_WRONG, &
            msg="Only UGRID file need the optional arguments meshname or varname", &
            ESMF_CONTEXT, rcToReturn=rc)
       return
    endif
#endif

    ! Initialize return code; assume failure until success is certain
    localrc = ESMF_RC_NOT_IMPL
    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    ! get global vm information
    !
    call ESMF_VMGetCurrent(vm, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    ! set up local pet info
    call ESMF_VMGet(vm, localPet=PetNo, petCount=PetCnt, rc=localrc)
    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    if (localfileformat == ESMF_FILEFORMAT_SHAPEFILE) then
       ! Read shapefile using GDAL - TRUE single-pass coordinate extraction
       ! One call: C reads features once with dynamic realloc, returns pointers.
       ! For connectivity/network topology, use ESMF_LineStreamCreateFromShapefile()

       call c_esmc_gdal_shpreadcoords(trim(filename)//C_NULL_CHAR, PetNo, PetCnt, &
            localcount, totaldims, cx_ptr, cy_ptr, localrc)
       if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return

       ! Set up dimensions
       totalpoints = localcount
       coordSys = ESMF_COORDSYS_SPH_RAD

       ! Map C pointers to Fortran arrays, copy, then free C memory
       allocate(coordX(localcount), coordY(localcount), imask(localcount))

       if (localcount > 0) then
          call c_f_pointer(cx_ptr, c_coordX, [localcount])
          call c_f_pointer(cy_ptr, c_coordY, [localcount])
          coordX(1:localcount) = c_coordX(1:localcount)
          coordY(1:localcount) = c_coordY(1:localcount)
       endif

       call c_esmc_gdal_shpfreecoords(cx_ptr, cy_ptr)

       ! Set all points as valid (mask = 1)
       imask(:) = 1
       
    endif
    
    ! create Line Stream

    locStream = ESMF_LineStreamCreate(name=name, localcount=localcount, indexflag=indexflagLocal,&
         coordSys = coordSys, rc=localrc)

    if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

!<>    ! Add coordinate keys based on coordSys
!<>    if ((coordSys == ESMF_COORDSYS_SPH_DEG) .or. (coordSys == ESMF_COORDSYS_SPH_RAD)) then 
!<>
!<>       call ESMF_LineStreamAddKey(locStream, 'ESMF:Lon',coordX, keyUnits=units, &
!<>            keyLongName='Longitude', &
!<>            datacopyflag=ESMF_DATACOPY_VALUE, rc=localrc)
!<>       if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
!<>            ESMF_CONTEXT, rcToReturn=rc)) return
!<>
!<>       call ESMF_LineStreamAddKey(locStream, 'ESMF:Lat',coordY, keyUnits=units, &
!<>            keyLongName='Latitude', &
!<>            datacopyflag=ESMF_DATACOPY_VALUE, rc=localrc)
!<>       if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
!<>            ESMF_CONTEXT, rcToReturn=rc)) return
!<>
!<>       !If 3D grid, add the height coordinates
!<>       if (totaldims == 3) then
!<>          if (localcount == 0) allocate(coordZ(localcount))
!<>          call ESMF_LineStreamAddKey(locStream, 'ESMF:Radius',coordZ, &
!<>               keyUnits='radius', &
!<>               keyLongName='Height', &
!<>               datacopyflag=ESMF_DATACOPY_VALUE, rc=localrc)
!<>          if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
!<>               ESMF_CONTEXT, rcToReturn=rc)) return
!<>          deallocate(coordZ)
!<>       endif
!<>
!<>    else if (coordSys == ESMF_COORDSYS_CART) then
!<>
!<>       call ESMF_LineStreamAddKey(locStream, 'ESMF:X',coordX, keyUnits=units, &
!<>            datacopyflag=ESMF_DATACOPY_VALUE, rc=localrc)
!<>       if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
!<>            ESMF_CONTEXT, rcToReturn=rc)) return
!<>
!<>       call ESMF_LineStreamAddKey(locStream, 'ESMF:Y',coordY, keyUnits=units, &
!<>            datacopyflag=ESMF_DATACOPY_VALUE, rc=localrc)
!<>       if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
!<>            ESMF_CONTEXT, rcToReturn=rc)) return
!<>
!<>       !If 3D grid, add the height coordinates
!<>       if (totaldims == 3) then
!<>          if (localcount == 0) allocate(coordZ(localcount))
!<>          call ESMF_LineStreamAddKey(locStream, 'ESMF:Z',coordZ, &
!<>               datacopyflag=ESMF_DATACOPY_VALUE, rc=localrc)
!<>          if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
!<>               ESMF_CONTEXT, rcToReturn=rc)) return
!<>          deallocate(coordZ)
!<>       endif
!<>    else
!<>       call ESMF_LogSetError(rcToCheck=ESMF_RC_ARG_BAD, &
!<>            msg="Unrecognized coordinate system.", &
!<>            ESMF_CONTEXT, rcToReturn=rc)
!<>       return
!<>    endif
!<>
!<>    !Add mask key
!<>    if (localfileformat .ne. ESMF_FILEFORMAT_SHAPEFILE) then
!<>       call ESMF_LineStreamAddKey(locStream, 'ESMF:Mask',imask,  &
!<>            keyLongName='Mask', &
!<>            datacopyflag=ESMF_DATACOPY_VALUE, rc=localrc)
!<>       if (ESMF_LogFoundError(localrc, ESMF_ERR_PASSTHRU, &
!<>            ESMF_CONTEXT, rcToReturn=rc)) return
!<>    endif

    ! local garbage collection
    deallocate(coordX, coordY, imask)

    ESMF_LineStreamCreateFromFile = locStream

    if (present(rc)) rc=ESMF_SUCCESS
    return

#else
    if (present(rc)) rc = ESMF_RC_LIB_NOT_PRESENT
#endif

  end function ESMF_LineStreamCreateFromFile

  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamCreate"
  !BOP
  ! !IROUTINE: ESMF_LineStreamCreate - Create a new LineStream from a local count

  ! !INTERFACE:
  ! Private name: call using ESMF_LineStreamCreate()
  function ESMF_LineStreamCreateFromLocal(localCount, keywordEnforcer, &
       indexflag, coordSys, name, rc)
    !
    ! !RETURN VALUE:
    type(ESMF_LineStream) :: ESMF_LineStreamCreateFromLocal

    !
    ! !ARGUMENTS:
    integer, intent(in)                             :: localCount
    type(ESMF_KeywordEnforcer), optional            :: keywordEnforcer ! must use keywords below
    type(ESMF_Index_Flag), intent(in), optional     :: indexflag
    type(ESMF_CoordSys_Flag), intent(in),  optional :: coordSys
    character (len=*), intent(in), optional         :: name
    integer, intent(out), optional                  :: rc
    !
    ! !DESCRIPTION:
    !     Allocates memory for a new {\tt ESMF\_LineStream} object, constructs its
    !     internal derived types.  The {\tt ESMF\_DistGrid} is set up, indicating
    !     how the LineStream is distributed. The assumed layout is one DE per PET.
    !
    !     The arguments are:
    !     \begin{description}
    !     \item[localCount]
    !          Number of grid cells to be distributed to this DE/PET.
    !     \item[{[indexflag]}]
    !          Flag that indicates how the DE-local indices are to be defined.
    !          Defaults to {\tt ESMF\_INDEX\_DELOCAL}, which indicates
    !          that the index range on each DE starts at 1. See Section~\ref{const:indexflag}
    !          for the full range of options. 
    !     \item[{[coordSys]}]
    !         The coordinate system of the location stream coordinate data.
    !         For a full list of options, please see Section~\ref{const:coordsys}.
    !         If not specified then defaults to ESMF\_COORDSYS\_SPH\_DEG.
    !     \item[{[name]}]
    !          Name of the location stream
    !     \item[{[rc]}]
    !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
    !   \end{description}
    !
    !EOP

    integer                                               :: localrc  ! Error status
    type(ESMF_VM)                                   :: vm       ! Virtual machine used
    integer, allocatable  :: countsPerPet(:)
    integer :: localPet, petCount
    integer :: i, currMin
    type(ESMF_DistGrid)                 :: distgrid
    integer, pointer :: deBLockList(:,:,:)   
    integer               :: minIndex(1), maxIndex(1)
    type(ESMF_Index_Flag)  :: indexflagLocal
    type(ESMF_CoordSys_Flag) :: coordSysLocal


    ! Initialize return code; assume failure until success is certain
    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    ! Set defaults
    if (present(indexflag)) then
       indexflagLocal=indexflag
    else
       indexflagLocal=ESMF_INDEX_DELOCAL
    endif

    if (present(coordSys)) then
       coordSysLocal=coordSys
    else
       coordSysLocal=ESMF_COORDSYS_SPH_DEG
    endif

    ! Get VM for this context
    call ESMF_VMGetCurrent(vm, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
         ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    ! Gather localCount for each Pet
    call ESMF_VMGet( vm, localPet = localPet,                        &
         petCount = petCount, rc=localrc )
    if (ESMF_LogFoundError(localrc, &
         ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    allocate(countsPerPet(petCount), stat=localrc)
    if (ESMF_LogFoundAllocError(localrc, msg="Allocating countsPerPet", &
         ESMF_CONTEXT, rcToReturn=rc)) return

    call ESMF_VMAllGather(vm, sendData=(/localCount/),               &
         recvData=countsPerPet, count=1, rc=localrc)
    if (ESMF_LogFoundError(localrc, &
         ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    ! Setup DistGrid
    !! define min and maxIndex
    minIndex(1)=1
    maxIndex(1)=sum(countsPerPet(:))

    !! setup deBlockList
    allocate(deBlockList(1,2,petCount), stat=localrc)
    if (ESMF_LogFoundAllocError(localrc, msg="Allocating deBlockList", &
         ESMF_CONTEXT, rcToReturn=rc)) return
    currMin=1
    do i=1,petCount
       deBlockList(1,1,i)=currMin    
       deBlockList(1,2,i)=currMin+countsPerPet(i)-1
       currMin=deBlockList(1,2,i)+1
    enddo

    !! Create DistGrid
    distgrid=ESMF_DistGridCreate(minIndex=minIndex, &
         maxIndex=maxIndex, &
         deBlockList=deBlockList, &
         indexflag=indexflagLocal, &
         rc=localrc)
    if (ESMF_LogFoundError(localrc, &
         ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    ! cleanup local allocations
    deallocate(countsPerPet)
    deallocate(deBlockList)

    ! Create LineStream using CreateFromDistGrid version
    ESMF_LineStreamCreateFromLocal=ESMF_LineStreamCreateFromDG(name=name, &
         distgrid=distgrid, &
         indexflag=indexflagLocal, &
         coordSys=coordSysLocal, &
         rc=localrc )
    if (ESMF_LogFoundError(localrc, &
         ESMF_ERR_PASSTHRU, &
         ESMF_CONTEXT, rcToReturn=rc)) return

    ! Set distgrid to be destroyed, since ESMF created it
    ESMF_LineStreamCreateFromLocal%lstypep%destroyDistgrid=.true.

    ! return successfully
    if (present(rc)) rc = ESMF_SUCCESS

  end function ESMF_LineStreamCreateFromLocal
  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
#undef  ESMF_METHOD
#define ESMF_METHOD "ESMF_LineStreamCreate"
  !BOP
  ! !IROUTINE: ESMF_LineStreamCreate - Create a new LineStream from a distgrid

  ! !INTERFACE:
  ! Private name: call using ESMF_LineStreamCreate()
  function ESMF_LineStreamCreateFromDG(distgrid, keywordEnforcer, &
       indexflag, coordSys, name, vm, rc )
    !
    ! !RETURN VALUE:
    type(ESMF_LineStream) :: ESMF_LineStreamCreateFromDG

    !
    ! !ARGUMENTS:
    type(ESMF_DistGrid),      intent(in)            :: distgrid
    type(ESMF_KeywordEnforcer), optional:: keywordEnforcer ! must use keywords below
    type(ESMF_Index_Flag),    intent(in),  optional :: indexflag    
    type(ESMF_CoordSys_Flag), intent(in),  optional :: coordSys
    character (len=*),        intent(in),  optional :: name
    type(ESMF_VM),            intent(in),  optional :: vm
    integer,                  intent(out), optional :: rc
    !
    ! !DESCRIPTION:
    !     Allocates memory for a new {\tt ESMF\_LineStream} object, constructs its
    !     internal derived types. 
    !
    !     The arguments are:
    !     \begin{description}
    !     \item[distgrid]
    !          Distgrid specifying size and distribution. Only 1D distgrids are allowed.
    !     \item[{[indexflag]}]
    !          Flag that indicates how the DE-local indices are to be defined.
    !          Defaults to {\tt ESMF\_INDEX\_DELOCAL}, which indicates
    !          that the index range on each DE starts at 1. See Section~\ref{const:indexflag}
    !          for the full range of options. 
    !     \item[{[coordSys]}]
    !         The coordinate system of the location stream coordinate data.
    !         For a full list of options, please see Section~\ref{const:coordsys}.
    !         If not specified then defaults to ESMF\_COORDSYS\_SPH\_DEG.
    !     \item[{[name]}]
    !          Name of the location stream
    !     \item[{[vm]}]
    !         If present, the LineStream object is created on the specified 
    !         {\tt ESMF\_VM} object. The default is to create on the VM of the 
    !         current context.
    !     \item[{[rc]}]
    !          Return code; equals {\tt ESMF\_SUCCESS} if there are no errors.
    !   \end{description}
    !
    !EOP

    integer                             :: localrc  ! Error status
    type (ESMF_LineStreamType), pointer  :: lstypep
    type(ESMF_LineStream)                :: locstream 
    integer                             :: dimCount 
    type(ESMF_Index_Flag)               :: indexflagLocal
    type(ESMF_CoordSys_Flag)            :: coordSysLocal
    type(ESMF_Pointer)                  :: vmThis
    logical                             :: actualFlag

    ! Initialize return code; assume failure until success is certain
    if (present(rc)) rc = ESMF_RC_NOT_IMPL

    ! Init check input types
    ESMF_INIT_CHECK_DEEP_SHORT(ESMF_DistGridGetInit,distgrid,rc)      

    ! Must make sure the local PET is associated with an actual member
    actualFlag = .true.
    if (present(vm)) then
       call ESMF_VMGetThis(vm, vmThis)
       if (vmThis == ESMF_NULL_POINTER) then
          actualFlag = .false.  ! local PET is not for an actual member of Array
       endif
    endif

    if (actualFlag) then
       ! only actual member PETs worry about the DistGrid

       ! Make sure DistGrid is 1D
       call ESMF_DistGridGet(distgrid, dimCount=dimCount, rc=localrc)
       if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return
       if (dimCount .ne. 1) then
          if (ESMF_LogFoundError(ESMF_RC_ARG_RANK, &
               msg=" - DistGrid must be 1D", &
               ESMF_CONTEXT, rcToReturn=rc)) return
       endif

    endif

    ! Initialize pointers
    nullify(lstypep)
    nullify(ESMF_LineStreamCreateFromDG%lstypep)

    ! allocate LineStream type
    allocate(lstypep, stat=localrc)
    if (ESMF_LogFoundAllocError(localrc, msg="Allocating LineStream type object", &
         ESMF_CONTEXT, rcToReturn=rc)) return

    ! Initialize key member variables
!    nullify(lstypep%keyNames)
!    nullify(lstypep%keyUnits)
!    nullify(lstypep%keyLongNames)
!    nullify(lstypep%keys)
!    nullify(lstypep%destroyKeys)

    ! Set defaults
    if (present(indexflag)) then
       indexflagLocal=indexflag
    else
       indexflagLocal=ESMF_INDEX_DELOCAL
    endif

    if (present(coordSys)) then
       coordSysLocal=coordSys
    else
       coordSysLocal=ESMF_COORDSYS_SPH_DEG
    endif

    ! Set some remaining info into the struct      
    lstypep%indexflag=indexflagLocal
    lstypep%coordSys=coordSysLocal
    lstypep%destroyDistgrid=.false.
!    lstypep%keyCount=0
!    lstypep%has_connectivity=.false.

    if (actualFlag) then
       ! only actual member PETs set distgrid
       lstypep%distgrid=distgrid

       ! create base object and set name
       call ESMF_BaseCreate(lstypep%base,"LineStream",name,0,rc=localrc)       
       if (ESMF_LogFoundError(localrc, &
            ESMF_ERR_PASSTHRU, &
            ESMF_CONTEXT, rcToReturn=rc)) return

       ! Set pointer to internal locstream type
       locstream%lstypep=>lstypep

       ! Set return value.
       ESMF_LineStreamCreateFromDG=locstream

       ! Add reference to this object into ESMF garbage collection table
       ! Only call this in those Create() methods that do not call other LSCreate()
       call c_ESMC_VMAddFObject(locstream, &
            ESMF_ID_LOCSTREAM%objectID)

    endif

    ! set init status to created
    ESMF_INIT_SET_CREATED(ESMF_LineStreamCreateFromDG)

    ! return successfully
    if (present(rc)) rc = ESMF_SUCCESS

  end function ESMF_LineStreamCreateFromDG
  !------------------------------------------------------------------------------

end module ESMF_LineStreamMod
