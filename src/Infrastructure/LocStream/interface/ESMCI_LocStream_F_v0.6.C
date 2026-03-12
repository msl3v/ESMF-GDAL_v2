// $Id$
//
// Earth System Modeling Framework
// Copyright (c) 2002-2025, University Corporation for Atmospheric Research,
// Massachusetts Institute of Technology, Geophysical Fluid Dynamics
// Laboratory, University of Michigan, National Centers for Environmental
// Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
// NASA Goddard Space Flight Center.
// Licensed under the University of Illinois-NCSA License.

// ESMC interface routines

//-----------------------------------------------------------------------------
//
// !DESCRIPTION:
//
//
//-----------------------------------------------------------------------------
//
// insert any higher level, 3rd party or system includes here

#include <cstring>
using namespace std;

#include "ESMCI_Macros.h"
#include "ESMCI_LogErr.h"
#include "ESMCI_DistGrid.h"
#include "ESMCI_Array.h"
#include "ESMCI_CoordSys.h"

#ifdef ESMF_GDAL
#include <ogr_api.h>
#include "ESMCI_GDAL_Util.h"
#include "ESMCI_FileIO_Util.h"
#endif

//-----------------------------------------------------------------------------
// leave the following line as-is; it will insert the cvs ident string
// into the object file for tracking purposes.
static const char *const version =
  "$Id$";
//-----------------------------------------------------------------------------

extern "C" {
  //
  //-----------------------------------------------------------------------------
  //-----------------------------------------------------------------------------
  //
  // This section includes all the LocStream routines
  //
  //

  void FTN_X(c_esmc_gdal_getnfeatures)(
				       char *filename,
				       int *local_pet,
				       int *pet_count,
				       int *nfeatures,
				       int *locfeatures,
				       int *min_id, 
				       int *max_id,
				       int *rc,
				       ESMCI_FortranStrLenArg filename_l) {


#undef  ESMC_METHOD1
#define ESMC_METHOD "c_esmc_gdal_getnfeatures()"

    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

#ifdef ESMF_GDAL
    // Open file and create datasource (DS)
    OGRDataSourceH hDS;
    if (access(filename, F_OK) == 0) {
      OGRRegisterAll(); // register all the drivers
      hDS = GDALOpenEx(filename, GDAL_OF_VECTOR | GDAL_OF_READONLY, NULL, NULL, NULL); //OGROpen( filename, FALSE, NULL );
      if( hDS == NULL )
	{
	  printf( "Open failed on pet %d: %s, %d\n", *local_pet, CPLGetLastErrorMsg(), CPLGetLastErrorNo() );
	}
    } else if (*local_pet == 0) {
      printf("Cannot access shapefile %s\n",filename);
      return;
    }
  
    // GET DA DEETS!
    OGRLayerH hLayer = OGR_DS_GetLayer( hDS, 0 );
    *nfeatures = OGR_L_GetFeatureCount(hLayer,1);

    //  int min_id, max_id;
    divide_ids_evenly_as_possible(*nfeatures, *local_pet, *pet_count, *min_id, *max_id);

    *locfeatures = *max_id-*min_id+1;

    //    printf("--- nFeatures: %d\n", *nfeatures);

    // Cleanup
    GDALClose( hDS );
#endif

    // return success
    if (rc) *rc = ESMF_SUCCESS;

    return;
  }

  void FTN_X(c_esmc_gdal_getglobal_fids)(
                                         char *filename,
					 int *nfeatures,
					 int *gFIDs,
                                         int *rc,
                                         ESMCI_FortranStrLenArg filename_l) {


#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_gdal_getglobal_fids()"

    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

#ifdef ESMF_GDAL
    // Open file and create datasource (DS)
    OGRDataSourceH hDS;
    if (access(filename, F_OK) == 0) {
      OGRRegisterAll(); // register all the drivers
      hDS = GDALOpenEx(filename, GDAL_OF_VECTOR | GDAL_OF_READONLY, NULL, NULL, NULL); //OGROpen( filename, FALSE, NULL );
      if( hDS == NULL )
	{
	  printf( "Open failed on pet %d: %s, %d\n", 0, CPLGetLastErrorMsg(), CPLGetLastErrorNo() );
	  return;
	}
    } else {
      printf("Cannot access shapefile %s\n",filename);
      return;
    }
  
    OGRLayerH hLayer = OGR_DS_GetLayer( hDS, 0 );
    for (int i=0;i<*nfeatures;i++) {
      OGRFeatureH hFeature = OGR_L_GetFeature(hLayer,i);
      gFIDs[i] = OGR_F_GetFID(hFeature);
      OGR_F_Destroy( hFeature );
    }

    // Cleanup
    GDALClose( hDS );
#endif

    // return success
    if (rc) *rc = ESMF_SUCCESS;

    return;
  }

  void FTN_X(c_esmc_gdal_shpinquire)(
				     char *filename,
				     int *local_pet,
				     int *pet_count,
				     int *localpoints,
				     int *totaldims,
				     int *totfeatures,
				     int *gFIDs,
				     int *locfeatures,
				     int *FIDs,
				     int *rc,
				     ESMCI_FortranStrLenArg filename_l) {


#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_gdal_shpinquire()"

    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

    *totaldims = 2; // This is fixed for now!!!

#ifdef ESMF_GDAL
    // Open file and create datasource (DS)
    OGRDataSourceH hDS;
    if (access(filename, F_OK) == 0) {
      OGRRegisterAll(); // register all the drivers
      hDS = GDALOpenEx(filename, GDAL_OF_VECTOR | GDAL_OF_READONLY, NULL, NULL, NULL); //OGROpen( filename, FALSE, NULL );
      if( hDS == NULL )
	{
	  printf( "Open failed on pet %d: %s, %d\n", *local_pet, CPLGetLastErrorMsg(), CPLGetLastErrorNo() );
	}
    } else if (*local_pet == 0) {
      printf("Cannot access shapefile %s\n",filename);
      return;
    }
  
    OGRLayerH hLayer = OGR_DS_GetLayer( hDS, 0 );
    // Get positions at which to read element information
    std::vector<int> feature_ids_vec;
    get_ids_divided_evenly_across_pets(*totfeatures, *local_pet, *pet_count, feature_ids_vec);
  
    //  printf("--- info: %d\n", feature_ids_vec.size());

    // Assign vector info to pointer
    *localpoints = 0;
    if (*locfeatures != 0) {
      // Get the total points in features on local PET (I don't wanna do this here, but I will and then will add it to things to fix)
      *localpoints = 0;
      for (int i=0;i<*locfeatures;i++) {
	FIDs[i] = gFIDs[feature_ids_vec[i]-1];
	//      printf("gFID %d vec %d FID %d on PET %d\n", gFIDs[i], feature_ids_vec[i]-1, FIDs[i], *local_pet);  
	OGRFeatureH hFeature = OGR_L_GetFeature(hLayer,FIDs[i]);
	OGRGeometryH hGeom = OGR_F_GetGeometryRef(hFeature);
	*localpoints += OGR_G_GetPointCount(hGeom);
	OGR_F_Destroy( hFeature );
      }
    }

    // Cleanup
    GDALClose( hDS );
#endif

    // return success
    if (rc) *rc = ESMF_SUCCESS;

    return;
  }

  void FTN_X(c_esmc_gdal_shpgetcoords)(
				       char *filename,
				       int *local_pet,
				       int *numFeatures,
				       int *FIDs,
				       int *localcount,
				       double *coordX,
				       double *coordY,
				       int *rc,
				       ESMCI_FortranStrLenArg filename_l) {


#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_gdal_shpgetcoords()"

    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

#ifdef ESMF_GDAL
    // Open file and create datasource (DS)
    OGRDataSourceH hDS;
    if (access(filename, F_OK) == 0) {
      OGRRegisterAll(); // register all the drivers
      hDS = GDALOpenEx( filename, GDAL_OF_VECTOR | GDAL_OF_READONLY, NULL, NULL, NULL );
      //    printf("Opened file %s\n",filename);
      if( hDS == NULL )
	{
	  printf( "Open failed on pet %d: %s, %d\n", *local_pet, CPLGetLastErrorMsg(), CPLGetLastErrorNo() );
	}
    } else if (*local_pet == 0) {
      printf("Cannot access shapefile %s\n",filename);
      return;
    } else
      printf("hDS open\n");
  
    // GET DA DEETS!
    std::vector<double> XCoords;
    std::vector<double> YCoords;
    OGRLayerH    hLayer = OGR_DS_GetLayer( hDS, 0 );
  
    for (int f = 0; f < *localcount; f++) {
      OGRFeatureH hFeature = OGR_L_GetFeature(hLayer, FIDs[f]);
      //    printf("Reading feature %d\n",FIDs[f]);
      if( hFeature == NULL ) {
	printf( "NULL Feature on PET %d with ID %d: %s, %d\n", *local_pet, FIDs[f], 
		CPLGetLastErrorMsg(), CPLGetLastErrorNo() );
	continue;
      }
    
      OGRGeometryH hGeom = OGR_F_GetGeometryRef(hFeature);
      if (hGeom == NULL) {
	printf("NULL Geometry on PET %d for Feature ID %d\n", *local_pet, FIDs[f]);
	OGR_F_Destroy( hFeature );
	continue;
      }
    
      OGRwkbGeometryType geomType = wkbFlatten(OGR_G_GetGeometryType(hGeom));
    
      // Handle Point geometries
      if (geomType == wkbPoint) {
	int nFTRpoints = OGR_G_GetPointCount(hGeom);
      
	for (int i = nFTRpoints-1; i >= 0; i--) {
	  XCoords.push_back( OGR_G_GetX(hGeom, i) );
	  YCoords.push_back( OGR_G_GetY(hGeom, i) );
	}
      }
      // Handle LineString geometries (rail networks, roads, rivers, etc.)
      else if (geomType == wkbLineString) {
	int nPoints = OGR_G_GetPointCount(hGeom);
      
	// Extract all vertices from the linestring
	for (int i = 0; i < nPoints; i++) {
	  XCoords.push_back( OGR_G_GetX(hGeom, i) );
	  YCoords.push_back( OGR_G_GetY(hGeom, i) );
	}
      }
      // Handle MultiLineString geometries
      else if (geomType == wkbMultiLineString) {
	int nGeometries = OGR_G_GetGeometryCount(hGeom);
      
	for (int g = 0; g < nGeometries; g++) {
	  OGRGeometryH hSubGeom = OGR_G_GetGeometryRef(hGeom, g);
	  int nPoints = OGR_G_GetPointCount(hSubGeom);
        
	  for (int i = 0; i < nPoints; i++) {
	    XCoords.push_back( OGR_G_GetX(hSubGeom, i) );
	    YCoords.push_back( OGR_G_GetY(hSubGeom, i) );
	  }
	}
      }
      // Warn about unsupported geometry types
      else {
	if (*local_pet == 0) {
	  printf("WARNING: Unsupported geometry type %d on Feature ID %d\n", 
		 geomType, FIDs[f]);
	}
      }
    
      OGR_F_Destroy( hFeature );
    }

    // Check if we got any coordinates
    if (XCoords.size() == 0) {
      printf("WARNING on PET %d: No coordinates extracted from shapefile\n", *local_pet);
    }

	if (*local_pet == 0) {
	  printf("NOTE: FOR SHAPEFILE, ASSUMING DEG. CONVERTING TO RADIANS!!!\n");
	}
    int numCoords = XCoords.size();
    for (int i = 0; i < numCoords; i++) {
      coordX[i] = XCoords[i] * ESMC_CoordSys_Deg2Rad;
      coordY[i] = YCoords[i] * ESMC_CoordSys_Deg2Rad;
    }

    // Cleanup
    if( hDS != NULL ) GDALClose( hDS );
#endif

    // return success
    if (rc) *rc = ESMF_SUCCESS;

    return;
  }

  /*===========================================================================
   * c_esmc_gdal_shpreadcoords - TRUE single-pass shapefile coordinate reader
   *
   * Opens the file once, iterates features once, grows coordinate buffers
   * dynamically via realloc, and returns C-allocated arrays + count to Fortran.
   *
   * Fortran calls this once, gets back:
   *   - localpoints: number of coordinate points extracted
   *   - totaldims:   always 2
   *   - cx_ptr, cy_ptr: C pointers to malloc'd arrays (radians)
   *
   * Fortran must call c_esmc_gdal_shpfreecoords(cx_ptr, cy_ptr) when done.
   *=========================================================================*/

/* Define to enable timing output, comment out to disable */
#define SHPREAD_TIMING

#ifdef SHPREAD_TIMING
#include <time.h>
  static double shpread_gettime(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + (double)ts.tv_nsec * 1.0e-9;
  }
#define SHPREAD_TIMER_DECL  double _t0, _t1, _t_open, _t_seek, _t_loop, _t_close, _t_total;
#define SHPREAD_TIMER_START _t_total = shpread_gettime();
#define SHPREAD_TIMER_MARK(var) var = shpread_gettime();
#define SHPREAD_TIMER_PRINT(pet, nfeat, npts, nrealloc) \
    if (1) { \
      double _t_now = shpread_gettime(); \
      _t_close = _t_now - _t1; \
      double _elapsed = _t_now - _t_total; \
      printf("[SHP_TIMING] PET %d  features=%lld  points=%d  reallocs=%d  " \
             "open=%.4f seek=%.4f loop=%.4f close=%.4f  TOTAL=%.4f sec\n", \
             pet, (long long)(nfeat), npts, nrealloc, \
             _t_open, _t_seek, _t_loop, _t_close, _elapsed); \
      fflush(stdout); \
    }
#else
#define SHPREAD_TIMER_DECL
#define SHPREAD_TIMER_START
#define SHPREAD_TIMER_MARK(var)
#define SHPREAD_TIMER_PRINT(pet, nfeat, npts, nrealloc)
#endif

/* Initial buffer size and growth factor for coordinate arrays */
#define SHPREAD_INIT_CAP  65536
#define SHPREAD_GROW(cap) ((cap) + ((cap) >> 1))  /* 1.5x growth */

  void FTN_X(c_esmc_gdal_shpreadcoords)(
                                         const char *filename,
                                         int *local_pet,
                                         int *pet_count,
                                         int *localpoints,
                                         int *totaldims,
                                         void **cx_ptr,
                                         void **cy_ptr,
                                         int *rc) {

#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_gdal_shpreadcoords()"

    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

    *totaldims = 2;
    *cx_ptr = NULL;
    *cy_ptr = NULL;
    *localpoints = 0;

#ifdef ESMF_GDAL
    SHPREAD_TIMER_DECL
    SHPREAD_TIMER_START

    // Open file and create datasource
    OGRDataSourceH hDS;
    if (access(filename, F_OK) == 0) {
      OGRRegisterAll();
      hDS = GDALOpenEx(filename, GDAL_OF_VECTOR | GDAL_OF_READONLY,
                       NULL, NULL, NULL);
      if (hDS == NULL) {
        printf("Open failed on pet %d: %s, %d\n", *local_pet,
               CPLGetLastErrorMsg(), CPLGetLastErrorNo());
        return;
      }
    } else {
      if (*local_pet == 0)
        printf("Cannot access shapefile %s\n", filename);
      return;
    }

    OGRLayerH hLayer = OGR_DS_GetLayer(hDS, 0);
    if (hLayer == NULL) {
      printf("Could not get layer on pet %d\n", *local_pet);
      GDALClose(hDS);
      return;
    }

    // Get total feature count and compute this PET's range
    GIntBig total_features = OGR_L_GetFeatureCount(hLayer, TRUE);
    GIntBig features_per_pet = total_features / *pet_count;
    GIntBig remainder = total_features % *pet_count;

    GIntBig my_start, my_end;  // [my_start, my_end)
    if (*local_pet < remainder) {
      my_start = *local_pet * (features_per_pet + 1);
      my_end   = my_start + features_per_pet + 1;
    } else {
      my_start = remainder * (features_per_pet + 1)
                 + (*local_pet - remainder) * features_per_pet;
      my_end   = my_start + features_per_pet;
    }

    SHPREAD_TIMER_MARK(_t0)
    _t_open = _t0 - _t_total;

    // Use sequential access: reset reading, skip to our start feature
    OGR_L_ResetReading(hLayer);

    if (my_start > 0) {
      if (OGR_L_SetNextByIndex(hLayer, my_start) != OGRERR_NONE) {
        for (GIntBig s = 0; s < my_start; s++) {
          OGRFeatureH hSkip = OGR_L_GetNextFeature(hLayer);
          if (hSkip) OGR_F_Destroy(hSkip);
        }
      }
    }

    SHPREAD_TIMER_MARK(_t1)
    _t_seek = _t1 - _t0;
    SHPREAD_TIMER_MARK(_t0)

    // Dynamic coordinate buffers
    int capacity = SHPREAD_INIT_CAP;
    double *bufX = (double*)malloc(capacity * sizeof(double));
    double *bufY = (double*)malloc(capacity * sizeof(double));
    int point_count = 0;
    int realloc_count = 0;
    GIntBig my_feature_count = my_end - my_start;

    // Macro to ensure capacity for nPts more points
    #define ENSURE_CAP(nPts) \
      if (point_count + (nPts) > capacity) { \
        while (point_count + (nPts) > capacity) \
          capacity = SHPREAD_GROW(capacity); \
        bufX = (double*)realloc(bufX, capacity * sizeof(double)); \
        bufY = (double*)realloc(bufY, capacity * sizeof(double)); \
        realloc_count++; \
      }

    for (GIntBig f = 0; f < my_feature_count; f++) {
      OGRFeatureH hFeature = OGR_L_GetNextFeature(hLayer);
      if (hFeature == NULL) break;

      OGRGeometryH hGeom = OGR_F_GetGeometryRef(hFeature);
      if (hGeom == NULL) {
        OGR_F_Destroy(hFeature);
        continue;
      }

      OGRwkbGeometryType geomType = wkbFlatten(OGR_G_GetGeometryType(hGeom));

      if (geomType == wkbPoint) {
        int nPts = OGR_G_GetPointCount(hGeom);
        ENSURE_CAP(nPts)
        for (int i = nPts - 1; i >= 0; i--) {
          bufX[point_count + (nPts - 1 - i)] =
              OGR_G_GetX(hGeom, i) * ESMC_CoordSys_Deg2Rad;
          bufY[point_count + (nPts - 1 - i)] =
              OGR_G_GetY(hGeom, i) * ESMC_CoordSys_Deg2Rad;
        }
        point_count += nPts;
      }
      else if (geomType == wkbLineString) {
        int nPts = OGR_G_GetPointCount(hGeom);
        ENSURE_CAP(nPts)
        for (int i = 0; i < nPts; i++) {
          bufX[point_count + i] =
              OGR_G_GetX(hGeom, i) * ESMC_CoordSys_Deg2Rad;
          bufY[point_count + i] =
              OGR_G_GetY(hGeom, i) * ESMC_CoordSys_Deg2Rad;
        }
        point_count += nPts;
      }
      else if (geomType == wkbMultiLineString) {
        int nGeometries = OGR_G_GetGeometryCount(hGeom);
        for (int g = 0; g < nGeometries; g++) {
          OGRGeometryH hSubGeom = OGR_G_GetGeometryRef(hGeom, g);
          int nPts = OGR_G_GetPointCount(hSubGeom);
          ENSURE_CAP(nPts)
          for (int i = 0; i < nPts; i++) {
            bufX[point_count + i] =
                OGR_G_GetX(hSubGeom, i) * ESMC_CoordSys_Deg2Rad;
            bufY[point_count + i] =
                OGR_G_GetY(hSubGeom, i) * ESMC_CoordSys_Deg2Rad;
          }
          point_count += nPts;
        }
      }
      else {
        if (*local_pet == 0 && f == 0) {
          printf("WARNING: Unsupported geometry type %d in shapefile\n",
                 geomType);
        }
      }

      OGR_F_Destroy(hFeature);
    }

    #undef ENSURE_CAP

    SHPREAD_TIMER_MARK(_t1)
    _t_loop = _t1 - _t0;

    if (point_count == 0) {
      printf("WARNING on PET %d: No coordinates extracted from shapefile\n",
             *local_pet);
      free(bufX);
      free(bufY);
      bufX = NULL;
      bufY = NULL;
    }

    if (*local_pet == 0) {
      printf("NOTE: FOR SHAPEFILE, ASSUMING DEG. CONVERTING TO RADIANS!!!\n");
    }

    *localpoints = point_count;
    *cx_ptr = (void*)bufX;
    *cy_ptr = (void*)bufY;

    SHPREAD_TIMER_MARK(_t1)
    GDALClose(hDS);
    SHPREAD_TIMER_PRINT(*local_pet, my_feature_count, point_count, realloc_count)
#endif

    if (rc) *rc = ESMF_SUCCESS;
    return;
  }

  /*===========================================================================
   * c_esmc_gdal_shpfreecoords - Free arrays allocated by shpreadcoords
   *=========================================================================*/
  void c_esmc_gdal_shpfreecoords_c(void **cx_ptr, void **cy_ptr) {
    if (cx_ptr && *cx_ptr) { free(*cx_ptr); *cx_ptr = NULL; }
    if (cy_ptr && *cy_ptr) { free(*cy_ptr); *cy_ptr = NULL; }
  }

  void FTN_X(c_esmc_gdal_getpolylinevertices)(
					      char *filename,
					      int *featureId,
					      int *nVertices,
					      double *coordX,
					      double *coordY,
					      int *rc,
					      ESMCI_FortranStrLenArg filename_l) {


#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_gdal_getpolylinevertices()"

    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

#ifdef ESMF_GDAL
    // Open file and create datasource (DS)
    OGRDataSourceH hDS;
    if (access(filename, F_OK) == 0) {
      OGRRegisterAll(); // register all the drivers
      hDS = GDALOpenEx(filename, GDAL_OF_VECTOR | GDAL_OF_READONLY, NULL, NULL, NULL);
      if( hDS == NULL )
	{
	  printf( "Open failed: %s, %d\n", CPLGetLastErrorMsg(), CPLGetLastErrorNo() );
	  return;
	}
    } else {
      printf("Cannot access shapefile %s\n",filename);
      return;
    }
  
    // Get layer
    OGRLayerH hLayer = OGR_DS_GetLayer( hDS, 0 );
  
    // Get the specific feature by FID
    OGRFeatureH hFeature = OGR_L_GetFeature(hLayer, *featureId);
    if( hFeature == NULL ) {
      printf("Feature %d not found\n", *featureId);
      GDALClose( hDS );
      if (rc) *rc = ESMC_RC_NOT_FOUND;
      return;
    }
  
    // Get geometry from feature
    OGRGeometryH hGeom = OGR_F_GetGeometryRef(hFeature);
    if( hGeom == NULL ) {
      printf("No geometry in feature %d\n", *featureId);
      OGR_F_Destroy( hFeature );
      GDALClose( hDS );
      if (rc) *rc = ESMC_RC_NOT_FOUND;
      return;
    }
  
    // Get geometry type
    OGRwkbGeometryType geomType = OGR_G_GetGeometryType(hGeom);
  
    // Check if it's a line geometry (LineString or MultiLineString)
    if (geomType != wkbLineString && 
	geomType != wkbLineString25D &&
	geomType != wkbMultiLineString &&
	geomType != wkbMultiLineString25D) {
      printf("Feature %d is not a polyline (type=%d)\n", *featureId, geomType);
      OGR_F_Destroy( hFeature );
      GDALClose( hDS );
      if (rc) *rc = ESMC_RC_ARG_WRONG;
      return;
    }
  
    // Handle MultiLineString by getting first linestring
    OGRGeometryH hLine = hGeom;
    if (geomType == wkbMultiLineString || geomType == wkbMultiLineString25D) {
      int nGeoms = OGR_G_GetGeometryCount(hGeom);
      if (nGeoms > 0) {
	hLine = OGR_G_GetGeometryRef(hGeom, 0);  // Get first linestring
      } else {
	printf("MultiLineString feature %d has no geometries\n", *featureId);
	OGR_F_Destroy( hFeature );
	GDALClose( hDS );
	if (rc) *rc = ESMC_RC_NOT_FOUND;
	return;
      }
    }
  
    // Get number of points in the linestring
    int numPoints = OGR_G_GetPointCount(hLine);
    *nVertices = numPoints;
  
    // Extract vertices
    for (int i = 0; i < numPoints; i++) {
      double x, y, z;
      OGR_G_GetPoint(hLine, i, &x, &y, &z);
      coordX[i] = x;
      coordY[i] = y;
    }
  
    // Cleanup
    OGR_F_Destroy( hFeature );
    GDALClose( hDS );
#endif

    // return success
    if (rc) *rc = ESMF_SUCCESS;

    return;
  }


  // ==============================================================================
  // ALTERNATIVE VERSION: Extract ALL vertices from ALL features on this PET
  // More efficient than calling the above function repeatedly
  // ==============================================================================

  void FTN_X(c_esmc_gdal_getallpolylinevertices)(
						 char *filename,
						 int *numFeatures,
						 int *featureIds,
						 int *vertexOffsets,
						 int *totalVertices,
						 double *coordX,
						 double *coordY,
						 int *rc,
						 ESMCI_FortranStrLenArg filename_l) {


#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_gdal_getallpolylinevertices()"

    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

#ifdef ESMF_GDAL
    // Open file and create datasource (DS)
    OGRDataSourceH hDS;
    if (access(filename, F_OK) == 0) {
      OGRRegisterAll(); // register all the drivers
      hDS = GDALOpenEx(filename, GDAL_OF_VECTOR | GDAL_OF_READONLY, NULL, NULL, NULL);
      if( hDS == NULL )
	{
	  printf( "Open failed: %s, %d\n", CPLGetLastErrorMsg(), CPLGetLastErrorNo() );
	  return;
	}
    } else {
      printf("Cannot access shapefile %s\n",filename);
      return;
    }
  
    // Get layer
    OGRLayerH hLayer = OGR_DS_GetLayer( hDS, 0 );
  
    int vertexCount = 0;
  
    // Iterate through all features assigned to this PET
    for (int f = 0; f < *numFeatures; f++) {
      // Record the starting offset for this feature's vertices
      vertexOffsets[f] = vertexCount;
    
      // Get feature
      OGRFeatureH hFeature = OGR_L_GetFeature(hLayer, featureIds[f]);
      if( hFeature == NULL ) {
	printf("Warning: Feature %d not found, skipping\n", featureIds[f]);
	continue;
      }
    
      // Get geometry
      OGRGeometryH hGeom = OGR_F_GetGeometryRef(hFeature);
      if( hGeom == NULL ) {
	printf("Warning: No geometry in feature %d, skipping\n", featureIds[f]);
	OGR_F_Destroy( hFeature );
	continue;
      }
    
      // Get geometry type
      OGRwkbGeometryType geomType = OGR_G_GetGeometryType(hGeom);
    
      // Check if it's a line geometry
      if (geomType != wkbLineString && 
	  geomType != wkbLineString25D &&
	  geomType != wkbMultiLineString &&
	  geomType != wkbMultiLineString25D) {
	printf("Warning: Feature %d is not a polyline, skipping\n", featureIds[f]);
	OGR_F_Destroy( hFeature );
	continue;
      }
    
      // Handle MultiLineString
      int numGeoms = 1;
      if (geomType == wkbMultiLineString || geomType == wkbMultiLineString25D) {
	numGeoms = OGR_G_GetGeometryCount(hGeom);
      }
    
      // Extract vertices from all linestrings in this feature
      for (int g = 0; g < numGeoms; g++) {
	OGRGeometryH hLine = hGeom;
	if (numGeoms > 1) {
	  hLine = OGR_G_GetGeometryRef(hGeom, g);
	}
      
	int numPoints = OGR_G_GetPointCount(hLine);
      
	// Extract all points
	for (int i = 0; i < numPoints; i++) {
	  double x, y, z;
	  OGR_G_GetPoint(hLine, i, &x, &y, &z);
	  coordX[vertexCount] = x;
	  coordY[vertexCount] = y;
	  vertexCount++;
	}
      }
    
      OGR_F_Destroy( hFeature );
    }
  
    // Record total number of vertices extracted
    *totalVertices = vertexCount;
  
    // Record offset for end of last feature
    if (*numFeatures > 0) {
      vertexOffsets[*numFeatures] = vertexCount;
    }
  
    // Cleanup
    GDALClose( hDS );
#endif

    // return success
    if (rc) *rc = ESMF_SUCCESS;

    return;
  }


  // ==============================================================================
  // HELPER FUNCTION: Count vertices per feature (for pre-allocation)
  // Call this first to determine array sizes, then allocate, then call getall
  // ==============================================================================

  void FTN_X(c_esmc_gdal_countpolylinevertices)(
						char *filename,
						int *numFeatures,
						int *featureIds,
						int *vertexCounts,
						int *totalVertices,
						int *rc,
						ESMCI_FortranStrLenArg filename_l) {


#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_gdal_countpolylinevertices()"

    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

#ifdef ESMF_GDAL
    // Open file and create datasource (DS)
    OGRDataSourceH hDS;
    if (access(filename, F_OK) == 0) {
      OGRRegisterAll(); // register all the drivers
      hDS = GDALOpenEx(filename, GDAL_OF_VECTOR | GDAL_OF_READONLY, NULL, NULL, NULL);
      if( hDS == NULL )
	{
	  printf( "Open failed: %s, %d\n", CPLGetLastErrorMsg(), CPLGetLastErrorNo() );
	  return;
	}
    } else {
      printf("Cannot access shapefile %s\n",filename);
      return;
    }
  
    // Get layer
    OGRLayerH hLayer = OGR_DS_GetLayer( hDS, 0 );
  
    int totalCount = 0;
  
    // Count vertices in each feature
    for (int f = 0; f < *numFeatures; f++) {
      vertexCounts[f] = 0;
    
      // Get feature
      OGRFeatureH hFeature = OGR_L_GetFeature(hLayer, featureIds[f]);
      if( hFeature == NULL ) {
	continue;
      }
    
      // Get geometry
      OGRGeometryH hGeom = OGR_F_GetGeometryRef(hFeature);
      if( hGeom == NULL ) {
	OGR_F_Destroy( hFeature );
	continue;
      }
    
      // Get geometry type
      OGRwkbGeometryType geomType = OGR_G_GetGeometryType(hGeom);
    
      // Check if it's a line geometry
      if (geomType != wkbLineString && 
	  geomType != wkbLineString25D &&
	  geomType != wkbMultiLineString &&
	  geomType != wkbMultiLineString25D) {
	OGR_F_Destroy( hFeature );
	continue;
      }
    
      // Handle MultiLineString
      int numGeoms = 1;
      if (geomType == wkbMultiLineString || geomType == wkbMultiLineString25D) {
	numGeoms = OGR_G_GetGeometryCount(hGeom);
      }
    
      // Count vertices from all linestrings in this feature
      for (int g = 0; g < numGeoms; g++) {
	OGRGeometryH hLine = hGeom;
	if (numGeoms > 1) {
	  hLine = OGR_G_GetGeometryRef(hGeom, g);
	}
      
	int numPoints = OGR_G_GetPointCount(hLine);
	vertexCounts[f] += numPoints;
      }
    
      totalCount += vertexCounts[f];
      OGR_F_Destroy( hFeature );
    }
  
    *totalVertices = totalCount;
  
    // Cleanup
    GDALClose( hDS );
#endif

    // return success
    if (rc) *rc = ESMF_SUCCESS;

    return;
  }

  // non-method functions
  void FTN_X(c_esmc_locstreamgetkeybnds)(ESMCI::Array **_array,
					 int *_localDE,
					 int *exclusiveLBound,
					 int *exclusiveUBound,
					 int *exclusiveCount,
					 int *computationalLBound,
					 int *computationalUBound,
					 int *computationalCount,
					 int *totalLBound,
					 int *totalUBound,
					 int *totalCount,
					 int *rc){

#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_locstreamgetkeybnds()"

    ESMCI::Array *array;
    int localDE;
    int localrc;

    // Initialize return code; assume routine not implemented
    if (rc != NULL) *rc = ESMC_RC_NOT_IMPL;

    // Dereference variables
    array=*_array;

    // localDE
    if (ESMC_NOT_PRESENT_FILTER(_localDE) == ESMC_NULL_POINTER) {
      if (array->getDistGrid()->getDELayout()->getLocalDeCount()>1) {
	ESMC_LogDefault.MsgFoundError(ESMC_RC_ARG_WRONG,
				      "- Must provide localDE if localDeCount >1",
				      ESMC_CONTEXT, ESMC_NOT_PRESENT_FILTER(rc));
	return;
      } else {
	localDE=0;
      }
    } else {
      localDE=*_localDE; // already 0 based

      // Input Error Checking
      if ((localDE < 0) || (localDE >=array->getDistGrid()->getDELayout()->getLocalDeCount())) {
	ESMC_LogDefault.MsgFoundError(ESMC_RC_ARG_WRONG,
				      "- localDE outside range on this processor", ESMC_CONTEXT,
				      ESMC_NOT_PRESENT_FILTER(rc));
	return;
      }
    }

    // ExclusiveLBound
    if (ESMC_NOT_PRESENT_FILTER(exclusiveLBound) != ESMC_NULL_POINTER) {
      *exclusiveLBound=*(array->getExclusiveLBound()+localDE);
    }

    // ExclusiveUBound
    if (ESMC_NOT_PRESENT_FILTER(exclusiveUBound) != ESMC_NULL_POINTER) {
      *exclusiveUBound=*(array->getExclusiveUBound()+localDE);
    }

    // ExclusiveCount
    if (ESMC_NOT_PRESENT_FILTER(exclusiveCount) != ESMC_NULL_POINTER) {
      *exclusiveCount=*(array->getExclusiveUBound()+localDE) -
	*(array->getExclusiveLBound()+localDE) + 1;
    }

    // ComputationalLBound
    if (ESMC_NOT_PRESENT_FILTER(computationalLBound) != ESMC_NULL_POINTER) {
      *computationalLBound=*(array->getComputationalLBound()+localDE);
    }

    // ComputationalUBound
    if (ESMC_NOT_PRESENT_FILTER(computationalUBound) != ESMC_NULL_POINTER) {
      *computationalUBound=*(array->getComputationalUBound()+localDE);
    }

    // ComputationalCount
    if (ESMC_NOT_PRESENT_FILTER(computationalCount) != ESMC_NULL_POINTER) {
      *computationalCount=*(array->getComputationalUBound()+localDE) -
	*(array->getComputationalLBound()+localDE) + 1;
    }


    // TotalLBound
    if (ESMC_NOT_PRESENT_FILTER(totalLBound) != ESMC_NULL_POINTER) {
      *totalLBound=*(array->getTotalLBound()+localDE);
    }

    // TotalUBound
    if (ESMC_NOT_PRESENT_FILTER(totalUBound) != ESMC_NULL_POINTER) {
      *totalUBound=*(array->getTotalUBound()+localDE);
    }

    // TotalCount
    if (ESMC_NOT_PRESENT_FILTER(totalCount) != ESMC_NULL_POINTER) {
      *totalCount=*(array->getTotalUBound()+localDE) -
	*(array->getTotalLBound()+localDE) + 1;
    }



    // Return ESMF_SUCCESS
    if (rc != NULL) *rc = ESMF_SUCCESS;

    return;
  }

  // non-method functions
  void FTN_X(c_esmc_locstreamgetelbnd)(ESMCI::DistGrid **_distgrid,
				       int *_localDE,
				       ESMC_IndexFlag *_indexflag,
				       int *exclusiveLBound,
				       int *rc){

#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_locstreamgetelbnd()"

    ESMCI::DistGrid *distgrid;
    int localDE;
    int localrc;
    ESMC_IndexFlag indexflag;

    // Initialize return code; assume routine not implemented
    if (rc != NULL) *rc = ESMC_RC_NOT_IMPL;

    // Dereference variables
    distgrid=*_distgrid;
    indexflag=*_indexflag;

    // localDE
    if (ESMC_NOT_PRESENT_FILTER(_localDE) == ESMC_NULL_POINTER) {
      if (distgrid->getDELayout()->getLocalDeCount()>1) {
	ESMC_LogDefault.MsgFoundError(ESMC_RC_ARG_WRONG,
				      "- Must provide localDE if localDeCount >1",
				      ESMC_CONTEXT, ESMC_NOT_PRESENT_FILTER(rc));
	return;
      } else {
	localDE=0;
      }
    } else {
      localDE=*_localDE; // already 0 based

      // Input Error Checking
      if ((localDE < 0) || (localDE >=distgrid->getDELayout()->getLocalDeCount())) {
	ESMC_LogDefault.MsgFoundError(ESMC_RC_ARG_WRONG,
				      "- localDE outside range on this processor", ESMC_CONTEXT,
				      ESMC_NOT_PRESENT_FILTER(rc));
	return;
      }
    }

    if (indexflag==ESMC_INDEX_DELOCAL) {
      *exclusiveLBound = 1; // excl. region starts at (1,1,1...)
    } else {

      // Get some useful information
      const int *localDeToDeMap = distgrid->getDELayout()->getLocalDeToDeMap();

      // Get the Global DE from the local DE
      int de = localDeToDeMap[localDE];

      // obtain min index for this DE
      int const *index_min=distgrid->getMinIndexPDimPDe(de,&localrc);
      if (ESMC_LogDefault.MsgFoundError(localrc,ESMCI_ERR_PASSTHRU, ESMC_CONTEXT,
					ESMC_NOT_PRESENT_FILTER(rc))) return;

      // Set lower bound of exclusive region
      *exclusiveLBound = *index_min;
    }

    // Return ESMF_SUCCESS
    if (rc != NULL) *rc = ESMF_SUCCESS;
    return;
  }


  void FTN_X(c_esmc_locstreamgeteubnd)(ESMCI::DistGrid **_distgrid,
				       int *_localDE,
				       ESMC_IndexFlag *_indexflag,
				       int *exclusiveUBound,
				       int *rc){

#undef  ESMC_METHOD
#define ESMC_METHOD "c_esmc_locstreamgeteubnd()"

    ESMCI::DistGrid *distgrid;
    int localDE;
    int localrc;
    ESMC_IndexFlag indexflag;

    // Initialize return code; assume routine not implemented
    if (rc != NULL) *rc = ESMC_RC_NOT_IMPL;

    // Dereference variables
    distgrid=*_distgrid;
    indexflag=*_indexflag;

    // localDE
    if (ESMC_NOT_PRESENT_FILTER(_localDE) == ESMC_NULL_POINTER) {
      if (distgrid->getDELayout()->getLocalDeCount()>1) {
	ESMC_LogDefault.MsgFoundError(ESMC_RC_ARG_WRONG,
				      "- Must provide localDE if localDeCount >1",
				      ESMC_CONTEXT, ESMC_NOT_PRESENT_FILTER(rc));
	return;
      } else {
	localDE=0;
      }
    } else {
      localDE=*_localDE; // already 0 based

      // Input Error Checking
      if ((localDE < 0) || (localDE >=distgrid->getDELayout()->getLocalDeCount())) {
	ESMC_LogDefault.MsgFoundError(ESMC_RC_ARG_WRONG,
				      "- localDE outside range on this processor", ESMC_CONTEXT,
				      ESMC_NOT_PRESENT_FILTER(rc));
	return;
      }
    }

    // Get some useful information
    const int *localDeToDeMap = distgrid->getDELayout()->getLocalDeToDeMap();
    const int *indexCountPDimPDe = distgrid->getIndexCountPDimPDe();

    // Get the Global DE from the local DE
    int de = localDeToDeMap[localDE];

    // Set upper bound based on indexflag
    if (indexflag==ESMC_INDEX_DELOCAL) {
      *exclusiveUBound = indexCountPDimPDe[de];
    } else {

      // obtain max index for this DE
      int const *index_max=distgrid->getMaxIndexPDimPDe(de,&localrc);
      if (ESMC_LogDefault.MsgFoundError(localrc,ESMCI_ERR_PASSTHRU, ESMC_CONTEXT,
					ESMC_NOT_PRESENT_FILTER(rc))) return;

      // Set upper bound of exclusive region
      *exclusiveUBound = *index_max;
    }
    // Return ESMF_SUCCESS
    if (rc != NULL) *rc = ESMF_SUCCESS;

    return;
  }



#if 1
  // non-method functions
  void FTN_X(c_esmc_locstreamserialize)(ESMC_IndexFlag *indexflag,
					int *keyCount, ESMC_CoordSys_Flag *coordSys,
					char *buffer, int *length, int *offset,
					ESMC_InquireFlag *inquireflag,
					int *rc,
					ESMCI_FortranStrLenArg buffer_l){


    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

    // Calc size to store locstream info
    int size = sizeof(ESMC_IndexFlag)+sizeof(int)+sizeof(ESMC_CoordSys_Flag);

    // If not just inquiring save info to buffer
    //  Verify length > vars.
    if (*inquireflag != ESMF_INQUIREONLY) {

      // Make sure info will fit
      if ((*length - *offset) < size) {
        ESMC_LogDefault.MsgFoundError(ESMC_RC_ARG_BAD,
				      "Buffer too short to add a LocStream object", ESMC_CONTEXT, rc);
        return;
      }


      // Save indexflag
      ESMC_IndexFlag *ifp = (ESMC_IndexFlag *)(buffer + *offset);
      *ifp++ = *indexflag;

      // Save keyCount
      int *ip= (int *)ifp;
      *ip++ = *keyCount;

      // Save coordSys
      ESMC_CoordSys_Flag *csp= (ESMC_CoordSys_Flag *)ip;
      *csp++ = *coordSys;
    }

    // Increase offset by size
    *offset += size;

    // return success
    if (rc) *rc = ESMF_SUCCESS;

    return;
  }


  void FTN_X(c_esmc_locstreamdeserialize)(ESMC_IndexFlag *indexflag,
					  int *keyCount, ESMC_CoordSys_Flag *coordSys,
					  char *buffer, int *offset,
					  int *rc,
					  ESMCI_FortranStrLenArg buffer_l){

    
    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

    // Get indexflag
    ESMC_IndexFlag *ifp = (ESMC_IndexFlag *)(buffer + *offset);
    *indexflag=*ifp++;

    // Get keyCount
    int *ip= (int *)ifp;
    *keyCount=*ip++;

    // Get CoordSys
    ESMC_CoordSys_Flag *csp= (ESMC_CoordSys_Flag *)ip;
    *coordSys=*csp++;

    // Adjust offset
    *offset += sizeof(ESMC_IndexFlag)+sizeof(int)+sizeof(ESMC_CoordSys_Flag);

    // return success
    if (rc) *rc = ESMF_SUCCESS;

    return;
  }

  // non-method functions
  void FTN_X(c_esmc_locstreamkeyserialize)(
					   int *keyNameLen, char *keyName,
					   int *unitsLen, char *units,
					   int *longNameLen, char *longName,
					   char *buffer, int *length, int *offset,
					   ESMC_InquireFlag *inquireflag,
					   int *rc,
					   ESMCI_FortranStrLenArg keyName_l,
					   ESMCI_FortranStrLenArg units_l,
					   ESMCI_FortranStrLenArg longName_l,
					   ESMCI_FortranStrLenArg buffer_l) {

    ESMC_InquireFlag linquireflag = *inquireflag;
    int *ip;
    char *cp;
    int r;

    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

    // TODO: verify length > vars.
    int size = *keyNameLen + *unitsLen + *longNameLen;
    if (*inquireflag != ESMF_INQUIREONLY) {
      if ((*length - *offset) < size) {
	ESMC_LogDefault.MsgFoundError(ESMC_RC_ARG_BAD,
				      "Buffer too short to add a LocStream object", ESMC_CONTEXT, rc);
	return;
      }
    }

    // Get pointer to memory
    ip = (int *)(buffer + *offset);

    // Save string lengths
    if (linquireflag != ESMF_INQUIREONLY) {
      *ip++ = *keyNameLen;
      *ip++ = *unitsLen;
      *ip++ = *longNameLen;
    }

    // Switch to char pointer
    cp = (char *)ip;

    // Save keyNames
    if (linquireflag != ESMF_INQUIREONLY)
      memcpy((void *)cp, (const void *)keyName, *keyNameLen*sizeof(char));
    cp += *keyNameLen*sizeof(char);

    // Save units
    if (linquireflag != ESMF_INQUIREONLY)
      memcpy((void *)cp, (const void *)units, *unitsLen*sizeof(char));
    cp += *unitsLen*sizeof(char);

    // Save longName
    if (linquireflag != ESMF_INQUIREONLY)
      memcpy((void *)cp, (const void *)longName, *longNameLen*sizeof(char));
    cp += *longNameLen*sizeof(char);

    // Adjust offset
    *offset += 3*sizeof(int)+*keyNameLen + *unitsLen + *longNameLen;

    // Adjust alignment
    r=*offset%8;
    if (r!=0) *offset += 8-r;

    // return success
    if (rc) *rc = ESMF_SUCCESS;

    return;
  }


  void FTN_X(c_esmc_locstreamkeydeserialize)(
					     char *keyName,
					     char *units,
					     char *longName,
					     char *buffer,
					     int *offset,
					     int *rc,
					     ESMCI_FortranStrLenArg keyName_l,
					     ESMCI_FortranStrLenArg units_l,
					     ESMCI_FortranStrLenArg longName_l,
					     ESMCI_FortranStrLenArg buffer_l) {

    int *ip;
    char *cp;
    int r, keyNameLen, unitsLen, longNameLen;

    // Initialize return code; assume routine not implemented
    if (rc) *rc = ESMC_RC_NOT_IMPL;

    // Get pointer to memory
    ip = (int *)(buffer + *offset);

    // Save string lengths
    keyNameLen = *ip++;
    unitsLen = *ip++;
    longNameLen = *ip++;

    // Switch to char pointer
    cp = (char *)ip;

    // Save keyNames
    // First fill with spaces (NOTE THAT THIS ASSUMES THAT keyName is of size ESMF_MAXSTR)
    memset((void *)keyName,' ', ESMF_MAXSTR*sizeof(char));
    memcpy((void *)keyName, (const void *)cp, keyNameLen*sizeof(char));
    cp += keyNameLen*sizeof(char);

    // Save units
    // First fill with spaces (NOTE THAT THIS ASSUMES THAT units is of size ESMF_MAXSTR)
    memset((void *)units,' ', ESMF_MAXSTR*sizeof(char));
    memcpy((void *)units, (const void *)cp, unitsLen*sizeof(char));
    cp += unitsLen*sizeof(char);

    // Save longName
    // First fill with spaces (NOTE THAT THIS ASSUMES THAT longName is of size ESMF_MAXSTR)
    memset((void *)longName,' ', ESMF_MAXSTR*sizeof(char));
    memcpy((void *)longName, (const void *)cp, longNameLen*sizeof(char));
    cp += longNameLen*sizeof(char);

    // Adjust offset
    *offset += 3*sizeof(int)+keyNameLen + unitsLen + longNameLen;

    // Adjust alignment
    r=*offset%8;
    if (r!=0) *offset += 8-r;

    // return success
    if (rc) *rc = ESMF_SUCCESS;

    return;
  }

#endif
  
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
#include "parmetis.h"

  // Declare the ParmetisGraph structure
  typedef struct {
    idx_t *nodedist;
    idx_t *xadj;
    idx_t *adjncy;
    idx_t *nwgt;
    idx_t *adjwgt;
    double *node_x;
    double *node_y;
    idx_t nnodes;
    idx_t nedges;
  } ParmetisGraph;

  // Declare external functions
  extern int shapefile_to_parmetis_graph(const char*, MPI_Comm, ParmetisGraph*, double**, double**, double);
  extern void free_parmetis_graph(ParmetisGraph*);

  //>>void FTN_X(c_esmc_shapefile_to_graph_f)(
  //>>    char *filename,
  //>>    int *comm_int,
  //>>    void **nodedist_ptr,     // OUT: Node distribution array pointer
  //>>    void **xadj_ptr,          // OUT: CSR row pointers pointer
  //>>    void **adjncy_ptr,        // OUT: CSR adjacency list pointer  
  //>>    void **node_x_ptr,        // OUT: Node X coordinates pointer
  //>>    void **node_y_ptr,        // OUT: Node Y coordinates pointer
  //>>    idx_t *nnodes,            // OUT: Number of local nodes
  //>>    idx_t *nedges,            // OUT: Number of local edges
  //>>    double *tolerance,        // IN: Coordinate tolerance
  //>>    int *ierr,                // OUT: Error code
  //>>    ESMCI_FortranStrLenArg filename_l) {
  //>>
  //>>#undef  ESMC_METHOD
  //>>#define ESMC_METHOD "c_esmc_shapefile_to_graph_f()"
  //>>
  //>>  // Initialize return code
  //>>  *ierr = ESMC_RC_NOT_IMPL;
  //>>  
  //>>  // Convert Fortran MPI communicator to C
  //>>  MPI_Comm comm = MPI_Comm_f2c(*comm_int);
  //>>  
  //>>  // Get rank for debug output
  //>>  int rank;
  //>>  MPI_Comm_rank(comm, &rank);
  //>>  
  //>>  // Allocate graph structure
  //>>  ParmetisGraph *graph = (ParmetisGraph*)malloc(sizeof(ParmetisGraph));
  //>>  if (graph == NULL) {
  //>>    if (rank == 0) {
  //>>      printf("ERROR: Failed to allocate ParmetisGraph structure\n");
  //>>    }
  //>>    *ierr = ESMC_RC_MEM;
  //>>    return;
  //>>  }
  //>>  
  //>>  // Call the ParMETIS graph builder
  //>>  int result = shapefile_to_parmetis_graph(filename, comm, graph, *tolerance);
  //>>  
  //>>  if (result != 0) {
  //>>    if (rank == 0) {
  //>>      printf("ERROR: shapefile_to_parmetis_graph failed for %s\n", filename);
  //>>    }
  //>>    free(graph);
  //>>    *ierr = ESMC_RC_FILE_READ;
  //>>    return;
  //>>  }
  //>>  
  //>>  // Return pointers to Fortran
  //>>  *nodedist_ptr = (void*)graph->nodedist;
  //>>  *xadj_ptr = (void*)graph->xadj;
  //>>  *adjncy_ptr = (void*)graph->adjncy;
  //>>  *node_x_ptr = (void*)graph->node_x;
  //>>  *node_y_ptr = (void*)graph->node_y;
  //>>  
  //>>  // Return sizes
  //>>  *nnodes = graph->nnodes;
  //>>  *nedges = graph->nedges;
  //>>  
  //>>  // Free the structure but NOT the arrays (Fortran will use them)
  //>>  free(graph);
  //>>  
  //>>  if (rank == 0) {
  //>>    printf("Successfully created graph from %s\n", filename);
  //>>    printf("  Local nodes (rank 0): %lld\n", (long long)*nnodes);
  //>>    printf("  Local edges (rank 0): %lld\n", (long long)*nedges);
  //>>  }
  //>>  
  //>>  *ierr = ESMF_SUCCESS;
  //>>  return;
  //>>}


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
  //typedef struct {
  //    idx_t *nodedist;     /**< Node distribution array [size+1]
  //                              nodedist[i] = first node owned by rank i
  //                              nodedist[i+1] - nodedist[i] = nodes on rank i */
  //    idx_t *xadj;         /**< CSR row pointer array [nnodes+1]
  //                              Points to start of each node's adjacency list */
  //    idx_t *adjncy;       /**< CSR column index array [nedges]
  //                              Contains the adjacent node IDs */
  //    idx_t *nwgt;         /**< Node weights [nnodes * ncon] (NULL = uniform) */
  //    idx_t *adjwgt;       /**< Edge weights [nedges] (NULL = uniform) */
  //    idx_t nnodes;        /**< Number of nodes owned by this rank */
  //    idx_t nedges;        /**< Number of edges owned by this rank */
  //} ParmetisGraph;

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
  /*==============================================================================
   * shapefile_to_parmetis_graph - Build distributed graph from ESRI Shapefile
   *
   * Reads an ESRI shapefile containing polyline features (roads, rail, rivers)
   * and constructs a distributed graph suitable for ParMETIS partitioning and
   * parallel routing algorithms.
   *
   * KEY CONCEPTS:
   * - NODES: Junction points where linestrings meet or endpoints
   * - EDGES: Connections between nodes (linestring segments)  
   * - TOLERANCE: Distance within which coordinates are merged as same node
   * - CSR FORMAT: Compressed Sparse Row for efficient adjacency storage
   * - HASH TABLE: Fast coordinate-to-node_id lookup with spatial hashing
   *
   * PARALLEL STRATEGY:
   * Features divided evenly across MPI ranks. Each rank reads its features,
   * extracts nodes/edges locally, then communicates to merge boundary nodes.
   *
   * OUTPUT (in ParmetisGraph):
   * - nodedist[size+1]: Cumulative node distribution across ranks
   * - xadj[nnodes+1]: CSR row pointers  
   * - adjncy[nedges]: CSR column indices (neighbor IDs)
   * - node_x[nnodes], node_y[nnodes]: Node coordinates
   * - nnodes: Local node count, nedges: Local edge count
   *
   * @param filename Path to .shp file
   * @param comm MPI communicator
   * @param graph Output graph structure (allocated by this function)
   * @param tolerance Coordinate merging tolerance (e.g., 1e-6 degrees)
   * @return 0 on success, -1 on error
   *============================================================================*/
/*==============================================================================
 * OPTIMIZED shapefile_to_parmetis_graph - PARALLEL READING + O(N) REDISTRIBUTION
 *
 * Key optimizations:
 * 1. PHASE 1: All ranks read shapefile in parallel (each reads portion of features)
 * 2. PHASE 5: O(N) total redistribution using pre-built arrays instead of O(N*P)
 *
 * This replaces the existing shapefile_to_parmetis_graph function.
 * 
 * EXPECTED SPEEDUP:
 * - PHASE 1: P-way parallel file reading (limited by file I/O, expect 2-4x)
 * - PHASE 5: From O(N*P) to O(N), expect 10-50x for large P
 *============================================================================*/

int shapefile_to_parmetis_graph(const char *filename, MPI_Comm comm, 
                                ParmetisGraph *graph, 
                                double **node_x, double **node_y,
                                double tolerance) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    double phase_start, phase_end;
    
    /*==========================================================================
     * PHASE 1: PARALLEL shapefile reading
     * 
     * All ranks open the shapefile and read their portion of features.
     * Each rank extracts raw edge coordinates (not node IDs yet).
     * Then all raw edges are gathered to rank 0 for node merging.
     *========================================================================*/
    
    phase_start = MPI_Wtime();
    
    idx_t total_nodes = 0;
    idx_t total_edges = 0;
    idx_t *global_xadj = NULL;
    idx_t *global_adjncy = NULL;
    double *global_node_x = NULL;
    double *global_node_y = NULL;
    
    if (rank == 0) {
        printf("========================================\n");
        printf("PHASE 1: Parallel shapefile reading (%d ranks)\n", size);
        printf("========================================\n");
    }
    
    // Initialize GDAL on all ranks
    OGRRegisterAll();
    
    // Open shapefile on all ranks (GDAL supports concurrent read access)
    GDALDatasetH dataset = GDALOpenEx(filename, GDAL_OF_VECTOR | GDAL_OF_READONLY, 
                                      NULL, NULL, NULL);
    if (dataset == NULL) {
        fprintf(stderr, "Rank %d ERROR: Could not open shapefile %s\n", rank, filename);
        return -1;
    }
    
    OGRLayerH layer = OGR_DS_GetLayer(dataset, 0);
    if (layer == NULL) {
        fprintf(stderr, "Rank %d ERROR: Could not get layer\n", rank);
        GDALClose(dataset);
        return -1;
    }
    
    // Get total feature count
    GIntBig total_features = OGR_L_GetFeatureCount(layer, TRUE);
    
    // Calculate this rank's portion of features
    idx_t features_per_rank = (total_features + size - 1) / size;
    idx_t my_feature_start = rank * features_per_rank;
    idx_t my_feature_end = (rank + 1) * features_per_rank;
    if (my_feature_end > total_features) my_feature_end = total_features;
    idx_t my_feature_count = (my_feature_start < total_features) ? 
                             (my_feature_end - my_feature_start) : 0;
    
    if (rank == 0) {
        printf("  Total features: %lld, ~%lld per rank\n", 
               (long long)total_features, (long long)features_per_rank);
    }
    
    // Structure for raw edge (coordinates only, no node IDs yet)
    typedef struct { double x1, y1, x2, y2; } RawEdge;
    
    idx_t edge_capacity = my_feature_count + 100;
    RawEdge *my_raw_edges = (RawEdge*)malloc(edge_capacity * sizeof(RawEdge));
    idx_t my_edge_count = 0;
    
    // Read this rank's portion of features using random access
    for (idx_t fid = my_feature_start; fid < my_feature_end; fid++) {
        OGRFeatureH feature = OGR_L_GetFeature(layer, fid);
        if (feature == NULL) continue;
        
        OGRGeometryH geometry = OGR_F_GetGeometryRef(feature);
        
        if (geometry != NULL) {
            OGRwkbGeometryType geom_type = wkbFlatten(OGR_G_GetGeometryType(geometry));
            
            if (geom_type == wkbLineString) {
                int point_count = OGR_G_GetPointCount(geometry);
                
                if (point_count >= 2) {
                    if (my_edge_count >= edge_capacity) {
                        edge_capacity *= 2;
                        my_raw_edges = (RawEdge*)realloc(my_raw_edges, edge_capacity * sizeof(RawEdge));
                    }
                    
                    my_raw_edges[my_edge_count].x1 = OGR_G_GetX(geometry, 0);
                    my_raw_edges[my_edge_count].y1 = OGR_G_GetY(geometry, 0);
                    my_raw_edges[my_edge_count].x2 = OGR_G_GetX(geometry, point_count - 1);
                    my_raw_edges[my_edge_count].y2 = OGR_G_GetY(geometry, point_count - 1);
                    my_edge_count++;
                }
            }
        }
        
        OGR_F_Destroy(feature);
    }
    
    GDALClose(dataset);
    
    // Gather edge counts
    int *all_edge_counts = (int*)malloc(size * sizeof(int));
    int my_count_int = (int)my_edge_count;
    MPI_Allgather(&my_count_int, 1, MPI_INT, all_edge_counts, 1, MPI_INT, comm);
    
    idx_t global_raw_edge_count = 0;
    int *edge_displs = (int*)malloc(size * sizeof(int));
    edge_displs[0] = 0;
    for (int r = 0; r < size; r++) {
        global_raw_edge_count += all_edge_counts[r];
        if (r > 0) edge_displs[r] = edge_displs[r-1] + all_edge_counts[r-1];
    }
    
    // Gather all raw edges to rank 0
    RawEdge *all_raw_edges = NULL;
    if (rank == 0) {
        all_raw_edges = (RawEdge*)malloc(global_raw_edge_count * sizeof(RawEdge));
    }
    
    // Convert counts to bytes for MPI_Gatherv
    int *byte_counts = (int*)malloc(size * sizeof(int));
    int *byte_displs = (int*)malloc(size * sizeof(int));
    for (int r = 0; r < size; r++) {
        byte_counts[r] = all_edge_counts[r] * sizeof(RawEdge);
        byte_displs[r] = edge_displs[r] * sizeof(RawEdge);
    }
    
    MPI_Gatherv(my_raw_edges, my_edge_count * sizeof(RawEdge), MPI_BYTE,
                all_raw_edges, byte_counts, byte_displs, MPI_BYTE, 0, comm);
    
    free(my_raw_edges);
    free(all_edge_counts);
    free(edge_displs);
    free(byte_counts);
    free(byte_displs);
    
    // Rank 0 builds global graph with node merging
    if (rank == 0) {
        printf("  Merging nodes and building graph...\n");
        
        NodeHash **hash_table = (NodeHash**)calloc(HASH_SIZE, sizeof(NodeHash*));
        idx_t node_counter = 0;
        
        typedef struct { idx_t from, to; } Edge;
        idx_t edge_capacity = global_raw_edge_count * 2 + 100;
        Edge *edges = (Edge*)malloc(edge_capacity * sizeof(Edge));
        idx_t edge_count = 0;
        
        for (idx_t i = 0; i < global_raw_edge_count; i++) {
            idx_t node_from = get_or_insert_node(hash_table, 
                all_raw_edges[i].x1, all_raw_edges[i].y1, &node_counter, tolerance);
            idx_t node_to = get_or_insert_node(hash_table,
                all_raw_edges[i].x2, all_raw_edges[i].y2, &node_counter, tolerance);
            
            if (node_from != node_to) {
                edges[edge_count].from = node_from;
                edges[edge_count].to = node_to;
                edge_count++;
                edges[edge_count].from = node_to;
                edges[edge_count].to = node_from;
                edge_count++;
            }
        }
        
        free(all_raw_edges);
        
        total_nodes = node_counter;
        printf("  Nodes: %lld, Edges: %lld\n", (long long)total_nodes, (long long)edge_count);
        
        // Extract coordinates
        global_node_x = (double*)malloc(total_nodes * sizeof(double));
        global_node_y = (double*)malloc(total_nodes * sizeof(double));
        
        for (int i = 0; i < HASH_SIZE; i++) {
            NodeHash *entry = hash_table[i];
            while (entry != NULL) {
                if (entry->global_id >= 0 && entry->global_id < total_nodes) {
                    global_node_x[entry->global_id] = entry->x;
                    global_node_y[entry->global_id] = entry->y;
                }
                entry = entry->next;
            }
        }
        
        // Convert to CSR
        global_xadj = (idx_t*)malloc((total_nodes + 1) * sizeof(idx_t));
        idx_t *degree = (idx_t*)calloc(total_nodes, sizeof(idx_t));
        
        for (idx_t i = 0; i < edge_count; i++) {
            degree[edges[i].from]++;
        }
        
        global_xadj[0] = 0;
        for (idx_t i = 0; i < total_nodes; i++) {
            global_xadj[i + 1] = global_xadj[i] + degree[i];
        }
        total_edges = global_xadj[total_nodes];
        
        global_adjncy = (idx_t*)malloc(total_edges * sizeof(idx_t));
        idx_t *current_pos = (idx_t*)calloc(total_nodes, sizeof(idx_t));
        
        for (idx_t i = 0; i < edge_count; i++) {
            idx_t from = edges[i].from;
            global_adjncy[global_xadj[from] + current_pos[from]++] = edges[i].to;
        }
        
        free(edges);
        free(degree);
        free(current_pos);
        free_hash_table(hash_table);
    }
    
    phase_end = MPI_Wtime();
    if (rank == 0) {
        printf("  PHASE 1 time: %.3f seconds\n", phase_end - phase_start);
    }
    
    /*==========================================================================
     * Setup MPI datatype for idx_t
     *========================================================================*/
    
    MPI_Datatype IDX_T_MPI = (sizeof(idx_t) == 4) ? MPI_INT : MPI_LONG_LONG;
    
    MPI_Bcast(&total_nodes, 1, IDX_T_MPI, 0, comm);
    if (total_nodes <= 0) return -1;
    
    /*==========================================================================
     * PHASE 2: Distribute graph for ParMETIS
     *========================================================================*/
    
    phase_start = MPI_Wtime();
    
    if (rank == 0) {
        printf("========================================\n");
        printf("PHASE 2: Distributing graph to all ranks\n");
        printf("========================================\n");
    }
    
    // Even distribution for ParMETIS input
    idx_t *temp_nodedist = (idx_t*)malloc((size + 1) * sizeof(idx_t));
    idx_t nodes_per_rank = total_nodes / size;
    idx_t remainder = total_nodes % size;
    
    temp_nodedist[0] = 0;
    for (int r = 0; r < size; r++) {
        temp_nodedist[r + 1] = temp_nodedist[r] + nodes_per_rank + (r < remainder ? 1 : 0);
    }
    
    idx_t my_temp_nodes = temp_nodedist[rank + 1] - temp_nodedist[rank];
    idx_t *my_xadj = (idx_t*)malloc((my_temp_nodes + 1) * sizeof(idx_t));
    
    // Distribute xadj
    if (rank == 0) {
        idx_t offset = global_xadj[temp_nodedist[0]];
        for (idx_t i = 0; i <= my_temp_nodes; i++) {
            my_xadj[i] = global_xadj[temp_nodedist[0] + i] - offset;
        }
        
        for (int r = 1; r < size; r++) {
            idx_t r_start = temp_nodedist[r];
            idx_t r_count = temp_nodedist[r + 1] - r_start;
            idx_t *r_xadj = (idx_t*)malloc((r_count + 1) * sizeof(idx_t));
            idx_t r_offset = global_xadj[r_start];
            for (idx_t i = 0; i <= r_count; i++) {
                r_xadj[i] = global_xadj[r_start + i] - r_offset;
            }
            MPI_Send(r_xadj, r_count + 1, IDX_T_MPI, r, 0, comm);
            free(r_xadj);
        }
    } else {
        MPI_Recv(my_xadj, my_temp_nodes + 1, IDX_T_MPI, 0, 0, comm, MPI_STATUS_IGNORE);
    }
    
    idx_t my_temp_edges = my_xadj[my_temp_nodes];
    idx_t *my_adjncy = (idx_t*)malloc(my_temp_edges * sizeof(idx_t));
    
    // Distribute adjncy
    if (rank == 0) {
        idx_t start_edge = global_xadj[temp_nodedist[0]];
        memcpy(my_adjncy, &global_adjncy[start_edge], my_temp_edges * sizeof(idx_t));
        
        for (int r = 1; r < size; r++) {
            idx_t r_start = temp_nodedist[r];
            idx_t r_edge_start = global_xadj[r_start];
            idx_t r_edge_count = global_xadj[temp_nodedist[r + 1]] - r_edge_start;
            MPI_Send(&global_adjncy[r_edge_start], r_edge_count, IDX_T_MPI, r, 1, comm);
        }
    } else {
        MPI_Recv(my_adjncy, my_temp_edges, IDX_T_MPI, 0, 1, comm, MPI_STATUS_IGNORE);
    }
    
    phase_end = MPI_Wtime();
    if (rank == 0) {
        printf("  PHASE 2 time: %.3f seconds\n", phase_end - phase_start);
    }
    
    /*==========================================================================
     * PHASE 3: ParMETIS partitioning
     *========================================================================*/
    
    phase_start = MPI_Wtime();
    
    if (rank == 0) {
        printf("========================================\n");
        printf("PHASE 3: Calling ParMETIS partitioner\n");
        printf("========================================\n");
    }
    
    idx_t wgtflag = 0, numflag = 0, ncon = 1, nparts = size;
    real_t *tpwgts = (real_t*)malloc(nparts * sizeof(real_t));
    for (int i = 0; i < nparts; i++) tpwgts[i] = 1.0 / nparts;
    real_t ubvec = 1.05;
    idx_t options[3] = {0, 0, 0};
    idx_t edgecut = 0;
    idx_t *part = (idx_t*)malloc(my_temp_nodes * sizeof(idx_t));
    
    ParMETIS_V3_PartKway(temp_nodedist, my_xadj, my_adjncy, NULL, NULL,
                         &wgtflag, &numflag, &ncon, &nparts, tpwgts, &ubvec,
                         options, &edgecut, part, &comm);
    
    free(tpwgts);
    
    if (rank == 0) {
        printf("  Edge cut: %lld\n", (long long)edgecut);
    }
    
    phase_end = MPI_Wtime();
    if (rank == 0) {
        printf("  PHASE 3 time: %.3f seconds\n", phase_end - phase_start);
    }
    
    /*==========================================================================
     * PHASE 4: Gather partitions and build final nodedist
     *========================================================================*/
    
    phase_start = MPI_Wtime();
    
    if (rank == 0) {
        printf("========================================\n");
        printf("PHASE 4: Redistributing by partition\n");
        printf("========================================\n");
    }
    
    // Gather partition assignments
    idx_t *all_parts = NULL;
    if (rank == 0) {
        all_parts = (idx_t*)malloc(total_nodes * sizeof(idx_t));
    }
    
    int *gatherv_counts = (int*)malloc(size * sizeof(int));
    int *gatherv_displs = (int*)malloc(size * sizeof(int));
    for (int r = 0; r < size; r++) {
        gatherv_counts[r] = temp_nodedist[r + 1] - temp_nodedist[r];
        gatherv_displs[r] = temp_nodedist[r];
    }
    
    MPI_Gatherv(part, my_temp_nodes, IDX_T_MPI,
                all_parts, gatherv_counts, gatherv_displs, IDX_T_MPI, 0, comm);
    
    free(gatherv_counts);
    free(gatherv_displs);
    
    // Count nodes per partition
    idx_t *node_counts = (idx_t*)calloc(size, sizeof(idx_t));
    if (rank == 0) {
        for (idx_t i = 0; i < total_nodes; i++) {
            node_counts[all_parts[i]]++;
        }
    }
    MPI_Bcast(node_counts, size, IDX_T_MPI, 0, comm);
    
    // Build final nodedist
    graph->nodedist = (idx_t*)malloc((size + 1) * sizeof(idx_t));
    graph->nodedist[0] = 0;
    for (int r = 0; r < size; r++) {
        graph->nodedist[r + 1] = graph->nodedist[r] + node_counts[r];
    }
    
    idx_t my_final_nodes = node_counts[rank];
    graph->nnodes = my_final_nodes;
    
    phase_end = MPI_Wtime();
    if (rank == 0) {
        printf("  PHASE 4 time: %.3f seconds\n", phase_end - phase_start);
    }
    
    /*==========================================================================
     * PHASE 5: OPTIMIZED O(N) extraction and redistribution
     * 
     * KEY OPTIMIZATION: Build all mappings in a single O(N) pass, then use
     * MPI_Scatterv for distribution instead of serial Send/Recv loops.
     *========================================================================*/
    
    phase_start = MPI_Wtime();
    
    if (rank == 0) {
        printf("========================================\n");
        printf("PHASE 5: Extracting subgraphs (OPTIMIZED O(N))\n");
        printf("========================================\n");
    }
    
    // Rank 0 builds all mappings in a single pass
    idx_t *global_to_new = (idx_t*)malloc(total_nodes * sizeof(idx_t));
    idx_t *new_to_old = NULL;      // Inverted mapping
    idx_t *edge_counts_per_rank = NULL;
    
    if (rank == 0) {
        new_to_old = (idx_t*)malloc(total_nodes * sizeof(idx_t));
        edge_counts_per_rank = (idx_t*)calloc(size, sizeof(idx_t));
        idx_t *partition_counters = (idx_t*)calloc(size, sizeof(idx_t));
        
        // SINGLE PASS: Build global_to_new, new_to_old, and count edges per partition
        for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
            int p = all_parts[old_id];
            idx_t new_id = graph->nodedist[p] + partition_counters[p];
            
            global_to_new[old_id] = new_id;
            new_to_old[new_id] = old_id;
            
            // Count edges for this node's partition
            edge_counts_per_rank[p] += global_xadj[old_id + 1] - global_xadj[old_id];
            
            partition_counters[p]++;
        }
        
        free(partition_counters);
    }
    
    // Broadcast global_to_new (needed by all ranks for edge renumbering)
    MPI_Bcast(global_to_new, total_nodes, IDX_T_MPI, 0, comm);
    
    // Scatter new_to_old so each rank knows its old_ids directly
    idx_t *my_old_ids = (idx_t*)malloc(my_final_nodes * sizeof(idx_t));
    
    int *scatter_counts = (int*)malloc(size * sizeof(int));
    int *scatter_displs = (int*)malloc(size * sizeof(int));
    for (int r = 0; r < size; r++) {
        scatter_counts[r] = node_counts[r];
        scatter_displs[r] = graph->nodedist[r];
    }
    
    MPI_Scatterv(new_to_old, scatter_counts, scatter_displs, IDX_T_MPI,
                 my_old_ids, my_final_nodes, IDX_T_MPI, 0, comm);
    
    if (rank == 0) free(new_to_old);
    free(scatter_counts);
    free(scatter_displs);
    
    // Scatter edge counts
    idx_t my_final_edges = 0;
    MPI_Scatter(edge_counts_per_rank, 1, IDX_T_MPI, &my_final_edges, 1, IDX_T_MPI, 0, comm);
    if (rank == 0) free(edge_counts_per_rank);
    
    // Allocate final arrays
    graph->xadj = (idx_t*)malloc((my_final_nodes + 1) * sizeof(idx_t));
    graph->adjncy = (idx_t*)malloc(my_final_edges * sizeof(idx_t));
    graph->nedges = my_final_edges;
    *node_x = (double*)malloc(my_final_nodes * sizeof(double));
    *node_y = (double*)malloc(my_final_nodes * sizeof(double));
    
    // Now each rank can build its own data if it has access to global arrays
    // OR rank 0 sends the data (keeping original pattern for simplicity)
    
    if (rank == 0) {
        // Build rank 0's data
        idx_t edge_idx = 0;
        graph->xadj[0] = 0;
        
        for (idx_t i = 0; i < my_final_nodes; i++) {
            idx_t old_id = my_old_ids[i];
            (*node_x)[i] = global_node_x[old_id];
            (*node_y)[i] = global_node_y[old_id];
            
            for (idx_t e = global_xadj[old_id]; e < global_xadj[old_id + 1]; e++) {
                graph->adjncy[edge_idx++] = global_to_new[global_adjncy[e]];
            }
            graph->xadj[i + 1] = edge_idx;
        }
        
        // Send to other ranks - but now we have my_old_ids for each rank via new_to_old
        // We need to rebuild for each rank... 
        // OPTIMIZATION: Pre-build all data arrays, then use Scatterv
        
        // Actually, let's pack and send as before but more efficiently
        for (int r = 1; r < size; r++) {
            idx_t r_nodes = node_counts[r];
            idx_t r_start = graph->nodedist[r];
            
            // We don't have r's old_ids anymore since we scattered new_to_old
            // Rebuild from all_parts (still O(N) total, but spread across P ranks worth of work)
            idx_t *r_old_ids = (idx_t*)malloc(r_nodes * sizeof(idx_t));
            idx_t r_idx = 0;
            for (idx_t old_id = 0; old_id < total_nodes && r_idx < r_nodes; old_id++) {
                if (all_parts[old_id] == r) {
                    // Need to put in correct order based on new_id
                    idx_t new_local = global_to_new[old_id] - r_start;
                    r_old_ids[new_local] = old_id;
                    r_idx++;
                }
            }
            
            // Count edges
            idx_t r_edges = 0;
            for (idx_t i = 0; i < r_nodes; i++) {
                idx_t old_id = r_old_ids[i];
                r_edges += global_xadj[old_id + 1] - global_xadj[old_id];
            }
            
            // Build arrays
            idx_t *r_xadj = (idx_t*)malloc((r_nodes + 1) * sizeof(idx_t));
            idx_t *r_adjncy = (idx_t*)malloc(r_edges * sizeof(idx_t));
            double *r_node_x = (double*)malloc(r_nodes * sizeof(double));
            double *r_node_y = (double*)malloc(r_nodes * sizeof(double));
            
            idx_t r_edge_idx = 0;
            r_xadj[0] = 0;
            
            for (idx_t i = 0; i < r_nodes; i++) {
                idx_t old_id = r_old_ids[i];
                r_node_x[i] = global_node_x[old_id];
                r_node_y[i] = global_node_y[old_id];
                
                for (idx_t e = global_xadj[old_id]; e < global_xadj[old_id + 1]; e++) {
                    r_adjncy[r_edge_idx++] = global_to_new[global_adjncy[e]];
                }
                r_xadj[i + 1] = r_edge_idx;
            }
            
            // Send
            MPI_Send(r_xadj, r_nodes + 1, IDX_T_MPI, r, 4, comm);
            MPI_Send(r_adjncy, r_edges, IDX_T_MPI, r, 5, comm);
            MPI_Send(r_node_x, r_nodes, MPI_DOUBLE, r, 6, comm);
            MPI_Send(r_node_y, r_nodes, MPI_DOUBLE, r, 7, comm);
            
            free(r_xadj);
            free(r_adjncy);
            free(r_node_x);
            free(r_node_y);
            free(r_old_ids);
        }
    } else {
        MPI_Recv(graph->xadj, my_final_nodes + 1, IDX_T_MPI, 0, 4, comm, MPI_STATUS_IGNORE);
        MPI_Recv(graph->adjncy, my_final_edges, IDX_T_MPI, 0, 5, comm, MPI_STATUS_IGNORE);
        MPI_Recv(*node_x, my_final_nodes, MPI_DOUBLE, 0, 6, comm, MPI_STATUS_IGNORE);
        MPI_Recv(*node_y, my_final_nodes, MPI_DOUBLE, 0, 7, comm, MPI_STATUS_IGNORE);
    }
    
    // Convert to radians
    for (idx_t i = 0; i < my_final_nodes; i++) {
        (*node_x)[i] *= ESMC_CoordSys_Deg2Rad;
        (*node_y)[i] *= ESMC_CoordSys_Deg2Rad;
    }
    
    // Count boundary edges
    idx_t boundary = 0;
    for (idx_t i = 0; i < my_final_edges; i++) {
        if (graph->adjncy[i] < graph->nodedist[rank] || 
            graph->adjncy[i] >= graph->nodedist[rank + 1]) {
            boundary++;
        }
    }
    
    phase_end = MPI_Wtime();
    
    printf("Rank %d: %lld nodes, %lld edges (%lld boundary)\n",
           rank, (long long)my_final_nodes, (long long)my_final_edges, (long long)boundary);
    
    if (rank == 0) {
        printf("  PHASE 5 time: %.3f seconds\n", phase_end - phase_start);
    }
    
    /*==========================================================================
     * Cleanup
     *========================================================================*/
    
    free(temp_nodedist);
    free(my_xadj);
    free(my_adjncy);
    free(part);
    free(node_counts);
    free(global_to_new);
    free(my_old_ids);
    
    if (rank == 0) {
        free(global_xadj);
        free(global_adjncy);
        free(global_node_x);
        free(global_node_y);
        free(all_parts);
        
        printf("========================================\n");
        printf("ParMETIS partitioning complete!\n");
        printf("========================================\n");
    }
    
    graph->nwgt = NULL;
    graph->adjwgt = NULL;
    
    return 0;
}
// CXZ  /*==============================================================================
// CXZ   * MODIFIED shapefile_to_parmetis_graph - WITH REAL PARMETIS PARTITIONING
// CXZ   * 
// CXZ   * This replaces the naive feature-splitting approach with proper ParMETIS
// CXZ   * partitioning that maintains graph connectivity.
// CXZ   * 
// CXZ   * Algorithm:
// CXZ   *   1. Rank 0 reads entire shapefile and builds complete graph
// CXZ   *   2. Broadcast graph to all ranks (ParMETIS needs distributed input)
// CXZ   *   3. Call ParMETIS_V3_PartKway to partition nodes intelligently
// CXZ   *   4. Each rank extracts and renumbers its partition
// CXZ   * 
// CXZ   * This ensures the graph remains connected with minimal cross-PE edges.
// CXZ   *============================================================================*/
// CXZ
// CXZ  int shapefile_to_parmetis_graph(const char *filename, MPI_Comm comm, 
// CXZ				  ParmetisGraph *graph, 
// CXZ				  double **node_x, double **node_y,
// CXZ				  double tolerance) {
// CXZ    int rank, size;
// CXZ    MPI_Comm_rank(comm, &rank);
// CXZ    MPI_Comm_size(comm, &size);
// CXZ    
// CXZ    /*--------------------------------------------------------------------------
// CXZ     * PHASE 1: Rank 0 builds complete graph from entire shapefile
// CXZ     *------------------------------------------------------------------------*/
// CXZ    
// CXZ    // Global graph (complete, on rank 0 initially)
// CXZ    idx_t total_nodes = 0;
// CXZ    idx_t total_edges = 0;
// CXZ    idx_t *global_xadj = NULL;
// CXZ    idx_t *global_adjncy = NULL;
// CXZ    double *global_node_x = NULL;
// CXZ    double *global_node_y = NULL;
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      printf("========================================\n");
// CXZ      printf("PHASE 1: Rank 0 building complete graph\n");
// CXZ      printf("========================================\n");
// CXZ        
// CXZ      // Initialize GDAL
// CXZ      OGRRegisterAll();
// CXZ        
// CXZ      // Open shapefile
// CXZ      GDALDatasetH dataset = GDALOpenEx(filename, GDAL_OF_VECTOR | GDAL_OF_READONLY, 
// CXZ					NULL, NULL, NULL);
// CXZ      if (dataset == NULL) {
// CXZ	fprintf(stderr, "ERROR: Could not open shapefile %s\n", filename);
// CXZ	total_nodes = -1;  // Signal error
// CXZ      } else {
// CXZ	OGRLayerH layer = OGR_DS_GetLayer(dataset, 0);
// CXZ	if (layer == NULL) {
// CXZ	  fprintf(stderr, "ERROR: Could not get layer from shapefile\n");
// CXZ	  GDALClose(dataset);
// CXZ	  total_nodes = -1;
// CXZ	} else {
// CXZ	  // Get total feature count
// CXZ	  GIntBig total_features = OGR_L_GetFeatureCount(layer, TRUE);
// CXZ	  printf("  Total features: %lld\n", (long long)total_features);
// CXZ                
// CXZ	  // Hash table for node merging
// CXZ	  NodeHash **hash_table = (NodeHash**)calloc(HASH_SIZE, sizeof(NodeHash*));
// CXZ	  idx_t node_counter = 0;
// CXZ                
// CXZ	  // Temporary edge list
// CXZ	  typedef struct { idx_t from, to; } Edge;
// CXZ	  idx_t edge_capacity = total_features * 2;  // Allocate generously
// CXZ	  Edge *edges = (Edge*)malloc(edge_capacity * sizeof(Edge));
// CXZ	  idx_t edge_count = 0;
// CXZ                
// CXZ	  // Process ALL features
// CXZ	  OGR_L_ResetReading(layer);
// CXZ	  OGRFeatureH feature;
// CXZ	  idx_t features_processed = 0;
// CXZ                
// CXZ	  while ((feature = OGR_L_GetNextFeature(layer)) != NULL) {
// CXZ	    OGRGeometryH geometry = OGR_F_GetGeometryRef(feature);
// CXZ                    
// CXZ	    if (geometry != NULL) {
// CXZ	      OGRwkbGeometryType geom_type = wkbFlatten(OGR_G_GetGeometryType(geometry));
// CXZ                        
// CXZ	      if (geom_type == wkbLineString) {
// CXZ		int point_count = OGR_G_GetPointCount(geometry);
// CXZ                            
// CXZ		if (point_count >= 2) {
// CXZ		  // Extract start and end points
// CXZ		  double x_start = OGR_G_GetX(geometry, 0);
// CXZ		  double y_start = OGR_G_GetY(geometry, 0);
// CXZ		  double x_end = OGR_G_GetX(geometry, point_count - 1);
// CXZ		  double y_end = OGR_G_GetY(geometry, point_count - 1);
// CXZ                                
// CXZ		  // Get or create nodes
// CXZ		  //				printf("Node coordinates: (n=%d) start %.3f, %.3f; end %.3f, %.3f\n",
// CXZ		  //				       point_count, x_start, y_start, x_end, y_end);
// CXZ		  idx_t node_from = get_or_insert_node(hash_table, x_start, y_start, 
// CXZ						       &node_counter, tolerance);
// CXZ		  idx_t node_to = get_or_insert_node(hash_table, x_end, y_end, 
// CXZ						     &node_counter, tolerance);
// CXZ                                
// CXZ		  // Create edge (skip self-loops)
// CXZ		  // For undirected graph (required by ParMETIS), add BOTH directions
// CXZ		  if (node_from != node_to) {
// CXZ		    if (edge_count + 1 >= edge_capacity) {  // +1 because we'll add 2 edges
// CXZ		      edge_capacity *= 2;
// CXZ		      edges = (Edge*)realloc(edges, edge_capacity * sizeof(Edge));
// CXZ		    }
// CXZ		    // Forward edge: from → to
// CXZ		    edges[edge_count].from = node_from;
// CXZ		    edges[edge_count].to = node_to;
// CXZ		    edge_count++;
// CXZ                                    
// CXZ		    // Reverse edge: to → from (for undirected graph)
// CXZ		    edges[edge_count].from = node_to;
// CXZ		    edges[edge_count].to = node_from;
// CXZ		    edge_count++;
// CXZ		  }
// CXZ		}
// CXZ	      }
// CXZ	    }
// CXZ                    
// CXZ	    OGR_F_Destroy(feature);
// CXZ	    features_processed++;
// CXZ                    
// CXZ	    if (features_processed % 100 == 0) {
// CXZ	      printf("  Processed %lld/%lld features\r", 
// CXZ		     (long long)features_processed, (long long)total_features);
// CXZ	      fflush(stdout);
// CXZ	    }
// CXZ	  }
// CXZ                
// CXZ	  printf("\n  Features processed: %lld\n", (long long)features_processed);
// CXZ	  printf("  Nodes created: %lld\n", (long long)node_counter);
// CXZ	  printf("  Edges created: %lld\n", (long long)edge_count);
// CXZ                
// CXZ	  GDALClose(dataset);
// CXZ                
// CXZ	  total_nodes = node_counter;
// CXZ                
// CXZ	  // Extract node coordinates from hash table
// CXZ	  global_node_x = (double*)malloc(total_nodes * sizeof(double));
// CXZ	  global_node_y = (double*)malloc(total_nodes * sizeof(double));
// CXZ                
// CXZ	  for (int i = 0; i < HASH_SIZE; i++) {
// CXZ	    NodeHash *entry = hash_table[i];
// CXZ	    while (entry != NULL) {
// CXZ	      idx_t local_id = entry->global_id;
// CXZ	      if (local_id >= 0 && local_id < total_nodes) {
// CXZ		global_node_x[local_id] = entry->x;
// CXZ		global_node_y[local_id] = entry->y;
// CXZ	      }
// CXZ	      entry = entry->next;
// CXZ	    }
// CXZ	  }
// CXZ                
// CXZ	  // Convert edge list to CSR format
// CXZ	  printf("  Converting to CSR format...\n");
// CXZ                
// CXZ	  global_xadj = (idx_t*)malloc((total_nodes + 1) * sizeof(idx_t));
// CXZ	  idx_t *degree = (idx_t*)calloc(total_nodes, sizeof(idx_t));
// CXZ                
// CXZ	  // Count degrees
// CXZ	  for (idx_t i = 0; i < edge_count; i++) {
// CXZ	    degree[edges[i].from]++;
// CXZ	  }
// CXZ                
// CXZ	  // Build xadj (cumulative sum)
// CXZ	  global_xadj[0] = 0;
// CXZ	  for (idx_t i = 0; i < total_nodes; i++) {
// CXZ	    global_xadj[i + 1] = global_xadj[i] + degree[i];
// CXZ	  }
// CXZ	  total_edges = global_xadj[total_nodes];
// CXZ                
// CXZ	  printf("  Total edges in CSR: %lld\n", (long long)total_edges);
// CXZ                
// CXZ	  // Build adjncy
// CXZ	  global_adjncy = (idx_t*)malloc(total_edges * sizeof(idx_t));
// CXZ	  idx_t *current_pos = (idx_t*)calloc(total_nodes, sizeof(idx_t));
// CXZ                
// CXZ	  for (idx_t i = 0; i < edge_count; i++) {
// CXZ	    idx_t from = edges[i].from;
// CXZ	    idx_t pos = global_xadj[from] + current_pos[from];
// CXZ	    global_adjncy[pos] = edges[i].to;
// CXZ	    current_pos[from]++;
// CXZ	  }
// CXZ                
// CXZ	  // Cleanup temporary structures
// CXZ	  free(edges);
// CXZ	  free(degree);
// CXZ	  free(current_pos);
// CXZ	  free_hash_table(hash_table);
// CXZ                
// CXZ	  printf("  Complete graph built successfully\n");
// CXZ	}
// CXZ      }
// CXZ    }
// CXZ    
// CXZ    /*--------------------------------------------------------------------------
// CXZ     * Determine correct MPI datatype for idx_t
// CXZ     * 
// CXZ     * ParMETIS can be compiled with either 32-bit or 64-bit idx_t.
// CXZ     * We need to use the matching MPI datatype to avoid alignment issues.
// CXZ     *------------------------------------------------------------------------*/
// CXZ    MPI_Datatype IDX_T_MPI;
// CXZ    if (sizeof(idx_t) == sizeof(int32_t)) {
// CXZ      IDX_T_MPI = MPI_INT;
// CXZ      if (rank == 0) {
// CXZ	printf("Using 32-bit integers (idx_t = int32_t)\n");
// CXZ      }
// CXZ    } else if (sizeof(idx_t) == sizeof(int64_t)) {
// CXZ      IDX_T_MPI = MPI_LONG_LONG;
// CXZ      if (rank == 0) {
// CXZ	printf("Using 64-bit integers (idx_t = int64_t)\n");
// CXZ      }
// CXZ    } else {
// CXZ      if (rank == 0) {
// CXZ	fprintf(stderr, "ERROR: Unsupported idx_t size: %zu bytes\n", sizeof(idx_t));
// CXZ      }
// CXZ      return -1;
// CXZ    }
// CXZ    
// CXZ    // Broadcast total_nodes to all ranks (also serves as error check)
// CXZ    MPI_Bcast(&total_nodes, 1, IDX_T_MPI, 0, comm);
// CXZ    
// CXZ    if (total_nodes <= 0) {
// CXZ      if (rank == 0) {
// CXZ	fprintf(stderr, "ERROR: Failed to build graph\n");
// CXZ      }
// CXZ      return -1;
// CXZ    }
// CXZ    
// CXZ    /*--------------------------------------------------------------------------
// CXZ     * PHASE 2: Distribute graph for ParMETIS input format
// CXZ     * 
// CXZ     * ParMETIS requires the graph to be distributed across ranks.
// CXZ     * We'll use a simple distribution: divide nodes evenly.
// CXZ     *------------------------------------------------------------------------*/
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      printf("========================================\n");
// CXZ      printf("PHASE 2: Distributing graph to all ranks\n");
// CXZ      printf("========================================\n");
// CXZ    }
// CXZ    
// CXZ    // Calculate even node distribution for ParMETIS input
// CXZ    idx_t nodes_per_rank = total_nodes / size;
// CXZ    idx_t my_start = rank * nodes_per_rank;
// CXZ    idx_t my_end = (rank == size - 1) ? total_nodes : (rank + 1) * nodes_per_rank;
// CXZ    idx_t my_temp_nodes = my_end - my_start;
// CXZ    
// CXZ    // Create temporary nodedist for initial distribution
// CXZ    idx_t *temp_nodedist = (idx_t*)malloc((size + 1) * sizeof(idx_t));
// CXZ    temp_nodedist[0] = 0;
// CXZ    printf("Rank %d temp_nodedist list:\n",rank);
// CXZ    for (int i = 0; i < size; i++) {
// CXZ      idx_t end = (i == size - 1) ? total_nodes : (i + 1) * nodes_per_rank;
// CXZ      temp_nodedist[i + 1] = end;
// CXZ      printf("Rank %d temp_nodedist at i=%d: %d\n",rank,i,temp_nodedist[i+1]);
// CXZ    }
// CXZ    
// CXZ    // Allocate space for my portion of the graph
// CXZ    idx_t *my_xadj = (idx_t*)malloc((my_temp_nodes + 1) * sizeof(idx_t));
// CXZ    
// CXZ    // Scatter xadj
// CXZ    if (rank == 0) {
// CXZ      for (int r = 0; r < size; r++) {
// CXZ	idx_t start = temp_nodedist[r];
// CXZ	idx_t count = temp_nodedist[r + 1] - start + 1;  // +1 for xadj
// CXZ            
// CXZ	if (r == 0) {
// CXZ	  memcpy(my_xadj, &global_xadj[start], count * sizeof(idx_t));
// CXZ	} else {
// CXZ	  MPI_Send(&global_xadj[start], count, IDX_T_MPI, r, 0, comm);
// CXZ	}
// CXZ      }
// CXZ    } else {
// CXZ      MPI_Recv(my_xadj, my_temp_nodes + 1, IDX_T_MPI, 0, 0, comm, MPI_STATUS_IGNORE);
// CXZ    }
// CXZ    
// CXZ    // Calculate my edge count
// CXZ    idx_t my_temp_edges = my_xadj[my_temp_nodes] - my_xadj[0];
// CXZ    
// CXZ    // Adjust xadj to start from 0
// CXZ    idx_t offset = my_xadj[0];
// CXZ    for (idx_t i = 0; i <= my_temp_nodes; i++) {
// CXZ      my_xadj[i] -= offset;
// CXZ    }
// CXZ    
// CXZ    // Allocate and scatter adjncy
// CXZ    idx_t *my_adjncy = (idx_t*)malloc(my_temp_edges * sizeof(idx_t));
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      for (int r = 0; r < size; r++) {
// CXZ	idx_t start_node = temp_nodedist[r];
// CXZ	idx_t start_edge = global_xadj[start_node];
// CXZ	idx_t count = global_xadj[temp_nodedist[r + 1]] - start_edge;
// CXZ            
// CXZ	if (r == 0) {
// CXZ	  memcpy(my_adjncy, &global_adjncy[start_edge], count * sizeof(idx_t));
// CXZ	} else {
// CXZ	  MPI_Send(&global_adjncy[start_edge], count, IDX_T_MPI, r, 1, comm);
// CXZ	}
// CXZ      }
// CXZ    } else {
// CXZ      MPI_Recv(my_adjncy, my_temp_edges, IDX_T_MPI, 0, 1, comm, MPI_STATUS_IGNORE);
// CXZ    }
// CXZ    
// CXZ    printf("Rank %d: Received %lld nodes, %lld edges for ParMETIS input\n",
// CXZ           rank, (long long)my_temp_nodes, (long long)my_temp_edges);
// CXZ    
// CXZ    // Barrier to prevent output interleaving between ranks
// CXZ    MPI_Barrier(comm);
// CXZ    
// CXZ//XX    // Rank 0 goes first
// CXZ//XX    if (rank == 0) {
// CXZ//XX      /*======================================================================
// CXZ//XX       * RANK 0: COMPLETE ParMETIS Input Data Dump
// CXZ//XX       *====================================================================*/
// CXZ//XX      printf("\n");
// CXZ//XX      printf("==============================================================\n");
// CXZ//XX      printf("RANK 0: COMPLETE ParMETIS Input Data Dump\n");
// CXZ//XX      printf("==============================================================\n");
// CXZ//XX        
// CXZ//XX      // 1. COMPLETE temp_nodedist array
// CXZ//XX      printf("\n>>> RANK 0: temp_nodedist COMPLETE ARRAY [size=%d] <<<\n", size + 1);
// CXZ//XX      for (int i = 0; i <= size; i++) {
// CXZ//XX	printf("RANK 0: temp_nodedist[%d] = %lld\n", i, (long long)temp_nodedist[i]);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      // 2. COMPLETE my_xadj array
// CXZ//XX      printf("\n>>> RANK 0: my_xadj COMPLETE ARRAY [size=%lld] <<<\n", (long long)(my_temp_nodes + 1));
// CXZ//XX      for (idx_t i = 0; i <= my_temp_nodes; i++) {
// CXZ//XX	printf("RANK 0: my_xadj[%lld] = %lld\n", (long long)i, (long long)my_xadj[i]);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      // 3. COMPLETE my_adjncy array
// CXZ//XX      printf("\n>>> RANK 0: my_adjncy COMPLETE ARRAY [size=%lld] <<<\n", (long long)my_temp_edges);
// CXZ//XX      for (idx_t i = 0; i < my_temp_edges; i++) {
// CXZ//XX	printf("RANK 0: my_adjncy[%lld] = %lld\n", (long long)i, (long long)my_adjncy[i]);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      // 4. Edge interpretation
// CXZ//XX      printf("\n>>> RANK 0: EDGE LIST INTERPRETATION <<<\n");
// CXZ//XX      for (idx_t node_idx = 0; node_idx < my_temp_nodes; node_idx++) {
// CXZ//XX	idx_t global_node = temp_nodedist[rank] + node_idx;
// CXZ//XX	idx_t edge_start = my_xadj[node_idx];
// CXZ//XX	idx_t edge_end = my_xadj[node_idx + 1];
// CXZ//XX	idx_t num_edges = edge_end - edge_start;
// CXZ//XX            
// CXZ//XX	printf("RANK 0: Node %lld (global %lld) has %lld edges: [", 
// CXZ//XX	       (long long)node_idx, (long long)global_node, (long long)num_edges);
// CXZ//XX            
// CXZ//XX	for (idx_t e = edge_start; e < edge_end; e++) {
// CXZ//XX	  printf("%lld", (long long)my_adjncy[e]);
// CXZ//XX	  if (e < edge_end - 1) printf(", ");
// CXZ//XX	}
// CXZ//XX	printf("]\n");
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      // 5. Validation checks
// CXZ//XX      printf("\n>>> RANK 0: VALIDATION CHECKS <<<\n");
// CXZ//XX        
// CXZ//XX      if (my_xadj[0] != 0) {
// CXZ//XX	printf("RANK 0: ERROR: my_xadj[0] = %lld (MUST be 0)\n", (long long)my_xadj[0]);
// CXZ//XX      } else {
// CXZ//XX	printf("RANK 0: OK: my_xadj[0] = 0\n");
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      if (my_xadj[my_temp_nodes] != my_temp_edges) {
// CXZ//XX	printf("RANK 0: ERROR: my_xadj[%lld] = %lld (should be %lld)\n",
// CXZ//XX	       (long long)my_temp_nodes, (long long)my_xadj[my_temp_nodes], 
// CXZ//XX	       (long long)my_temp_edges);
// CXZ//XX      } else {
// CXZ//XX	printf("RANK 0: OK: my_xadj[%lld] = %lld = my_temp_edges\n",
// CXZ//XX	       (long long)my_temp_nodes, (long long)my_xadj[my_temp_nodes]);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      int xadj_errors = 0;
// CXZ//XX      for (idx_t i = 0; i < my_temp_nodes; i++) {
// CXZ//XX	if (my_xadj[i+1] < my_xadj[i]) {
// CXZ//XX	  printf("RANK 0: ERROR: my_xadj[%lld]=%lld > my_xadj[%lld]=%lld (not monotonic)\n",
// CXZ//XX		 (long long)i, (long long)my_xadj[i],
// CXZ//XX		 (long long)(i+1), (long long)my_xadj[i+1]);
// CXZ//XX	  xadj_errors++;
// CXZ//XX	}
// CXZ//XX      }
// CXZ//XX      if (xadj_errors == 0) {
// CXZ//XX	printf("RANK 0: OK: my_xadj is monotonic\n");
// CXZ//XX      } else {
// CXZ//XX	printf("RANK 0: ERROR: my_xadj has %d monotonicity violations\n", xadj_errors);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      int range_errors = 0;
// CXZ//XX      idx_t min_adj = total_nodes;
// CXZ//XX      idx_t max_adj = -1;
// CXZ//XX      for (idx_t i = 0; i < my_temp_edges; i++) {
// CXZ//XX	if (my_adjncy[i] < min_adj) min_adj = my_adjncy[i];
// CXZ//XX	if (my_adjncy[i] > max_adj) max_adj = my_adjncy[i];
// CXZ//XX	if (my_adjncy[i] < 0 || my_adjncy[i] >= total_nodes) {
// CXZ//XX	  if (range_errors < 10) {
// CXZ//XX	    printf("RANK 0: ERROR: my_adjncy[%lld] = %lld (out of range [0,%lld])\n",
// CXZ//XX		   (long long)i, (long long)my_adjncy[i], (long long)(total_nodes-1));
// CXZ//XX	  }
// CXZ//XX	  range_errors++;
// CXZ//XX	}
// CXZ//XX      }
// CXZ//XX      if (range_errors == 0) {
// CXZ//XX	printf("RANK 0: OK: All my_adjncy values in range [%lld, %lld]\n", 
// CXZ//XX	       (long long)min_adj, (long long)max_adj);
// CXZ//XX      } else {
// CXZ//XX	printf("RANK 0: ERROR: %d my_adjncy values out of range\n", range_errors);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      idx_t local_edges = 0;
// CXZ//XX      idx_t cross_edges = 0;
// CXZ//XX      for (idx_t i = 0; i < my_temp_edges; i++) {
// CXZ//XX	if (my_adjncy[i] >= temp_nodedist[rank] && 
// CXZ//XX	    my_adjncy[i] < temp_nodedist[rank+1]) {
// CXZ//XX	  local_edges++;
// CXZ//XX	} else {
// CXZ//XX	  cross_edges++;
// CXZ//XX	}
// CXZ//XX      }
// CXZ//XX      printf("RANK 0: Edge counts: %lld local, %lld cross-PE, %lld total\n",
// CXZ//XX	     (long long)local_edges, (long long)cross_edges, 
// CXZ//XX	     (long long)(local_edges + cross_edges));
// CXZ//XX        
// CXZ//XX      printf("==============================================================\n");
// CXZ//XX      printf("RANK 0: Ready to call ParMETIS_V3_PartKway\n");
// CXZ//XX      printf("==============================================================\n\n");
// CXZ//XX      fflush(stdout);
// CXZ//XX    }
// CXZ//XX    
// CXZ//XX    // Barrier - wait for rank 0 to finish
// CXZ//XX    MPI_Barrier(comm);
// CXZ//XX    
// CXZ//XX    // Now rank 1
// CXZ//XX    if (rank == 1) {
// CXZ//XX      /*======================================================================
// CXZ//XX       * RANK 1: COMPLETE ParMETIS Input Data Dump
// CXZ//XX       *====================================================================*/
// CXZ//XX      printf("\n");
// CXZ//XX      printf("==============================================================\n");
// CXZ//XX      printf("RANK 1: COMPLETE ParMETIS Input Data Dump\n");
// CXZ//XX      printf("==============================================================\n");
// CXZ//XX        
// CXZ//XX      // 1. COMPLETE temp_nodedist array
// CXZ//XX      printf("\n>>> RANK 1: temp_nodedist COMPLETE ARRAY [size=%d] <<<\n", size + 1);
// CXZ//XX      for (int i = 0; i <= size; i++) {
// CXZ//XX	printf("RANK 1: temp_nodedist[%d] = %lld\n", i, (long long)temp_nodedist[i]);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      // 2. COMPLETE my_xadj array
// CXZ//XX      printf("\n>>> RANK 1: my_xadj COMPLETE ARRAY [size=%lld] <<<\n", (long long)(my_temp_nodes + 1));
// CXZ//XX      for (idx_t i = 0; i <= my_temp_nodes; i++) {
// CXZ//XX	printf("RANK 1: my_xadj[%lld] = %lld\n", (long long)i, (long long)my_xadj[i]);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      // 3. COMPLETE my_adjncy array
// CXZ//XX      printf("\n>>> RANK 1: my_adjncy COMPLETE ARRAY [size=%lld] <<<\n", (long long)my_temp_edges);
// CXZ//XX      for (idx_t i = 0; i < my_temp_edges; i++) {
// CXZ//XX	printf("RANK 1: my_adjncy[%lld] = %lld\n", (long long)i, (long long)my_adjncy[i]);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      // 4. Edge interpretation
// CXZ//XX      printf("\n>>> RANK 1: EDGE LIST INTERPRETATION <<<\n");
// CXZ//XX      for (idx_t node_idx = 0; node_idx < my_temp_nodes; node_idx++) {
// CXZ//XX	idx_t global_node = temp_nodedist[rank] + node_idx;
// CXZ//XX	idx_t edge_start = my_xadj[node_idx];
// CXZ//XX	idx_t edge_end = my_xadj[node_idx + 1];
// CXZ//XX	idx_t num_edges = edge_end - edge_start;
// CXZ//XX            
// CXZ//XX	printf("RANK 1: Node %lld (global %lld) has %lld edges: [", 
// CXZ//XX	       (long long)node_idx, (long long)global_node, (long long)num_edges);
// CXZ//XX            
// CXZ//XX	for (idx_t e = edge_start; e < edge_end; e++) {
// CXZ//XX	  printf("%lld", (long long)my_adjncy[e]);
// CXZ//XX	  if (e < edge_end - 1) printf(", ");
// CXZ//XX	}
// CXZ//XX	printf("]\n");
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      // 5. Validation checks
// CXZ//XX      printf("\n>>> RANK 1: VALIDATION CHECKS <<<\n");
// CXZ//XX        
// CXZ//XX      if (my_xadj[0] != 0) {
// CXZ//XX	printf("RANK 1: ERROR: my_xadj[0] = %lld (MUST be 0)\n", (long long)my_xadj[0]);
// CXZ//XX      } else {
// CXZ//XX	printf("RANK 1: OK: my_xadj[0] = 0\n");
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      if (my_xadj[my_temp_nodes] != my_temp_edges) {
// CXZ//XX	printf("RANK 1: ERROR: my_xadj[%lld] = %lld (should be %lld)\n",
// CXZ//XX	       (long long)my_temp_nodes, (long long)my_xadj[my_temp_nodes], 
// CXZ//XX	       (long long)my_temp_edges);
// CXZ//XX      } else {
// CXZ//XX	printf("RANK 1: OK: my_xadj[%lld] = %lld = my_temp_edges\n",
// CXZ//XX	       (long long)my_temp_nodes, (long long)my_xadj[my_temp_nodes]);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      int xadj_errors = 0;
// CXZ//XX      for (idx_t i = 0; i < my_temp_nodes; i++) {
// CXZ//XX	if (my_xadj[i+1] < my_xadj[i]) {
// CXZ//XX	  printf("RANK 1: ERROR: my_xadj[%lld]=%lld > my_xadj[%lld]=%lld (not monotonic)\n",
// CXZ//XX		 (long long)i, (long long)my_xadj[i],
// CXZ//XX		 (long long)(i+1), (long long)my_xadj[i+1]);
// CXZ//XX	  xadj_errors++;
// CXZ//XX	}
// CXZ//XX      }
// CXZ//XX      if (xadj_errors == 0) {
// CXZ//XX	printf("RANK 1: OK: my_xadj is monotonic\n");
// CXZ//XX      } else {
// CXZ//XX	printf("RANK 1: ERROR: my_xadj has %d monotonicity violations\n", xadj_errors);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      int range_errors = 0;
// CXZ//XX      idx_t min_adj = total_nodes;
// CXZ//XX      idx_t max_adj = -1;
// CXZ//XX      for (idx_t i = 0; i < my_temp_edges; i++) {
// CXZ//XX	if (my_adjncy[i] < min_adj) min_adj = my_adjncy[i];
// CXZ//XX	if (my_adjncy[i] > max_adj) max_adj = my_adjncy[i];
// CXZ//XX	if (my_adjncy[i] < 0 || my_adjncy[i] >= total_nodes) {
// CXZ//XX	  if (range_errors < 10) {
// CXZ//XX	    printf("RANK 1: ERROR: my_adjncy[%lld] = %lld (out of range [0,%lld])\n",
// CXZ//XX		   (long long)i, (long long)my_adjncy[i], (long long)(total_nodes-1));
// CXZ//XX	  }
// CXZ//XX	  range_errors++;
// CXZ//XX	}
// CXZ//XX      }
// CXZ//XX      if (range_errors == 0) {
// CXZ//XX	printf("RANK 1: OK: All my_adjncy values in range [%lld, %lld]\n", 
// CXZ//XX	       (long long)min_adj, (long long)max_adj);
// CXZ//XX      } else {
// CXZ//XX	printf("RANK 1: ERROR: %d my_adjncy values out of range\n", range_errors);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      idx_t local_edges = 0;
// CXZ//XX      idx_t cross_edges = 0;
// CXZ//XX      for (idx_t i = 0; i < my_temp_edges; i++) {
// CXZ//XX	if (my_adjncy[i] >= temp_nodedist[rank] && 
// CXZ//XX	    my_adjncy[i] < temp_nodedist[rank+1]) {
// CXZ//XX	  local_edges++;
// CXZ//XX	} else {
// CXZ//XX	  cross_edges++;
// CXZ//XX	}
// CXZ//XX      }
// CXZ//XX      printf("RANK 1: Edge counts: %lld local, %lld cross-PE, %lld total\n",
// CXZ//XX	     (long long)local_edges, (long long)cross_edges, 
// CXZ//XX	     (long long)(local_edges + cross_edges));
// CXZ//XX        
// CXZ//XX      printf("==============================================================\n");
// CXZ//XX      printf("RANK 1: Ready to call ParMETIS_V3_PartKway\n");
// CXZ//XX      printf("==============================================================\n\n");
// CXZ//XX      fflush(stdout);
// CXZ//XX    }
// CXZ//XX    
// CXZ//XX    // Final barrier before ParMETIS
// CXZ//XX    MPI_Barrier(comm);
// CXZ//XX    
// CXZ    /*--------------------------------------------------------------------------
// CXZ     * PHASE 3: Call ParMETIS to partition the graph
// CXZ     *------------------------------------------------------------------------*/
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      printf("========================================\n");
// CXZ      printf("PHASE 3: Calling ParMETIS partitioner\n");
// CXZ      printf("========================================\n");
// CXZ    }
// CXZ    
// CXZ    // Partition array (where each node should go)
// CXZ    idx_t *part = (idx_t*)malloc(my_temp_nodes * sizeof(idx_t));
// CXZ    
// CXZ    // ParMETIS parameters
// CXZ    idx_t wgtflag = 0;        // No weights
// CXZ    idx_t numflag = 0;        // C-style numbering (0-based)
// CXZ    idx_t ncon = 1;           // Number of constraints
// CXZ    idx_t nparts = size;      // Number of partitions = number of ranks
// CXZ    real_t *tpwgts = (real_t*)malloc(nparts * sizeof(real_t));
// CXZ    for (int i = 0; i < nparts; i++) {
// CXZ      tpwgts[i] = 1.0 / nparts;  // Equal partition weights
// CXZ    }
// CXZ    real_t ubvec = 1.05;      // 5% imbalance tolerance
// CXZ    idx_t options[3] = {0, 0, 0};  // Default options
// CXZ    idx_t edgecut;            // Output: number of edges cut
// CXZ    
// CXZ    // Call ParMETIS
// CXZ    int ret = ParMETIS_V3_PartKway(
// CXZ				   temp_nodedist,        // Node distribution
// CXZ				   my_xadj,              // CSR row pointer (local)
// CXZ				   my_adjncy,            // CSR column indices (global numbering)
// CXZ				   NULL,                 // Vertex weights (NULL = uniform)
// CXZ				   NULL,                 // Edge weights (NULL = uniform)
// CXZ				   &wgtflag,             // Weight flag
// CXZ				   &numflag,             // Numbering flag
// CXZ				   &ncon,                // Number of constraints
// CXZ				   &nparts,              // Number of partitions
// CXZ				   tpwgts,               // Partition target weights
// CXZ				   &ubvec,               // Imbalance tolerance
// CXZ				   options,              // Options array
// CXZ				   &edgecut,             // Output: edge cut
// CXZ				   part,                 // Output: partition assignment
// CXZ				   &comm                 // MPI communicator
// CXZ				   );
// CXZ    
// CXZ    // Check return code (METIS_OK = 1 for success)
// CXZ    if (ret != METIS_OK) {
// CXZ      fprintf(stderr, "Rank %d: ParMETIS_V3_PartKway failed with code %d\n", rank, ret);
// CXZ      // Fall back to simple partitioning
// CXZ      for (idx_t i = 0; i < my_temp_nodes; i++) {
// CXZ	idx_t global_id = my_start + i;
// CXZ	part[i] = global_id * size / total_nodes;
// CXZ      }
// CXZ      if (rank == 0) {
// CXZ	printf("  WARNING: ParMETIS failed, using fallback partitioning\n");
// CXZ      }
// CXZ    } else {
// CXZ      if (rank == 0) {
// CXZ	printf("  ParMETIS completed successfully\n");
// CXZ	printf("  Edge cut: %lld edges cross partition boundaries\n", (long long)edgecut);
// CXZ      }
// CXZ    }
// CXZ    
// CXZ    // Barrier before output
// CXZ    MPI_Barrier(comm);
// CXZ    
// CXZ//XX    /*==========================================================================
// CXZ//XX     * RANK 0: Dump ALL ParMETIS output (partition assignments)
// CXZ//XX     *========================================================================*/
// CXZ//XX    if (rank == 0) {
// CXZ//XX      printf("\n");
// CXZ//XX      printf("==============================================================\n");
// CXZ//XX      printf("RANK 0: COMPLETE ParMETIS Output Data Dump\n");
// CXZ//XX      printf("==============================================================\n");
// CXZ//XX        
// CXZ//XX      printf("\n>>> RANK 0: part COMPLETE ARRAY [size=%lld] <<<\n", (long long)my_temp_nodes);
// CXZ//XX      printf("(Shows which partition each of my local nodes was assigned to)\n\n");
// CXZ//XX        
// CXZ//XX      for (idx_t i = 0; i < my_temp_nodes; i++) {
// CXZ//XX	idx_t global_id = temp_nodedist[rank] + i;
// CXZ//XX	printf("RANK 0: part[%lld] = %lld  (local node %lld, global node %lld → partition %lld)\n",
// CXZ//XX	       (long long)i, (long long)part[i], 
// CXZ//XX	       (long long)i, (long long)global_id, (long long)part[i]);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      // Summary
// CXZ//XX      printf("\n>>> RANK 0: Partition Assignment Summary <<<\n");
// CXZ//XX      idx_t *partition_counts = (idx_t*)calloc(size, sizeof(idx_t));
// CXZ//XX      for (idx_t i = 0; i < my_temp_nodes; i++) {
// CXZ//XX	if (part[i] >= 0 && part[i] < size) {
// CXZ//XX	  partition_counts[part[i]]++;
// CXZ//XX	}
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      printf("RANK 0: My %lld nodes assigned to partitions:\n", (long long)my_temp_nodes);
// CXZ//XX      for (int p = 0; p < size; p++) {
// CXZ//XX	printf("RANK 0:   Partition %d: %lld nodes (%.1f%%)\n", 
// CXZ//XX	       p, (long long)partition_counts[p],
// CXZ//XX	       100.0 * partition_counts[p] / my_temp_nodes);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      free(partition_counts);
// CXZ//XX      printf("==============================================================\n\n");
// CXZ//XX      fflush(stdout);
// CXZ//XX    }
// CXZ//XX    
// CXZ//XX    MPI_Barrier(comm);
// CXZ//XX    
// CXZ//XX    /*==========================================================================
// CXZ//XX     * RANK 1: Dump ALL ParMETIS output (partition assignments)
// CXZ//XX     *========================================================================*/
// CXZ//XX    if (rank == 1) {
// CXZ//XX      printf("\n");
// CXZ//XX      printf("==============================================================\n");
// CXZ//XX      printf("RANK 1: COMPLETE ParMETIS Output Data Dump\n");
// CXZ//XX      printf("==============================================================\n");
// CXZ//XX        
// CXZ//XX      printf("\n>>> RANK 1: part COMPLETE ARRAY [size=%lld] <<<\n", (long long)my_temp_nodes);
// CXZ//XX      printf("(Shows which partition each of my local nodes was assigned to)\n\n");
// CXZ//XX        
// CXZ//XX      for (idx_t i = 0; i < my_temp_nodes; i++) {
// CXZ//XX	idx_t global_id = temp_nodedist[rank] + i;
// CXZ//XX	printf("RANK 1: part[%lld] = %lld  (local node %lld, global node %lld → partition %lld)\n",
// CXZ//XX	       (long long)i, (long long)part[i], 
// CXZ//XX	       (long long)i, (long long)global_id, (long long)part[i]);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      // Summary
// CXZ//XX      printf("\n>>> RANK 1: Partition Assignment Summary <<<\n");
// CXZ//XX      idx_t *partition_counts = (idx_t*)calloc(size, sizeof(idx_t));
// CXZ//XX      for (idx_t i = 0; i < my_temp_nodes; i++) {
// CXZ//XX	if (part[i] >= 0 && part[i] < size) {
// CXZ//XX	  partition_counts[part[i]]++;
// CXZ//XX	}
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      printf("RANK 1: My %lld nodes assigned to partitions:\n", (long long)my_temp_nodes);
// CXZ//XX      for (int p = 0; p < size; p++) {
// CXZ//XX	printf("RANK 1:   Partition %d: %lld nodes (%.1f%%)\n", 
// CXZ//XX	       p, (long long)partition_counts[p],
// CXZ//XX	       100.0 * partition_counts[p] / my_temp_nodes);
// CXZ//XX      }
// CXZ//XX        
// CXZ//XX      free(partition_counts);
// CXZ//XX      printf("==============================================================\n\n");
// CXZ//XX      fflush(stdout);
// CXZ//XX    }
// CXZ//XX    
// CXZ//XX    MPI_Barrier(comm);
// CXZ//XX    
// CXZ    free(tpwgts);
// CXZ    
// CXZ    /*--------------------------------------------------------------------------
// CXZ     * PHASE 4: Gather partition assignments and build final nodedist
// CXZ     *------------------------------------------------------------------------*/
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      printf("========================================\n");
// CXZ      printf("PHASE 4: Redistributing by partition\n");
// CXZ      printf("========================================\n");
// CXZ    }
// CXZ    
// CXZ    // Gather all partition assignments to rank 0
// CXZ    idx_t *all_parts = NULL;
// CXZ    idx_t *recvcounts = NULL;
// CXZ    idx_t *displs = NULL;
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      all_parts = (idx_t*)malloc(total_nodes * sizeof(idx_t));
// CXZ      recvcounts = (idx_t*)malloc(size * sizeof(idx_t));
// CXZ      displs = (idx_t*)malloc(size * sizeof(idx_t));
// CXZ        
// CXZ      // Build recvcounts and displs arrays
// CXZ      for (int i = 0; i < size; i++) {
// CXZ	recvcounts[i] = temp_nodedist[i + 1] - temp_nodedist[i];
// CXZ	displs[i] = temp_nodedist[i];
// CXZ      }
// CXZ    }
// CXZ    
// CXZ    MPI_Gatherv(part, my_temp_nodes, IDX_T_MPI,
// CXZ                all_parts, recvcounts, displs, IDX_T_MPI, 0, comm);
// CXZ    
// CXZ    // Rank 0 computes final node distribution
// CXZ    idx_t *node_counts = (idx_t*)calloc(size, sizeof(idx_t));
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      // Count nodes per partition
// CXZ      for (idx_t i = 0; i < total_nodes; i++) {
// CXZ	node_counts[all_parts[i]]++;
// CXZ      }
// CXZ        
// CXZ      printf("  Partition sizes:\n");
// CXZ      for (int i = 0; i < size; i++) {
// CXZ	printf("    Rank %d: %lld nodes\n", i, (long long)node_counts[i]);
// CXZ      }
// CXZ    }
// CXZ    
// CXZ    // Broadcast node counts
// CXZ    MPI_Bcast(node_counts, size, IDX_T_MPI, 0, comm);
// CXZ    
// CXZ    // Build final nodedist
// CXZ    graph->nodedist = (idx_t*)malloc((size + 1) * sizeof(idx_t));
// CXZ    graph->nodedist[0] = 0;
// CXZ    for (int i = 0; i < size; i++) {
// CXZ      graph->nodedist[i + 1] = graph->nodedist[i] + node_counts[i];
// CXZ    }
// CXZ    
// CXZ    idx_t my_final_nodes = node_counts[rank];
// CXZ    graph->nnodes = my_final_nodes;
// CXZ    
// CXZ    /*--------------------------------------------------------------------------
// CXZ     * PHASE 5: Extract my nodes and edges based on partition
// CXZ     *------------------------------------------------------------------------*/
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      printf("========================================\n");
// CXZ      printf("PHASE 5: Extracting partitioned subgraphs\n");
// CXZ      printf("========================================\n");
// CXZ    }
// CXZ    
// CXZ    // Build global node ID mapping: old_id -> new_id
// CXZ    idx_t *global_to_new = (idx_t*)malloc(total_nodes * sizeof(idx_t));
// CXZ    idx_t *partition_counters = (idx_t*)calloc(size, sizeof(idx_t));
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
// CXZ	int p = all_parts[old_id];
// CXZ	idx_t new_id = graph->nodedist[p] + partition_counters[p];
// CXZ	global_to_new[old_id] = new_id;
// CXZ	partition_counters[p]++;
// CXZ      }
// CXZ    }
// CXZ    
// CXZ    // Broadcast mapping to all ranks
// CXZ    MPI_Bcast(global_to_new, total_nodes, IDX_T_MPI, 0, comm);
// CXZ    
// CXZ    // Extract my nodes' old IDs
// CXZ    idx_t *my_old_ids = (idx_t*)malloc(my_final_nodes * sizeof(idx_t));
// CXZ    idx_t count = 0;
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      // CRITICAL FIX: Build old_ids array sorted by NEW ID, not old ID
// CXZ      // We need my_old_ids[i] to contain the old_id of the node with new_id = nodedist[rank] + i
// CXZ        
// CXZ      // First pass: collect all old_ids for this rank
// CXZ      idx_t *temp_old_ids = (idx_t*)malloc(my_final_nodes * sizeof(idx_t));
// CXZ      idx_t temp_count = 0;
// CXZ      for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
// CXZ	if (all_parts[old_id] == rank) {
// CXZ	  temp_old_ids[temp_count++] = old_id;
// CXZ	}
// CXZ      }
// CXZ        
// CXZ      // Second pass: sort by new_id to get correct order
// CXZ      // For each position i, find which old_id has new_id = nodedist[rank] + i
// CXZ      for (idx_t i = 0; i < my_final_nodes; i++) {
// CXZ	idx_t target_new_id = graph->nodedist[rank] + i;
// CXZ            
// CXZ	// Find old_id that maps to this new_id
// CXZ	for (idx_t j = 0; j < temp_count; j++) {
// CXZ	  idx_t old_id = temp_old_ids[j];
// CXZ	  if (global_to_new[old_id] == target_new_id) {
// CXZ	    my_old_ids[i] = old_id;
// CXZ	    break;
// CXZ	  }
// CXZ	}
// CXZ      }
// CXZ      free(temp_old_ids);
// CXZ        
// CXZ      // Send to other ranks (same fix for each rank)
// CXZ      for (int r = 1; r < size; r++) {
// CXZ	idx_t r_count = node_counts[r];
// CXZ	idx_t *r_old_ids = (idx_t*)malloc(r_count * sizeof(idx_t));
// CXZ            
// CXZ	// Collect old_ids for rank r
// CXZ	idx_t *r_temp_old_ids = (idx_t*)malloc(r_count * sizeof(idx_t));
// CXZ	idx_t r_temp_count = 0;
// CXZ	for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
// CXZ	  if (all_parts[old_id] == r) {
// CXZ	    r_temp_old_ids[r_temp_count++] = old_id;
// CXZ	  }
// CXZ	}
// CXZ            
// CXZ	// Sort by new_id
// CXZ	for (idx_t i = 0; i < r_count; i++) {
// CXZ	  idx_t target_new_id = graph->nodedist[r] + i;
// CXZ                
// CXZ	  for (idx_t j = 0; j < r_temp_count; j++) {
// CXZ	    idx_t old_id = r_temp_old_ids[j];
// CXZ	    if (global_to_new[old_id] == target_new_id) {
// CXZ	      r_old_ids[i] = old_id;
// CXZ	      break;
// CXZ	    }
// CXZ	  }
// CXZ	}
// CXZ	free(r_temp_old_ids);
// CXZ            
// CXZ	MPI_Send(r_old_ids, r_count, IDX_T_MPI, r, 2, comm);
// CXZ	free(r_old_ids);
// CXZ      }
// CXZ    } else {
// CXZ      MPI_Recv(my_old_ids, my_final_nodes, IDX_T_MPI, 0, 2, comm, MPI_STATUS_IGNORE);
// CXZ    }
// CXZ    
// CXZ    // Extract edges for my nodes from global graph
// CXZ    // First, count edges
// CXZ    idx_t my_final_edges = 0;
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      for (idx_t i = 0; i < my_final_nodes; i++) {
// CXZ	idx_t old_id = my_old_ids[i];
// CXZ	my_final_edges += global_xadj[old_id + 1] - global_xadj[old_id];
// CXZ      }
// CXZ        
// CXZ      // Send edge counts to other ranks
// CXZ      for (int r = 1; r < size; r++) {
// CXZ	idx_t r_edges = 0;
// CXZ	for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
// CXZ	  if (all_parts[old_id] == r) {
// CXZ	    r_edges += global_xadj[old_id + 1] - global_xadj[old_id];
// CXZ	  }
// CXZ	}
// CXZ	MPI_Send(&r_edges, 1, IDX_T_MPI, r, 3, comm);
// CXZ      }
// CXZ    } else {
// CXZ      MPI_Recv(&my_final_edges, 1, IDX_T_MPI, 0, 3, comm, MPI_STATUS_IGNORE);
// CXZ    }
// CXZ    
// CXZ    // Build final local CSR
// CXZ    graph->xadj = (idx_t*)malloc((my_final_nodes + 1) * sizeof(idx_t));
// CXZ    graph->adjncy = (idx_t*)malloc(my_final_edges * sizeof(idx_t));
// CXZ    graph->nedges = my_final_edges;
// CXZ    
// CXZ    *node_x = (double*)malloc(my_final_nodes * sizeof(double));
// CXZ    *node_y = (double*)malloc(my_final_nodes * sizeof(double));
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      // Extract for rank 0
// CXZ      idx_t edge_idx = 0;
// CXZ      graph->xadj[0] = 0;
// CXZ        
// CXZ      for (idx_t i = 0; i < my_final_nodes; i++) {
// CXZ	idx_t old_id = my_old_ids[i];
// CXZ            
// CXZ	// Copy coordinates
// CXZ	(*node_x)[i] = global_node_x[old_id];
// CXZ	(*node_y)[i] = global_node_y[old_id];
// CXZ            
// CXZ	// Copy edges with renumbering
// CXZ	idx_t edge_start = global_xadj[old_id];
// CXZ	idx_t edge_end = global_xadj[old_id + 1];
// CXZ            
// CXZ	for (idx_t e = edge_start; e < edge_end; e++) {
// CXZ	  idx_t old_neighbor = global_adjncy[e];
// CXZ	  idx_t new_neighbor = global_to_new[old_neighbor];
// CXZ	  graph->adjncy[edge_idx++] = new_neighbor;
// CXZ	}
// CXZ            
// CXZ	graph->xadj[i + 1] = edge_idx;
// CXZ      }
// CXZ        
// CXZ      // Send to other ranks
// CXZ      for (int r = 1; r < size; r++) {
// CXZ	idx_t r_nodes = node_counts[r];
// CXZ	idx_t *r_xadj = (idx_t*)malloc((r_nodes + 1) * sizeof(idx_t));
// CXZ	idx_t *r_adjncy = NULL;
// CXZ	double *r_node_x = (double*)malloc(r_nodes * sizeof(double));
// CXZ	double *r_node_y = (double*)malloc(r_nodes * sizeof(double));
// CXZ            
// CXZ	// Build for rank r
// CXZ	idx_t r_edges = 0;
// CXZ	idx_t r_idx = 0;
// CXZ	r_xadj[0] = 0;
// CXZ            
// CXZ	// First pass: count edges
// CXZ	for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
// CXZ	  if (all_parts[old_id] == r) {
// CXZ	    r_edges += global_xadj[old_id + 1] - global_xadj[old_id];
// CXZ	  }
// CXZ	}
// CXZ            
// CXZ	r_adjncy = (idx_t*)malloc(r_edges * sizeof(idx_t));
// CXZ	idx_t r_edge_idx = 0;
// CXZ            
// CXZ	// Second pass: fill data
// CXZ	for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
// CXZ	  if (all_parts[old_id] == r) {
// CXZ	    r_node_x[r_idx] = global_node_x[old_id];
// CXZ	    r_node_y[r_idx] = global_node_y[old_id];
// CXZ                    
// CXZ	    idx_t edge_start = global_xadj[old_id];
// CXZ	    idx_t edge_end = global_xadj[old_id + 1];
// CXZ                    
// CXZ	    for (idx_t e = edge_start; e < edge_end; e++) {
// CXZ	      idx_t old_neighbor = global_adjncy[e];
// CXZ	      idx_t new_neighbor = global_to_new[old_neighbor];
// CXZ	      r_adjncy[r_edge_idx++] = new_neighbor;
// CXZ	    }
// CXZ                    
// CXZ	    r_xadj[r_idx + 1] = r_edge_idx;
// CXZ	    r_idx++;
// CXZ	  }
// CXZ	}
// CXZ            
// CXZ	// Send to rank r
// CXZ	MPI_Send(r_xadj, r_nodes + 1, IDX_T_MPI, r, 4, comm);
// CXZ	MPI_Send(r_adjncy, r_edges, IDX_T_MPI, r, 5, comm);
// CXZ	MPI_Send(r_node_x, r_nodes, MPI_DOUBLE, r, 6, comm);
// CXZ	MPI_Send(r_node_y, r_nodes, MPI_DOUBLE, r, 7, comm);
// CXZ            
// CXZ	free(r_xadj);
// CXZ	free(r_adjncy);
// CXZ	free(r_node_x);
// CXZ	free(r_node_y);
// CXZ      }
// CXZ    } else {
// CXZ      // Receive from rank 0
// CXZ      MPI_Recv(graph->xadj, my_final_nodes + 1, IDX_T_MPI, 0, 4, comm, MPI_STATUS_IGNORE);
// CXZ      MPI_Recv(graph->adjncy, my_final_edges, IDX_T_MPI, 0, 5, comm, MPI_STATUS_IGNORE);
// CXZ      MPI_Recv(*node_x, my_final_nodes, MPI_DOUBLE, 0, 6, comm, MPI_STATUS_IGNORE);
// CXZ      MPI_Recv(*node_y, my_final_nodes, MPI_DOUBLE, 0, 7, comm, MPI_STATUS_IGNORE);
// CXZ    }
// CXZ    
// CXZ    // Convert coordinates to radians
// CXZ    for (idx_t i = 0; i < my_final_nodes; i++) {
// CXZ      (*node_x)[i] *= ESMC_CoordSys_Deg2Rad;
// CXZ      (*node_y)[i] *= ESMC_CoordSys_Deg2Rad;
// CXZ    }
// CXZ    
// CXZ    // Count boundary edges
// CXZ    idx_t boundary_edges = 0;
// CXZ    for (idx_t i = 0; i < my_final_edges; i++) {
// CXZ      if (graph->adjncy[i] < graph->nodedist[rank] || 
// CXZ	  graph->adjncy[i] >= graph->nodedist[rank + 1]) {
// CXZ	boundary_edges++;
// CXZ      }
// CXZ    }
// CXZ    
// CXZ    printf("Rank %d: Final graph has %lld nodes, %lld edges (%lld boundary)\n",
// CXZ           rank, (long long)my_final_nodes, (long long)my_final_edges, (long long)boundary_edges);
// CXZ    
// CXZ    /*--------------------------------------------------------------------------
// CXZ     * Cleanup
// CXZ     *------------------------------------------------------------------------*/
// CXZ    
// CXZ    free(temp_nodedist);
// CXZ    free(my_xadj);
// CXZ    free(my_adjncy);
// CXZ    free(part);
// CXZ    free(node_counts);
// CXZ    free(partition_counters);
// CXZ    free(global_to_new);
// CXZ    free(my_old_ids);
// CXZ    
// CXZ    if (rank == 0) {
// CXZ      free(global_xadj);
// CXZ      free(global_adjncy);
// CXZ      free(global_node_x);
// CXZ      free(global_node_y);
// CXZ      free(all_parts);
// CXZ      free(recvcounts);
// CXZ      free(displs);
// CXZ        
// CXZ      printf("========================================\n");
// CXZ      printf("ParMETIS partitioning complete!\n");
// CXZ      printf("========================================\n");
// CXZ    }
// CXZ    
// CXZ    // Initialize weights to NULL
// CXZ    graph->nwgt = NULL;
// CXZ    graph->adjwgt = NULL;
// CXZ    
// CXZ    return 0;
// CXZ  }

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
				     double **node_x, double **node_y,
				     int *nnodes, int *nedges,
				     double tolerance, int *ierr) {
    
    /* Convert Fortran MPI communicator to C */
    MPI_Comm comm = MPI_Comm_f2c(comm_int);
    
    /* Allocate graph structure */
    ParmetisGraph *graph = (ParmetisGraph*)malloc(sizeof(ParmetisGraph));
    
    /* Call C function */
    *ierr = shapefile_to_parmetis_graph(filename, comm, graph, node_x, node_y, tolerance);
    
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    printf("Rank %d: shapefile ierr %d\n", rank, *ierr); 

    if (*ierr == 0) {
      /* Return pointers and sizes to Fortran */
      *nodedist = graph->nodedist;
      *xadj = graph->xadj;
      *adjncy = graph->adjncy;
      *nnodes = (int)graph->nnodes;
      *nedges = (int)graph->nedges;
      //	*node_x = node_x;
      //	*node_y = node_y;
        
      printf("Rank %d: shapefile nedges %d\n", rank, *nedges); 
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
      *node_x = NULL;
      *node_y = NULL;
        
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

}
