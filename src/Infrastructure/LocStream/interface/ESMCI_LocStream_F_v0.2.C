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

  printf("--- nFeatures: %d\n", *nfeatures);

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

  printf("NOTE: ASSUMING DEG. CONVERTING TO RADIANS!!!\n");
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
 * MODIFIED shapefile_to_parmetis_graph - WITH REAL PARMETIS PARTITIONING
 * 
 * This replaces the naive feature-splitting approach with proper ParMETIS
 * partitioning that maintains graph connectivity.
 * 
 * Algorithm:
 *   1. Rank 0 reads entire shapefile and builds complete graph
 *   2. Broadcast graph to all ranks (ParMETIS needs distributed input)
 *   3. Call ParMETIS_V3_PartKway to partition nodes intelligently
 *   4. Each rank extracts and renumbers its partition
 * 
 * This ensures the graph remains connected with minimal cross-PE edges.
 *============================================================================*/

int shapefile_to_parmetis_graph(const char *filename, MPI_Comm comm, 
                                 ParmetisGraph *graph, 
                                 double **node_x, double **node_y,
                                 double tolerance) {
    int rank, size;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    
    /*--------------------------------------------------------------------------
     * PHASE 1: Rank 0 builds complete graph from entire shapefile
     *------------------------------------------------------------------------*/
    
    // Global graph (complete, on rank 0 initially)
    idx_t total_nodes = 0;
    idx_t total_edges = 0;
    idx_t *global_xadj = NULL;
    idx_t *global_adjncy = NULL;
    double *global_node_x = NULL;
    double *global_node_y = NULL;
    
    if (rank == 0) {
        printf("========================================\n");
        printf("PHASE 1: Rank 0 building complete graph\n");
        printf("========================================\n");
        
        // Initialize GDAL
        OGRRegisterAll();
        
        // Open shapefile
        GDALDatasetH dataset = GDALOpenEx(filename, GDAL_OF_VECTOR | GDAL_OF_READONLY, 
                                          NULL, NULL, NULL);
        if (dataset == NULL) {
            fprintf(stderr, "ERROR: Could not open shapefile %s\n", filename);
            total_nodes = -1;  // Signal error
        } else {
            OGRLayerH layer = OGR_DS_GetLayer(dataset, 0);
            if (layer == NULL) {
                fprintf(stderr, "ERROR: Could not get layer from shapefile\n");
                GDALClose(dataset);
                total_nodes = -1;
            } else {
                // Get total feature count
                GIntBig total_features = OGR_L_GetFeatureCount(layer, TRUE);
                printf("  Total features: %lld\n", (long long)total_features);
                
                // Hash table for node merging
                NodeHash **hash_table = (NodeHash**)calloc(HASH_SIZE, sizeof(NodeHash*));
                idx_t node_counter = 0;
                
                // Temporary edge list
                typedef struct { idx_t from, to; } Edge;
                idx_t edge_capacity = total_features * 2;  // Allocate generously
                Edge *edges = (Edge*)malloc(edge_capacity * sizeof(Edge));
                idx_t edge_count = 0;
                
                // Process ALL features
                OGR_L_ResetReading(layer);
                OGRFeatureH feature;
                idx_t features_processed = 0;
                
                while ((feature = OGR_L_GetNextFeature(layer)) != NULL) {
                    OGRGeometryH geometry = OGR_F_GetGeometryRef(feature);
                    
                    if (geometry != NULL) {
                        OGRwkbGeometryType geom_type = wkbFlatten(OGR_G_GetGeometryType(geometry));
                        
                        if (geom_type == wkbLineString) {
                            int point_count = OGR_G_GetPointCount(geometry);
                            
                            if (point_count >= 2) {
                                // Extract start and end points
                                double x_start = OGR_G_GetX(geometry, 0);
                                double y_start = OGR_G_GetY(geometry, 0);
                                double x_end = OGR_G_GetX(geometry, point_count - 1);
                                double y_end = OGR_G_GetY(geometry, point_count - 1);
                                
                                // Get or create nodes
                                idx_t node_from = get_or_insert_node(hash_table, x_start, y_start, 
                                                                     &node_counter, tolerance);
                                idx_t node_to = get_or_insert_node(hash_table, x_end, y_end, 
                                                                   &node_counter, tolerance);
                                
                                // Create edge (skip self-loops)
                                if (node_from != node_to) {
                                    if (edge_count >= edge_capacity) {
                                        edge_capacity *= 2;
                                        edges = (Edge*)realloc(edges, edge_capacity * sizeof(Edge));
                                    }
                                    edges[edge_count].from = node_from;
                                    edges[edge_count].to = node_to;
                                    edge_count++;
                                }
                            }
                        }
                    }
                    
                    OGR_F_Destroy(feature);
                    features_processed++;
                    
                    if (features_processed % 100 == 0) {
                        printf("  Processed %lld/%lld features\r", 
                               (long long)features_processed, (long long)total_features);
                        fflush(stdout);
                    }
                }
                
                printf("\n  Features processed: %lld\n", (long long)features_processed);
                printf("  Nodes created: %lld\n", (long long)node_counter);
                printf("  Edges created: %lld\n", (long long)edge_count);
                
                GDALClose(dataset);
                
                total_nodes = node_counter;
                
                // Extract node coordinates from hash table
                global_node_x = (double*)malloc(total_nodes * sizeof(double));
                global_node_y = (double*)malloc(total_nodes * sizeof(double));
                
                for (int i = 0; i < HASH_SIZE; i++) {
                    NodeHash *entry = hash_table[i];
                    while (entry != NULL) {
                        idx_t local_id = entry->global_id;
                        if (local_id >= 0 && local_id < total_nodes) {
                            global_node_x[local_id] = entry->x;
                            global_node_y[local_id] = entry->y;
                        }
                        entry = entry->next;
                    }
                }
                
                // Convert edge list to CSR format
                printf("  Converting to CSR format...\n");
                
                global_xadj = (idx_t*)malloc((total_nodes + 1) * sizeof(idx_t));
                idx_t *degree = (idx_t*)calloc(total_nodes, sizeof(idx_t));
                
                // Count degrees
                for (idx_t i = 0; i < edge_count; i++) {
                    degree[edges[i].from]++;
                }
                
                // Build xadj (cumulative sum)
                global_xadj[0] = 0;
                for (idx_t i = 0; i < total_nodes; i++) {
                    global_xadj[i + 1] = global_xadj[i] + degree[i];
                }
                total_edges = global_xadj[total_nodes];
                
                printf("  Total edges in CSR: %lld\n", (long long)total_edges);
                
                // Build adjncy
                global_adjncy = (idx_t*)malloc(total_edges * sizeof(idx_t));
                idx_t *current_pos = (idx_t*)calloc(total_nodes, sizeof(idx_t));
                
                for (idx_t i = 0; i < edge_count; i++) {
                    idx_t from = edges[i].from;
                    idx_t pos = global_xadj[from] + current_pos[from];
                    global_adjncy[pos] = edges[i].to;
                    current_pos[from]++;
                }
                
                // Cleanup temporary structures
                free(edges);
                free(degree);
                free(current_pos);
                free_hash_table(hash_table);
                
                printf("  Complete graph built successfully\n");
            }
        }
    }
    
    /*--------------------------------------------------------------------------
     * Determine correct MPI datatype for idx_t
     * 
     * ParMETIS can be compiled with either 32-bit or 64-bit idx_t.
     * We need to use the matching MPI datatype to avoid alignment issues.
     *------------------------------------------------------------------------*/
    MPI_Datatype IDX_T_MPI;
    if (sizeof(idx_t) == sizeof(int32_t)) {
        IDX_T_MPI = MPI_INT;
        if (rank == 0) {
            printf("Using 32-bit integers (idx_t = int32_t)\n");
        }
    } else if (sizeof(idx_t) == sizeof(int64_t)) {
        IDX_T_MPI = MPI_LONG_LONG;
        if (rank == 0) {
            printf("Using 64-bit integers (idx_t = int64_t)\n");
        }
    } else {
        if (rank == 0) {
            fprintf(stderr, "ERROR: Unsupported idx_t size: %zu bytes\n", sizeof(idx_t));
        }
        return -1;
    }
    
    // Broadcast total_nodes to all ranks (also serves as error check)
    MPI_Bcast(&total_nodes, 1, IDX_T_MPI, 0, comm);
    
    if (total_nodes <= 0) {
        if (rank == 0) {
            fprintf(stderr, "ERROR: Failed to build graph\n");
        }
        return -1;
    }
    
    /*--------------------------------------------------------------------------
     * PHASE 2: Distribute graph for ParMETIS input format
     * 
     * ParMETIS requires the graph to be distributed across ranks.
     * We'll use a simple distribution: divide nodes evenly.
     *------------------------------------------------------------------------*/
    
    if (rank == 0) {
        printf("========================================\n");
        printf("PHASE 2: Distributing graph to all ranks\n");
        printf("========================================\n");
    }
    
    // Calculate even node distribution for ParMETIS input
    idx_t nodes_per_rank = total_nodes / size;
    idx_t my_start = rank * nodes_per_rank;
    idx_t my_end = (rank == size - 1) ? total_nodes : (rank + 1) * nodes_per_rank;
    idx_t my_temp_nodes = my_end - my_start;
    
    // Create temporary nodedist for initial distribution
    idx_t *temp_nodedist = (idx_t*)malloc((size + 1) * sizeof(idx_t));
    temp_nodedist[0] = 0;
    printf("Rank %d temp_nodedist list:\n",rank);
    for (int i = 0; i < size; i++) {
        idx_t end = (i == size - 1) ? total_nodes : (i + 1) * nodes_per_rank;
        temp_nodedist[i + 1] = end;
	printf("Rank %d temp_nodedist at i=%d: %d\n",rank,i,temp_nodedist[i+1]);
    }
    
    // Allocate space for my portion of the graph
    idx_t *my_xadj = (idx_t*)malloc((my_temp_nodes + 1) * sizeof(idx_t));
    
    // Scatter xadj
    if (rank == 0) {
        for (int r = 0; r < size; r++) {
            idx_t start = temp_nodedist[r];
            idx_t count = temp_nodedist[r + 1] - start + 1;  // +1 for xadj
            
            if (r == 0) {
                memcpy(my_xadj, &global_xadj[start], count * sizeof(idx_t));
            } else {
                MPI_Send(&global_xadj[start], count, IDX_T_MPI, r, 0, comm);
            }
        }
    } else {
        MPI_Recv(my_xadj, my_temp_nodes + 1, IDX_T_MPI, 0, 0, comm, MPI_STATUS_IGNORE);
    }
    
    // Calculate my edge count
    idx_t my_temp_edges = my_xadj[my_temp_nodes] - my_xadj[0];
    
    // Adjust xadj to start from 0
    idx_t offset = my_xadj[0];
    for (idx_t i = 0; i <= my_temp_nodes; i++) {
        my_xadj[i] -= offset;
    }
    
    // Allocate and scatter adjncy
    idx_t *my_adjncy = (idx_t*)malloc(my_temp_edges * sizeof(idx_t));
    
    if (rank == 0) {
        for (int r = 0; r < size; r++) {
            idx_t start_node = temp_nodedist[r];
            idx_t start_edge = global_xadj[start_node];
            idx_t count = global_xadj[temp_nodedist[r + 1]] - start_edge;
            
            if (r == 0) {
                memcpy(my_adjncy, &global_adjncy[start_edge], count * sizeof(idx_t));
            } else {
                MPI_Send(&global_adjncy[start_edge], count, IDX_T_MPI, r, 1, comm);
            }
        }
    } else {
        MPI_Recv(my_adjncy, my_temp_edges, IDX_T_MPI, 0, 1, comm, MPI_STATUS_IGNORE);
    }
    
    printf("Rank %d: Received %lld nodes, %lld edges for ParMETIS input\n",
           rank, (long long)my_temp_nodes, (long long)my_temp_edges);
    
    // Barrier to prevent output interleaving between ranks
    MPI_Barrier(comm);
    
    // Rank 0 goes first
    if (rank == 0) {
        /*======================================================================
         * RANK 0: COMPLETE ParMETIS Input Data Dump
         *====================================================================*/
        printf("\n");
        printf("==============================================================\n");
        printf("RANK 0: COMPLETE ParMETIS Input Data Dump\n");
        printf("==============================================================\n");
        
        // 1. COMPLETE temp_nodedist array
        printf("\n>>> RANK 0: temp_nodedist COMPLETE ARRAY [size=%d] <<<\n", size + 1);
        for (int i = 0; i <= size; i++) {
            printf("RANK 0: temp_nodedist[%d] = %lld\n", i, (long long)temp_nodedist[i]);
        }
        
        // 2. COMPLETE my_xadj array
        printf("\n>>> RANK 0: my_xadj COMPLETE ARRAY [size=%lld] <<<\n", (long long)(my_temp_nodes + 1));
        for (idx_t i = 0; i <= my_temp_nodes; i++) {
            printf("RANK 0: my_xadj[%lld] = %lld\n", (long long)i, (long long)my_xadj[i]);
        }
        
        // 3. COMPLETE my_adjncy array
        printf("\n>>> RANK 0: my_adjncy COMPLETE ARRAY [size=%lld] <<<\n", (long long)my_temp_edges);
        for (idx_t i = 0; i < my_temp_edges; i++) {
            printf("RANK 0: my_adjncy[%lld] = %lld\n", (long long)i, (long long)my_adjncy[i]);
        }
        
        // 4. Edge interpretation
        printf("\n>>> RANK 0: EDGE LIST INTERPRETATION <<<\n");
        for (idx_t node_idx = 0; node_idx < my_temp_nodes; node_idx++) {
            idx_t global_node = temp_nodedist[rank] + node_idx;
            idx_t edge_start = my_xadj[node_idx];
            idx_t edge_end = my_xadj[node_idx + 1];
            idx_t num_edges = edge_end - edge_start;
            
            printf("RANK 0: Node %lld (global %lld) has %lld edges: [", 
                   (long long)node_idx, (long long)global_node, (long long)num_edges);
            
            for (idx_t e = edge_start; e < edge_end; e++) {
                printf("%lld", (long long)my_adjncy[e]);
                if (e < edge_end - 1) printf(", ");
            }
            printf("]\n");
        }
        
        // 5. Validation checks
        printf("\n>>> RANK 0: VALIDATION CHECKS <<<\n");
        
        if (my_xadj[0] != 0) {
            printf("RANK 0: ERROR: my_xadj[0] = %lld (MUST be 0)\n", (long long)my_xadj[0]);
        } else {
            printf("RANK 0: OK: my_xadj[0] = 0\n");
        }
        
        if (my_xadj[my_temp_nodes] != my_temp_edges) {
            printf("RANK 0: ERROR: my_xadj[%lld] = %lld (should be %lld)\n",
                   (long long)my_temp_nodes, (long long)my_xadj[my_temp_nodes], 
                   (long long)my_temp_edges);
        } else {
            printf("RANK 0: OK: my_xadj[%lld] = %lld = my_temp_edges\n",
                   (long long)my_temp_nodes, (long long)my_xadj[my_temp_nodes]);
        }
        
        int xadj_errors = 0;
        for (idx_t i = 0; i < my_temp_nodes; i++) {
            if (my_xadj[i+1] < my_xadj[i]) {
                printf("RANK 0: ERROR: my_xadj[%lld]=%lld > my_xadj[%lld]=%lld (not monotonic)\n",
                       (long long)i, (long long)my_xadj[i],
                       (long long)(i+1), (long long)my_xadj[i+1]);
                xadj_errors++;
            }
        }
        if (xadj_errors == 0) {
            printf("RANK 0: OK: my_xadj is monotonic\n");
        } else {
            printf("RANK 0: ERROR: my_xadj has %d monotonicity violations\n", xadj_errors);
        }
        
        int range_errors = 0;
        idx_t min_adj = total_nodes;
        idx_t max_adj = -1;
        for (idx_t i = 0; i < my_temp_edges; i++) {
            if (my_adjncy[i] < min_adj) min_adj = my_adjncy[i];
            if (my_adjncy[i] > max_adj) max_adj = my_adjncy[i];
            if (my_adjncy[i] < 0 || my_adjncy[i] >= total_nodes) {
                if (range_errors < 10) {
                    printf("RANK 0: ERROR: my_adjncy[%lld] = %lld (out of range [0,%lld])\n",
                           (long long)i, (long long)my_adjncy[i], (long long)(total_nodes-1));
                }
                range_errors++;
            }
        }
        if (range_errors == 0) {
            printf("RANK 0: OK: All my_adjncy values in range [%lld, %lld]\n", 
                   (long long)min_adj, (long long)max_adj);
        } else {
            printf("RANK 0: ERROR: %d my_adjncy values out of range\n", range_errors);
        }
        
        idx_t local_edges = 0;
        idx_t cross_edges = 0;
        for (idx_t i = 0; i < my_temp_edges; i++) {
            if (my_adjncy[i] >= temp_nodedist[rank] && 
                my_adjncy[i] < temp_nodedist[rank+1]) {
                local_edges++;
            } else {
                cross_edges++;
            }
        }
        printf("RANK 0: Edge counts: %lld local, %lld cross-PE, %lld total\n",
               (long long)local_edges, (long long)cross_edges, 
               (long long)(local_edges + cross_edges));
        
        printf("==============================================================\n");
        printf("RANK 0: Ready to call ParMETIS_V3_PartKway\n");
        printf("==============================================================\n\n");
        fflush(stdout);
    }
    
    // Barrier - wait for rank 0 to finish
    MPI_Barrier(comm);
    
    // Now rank 1
    if (rank == 1) {
        /*======================================================================
         * RANK 1: COMPLETE ParMETIS Input Data Dump
         *====================================================================*/
        printf("\n");
        printf("==============================================================\n");
        printf("RANK 1: COMPLETE ParMETIS Input Data Dump\n");
        printf("==============================================================\n");
        
        // 1. COMPLETE temp_nodedist array
        printf("\n>>> RANK 1: temp_nodedist COMPLETE ARRAY [size=%d] <<<\n", size + 1);
        for (int i = 0; i <= size; i++) {
            printf("RANK 1: temp_nodedist[%d] = %lld\n", i, (long long)temp_nodedist[i]);
        }
        
        // 2. COMPLETE my_xadj array
        printf("\n>>> RANK 1: my_xadj COMPLETE ARRAY [size=%lld] <<<\n", (long long)(my_temp_nodes + 1));
        for (idx_t i = 0; i <= my_temp_nodes; i++) {
            printf("RANK 1: my_xadj[%lld] = %lld\n", (long long)i, (long long)my_xadj[i]);
        }
        
        // 3. COMPLETE my_adjncy array
        printf("\n>>> RANK 1: my_adjncy COMPLETE ARRAY [size=%lld] <<<\n", (long long)my_temp_edges);
        for (idx_t i = 0; i < my_temp_edges; i++) {
            printf("RANK 1: my_adjncy[%lld] = %lld\n", (long long)i, (long long)my_adjncy[i]);
        }
        
        // 4. Edge interpretation
        printf("\n>>> RANK 1: EDGE LIST INTERPRETATION <<<\n");
        for (idx_t node_idx = 0; node_idx < my_temp_nodes; node_idx++) {
            idx_t global_node = temp_nodedist[rank] + node_idx;
            idx_t edge_start = my_xadj[node_idx];
            idx_t edge_end = my_xadj[node_idx + 1];
            idx_t num_edges = edge_end - edge_start;
            
            printf("RANK 1: Node %lld (global %lld) has %lld edges: [", 
                   (long long)node_idx, (long long)global_node, (long long)num_edges);
            
            for (idx_t e = edge_start; e < edge_end; e++) {
                printf("%lld", (long long)my_adjncy[e]);
                if (e < edge_end - 1) printf(", ");
            }
            printf("]\n");
        }
        
        // 5. Validation checks
        printf("\n>>> RANK 1: VALIDATION CHECKS <<<\n");
        
        if (my_xadj[0] != 0) {
            printf("RANK 1: ERROR: my_xadj[0] = %lld (MUST be 0)\n", (long long)my_xadj[0]);
        } else {
            printf("RANK 1: OK: my_xadj[0] = 0\n");
        }
        
        if (my_xadj[my_temp_nodes] != my_temp_edges) {
            printf("RANK 1: ERROR: my_xadj[%lld] = %lld (should be %lld)\n",
                   (long long)my_temp_nodes, (long long)my_xadj[my_temp_nodes], 
                   (long long)my_temp_edges);
        } else {
            printf("RANK 1: OK: my_xadj[%lld] = %lld = my_temp_edges\n",
                   (long long)my_temp_nodes, (long long)my_xadj[my_temp_nodes]);
        }
        
        int xadj_errors = 0;
        for (idx_t i = 0; i < my_temp_nodes; i++) {
            if (my_xadj[i+1] < my_xadj[i]) {
                printf("RANK 1: ERROR: my_xadj[%lld]=%lld > my_xadj[%lld]=%lld (not monotonic)\n",
                       (long long)i, (long long)my_xadj[i],
                       (long long)(i+1), (long long)my_xadj[i+1]);
                xadj_errors++;
            }
        }
        if (xadj_errors == 0) {
            printf("RANK 1: OK: my_xadj is monotonic\n");
        } else {
            printf("RANK 1: ERROR: my_xadj has %d monotonicity violations\n", xadj_errors);
        }
        
        int range_errors = 0;
        idx_t min_adj = total_nodes;
        idx_t max_adj = -1;
        for (idx_t i = 0; i < my_temp_edges; i++) {
            if (my_adjncy[i] < min_adj) min_adj = my_adjncy[i];
            if (my_adjncy[i] > max_adj) max_adj = my_adjncy[i];
            if (my_adjncy[i] < 0 || my_adjncy[i] >= total_nodes) {
                if (range_errors < 10) {
                    printf("RANK 1: ERROR: my_adjncy[%lld] = %lld (out of range [0,%lld])\n",
                           (long long)i, (long long)my_adjncy[i], (long long)(total_nodes-1));
                }
                range_errors++;
            }
        }
        if (range_errors == 0) {
            printf("RANK 1: OK: All my_adjncy values in range [%lld, %lld]\n", 
                   (long long)min_adj, (long long)max_adj);
        } else {
            printf("RANK 1: ERROR: %d my_adjncy values out of range\n", range_errors);
        }
        
        idx_t local_edges = 0;
        idx_t cross_edges = 0;
        for (idx_t i = 0; i < my_temp_edges; i++) {
            if (my_adjncy[i] >= temp_nodedist[rank] && 
                my_adjncy[i] < temp_nodedist[rank+1]) {
                local_edges++;
            } else {
                cross_edges++;
            }
        }
        printf("RANK 1: Edge counts: %lld local, %lld cross-PE, %lld total\n",
               (long long)local_edges, (long long)cross_edges, 
               (long long)(local_edges + cross_edges));
        
        printf("==============================================================\n");
        printf("RANK 1: Ready to call ParMETIS_V3_PartKway\n");
        printf("==============================================================\n\n");
        fflush(stdout);
    }
    
    // Final barrier before ParMETIS
    MPI_Barrier(comm);
    
    /*--------------------------------------------------------------------------
     * PHASE 3: Call ParMETIS to partition the graph
     *------------------------------------------------------------------------*/
    
    if (rank == 0) {
        printf("========================================\n");
        printf("PHASE 3: Calling ParMETIS partitioner\n");
        printf("========================================\n");
    }
    
    // Partition array (where each node should go)
    idx_t *part = (idx_t*)malloc(my_temp_nodes * sizeof(idx_t));
    
    // ParMETIS parameters
    idx_t wgtflag = 0;        // No weights
    idx_t numflag = 0;        // C-style numbering (0-based)
    idx_t ncon = 1;           // Number of constraints
    idx_t nparts = size;      // Number of partitions = number of ranks
    real_t *tpwgts = (real_t*)malloc(nparts * sizeof(real_t));
    for (int i = 0; i < nparts; i++) {
        tpwgts[i] = 1.0 / nparts;  // Equal partition weights
    }
    real_t ubvec = 1.05;      // 5% imbalance tolerance
    idx_t options[3] = {0, 0, 0};  // Default options
    idx_t edgecut;            // Output: number of edges cut
    
    // Call ParMETIS
    int ret = ParMETIS_V3_PartKway(
        temp_nodedist,        // Node distribution
        my_xadj,              // CSR row pointer (local)
        my_adjncy,            // CSR column indices (global numbering)
        NULL,                 // Vertex weights (NULL = uniform)
        NULL,                 // Edge weights (NULL = uniform)
        &wgtflag,             // Weight flag
        &numflag,             // Numbering flag
        &ncon,                // Number of constraints
        &nparts,              // Number of partitions
        tpwgts,               // Partition target weights
        &ubvec,               // Imbalance tolerance
        options,              // Options array
        &edgecut,             // Output: edge cut
        part,                 // Output: partition assignment
        &comm                 // MPI communicator
    );
    
    // Check return code (METIS_OK = 1 for success)
    if (ret != METIS_OK) {
        fprintf(stderr, "Rank %d: ParMETIS_V3_PartKway failed with code %d\n", rank, ret);
        // Fall back to simple partitioning
        for (idx_t i = 0; i < my_temp_nodes; i++) {
            idx_t global_id = my_start + i;
            part[i] = global_id * size / total_nodes;
        }
        if (rank == 0) {
            printf("  WARNING: ParMETIS failed, using fallback partitioning\n");
        }
    } else {
        if (rank == 0) {
            printf("  ParMETIS completed successfully\n");
            printf("  Edge cut: %lld edges cross partition boundaries\n", (long long)edgecut);
        }
    }
    
    // Barrier before output
    MPI_Barrier(comm);
    
    /*==========================================================================
     * RANK 0: Dump ALL ParMETIS output (partition assignments)
     *========================================================================*/
    if (rank == 0) {
        printf("\n");
        printf("==============================================================\n");
        printf("RANK 0: COMPLETE ParMETIS Output Data Dump\n");
        printf("==============================================================\n");
        
        printf("\n>>> RANK 0: part COMPLETE ARRAY [size=%lld] <<<\n", (long long)my_temp_nodes);
        printf("(Shows which partition each of my local nodes was assigned to)\n\n");
        
        for (idx_t i = 0; i < my_temp_nodes; i++) {
            idx_t global_id = temp_nodedist[rank] + i;
            printf("RANK 0: part[%lld] = %lld  (local node %lld, global node %lld → partition %lld)\n",
                   (long long)i, (long long)part[i], 
                   (long long)i, (long long)global_id, (long long)part[i]);
        }
        
        // Summary
        printf("\n>>> RANK 0: Partition Assignment Summary <<<\n");
        idx_t *partition_counts = (idx_t*)calloc(size, sizeof(idx_t));
        for (idx_t i = 0; i < my_temp_nodes; i++) {
            if (part[i] >= 0 && part[i] < size) {
                partition_counts[part[i]]++;
            }
        }
        
        printf("RANK 0: My %lld nodes assigned to partitions:\n", (long long)my_temp_nodes);
        for (int p = 0; p < size; p++) {
            printf("RANK 0:   Partition %d: %lld nodes (%.1f%%)\n", 
                   p, (long long)partition_counts[p],
                   100.0 * partition_counts[p] / my_temp_nodes);
        }
        
        free(partition_counts);
        printf("==============================================================\n\n");
        fflush(stdout);
    }
    
    MPI_Barrier(comm);
    
    /*==========================================================================
     * RANK 1: Dump ALL ParMETIS output (partition assignments)
     *========================================================================*/
    if (rank == 1) {
        printf("\n");
        printf("==============================================================\n");
        printf("RANK 1: COMPLETE ParMETIS Output Data Dump\n");
        printf("==============================================================\n");
        
        printf("\n>>> RANK 1: part COMPLETE ARRAY [size=%lld] <<<\n", (long long)my_temp_nodes);
        printf("(Shows which partition each of my local nodes was assigned to)\n\n");
        
        for (idx_t i = 0; i < my_temp_nodes; i++) {
            idx_t global_id = temp_nodedist[rank] + i;
            printf("RANK 1: part[%lld] = %lld  (local node %lld, global node %lld → partition %lld)\n",
                   (long long)i, (long long)part[i], 
                   (long long)i, (long long)global_id, (long long)part[i]);
        }
        
        // Summary
        printf("\n>>> RANK 1: Partition Assignment Summary <<<\n");
        idx_t *partition_counts = (idx_t*)calloc(size, sizeof(idx_t));
        for (idx_t i = 0; i < my_temp_nodes; i++) {
            if (part[i] >= 0 && part[i] < size) {
                partition_counts[part[i]]++;
            }
        }
        
        printf("RANK 1: My %lld nodes assigned to partitions:\n", (long long)my_temp_nodes);
        for (int p = 0; p < size; p++) {
            printf("RANK 1:   Partition %d: %lld nodes (%.1f%%)\n", 
                   p, (long long)partition_counts[p],
                   100.0 * partition_counts[p] / my_temp_nodes);
        }
        
        free(partition_counts);
        printf("==============================================================\n\n");
        fflush(stdout);
    }
    
    MPI_Barrier(comm);
    
    free(tpwgts);
    
    /*--------------------------------------------------------------------------
     * PHASE 4: Gather partition assignments and build final nodedist
     *------------------------------------------------------------------------*/
    
    if (rank == 0) {
        printf("========================================\n");
        printf("PHASE 4: Redistributing by partition\n");
        printf("========================================\n");
    }
    
    // Gather all partition assignments to rank 0
    idx_t *all_parts = NULL;
    idx_t *recvcounts = NULL;
    idx_t *displs = NULL;
    
    if (rank == 0) {
        all_parts = (idx_t*)malloc(total_nodes * sizeof(idx_t));
        recvcounts = (idx_t*)malloc(size * sizeof(idx_t));
        displs = (idx_t*)malloc(size * sizeof(idx_t));
        
        // Build recvcounts and displs arrays
        for (int i = 0; i < size; i++) {
            recvcounts[i] = temp_nodedist[i + 1] - temp_nodedist[i];
            displs[i] = temp_nodedist[i];
        }
    }
    
    MPI_Gatherv(part, my_temp_nodes, IDX_T_MPI,
                all_parts, recvcounts, displs, IDX_T_MPI, 0, comm);
    
    // Rank 0 computes final node distribution
    idx_t *node_counts = (idx_t*)calloc(size, sizeof(idx_t));
    
    if (rank == 0) {
        // Count nodes per partition
        for (idx_t i = 0; i < total_nodes; i++) {
            node_counts[all_parts[i]]++;
        }
        
        printf("  Partition sizes:\n");
        for (int i = 0; i < size; i++) {
            printf("    Rank %d: %lld nodes\n", i, (long long)node_counts[i]);
        }
    }
    
    // Broadcast node counts
    MPI_Bcast(node_counts, size, IDX_T_MPI, 0, comm);
    
    // Build final nodedist
    graph->nodedist = (idx_t*)malloc((size + 1) * sizeof(idx_t));
    graph->nodedist[0] = 0;
    for (int i = 0; i < size; i++) {
        graph->nodedist[i + 1] = graph->nodedist[i] + node_counts[i];
    }
    
    idx_t my_final_nodes = node_counts[rank];
    graph->nnodes = my_final_nodes;
    
    /*--------------------------------------------------------------------------
     * PHASE 5: Extract my nodes and edges based on partition
     *------------------------------------------------------------------------*/
    
    if (rank == 0) {
        printf("========================================\n");
        printf("PHASE 5: Extracting partitioned subgraphs\n");
        printf("========================================\n");
    }
    
    // Build global node ID mapping: old_id -> new_id
    idx_t *global_to_new = (idx_t*)malloc(total_nodes * sizeof(idx_t));
    idx_t *partition_counters = (idx_t*)calloc(size, sizeof(idx_t));
    
    if (rank == 0) {
        for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
            int p = all_parts[old_id];
            idx_t new_id = graph->nodedist[p] + partition_counters[p];
            global_to_new[old_id] = new_id;
            partition_counters[p]++;
        }
    }
    
    // Broadcast mapping to all ranks
    MPI_Bcast(global_to_new, total_nodes, IDX_T_MPI, 0, comm);
    
    // Extract my nodes' old IDs
    idx_t *my_old_ids = (idx_t*)malloc(my_final_nodes * sizeof(idx_t));
    idx_t count = 0;
    
    if (rank == 0) {
        // CRITICAL FIX: Build old_ids array sorted by NEW ID, not old ID
        // We need my_old_ids[i] to contain the old_id of the node with new_id = nodedist[rank] + i
        
        // First pass: collect all old_ids for this rank
        idx_t *temp_old_ids = (idx_t*)malloc(my_final_nodes * sizeof(idx_t));
        idx_t temp_count = 0;
        for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
            if (all_parts[old_id] == rank) {
                temp_old_ids[temp_count++] = old_id;
            }
        }
        
        // Second pass: sort by new_id to get correct order
        // For each position i, find which old_id has new_id = nodedist[rank] + i
        for (idx_t i = 0; i < my_final_nodes; i++) {
            idx_t target_new_id = graph->nodedist[rank] + i;
            
            // Find old_id that maps to this new_id
            for (idx_t j = 0; j < temp_count; j++) {
                idx_t old_id = temp_old_ids[j];
                if (global_to_new[old_id] == target_new_id) {
                    my_old_ids[i] = old_id;
                    break;
                }
            }
        }
        free(temp_old_ids);
        
        // Send to other ranks (same fix for each rank)
        for (int r = 1; r < size; r++) {
            idx_t r_count = node_counts[r];
            idx_t *r_old_ids = (idx_t*)malloc(r_count * sizeof(idx_t));
            
            // Collect old_ids for rank r
            idx_t *r_temp_old_ids = (idx_t*)malloc(r_count * sizeof(idx_t));
            idx_t r_temp_count = 0;
            for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
                if (all_parts[old_id] == r) {
                    r_temp_old_ids[r_temp_count++] = old_id;
                }
            }
            
            // Sort by new_id
            for (idx_t i = 0; i < r_count; i++) {
                idx_t target_new_id = graph->nodedist[r] + i;
                
                for (idx_t j = 0; j < r_temp_count; j++) {
                    idx_t old_id = r_temp_old_ids[j];
                    if (global_to_new[old_id] == target_new_id) {
                        r_old_ids[i] = old_id;
                        break;
                    }
                }
            }
            free(r_temp_old_ids);
            
            MPI_Send(r_old_ids, r_count, IDX_T_MPI, r, 2, comm);
            free(r_old_ids);
        }
    } else {
        MPI_Recv(my_old_ids, my_final_nodes, IDX_T_MPI, 0, 2, comm, MPI_STATUS_IGNORE);
    }
    
    // Extract edges for my nodes from global graph
    // First, count edges
    idx_t my_final_edges = 0;
    
    if (rank == 0) {
        for (idx_t i = 0; i < my_final_nodes; i++) {
            idx_t old_id = my_old_ids[i];
            my_final_edges += global_xadj[old_id + 1] - global_xadj[old_id];
        }
        
        // Send edge counts to other ranks
        for (int r = 1; r < size; r++) {
            idx_t r_edges = 0;
            for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
                if (all_parts[old_id] == r) {
                    r_edges += global_xadj[old_id + 1] - global_xadj[old_id];
                }
            }
            MPI_Send(&r_edges, 1, IDX_T_MPI, r, 3, comm);
        }
    } else {
        MPI_Recv(&my_final_edges, 1, IDX_T_MPI, 0, 3, comm, MPI_STATUS_IGNORE);
    }
    
    // Build final local CSR
    graph->xadj = (idx_t*)malloc((my_final_nodes + 1) * sizeof(idx_t));
    graph->adjncy = (idx_t*)malloc(my_final_edges * sizeof(idx_t));
    graph->nedges = my_final_edges;
    
    *node_x = (double*)malloc(my_final_nodes * sizeof(double));
    *node_y = (double*)malloc(my_final_nodes * sizeof(double));
    
    if (rank == 0) {
        // Extract for rank 0
        idx_t edge_idx = 0;
        graph->xadj[0] = 0;
        
        for (idx_t i = 0; i < my_final_nodes; i++) {
            idx_t old_id = my_old_ids[i];
            
            // Copy coordinates
            (*node_x)[i] = global_node_x[old_id];
            (*node_y)[i] = global_node_y[old_id];
            
            // Copy edges with renumbering
            idx_t edge_start = global_xadj[old_id];
            idx_t edge_end = global_xadj[old_id + 1];
            
            for (idx_t e = edge_start; e < edge_end; e++) {
                idx_t old_neighbor = global_adjncy[e];
                idx_t new_neighbor = global_to_new[old_neighbor];
                graph->adjncy[edge_idx++] = new_neighbor;
            }
            
            graph->xadj[i + 1] = edge_idx;
        }
        
        // Send to other ranks
        for (int r = 1; r < size; r++) {
            idx_t r_nodes = node_counts[r];
            idx_t *r_xadj = (idx_t*)malloc((r_nodes + 1) * sizeof(idx_t));
            idx_t *r_adjncy = NULL;
            double *r_node_x = (double*)malloc(r_nodes * sizeof(double));
            double *r_node_y = (double*)malloc(r_nodes * sizeof(double));
            
            // Build for rank r
            idx_t r_edges = 0;
            idx_t r_idx = 0;
            r_xadj[0] = 0;
            
            // First pass: count edges
            for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
                if (all_parts[old_id] == r) {
                    r_edges += global_xadj[old_id + 1] - global_xadj[old_id];
                }
            }
            
            r_adjncy = (idx_t*)malloc(r_edges * sizeof(idx_t));
            idx_t r_edge_idx = 0;
            
            // Second pass: fill data
            for (idx_t old_id = 0; old_id < total_nodes; old_id++) {
                if (all_parts[old_id] == r) {
                    r_node_x[r_idx] = global_node_x[old_id];
                    r_node_y[r_idx] = global_node_y[old_id];
                    
                    idx_t edge_start = global_xadj[old_id];
                    idx_t edge_end = global_xadj[old_id + 1];
                    
                    for (idx_t e = edge_start; e < edge_end; e++) {
                        idx_t old_neighbor = global_adjncy[e];
                        idx_t new_neighbor = global_to_new[old_neighbor];
                        r_adjncy[r_edge_idx++] = new_neighbor;
                    }
                    
                    r_xadj[r_idx + 1] = r_edge_idx;
                    r_idx++;
                }
            }
            
            // Send to rank r
            MPI_Send(r_xadj, r_nodes + 1, IDX_T_MPI, r, 4, comm);
            MPI_Send(r_adjncy, r_edges, IDX_T_MPI, r, 5, comm);
            MPI_Send(r_node_x, r_nodes, MPI_DOUBLE, r, 6, comm);
            MPI_Send(r_node_y, r_nodes, MPI_DOUBLE, r, 7, comm);
            
            free(r_xadj);
            free(r_adjncy);
            free(r_node_x);
            free(r_node_y);
        }
    } else {
        // Receive from rank 0
        MPI_Recv(graph->xadj, my_final_nodes + 1, IDX_T_MPI, 0, 4, comm, MPI_STATUS_IGNORE);
        MPI_Recv(graph->adjncy, my_final_edges, IDX_T_MPI, 0, 5, comm, MPI_STATUS_IGNORE);
        MPI_Recv(*node_x, my_final_nodes, MPI_DOUBLE, 0, 6, comm, MPI_STATUS_IGNORE);
        MPI_Recv(*node_y, my_final_nodes, MPI_DOUBLE, 0, 7, comm, MPI_STATUS_IGNORE);
    }
    
    // Convert coordinates to radians
    for (idx_t i = 0; i < my_final_nodes; i++) {
        (*node_x)[i] *= ESMC_CoordSys_Deg2Rad;
        (*node_y)[i] *= ESMC_CoordSys_Deg2Rad;
    }
    
    // Count boundary edges
    idx_t boundary_edges = 0;
    for (idx_t i = 0; i < my_final_edges; i++) {
        if (graph->adjncy[i] < graph->nodedist[rank] || 
            graph->adjncy[i] >= graph->nodedist[rank + 1]) {
            boundary_edges++;
        }
    }
    
    printf("Rank %d: Final graph has %lld nodes, %lld edges (%lld boundary)\n",
           rank, (long long)my_final_nodes, (long long)my_final_edges, (long long)boundary_edges);
    
    /*--------------------------------------------------------------------------
     * Cleanup
     *------------------------------------------------------------------------*/
    
    free(temp_nodedist);
    free(my_xadj);
    free(my_adjncy);
    free(part);
    free(node_counts);
    free(partition_counters);
    free(global_to_new);
    free(my_old_ids);
    
    if (rank == 0) {
        free(global_xadj);
        free(global_adjncy);
        free(global_node_x);
        free(global_node_y);
        free(all_parts);
        free(recvcounts);
        free(displs);
        
        printf("========================================\n");
        printf("ParMETIS partitioning complete!\n");
        printf("========================================\n");
    }
    
    // Initialize weights to NULL
    graph->nwgt = NULL;
    graph->adjwgt = NULL;
    
    return 0;
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
