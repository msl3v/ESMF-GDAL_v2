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
int shapefile_to_parmetis_graph(const char *filename, MPI_Comm comm, 
                                 ParmetisGraph *graph, 
                                 double **node_x, double **node_y,
				 double tolerance) {
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
    
    if (rank == 0) {
        printf("DEBUG: Shapefile opened successfully\n");
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
    
    if (rank == 0) {
        printf("DEBUG: Layer 0 retrieved successfully\n");
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
    
    printf("DEBUG Rank %d: Processing features %lld to %lld (%lld features)\n",
           rank, (long long)start_feature, (long long)end_feature-1, (long long)local_feature_count);
    
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
    printf(" local_nodes rank%d=%lld", rank, (long long)local_nodes);
//    MPI_Gather(&local_nodes, 1, MPI_INT64_T, node_counts, 1, MPI_INT64_T, 0, comm);
    MPI_Gather(&local_nodes, sizeof(idx_t), MPI_BYTE, node_counts, sizeof(idx_t), MPI_BYTE, 0, comm);
    
    if (rank == 0) {
        printf("DEBUG Rank 0: Gathered node counts:");
        for (int i = 0; i < size; i++) {
            printf(" rank%d=%lld", i, (long long)node_counts[i]);
        }
        printf("\n");
    }
    
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
    
    printf("DEBUG Rank %d: nodedist[%d]=%lld, nodedist[%d]=%lld (owns %lld nodes)\n",
           rank, rank, (long long)graph->nodedist[rank], 
           rank+1, (long long)graph->nodedist[rank+1],
           (long long)(graph->nodedist[rank+1] - graph->nodedist[rank]));
    
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
    
    // Allocate and populate coordinate arrays from hash table
    double *local_node_x = (double*)malloc(local_nodes * sizeof(double));
    double *local_node_y = (double*)malloc(local_nodes * sizeof(double));
    
    // Extract coordinates from hash table
    for (int i = 0; i < HASH_SIZE; i++) {
      NodeHash *entry = hash_table[i];
      while (entry != NULL) {
        idx_t local_id = entry->global_id;
        if (local_id >= 0 && local_id < local_nodes) {
	  local_node_x[local_id] = entry->x * ESMC_CoordSys_Deg2Rad; // Assume deg.
	  local_node_y[local_id] = entry->y * ESMC_CoordSys_Deg2Rad; // Assume deg. need rads
        }
        entry = entry->next;
      }
    }

    *node_x = local_node_x;
    *node_y = local_node_y;

    printf("DEBUG Rank %d: local_nodes=%lld, node_counter=%lld\n",
           rank, (long long)local_nodes, (long long)node_counter);
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
