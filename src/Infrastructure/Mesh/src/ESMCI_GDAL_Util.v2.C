// $Id$
//
// Earth System Modeling Framework
// Copyright (c) 2002-2023, University Corporation for Atmospheric Research,
// Massachusetts Institute of Technology, Geophysical Fluid Dynamics
// Laboratory, University of Michigan, National Centers for Environmental
// Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
// NASA Goddard Space Flight Center.
// Licensed under the University of Illinois-NCSA License.
//
//==============================================================================

//==============================================================================
//
// This file contains the Fortran interface code to link F90 and C++.
//
//------------------------------------------------------------------------------
// INCLUDES
//------------------------------------------------------------------------------

#include <string>
#include <ostream>
#include <vector>
#include <iterator>
#include <algorithm>
#include <unordered_set>

#include "ESMCI_Macros.h"
#include "ESMCI_F90Interface.h"
#include "ESMCI_LogErr.h"
#include "ESMCI_VM.h"
#include "ESMCI_CoordSys.h"
#include "ESMCI_Array.h"
#include "ESMC_Util.h"

#include "ESMCI_TraceMacros.h"  // for profiling

#include "Mesh/include/ESMCI_Mesh.h"
#include "Mesh/include/Legacy/ESMCI_MeshRead.h"
#include "Mesh/include/Regridding/ESMCI_MeshRegrid.h" //only for the conservative flag in add_elements
#include "Mesh/include/Legacy/ESMCI_MeshVTK.h"
#include "Mesh/include/Legacy/ESMCI_ParEnv.h"
#include "Mesh/include/Legacy/ESMCI_MeshUtils.h"
#include "Mesh/include/Legacy/ESMCI_GlobalIds.h"
#include "Mesh/include/ESMCI_MeshRedist.h"
#include "Mesh/include/ESMCI_MeshDual.h"
#include "Mesh/include/ESMCI_Mesh_Glue.h"
#include "Mesh/include/ESMCI_FileIO_Util.h"
#include "Mesh/include/ESMCI_GDAL_Util.h"

// These internal functions can only be used if GDAL is available
#ifdef ESMF_GDAL

// TODO: SWITCH THIS TO SHAPELIB, WHEN WE KNOW WHAT IT'S CALLED
#include <ogr_api.h>
#include <gdal.h>
#include <ogr_srs_api.h>

//-----------------------------------------------------------------------------
 // leave the following line as-is; it will insert the cvs ident string
 // into the object file for tracking purposes.
 static const char *const version = "$Id$";
//-----------------------------------------------------------------------------
using namespace ESMCI;

int elm       = 0;
int totpoints = 0;

//int numNodes, numElems;
//double *nodeXCoords,*nodeYCoords;

// Local Routines
int processPolygon(OGRGeometryH fGeom, std::vector<double> &XCoords, std::vector<double> &YCoords, std::vector<int> 
		   &elemConn, std::vector<int> &numelemConn, std::vector<int> &elemNodeIDs, int totpoints, int *nPpoints);
int processMultiPolygon(OGRGeometryH hGeom, std::vector<double> &mXCoords, std::vector<double> &mYCoords, std::vector<int> &melemConn, 
			std::vector<int> &nelemConn, std::vector<int> &elemNodeIDs, int totpoints, int *nMPpoints);
int processLineString(OGRGeometryH fGeom, std::vector<double> &XCoords, std::vector<double> &YCoords, std::vector<int> &elemNodeIDs, int *nPpoints);
int processMultiLineString(OGRGeometryH hGeom, std::vector<double> &mXCoords, std::vector<double> &mYCoords, std::vector<int> &elemNodeIDs, int *nMPpoints);
bool valueinarray(int val, int *arr, int n);
int countMultiPolygon(OGRGeometryH hGeom, int *nMPoints);

// Get the dimension of the mesh in the SHP file
// (This dimension is both the pdim and orig_sdim of the mesh)
void ESMCI_GDAL_SHP_get_dim_from_file(OGRDataSourceH hDS, char *filename, int &dim) {
#undef ESMC_METHOD
#define ESMC_METHOD "get_dim_from_SHP_file()"

  // value fixed for now. Just assume horizontal dims.
  // Geometries will be 'flattened' if they are 3D.
  // -- MSL 5/31/2023
  dim = 2; 

  return;
}

void ESMCI_GDAL_SHP_get_feature_info(OGRDataSourceH hDS, int *nFeatures, int *&FeatureIDs) {
  OGRLayerH hLayer;
  OGRFeatureH hFeature;

  // Access the layer (associate the handle)
  // Assume that index 0 is the layer we want.
  hLayer = OGR_DS_GetLayer( hDS, 0 );

  // Get the number of elements
  *nFeatures = OGR_L_GetFeatureCount(hLayer,1);

  // -- elemIDs
  FeatureIDs   = (int *)malloc(*nFeatures*sizeof(int));

  for (int i=0;i<*nFeatures;i++) {
//    hFeature = OGR_L_GetNextFeature(hLayer);
    hFeature = OGR_L_GetFeature(hLayer,i);
    FeatureIDs[i] = OGR_F_GetFID(hFeature)+1; // IDs can't be zero in meshes
    OGR_F_Destroy( hFeature );
  }


  return;
}

void ESMCI_GDAL_process_shapefile_serial(
// inputs
		       OGRDataSourceH hDS, 
// outputs
		       double *&nodeCoords, 
		       int *&nodeIDs, 
		       int *&elemIDs, 
		       int *&elemConn,
		       int *&numElemConn, 
		       int *totNumElemConn, 
		       int *numNodes, 
		       int *numElems) {
}

void ESMCI_GDAL_process_shapefile_distributed(
// inputs
		       OGRDataSourceH hDS,
		       int *nFeatures,
		       int *&featureIDs,
		       int *&globFeatureIDs,
// outputs
		       double *&nodeCoords,
		       std::vector<int> &nodeIDs,
		       std::vector<int> &elemIDs,
		       std::vector<int> &elemConn,
		       std::vector<double> &elemCoords,
		       std::vector<int> &numElemConn,
		       int *totNumElemConn,
		       int *nNodes,
		       int *nElems) {

  std::vector<double> XCoords;
  std::vector<double> YCoords;

  int pet_rank;// = 0;
//#ifdef ESMF_COMM
  MPI_Comm_rank(MPI_COMM_WORLD, &pet_rank);
//#endif

  printf("<<>> %d nfeatures in distrib at pet %d\n",*nFeatures,pet_rank);
  for (int i = 0; i < *nFeatures; i++) {
    printf("<<>> PE %d: featureIDs[%d]=%d\n",pet_rank,i,featureIDs[i]);
  }

  OGRRegisterAll(); // register all the drivers

  OGRLayerH hLayer;
  OGRFeatureH hFeature;
  OGRGeometryH hGeom = nullptr;

  // Step 1: Access the layer (assume index 0 is the layer we want)
  hLayer = OGR_DS_GetLayer( hDS, 0 );
  if (!hLayer) return; // FIX #5 variant: guard against null layer

  // Step 2: Get number of elements.
  // Use the caller-supplied nFeatures to avoid a redundant OGR_L_GetFeatureCount
  // call (fix #9). *nElems is set here for the caller's benefit.
  *nElems = *nFeatures;

  // FIX #10: removed the getLayerInfo() call — its results (nPoints, nGeom)
  // were computed but never used, wasting a full O(F) sequential pass.

  // Step 3: Initialise counters and rewind layer once.
  totpoints = 0;
  int localpoints = 0;
  *totNumElemConn = 0;

  OGR_L_ResetReading(hLayer);

  // FIX #8: Build a hash-set of this PET's feature IDs so membership
  // lookup is O(1) instead of O(N) per feature.
  std::unordered_set<int> localFeatureSet;
  if (featureIDs && *nFeatures > 0) {
    printf("<<>> here %d\n",pet_rank);
    localFeatureSet.insert(featureIDs, featureIDs + *nFeatures);
  }

  // FIX #11: Reserve coordinate vectors to avoid repeated reallocations.
  // Use nFeatures as a conservative lower-bound hint.
  XCoords.reserve(*nFeatures);
  YCoords.reserve(*nFeatures);

  // Loop through ALL features in the layer (distributed by featureID filter).
  for (int i = 0; i < *nElems; i++) {

    hFeature = OGR_L_GetFeature(hLayer, i);
    if (!hFeature) continue; // FIX #1 variant: guard null feature

    int FID = OGR_F_GetFID(hFeature) + 1;

    hGeom = OGR_F_GetGeometryRef(hFeature);
    if (!hGeom) {
      OGR_F_Destroy(hFeature);
      continue;
    }

    OGRwkbGeometryType geomType = wkbFlatten(OGR_G_GetGeometryType(hGeom));

    // FIX #1: Initialize nFTRpoints to 0 and resolve geometry type once.
    // Skip unknown/unhandled types rather than using garbage nFTRpoints.
    int nFTRpoints = 0;
    OGRGeometryH fGeom = nullptr; // FIX #3: always initialised

    if (geomType == wkbPolygon) {
      fGeom      = OGR_G_GetGeometryRef(hGeom, 0);
      nFTRpoints = fGeom ? OGR_G_GetPointCount(fGeom) - 1 : 0;
    } else if (geomType == wkbMultiPolygon) {
      countMultiPolygon(hGeom, &nFTRpoints);
    } else if (geomType == wkbLineString) {
      nFTRpoints = OGR_G_GetPointCount(hGeom);
    } else if (geomType == wkbPoint) {
      nFTRpoints = OGR_G_GetPointCount(hGeom);
    } else {
      // Unhandled geometry type — skip cleanly instead of using
      // an uninitialised nFTRpoints value (fix #1).
      OGR_F_Destroy(hFeature);
      continue;
    }

    // FIX #8: O(1) set lookup replaces O(N) valueinarray scan.
    if (localFeatureSet.count(FID)) {

      localpoints += nFTRpoints;

      // Get element centroid coordinates (2D assumed).
      // FIX #2: Create Cpoint only when needed and destroy it afterwards.
      OGRGeometryH Cpoint = OGR_G_CreateGeometry(wkbPoint);
      OGR_G_Centroid(hGeom, Cpoint);
      elemCoords.push_back(OGR_G_GetX(Cpoint, 0));
      elemCoords.push_back(OGR_G_GetY(Cpoint, 0));
      OGR_G_DestroyGeometry(Cpoint); // FIX #2: always destroy

      if (geomType == wkbPolygon) {
        // ADD POLYGON
        processPolygon(fGeom, XCoords, YCoords, elemConn, numElemConn,
                       nodeIDs, totpoints, &nFTRpoints);
      } else if (geomType == wkbMultiPolygon) {
        // BREAK DOWN MULTIPOLYGON AND ADD SUB-POLYGONS
        processMultiPolygon(hGeom, XCoords, YCoords, elemConn, numElemConn,
                            nodeIDs, totpoints, &nFTRpoints);
      } else if (geomType == wkbLineString) {
        // ADD LineString
        processLineString(hGeom, XCoords, YCoords, nodeIDs, &nFTRpoints);
      } else if (geomType == wkbMultiLineString) {
        // BREAK DOWN MULTILINESTRING AND ADD SUB-LineStrings
        processMultiLineString(hGeom, XCoords, YCoords, nodeIDs, &nFTRpoints);
      } else if (geomType == wkbPoint) {
        for (int p = nFTRpoints - 1; p >= 0; p--) {
          XCoords.push_back(OGR_G_GetX(hGeom, p));
          YCoords.push_back(OGR_G_GetY(hGeom, p));
          nodeIDs.push_back(totpoints + p + 1);
        }
      }

      // FIX #12: accumulate totNumElemConn incrementally rather than
      // in a separate loop after the fact.
      if (!numElemConn.empty()) {
        *totNumElemConn += numElemConn.back();
      }
    }

    // FIX #6: OGR_L_ResetReading removed from inside the loop — it was
    // resetting the sequential cursor on every iteration unnecessarily.
    totpoints += nFTRpoints;

    OGR_F_Destroy(hFeature);
  }

  if (localpoints <= 0 || nodeIDs.empty()) { 
    printf("<<>> %d LOCALPOINTS/%d nodeIDs on PE %d\n",localpoints, nodeIDs.size(),pet_rank);
    return; 
  }

  // FIX #4: Allocate nodeCoords for localpoints only, not totpoints.
  // totpoints counts ALL features seen across all PETs; only localpoints
  // were actually filled into XCoords/YCoords on this PET.
  nodeCoords = new double[2 * localpoints];

  int j = 0;
  for (int i = 0; i < localpoints; i++) {
    nodeCoords[j]     = XCoords[i];
    nodeCoords[j + 1] = YCoords[i];
    j += 2;
  }

  *nNodes = localpoints;
  // totNumElemConn already accumulated incrementally above (fix #12).

  // ---- DIAGNOSTICS (shapefile_distributed) ----
  // These mirror the inputs that ESMCI_meshaddelements will receive so that
  // the source of the SIGSEGV at Mesh_Glue.C:810 can be isolated.

  printf("DIAG_SHP %d# === shapefile_distributed output summary ===\n", pet_rank);
  printf("DIAG_SHP %d# nFeatures(local)=%d  nNodes=%d  nElems=%d  totNumElemConn=%d\n",
         pet_rank, *nElems, *nNodes, *nElems, *totNumElemConn);

  // Sanity: numElemConn size must equal nElems
  printf("DIAG_SHP %d# numElemConn.size()=%zu  (expected nElems=%d)\n",
         pet_rank, numElemConn.size(), *nElems);
  if ((int)numElemConn.size() != *nElems) {
    printf("DIAG_SHP %d# *** MISMATCH: numElemConn.size() != nElems — "
           "this will corrupt the elemType array passed to meshaddelements ***\n", pet_rank);
  }

  // Sanity: elemConn.size() must equal totNumElemConn
  printf("DIAG_SHP %d# elemConn.size()=%zu  totNumElemConn=%d\n",
         pet_rank, elemConn.size(), *totNumElemConn);
  if ((int)elemConn.size() != *totNumElemConn) {
    printf("DIAG_SHP %d# *** MISMATCH: elemConn.size() != totNumElemConn ***\n", pet_rank);
  }

  // Sanity: elemCoords.size() must be 2 * nElems (2D centroids)
  printf("DIAG_SHP %d# elemCoords.size()=%zu  (expected 2*nElems=%d)\n",
         pet_rank, elemCoords.size(), 2 * (*nElems));
  if ((int)elemCoords.size() != 2 * (*nElems)) {
    printf("DIAG_SHP %d# *** MISMATCH: elemCoords.size() != 2*nElems ***\n", pet_rank);
  }

  // Sanity: nodeIDs.size() must equal nNodes
  printf("DIAG_SHP %d# nodeIDs.size()=%zu  (expected nNodes=%d)\n",
         pet_rank, nodeIDs.size(), *nNodes);
  if ((int)nodeIDs.size() != *nNodes) {
    printf("DIAG_SHP %d# *** MISMATCH: nodeIDs.size() != nNodes ***\n", pet_rank);
  }

  // Sanity: nodeCoords pointer
  printf("DIAG_SHP %d# nodeCoords ptr=%s\n",
         pet_rank, (nodeCoords != NULL) ? "non-NULL" : "NULL *** BAD ***");

  // Per-element detail: elemType (numElemConn[i]) and connectivity range
  {
    int conn_offset = 0;
    for (int i = 0; i < (int)numElemConn.size() && i < *nElems; i++) {
      int nconn = numElemConn[i];
      printf("DIAG_SHP %d# elem[%d]: numElemConn=%d  conn_offset=%d",
             pet_rank, i, nconn, conn_offset);
      if (nconn <= 0) {
        printf("  *** BAD numElemConn — will cause sum loop at Glue:810 to under/overflow ***");
      }
      // Print the connectivity entries for this element
      printf("  conn=[");
      for (int c = 0; c < nconn && (conn_offset + c) < (int)elemConn.size(); c++) {
        printf("%d", elemConn[conn_offset + c]);
        if (c < nconn - 1) printf(",");
      }
      printf("]\n");
      conn_offset += nconn;
    }
    if (conn_offset != *totNumElemConn) {
      printf("DIAG_SHP %d# *** sum(numElemConn)=%d != totNumElemConn=%d ***\n",
             pet_rank, conn_offset, *totNumElemConn);
    }
  }

  // elemIDs
  printf("DIAG_SHP %d# elemIDs (first %d): [", pet_rank, (int)elemIDs.size());
  for (int i = 0; i < (int)elemIDs.size(); i++) {
    printf("%d", elemIDs[i]);
    if (i < (int)elemIDs.size()-1) printf(",");
  }
  printf("]\n");

  // nodeIDs (first up to 20)
  {
    int show = (int)nodeIDs.size() < 20 ? (int)nodeIDs.size() : 20;
    printf("DIAG_SHP %d# nodeIDs (first %d of %zu): [", pet_rank, show, nodeIDs.size());
    for (int i = 0; i < show; i++) {
      printf("%d", nodeIDs[i]);
      if (i < show-1) printf(",");
    }
    if ((int)nodeIDs.size() > show) printf(",...");
    printf("]\n");
  }
  printf("DIAG_SHP %d# === end shapefile_distributed output summary ===\n", pet_rank);
  // ---- END DIAGNOSTICS ----

  return;
}

int processPolygon(
  // Just return specs of a given geometry handle. Should be a linestring or polygon.
  // We don't add to the global spec arrays here because we may need to wedge in
  // polybreaks for ESMF in the case we're in a multipolygon loop.
   OGRGeometryH fGeom,  // Geometry handle for this polygon/linestring
   std::vector<double> &XCoords, // X coords in this polygon/linestring.
   std::vector<double> &YCoords, // X coords in this polygon/linestring.
   std::vector<int> &elemConn,
   std::vector<int> &numelemConn,
   std::vector<int> &elemNodeIDs,
   int    totpoints,
   int    *nPpoints)   // Number of points in this polygon/linestring
{

  // Points of a polygon should start and end the same.
  // This will determine the element params and serve
  // as a check.
  if( OGR_G_GetX(fGeom,0) != OGR_G_GetX(fGeom,*nPpoints) &&
      OGR_G_GetY(fGeom,0) != OGR_G_GetY(fGeom,*nPpoints) )
    {
      printf( "Incomplete polygon. (X,Y)i != (X,Y)e\n" );
      return 1;
    }
  
  // Set elemConns & polyCoords
  for (int i = 0; i < *nPpoints; i++) {
    elemConn.push_back(totpoints+i+1); 
    elemNodeIDs.push_back(totpoints+i+1);
  }
  numelemConn.push_back(*nPpoints);

  // -- Set coords: 
  //    this is done this way because it appears GDAL reads these clockwise
  //    and we need it to be counterclockwise. So the loop is reversed.
  //    Haven't found a way to test for this using GDAL, so this is
  //    a hardwire
  for (int i = *nPpoints-1; i >=0; i--) {
//    printf("<<>> <<>> X,Y: %.2f %.2f\n",OGR_G_GetX(fGeom, i),OGR_G_GetY(fGeom, i));
    XCoords.push_back( OGR_G_GetX(fGeom, i) );
    YCoords.push_back( OGR_G_GetY(fGeom, i) );
  }

// Success (currently, there's no error checking)
  return 0;
}

int processMultiPolygon(
  OGRGeometryH hGeom, 
  std::vector<double> &mXCoords, 
  std::vector<double> &mYCoords, 
  std::vector<int> &elemConn, 
  std::vector<int> &nelemConn, 
  std::vector<int> &elemNodeIDs,
  int    totpoints,
  int    *nPpoints)
{ 
  int j;
  int nGeom;
  //int nPpoints,
  int nPconn;
  int nMPp = 0;

  OGRGeometryH fGeom,mGeom;

//  nMPpoints = 0;

  // get number of geometries in MP
  nGeom = OGR_G_GetGeometryCount(hGeom);
  
  // Loop over nGeom in the multipolygon
  int localtotal = totpoints;
  for (j = 0; j < nGeom; j++) {
    fGeom = OGR_G_GetGeometryRef(hGeom,j);
    mGeom = OGR_G_GetGeometryRef(fGeom,0);
    
    // ADD POLYGON
    if (wkbFlatten(OGR_G_GetGeometryType(fGeom)) == wkbPolygon) {
      int nPpoints = OGR_G_GetPointCount(mGeom)-1;
      processPolygon(mGeom, mXCoords, mYCoords, elemConn, nelemConn, elemNodeIDs, localtotal, &nPpoints);
      localtotal += nPpoints; // Updated number of total points in MP
    } // Or RECURSE INTO SUB MULTIPOLYGON
    else if(wkbFlatten(OGR_G_GetGeometryType(fGeom)) == wkbMultiPolygon) {
      // FIX #5: Replaced exit(1) with error return to allow ESMF/MPI
      // to shut down gracefully. Nested multipolygons are not yet supported.
      ESMC_LogDefault.MsgFoundError(ESMC_RC_NOT_IMPL,
        "Nested MultiPolygon geometry (MultiPolygon of MultiPolygons) is not supported.",
        ESMC_CONTEXT, NULL);
      return 1;
    }

    // Now, we need to expand the arrays to fit the new polygon info, and
    // include a polybreak.
    //
    // -- append a polybreak to the element connectivity vectors
    if (j < nGeom-1) {
      elemConn.push_back(MESH_POLYBREAK_IND); // THIS SHOULD BE MESH_POLYBREAK_IND
      nelemConn.back() = nelemConn.back()+1;
    }

    // Need to collapse the number of connectivity entries, since this is a multipolygon.
    // processPolygon() treats them like individual entries.  
    if ( j > 0 ) {
      nelemConn[nelemConn.size()-2] += nelemConn.back();
      nelemConn.pop_back();
    }

    // -- append the coordinate vectors
    //mXCoords.insert(mXCoords.end(),XCoords.begin(),XCoords.end());
    //mYCoords.insert(mYCoords.end(),YCoords.begin(),YCoords.end());

  }

  return 0;
}

int processLineString(
  // Just return specs of a given geometry handle. Should be a linestring or polygon.
  // We don't add to the global spec arrays here because we may need to wedge in
  // polybreaks for ESMF in the case we're in a multipolygon loop.
   OGRGeometryH fGeom,  // Geometry handle for this polygon/linestring
   std::vector<double> &polyXCoords, // X coords in this polygon/linestring.
   std::vector<double> &polyYCoords, // X coords in this polygon/linestring.
   std::vector<int> &elemNodeIDs,
   int    *nPpoints)   // Number of points in this polygon/linestring
{

  *nPpoints    = OGR_G_GetPointCount(fGeom);

  // -- Set coords: 
  //    this is done this way because it appears GDAL reads these clockwise
  //    and we need it to be counterclockwise. So the loop is reversed.
  //    Haven't found a way to test for this using GDAL, so this is
  //    a hardwire
  for (int i = *nPpoints-1; i >=0; i--) {
    polyXCoords.push_back( OGR_G_GetX(fGeom, i) );
    polyYCoords.push_back( OGR_G_GetY(fGeom, i) );
    elemNodeIDs.push_back(totpoints+i+1);
  }

// Success (currently, there's no error checking)
  return 0;
}

int processMultiLineString(
  OGRGeometryH hGeom, 
  std::vector<double> &mXCoords, 
  std::vector<double> &mYCoords, 
  std::vector<int> &elemNodeIDs,
  int    *nMPpoints)
{ 
  int j;
  int nGeom;
  int nPpoints, nPconn;
  int nMPp = 0;

//  double *XCoords, *YCoords;
  std::vector<double> XCoords;
  std::vector<double> YCoords;

  OGRGeometryH fGeom,mGeom;

//  nMPpoints = 0;

  // get number of geometries in MP
  nGeom = OGR_G_GetGeometryCount(hGeom);
  
  // Loop over nGeom in the multilinestring
  for (j = 0; j < nGeom; j++) {
    fGeom = OGR_G_GetGeometryRef(hGeom,j);
    
    // ADD LINESTRING
    if (wkbFlatten(OGR_G_GetGeometryType(fGeom)) == wkbLineString) {
      processLineString(fGeom, mXCoords, mYCoords, elemNodeIDs, &nPpoints);
      nMPp += nPpoints; // Updated number of total points in MP
    } // Or RECURSE INTO SUB MULTILINESTRING
    else if(wkbFlatten(OGR_G_GetGeometryType(fGeom)) == wkbMultiLineString) {
      // FIX #5: Replaced exit(1) with error return to allow ESMF/MPI
      // to shut down gracefully. Nested multilinestrings are not yet supported.
      ESMC_LogDefault.MsgFoundError(ESMC_RC_NOT_IMPL,
        "Nested MultiLineString geometry (MultiLineString of MultiLineStrings) is not supported.",
        ESMC_CONTEXT, NULL);
      return 1;
    }

    XCoords.clear();
    YCoords.clear();

  }
  *nMPpoints = nMPp;

  return 0;
}

int getLayerInfo( OGRLayerH hL, int *nPoints, int *nGeom) 
{
  int nGeom_tmp = 0;
  int nPnts_tmp = 0;
  OGRFeatureH hFeature;

  // Rewind to the beginning, just in case
  OGR_L_ResetReading(hL);

  // Loop through features in layer and establish extents for allocation
  while( (hFeature = OGR_L_GetNextFeature(hL)) != NULL ) {
    OGRGeometryH hGeom1,hGeom2,hGeom3;
    
    // Get geometry handles
    hGeom1 = OGR_F_GetGeometryRef(hFeature); // looks like this should be a polygon
    
    // Process POLYGON
    if (wkbFlatten(OGR_G_GetGeometryType(hGeom1)) == wkbPolygon) {
      hGeom2 = OGR_G_GetGeometryRef(hGeom1,0);  // and this should be linestring
      nGeom_tmp++;
      *nGeom=nGeom_tmp;
      nPnts_tmp = OGR_G_GetPointCount(hGeom2);
      *nPoints += nPnts_tmp-1;
    } // or process a MULTIPOLYGON
    else if(wkbFlatten(OGR_G_GetGeometryType(hGeom1)) == wkbMultiPolygon) {
      int nGeomMP = OGR_G_GetGeometryCount(hGeom1);
      for (int j = 0; j < nGeomMP; j++) {
	hGeom2 = OGR_G_GetGeometryRef(hGeom1,j); // this should be a polygon
	// Process sub-polygon
	if (wkbFlatten(OGR_G_GetGeometryType(hGeom2)) == wkbPolygon) {
	  hGeom3 = OGR_G_GetGeometryRef(hGeom2,0);  // and this should be linestring
	  nGeom_tmp++;
	  *nGeom=nGeom_tmp;
	  nPnts_tmp = OGR_G_GetPointCount(hGeom3);
	  *nPoints += nPnts_tmp-1; // subtract one so we don't repeat the point that closes the ring
	} // Or RECURSE INTO SUB MULTIPOLYGON
	else if(wkbFlatten(OGR_G_GetGeometryType(hGeom2)) == wkbMultiPolygon) {
	  // FIX #5: Replaced exit(1) with error return.
	  ESMC_LogDefault.MsgFoundError(ESMC_RC_NOT_IMPL,
            "Nested MultiPolygon geometry (MultiPolygon of MultiPolygons) is not supported.",
            ESMC_CONTEXT, NULL);
	  OGR_F_Destroy(hFeature);
	  return 1;
	}
      }
    }
    // Process Linestring
    if (wkbFlatten(OGR_G_GetGeometryType(hGeom1)) == wkbLineString) {
      nGeom_tmp++;
      *nGeom=nGeom_tmp;
      nPnts_tmp = OGR_G_GetPointCount(hGeom1);
      *nPoints += nPnts_tmp;
    } // or process a MULTILINESTRING
    else if(wkbFlatten(OGR_G_GetGeometryType(hGeom1)) == wkbMultiLineString) {
      int nGeomMP = OGR_G_GetGeometryCount(hGeom1);
      for (int j = 0; j < nGeomMP; j++) {
	hGeom2 = OGR_G_GetGeometryRef(hGeom1,j); // this should be a linestring
	// Process sub-string
	if (wkbFlatten(OGR_G_GetGeometryType(hGeom2)) == wkbLineString) {
	  nGeom_tmp++;
	  *nGeom=nGeom_tmp;
	  nPnts_tmp = OGR_G_GetPointCount(hGeom2);
	  *nPoints += nPnts_tmp;
	} // Or RECURSE INTO SUB MULTISTRING
	else if(wkbFlatten(OGR_G_GetGeometryType(hGeom2)) == wkbMultiLineString) {
	  // FIX #5: Replaced exit(1) with error return.
	  ESMC_LogDefault.MsgFoundError(ESMC_RC_NOT_IMPL,
            "Nested MultiLineString geometry (MultiLineString of MultiLineStrings) is not supported.",
            ESMC_CONTEXT, NULL);
	  OGR_F_Destroy(hFeature);
	  return 1;
	}
      }
    }
    OGR_F_Destroy(hFeature);
  }
  return 0;
}

bool valueinarray(int val, int *arr, int n) {
    for(size_t i = 0; i < n; i++) {
        if(arr[i] == val)
            return true;
    }
    return false;
}

int countMultiPolygon(
  OGRGeometryH hGeom, 
  int    *nMPpoints )
  { 
  OGRGeometryH fGeom,mGeom;

  // get number of geometries in MP
  int nGeom = OGR_G_GetGeometryCount(hGeom);
  
  // Loop over nGeom in the multipolygon
  for (int j = 0; j < nGeom; j++) {
    fGeom = OGR_G_GetGeometryRef(hGeom,j);
    mGeom = OGR_G_GetGeometryRef(fGeom,0);
    
    // Get points from member polygon
    if (wkbFlatten(OGR_G_GetGeometryType(fGeom)) == wkbPolygon) {
      int nPpoints = OGR_G_GetPointCount(mGeom)-1;
      *nMPpoints += nPpoints; // Updated number of total points in MP
    } // Or RECURSE INTO SUB MULTIPOLYGON
    else if(wkbFlatten(OGR_G_GetGeometryType(fGeom)) == wkbMultiPolygon) {
      // FIX #5: Replaced exit(1) with error return.
      ESMC_LogDefault.MsgFoundError(ESMC_RC_NOT_IMPL,
        "Nested MultiPolygon geometry (MultiPolygon of MultiPolygons) is not supported.",
        ESMC_CONTEXT, NULL);
      return 1;
    }
  }

  return 0;
}

struct NODE_INFO {
  int node_id;
  int local_elem_conn_pos;

  bool operator< (const NODE_INFO &rhs) const {
    return node_id < rhs.node_id;
  }

};

// Note that local_elem_conn is base-1 (as expected by mesh create routines)
void convert_global_elem_conn_to_local_elem_info(int num_local_elem, int tot_num_elem_conn, int *num_elem_conn, int *global_elem_conn, int*& local_elem_conn) {
  // Init output
  local_elem_conn=NULL;

  // If nothing to do leave
  if (tot_num_elem_conn < 1) return;

  // Allocate conversion list
  NODE_INFO *convert_list=new NODE_INFO[tot_num_elem_conn];

  // Copy global elem connection info into conversion list
  int num_node_conn=0; // Number of connections that are nodes (vs. polybreak)
  for (int i=0; i<tot_num_elem_conn; i++) {

    // Skip polygon break entries 
    if (global_elem_conn[i] == MESH_POLYBREAK_IND) continue;
    
    // Add node entiries to conversion list
    convert_list[num_node_conn].node_id=global_elem_conn[i];
    convert_list[num_node_conn].local_elem_conn_pos=i;
    num_node_conn++;
  }

  // Sort list by node_id, to make it easy to find unique node_ids
  std::sort(convert_list,convert_list+num_node_conn);

  // Count number of unique node ids in  convert_list
  int num_unique_node_ids=1;                 // There has to be at least 1, 
  int prev_node_id=convert_list[0].node_id;  // because we leave if < 1 above
  for (int i=1; i<num_node_conn; i++) {

    // If not the same as the last one count a new one
    if (convert_list[i].node_id != prev_node_id) {
      num_unique_node_ids++;
      prev_node_id=convert_list[i].node_id;
    }
  }

  // Allocate local elem conn
  local_elem_conn=new int[tot_num_elem_conn];

  // Set to polybreak value so that it's in the correct places
  // after the code below fills in the node connection values
  for (int i=0; i<tot_num_elem_conn; i++) {
    local_elem_conn[i]=MESH_POLYBREAK_IND;
  }
  
  // Translate convert_list to node_ids and local_elem_conn
  int node_ids_pos=0;                             // There has to be at least 1, 
  local_elem_conn[convert_list[0].local_elem_conn_pos]=node_ids_pos+1; // +1 to make base-1
  for (int i=1; i<num_node_conn; i++) {

    // Add an entry for this in local_elem_conn
    local_elem_conn[convert_list[i].local_elem_conn_pos]=node_ids_pos+1; // +1 to make base-1
  }


  // Get rid of conversion list
  delete [] convert_list;    
}

#endif // ifdef ESMF_GDAL

