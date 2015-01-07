// $Id$
//
// Earth System Modeling Framework
// Copyright 2002-2014, University Corporation for Atmospheric Research, 
// Massachusetts Institute of Technology, Geophysical Fluid Dynamics 
// Laboratory, University of Michigan, National Centers for Environmental 
// Prediction, Los Alamos National Laboratory, Argonne National Laboratory, 
// NASA Goddard Space Flight Center.
// Licensed under the University of Illinois-NCSA License.
//-----------------------------------------------------------------------------

//-----------------------------------------------------------------------------
// This file is part of the pure C public ESMC API
//-----------------------------------------------------------------------------

//-------------------------------------------------------------------------
// (all lines below between the !BOP and !EOP markers will be included in
//  the automated document processing.)
//-------------------------------------------------------------------------
// these lines prevent this file from being read more than once if it
// ends up being included multiple times

#ifndef ESMC_Interface_H
#define ESMC_Interface_H

//-----------------------------------------------------------------------------
//BOPI
//
// !DESCRIPTION:
//
// The code in this file provides a types used in the public ESMC interfaces of
// ESMF.
//
//EOPI
//-----------------------------------------------------------------------------


#ifdef __cplusplus
extern "C" {
#endif

// Class declaration type
typedef struct{
  char shallowMem[80];
}ESMC_InterfaceInt;

// Class API
int ESMC_InterfaceIntSet(ESMC_InterfaceInt *interfaceIntArg,
  int *arrayArg, int lenArg);

#ifdef __cplusplus
} // extern "C"
#endif

#endif  // ESMC_Interface_H
