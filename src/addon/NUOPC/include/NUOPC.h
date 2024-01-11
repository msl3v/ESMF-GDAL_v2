// $Id$
//
// Earth System Modeling Framework
// Copyright (c) 2002-2024, University Corporation for Atmospheric Research,
// Massachusetts Institute of Technology, Geophysical Fluid Dynamics
// Laboratory, University of Michigan, National Centers for Environmental
// Prediction, Los Alamos National Laboratory, Argonne National Laboratory,
// NASA Goddard Space Flight Center.
// Licensed under the University of Illinois-NCSA License.
//-----------------------------------------------------------------------------

#include "ESMC.h"

//-----------------------------------------------------------------------------
// This file is part of the pure C public NUOPC API
//-----------------------------------------------------------------------------

#ifndef NUOPC_H
#define NUOPC_H

// TODO: access these string constants via ISO_C interop. from Fortran definition
const char *label_InternalState = "ModelBase_InternalState";
const char *label_Advance = "ModelBase_Advance";
const char *label_AdvanceClock = "ModelBase_AdvanceClock";
const char *label_CheckImport = "ModelBase_CheckImport";
const char *label_SetRunClock = "ModelBase_SetRunClock";
const char *label_TimestampExport = "ModelBase_TimestampExport";
const char *label_Finalize = "ModelBase_Finalize";
const char *label_Advertise = "ModelBase_Advertise";
const char *label_ModifyAdvertised = "ModelBase_ModifyAdvertised";
const char *label_RealizeProvided = "ModelBase_RealizeProvided";
const char *label_AcceptTransfer = "ModelBase_AcceptTransfer";
const char *label_RealizeAccepted = "ModelBase_RealizeAccepted";
const char *label_SetClock = "ModelBase_SetClock";
const char *label_DataInitialize = "ModelBase_DataInitialize";

int NUOPC_CompDerive(
  ESMC_GridComp     ,                           // in
  void (*userRoutine)(ESMC_GridComp, int *)     // in
);

int NUOPC_CompSpecialize(
  ESMC_GridComp,                                // in
  char *,                                       // in
  void (*specLabel)(ESMC_GridComp, int *)       // in
);

void NUOPC_ModelSetServices(ESMC_GridComp, int *);

#endif  // NUOPC_H
