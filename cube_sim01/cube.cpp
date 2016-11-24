// Copyright (C) 2008-2016 Johan Jansson and Niclas Jansson as main authors.

// Licensed under the GNU LGPL Version 2.1.
//
// This program solves the incompressible Navier-Stokes equations using
// least-squares-Galerkin stabilized FEM with a Schur-preconditioned fixed-point
// iteration between the momentunm and pressure equations and a do-nothing adaptive
// method.

#include <dolfin.h>
#include <dolfin/main/MPI.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/fem/UFC.h>
#ifdef ENABLE_UFL 
#include "../ufc2/NSEMomentum3D.h"
#include "../ufc2/NSEContinuity3D.h"
#include "../ufc2/NSEDualMomentum3D.h"
#include "../ufc2/NSEDualContinuity3D.h"
#include "../ufc2/NSEErrRepMomentum3D.h"
#include "../ufc2/NSEErrRepContinuity3D.h"
#include "../ufc2/Drag3D.h"
#include "../ufc2/NSEH1.h"
#include "../ufc2/NSEH12.h"
#include "../ufc2/NSEH1Momentum3D.h"
#include "../ufc2/NSEH1Continuity3D.h"
#include "../ufc2/NSEH1MomentumGlobal3D.h"
#include "../ufc2/NSEH1ContinuityGlobal3D.h"
#include "../ufc2/NSEMomentumResidual3D.h"
#include "../ufc2/NSEContinuityResidual3D.h"
#include "../ufc2/NSEMomentumResidualGlobal3D.h"
#include "../ufc2/NSEContinuityResidualGlobal3D.h"
#include "../ufc2/NSEErrEst.h"
#include "../ufc2/NSEErrEstGlobal.h"
#else
#include "../ufc2/NSEMomentum3D.h"
#include "../ufc2/NSEContinuity3D.h"
#endif

#include "../dolfin/NodeNormal.h"
#include "../dolfin/SpaceTimeFunction.h"
#include "../dolfin/SlipBC.h"

#include <ostream>
#include <iomanip>
#include <cstring>
#include <sstream>
#include <string>
#include <algorithm>
#include <map>
#include <mpi.h>

#include "../NSESolver.h"

using namespace dolfin;

int main(int argc, char* argv[])
{
  // Create mesh
  Mesh mesh("mesh.bin");
  //Mesh mesh("mesh.xml");

  NSESolver solver(mesh);
  solver.run();

  return 0;
}
