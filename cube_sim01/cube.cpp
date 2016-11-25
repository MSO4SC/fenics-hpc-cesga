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

// Inflow velocity
class Inflow : public Function
{
public:

  Inflow(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;

    if(x[0] <= xmax - bmarg)
      values[0] = 1.0;
  }
};

// Inflow velocity for the dual problem
class DualInflow : public Function
{
public:

  DualInflow(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;

    if(fabs(x[0] - 0.0) < robj &&
       fabs(x[1] - 0.0) < robj &&
       fabs(x[2] - 0.0) < robj)
    {
      values[0] = 1.0;
    }
  }
};

// Outflow pressure
class Outflow : public Function
{
public:

  Outflow(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = 0.0;
  }

};

// Sub domain for Dirichlet boundary condition
class AllBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return on_boundary;
  }
};

// Sub domain for Dirichlet boundary condition
class DirichletBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return x[0] <= xmax - bmarg && on_boundary;
  }
};

// Sub domain for Dirichlet boundary condition
class SlipBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return x[0] >= xmin + bmarg && x[0] <= xmax - bmarg && on_boundary;
  }
};

// Sub domain for Dirichlet boundary condition
class InflowBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return x[0] <= xmin + bmarg && on_boundary;
  }
};

// Sub domain for Dirichlet boundary condition
class OutflowBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return x[0] >= xmax - bmarg && on_boundary;
  }
};

int main(int argc, char* argv[])
{
  // Create mesh
  Mesh mesh("mesh.bin");
  
  // Create boundary conditions
  InflowBoundary iboundary;
  Inflow inflow(mesh);

  AllBoundary aboundary;
  DualInflow dinflow(mesh);

  OutflowBoundary oboundary;
  Outflow outflow(mesh);

  SlipBoundary sboundary;

#warning "unused variables"
  DirichletBoundary dboundary;

  ThetaDrag thetadrag(mesh);
  ThetaLift thetalift(mesh);
  SlipMarker sm(mesh);

  DirichletBCList dbcs_m;
  dbcs_m.push_back(std::make_pair(&iboundary,&inflow));

  SlipBCList sbcs_m;
  sbcs_m.push_back(&sboundary);

  DirichletBCList dbcs_c;
  dbcs_c.push_back(std::make_pair(&oboundary,&outflow));

  DirichletBCList dbcs_dm;
  dbcs_dm.push_back(std::make_pair(&aboundary,&dinflow));

  NSESolver solver(
      mesh,
      dbcs_m,
      sbcs_m,
      dbcs_c,
      dbcs_dm
      );
  solver.run();

  return 0;
}
