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

int main(int argc, char* argv[])
{
  // Create mesh
  Mesh mesh("mesh.bin");
  //Mesh mesh("mesh.xml");
  
      NodeNormal *nn;
      nn = new NodeNormal(mesh);

      // Create boundary conditions
      Inflow *inflow;
      DualInflow *dinflow;
      Outflow *outflow;
      ThetaDrag *thetadrag;
      ThetaLift *thetalift;
      SlipMarker *sm;
      DirichletBoundary *dboundary;
      OutflowBoundary *oboundary;
      InflowBoundary *iboundary;
      SlipBoundary *sboundary;
      AllBoundary *aboundary;
      DirichletBC *bc_in_m;
      DirichletBC *bc_c;
      SlipBC *slipbc_m;
      DirichletBC *bc_dm;

      Array<BoundaryCondition*> bcs_m;
      Array<BoundaryCondition*> bcs_dm;

      inflow = new Inflow(mesh);
      dinflow = new DualInflow(mesh);
      outflow = new Outflow(mesh);
      thetadrag = new ThetaDrag(mesh);
      thetalift = new ThetaLift(mesh);
      sm = new SlipMarker(mesh);
      dboundary = new DirichletBoundary;
      oboundary = new OutflowBoundary;
      iboundary = new InflowBoundary;
      sboundary = new SlipBoundary;
      aboundary = new AllBoundary;
      bc_in_m = new DirichletBC(*inflow, mesh, *iboundary);
      bc_c = new DirichletBC(*outflow, mesh, *oboundary);
      slipbc_m = new SlipBC(mesh, *sboundary, *nn);
      bc_dm = new DirichletBC(*dinflow, mesh, *aboundary);

      bcs_m.push_back(bc_in_m);
      bcs_m.push_back(slipbc_m);

      bcs_dm.push_back(bc_dm);

  NSESolver solver(
      mesh,
      bcs_m,
      bcs_dm,
      bc_c
      );
  solver.run();

  return 0;
}
