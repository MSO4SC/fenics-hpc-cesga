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

// Temporary defines to avoid the need for c++11
#define constexpr const
#define nullptr NULL


#include "../NSESolver.h"

using namespace dolfin;

constexpr real bmarg = 1.0e-5 + DOLFIN_EPS;
constexpr real robj = 1. - bmarg;

constexpr real xmin = -10.0;
constexpr real xmax = 30.0;
constexpr real ymin = -10.0;
constexpr real ymax = 10.0;
constexpr real zmin = -10.0;
constexpr real zmax = 10.0;

constexpr real xbody = 0;
constexpr real ybody = 0;
constexpr real zbody = 0;

constexpr real tFinal = 5.0;

constexpr real Uin = 1;
constexpr real Uin_dual = 1;

// Inflow velocity
class Inflow : public Function
{
public:
  Inflow(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = Uin;
    values[1] = 0.0;
    values[2] = 0.0;
  }
};

// Sub domain for Inflow boundary condition
class InflowBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return x[0] <= xmin + bmarg && on_boundary;
  }
};

// Inflow velocity for the dual problem
class BodyDualInflow : public Function
{
public:
  BodyDualInflow(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = Uin_dual;
    values[1] = 0.0;
    values[2] = 0.0;
  }
};

// Sub domain for Dirichlet boundary condition
class Body : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return on_boundary and
        fabs(x[0] - xbody) < robj and
        fabs(x[1] - ybody) < robj and
        fabs(x[2] - zbody) < robj;
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

// Sub domain for Outflow boundary condition
class OutflowBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return x[0] >= xmax - bmarg && on_boundary;
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

// Sub domain for Slip boundary condition
class SlipBoundary : public SubDomain
{
  bool inside(const real* x, bool on_boundary) const
  {
    return x[0] >= xmin + bmarg && x[0] <= xmax - bmarg && on_boundary;
  }
};

// Marker and orientation function for drag computation
class ThetaDrag : public Function
{
public:
  ThetaDrag(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;

    if(fabs(x[0] - xbody) < robj &&
       fabs(x[1] - ybody) < robj &&
       fabs(x[2] - zbody) < robj)
    {
      values[0] = 1.0;
    }
  }

  uint rank() const
  {
    return 1;
  }

  uint dim(uint i) const
  {
    return 3;
  }
};

// Marker and orientation function for lift computation
class ThetaLift : public Function
{
public:

  ThetaLift(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;

    if(fabs(x[0] - xbody) < robj &&
       fabs(x[1] - ybody) < robj &&
       fabs(x[2] - zbody) < robj)
    {
      values[1] = 1.0;
    }
  }

  uint rank() const
  {
    return 1;
  }

  uint dim(uint i) const
  {
    return 3;
  }
};

// Marker function for weak boundary condition
class SlipMarker : public Function
{
public:

  SlipMarker(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = 0.0;

    if(x[0] >= xmin + bmarg && x[0] <= xmax - bmarg)
      values[0] = 0.0;
  }
};

// Dual volume source for momentum
class PsiMomentum : public Function
{
public:

  PsiMomentum(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;
  }

  uint rank() const
  {
    return 1;
  }

  uint dim(uint i) const
  {
    return 3;
  }
};

// Dual volume source for continuity
class PsiContinuity : public Function
{
public:

  PsiContinuity(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = 0.0;
  }

  uint rank() const
  {
    return 0;
  }

  uint dim(uint i) const
  {
    return 0;
  }
};

// Dual boundary source for momentum
class BPsiMomentum : public Function
{
public:

  BPsiMomentum(Mesh& mesh) : Function(mesh) {}

  void eval(real* values, const real* x) const
  {
    values[0] = 0.0;
    values[1] = 0.0;
    values[2] = 0.0;
    if(x[0] >= xmin + bmarg && x[0] <= xmax - bmarg &&
       x[1] >= ymin + bmarg && x[1] <= ymax - bmarg &&
       x[2] >= zmin + bmarg && x[2] <= zmax - bmarg)
    {
      values[0] = 0.0;
    }
  }

  uint rank() const
  {
    return 1;
  }

  uint dim(uint i) const
  {
    return 3;
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

  Body body;
  BodyDualInflow bodyDualInflow(mesh);

  OutflowBoundary oboundary;
  Outflow outflow(mesh);

  SlipBoundary sboundary;

#warning "unused variables"
  DirichletBoundary dboundary;

  ThetaDrag thetadrag(mesh);
  ThetaLift thetalift(mesh);
  SlipMarker sm(mesh);

  PsiMomentum psim(mesh);
  PsiContinuity psic(mesh);
  BPsiMomentum bpsim(mesh);

  DirichletBCList dbcs_m;
  dbcs_m.push_back(std::make_pair(&iboundary,&inflow));

  SlipBCList sbcs_m;
  sbcs_m.push_back(&sboundary);

  DirichletBCList dbcs_c;
  dbcs_c.push_back(std::make_pair(&oboundary,&outflow));

  DirichletBCList dbcs_dm;
#warning "order is critical here"
  dbcs_dm.push_back(std::make_pair(&aboundary,&dinflow));
  dbcs_dm.push_back(std::make_pair(&body,&bodyDualInflow));

  NSESolver solver(
      mesh,
      dbcs_m,
      sbcs_m,
      dbcs_c,
      dbcs_dm,
      &thetadrag,
      &thetalift,
      &sm,
      &psim,
      &psic,
      &bpsim
      );
  solver.setT(tFinal);
  solver.run();

  exit(0);

  return 0;
}
