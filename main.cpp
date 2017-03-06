// Copyright (C) 2008-2016 Johan Jansson and Niclas Jansson as main authors.
// Modified by Niyazi Cem Degirmenci
// Modified by Daniel Castanon Quiroz
// Last changed: 2017-01-06


// Licensed under the GNU LGPL Version 2.1.
//
// This program solves the incompressible Navier-Stokes equations using
// least-squares-Galerkin stabilized FEM with a Schur-preconditioned fixed-point
// iteration between the momentunm and pressure equations and a do-nothing adaptive
// method.

// ** Added Functionality**: Two phase flow and ALE.


#include <dolfin.h>
#include <dolfin/main/MPI.h>
#include <dolfin/config/dolfin_config.h>
#include <dolfin/fem/UFC.h>
#ifdef ENABLE_UFL 
#include "ufc2/NSEMomentum3D.h"
#include "ufc2/NSEDensity3D.h"
#include "ufc2/NSEContinuity3D.h"
#include "ufc2/NSEResidualSC3D.h"
#include "ufc2/NSEDualMomentum3D.h"
#include "ufc2/NSEDualContinuity3D.h"
#include "ufc2/NSEErrRepMomentum3D.h"
#include "ufc2/NSEErrRepContinuity3D.h"
#include "ufc2/Drag3D.h"
#include "ufc2/ProjectDensity3D.h"
#else
#include "ufc2/NSEMomentum3D.h"
#include "ufc2/NSEContinuity3D.h"
#endif

#include "dolfin/NodeNormal.h"
#include "dolfin/SpaceTimeFunction.h"
#include "dolfin/SlipBC.h"


#include <dolfin/fem/UFC.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/mesh/Vertex.h>
#include "dolfin/LaplacianSmoother.h"
#include <dolfin/mesh/Mesh.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/MeshFunction.h>
#include <dolfin/mesh/MeshData.h>


#include "ufc2/L2ProjPfromM.h"
#include "ufc2/L2ProjUfromM.h"


#include <ostream>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <sstream>
#include <string>
#include <algorithm>
#include <map>

#include <mpi.h>

using namespace dolfin;

//////////////////////////////////////////////////////////////////////////////////////////////
/*********************************************************************************************
//Description of the Program:

This code solves two phase flow with ALE when an solid body is immersed with an initial velocity.

1) The solid body moves accordingly with Newton's Second law only in the Y-Direction. 
So it suffices to compute the lift force each time, and move the solid body accordingly.

2)  Slip Boundaries  are imposed weakly in the  ufc/NSEbase.ufl  file
    and are equal to the mesh velocity in all boundaries (interior and exterior).

3) The immersed solid is a cube centered at the origin and with ratio 0.5  (in each axis, not diagonally)

4) The whole domain is a cube  centered at the origin and with ratio 8.0  (in each axis, not diagonally)

5) **Important: the computional domain needs to be big to have zero slip conditions at the outside boundary.
If not is possible that the whole mesh moves, i.e. the LaplacianSmoother moves all  vertices, and the
fluid velocity is not zero at the outside boundary.

6) In the class InitialDensity we set how the densities are at t=0.

7) Note about meshes:
   in the directories meshes/one_phase and meshes/two_phase the same computational domain is meshed.
   However, the one in the two_phase directory has 500,000 cells, and is meshed properly for this problem.
   The meshes in meshes/one_phase are for development and quick experimentation. In both cases the gmsh
   file  (.geo) is also distributed.


/******************************************************************************************/
// ** User Parameters**
//Coordinates of the full domain
const real xmin = -4.0;
const real xmax =  4.0;
const real ymin = -4.0;
const real ymax =  4.0;
const real zmin = -4.0;
const real zmax =  4.0;

const real rhoMin = 1.0e-3;//min density: if this is  modified, it  needs to be modified in ufc2/NSEbase.ufl as well
const real rhoMax = 1.0;   //max density: if this is  modified, it  needs to be modified in ufc2/NSEbase.ufl as well

real obstacle_yvel   = 0.0;//Initial velocity  of the obstacle (immersed body) in Y
const real obstacle_rho   =   0.55;//density   of the obstacle
const real obstacle_vol   =   std::pow(0.5,3);//volume of the obstacle
const real gravity = -9.81;


const  real T = 1.0;//Final time
const  real primal_T =T;
const  int  no_samples = 200;//total of *.bin samples in [0,T]
const  real nu = 1.0e-5;//viscosity in Momentum Equation
const  real dual_T = T/2.0;  /*Dual functionality does not work*/


/////////////////////////////////////////////////////////////////////////////////////////////
//Parameters
const real obstacle_mass  =   obstacle_vol*obstacle_rho;//mass  of the obstacle
real bmarg = 1.0e-3 + DOLFIN_EPS;



//Compute XX (coordinates as function) given mesh
void computeX(Function& XX, Mesh& mesh, Form* aM, uint index);

//deform Mesh given XX
void deform(Function& XX, Mesh& mesh, Form* CG1vform, uint CG1vindex);


//set the mesh velocity and the boundary equal to (vx,vy,vz)
void setMeshVelocity(Mesh &mesh,
                     MeshFunction<bool> &MeshBoundaryMarker,
		     Function &MeshVelBoundaryFunction,
		     Vector &MeshVelBoundaryFunctionx,
		     Form* &CG1vform,
		     uint CG1vindex,
		     BoundaryMesh &boundary,
		     MeshFunction<uint>* &vertex_map,
		     real vx, real vy, real vz);

//refresh the normal and tangential vectors of the mesh
void ComputeTangentialVectors(Mesh& mesh,  Vector& tau_1, 
			      Vector& tau_2, Vector& normal,

			      Form& form, NodeNormal& node_normal);

//////////////////////////////////////////////////////////////////////////////////////////////
//For adaptivity(not used in this test)
void merge(real *a,real *b,real *res,int an,int bn);
void ComputeLargestIndicators_cell(Mesh& mesh, Vector& e_indx, std::vector<int>& cells,
				   real percentage);
void ComputeLargestIndicators_eind(Mesh& mesh, Vector& e_indx, std::vector<int>& cells,
				   real percentage);
void ComputeRefinementMarkers(Mesh& mesh, real percentage, Vector& e_indx,
			      MeshFunction<bool>& cell_refinement_marker);
//////////////////////////////////////////////////////////////////////////////////////////////
//Computes the mean of a function (not used in this test)
void ComputeMean(Mesh& mesh, Function& vmean, Function& v, 
		 Form& form, uint* indices, uint* c_indices);


//-----------------------------------------------------------------------------

int main(int argc, char **argv)
{

  dolfin_init(argc,argv);

  dolfin_set("output destination","silent");
  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");


  // Function for no-slip boundary condition for velocity
  class Noslip : public Function
  {
  public:

    Noslip(Mesh& mesh) : Function(mesh) {}

    void eval(real* values, const real* x) const
    {
      values[0] = 0.0;
      values[1] = 0.0;
      values[2] = 0.0;

    }
  };

  // Function for no-slip boundary condition for velocity
  class DualNoslip : public Function
  {
  public:

    DualNoslip(Mesh& mesh) : Function(mesh) {}

    void eval(real* values, const real* x) const
    {
      values[0] = 0.0;
      values[1] = 0.0;
      values[2] = 0.0;
    }
  };

  // Function for no-slip boundary condition for velocity
  class Outflow : public Function
  {
  public:

    Outflow(Mesh& mesh) : Function(mesh) {}

    void eval(real* values, const real* x) const
    {
      values[0] = 0.0;
    }

  };

  // Function for no-slip boundary condition for velocity
  class DensityOutflow : public Function
  {
  public:

    DensityOutflow(Mesh& mesh) : Function(mesh) {}

    void eval(real* values, const real* x) const
    {
      values[0] = rhoMin;
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
      return on_boundary && false;//Slip Boundaries are imposed weakly in ufc/NSEbase.ufl 
                                  //and are equal to the mesh velocity uw
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
      return x[1] >= ymax - bmarg && on_boundary;
    }
  };

  // Function for no-slip boundary condition for velocity                                                                              
  //Used for weak boundary conditions
  class SlipMarker : public Function
  {
  public:

    SlipMarker(Mesh& mesh) : Function(mesh) {}

    void eval(real* values, const real* x) const
    {
      values[0] = 1.0;
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

      double Xmin=xmin+ bmarg ,Xmax=xmax- bmarg;
      double Ymin=ymin+ bmarg ,Ymax=ymax- bmarg;
      double Zmin=zmin+ bmarg ,Zmax=ymax- bmarg;

      if ( (x[0]> Xmin) &&  (x[0]  < Xmax) && (x[1]  > Ymin) &&  (x[1] < Ymax) && (x[2]> Zmin) &&  (x[2] < Zmax))
	      values[0] = 1.0;

      
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

      double Xmin=xmin+ bmarg ,Xmax=xmax- bmarg;
      double Ymin=ymin+ bmarg ,Ymax=ymax- bmarg;
      double Zmin=zmin+ bmarg ,Zmax=ymax- bmarg;

      if ( (x[0]> Xmin) &&  (x[0]  < Xmax) && (x[1]  > Ymin) &&  (x[1] < Ymax) && (x[2]> Zmin) &&  (x[2] < Zmax))
	      values[1] = 1.0;


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


  class Force : public Function
  {
  public:

    Force(Mesh& mesh) : Function(mesh) {}

    void eval(real* values, const real* x) const
    {
      values[0] = 0.0;
      values[1] = gravity;
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

  class InitialDensity : public Function
  {
  public:

    InitialDensity(Mesh& mesh) : Function(mesh) {}
    void eval(real* values, const real* x) const
    {

      values[0] = rhoMin;

      if(x[1] <= 0)
	values[0] = rhoMax;
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

  // Function for no-slip boundary condition for velocity
  //For dual computations
  class PsiMomentum : public Function
  {
  public:

    PsiMomentum(Mesh& mesh) : Function(mesh) {}

    void eval(real* values, const real* x) const
    {
      values[0] = 0.0;
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

  // Function for no-slip boundary condition for velocity
  //For dual computations
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

  // Function for no-slip boundary condition for velocity
  //For dual computations

  class BPsiContinuity : public Function
  {
  public:

    BPsiContinuity(Mesh& mesh) : Function(mesh) {}
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




  // Create mesh
  Mesh mesh("mesh.bin");

  Assembler assembler(mesh);

  message("Running on %d %s", dolfin::MPI::numProcesses(), 
	  (dolfin::MPI::numProcesses() > 1 ? "nodes" : "node"));
  message("Global number of vertices: %d", 
	  (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numVertices() : mesh.numVertices()));
  message("Global number of cells: %d", 
	  (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numCells() : mesh.numCells()));


  message("Global number of vertices: %d", 
	  (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numVertices() : mesh.numVertices()));
  message("Global number of cells: %d", 
	  (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numCells() : mesh.numCells()));

  MeshSize h(mesh);
  FacetNormal n(mesh);
  NodeNormal nn(mesh);
  CellVolume cv(mesh);

  // Create boundary condition
  Noslip noslip(mesh);
  DualNoslip dnoslip(mesh);
  DensityOutflow density_outflow(mesh);
  Outflow outflow(mesh);
  ThetaLift thetaLift(mesh);
  ThetaDrag thetaDrag(mesh);

  DirichletBoundary dboundary;
  OutflowBoundary oboundary;
  InflowBoundary iboundary;
  SlipBoundary sboundary;
  AllBoundary aboundary;
  DirichletBC bc_m(noslip, mesh, dboundary);
  DirichletBC bc_in_m(noslip, mesh, iboundary);
  DirichletBC bc_c(outflow, mesh, oboundary);
  DirichletBC bc_rho(density_outflow, mesh, oboundary);
  SlipBC slipbc_m(mesh, sboundary, nn);
  SlipMarker sm(mesh);


  //DirichletBC slipbc_m(noslip, mesh, sboundary);
  DirichletBC bc_dm(dnoslip, mesh, aboundary);

  Array<BoundaryCondition*> bcs_m;
  //bcs_m.push_back(&bc_in_m);
  bcs_m.push_back(&slipbc_m);

  Array<BoundaryCondition*> bcs_rho;
  bcs_rho.push_back(&bc_rho);

  Array<BoundaryCondition*> bcs_dm;
  bcs_dm.push_back(&bc_dm);


  uint *c_indices = 0;
  uint *indices = 0;

  real hmin = h.min();
  cout << "hmin: " << hmin << endl;
  
  InitialDensity rho_i0(mesh);
  Force f(mesh);
  PsiMomentum psim(mesh);
  PsiContinuity psic(mesh);
  BPsiContinuity bpsic(mesh);

  real k = 0.01*hmin;
  //real c1 = 0.2;
  real c1 = 0.02;
  real c2 = 0.01;
  real c3 = 0.0;
  //real c1rho = 0.2;
  real c1rho = 0.02;
  //real c2rho = 0.1;
  real c2rho = 0.05;
  real c3rho = 0.0;
  real c1p = 1.0;


  cout << "Stabilization: " << " c1: " << c1 << " c2: " << c2 << " c1rho: " << c1rho << " c2rho: " << c2rho << endl;

  real rdtol = 1e-3;
  real rdtol_u = 1e-2;
  int maxit = 2000;

  Function u;
  Function u0;
  Function rho;
  Function rho0;
  Function p;
  Function p0;
  Function nuf(mesh, nu);
  Function kf(mesh, k);
  Function k0f(mesh, 0.0);
  Function c1f(mesh, c1);
  Function c2f(mesh, c2);
  Function c3f(mesh, c3);
  Function c1pf(mesh, c1p);
  Function c1rhof(mesh, c1rho);
  Function c2rhof(mesh, c2rho);
  Function c3rhof(mesh, c3rho);
  Function hminf(mesh, hmin);
  Function umean;

  /* for moving the domain */
  /*-------------------------------------------------------*/
  Function p_n1(mesh);
  Function u_n1(mesh);
  Function pu_n1(mesh);

  Vector pu_n1x;
  Vector p_n1x;
  Vector u_n1x;

  L2ProjUfromMBilinearForm bipuf;
  L2ProjUfromMLinearForm Lpuf(pu_n1);
  
  L2ProjPfromMBilinearForm bippf;
  L2ProjPfromMLinearForm Lppf(pu_n1);

  
  Form *CG1form = &Lppf;
  int  CG1index = 0;
  
  Form *CG1vform = &Lpuf;
  int CG1vindex = 0;

  Form *CG1CG1vform = &Lpuf;
  int CG1CG1vindex = 1;

  p_n1.init(mesh,p_n1x, *CG1form, CG1index);
  u_n1.init(mesh,u_n1x, *CG1vform, CG1vindex);
  pu_n1.init(mesh, pu_n1x, *CG1CG1vform, CG1CG1vindex);
   // initial conditions 
  p_n1.vector() = 0;
  u_n1.vector() = 0;


  MeshFunction<bool> MeshBoundaryMarker(mesh, 0);

  Function MeshVelBoundaryFunction;
  Vector MeshVelBoundaryFunctionx; /* doesn't change after init */

  //Mesh Velocity
  //DCQ
  Function w,w0;

  Vector wx,w0x;   
  Function X,Xtmp;
  Vector Xx;
  Vector Xtmpx;

  MeshVelBoundaryFunction.init(mesh, MeshVelBoundaryFunctionx, Lpuf, 0);
  w.init(mesh, wx, Lpuf, 0);
  w0.init(mesh, w0x, Lpuf, 0);
  X.init(mesh, Xx, Lpuf, 0);
  Xtmp.init(mesh, Xtmpx, Lpuf, 0);

  wx = 0;
  w0x = 0;
  MeshVelBoundaryFunctionx = 0;

  //Moving Domains Ends
  /*-------------------------------------------------------*/



  Function R_sc;

  Function rd_u;
  Function rd_rho;
  Function rd_p;

  Function up;
  Function rhop;
  Function pp;
  Function up0;

  Function ei_m;
  Function ei_c;
  Function eij_m;
  Function eij_c;
  Function eif;

  Function tau_1, tau_2, normal;

  Vector ei;

  // Declare primal and dual forms
  Form *a_m, *L_m, *a_rho, *L_rho, *a_c, *L_c;

  NSEMomentum3DBilinearForm ap_m(u, rho, nuf, h, kf, c1f, c2f, c3f, u0, rho0, f, R_sc,w0,w,sm,normal);
  NSEMomentum3DLinearForm Lp_m(u, rho, p, nuf, h, kf, c1f, c2f, c3f, u0, rho0, f, R_sc,w0,w,sm,normal);

  NSEDensity3DBilinearForm ap_rho(u, rho, h, kf, c1rhof, c2rhof, c3rhof, u0, rho0, f, R_sc,w0,w);
  NSEDensity3DLinearForm Lp_rho(u, rho, h, kf, c1rhof, c2rhof, c3rhof, u0, rho0, f, R_sc,w0,w);

  NSEContinuity3DBilinearForm ap_c(h, c1pf);
  NSEContinuity3DLinearForm Lp_c(u, rho, p, h, c1pf, u0, rho0, p0, f);

  NSEResidualSC3DLinearForm LR_sc(u, rho, p, kf, u0, rho0, f, cv);
  //DCQ: No dual
  /*
  NSEDualMomentum3DBilinearForm ad_m(up, p, nuf, h, kf, c1f, u0, hminf, umean);
  NSEDualMomentum3DLinearForm Ld_m(up, p, h, kf, u0, hminf, umean, psim);

  NSEDualContinuity3DBilinearForm ad_c(h);
  NSEDualContinuity3DLinearForm Ld_c(up, psic, bpsic);
  */

  NSEErrRepMomentum3DLinearForm Lrep_m(up, pp, nuf, kf, up0, u);
  NSEErrRepContinuity3DLinearForm Lrep_c(up, up0, p);

  Drag3DFunctional Md(normal, thetaDrag, p);
  Drag3DFunctional Ml(normal, thetaLift, p);

  ProjectDensity3DBilinearForm aP_rho(h, c2f);
  ProjectDensity3DLinearForm LP_rho(rho_i0);

  // Initialize functions
  u.init(mesh, ap_m, 0);
  u0.init(mesh, ap_m, 0);
  rho.init(mesh, ap_rho, 0);
  rho0.init(mesh, ap_rho, 0);
  p.init(mesh, ap_c, 0);
  p0.init(mesh, ap_c, 0);


  
  R_sc.init(mesh, LR_sc, 0);

  rd_u.init(mesh, ap_m, 0);
  rd_rho.init(mesh, ap_rho, 0);
  rd_p.init(mesh, ap_c, 0);

  //umean.init(mesh, ap_m, 12);

  up.init(mesh, ap_m, 0);
  up0.init(mesh, ap_m, 0);
  rhop.init(mesh, ap_rho, 0);
  pp.init(mesh, ap_c, 0);

  ei_m.init(mesh, Lrep_m, 0);
  eij_m.init(mesh, Lrep_m, 0);
  ei_c.init(mesh, Lrep_c, 0);
  eij_c.init(mesh, Lrep_c, 0);
  eif.init(mesh, ei, Lrep_c, 0);

  normal.init(mesh, ap_m, 0);
  tau_1.init(mesh, ap_m, 0);
  tau_2.init(mesh, ap_m, 0);

  p.vector() = 1.0;
  p0.vector() = 1.0;


  ei_m.vector() = 0.0;
  ei_c.vector() = 0.0;

  ComputeTangentialVectors(mesh, (Vector&)tau_1.vector(), (Vector&)tau_2.vector(), (Vector&)normal.vector(), ap_m, nn);

  dolfin_set("PDE linear solver", "iterative");


  // Declare PDE solvers

  Function U;
  Function Rho;
  Function P;

  LinearPDE *pde_m, *pde_rho, *pde_c;

  //DCQ
  LinearPDE pdep_m(ap_m, Lp_m, mesh, bcs_m);
  //LinearPDE pdep_m(ap_m, Lp_m, mesh);
  
 LinearPDE pdep_rho(ap_rho, Lp_rho, mesh, bcs_rho);
  LinearPDE pdep_c(ap_c, Lp_c, mesh, bc_c, cg);

  /*
  //DCQ: No Dual
  LinearPDE pded_m(ad_m, Ld_m, mesh, bcs_dm);
  LinearPDE pded_c(ad_c, Ld_c, mesh, bc_c);
  */

  File file_u("velocity.bin");
  File file_rho("density.bin");
  File file_p("pressure.bin");
  File file_du("dvelocity.bin");
  File file_dp("dpressure.bin");
  File file_m("mesh2.bin");


  file_m << mesh;

  int iteration0 = 0;

  real t = 0; real s = 0;

  SpaceTimeFunction* Up = 0;
  SpaceTimeFunction* Pp = 0;
  /*
  //DCQ 
  //No dual
  Up = new SpaceTimeFunction(mesh, up);
  std::vector<std::string> uprimal_fnames;

  //Up->util_fileList("velocity_v", no_samples, uprimal_fnames);
  //Up->util_addFiles(uprimal_fnames, T);

  Pp = new SpaceTimeFunction(mesh, pp);
  std::vector<std::string> pprimal_fnames;
  //Pp->util_fileList("pressure_v", no_samples, pprimal_fnames);
  //Pp->util_addFiles(pprimal_fnames, T);
  */

  std::string solver = "primal";

  bool coeffchanged = true;

  for(int solver_idx = 0; solver_idx < 2; solver_idx++)
  {
    if(solver_idx == 1)
      solver = "dual";
    
    if(solver == "primal")
    {
      cout << "Starting primal solver" << endl;
      a_m = &ap_m; L_m = &Lp_m; a_rho = &ap_rho; L_rho = &Lp_rho; a_c = &ap_c; L_c = &Lp_c;
      pde_m = &pdep_m; pde_rho = &pdep_rho; pde_c = &pdep_c;
    }
    else
    {
      cout << "Starting dual solver" << endl;
      break;
      //DCQ: No dual
      /*
      a_m = &ad_m; L_m = &Ld_m; a_c = &ad_c; L_c = &Ld_c;
      pde_m = &pded_m; pde_c = &pded_c;
      T = dual_T;
      */
    }

    u.vector() = 0.0;
    u0.vector() = 0.0;
    rho.vector() = 1.0;
    rho0.vector() = 1.0;
    p.vector() = 0.0;
    p0.vector() = 0.0;
    up0.vector() = 0.0;
    rhop.vector() = 0.0;


    
    LinearPDE pde_P_rho(aP_rho, LP_rho, mesh);
    pde_P_rho.solve(Rho);
    rho.vector() = Rho.vector();
    rho0.vector() = rho.vector();

    int stepcounter = 0;
    int sample = 0;
    t = 0;
    /*--------------------------------------------------------------------------------------*/
    /* Moving boundary init */
    BoundaryMesh boundary(mesh);
    MeshFunction<uint>* vertex_map = boundary.data().meshFunction("vertex map");

    int d = mesh.topology().dim();
    UFC ufc(CG1vform->form(), mesh, CG1vform->dofMaps());
    Cell c(mesh, 0);
    uint local_dim = c.numEntities(0);
    uint *idx  = new uint[d * local_dim];
    uint *id  = new uint[d * local_dim];
    real *XX_block = new real[d * local_dim];
    //DCQ
    //Mesh Velocity
    File MeshVelFile("meshvelocity.bin");
    /*--------------------------------------------------------------------------------------*/
    /*To move the obstacle accordingly */
    real lift = 0.0;
    real drag = 0.0;
    /*--------------------------------------------------------------------------------------*/


    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //Init Obstacle/Mesh velocity
    {
      cout << "Initial Obstacle vel: " << obstacle_yvel << " t = " << t << endl;
      cout << "Init Mesh Velocity" << endl;
      setMeshVelocity(mesh,MeshBoundaryMarker, MeshVelBoundaryFunction,
		      MeshVelBoundaryFunctionx,CG1vform,CG1vindex,boundary, vertex_map,
		      0.0, obstacle_yvel, 0.0);
      //smooth but do not move mesh
      NodeNormal nodenormal(mesh);
      LaplacianSmoother *lsmoother= new LaplacianSmoother(mesh,  MeshBoundaryMarker, &MeshVelBoundaryFunctionx, nodenormal);
      cout << "Calling smoother at t: "<< t <<", stepcounter: "<<stepcounter<<endl;
      w0x=wx;
      w0x.apply();
      w0.sync_ghosts();
      lsmoother->smooth(wx, true, false);            
      delete lsmoother;
    }
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    NodeNormal nodenormal(mesh);
    LaplacianSmoother *lsmoother;
    // Time-stepping
    while(t <= T)
    {


      cout << "Starting timestepping" << endl;
  

      s = primal_T - t;
      
      if(solver == "dual")
      {
	/*
	//DCQ: No Dual
	cout << "eval dual" << endl;
	Up->eval(s);
	Pp->eval(s);
	cout << "eval dual done" << endl;
	*/
      }
      
      if(t <= 0.01 * T)
      {
	k = 0.02*hmin;
      }
      else
      {
	k = 0.02*hmin;
      }

      //k = 0.02*hmin;
      
      //nuf.init(mesh, nu);
      kf.init(mesh, k);
      //c1f.init(mesh, c1);
      real stimer = time();


      
      // Fixed-point iteration
      for(int i = 0; i < maxit; i++)
      {
	cout << "Fixed-point iteration" << endl;
	real rd_u_norm, rd_rho_norm, rd_p_norm;
	//ComputeMean(mesh, umean, u, ap_m, indices, c_indices);

	// Compute residual for shock-capturing
	assembler.assemble(R_sc.vector(), LR_sc);

	real timer = time();

	cout << "Solving momentum" << endl;
	pde_m->solve(U);
	rd_u.vector() = U.vector();
	rd_u.vector() -= u.vector();
	rd_u_norm = rd_u.vector().norm(l2);
	u.vector() = U.vector();
	
	cout << "Solving density" << endl;
	pde_rho->solve(Rho);
	rd_rho.vector() = Rho.vector();
	rd_rho.vector() -= rho.vector();
	rd_rho_norm = rd_rho.vector().norm(l2);
	rho.vector() = Rho.vector();

	cout << "Solving continuity" << endl;
	pde_c->solve(P);
	p.vector() = P.vector();
	rd_p.vector() = p.vector();
	rd_p.vector() -= p0.vector();
	rd_p_norm = rd_p.vector().norm(l2);

	cout << "Iteration info: " << "Unorm: " << U.vector().norm(linf) << " Rhonorm: " << Rho.vector().norm(linf) << " Pnorm: " << P.vector().norm(linf) << " Uincr: " <<  rd_u_norm / u.vector().norm(l2) << " Rhoincr: " << rd_rho_norm / rho.vector().norm(l2) << " Pincr: " <<  rd_p_norm / p.vector().norm(l2) << " k: " << k << " step: " << stepcounter << " t: " << t << " timer: " << time() - timer << endl;
	cout << "iteration: " << i << endl;
	iteration0 = i;
	p0.vector() = p.vector();
 	if(rd_u_norm / u.vector().norm(l2) <= rdtol_u &&
	   rd_rho_norm / rho.vector().norm(l2) <= rdtol &&
	   rd_p_norm / p.vector().norm(l2) <= rdtol ||
	   i >= maxit - 1)
	{
	  cout << "Step info: " << "Unorm: " << U.vector().norm(linf) << " Rhonorm: " << Rho.vector().norm(linf) << " Pnorm: " << P.vector().norm(linf) << " Uincr: " <<  rd_u_norm / u.vector().norm(l2) << " Rhoincr: " << rd_rho_norm / rho.vector().norm(l2) << " Pincr: " <<  rd_p_norm / p.vector().norm(l2) << " k: " << k << " step: " << stepcounter << " iters: " << iteration0 + 1 << " t: " << t << " timer: " << time() - stimer << endl;
	  break;
	}
      }

      u0.vector() = u.vector();
      rho0.vector() = rho.vector();
      up0.vector() = up.vector();
      
      if(solver == "dual")
      {
	cout << "errest" << endl;
	assembler.assemble(eij_m.vector(), Lrep_m);
	ei_m.vector() += eij_m.vector();
	assembler.assemble(eij_c.vector(), Lrep_c);
	ei_c.vector() += eij_c.vector();
	cout << "errest done: " << ei_m.vector().norm(linf) << " " << ei_c.vector().norm(linf) << endl;
      }

      if(solver == "primal")
      {
	lift = assembler.assemble(Ml);
	drag = assembler.assemble(Md);
	cout << "drag: " << drag << " t = " << t << endl;
	cout << "lift: " << lift << " t = " << t << endl;


      }
      
      
      //DCQ      
      //Compute  Total Force in Y-direction on the obstacle
      real  yforce= lift + obstacle_mass*gravity;
      obstacle_yvel  = obstacle_yvel + k*yforce/obstacle_mass; //Newton's 2nd Law and Euler'scheme
      if(dolfin::MPI::processNumber() == 0){
	//setprecision only works if everything is from stdlib
	std::cout <<std::setprecision(6) << "Total YForce: " << yforce << " t = "          << t << std::endl;
	std::cout <<std::setprecision(6) << "Obstacle vel: " << obstacle_yvel << " t = "   << t << std::endl;
      }

      /* fill the boundary conditions*/
      
      cout << "Updating Velocity Mesh" << endl;
      setMeshVelocity(mesh,MeshBoundaryMarker, MeshVelBoundaryFunction,
		      MeshVelBoundaryFunctionx,CG1vform,CG1vindex,boundary, vertex_map,
		      0.0, obstacle_yvel, 0.0);
    

      //save previous velocity
      w0x=wx;
      w0x.apply();
      w0.sync_ghosts();
      //smooth and move mesh
      cout << "Calling smoother at t: "<< t <<", stepcounter: "<<stepcounter<<endl;
      if(stepcounter == 0)
	{
	  lsmoother = new LaplacianSmoother(mesh,  MeshBoundaryMarker, &MeshVelBoundaryFunctionx, nodenormal);
	  lsmoother->smooth(wx, true, false);            
	}
      else
	lsmoother->smooth(wx, false, false);    


      computeX(X, mesh, CG1vform, CG1vindex );
      Xtmpx = wx;
      Xtmpx *= k; 
      Xtmpx += Xx;
      Xtmpx.apply();
      Xtmp.sync_ghosts();
      deform(Xtmp, mesh, CG1vform, CG1vindex);	
      nn.__compute_normal(mesh);
      ComputeTangentialVectors(mesh, (Vector&)tau_1.vector(), (Vector&)tau_2.vector(), (Vector&)normal.vector(), ap_m, nn);

      
      //save data
      if(stepcounter == 0 || t > T*(real(sample)/real(no_samples)))
      {
	if(solver == "primal")
	{
	  file_u << u;
	  file_rho << rho;
	  file_p << p;

	  cout << "Saving data at t: "<< t <<", sample: "<<sample<<", stepcounter: "<<stepcounter<<endl;

	  
	  // Record primal solution
	  std::stringstream number;
	  number << std::setfill('0') << std::setw(6) << sample;
	  
	  // Save primal velocity
	  std::stringstream ufilename;
	  ufilename << "velocity_v" << number.str() << ".bin" << std::ends;
	  File velbinfile(ufilename.str());
	  velbinfile << u.vector();
	  
	  // Save primal pressure
	  std::stringstream pfilename;
	  pfilename << "pressure_v" << number.str() << ".bin" << std::ends;
	  File pbinfile(pfilename.str());
	  pbinfile << p.vector();

	  //Mesh and mesh velocity 
	  std::stringstream filename;
	  filename << "mesh" <<number.str() << ".bin" << std::ends;
	  File meshfile2(filename.str());
	  meshfile2 << mesh;  
	  MeshVelFile << w;
	
	}
	else
	  {
	  file_du << u;
	  file_dp << p;
	}
	
	sample++;
      }

      
      t += k;
      stepcounter++;
    }
    //Free resouces of Moving BD
    delete lsmoother;
    delete idx; 
    delete id; 
    delete XX_block;


    
    cout << "Solver done" << endl;

    if(solver == "dual")
    {
      cout << "Preparing adaptivity" << endl;

      // Adaptive error control
      if(!ParameterSystem::parameters.defined("adapt_algorithm"))
	dolfin_add("adapt_algorithm", "rivara");
      dolfin_set("adapt_algorithm", "rivara");
      if(!ParameterSystem::parameters.defined("output_format"))
	dolfin_add("output_format", "binary");
      dolfin_set("output_format", "binary");
      MeshFunction<bool> cell_marker;
      ei = 0.0;
      ei += ei_m.vector();
      ei += ei_c.vector();
      
      File file_ei("ei.bin");
      file_ei << eif;

      ComputeRefinementMarkers(mesh, 10.0, ei, cell_marker);
      RivaraRefinement::refine(mesh, cell_marker);

      
      File file_rm("rmesh.bin");
      file_rm << mesh;
    }
  }
  

  return 0;
}



void computeX(Function& XX, Mesh& mesh, Form* aM, uint index)
{
  // Copy mesh coordinates into X array/function
  message("entered computeX");
  int d = mesh.topology().dim();
  UFC ufc(aM->form(), mesh, aM->dofMaps());
  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[d * local_dim];
  uint *id  = new uint[d * local_dim];
  real *XX_block = new real[d * local_dim];  
  
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    ufc.update(*cell, mesh.distdata());
    (aM->dofMaps())[index].tabulate_dofs(idx, ufc.cell, cell->index());
    
    uint ii = 0; 
    uint jj = 0;    
    for(uint i = 0; i < d; i++) 
    {    
      for(VertexIterator v(*cell); !v.end(); ++v, ii++) 
      {    
        if (!mesh.distdata().is_ghost(v->index(), 0))  
        {    
          XX_block[jj] = v->x()[i];
          id[jj++] = idx[ii];
        }    
      }    
    }    
    XX.vector().set(XX_block, jj, id); 
  }
  XX.vector().apply();
  XX.sync_ghosts();
  
  delete[] XX_block;
  delete[] idx; 
  delete[] id;
}


void deform(Function& XX, Mesh& mesh, Form* CG1vform, uint CG1vindex)
{
  message("entered deform");
  MeshGeometry& geometry = mesh.geometry();
  
  uint d = mesh.topology().dim();
  uint N = mesh.numVertices();
  if(dolfin::MPI::numProcesses() > 1) 
    N = mesh.distdata().global_numVertices();
  
  UFC ufc(CG1vform->form(), mesh, CG1vform->dofMaps());
  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[d * local_dim];
  uint *id  = new uint[d * local_dim];
  real *XX_block = new real[d * local_dim];  
  
  // Update the mesh
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    ufc.update(*cell, mesh.distdata());
    (CG1vform->dofMaps())[CG1vindex].tabulate_dofs(idx, ufc.cell, cell->index());

    XX.vector().get(XX_block, d * local_dim, idx);

    uint j = 0; 
    for(VertexIterator v(*cell); !v.end(); ++v) 
    {    
      Vertex& vertex = *v;

      {    
        for(unsigned int i = 0; i < d; i++) 
        {    
          geometry.x(vertex.index(), i) = XX_block[i * local_dim + j];
        }    
      }    
      j++; 
    }    
  }

  delete[] XX_block;
  delete[] idx; 
  delete[] id;

  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
}



void ComputeMean(Mesh& mesh, Function& vmean, Function& v, 
		 Form& form, uint* indices, uint* c_indices)
{

  Cell cell_tmp(mesh, 0);
  uint nsd = mesh.topology().dim(); 
  uint local_dim = cell_tmp.numEntities(0);
  UFC ufc(form.form(), mesh, form.dofMaps()); 
  real *v_block = new real[nsd * local_dim * mesh.numCells()];
  real *vmean_block = new real[nsd*mesh.numCells()];
  
  if(!indices) {
    indices = new uint[3 * local_dim * mesh.numCells()];

    uint *ip = &indices[0];
  
    for (CellIterator cell(mesh); !cell.end(); ++cell)
    {
      
      ufc.update(*cell, mesh.distdata());
            
      (form.dofMaps())[0].tabulate_dofs(ip, ufc.cell, cell->index());
            
      ip += 3 * local_dim;
      
    }
  }
  if(!c_indices) {
    c_indices = new uint[nsd * mesh.numCells()];

    uint *cip = &c_indices[0];
    for(CellIterator c(mesh); !c.end(); ++c) {
      ufc.update(*c, mesh.distdata());
      (form.dofMaps())[12].tabulate_dofs(cip, ufc.cell, c->index());
      
      cip += nsd;
    }
  }

  v.vector().get(v_block, nsd * local_dim * mesh.numCells(), indices);

  uint mi = 0;
  real cellmean = 0.0;  
  uint ri = 0;
  for (CellIterator c(mesh); !c.end(); ++c)
  {

    for (uint i = 0; i < nsd; i++) {
      cellmean = 0.0;
      for (VertexIterator n(*c); !n.end(); ++n)
	cellmean += v_block[ri++];
      cellmean /= c->numEntities(0);
      vmean_block[mi++] = cellmean;
    }    
  }
  vmean.vector().set(vmean_block,nsd * mesh.numCells(), c_indices);
  vmean.vector().apply();
  
  delete[] v_block;
  delete[] vmean_block;

}

void ComputeTangentialVectors(Mesh& mesh,  Vector& tau_1, 
			      Vector& tau_2, Vector& normal,
			      Form& form, NodeNormal& node_normal)
{
  UFC ufc(form.form(), mesh, form.dofMaps());
  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[3 * local_dim];
  uint *id  = new uint[3 * local_dim];
  real *tau_1_block = new real[3 * local_dim];  
  real *tau_2_block = new real[3 * local_dim];  
  real *normal_block = new real[3 * local_dim];

  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    
    ufc.update(*cell, mesh.distdata());
    
    (form.dofMaps())[1].tabulate_dofs(idx, ufc.cell, cell->index());
    
    uint ii = 0;
    uint jj = 0;    
    for(uint i = 0; i < 3; i++) 
    {
      for(VertexIterator v(*cell); !v.end(); ++v, ii++) 
      {
	if (!mesh.distdata().is_ghost(v->index(), 0)) 
	{
	  tau_1_block[jj] = node_normal.tau_1[i].get(*v);
	  tau_2_block[jj] = node_normal.tau_2[i].get(*v);
	  normal_block[jj] = node_normal.normal[i].get(*v);
	  id[jj++] = idx[ii];
	}
      }
    }

    tau_1.set(tau_1_block, jj, id);
    tau_2.set(tau_2_block, jj, id);
    normal.set(normal_block, jj, id);
  }

  tau_1.apply();
  tau_2.apply();
  normal.apply();
  delete[] tau_1_block;
  delete[] tau_2_block;
  delete[] normal_block;
  delete[] idx;
  delete[] id;

}

// Comparison operator for index/value pairs
struct less_pair : public std::binary_function<std::pair<int, real>,
					       std::pair<int, real>, bool>
{
  bool operator()(std::pair<int, real> x, std::pair<int, real> y)
  {
    return x.second < y.second;
  }
};


void merge(real *a,real *b,real *res,int an,int bn)
{
  real *ap,*bp,*rp;
  ap=a;
  bp=b;
  rp=res;

  while(ap<a+an && bp<b+bn){ 
    if(*ap <= *bp){
      *rp=*ap;
      ap++;
      rp++;
    }
    else { 
      *rp=*bp;
      rp++;
      bp++;
    }
  }
  if(ap<a+an){
    do
      *rp=*ap;
    while(++rp && ++ap<a+an);
  }
  else{
    do
      *rp=*bp;
    while(++rp && ++bp<b+bn);
  }
}

void ComputeLargestIndicators_cell(Mesh& mesh, Vector& e_indx, std::vector<int>& cells,
						  real percentage)
{
  int N = mesh.numCells();
  int M = std::min((int)(N), 
		   (int)((real) 
			 (dolfin::MPI::numProcesses() > 1 ? 
			  mesh.distdata().global_numCells() : mesh.numCells()) * percentage * 0.01));
  
  if(dolfin::MPI::processNumber() == 1)
    dolfin_set("output destination","terminal");
  message("Computing largest indicators");
  message("percentage: %f", percentage);
  message("N: %d", N);
  message("M: %d", M);
  dolfin_set("output destination","silent");


  std::vector<std::pair<int, real> > indicators(N);
  real eind;
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    int id = (*cell).index();
    std::pair<int, real> p;
    p.first = id;
    uint ci = id;    
    if(dolfin::MPI::numProcesses() > 1)
      ci = mesh.distdata().get_cell_global(ci);
    e_indx.get(&eind, 1, &ci);      
    p.second = eind;    
    indicators[id] = p;
  }

  less_pair comp;
  std::sort(indicators.begin(), indicators.end(), comp);


  real *local_eind = new real[M];
  for(int i = 0; i < M; i++)
  {
    std::pair<int, real> p = indicators[N - 1 - i];
    local_eind[M - 1 - i] = p.second;
  }


  /*
   *  FIXME reduce memory usage
   *  merge only half of the recived data
   */

  uint M_max, M_tot;
  MPI_Allreduce(&M, &M_max, 1, MPI_UNSIGNED, MPI_MAX, dolfin::MPI::DOLFIN_COMM);
  MPI_Allreduce(&M, &M_tot, 1, MPI_UNSIGNED, MPI_SUM, dolfin::MPI::DOLFIN_COMM);

  double *recv_eind = new double[M_max];
  double *global_eind = new double[M_tot];
  double *work = new double[M_tot];

  //  std::vector<double> global_eind;

  MPI_Status status;
  uint src,dest;
  uint rank =  dolfin::MPI::processNumber();
  uint size =  dolfin::MPI::numProcesses();
  uint nm = M;
  int num_recv;
  //  global_eind.insert(global_eind.begin(), local_eind, local_eind + M);
  std::memcpy(global_eind, local_eind, M*sizeof(real));

  for(uint i = 1; i < size; i++) {
    src =(rank - i + size) % size;
    dest = (rank + i) % size;

    MPI_Sendrecv(local_eind, M, MPI_DOUBLE, dest, 0, 
		 recv_eind, M_max, MPI_DOUBLE, src, 0, dolfin::MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_DOUBLE,&num_recv);
    //global_eind.insert(global_eind.end(), recv_eind, recv_eind + num_recv);
    merge(recv_eind, global_eind, work, num_recv, nm);
    std::memcpy(global_eind, work, M_tot * sizeof(real));
    nm += num_recv;
    
  }

  //  std::sort(global_eind.begin(), global_eind.end());
  cells.clear();
  int MM = (int)((real) (dolfin::MPI::numProcesses() > 1 ? 
			 mesh.distdata().global_numCells() : mesh.numCells()) * percentage * 0.01);
  int i = 0;
  for(int j = 0; j < MM; j++) {
    if( local_eind[M - 1 - i] >= global_eind[M_tot - 1 - j] ) {
      std::pair<int, real> p = indicators[N - 1 - i];
      cells.push_back(p.first);
      if( (i++) >= std::min(N, MM)) break;    
    }
  }

  dolfin_set("output destination", "terminal");
  message("%d marked cells on cpu %d", cells.size(), dolfin::MPI::processNumber());
  dolfin_set("output destination", "silent");

  
  delete[] local_eind;
  delete[] recv_eind;
  delete[] global_eind;
  delete[] work;
}


void ComputeLargestIndicators_eind(Mesh& mesh, Vector& e_indx, std::vector<int>& cells,
						  real percentage)
{
  int N = mesh.numCells();
  real eind, sum_e, sum_e_local, max_e, max_e_local, min_e, min_e_local;
  sum_e = sum_e_local = max_e_local = 0.0;
  min_e_local = 1e6;
  
  std::vector<std::pair<int, real> > indicators(N);

  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    int id = (*cell).index();
    std::pair<int, real> p;
    p.first = id;
    uint ci = id;    
    if(dolfin::MPI::numProcesses() > 1)
      ci = mesh.distdata().get_cell_global(ci);
    e_indx.get(&eind, 1, &ci);      
    // Take absolute value
    eind = abs(eind);
    p.second = eind;    
    indicators[id] = p;
    max_e_local = std::max(max_e_local, eind);
    min_e_local = std::min(min_e_local, eind);
    sum_e_local += p.second;
  }

  less_pair comp;
  std::sort(indicators.begin(), indicators.end(), comp);

  MPI_Allreduce(&sum_e_local, &sum_e, 1, MPI_DOUBLE,
		MPI_SUM, dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&max_e_local, &max_e, 1, MPI_DOUBLE, 
		MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&min_e_local, &min_e, 1, MPI_DOUBLE, 
		MPI_MIN, dolfin::MPI::DOLFIN_COMM);

  real threshold = (percentage * 0.01 * sum_e);
  real cutoff = (max_e + min_e) / 2.0;
  real acc_local, acc;
  acc_local = acc = 0.0;

  int iter = 0;
  while ( (fabs(acc - threshold) / threshold )  > 1e-2  && (iter++) < 10)
  {
    cutoff = (max_e + min_e) / 2.0;
    acc = acc_local = 0.0;
    cells.clear();

    for (int i = 0; i < N; i++) 
    {
      std::pair<int, real> p = indicators[N - 1 - i];

      cells.push_back(p.first);
      acc_local += p.second;

      if ( p.second < cutoff )
	break;     
    }

    MPI_Allreduce(&acc_local, &acc, 1, MPI_DOUBLE, 
		  MPI_SUM, dolfin::MPI::DOLFIN_COMM);
        
    ( acc > threshold ? (min_e = cutoff ) : (max_e = cutoff));    
  }
}

void ComputeRefinementMarkers(Mesh& mesh, real percentage, Vector& e_indx,
			      MeshFunction<bool>& cell_refinement_marker)
{

  real error = 0.0;
  //ComputeError(error);

  //message("err: %g", error);
  
  std::vector<int> cells;
  ComputeLargestIndicators_cell(mesh, e_indx, cells, percentage);
    
  cell_refinement_marker.init(mesh, mesh.topology().dim());
  cell_refinement_marker = false;
    
  int M = cells.size();
       
  for(int i = 0; i < M; i++)
  {
    cell_refinement_marker.set(cells[i], true);
  }

}



void setMeshVelocity(Mesh &mesh,
                     MeshFunction<bool> &MeshBoundaryMarker,
		     Function &MeshVelBoundaryFunction,
		     Vector &MeshVelBoundaryFunctionx,
		     Form* &CG1vform, 
		     uint CG1vindex,
		     BoundaryMesh &boundary,
		     MeshFunction<uint>* &vertex_map,
		     real vx, real vy, real vz)
{

  int d = mesh.topology().dim();
  UFC ufc(CG1vform->form(), mesh, CG1vform->dofMaps());
  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[d * local_dim];
  uint *id  = new uint[d * local_dim];
  real *XX_block = new real[d * local_dim];

  for (VertexIterator bv(boundary); !bv.end(); ++bv)
    {
      Vertex mv(mesh, (*vertex_map)(*bv));      
      real totalw[3];

      totalw[0]  = vx;
      totalw[1]  = vy;
      totalw[2]  = vz;
      //area that is surrouding  the boundary that moves
      double Xmin=xmin+bmarg ,Xmax=xmax-bmarg;
      double Ymin=ymin+bmarg ,Ymax=ymax-bmarg;
      double Zmin=zmin+bmarg ,Zmax=zmax-bmarg;
	  
      double xx= mv.x()[0];
      double yy= mv.x()[1];
      double zz= mv.x()[2];
	  
      if ( (xx> Xmin) &&  (xx < Xmax) && (yy > Ymin) &&  (yy < Ymax) && (zz> Zmin) &&  (zz < Zmax))
	{
	  MeshBoundaryMarker.set(mv, true);
	  
	  for (CellIterator cv(mv); !cv.end(); ++cv)
	    {
	      ufc.update(*cv, mesh.distdata());
	      (CG1vform->dofMaps())[CG1vindex].tabulate_dofs(idx, ufc.cell, cv->index());
	      uint ii = 0;
	      uint jj = 0;
	      for (uint i = 0; i < d; i++)
		for (VertexIterator v(*cv); !v.end(); ++v, ii++)
		  {
			if (mv.index() == v->index())
			  {
			    XX_block[jj] = totalw[i];
			    id[jj++] = idx[ii];
			  }
		  }
	      
		  MeshVelBoundaryFunctionx.set(XX_block,jj,id);
	    }
	  
	}	
    }
  
  
      MeshVelBoundaryFunctionx.apply();
      MeshVelBoundaryFunction.sync_ghosts();
      
      //Free resouces of Moving BD
      delete idx; 
      delete id; 
      delete XX_block;
      
}      
