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

using namespace dolfin;

real bmarg = 1.0e-5 + DOLFIN_EPS;

std::string simcase = "cube";

real xmin;
real xmax;
real ymin;
real ymax;
real zmin;
real zmax;

real robj;


real adapt_percent = 5.;

real T;

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
  
  if(dolfin::MPI::processNumber() == 0)
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
//-----------------------------------------------------------------------------

// Utility function for mesh rotation
void computeX(Function& XX, Form* aM, Mesh& mesh)
{
  // Copy mesh coordinates into X array/function                                                                                                                                                            
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
    (aM->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());

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

// Utility function for mesh rotation
void Rotate(Function& XX, Form* aM, Mesh& mesh, double theta)
{


  MeshGeometry& geometry = mesh.geometry();

  uint d = mesh.topology().dim();
  uint N = mesh.numVertices();
  if(dolfin::MPI::numProcesses() > 1)
    N = mesh.distdata().global_numVertices();
  UFC ufc(aM->form(), mesh, aM->dofMaps());
  Cell c(mesh, 0);
  uint local_dim = c.numEntities(0);
  uint *idx  = new uint[d * local_dim];
  uint *id  = new uint[d * local_dim];
  real *XX_block = new real[d * local_dim];

  // Update the mesh                                                                                                                                                                                        
  for (CellIterator cell(mesh); !cell.end(); ++cell)
  {
    ufc.update(*cell, mesh.distdata());
    (aM->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());

    XX.vector().get(XX_block, d * local_dim, idx);

    std::vector<double> xx, yy, zz;
    uint jj = 0;
    for(VertexIterator v(*cell); !v.end(); ++v)
    {
      for(unsigned int i = 0; i < d; i++)
      {
        if (i==0)
          xx.push_back(XX_block[i * local_dim + jj]);
        if (i==1)
          yy.push_back(XX_block[i * local_dim + jj]);
      }
      jj++;
    }
    uint j = 0;
    for(VertexIterator v(*cell); !v.end(); ++v)
    {
      Vertex& vertex = *v;

      real theta2;

      Point cp(0.35, -0.06, vertex.point()[2]);
      Point pdiff = vertex.point() - cp;

      real r = pdiff.norm();
      if(r < 0.5)
        theta2 = (0.0 - theta)*2*DOLFIN_PI/360.0;
      else if(r >= 0.5 && r < 1.0)
        theta2 = (0.0 - theta)*2*DOLFIN_PI/360.0*(1.0 - r) / (1.0 - 0.5);
      else
        theta2 = 0.0*2*DOLFIN_PI/360.0;

      for(unsigned int i = 0; i < d; i++)
      {
        if (i==0)
          XX_block[i * local_dim + j] = xx[j]*cos(theta2) - yy[j]*sin(theta2);
        if (i==1)
          XX_block[i * local_dim + j] = xx[j]*sin(theta2) + yy[j]*cos(theta2);

        geometry.x(vertex.index(), i) = XX_block[i * local_dim + j];
      }
      j++;
    }
  }

  delete[] XX_block;
  delete[] idx;
  delete[] id;

  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
}


int main(int argc, char* argv[])
{

  real theta = 0.;

  // Parse command-line arguments
  if(argc >= 2)
  {
    simcase = std::string(argv[1]);

  }

  if(argc >= 3)
  {
    theta = atof(argv[2]);
  }

  xmin = -10.0;
  xmax = 30.0;
  ymin = -10.0;
  ymax = 10.0;
  zmin = -10.0;
  zmax = 10.0;
  robj = 1. - bmarg;
  T = 10.0;

  real primal_T = T;
  real dual_T = 1. * T / 2;

  dolfin_set("output destination","silent");
  if(dolfin::MPI::processNumber() == 0)
    dolfin_set("output destination","terminal");


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

      if(fabs(x[0] - 0.0) < robj &&
         fabs(x[1] - 0.0) < robj &&
         fabs(x[2] - 0.0) < robj)
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

      if(fabs(x[0] - 0.0) < robj &&
         fabs(x[1] - 0.0) < robj &&
         fabs(x[2] - 0.0) < robj)
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


  // Create mesh
  Mesh mesh("mesh.bin");
  //Mesh mesh("mesh.xml");

  Assembler assembler(mesh);

  message("Running on %d %s", dolfin::MPI::numProcesses(), 
          (dolfin::MPI::numProcesses() > 1 ? "nodes" : "node"));
  message("Global number of vertices: %d", 
          (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numVertices() : mesh.numVertices()));
  message("Global number of cells: %d", 
          (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numCells() : mesh.numCells()));

//    for(int i = 0; i < 1; i++)
//      mesh.refine();

  message("Global number of vertices: %d", 
          (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numVertices() : mesh.numVertices()));
  message("Global number of cells: %d", 
          (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numCells() : mesh.numCells()));

  MeshSize h(mesh);
  FacetNormal n(mesh);
  NodeNormal nn(mesh);
  CellVolume cv(mesh);

  // Create boundary conditions
  Inflow inflow(mesh);
  DualInflow dinflow(mesh);
  Outflow outflow(mesh);
  ThetaDrag thetadrag(mesh);
  ThetaLift thetalift(mesh);
  SlipMarker sm(mesh);
  DirichletBoundary dboundary;
  OutflowBoundary oboundary;
  InflowBoundary iboundary;
  SlipBoundary sboundary;
  AllBoundary aboundary;
  DirichletBC bc_in_m(inflow, mesh, iboundary);
  DirichletBC bc_c(outflow, mesh, oboundary);
  SlipBC slipbc_m(mesh, sboundary, nn);
  DirichletBC bc_dm(dinflow, mesh, aboundary);

  Array<BoundaryCondition*> bcs_m;
  bcs_m.push_back(&bc_in_m);
  bcs_m.push_back(&slipbc_m);

  Array<BoundaryCondition*> bcs_dm;
  bcs_dm.push_back(&bc_dm);


  uint *c_indices = 0;
  uint *indices = 0;

  real hmin = h.min();
  cout << "hmin: " << hmin << endl;
  
  PsiMomentum psim(mesh);
  PsiContinuity psic(mesh);
  BPsiMomentum bpsim(mesh);

  real k = 0.1*hmin;
  real nu = 0.0;
  real c1 = 0.1;
  real c2 = 0.0;
  real c3 = 0.0;

  real rdtol = 5e-2;
  int maxit = 10;

  cout << "c1: " << c1 << endl;
  cout << "c2: " << c2 << endl;

  // Declare all needed functions
  Function u;
  Function u0;
  Function p;
  Function p0;
  Function nuf(mesh, nu);
  Function kf(mesh, k);
  Function k0f(mesh, 0.0);
  Function c1f(mesh, c1);
  Function c2f(mesh, c2);
  Function hminf(mesh, c3*hmin);
  Function umean;
  Function dtu;

  Function X;

  Function rd_u;
  Function rd_p;

  Function up;
  Function pp;
  Function up0;
  Function dtup;

  Function Rm, Rmtot;
  Function Rc, Rctot;
  Function wm, wmtot;
  Function wc, wctot;

  Function ei_m;
  Function ei_c;
  Function eij_m;
  Function eij_c;
  Function eif;

  Function tau_1, tau_2, normal;

  Vector ei;

  // Declare primal and dual forms
  Form *a_m, *L_m, *a_c, *L_c;

  NSEMomentum3DBilinearForm ap_m(u, p, nuf, h, kf, c1f, u0, normal, sm);
  NSEMomentum3DLinearForm Lp_m(u, p, nuf, h, kf, c1f, u0, normal, sm);

  NSEContinuity3DBilinearForm ap_c(h, kf, c1f, hminf);
  NSEContinuity3DLinearForm Lp_c(u, p, h, kf, c1f, u0, p0, hminf);

  //NSEDualMomentum3DBilinearForm ad_m(u, up, nuf, h, kf, c1f, c2f, u0);
  NSEDualMomentum3DBilinearForm ad_m(u, up, p, nuf, h, kf, c1f, c2f, u0);
  NSEDualMomentum3DLinearForm Ld_m(u, up, p, nuf, h, kf, c1f, c2f, u0, psim, bpsim);

  NSEDualContinuity3DBilinearForm ad_c(u, h, kf, c1f, c2f, hminf);
  //NSEDualContinuity3DLinearForm Ld_c(u, up, p, h, kf, c1f, u0, p0, hminf, psic);
  NSEDualContinuity3DLinearForm Ld_c(u, p, h, kf, c1f, c2f, u0, p0, hminf, psic);

  NSEErrRepMomentum3DLinearForm Lrep_m(up, pp, nuf, dtup, u);
  NSEErrRepContinuity3DLinearForm Lrep_c(up, p);

  Drag3DFunctional Md(normal, thetadrag, p);
  Drag3DFunctional Ml(normal, thetalift, p);

  NSEH1Functional MH1(u, p, h);
  NSEH12Functional MH12(u, p);
  NSEErrEstFunctional MHerrest(h, cv, Rm, Rc, wm, wc);
  NSEErrEstGlobalFunctional MHerrestg(h, cv, Rm, Rc, wm, wc);

  NSEMomentumResidual3DLinearForm LRm(u, p, kf, u0);
  NSEMomentumResidual3DFunctional MRm(h, cv, Rm);
  NSEContinuityResidual3DLinearForm LRc(u, u0);
  NSEContinuityResidual3DFunctional MRc(h, cv, Rc);
  NSEMomentumResidualGlobal3DFunctional MRgm(u, p, h, kf, u0);
  NSEContinuityResidualGlobal3DFunctional MRgc(u, h, u0);

  NSEH1Momentum3DLinearForm Lwm(u, kf, u0);
  NSEH1Momentum3DFunctional Mwm(h, cv, wm);
  NSEH1Continuity3DLinearForm Lwc(p);
  NSEH1Continuity3DFunctional Mwc(h, cv, wc);

  NSEH1MomentumGlobal3DFunctional Mgwm(u, h, kf, u0);
  NSEH1ContinuityGlobal3DFunctional Mgwc(p, h);

  // Initialize functions
  u.init(mesh, ap_m, 0);
  u0.init(mesh, ap_m, 0);
  p.init(mesh, ap_c, 0);
  p0.init(mesh, ap_c, 0);

  X.init(mesh, ap_m, 0);

  rd_u.init(mesh, ap_m, 0);
  rd_p.init(mesh, ap_c, 0);

  up.init(mesh, ap_m, 0);
  up0.init(mesh, ap_m, 0);
  pp.init(mesh, ap_c, 0);
  dtup.init(mesh, ap_m, 0);

  dtu.init(mesh, ap_m, 0);

  Rm.init(mesh, LRm, 0);
  Rc.init(mesh, LRc, 0);
  wm.init(mesh, Lwm, 0);
  wc.init(mesh, Lwc, 0);

  Rmtot.init(mesh, LRm, 0);
  Rctot.init(mesh, LRc, 0);
  wmtot.init(mesh, Lwm, 0);
  wctot.init(mesh, Lwc, 0);

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

  LinearPDE *pde_m, *pde_c;

  LinearPDE pdep_m(ap_m, Lp_m, mesh, bcs_m);
  LinearPDE pdep_c(ap_c, Lp_c, mesh, bc_c, cg);

  LinearPDE pded_m(ad_m, Ld_m, mesh, bcs_dm);
  LinearPDE pded_c(ad_c, Ld_c, mesh, bc_c);

  Function U;
  Function P;

  computeX(X, &ap_m, mesh);
  Rotate(X, &ap_m, mesh, theta);

  File file_u("velocity.bin");
  File file_p("pressure.bin");
  File file_du("dvelocity.bin");
  File file_dp("dpressure.bin");
  File file_m("mesh.bin");

  //file_m << mesh;

  File meshfile("mesh_out.bin");
  meshfile << mesh;

  int iteration0 = 0;
  int stabcounter = 0;

  real t = 0; real s = 0;


  int no_samples = 200;

  SpaceTimeFunction* Up = 0;
  SpaceTimeFunction* dtUp = 0;
  SpaceTimeFunction* Pp = 0;
  SpaceTimeFunction* Rmp = 0;
  SpaceTimeFunction* Rcp = 0;

  Up = new SpaceTimeFunction(mesh, up);
  Up->setBasename("velocity_v");

  dtUp = new SpaceTimeFunction(mesh, dtup);
  dtUp->setBasename("dtvelocity_v");

  Pp = new SpaceTimeFunction(mesh, pp);
  Pp->setBasename("pressure_v");

  Rmp = new SpaceTimeFunction(mesh, Rm);
  Rmp->setBasename("Rm_v");

  Rcp = new SpaceTimeFunction(mesh, Rc);
  Rcp->setBasename("Rc_v");

  std::string solver = "primal";

  bool coeffchanged = true;

  real int_errest_gstcs = 0;

  for(int solver_idx = 0; solver_idx < 2; solver_idx++)
  {
    if(solver_idx == 1)
      solver = "dual";
    
    if(solver == "dual")
    {
      c1 = 4.0;
      c1f.init(mesh, c1);
    }

    if(solver == "primal")
    {
      cout << "Starting primal solver" << endl;
      a_m = &ap_m; L_m = &Lp_m; a_c = &ap_c; L_c = &Lp_c;
      pde_m = &pdep_m; pde_c = &pdep_c;
    }
    else
    {
      cout << "Starting dual solver" << endl;
      a_m = &ad_m; L_m = &Ld_m; a_c = &ad_c; L_c = &Ld_c;
      pde_m = &pded_m; pde_c = &pded_c;
      T = dual_T;
    }

    u.vector() = 0.0;
    u0.vector() = 0.0;
    p.vector() = 0.0;
    p0.vector() = 0.0;
    up0.vector() = 0.0;

    Rmtot.vector() = 0.0; Rctot.vector() = 0.0; wmtot.vector() = 0.0; wctot.vector() = 0.0;
    
    int stepcounter = 0;
    int sample = 0;
    t = 0;

    real tot_drag = 0;
    real tot_lift = 0;
    int n_mean = 0;

    real tot_H1dualm = 0;
    real tot_H1dualc = 0;
    real tot_H1dualgm = 0;
    real tot_H1dualgc = 0;
    real tot_H1dualgstm = 0;
    real tot_H1dualgstc = 0;
    real tot_H1primal = 0;
    real tot_H1primal2 = 0;

    real tot_Rm = 0;
    real tot_Rc = 0;
    real tot_Rgm = 0;
    real tot_Rgc = 0;
    real tot_Rgstm = 0;
    real tot_Rgstc = 0;

    real int_errest_cs = 0;
    real int_errest_gcs = 0;
    
    k = 0.1*hmin;
    kf.init(mesh, k);

    // Time-stepping
    while(t <= T)
    {
      real stimer = time();

      s = primal_T - t;
      
      if(solver == "dual")
      {
        cout << "eval dual" << endl;
        Up->eval(s);
        dtUp->eval(s);
        Pp->eval(s);
        Rmp->eval(s);
        Rcp->eval(s);
        cout << "eval dual done" << endl;
      }
      
      real umax = u.vector().norm(linf);
      if(solver == "dual")
        umax = up.vector().norm(linf);

      if(stepcounter >= 100)
      {
        if(iteration0 > 5)
          stabcounter = 10;
        k = 4.0*hmin/std::max(1., umax);
        if(hmin >= 0.1)
        {
          k = 0.5*hmin/std::max(1., umax);
        }
        if(stabcounter > 0)
          k /= 4.;
        kf.init(mesh, k);
      }

      // Fixed-point iteration
      for(int i = 0; i < maxit; i++)
      {
        cout << "Solving momentum" << endl;
        real timer = time();
        pde_m->solve(U);
        rd_u.vector() = U.vector();
        rd_u.vector() -= u.vector();
        real rd_u_norm = rd_u.vector().norm(l2);
        u.vector() = U.vector();
        
        cout << "Solving continuity" << endl;
        pde_c->solve(P);
        p.vector() = P.vector();
        rd_p.vector() = p.vector();
        rd_p.vector() -= p0.vector();
        real rd_p_norm = rd_p.vector().norm(l2);
        
        cout << "Iteration info: " <<
          "Unorm: " << U.vector().norm(linf) <<
          " Pnorm: " << P.vector().norm(linf) <<
          " Uincr: " <<  rd_u_norm <<
          " Pincr: " <<  rd_p_norm <<
          " k: " << k <<
          " step: " << stepcounter <<
          " t: " << t <<
          " timer: " << time() - timer << endl;
        cout << "iteration: " << i << endl;
        iteration0 = i;
        if(rd_u_norm / u.vector().norm(l2) <= rdtol && rd_p_norm / p.vector().norm(l2) <= rdtol)
        {
          cout << "Step info: " <<
            "Unorm: " << U.vector().norm(linf) <<
            " Pnorm: " << P.vector().norm(linf) <<
            " Uincr: " <<  rd_u_norm / u.vector().norm(l2) <<
            " Pincr: " <<  rd_p_norm / p.vector().norm(l2) <<
            " k: " << k <<
            " step: " << stepcounter <<
            " iters: " << iteration0 + 1 <<
            " t: " << t <<
            " timer: " << time() - stimer << endl;
          break;
        }
        p0.vector() = p.vector();
      }
      
      if(solver == "dual")
      {
        cout << "errest" << endl;
        assembler.assemble(eij_m.vector(), Lrep_m);
        ei_m.vector().axpy(k, eij_m.vector());
        assembler.assemble(eij_c.vector(), Lrep_c);
        ei_c.vector().axpy(k, eij_c.vector());
        cout << "errest done: " << ei_m.vector().norm(linf) << " " << ei_c.vector().norm(linf) << endl;

        assembler.assemble(wm.vector(), Lwm);
        real H1dualm = assembler.assemble(Mwm);
        tot_H1dualm += k*H1dualm;
        real H1dualgm = sqrt(assembler.assemble(Mgwm));
        tot_H1dualgm += k*H1dualgm;
        real H1dualgstm = assembler.assemble(Mgwm);
        tot_H1dualgstm += k*H1dualgstm;
        wmtot.vector().axpy(k, wm.vector());
        assembler.assemble(wc.vector(), Lwc);
        real H1dualc = assembler.assemble(Mwc);
        tot_H1dualc += k*H1dualc;
        real H1dualgc = sqrt(assembler.assemble(Mgwc));
        tot_H1dualgc += k*H1dualgc;
        real H1dualgstc = assembler.assemble(Mgwc);
        tot_H1dualgstc += k*H1dualgstc;
        wctot.vector().axpy(k, wc.vector());

        real Rmi = assembler.assemble(MRm);
        tot_Rm += k*Rmi;
        real Rgmi = 0.0;
        Rmtot.vector().axpy(k, Rm.vector());
        real Rci = assembler.assemble(MRc);
        tot_Rc += k*Rci;
        real Rgci = 0.0;
        Rctot.vector().axpy(k, Rc.vector());

        real errest_cs = assembler.assemble(MHerrest);
        int_errest_cs += k*errest_cs;
        real errest_gcs = sqrt(assembler.assemble(MHerrestg));
        int_errest_gcs += k*errest_gcs;
        real errest_gstcs = assembler.assemble(MHerrestg);
        int_errest_gstcs += k*errest_gstcs;
        n_mean++;

        cout << "step dual t: " << t <<
          " dualm: " << H1dualm <<
          " dualc: " << H1dualc <<
          " dualgm: " << H1dualgm <<
          " dualgc: " << H1dualgc <<
          " Rm: " << Rmi <<
          " Rc: " << Rci <<
          " Rgm: " << Rgmi <<
          " Rgc: " << Rgci <<
          " errest_cs: " << errest_cs <<
          " errest_gcs: " << errest_gcs <<
          endl;
      }

      if(solver == "primal")
      {
        real drag = 0.0, lift = 0.0;
        drag = assembler.assemble(Md);
        cout << "drag: " << drag << " t = " << t << endl;
        lift = assembler.assemble(Ml);
        cout << "lift: " << lift << " t = " << t << endl;

        assembler.assemble(Rm.vector(), LRm);
        real Rmi = assembler.assemble(MRm);
        cout << "step primal Rm: " << Rmi << endl;
        assembler.assemble(Rc.vector(), LRc);
        real Rci = assembler.assemble(MRc);
        cout << "step primal Rc: " << Rci << endl;

        if(t >= dual_T)
        {
          // Output drag and lift, together with other diagnostics

          tot_drag = (drag + n_mean*tot_drag) / (n_mean + 1);
          cout << "step t: " << t << " drag: " << drag << " lift: " << lift << endl;
          real H1primal = assembler.assemble(MH1);
          tot_H1primal = (H1primal + n_mean*tot_H1primal) / (n_mean + 1);
          real H1primal2 = assembler.assemble(MH12);
          tot_H1primal2 = (H1primal2 + n_mean*tot_H1primal2) / (n_mean + 1);
          cout << "step H1 primal: " << tot_H1primal << endl;
          cout << "step H1 primal2: " << tot_H1primal2 << endl;
          n_mean++;

          real Rgstmi = assembler.assemble(MRgm);
          tot_Rgstm += k*Rgstmi;
          real Rgstci = assembler.assemble(MRgc);
          tot_Rgstc += k*Rgstci;
        }
      }

      if(stepcounter == 0 || t > T*(real(sample)/real(no_samples)))
      {
        if(solver == "primal")
        {
          file_u << u; file_p << p;

          // Record primal solution
          std::stringstream number;
          number << std::setfill('0') << std::setw(6) << sample;

          // Save primal velocity
          up.vector() = u.vector(); up.vector() += u0.vector(); up.vector() /= 2.;
          File ubinfile(Up->getNewFilename(t));
          ubinfile << up.vector();

          // Save primal velocity time-derivative
          dtu.vector() = u.vector(); dtu.vector() -= u0.vector(); dtu.vector() /= k;
          File dtubinfile(dtUp->getNewFilename(t));
          dtubinfile << dtu.vector();

          // Save primal pressure
          File pbinfile(Pp->getNewFilename(t));
          pbinfile << p.vector();

          // Save primal residuals
          File Rmbinfile(Rmp->getNewFilename(t));
          Rmbinfile << Rm.vector();

          File Rcbinfile(Rcp->getNewFilename(t));
          Rcbinfile << Rc.vector();
        }
        else
        {
          file_du << u; file_dp << p;
        }
        
        sample++;
      }
      
      u0.vector() = u.vector();
      up0.vector() = up.vector();

      t += k;
      stepcounter++;
      if(stabcounter > 0)
        stabcounter--;
    }
    

    cout << "Solver done" << endl;

    if(solver == "primal")
    {
      cout << "mean drag: " << tot_drag << endl;
      cout << "total H1primal: " << sqrt(tot_H1primal) << endl;
      cout << "total H1primal2: " << sqrt(tot_H1primal2) << endl;
      cout << "total Rgstm: " << sqrt(tot_Rgstm) << endl;
      cout << "total Rgstc: " << sqrt(tot_Rgstc) << endl;

      int_errest_gstcs = (sqrt(tot_Rgstm) + sqrt(tot_Rgstc));
    }

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

      
      eif.vector() = 0.0;
      eif.vector() += ei_m.vector(); eif.vector() += ei_c.vector();
      //eif.vector() /= dual_T;
      
      NSEErrRepMomentum3DFunctional M_ei(eif, cv);
      real errest = fabs(assembler.assemble(M_ei));

      //int_errest_cs = sqrt(int_errest_cs);

      int_errest_gstcs *= (sqrt(tot_H1dualgstm) + sqrt(tot_H1dualgstc));

      cout << "error estimate: " << errest / dual_T << endl;
      cout << "error estimate cs: " << int_errest_cs / dual_T << endl;
      cout << "error estimate gcs: " << int_errest_gcs / dual_T << endl;
      cout << "error estimate gstcs: " << int_errest_gstcs / dual_T << endl;
      cout << "total H1dualm: " << tot_H1dualm << endl;
      cout << "total H1dualc: " << tot_H1dualc << endl;
      cout << "total H1dualgm: " << tot_H1dualgm << endl;
      cout << "total H1dualgc: " << tot_H1dualgc << endl;
      cout << "total H1dualgstm: " << sqrt(tot_H1dualgstm) << endl;
      cout << "total H1dualgstc: " << sqrt(tot_H1dualgstc) << endl;
      cout << "total Rm: " << tot_Rm << endl;
      cout << "total Rc: " << tot_Rc << endl;
      cout << "total Rgm: " << tot_Rgm << endl;
      cout << "total Rgc: " << tot_Rgc << endl;

      File file_ei("ei.bin");
      file_ei << eif;

      File file_Rmtot("Rmtot.bin");
      file_Rmtot << Rmtot;
      File file_Rctot("Rctot.bin");
      file_Rctot << Rctot;
      File file_wmtot("wmtot.bin");
      file_wmtot << wmtot;
      File file_wctot("wctot.bin");
      file_wctot << wctot;

      MeshFunction<real> eimf;
      eimf.init(mesh, mesh.topology().dim());
  
      // Initialize eimf - assumption on dofmap for DG0
      for (CellIterator c(mesh); !c.end(); ++c)
      {
        eimf.set(*c, eif.vector()[c->index()]);
      }

      MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
  
      cout << "Output eimf: " << endl;
      File file_eimf("eimf.bin");
      file_eimf << eimf;
      

      ComputeRefinementMarkers(mesh, adapt_percent, ei, cell_marker);

      if(MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");
      message("Adaptive refinement");
      message("cells before: %d",
              (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numCells() : mesh.numCells()));
      message("vertices before: %d",
              (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numVertices() : mesh.numVertices()));
      dolfin_set("output destination","silent");

      RivaraRefinement::refine(mesh, cell_marker);

      if(MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");
      message("cells after: %d",
              (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numCells() : mesh.numCells()));
      message("vertices after: %d",
              (dolfin::MPI::numProcesses() > 1 ? mesh.distdata().global_numVertices() : mesh.numVertices()));
      dolfin_set("output destination","silent");
      
      File file_rm("rmesh.bin");
      file_rm << mesh;
    }
  }
  

  return 0;
}
