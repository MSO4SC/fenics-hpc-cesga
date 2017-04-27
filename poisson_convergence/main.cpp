// Copyright (C) 2006-2007 Anders Logg.
// 2017 Niyazi Cem Degirmenci.
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2006-02-07
// Last changed: 2017-04-27
//
// This demo program solves Poisson's equation
//
//     - div grad u(x, y) = f(x, y)
//
// on the unit square with source f given by
//
//     f(x, y) = sin(pi*x)*sin(pi*y) in 2D
//     f(x, y, z) = sin(pi*x)*sin(pi*y)*sin(pi*z) in 3D
//
// and boundary conditions given by
//
//     u(x, y)     = 0               for x on boundary

#include <dolfin/main/init.h>
#include <dolfin/config/dolfin_config.h>
#include "ufc2/Poisson2D.h"
#include "ufc2/Poisson3D.h"

#include <dolfin/common/constants.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/DirichletBC.h>

#include <dolfin/function/Function.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/la/KrylovSolver.h>

using namespace dolfin;

int main(int argc, char **argv)
{
 dolfin_init(argc,argv);

// std::cout << argc << std::endl; 
// std::cout << sin(0.5) << std::endl; 
// std::cout << sin(DOLFIN_PI*0.5) << std::endl; 
dolfin_set("output destination","silent");
if(dolfin::MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");

{
  // Source term
  class Source : public Function
  {
  public:
    
    Source(Mesh& mesh) : Function(mesh) {
	geom_dim = mesh.geometry().dim();
	}
    
    void eval(real * value, const real* x) const
    {
     if (geom_dim == 2)
     	value[0] = 2*pow(DOLFIN_PI, 2)*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1]);

     if (geom_dim == 3)
        value[0] = 3*pow(DOLFIN_PI, 2)*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1])*sin(DOLFIN_PI*x[2]);
    }

    uint rank() const
    {
      return 0;
    }

    uint dim(uint i) const
    {
      return 1;
  }
    uint geom_dim;

  };

  // Neumann boundary condition
  class Flux : public Function
  {
  public:

    Flux(Mesh& mesh) : Function(mesh) {}

    void eval(real * value, const real* x) const
    {
        value[0] = 0.0;
    }

    uint rank() const
    {
      return 0;
    }

    uint dim(uint i) const
    {
      return 1;
    }

  };

  // Sub domain for Dirichlet boundary condition
  class DirichletBoundary : public SubDomain
  {
    bool inside(const real* x, bool on_boundary) const
    {
      return on_boundary; 
    }
  };

  // Create mesh
  //UnitSquare mesh(32, 32);a

  real time_v=time();

  Mesh mesh(argv[1]);
   
  message("Input file reading (%s) took : %f ", argv[1] ,time()- time_v);
  message("Total number of vertices : %d", mesh.distdata().global_numVertices());

  time_v = time();
 // hmin hmax say
 real hmin = 1.;
 real hmax = -1.;
 for (CellIterator c(mesh); !c.end(); ++c)
 {
	if(c->diameter() < hmin)
		hmin = c->diameter();
	if(c->diameter() > hmax)
		hmax = c->diameter();
 }
 real hmin_global, hmax_global;

  MPI_Allreduce(&hmin, &hmin_global, 1, MPI_DOUBLE,
                MPI_MIN, dolfin::MPI::DOLFIN_COMM);  
  MPI_Allreduce(&hmax, &hmax_global, 1, MPI_DOUBLE,
                MPI_MAX, dolfin::MPI::DOLFIN_COMM);  

 message("hmin hmax computation took : %f" ,time()- time_v);

 message("hmin : %f hmax : %f", hmin_global, hmax_global);



  // Create functions
  Source f(mesh);
  Flux g(mesh);

  Form *a, *L;

  // Create boundary condition
  Function u0(mesh, 0.0);
  DirichletBoundary boundary;
  DirichletBC bc(u0, mesh, boundary);
  

  // Define PDE
 
  if (mesh.topology().dim() == 2) {

   a = new Poisson2DBilinearForm();
   L = new Poisson2DLinearForm(f,g);
  
  }
  else if (mesh.topology().dim() == 3) {  

   a = new Poisson3DBilinearForm();
   L = new Poisson3DLinearForm(f,g);

  }
  else {
	error("dolfin-HPC is designed to work only on triangular or tetrahedral meshes");
  }
 
  Matrix A;
  Vector b; 



  Assembler ass(mesh); 


  time_v = time();
  ass.assemble(A, *a, true);
  message("Assembling A took : %f" ,time()- time_v);

  time_v = time();
  ass.assemble(b, *L, true);
  message("Assembling b took : %f" ,time()- time_v);
 

  time_v = time();
  bc.apply(A,b,*a);
  message("Applying BCs took : %f" ,time()- time_v);


  Function u(mesh);
  Vector x;
  u.init(mesh, x,  *a, 1);

  KrylovSolver solver(gmres);
  time_v = time();
  solver.solve(A,x,b);
  message("Solving took : %f" ,time()- time_v);
  x.apply(); 


  real value[1];
  real coords[] = {49e-2, 49e-2, 49e-2};

  u.eval(value, coords);
  real value_global; 
  MPI_Allreduce(&value[0], &value_global, 1, MPI_DOUBLE,
                MPI_MIN, dolfin::MPI::DOLFIN_COMM);  


  
  if (mesh.topology().dim() == 2)
  	message("Error in (%g, %g) : %g", coords[0], coords[1],  fabs(value_global - sin(DOLFIN_PI*coords[0])*sin(DOLFIN_PI*coords[1]) ));
  else
  	message("Error in (%g, %g, %g) : %g", coords[0], coords[1], coords[2], fabs(value_global - (sin(DOLFIN_PI*coords[0])*sin(DOLFIN_PI*coords[1])*sin(DOLFIN_PI*coords[2])    ) ));


   
  time_v = time();
  // Save solution to file
  File file("output.pvd");
  file << u;
  message("Writing pvd file took : %f" ,time()- time_v);



  time_v = time();
  File file2("output.bin");
  file2 << u;
  message("Writing bin file took : %f" ,time()- time_v);

  delete a;
  delete L;
}
dolfin_finalize();
return 0;
}
