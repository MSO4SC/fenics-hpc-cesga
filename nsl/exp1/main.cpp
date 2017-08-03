// Laura Saavedra, N. Cem Degirmenci 2017
#include <dolfin/main/init.h>
#include <dolfin/config/dolfin_config.h>
#include "ufc2/NavierStokes2D.h"
#include "ufc2/L2ProjPfromM2D.h"
#include "ufc2/L2ProjUfromM2D.h"
#include "ufc2/L2error.h"

#include <dolfin/common/constants.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/DirichletBC.h>

#include <dolfin/function/Function.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/la/KrylovSolver.h>
#include <dolfin/mesh/RivaraRefinement.h>


using namespace dolfin;

int main(int argc, char **argv)
{
 dolfin_init(argc, argv);

 dolfin_set("output destination","silent");
 if(dolfin::MPI::processNumber() == 0)
     dolfin_set("output destination","terminal");

 {
 class UIN : public Function
  {
  public:
    
    UIN(Mesh& mesh) : Function(mesh) {
        geom_dim = mesh.geometry().dim();
	std::cout << " geometrical dimension :" << geom_dim << std::endl ;
        }   
    
    void eval(real * value, const real* x) const
    {   
     if (geom_dim == 2) {
        value[0] = x[1]*x[1];
        value[1] = x[0]*x[0];
     }
     if (geom_dim == 3){
        value[0] = 3*pow(DOLFIN_PI, 2)*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1])*sin(DOLFIN_PI*x[2]);
      }
    }   

    uint rank() const
    {   
      return 1;
    }   

    uint dim(uint i) const
    {   
      return 2;
  }
    uint geom_dim;

  };  

 class F : public Function
  {
  public:
    
    F(Mesh& mesh) : Function(mesh) {
        geom_dim = mesh.geometry().dim();
        }   
    
    void eval(real * value, const real* x) const
    {   
     if (geom_dim == 2) {
        value[0] = 2*x[1]*x[0]*x[0] + 2;
        value[1] = 2*x[0]*x[1]*x[1] + 2;
     }
     if (geom_dim == 3){
        value[0] = 3*pow(DOLFIN_PI, 2)*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1])*sin(DOLFIN_PI*x[2]);
      }
    }   

    uint rank() const
    {   
      return 1;
    }   

    uint dim(uint i) const
    {   
      return 2;
  }
    uint geom_dim;

  };  


 
 class D : public Function
  {
  public:
    
    D(Mesh& mesh) : Function(mesh) {
        geom_dim = mesh.geometry().dim();
        }   
    
    void eval(real * value, const real* x) const
    {   
     if (geom_dim == 2) {
        value[0] = .2*pow(cell().diameter(), 1.5);
     }
     if (geom_dim == 3){
        value[0] = 3*pow(DOLFIN_PI, 2)*sin(DOLFIN_PI*x[0])*sin(DOLFIN_PI*x[1])*sin(DOLFIN_PI*x[2]);
      }
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


 
 class Gamma : public Function
  {
  public:
    
    Gamma(Mesh& mesh) : Function(mesh) {
        geom_dim = mesh.geometry().dim();
        }   
    
    void eval(real * value, const real* x) const
    {   
     if (geom_dim == 2) {
        value[0] = 1e3/(cell().diameter());
     }
     if (geom_dim == 3){
        value[0] =  1e3/(cell().diameter());
      }
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


 //Mesh mesh("mesh2D.xml");
 message("reading the mesh :%s ", argv[1]);
 Mesh mesh(argv[1]);



real  time_v = time();
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






/*

  MeshFunction<bool> cell_refinement_marker(mesh);
  cell_refinement_marker.init(mesh.topology().dim());
  cell_refinement_marker = true;
  RivaraRefinement::refine(mesh, cell_refinement_marker);

*/
  real dt = 0.002; 
  real t = 0;  
  real T = 1;  
  real nu = 1;  
  real gamma = 1e5;
  int max_nonlinear = 20; 
  real nonlinear_TOL = 1e-8;


  UIN uin(mesh);
  F   f(mesh);
  D   d(mesh);
  Function fk(mesh, dt);
  Function fnu(mesh, nu);
  Function fgamma(mesh, gamma);
 // Gamma fgamma(mesh); 

  Function up, u0, p, up_prev; 
  Vector   upx, u0x, px, up_prevx;

  Function soln; 
  Vector   solnx;


  Form *a, *L, *aUP, *LUP, *L2error;
  a = new NavierStokes2DBilinearForm( fk, fnu, d, fgamma, u0, soln);
  L = new NavierStokes2DLinearForm(uin, f, fk, fnu, d, fgamma, u0, soln );

  aUP = new L2ProjUfromM2DBilinearForm();
  LUP = new L2ProjUfromM2DLinearForm(soln);

  L2error = new L2errorFunctional(uin, u0); 

  up.init(mesh, upx, *L, 1); 
  
  up_prev.init(mesh, up_prevx, *L, 1); 

  u0.init(mesh, u0x, *L, 1); 
  soln.init(mesh, solnx, *L, 0); 


  u0x= 0;
  upx= 0;
  
  KrylovSolver solver(bicgstab);
  Assembler    ass(mesh);
 Matrix A; 
  Vector b;
 
  ass.assemble(A, *a, true);
  ass.assemble(b, *L, true);
 


  KrylovSolver solverU(gmres);
  Assembler    assU(mesh);
 
  Matrix AU;
  Vector bU;
 
  ass.assemble(AU, *aUP, true);
  ass.assemble(bU, *LUP, true);
 

 time_v = time();

 for (t = 0; t < T; t+=dt)
 {
        for (int i = 0; i < max_nonlinear; i++) {
        //a = new NavierStokes2DBilinearForm(fk, fnu, d, fgamma, up);
        //L = new NavierStokes2DLinearForm(uin, f, fk, fnu, d, fgamma, u0, up );


        up_prevx = upx;

        ass.assemble(A, *a, false);
        ass.assemble(b, *L, false);
        solver.solve(A, solnx, b);

        ass.assemble(bU, *LUP, false);
        solver.solve(AU, upx, bU);

        up_prevx -= upx;
        message("xxxx time step :%f, iteration :%d, diff.l2 :%g",  t, i,  up_prevx.norm());
        if (up_prevx.norm() < nonlinear_TOL)
                break;

        }
        u0.vector() = up.vector();

	real l2e = ass.assemble(*L2error);
        message("xxxx time step :%f, l2 error :%g",  t, sqrt(l2e));
 }

message("main loop took : %f" ,time()- time_v);

 }

dolfin_finalize();



}
