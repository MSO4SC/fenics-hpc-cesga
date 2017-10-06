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
#include "ufc2/compressible2d.h"
#include "ufc2/compressible2d_p.h"
#include "ufc2/compressible2d_mach.h"
#include "ufc2/compressible2d_velo.h"
#include "ufc2/compressible2d_sound.h"
#include "ufc2/L2proj2D.h"
#include "ufc2/L2proj2D_densi.h"
#include "ufc2/L2proj2D_energ.h"
#include "ufc2/L2proj2D_momen.h"
#include "ufc2/DirichletBC_CG1v_CG1.h"

#include <dolfin/common/constants.h>
#include <dolfin/fem/Assembler.h>
#include <dolfin/fem/DirichletBC.h>

#include <dolfin/function/Function.h>
#include <dolfin/parameter/parameters.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>
#include <dolfin/la/KrylovSolver.h>
#include <dolfin/fem/UFC.h>

using namespace dolfin;

real XMIN = 0.; 
real XMAX = 4.; 
real YMIN = 0.; 
real YMAX = 1.; 

real eps = 1e-5;
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
  class CellDiameter : public Function
  {
  public:
    
    CellDiameter(Mesh& mesh) : Function(mesh) {
	}
    
    void eval(real * value, const real* x) const
    {
     	value[0] = cell().diameter();

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


  class Phi_init : public Function
  {
  public:
    
    Phi_init(Mesh& mesh) : Function(mesh) {
	}
    
    void eval(real * value, const real* x) const
    {
     	value[0] =  1.0;
     	value[1] = -1.0; 
     	value[2] =  1.0; 
     	value[3] =  2.0; 

    }

    uint rank() const
    {
      return 1;
    }

    uint dim(uint i) const
    {
      return 4;
    }

  };

  // Sub domain for Dirichlet boundary condition
  class IN_boundary : public SubDomain
  {
    bool inside(const real* x, bool on_boundary) const
    {
      return on_boundary && ((fabs(x[0] - XMIN) < DOLFIN_EPS) || (fabs(x[1] - YMAX) < DOLFIN_EPS)  ); 
    }
  };

  class OUT_boundary: public SubDomain
  {
    bool inside(const real* x, bool on_boundary) const
    {
      return on_boundary && ((fabs(x[0] - XMAX) < DOLFIN_EPS)); 
    }
  };
 
  class WALL_boundary: public SubDomain
  {
    bool inside(const real* x, bool on_boundary) const
    {
      return on_boundary && ((fabs(x[1] - YMIN) < DOLFIN_EPS)); 
    }
  };

  // Create mesh
  //UnitSquare mesh(32, 32);a

  real time_v=time();

  Mesh mesh(argv[1]);
  
   
  message("Input file reading (%s) took : %f ", argv[1] ,time()- time_v);
  message("Total number of vertices : %d", mesh.distdata().global_numVertices());
  File fout("output_mesh.bin"); 
  fout << mesh;
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

  Assembler ass(mesh);
  KrylovSolver gmres_s(bicgstab); 

 // Line 50 in compressible-dim4.py
  Function Phi_0(mesh); 
  Vector Phi_0x; 
  Phi_init phi_init(mesh);  

  Matrix A_phi; 
  Vector b_phi; 
  

  Form *L2WWBilinear = new L2proj2DBilinearForm();
  Form *L2WWLinear =   new L2proj2DLinearForm(phi_init);

  Phi_0.init(mesh, Phi_0x, *L2WWBilinear, 0); 

  ass.assemble(A_phi, *L2WWBilinear, true);
  ass.assemble(b_phi, *L2WWLinear, true);
 
 

  gmres_s.solve(A_phi, Phi_0x, b_phi);

 // ------------------------------
 // BC's 
 
 Function u_0_in(mesh, 1.0);
 Function u_1_in(mesh, -1.0);
 Function E_in  (mesh, 2.0);
 Function rho_out(mesh, 1.0);
 Function u_1_wall(mesh, 0.0);

 IN_boundary in_boundary;
 OUT_boundary out_boundary; 
 WALL_boundary wall_boundary; 

 DirichletBC_CG1v_CG1 bc_0(mesh, in_boundary , 0, &u_0_in); 
 DirichletBC_CG1v_CG1 bc_1(mesh, in_boundary , 1, &u_1_in); 
 DirichletBC_CG1v_CG1 bc_2(mesh, in_boundary , 3, &E_in);
 DirichletBC_CG1v_CG1 bc_3(mesh, out_boundary, 2, &rho_out);
 DirichletBC_CG1v_CG1 bc_4(mesh, wall_boundary    , 1, &u_1_wall);
  


// ----------------------------------
// Files
//

File file_momen("momen.bin"); 
File file_densi("densi.bin"); 
File file_energ("energ.bin");
File file_mach ("mach.bin");
File file_veloc("veloc.bin");
File file_sound("sound.bin");



// ----------------------------------
// Vectors, Functions
CellDiameter h_el(mesh); 
Function R_gases(mesh, 2.0);
Function cv(mesh, 2.0); 
real dt = 5e-4;
Function delta_t(mesh, dt);  


Function momen, energ, densi, velo, mach, sound; 
Vector   momenX, energX, densiX, veloX, machX, soundX; 

Function pr;
Vector prX; 

// Forms
// ----------------------------------

  Form *c2DpBilinear = new compressible2d_pBilinearForm();
  Form *c2DpLinear =   new compressible2d_pLinearForm(Phi_0, R_gases, cv );
 

  Form *L2momenBilinear = new L2proj2D_momenBilinearForm();
  Form *L2momenLinear =   new L2proj2D_momenLinearForm(Phi_0);
  Form *L2energBilinear = new L2proj2D_energBilinearForm();
  Form *L2energLinear =   new L2proj2D_energLinearForm(Phi_0);
  Form *L2densiBilinear = new L2proj2D_densiBilinearForm();
  Form *L2densiLinear =   new L2proj2D_densiLinearForm(Phi_0);

  Form *c2DBilinear = new compressible2dBilinearForm();
  Form *c2DLinear =   new compressible2dLinearForm(Phi_0, h_el, R_gases, cv, delta_t, pr);
  Form *c2DveloBilinear = new compressible2d_veloBilinearForm();
  Form *c2DveloLinear =   new compressible2d_veloLinearForm(Phi_0);
  Form *c2DmachBilinear = new compressible2d_machBilinearForm();
  Form *c2DmachLinear =   new compressible2d_machLinearForm(Phi_0, R_gases, cv, pr);
  Form *c2DsoundBilinear = new compressible2d_soundBilinearForm();
  Form *c2DsoundLinear =   new compressible2d_soundLinearForm(Phi_0, R_gases, cv, pr);



// Matrices, vectors 
// ----------------------------------
 
  Matrix A_p; 
  Vector b_p; 

  Matrix A_L2momen; 
  Vector b_L2momen; 

  Matrix A_L2energ; 
  Vector b_L2energ; 

  Matrix A_L2densi; 
  Vector b_L2densi; 

  Matrix A; 
  Vector b; 

  Matrix A_velo;
  Vector b_velo; 

  Matrix A_mach;
  Vector b_mach; 

  Matrix A_sound; 
  Vector b_sound; 
  
// inits 
// ----------------------------------
  pr.init(mesh, prX, *c2DpBilinear,0);

  momen.init(mesh, momenX, *L2momenBilinear, 0); 
  energ.init(mesh, energX, *L2energBilinear, 0); 
  densi.init(mesh, densiX, *L2densiBilinear, 0); 

  velo.init( mesh, veloX , *c2DveloBilinear, 0); 
  mach.init( mesh, machX , *c2DmachBilinear, 0); 
  sound.init(mesh, soundX, *c2DsoundBilinear, 0); 


  ass.assemble(A_L2momen, *L2momenBilinear, true);
  ass.assemble(b_L2momen, *L2momenLinear, true);

  ass.assemble(A_L2energ, *L2energBilinear, true);
  ass.assemble(b_L2energ, *L2energLinear, true);
  
  ass.assemble(A_L2densi, *L2densiBilinear, true);
  ass.assemble(b_L2densi, *L2densiLinear, true);
  


  ass.assemble(A_p, *c2DpBilinear, true);
  ass.assemble(b_p, *c2DpLinear, true);

  ass.assemble(A, *c2DBilinear, true);
  ass.assemble(b, *c2DLinear, true);


  ass.assemble(A_velo, *c2DveloBilinear, true);
  ass.assemble(b_velo, *c2DveloBilinear, true);
  
  ass.assemble(A_mach, *c2DmachBilinear , true);
  ass.assemble(b_mach, *c2DmachBilinear , true);
  
  ass.assemble(A_sound, *c2DsoundBilinear , true);
  ass.assemble(b_sound, *c2DsoundBilinear , true);
  


// time stepping loop 
//-------------------------------------
//

real t = 0; 
for (int stepcounter=0; stepcounter< 300; stepcounter++)
{ 
	t+=dt; 
	warning("solving main equation for t= %g", t); 

  	ass.assemble(b_p, *c2DpLinear, false);
 	gmres_s.solve(A_p, prX, b_p);

  	// truncate pr 
  	{

    	int d = mesh.topology().dim();
    
    		{
      		UFC ufc(c2DpBilinear->form(), mesh, c2DpBilinear->dofMaps());
      		Cell c(mesh, 0);
      		uint local_dim = c.numEntities(0);
      		uint *idx  = new uint[local_dim];
      		uint *id  = new uint[local_dim];
      		real *rho_block = new real[local_dim];  

		real rho_val;

      		for (CellIterator cell(mesh); !cell.end(); ++cell)
      		{
			ufc.update(*cell, mesh.distdata());
			(c2DpBilinear->dofMaps())[0].tabulate_dofs(idx, ufc.cell, cell->index());
	
			pr.vector().get(rho_block, local_dim, idx);
		
			uint ii = 0;
			uint jj = 0;    

			for(VertexIterator v(*cell); !v.end(); ++v, ii++) 
        	 	{    
			        if (!mesh.distdata().is_ghost(v->index(), 0)) 
          			{    
            			 	if (rho_block[jj] < 0.0) 
						rho_block[jj]=1e-6;  			
            				id[jj++] = idx[ii];
          			}    
        		}  

			pr.vector().set(rho_block, jj, id);

      		}
      
        	pr.vector().apply();
      		delete[] rho_block;
      		delete[] idx;
      		delete[] id;
    		}
  	}



  	ass.assemble(b_velo, *c2DveloLinear, false);
  	ass.assemble(b_mach, *c2DmachLinear, false);
  	ass.assemble(b_sound, *c2DsoundLinear, false);

 	gmres_s.solve(A_velo, veloX, b_velo);
 	gmres_s.solve(A_mach, machX, b_mach);
 	gmres_s.solve(A_sound, soundX, b_sound);

  	ass.assemble(b, *c2DLinear, false);
        bc_0.apply(A, b, *c2DBilinear); 	
        bc_1.apply(A, b, *c2DBilinear); 	
        bc_2.apply(A, b, *c2DBilinear); 	
        bc_3.apply(A, b, *c2DBilinear); 	
        bc_4.apply(A, b, *c2DBilinear); 	


 	gmres_s.solve(A, Phi_0x, b);


  	ass.assemble(b_L2momen, *L2momenLinear, false);
  	ass.assemble(b_L2energ, *L2energLinear, false);
 	ass.assemble(b_L2densi, *L2densiLinear, false);

 	gmres_s.solve(A_L2momen, momenX, b_L2momen);
 	gmres_s.solve(A_L2energ, energX, b_L2energ);
 	gmres_s.solve(A_L2densi, densiX, b_L2densi);

	if (!(stepcounter %10)) {
	warning("saving at main equation for t= %g", t); 
	file_momen << momen;
	file_densi << densi; 
	file_energ << energ; 
	file_mach  << mach; 
	file_veloc << velo; 
	file_sound << sound; 
	}

	

}
// Cleanup 
// ----------------------------------
 delete L2WWLinear; 
 delete L2WWBilinear; 

 delete L2momenBilinear;
 delete L2momenLinear  ;
 delete L2energBilinear; 
 delete L2energLinear  ;
 delete L2densiBilinear; 
 delete L2densiLinear  ;

  delete c2DBilinear; 
  delete c2DLinear;
  delete c2DveloBilinear;
  delete c2DveloLinear;
  delete c2DmachBilinear;
  delete c2DmachLinear;
  delete c2DsoundBilinear;
  delete c2DsoundLinear;

}
dolfin_finalize();
return 0;
}
