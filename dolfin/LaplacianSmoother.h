#ifndef __LAPLACIAN_SMOOTHER_H
#define __LAPLACIAN_SMOOTHER_H

#include <dolfin.h>
#include <cstring>
#include <dolfin/fem/UFC.h>
#include "MeshBC.h"
#include "NodeNormal.h"
#include "SlipBC.h"

#if HAVE_SUNPERF_H
#include <sunperf.h>
#elif HAVE_SCSL_CBLAS_H
#include <cmplrs/cblas.h>
#elif HAVE_GSL_CBLAS_H
extern "C" {
#include <gsl_cblas.h>
}
#elif HAVE_CBLAS_H
extern "C" {
#include <cblas.h>
}
#endif

#define RM(row,col,nrow) ((row) + ((nrow)*(col)))
namespace dolfin 
{
//-----------------------------------------------------------------------------
  class LaplacianSmoother
  {
  public:
    
    // Sub domain for Dirichlet boundary condition
    class DirichletBoundary : public SubDomain
    {
    public:
      bool inside(const real* x, bool on_boundary) const
      {
	if(on_boundary)
	  return true;
	else
	  return false;
      }
    };

    class MySource : public Function
    {
    public:
    MySource(Mesh& mesh) : Function(mesh)
      {
      }
      
  void eval(real* values, const real* x, int& xx) const
	{eval(values,x);}
      void eval(real* values, const real* x) const
      {
	int d = cell().dim();

	for(int i = 0; i < d; i++)
	{
	  values[i] = 0.0;
	}
      }
    };
    
    class MyBC : public Function
    {
    public:
    MyBC(Mesh& mesh) :
      Function(mesh)
      {
      }
      
  void eval(real* values, const real* x, int& xx) const
	{eval(values,x);}
      void eval(real* values, const real* x) const
      {
	int d = mesh().topology().dim();
	
	for(int i = 0; i < d; i++)
	{
	  values[i] = 0.0;
	}
      }
    };

    /// Constructor
    LaplacianSmoother(Mesh& mesh, MeshFunction<bool>& masked_vertices,
		      Vector* node_values, NodeNormal& node_normal);
    ~LaplacianSmoother();

    void smooth(Vector& motionx,
		bool reset, bool with_slip_bc=false);

    Mesh& mesh;
    Matrix* storeA;
    Vector* storeb;
    Vector* node_values;
    //Assembler* assembler;

    MeshFunction<bool>& masked_vertices;
    NodeNormal& node_normal;
    DirichletBoundary dboundary;
    MySource f;
    MyBC bcf;
    DirichletBC bc0;
    SlipBC bcslip;
    MeshBC bc1;
    Form* a;
    Form* L;

    

  };
//-----------------------------------------------------------------------------
}

#endif
