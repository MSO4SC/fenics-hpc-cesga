#include "dolfin/LaplacianSmoother.h"
#include "ufc2/Laplacian2D.h"
#include "ufc2/Laplacian3D.h"

using namespace dolfin;

//-----------------------------------------------------------------------------
LaplacianSmoother::LaplacianSmoother(Mesh& mesh,
				     MeshFunction<bool>& masked_vertices,
				     Vector* node_values, NodeNormal& node_normal) :
  mesh(mesh), node_values(node_values), masked_vertices(masked_vertices),
  f(mesh), bcf(mesh), node_normal(node_normal),
  bcslip(mesh, dboundary, node_normal),
  bc0(bcf, mesh, dboundary),
  bc1(mesh, dboundary, masked_vertices, false, node_values)
{

  int d = mesh.topology().dim();

  cout << "LaplacianSmoother ctor" << endl;
  storeA = new Matrix;
  storeb = new Vector;

  if(d == 2)
  {
    a = new Laplacian2DBilinearForm;
    L = new Laplacian2DLinearForm(f);
  }
  else if(d == 3)
  {
    a = new Laplacian3DBilinearForm;
    L = new Laplacian3DLinearForm(f);
  }
  
}

LaplacianSmoother::~LaplacianSmoother(){


  delete L;
  delete a;
  delete storeb;
  delete storeA;
}

void LaplacianSmoother::smooth(Vector& motionx,
			       bool reset, bool with_slip_bc)
{
  cout << "laplacian smooth" << endl;
  int d = mesh.topology().dim();

  real timer = time();

  node_normal.__compute_normal(mesh);

  Matrix& A = *storeA;
  Vector& b = *storeb;

  Assembler assembler(mesh);

  assembler.assemble(A, *a, reset);
  assembler.assemble(b, *L, reset);

  message("LaplacianSmoother timer smooth assemble: %g", time() - timer);

  real timer2 = time();
  message("LaplacianSmoother timer smooth bc: %g", time() - timer2);
  timer2 = time();

  if(!with_slip_bc)
  {
    bc0.apply(A, b, *a);
    bc1.apply(A, b, *a);
    KrylovSolver ksolver(cg);
    ksolver.solve(A, motionx, b);
  }
  else
  {
    bcslip.update();
    bcslip.apply(A, b, *a);
    bc1.apply(A, b, *a);
    KrylovSolver ksolver(bicgstab);
    ksolver.solve(A, motionx, b);
  }
  message("LaplacianSmoother timer smooth linear solve: %g", time() - timer2);
  message("LaplacianSmoother timer smooth: %g", time() - timer);
}
//-----------------------------------------------------------------------------
