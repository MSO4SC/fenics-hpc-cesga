// Copyright (C) 2017 Niyazi Cem Degirmenci
// Licensed under the GNU LGPL Version 2.1.
//
// First added:  2017-04-27
// Last changed:  2017-04-27
//
//
// refines the mesh given by first argument and saves into second argument

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
#include <dolfin/mesh/RivaraRefinement.h>

using namespace dolfin;

int main(int argc, char **argv)
{
 dolfin_init(argc,argv);

dolfin_set("output destination","silent");
if(dolfin::MPI::processNumber() == 0)
        dolfin_set("output destination","terminal");

{

  Mesh mesh(argv[1]);
  MeshFunction<bool> cell_refinement_marker(mesh);
  cell_refinement_marker.init(mesh.topology().dim());
  cell_refinement_marker = true;
    

  RivaraRefinement::refine(mesh, cell_refinement_marker);

  File fout(argv[2]);
  fout << mesh;

   
}
dolfin_finalize();
return 0;
}
