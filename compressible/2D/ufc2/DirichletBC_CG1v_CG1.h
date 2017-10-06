// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Existing code for Dirichlet BC is used
//
// Modified by Niclas Jansson, 2008-2009.
// Modified by N. Cem Degirmenci, 2017
// First added:  2007-05-01
// Last changed: 2017-08-04


#ifndef __DIRICHLETBCCG1VCG1_H
#define __DIRICHLETBCCG1VCG1_H

#include <dolfin/fem/BoundaryCondition.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/fem/SubSystem.h>
#include <dolfin/la/Matrix.h>
#include <dolfin/la/Vector.h>

#include <ufc.h>

#include <set>

namespace dolfin
{

class DofMap;
class Form;
class Function;
class GenericMatrix;
class GenericVector;
class Mesh;
class SubDomain;

class DirichletBC_CG1v_CG1: public BoundaryCondition

{

public:

  /// Create boundary condition for sub domain
  DirichletBC_CG1v_CG1(Mesh& mesh, SubDomain& sub_domain, uint subsys, Function* value_fun );

  /// Destructor
  ~DirichletBC_CG1v_CG1();


  /// Apply boundary condition to linear system
  void apply(GenericMatrix& A, GenericVector& b, const Form& form);

  /// Apply boundary condition to linear system
  void apply(GenericMatrix& A, GenericVector& b, const DofMap& dof_map,
             const ufc::form& ufc_form);

  /// Apply boundary condition to linear system for a nonlinear problem
  void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
             const Form& form);

  /// Apply boundary condition to linear system for a nonlinear problem
  void apply(GenericMatrix& A, GenericVector& b, const GenericVector& x,
             const DofMap& dof_map, const ufc::form& ufc_form);

private:

  void applyDirichletBC_CG1v_CG1(Matrix& A, Matrix& As, Vector&, Mesh& mesh, uint node,
                   Array<uint>& nodes, real imp);

  // Do: A(row,col) = value   using setblock not setvalue
  void Aset(Matrix& A, uint row, uint col, real value);

  // Do: b(row) = value   using setblock not setvalue
  void bset(Vector& b, uint row, real value);

  // Initialize sub domain markers
  void init(SubDomain& sub_domain);

  void apply(GenericMatrix& A, GenericVector& b, DofMap const& dof_map,
             Form const& form);

  // The mesh
  Mesh& mesh;

  // Sub domain markers (if any)
  MeshFunction<uint>* sub_domains;

  // The sub domain
  uint sub_domain;

  // True if sub domain markers are created locally
  bool sub_domains_local;

  // Sub system
  SubSystem sub_system; // 0 pressure 1 ux 2 uy 3 uz

  uint sub_sys_num;

  // User defined sub domain
  SubDomain* user_sub_domain;


  Matrix* As;

  int N_local;
  int N_offset;
  std::set<uint> off_proc_rows;

  real *row_block;
  real *zero_block;
  uint *a1_indices_array;

  BoundaryMesh* boundary;
  MeshFunction<uint> *cell_map;
  MeshFunction<uint> *vertex_map;
  
  Function* value_fun;
};

//--- INLINE ------------------------------------------------------------------

inline void DirichletBC_CG1v_CG1::Aset(Matrix& A, uint row, uint col, real value)
{
  A.set(&value, 1, &row, 1, &col);
}

inline void DirichletBC_CG1v_CG1::bset(Vector& b, uint row, real value)
{
  b.set(&value, 1, &row);
}

}

#endif
