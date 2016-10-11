//#include <dolfin/common/common_includes.h>
// Copyright (C) 2007 Murtazo Nazarov
// Licensed under the GNU LGPL Version 2.1.
//
// Modified by Niclas Jansson, 2008-2009.
//
// First added:  2007-05-01
// Last changed: 2009-12-30

#include "dolfin/NodeNormal.h"
#include <dolfin/main/MPI.h>
#include <dolfin/math/basic.h>
#include <dolfin/mesh/BoundaryMesh.h>
#include <dolfin/mesh/Cell.h>
#include <dolfin/mesh/Facet.h>
#include <dolfin/mesh/MeshData.h>
#include <dolfin/mesh/Vertex.h>

#include <map>

#define B(row,col,nrow) ((row) + ((nrow)*(col)))

namespace dolfin
{

//-----------------------------------------------------------------------------
NodeNormal::NodeNormal(NodeNormal& node_normal) :
        normal(NULL),
        tau(NULL),
        tau_1(NULL),
        tau_2(NULL),
        mesh(node_normal.mesh),
        alpha_max_(DOLFIN_PI / 4.0),
        weighting_(none)
{
  *this = node_normal;
}

//-----------------------------------------------------------------------------
NodeNormal::NodeNormal(Mesh& mesh) :
        normal(NULL),
        tau(NULL),
        tau_1(NULL),
        tau_2(NULL),
        mesh(mesh),
        alpha_max_(DOLFIN_PI / 4.0),
        weighting_(none)
{
  uint const nsdim = mesh.topology().dim();

  for (uint d = 0; d < nsdim; ++d)
  {
    basis_.push_back(new MeshFunction<real> [nsdim]);
    for (uint i = 0; i < nsdim; ++i)
    {
      basis_.back()[i].init(mesh, 0);
    }
  }
  node_type.init(mesh, 0);

  // backward compatbility
  normal = basis_[0];
  if (nsdim == 2)
  {
    tau = basis_[1];
  }
  if (nsdim == 3)
  {
    tau_1 = basis_[1];
    tau_2 = basis_[2];
  }

  __compute_normal(mesh);
}

//-----------------------------------------------------------------------------
NodeNormal::~NodeNormal()
{
  clear();
}

//-----------------------------------------------------------------------------
void
NodeNormal::clear()
{
  while (!basis_.empty())
  {
    delete[] basis_.back();
    basis_.pop_back();
  }
}

//-----------------------------------------------------------------------------
NodeNormal&
NodeNormal::operator=(NodeNormal& node_normal)
{
  clear();

  normal = new MeshFunction<real> [mesh.topology().dim()];
  tau = new MeshFunction<real> [mesh.topology().dim()];
  tau_1 = new MeshFunction<real> [mesh.topology().dim()];
  tau_2 = new MeshFunction<real> [mesh.topology().dim()];

  uint const nsdim = mesh.topology().dim();

  node_type.init(mesh, 0);
  for (uint d = 0; d < nsdim; ++d)
  {
    basis_.push_back(new MeshFunction<real> [nsdim]);
    for (uint i = 0; i < nsdim; ++i)
    {
      basis_.back()[i].init(mesh, 0);
    }

  }
  for (VertexIterator v(mesh); !v.end(); ++v)
  {
    node_type.set(*v, node_normal.node_type.get(*v));
    for (uint d = 0; d < nsdim; ++d)
    {
      for (uint i = 0; i < nsdim; ++i)
      {
        basis_[d][i].set(*v, node_normal.basis()[d][i].get(*v));
      }
    }
  }

  // backward compatbility
  normal = basis_[0];
  if (nsdim == 2)
  {
    tau = basis_[1];
  }
  if (nsdim == 3)
  {
    tau_1 = basis_[1];
    tau_2 = basis_[2];
  }

  return *this;
}

//-----------------------------------------------------------------------------
void
NodeNormal::__compute_normal(Mesh& mesh)
{
  message("BoundaryNormals: Compute normals");
  mesh.renumber();
  uint rank = dolfin::MPI::processNumber();
  uint pe_size = dolfin::MPI::numProcesses();
  Array<real> *send_buff_type_basis = new Array<real> [pe_size];
  Array<uint> *send_buff_index = new Array<uint> [pe_size];
  _map<uint, bool> used_shared;

  for (MeshSharedIterator s(mesh.distdata(), 0); !s.end(); ++s)
  {
    used_shared[s.index()] = false;
  }

  // Iterate over all cells in the boundary mesh
  BoundaryMesh boundary(mesh);

  //  uint const boundary_nsdim = boundary.topology().dim();
  MeshFunction<uint>* cell_map = boundary.data().meshFunction("cell map");
  MeshFunction<uint>* vertex_map = boundary.data().meshFunction("vertex map");
  uint const nsdim = mesh.topology().dim();
  //
  if (dolfin::MPI::numProcesses() > 1)
  {
    __cache_shared_area(mesh, boundary);
  }

  //-------------------------------------------------------------------------
  real * basis_vec = new real[nsdim * nsdim];

  real * n_block = NULL;
  real * w_block = NULL;

  uint NbNeighCells = 0;
  // Computation of normals to the boundary vertices

  if (boundary.numCells())
  {

    for (VertexIterator boundary_vertex(boundary); !boundary_vertex.end();
        ++boundary_vertex)
    {
      uint boundary_id = vertex_map->get(*boundary_vertex);
      uint global_id = mesh.distdata().get_global(boundary_id, 0);
      //std_out << "global_id = " << global_id << std::endl;

      Point p = boundary_vertex->point();

      bool const vertex_is_shared = mesh.distdata().is_shared(boundary_id, 0);
      bool const vertex_is_ghosted = mesh.distdata().is_ghost(boundary_id, 0);

      // Reset basis to zero
      std::fill_n(basis_vec, nsdim * nsdim, 0.);

      NbNeighCells = 0;
      // Get number of neighbouring cells
      // Add storage for shared vertices cell normals
      if (dolfin::MPI::numProcesses() > 1 && vertex_is_shared)
      {
        NbNeighCells = num_neigh_cells_[global_id];
      }
      else
      {
        //uint nbnc = boundary_vertex->numEntities(boundary_nsdim);
        for (CellIterator boundary_cell(*boundary_vertex); !boundary_cell.end();
            ++boundary_cell)
        {
          ++NbNeighCells;
        }
        // FIXME: there seems to be a bug such that calling numEntities does not initialize the connectivities
        /*
         if (NbNeighCells != nbnc)
         {
         dolfin::error(
         "Different number of neighbouring facets returned by methods.");
         }
         */
      }

      if (NbNeighCells == 0)
      {
        dolfin::error("This vertex was found to have zero neighbouring facets");
      }

      // contains of all normals from the facets
      n_block = new real[NbNeighCells * nsdim];
      // contains of all weights from the facets
      w_block = new real[NbNeighCells];

      // sum of all areas of the elements
      real sum_weights = 0.0;
      uint neighcell_count = 0;

      for (CellIterator boundary_cell(*boundary_vertex); !boundary_cell.end();
          ++boundary_cell)
      {
        // Create mesh facet corresponding to boundary cell
        Facet mesh_facet(mesh, cell_map->get(*boundary_cell));
        dolfin_assert(mesh_facet.numEntities(nsdim) == 1);

        // Get cell to which facet belongs (pick first, there is only one)
        Cell mesh_cell(mesh, mesh_facet.entities(nsdim)[0]);
        // Get local index of facet with respect to the cell
        uint local_facet = mesh_cell.index(mesh_facet);

        // Compute area of facet
        real facet_normal_weight = 1.0;
        if (weighting_ == facet)
        {
          // Compute the measure of the facet
          facet_normal_weight = boundary_cell->volume();
        }
        else if (weighting_ == cell)
        {
          // Compute the measure of the cell
          facet_normal_weight = mesh_cell.volume();
        }

        w_block[neighcell_count] = facet_normal_weight;
        // Get sum of the length/area of all neighbouring cells
        sum_weights += facet_normal_weight;

        // Save facet normals to block
        for (uint d = 0; d < nsdim; ++d)
        {
          real nd = mesh_cell.normal(local_facet, d);
          basis_vec[d] += facet_normal_weight * nd;
          n_block[B(neighcell_count, d, NbNeighCells)] = nd;
        }
        ++neighcell_count;
      }

      // Add neighbouring bound.cell/facet contributions from other processes
      if (dolfin::MPI::numProcesses() > 1 && vertex_is_shared)
      {
        // Add normals computed at vertex from other processes
        uint tabulated_idx = shared_offsetidx_[global_id];
        for (uint d = 0; d < nsdim; ++d)
        {
          basis_vec[d] += shared_vertexnormals_[tabulated_idx + d];
        }
        Array<real> tmp = shared_facetweights_block_[global_id];
        Array<real>::iterator iter;

        int iirow = neighcell_count;
        for (iter = tmp.begin(); iter != tmp.end(); iter++)
        {
          w_block[iirow] = *iter;
          sum_weights += w_block[iirow];
          ++iirow;
        }

        tmp = shared_facetnormals_block_[global_id];
        for (iter = tmp.begin(); iter != tmp.end(); iter += nsdim)
        {
          for (uint d = 0; d < nsdim; ++d)
          {
            n_block[B(neighcell_count, d, NbNeighCells)] = *(iter + d);
          }
          ++neighcell_count;
        }
      }

      // Now, we have:
      //    basis_vec  : Vertex normal (*not* unit !) nsdim components
      //    n_block   : NbNeighCells facet unit normals with nsdim components
      //      w_block   : NbNeighCells facet weights
      //    sum_weights : sum of the facet weights
      Array<real> tangent_weights;
      uint vertex_type = 0; // if at the end the value is still zero then something is wrong

      real alpha_max1 = DOLFIN_PI/3.0;
      real alpha_max2 = DOLFIN_PI/4.0;
      real alpha_max3 = DOLFIN_PI/6.0;

      real alpha_max;

      if(fabs(p[1]) >= 0.5)
	alpha_max = alpha_max1;
      else
	alpha_max = alpha_max2;



      real cosalpha_max = std::cos(alpha_max);
      real cosalpha = 0.0;
      //real sinalpha_max = std::sin(alpha_max);

      // In 2D there are only and exactly only two facets but let us not
      // speculate maybe we live in the 5th dimension
      real * nref = new real[nsdim];

      // Let first initialize the list of offset to jump from one normal to another
      uint * normals_offsets = new uint[NbNeighCells];
      // We start by jumping across all the normals starting from 0 since
      // facet normal 0 is part of the computation of nS1
      // (but it's clear that normal 0 belongs to the surface it generates...)
      std::fill_n(normals_offsets, NbNeighCells, 1);
      normals_offsets[0] = 0;

      // DEBUG: map facets to surfaces
      std::map<uint, uint> surfaces;

      // add storage for normals to surfaces
      Array<real *> surface_normals;
      Array<real> surface_totalweights;

      // Has to be integer such that comparison to zero holds
      int remaining_normals_count = NbNeighCells;
      uint curr_surface = 0;
      uint curr_facet = 0;
      //      std_out << "Nb Facets = " << remaining_normals_count << std::endl;
      while (remaining_normals_count > 0)
      {
        ++curr_surface;

        // set new reference normal to first remaining normal.
        // Murtazo seems to use the second surface's averaged normal
        // to detect a third plane, should this be fixed ?
        uint nref_idx = normals_offsets[0];
        for (uint d = 0; d < nsdim; ++d)
        {
          nref[d] = n_block[B(nref_idx,d,NbNeighCells)];
        }

        // init storage for surface normal, deleted with cleanup of surface_normals Array
        real * nSx = new real[nsdim];
        real wSx = 0.;
        std::fill_n(nSx, nsdim, 0);

        // loop through remaining normals indexes
        // we read the sequence of offsets to jump through the remaining normals
        // the first sequence is { 0, 1, 1, 1, ...} with NbNeighCells elements
        curr_facet = normals_offsets[0];
        uint offset_to_update = 0;

        for (uint curr_offset = 0;
            curr_facet < NbNeighCells && curr_offset < NbNeighCells;
            curr_facet += normals_offsets[++curr_offset])
        {
          // then let us loop on the other facet normals to compute the scalar product
          cosalpha = 0.0;
          for (uint d = 0; d < nsdim; ++d)
          {
            cosalpha += nref[d] * n_block[B(curr_facet, d, NbNeighCells)];
          }
          if (cosalpha > cosalpha_max)
          {
            surfaces.insert(std::pair<uint, uint>(curr_facet, curr_surface));

            // add contribution to surface normal
            for (uint d = 0; d < nsdim; ++d)
            {
              nSx[d] += w_block[curr_facet]
                  * n_block[B(curr_facet,d,NbNeighCells)];
              wSx += w_block[curr_facet];
            }

            // eliminate from count of remaining normals to discriminate
            --remaining_normals_count;

            // update offset value
            ++normals_offsets[offset_to_update];

          }
          else
          {
            // found normal not belonging to the same plane
            // increase offset position then set offset value to one
            normals_offsets[offset_to_update] += normals_offsets[curr_offset]
                - 1;
            normals_offsets[++offset_to_update] = 1;
          }
        }

        // add surface normal to the list of surface normals
        surface_normals.push_back(nSx);

        // next loop we add a new surface and discriminate again across
        // the remaining normals

      }
      vertex_type = curr_surface;
      delete nref;
      delete normals_offsets;

      // n_k    = sum_{i=1}^k nS_i
      // tau1_k = |n_{k-1}|^2 nS_k - (n_{k-1} . nS_k ) n_{k-1}
      // tau2_k = n_k ^ tau1_k

      // In 2d tau2_k = ez = (0,0,1)

      // Now we have a nice list of surface normals...
      // ... and apparently we shamelessly take them as averaged normals and tangents
      // But why would I divide by the sum if I am merely normalizing them right after?
      uint const nb_surfaces = surface_normals.size();

      // Taken from V. John, J. Comp. and Appl. Math. 2002
      // Let's fit everything in basis_vec of size nsdim*nsdim
      // normalize the first vector in basis_vec (which is the normal)
      real normal_nrm = 0.0;
      for (uint d = 0; d < nsdim; ++d)
      {
        normal_nrm += basis_vec[d] * basis_vec[d];
      }
      normal_nrm = std::sqrt(normal_nrm);
      for (uint in = 0; in < nsdim; ++in)
      {
        basis_vec[in] /= normal_nrm;
      }
      // Compute unit tangents from unit normal
      if (nsdim == 2)
      {
        // 2D case, make orthonormal
        basis_vec[2] = -basis_vec[1];
        basis_vec[3] = basis_vec[0];
      }
      else
      {
        // 3D case
        if (vertex_type == 1)
        {
          real norm_inv = 0.0;
          if (std::fabs(basis_vec[0]) >= 0.5 || std::fabs(basis_vec[1]) >= 0.5)
          {
            norm_inv = 1.
                / std::sqrt(
                    basis_vec[0] * basis_vec[0] + basis_vec[1] * basis_vec[1]);
            // t11 = n2/n
            basis_vec[3] = basis_vec[1] * norm_inv;
            // t12 = -n1/n
            basis_vec[4] = -basis_vec[0] * norm_inv;
            // t13 = 0
            basis_vec[5] = 0.0;
            // t21 = -t12*n3
            basis_vec[6] = -basis_vec[4] * basis_vec[2];
            // t22 = t11*n3
            basis_vec[7] = basis_vec[3] * basis_vec[2];
            // t23 = t12*n1 - t11*n2
            basis_vec[8] = basis_vec[4] * basis_vec[0]
                - basis_vec[3] * basis_vec[1];
          }
          else
          {
            norm_inv = 1.
                / std::sqrt(
                    basis_vec[1] * basis_vec[1] + basis_vec[2] * basis_vec[2]);
            // t11 = 0
            basis_vec[3] = 0.0;
            // t12 = -n3/n
            basis_vec[4] = -basis_vec[2] * norm_inv;
            // t13 = n2/n
            basis_vec[5] = basis_vec[1] * norm_inv;
            // t21 = t13*n2 - t12*n3
            basis_vec[6] = basis_vec[5] * basis_vec[1]
                - basis_vec[4] * basis_vec[2];
            // t22 = -t13*n1
            basis_vec[7] = -basis_vec[5] * basis_vec[0];
            // t23 = t12*n1
            basis_vec[8] = basis_vec[4] * basis_vec[0];
          }
        }
        else
        {
          uint const Sk = vertex_type - 1;
          basis_vec[6] = basis_vec[1] * surface_normals[Sk][2] - surface_normals[Sk][1] * basis_vec[2];
          basis_vec[7] = basis_vec[2] * surface_normals[Sk][0] - surface_normals[Sk][2] * basis_vec[0];
          basis_vec[8] = basis_vec[0] * surface_normals[Sk][1] - surface_normals[Sk][0] * basis_vec[1];

          real tau2_norm = std::sqrt(basis_vec[6]*basis_vec[6]+basis_vec[7]*basis_vec[7]+basis_vec[8]*basis_vec[8]);

          basis_vec[6] /= tau2_norm;
          basis_vec[7] /= tau2_norm;
          basis_vec[8] /= tau2_norm;


          basis_vec[3] = basis_vec[7] * basis_vec[2] - basis_vec[8] * basis_vec[1];
          basis_vec[4] = basis_vec[8] * basis_vec[0] - basis_vec[6] * basis_vec[2];
          basis_vec[5] = basis_vec[6] * basis_vec[1] - basis_vec[7] * basis_vec[0];
        }

      }

      //---
      // Prepare data structures
      if (vertex_is_ghosted)
      {
        uint owner = mesh.distdata().get_owner(boundary_id, 0);
        send_buff_type_basis[owner].push_back((double) vertex_type);
        for (uint basisidx = 0; basisidx < nsdim * nsdim; ++basisidx)
        {
          send_buff_type_basis[owner].push_back(basis_vec[basisidx]);
        }
        send_buff_index[owner].push_back(global_id);
      }
      else
      {
        node_type.set(boundary_id, vertex_type);
        if (vertex_type == 0)
        {
          error("Surface multiplicity is equal to zero");
        }
        for (uint basisvec = 0; basisvec < nsdim; ++basisvec)
        {
          for (uint in = 0; in < nsdim; ++in)
          {
            basis_[basisvec][in].set(boundary_id,
                basis_vec[basisvec * nsdim + in]);
          }
        }
        used_shared[boundary_id] = true;
      }

      //

      while (!surface_normals.empty())
      {
        delete[] surface_normals.back();
        surface_normals.pop_back();
      }

      delete[] w_block;
      delete[] n_block;
    }
  }

  // Synchronize basis vectors, node_types across processes
  if (dolfin::MPI::numProcesses() > 1)
  {
    MPI_Status status;
    uint src = 0;
    uint dest = 0;
    int recv_count = 0;
    int recv_size = 0;
    int send_size = 0;
    int recv_count_data = 0;

    // Collect data size
    for (uint i = 0; i < pe_size; i++)
    {
      send_size = send_buff_index[i].size();
      MPI_Reduce(&send_size, &recv_count, 1, MPI_INT, MPI_SUM, i,
          dolfin::MPI::DOLFIN_COMM);

      send_size = send_buff_type_basis[i].size();
      MPI_Reduce(&send_size, &recv_count_data, 1, MPI_INT, MPI_SUM, i,
          dolfin::MPI::DOLFIN_COMM);
    }

    // Storage is (node_type (size = 1), basis (size = nsdim*nsdim))
    uint data_alignment = 1 + nsdim * nsdim;
    uint *recv_index = new uint[recv_count];
    real *recv_type = new real[recv_count_data];

    for (uint i = 1; i < pe_size; i++)
    {
      src = (rank - i + pe_size) % pe_size;
      dest = (rank + i) % pe_size;

      MPI_Sendrecv(&send_buff_index[dest][0], send_buff_index[dest].size(),
          MPI_UNSIGNED, dest, 0, recv_index, recv_count, MPI_UNSIGNED, src, 0,
          dolfin::MPI::DOLFIN_COMM, &status);

      MPI_Sendrecv(&send_buff_type_basis[dest][0],
          send_buff_type_basis[dest].size(), MPI_DOUBLE, dest, 1, recv_type,
          recv_count_data, MPI_DOUBLE, src, 1, dolfin::MPI::DOLFIN_COMM,
          &status);
      MPI_Get_count(&status, MPI_DOUBLE, &recv_size);
      // Insert check if value assigned
      uint idx = 0;
      // Data alignment is n_tau + 1
      for (int j = 0; j < recv_size; j += data_alignment, ++idx)
      {
        uint index = mesh.distdata().get_local(recv_index[idx], 0);
        if (!used_shared[index])
        {
          node_type.set(index, (uint) recv_type[j]);

          uint offset = j + 1;
          for (uint basisvec = 0; basisvec < nsdim; ++basisvec)
          {
            for (uint in = 0; in < nsdim; ++in)
            {
              basis_[basisvec][in].set(index,
                  recv_type[offset + basisvec * nsdim + in]);
            }
          }
          used_shared[index] = true;
        }
      }
    }

    delete[] recv_index;
    delete[] recv_type;

  }

  //--- Cleanup
  for (uint i = 0; i < pe_size; ++i)
  {
    send_buff_type_basis[i].clear();
    send_buff_index[i].clear();
  }
  delete[] send_buff_type_basis;
  delete[] send_buff_index;

  num_neigh_cells_.clear();
  shared_vertexnormals_.clear();
  shared_offsetidx_.clear();

  for (std::map<uint, Array<real> >::iterator it =
      shared_facetnormals_block_.begin();
      it != shared_facetnormals_block_.end(); ++it)
  {
    it->second.clear();
  }
  for (std::map<uint, Array<real> >::iterator it =
      shared_facetweights_block_.begin();
      it != shared_facetweights_block_.end(); ++it)
  {
    it->second.clear();
  }
  shared_facetnormals_block_.clear();
  shared_facetweights_block_.clear();
}

//-----------------------------------------------------------------------------
void
NodeNormal::__cache_shared_area(Mesh& mesh, BoundaryMesh& boundary)
{
  uint const nsdim = mesh.topology().dim();

  int rank = dolfin::MPI::processNumber();
  int pe_size = dolfin::MPI::numProcesses();

  MeshFunction<uint>* cell_map = boundary.data().meshFunction("cell map");
  MeshFunction<uint>* vertex_map = boundary.data().meshFunction("vertex map");

  // Send buff for global indices of shared vertices
  Array<uint> sendbuff_global_vert_indices;
  // Send buff for normal contribution of shared vertices
  Array<real> sendbuff_vertexnormals;
  // Send buff for facets normals associated with shared vertices
  Array<real> sendbuff_facetnormals;
  // Send buff for weights of shared vertices
  Array<real> sendbuff_facetweights;

  // Send buff for offset indices packed by triplet for each shared vertex
  // (NbNeighbouringCells, FacetNormalsOffset, FacetWeightsOffset )
  Array<uint> sendbuff_offset_indices;
  vertex_offset_ = 0;
  facetnormals_offset_ = 0;
  facetweights_offset_ = 0;

  // Map global id of boundary vertices to boundary marker
  std::map<uint, bool> GlobalIdOnBoundary;

  real * const normal = new real[nsdim];

  uint SharedVertexCount = 0;
  uint SharedMeshFacetCount = 0;

  // Computation of normals to the boundary vertices shared between processes
  if (boundary.numCells())
  {
    uint NbNeighCells = 0;
    for (VertexIterator boundary_vertex(boundary); !boundary_vertex.end();
        ++boundary_vertex, SharedMeshFacetCount += NbNeighCells)
    {
      uint boundary_id = vertex_map->get(*boundary_vertex);
      uint global_id = mesh.distdata().get_global(boundary_id, 0);

      // Mark vertex as part of the boundary
      GlobalIdOnBoundary[global_id] = true;

      // Reset normal to zero
      std::fill_n(normal, nsdim, 0.);

      // Cache number of neighbouring cells for each shared vertex
      NbNeighCells = 0;

      bool const vertex_is_shared = mesh.distdata().is_shared(boundary_id, 0);

      if (vertex_is_shared)
      {
        ++SharedVertexCount;
        sendbuff_global_vert_indices.push_back(global_id);

        //
        for (CellIterator boundary_cell(*boundary_vertex); !boundary_cell.end();
            ++boundary_cell)
        {
          // Create mesh facet corresponding to boundary cell
          Facet mesh_facet(mesh, cell_map->get(*boundary_cell));
          dolfin_assert(mesh_facet.numEntities(nsdim) == 1);

          // Get cell to which facet belongs (pick first, there is only one)
          Cell mesh_cell(mesh, mesh_facet.entities(nsdim)[0]);

          // Get local index of facet with respect to the cell
          uint local_facet = mesh_cell.index(mesh_facet);

          // Compute area of facet
          real facet_normal_weight = 1.0;
          if (weighting_ == facet)
          {
            // Compute the measure of the facet
            facet_normal_weight = boundary_cell->volume();
          }
          else if (weighting_ == cell)
          {
            // Compute the measure of the cell
            facet_normal_weight = mesh_cell.volume();
          }
          sendbuff_facetweights.push_back(facet_normal_weight);

          // Compute a normal to the boundary
          for (uint d = 0; d < nsdim; ++d)
          {
            real nd = mesh_cell.normal(local_facet, d);
            normal[d] += facet_normal_weight * nd;
            sendbuff_facetnormals.push_back(nd);
          }

          ++NbNeighCells;
        }
        num_neigh_cells_[global_id] = NbNeighCells;

        for (uint d = 0; d < nsdim; ++d)
        {
          sendbuff_vertexnormals.push_back(normal[d]);
        }

        sendbuff_offset_indices.push_back(NbNeighCells);
        sendbuff_offset_indices.push_back(facetnormals_offset_);
        sendbuff_offset_indices.push_back(facetweights_offset_);

        // Init datastructures for shared data
        shared_offsetidx_[global_id] = vertex_offset_;

        // Update shared vertex offset
        vertex_offset_ += nsdim;
        // Padding for NbNeighCells FacetNormals
        facetnormals_offset_ += nsdim * NbNeighCells;
        // Padding for NbNeighCells FacetWeights
        facetweights_offset_ += NbNeighCells;
      }
    }
  }

  delete[] normal;

  // Exchange values
  MPI_Status status;
  uint src;
  uint dest;

  int sh_vertidx_count = sendbuff_global_vert_indices.size();
  int sh_vertexnormals_count = sendbuff_vertexnormals.size();
  int sh_facetnormals_count = sendbuff_facetnormals.size();
  int sh_facetweights_count = sendbuff_facetweights.size();

  dolfin_assert(sh_vertidx_count == SharedVertexCount); dolfin_assert(sh_vertexnormals_count == SharedVertexCount*nsdim); dolfin_assert(sh_facetnormals_count == SharedMeshFacetCount*nsdim); dolfin_assert(sh_facetweights_count == SharedMeshFacetCount);

  int recv_size_vertidx;
  int recv_size_vertexnormals;
  int recv_size_facetnormals;
  int recv_size_facetweights;

  int recv_count;

  MPI_Barrier(dolfin::MPI::DOLFIN_COMM);
  MPI_Allreduce(&sh_vertidx_count, &recv_size_vertidx, 1, MPI_INT, MPI_MAX,
      dolfin::MPI::DOLFIN_COMM);
  MPI_Allreduce(&sh_vertexnormals_count, &recv_size_vertexnormals, 1, MPI_INT,
      MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&sh_facetweights_count, &recv_size_facetweights, 1, MPI_INT,
      MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  MPI_Allreduce(&sh_facetnormals_count, &recv_size_facetnormals, 1, MPI_INT,
      MPI_MAX, dolfin::MPI::DOLFIN_COMM);

  uint *recv_vertidx = new uint[recv_size_vertidx];
  uint *recv_offsetidx = new uint[offsetidx_padding_ * recv_size_vertidx];

  real *recv_vertexnormals = new real[recv_size_vertexnormals];
  real *recv_facetnormals = new real[recv_size_facetnormals];
  real *recv_facetweights = new real[recv_size_facetweights];

  // Init and fill BoundaryNormals storage
  shared_vertexnormals_.resize(SharedVertexCount * nsdim, 0.);

  // For each process
  for (int proc = 1; proc < pe_size; ++proc)
  {
    src = (rank - proc + pe_size) % pe_size;
    dest = (rank + proc) % pe_size;

    MPI_Sendrecv(&sendbuff_global_vert_indices[0], sh_vertidx_count,
        MPI_UNSIGNED, src, 1, recv_vertidx, recv_size_vertidx, MPI_UNSIGNED,
        dest, 1, dolfin::MPI::DOLFIN_COMM, &status);
    MPI_Get_count(&status, MPI_UNSIGNED, &recv_count);

    MPI_Sendrecv(&sendbuff_vertexnormals[0], sh_vertexnormals_count, MPI_DOUBLE,
        src, 1, recv_vertexnormals, recv_size_vertexnormals, MPI_DOUBLE, dest,
        1, dolfin::MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&sendbuff_facetnormals[0], sh_facetnormals_count, MPI_DOUBLE,
        src, 1, recv_facetnormals, recv_size_facetnormals, MPI_DOUBLE, dest, 1,
        dolfin::MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&sendbuff_facetweights[0], sh_facetweights_count, MPI_DOUBLE,
        src, 1, recv_facetweights, recv_size_facetweights, MPI_DOUBLE, dest, 1,
        dolfin::MPI::DOLFIN_COMM, &status);

    MPI_Sendrecv(&sendbuff_offset_indices[0],
        (offsetidx_padding_ * sh_vertidx_count), MPI_UNSIGNED, src, 1,
        recv_offsetidx, offsetidx_padding_ * recv_size_vertidx, MPI_UNSIGNED,
        dest, 1, dolfin::MPI::DOLFIN_COMM, &status);

    // Index for vertex
    uint vertidx = 0;
    // Index for boundary cells == global mesh facets
    uint facetidx = 0;
    // Number of neighbouring boundary cells
    uint nbneighcells = 0;
    // Offsets
    uint facetnoffset = 0;
    uint weightoffset = 0;
    for (int i = 0; i < recv_count; i++)
    {
      uint glb_index = recv_vertidx[i];
      if (mesh.distdata().have_global(glb_index, 0)
          && GlobalIdOnBoundary.count(glb_index) > 0
          && GlobalIdOnBoundary[glb_index])
      {
        // Compute shared normals
        for (uint d = 0; d < nsdim; ++d)
        {
          shared_vertexnormals_[shared_offsetidx_[glb_index] + d] +=
              recv_vertexnormals[vertidx + d];
        }

        // Get alignment and offsets
        nbneighcells = recv_offsetidx[facetidx];
        facetnoffset = recv_offsetidx[facetidx + 1];
        weightoffset = recv_offsetidx[facetidx + 2];

        // Update number of neighbouring boundary cells
        num_neigh_cells_[glb_index] += nbneighcells;

        // Add corresponding facet normals
        for (uint k = 0; k < nsdim * nbneighcells; ++k)
        {
          shared_facetnormals_block_[glb_index].push_back(
              recv_facetnormals[facetnoffset++]);
        }

        // Add corresponding facet weights
        for (uint k = 0; k < nbneighcells; ++k)
        {
          shared_facetweights_block_[glb_index].push_back(
              recv_facetweights[weightoffset++]);
        }
      }

      // Increase data offsets
      vertidx += nsdim;
      facetidx += offsetidx_padding_;
    }

  }
  delete[] recv_vertidx;
  delete[] recv_vertexnormals;
  delete[] recv_facetnormals;
  delete[] recv_offsetidx;
  delete[] recv_facetweights;
}

//-----------------------------------------------------------------------------

}

