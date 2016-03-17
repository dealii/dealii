// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2016 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>

#include <deal.II/lac/parallel_vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/matrix_free/fe_evaluation.h>

#include <algorithm>

DEAL_II_NAMESPACE_OPEN


template<int dim, typename Number>
MGTransferMatrixFree<dim,Number>::MGTransferMatrixFree ()
  :
  fe_degree(0),
  element_is_continuous(false),
  n_components(0),
  n_child_cell_dofs(0)
{}



template<int dim, typename Number>
MGTransferMatrixFree<dim,Number>::MGTransferMatrixFree (const MGConstrainedDoFs &mg_c)
  :
  fe_degree(0),
  element_is_continuous(false),
  n_components(0),
  n_child_cell_dofs(0)
{
  this->mg_constrained_dofs = &mg_c;
}



template <int dim, typename Number>
MGTransferMatrixFree<dim,Number>::~MGTransferMatrixFree ()
{}



template <int dim, typename Number>
void MGTransferMatrixFree<dim,Number>::initialize_constraints
(const MGConstrainedDoFs &mg_c)
{
  this->mg_constrained_dofs = &mg_c;
}



template <int dim, typename Number>
void MGTransferMatrixFree<dim,Number>::clear ()
{
  this->MGLevelGlobalTransfer<parallel::distributed::Vector<Number> >::clear();
  fe_degree = 0;
  element_is_continuous = false;
  n_components = 0;
  n_child_cell_dofs = 0;
  level_dof_indices.clear();
  parent_child_connect.clear();
  n_owned_level_cells.clear();
  shape_info = internal::MatrixFreeFunctions::ShapeInfo<Number>();
  evaluation_data.clear();
  weights_on_refined.clear();
}



namespace
{
  // given the collection of child cells in lexicographic ordering as seen
  // from the parent, compute the first index of the given child
  template <int dim>
  unsigned int
  compute_shift_within_children(const unsigned int child,
                                const unsigned int fe_shift_1d,
                                const unsigned int fe_degree)
  {
    // we put the degrees of freedom of all child cells in
    // lexicographic ordering
    unsigned int c_tensor_index[dim];
    unsigned int tmp = child;
    for (unsigned int d=0; d<dim; ++d)
      {
        c_tensor_index[d] = tmp % 2;
        tmp /= 2;
      }
    const unsigned int n_child_dofs_1d = fe_degree + 1 + fe_shift_1d;
    unsigned int factor = 1;
    unsigned int shift = fe_shift_1d * c_tensor_index[0];
    for (unsigned int d=1; d<dim; ++d)
      {
        factor *= n_child_dofs_1d;
        shift = shift + factor * fe_shift_1d * c_tensor_index[d];
      }
    return shift;
  }



  // puts the indices on the given child cell in lexicographic ordering with
  // respect to the collection of all child cells as seen from the parent
  template <int dim>
  void add_child_indices(const unsigned int child,
                         const unsigned int fe_shift_1d,
                         const unsigned int fe_degree,
                         const std::vector<unsigned int> &lexicographic_numbering,
                         const std::vector<types::global_dof_index> &local_dof_indices,
                         types::global_dof_index *target_indices)
  {
    const unsigned int n_child_dofs_1d = fe_degree + 1 + fe_shift_1d;
    const unsigned int shift =
      compute_shift_within_children<dim>(child, fe_shift_1d, fe_degree);
    const unsigned int n_components =
      local_dof_indices.size()/Utilities::fixed_power<dim>(fe_degree+1);
    types::global_dof_index *indices = target_indices + shift;
    const unsigned int n_scalar_cell_dofs = Utilities::fixed_power<dim>(n_child_dofs_1d);
    for (unsigned int c=0, m=0; c<n_components; ++c)
      for (unsigned int k=0; k<(dim>2 ? (fe_degree+1) : 1); ++k)
        for (unsigned int j=0; j<(dim>1 ? (fe_degree+1) : 1); ++j)
          for (unsigned int i=0; i<(fe_degree+1); ++i, ++m)
            {
              const unsigned int index = c*n_scalar_cell_dofs+k*n_child_dofs_1d*
                                         n_child_dofs_1d+j*n_child_dofs_1d+i;
              Assert(indices[index] == numbers::invalid_dof_index ||
                     indices[index] == local_dof_indices[lexicographic_numbering[m]],
                     ExcInternalError());
              indices[index] = local_dof_indices[lexicographic_numbering[m]];
            }
  }



  // initialize the vectors needed for the transfer (and merge with the
  // content in copy_indices_global_mine)
  template <typename Number>
  void
  reinit_ghosted_vector(const IndexSet &locally_owned,
                        std::vector<types::global_dof_index> &ghosted_level_dofs,
                        const MPI_Comm &communicator,
                        parallel::distributed::Vector<Number> &ghosted_level_vector,
                        std::vector<std::pair<unsigned int,unsigned int> > &copy_indices_global_mine)
  {
    std::sort(ghosted_level_dofs.begin(), ghosted_level_dofs.end());
    IndexSet ghosted_dofs(locally_owned.size());
    ghosted_dofs.add_indices(ghosted_level_dofs.begin(),
                             std::unique(ghosted_level_dofs.begin(),
                                         ghosted_level_dofs.end()));
    ghosted_dofs.compress();

    // Add possible ghosts from the previous content in the vector
    if (ghosted_level_vector.size() == locally_owned.size())
      {
        // shift the local number of the copy indices according to the new
        // partitioner that we are going to use for the vector
        const std_cxx11::shared_ptr<const Utilities::MPI::Partitioner> part
          = ghosted_level_vector.get_partitioner();
        ghosted_dofs.add_indices(part->ghost_indices());
        for (unsigned int i=0; i<copy_indices_global_mine.size(); ++i)
          copy_indices_global_mine[i].second =
            locally_owned.n_elements() +
            ghosted_dofs.index_within_set(part->local_to_global(copy_indices_global_mine[i].second));
      }
    ghosted_level_vector.reinit(locally_owned, ghosted_dofs, communicator);
  }

  // Transform the ghost indices to local index space for the vector
  void
  copy_indices_to_mpi_local_numbers(const Utilities::MPI::Partitioner &part,
                                    const std::vector<types::global_dof_index> &mine,
                                    const std::vector<types::global_dof_index> &remote,
                                    std::vector<unsigned int> &localized_indices)
  {
    localized_indices.resize(mine.size()+remote.size(),
                             numbers::invalid_unsigned_int);
    for (unsigned int i=0; i<mine.size(); ++i)
      if (mine[i] != numbers::invalid_dof_index)
        localized_indices[i] = part.global_to_local(mine[i]);

    for (unsigned int i=0; i<remote.size(); ++i)
      if (remote[i] != numbers::invalid_dof_index)
        localized_indices[i+mine.size()] = part.global_to_local(remote[i]);
  }
}



template <int dim, typename Number>
void MGTransferMatrixFree<dim,Number>::build
(const DoFHandler<dim,dim>  &mg_dof)
{
  this->fill_and_communicate_copy_indices(mg_dof);

  // we collect all child DoFs of a mother cell together. For faster
  // tensorized operations, we align the degrees of freedom
  // lexicographically. We distinguish FE_Q elements and FE_DGQ elements

  const Triangulation<dim> &tria = mg_dof.get_triangulation();

  // ---------------------------- 1. Extract 1D info about the finite element
  // step 1.1: create a 1D copy of the finite element from FETools where we
  // substitute the template argument
  AssertDimension(mg_dof.get_fe().n_base_elements(), 1);
  std::string fe_name = mg_dof.get_fe().base_element(0).get_name();
  {
    const std::size_t template_starts = fe_name.find_first_of('<');
    Assert (fe_name[template_starts+1] == (dim==1?'1':(dim==2?'2':'3')),
            ExcInternalError());
    fe_name[template_starts+1] = '1';
  }
  std_cxx11::shared_ptr<FiniteElement<1> > fe_1d
  (FETools::get_fe_from_name<1>(fe_name));
  const FiniteElement<1> &fe = *fe_1d;
  unsigned int n_child_dofs_1d = numbers::invalid_unsigned_int;

  {
    // currently, we have only FE_Q and FE_DGQ type elements implemented
    n_components = mg_dof.get_fe().element_multiplicity(0);
    AssertDimension(Utilities::fixed_power<dim>(fe.dofs_per_cell)*n_components,
                    mg_dof.get_fe().dofs_per_cell);
    AssertDimension(fe.degree, mg_dof.get_fe().degree);
    fe_degree = fe.degree;
    element_is_continuous = fe.dofs_per_vertex > 0;
    Assert(fe.dofs_per_vertex < 2, ExcNotImplemented());

    // step 1.2: get renumbering of 1D basis functions to lexicographic
    // numbers. The distinction according to fe.dofs_per_vertex is to support
    // both continuous and discontinuous bases.
    std::vector<unsigned int> renumbering(fe.dofs_per_cell);
    {
      AssertIndexRange(fe.dofs_per_vertex, 2);
      renumbering[0] = 0;
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        renumbering[i+fe.dofs_per_vertex] =
          GeometryInfo<1>::vertices_per_cell*fe.dofs_per_vertex + i;
      if (fe.dofs_per_vertex > 0)
        renumbering[fe.dofs_per_cell-fe.dofs_per_vertex] = fe.dofs_per_vertex;
    }

    // step 1.3: create a 1D quadrature formula from the finite element that
    // collects the support points of the basis functions on the two children.
    std::vector<Point<1> > basic_support_points = fe.get_unit_support_points();
    Assert(fe.dofs_per_vertex == 0 || fe.dofs_per_vertex == 1,
           ExcNotImplemented());
    std::vector<Point<1> > points_refined(fe.dofs_per_vertex > 0 ?
                                          (2 * fe.dofs_per_cell - 1) :
                                          (2 * fe.dofs_per_cell));
    const unsigned int shift = fe.dofs_per_cell - fe.dofs_per_vertex;
    for (unsigned int c=0; c<GeometryInfo<1>::max_children_per_cell; ++c)
      for (unsigned int j=0; j<basic_support_points.size(); ++j)
        points_refined[shift*c+j][0] =
          c*0.5 + 0.5 * basic_support_points[renumbering[j]][0];

    n_child_dofs_1d = points_refined.size();
    n_child_cell_dofs = n_components*Utilities::fixed_power<dim>(n_child_dofs_1d);

    // step 1.4: evaluate the polynomials and store the data in ShapeInfo
    const Quadrature<1> quadrature(points_refined);
    shape_info.reinit(quadrature, mg_dof.get_fe(), 0);

    for (unsigned int c=0; c<GeometryInfo<1>::max_children_per_cell; ++c)
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
          Assert(std::abs(shape_info.shape_values[i*n_child_dofs_1d+j+c*shift][0] -
                          fe.get_prolongation_matrix(c)(renumbering[j],renumbering[i]))
                 < std::max(2.*(double)std::numeric_limits<Number>::epsilon(),1e-12),
                 ExcInternalError());
  }

  // -------------- 2. Extract and match dof indices between child and parent
  const unsigned int n_levels = tria.n_global_levels();
  level_dof_indices.resize(n_levels);
  parent_child_connect.resize(n_levels-1);
  n_owned_level_cells.resize(n_levels-1);
  std::vector<std::vector<unsigned int> > coarse_level_indices(n_levels-1);
  for (unsigned int level=0; level<std::min(tria.n_levels(),n_levels-1); ++level)
    coarse_level_indices[level].resize(tria.n_raw_cells(level),
                                       numbers::invalid_unsigned_int);
  std::vector<types::global_dof_index> local_dof_indices(mg_dof.get_fe().dofs_per_cell);
  dirichlet_indices.resize(n_levels-1);

  // We use the vectors stored ghosted_level_vector in the base class for
  // keeping ghosted transfer indices. To avoid keeping two very similar
  // vectors, we merge them here.
  if (this->ghosted_level_vector.max_level() != n_levels-1)
    this->ghosted_level_vector.resize(0, n_levels-1);

  for (unsigned int level=n_levels-1; level > 0; --level)
    {
      unsigned int counter = 0;
      std::vector<types::global_dof_index> global_level_dof_indices;
      std::vector<types::global_dof_index> global_level_dof_indices_remote;
      std::vector<types::global_dof_index> ghosted_level_dofs;
      std::vector<types::global_dof_index> global_level_dof_indices_l0;
      std::vector<types::global_dof_index> ghosted_level_dofs_l0;

      // step 2.1: loop over the cells on the coarse side
      typename DoFHandler<dim>::cell_iterator cell, endc = mg_dof.end(level-1);
      for (cell = mg_dof.begin(level-1); cell != endc; ++cell)
        {
          // need to look into a cell if it has children and it is locally owned
          if (!cell->has_children())
            continue;

          bool consider_cell = false;
          if (tria.locally_owned_subdomain()==numbers::invalid_subdomain_id
              || cell->level_subdomain_id()==tria.locally_owned_subdomain()
             )
            consider_cell = true;

          // due to the particular way we store DoF indices (via children), we
          // also need to add the DoF indices for coarse cells where we own at
          // least one child
          bool cell_is_remote = !consider_cell;
          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            if (cell->child(c)->level_subdomain_id()==tria.locally_owned_subdomain())
              {
                consider_cell = true;
                break;
              }

          if (!consider_cell)
            continue;

          // step 2.2: loop through children and append the dof indices to the
          // appropriate list. We need separate lists for the owned coarse
          // cell case (which will be part of restriction/prolongation between
          // level-1 and level) and the remote case (which needs to store DoF
          // indices for the operations between level and level+1).
          AssertDimension(cell->n_children(),
                          GeometryInfo<dim>::max_children_per_cell);
          std::vector<types::global_dof_index> &next_indices =
            cell_is_remote ? global_level_dof_indices_remote : global_level_dof_indices;
          const std::size_t start_index = next_indices.size();
          next_indices.resize(start_index + n_child_cell_dofs,
                              numbers::invalid_dof_index);
          for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
            {
              if (cell_is_remote && cell->child(c)->level_subdomain_id() !=
                  tria.locally_owned_subdomain())
                continue;
              cell->child(c)->get_mg_dof_indices(local_dof_indices);

              const IndexSet &owned_level_dofs = mg_dof.locally_owned_mg_dofs(level);
              for (unsigned int i=0; i<local_dof_indices.size(); ++i)
                if (!owned_level_dofs.is_element(local_dof_indices[i]))
                  ghosted_level_dofs.push_back(local_dof_indices[i]);

              add_child_indices<dim>(c, fe.dofs_per_cell - fe.dofs_per_vertex,
                                     fe.degree, shape_info.lexicographic_numbering,
                                     local_dof_indices,
                                     &next_indices[start_index]);

              // step 2.3 store the connectivity to the parent
              if (cell->child(c)->has_children() &&
                  (tria.locally_owned_subdomain()==numbers::invalid_subdomain_id
                   || cell->child(c)->level_subdomain_id()==tria.locally_owned_subdomain()
                  ))
                {
                  const unsigned int child_index = coarse_level_indices[level][cell->child(c)->index()];
                  AssertIndexRange(child_index, parent_child_connect[level].size());
                  unsigned int parent_index = counter;
                  // remote cells, i.e., cells where we work on a further
                  // level but are not treated on the current level, need to
                  // be placed at the end of the list; however, we do not yet
                  // know the exact position in the array, so shift their
                  // parent index by the number of cells so we can set the
                  // correct number after the end of this loop
                  if (cell_is_remote)
                    parent_index = start_index/n_child_cell_dofs + tria.n_cells(level);
                  parent_child_connect[level][child_index] =
                    std::make_pair(parent_index, c);
                  AssertIndexRange(mg_dof.get_fe().dofs_per_cell,
                                   static_cast<unsigned short>(-1));

                  // set Dirichlet boundary conditions (as a list of
                  // constrained DoFs) for the child
                  if (this->mg_constrained_dofs != 0)
                    for (unsigned int i=0; i<mg_dof.get_fe().dofs_per_cell; ++i)
                      if (this->mg_constrained_dofs->is_boundary_index(level, local_dof_indices[shape_info.lexicographic_numbering[i]]))
                        dirichlet_indices[level][child_index].push_back(i);
                }
            }
          if (!cell_is_remote)
            {
              AssertIndexRange(static_cast<unsigned int>(cell->index()),
                               coarse_level_indices[level-1].size());
              coarse_level_indices[level-1][cell->index()] = counter++;
            }

          // step 2.4: include indices for the coarsest cells. we still insert
          // the indices as if they were from a child in order to use the same
          // code (the coarsest level does not matter much in terms of memory,
          // so we gain in code simplicity)
          if (level == 1 && !cell_is_remote)
            {
              cell->get_mg_dof_indices(local_dof_indices);

              const IndexSet &owned_level_dofs_l0 = mg_dof.locally_owned_mg_dofs(0);
              for (unsigned int i=0; i<local_dof_indices.size(); ++i)
                if (!owned_level_dofs_l0.is_element(local_dof_indices[i]))
                  ghosted_level_dofs_l0.push_back(local_dof_indices[i]);

              const std::size_t start_index = global_level_dof_indices_l0.size();
              global_level_dof_indices_l0.resize(start_index+n_child_cell_dofs,
                                                 numbers::invalid_dof_index);
              add_child_indices<dim>(0, fe.dofs_per_cell - fe.dofs_per_vertex,
                                     fe.degree, shape_info.lexicographic_numbering,
                                     local_dof_indices,
                                     &global_level_dof_indices_l0[start_index]);

              dirichlet_indices[0].push_back(std::vector<unsigned short>());
              if (this->mg_constrained_dofs != 0)
                for (unsigned int i=0; i<mg_dof.get_fe().dofs_per_cell; ++i)
                  if (this->mg_constrained_dofs->is_boundary_index(0, local_dof_indices[shape_info.lexicographic_numbering[i]]))
                    dirichlet_indices[0].back().push_back(i);
            }
        }

      // step 2.5: store information about the current level and prepare the
      // Dirichlet indices and parent-child relationship for the next coarser
      // level
      AssertDimension(counter*n_child_cell_dofs, global_level_dof_indices.size());
      n_owned_level_cells[level-1] = counter;
      dirichlet_indices[level-1].resize(counter);
      parent_child_connect[level-1].
      resize(counter, std::make_pair(numbers::invalid_unsigned_int,
                                     numbers::invalid_unsigned_int));

      // step 2.6: put the cells with remotely owned parent to the end of the
      // list (these are needed for the transfer from level to level+1 but not
      // for the transfer from level-1 to level).
      if (level < n_levels-1)
        for (std::vector<std::pair<unsigned int,unsigned int> >::iterator
             i=parent_child_connect[level].begin(); i!=parent_child_connect[level].end(); ++i)
          if (i->first >= tria.n_cells(level))
            {
              i->first -= tria.n_cells(level);
              i->first += counter;
            }

      // step 2.7: Initialize the ghosted vector
      const parallel::Triangulation<dim,dim> *ptria =
        (dynamic_cast<const parallel::Triangulation<dim,dim>*> (&tria));
      const MPI_Comm communicator =
        ptria != 0 ? ptria->get_communicator() : MPI_COMM_SELF;

      reinit_ghosted_vector(mg_dof.locally_owned_mg_dofs(level),
                            ghosted_level_dofs, communicator,
                            this->ghosted_level_vector[level],
                            this->copy_indices_global_mine[level]);

      copy_indices_to_mpi_local_numbers(*this->ghosted_level_vector[level].get_partitioner(),
                                        global_level_dof_indices,
                                        global_level_dof_indices_remote,
                                        level_dof_indices[level]);

      // step 2.8: Initialize the ghosted vector for level 0
      if (level == 1)
        {
          for (unsigned int i = 0; i<parent_child_connect[0].size(); ++i)
            parent_child_connect[0][i] = std::make_pair(i, 0U);

          reinit_ghosted_vector(mg_dof.locally_owned_mg_dofs(0),
                                ghosted_level_dofs_l0, communicator,
                                this->ghosted_level_vector[0],
                                this->copy_indices_global_mine[0]);

          copy_indices_to_mpi_local_numbers(*this->ghosted_level_vector[0].get_partitioner(),
                                            global_level_dof_indices_l0,
                                            std::vector<types::global_dof_index>(),
                                            level_dof_indices[0]);
        }
    }

  // ------------------------ 3. compute weights to make restriction additive
  //
  // get the valence of the individual components and compute the weights as
  // the inverse of the valence
  weights_on_refined.resize(n_levels);
  for (unsigned int level = 1; level<n_levels; ++level)
    {
      this->ghosted_level_vector[level] = 0;
      for (unsigned int c=0; c<n_owned_level_cells[level-1]; ++c)
        for (unsigned int j=0; j<n_child_cell_dofs; ++j)
          this->ghosted_level_vector[level].local_element(level_dof_indices[level][n_child_cell_dofs*c+j]) += Number(1.);
      this->ghosted_level_vector[level].compress(VectorOperation::add);
      this->ghosted_level_vector[level].update_ghost_values();

      const unsigned int vec_size = VectorizedArray<Number>::n_array_elements;
      std::vector<unsigned int> degree_to_3 (n_child_dofs_1d);
      degree_to_3[0] = 0;
      for (unsigned int i=1; i<n_child_dofs_1d-1; ++i)
        degree_to_3[i] = 1;
      degree_to_3.back() = 2;

      // we only store 3^dim weights because all dofs on a line have the same
      // valence, and all dofs on a quad have the same valence.
      weights_on_refined[level].resize(((n_owned_level_cells[level-1]+vec_size-1)/vec_size)*Utilities::fixed_power<dim>(3));
      for (unsigned int c=0; c<n_owned_level_cells[level-1]; ++c)
        {
          const unsigned int comp = c/vec_size;
          const unsigned int v = c%vec_size;

          for (unsigned int k=0, m=0; k<(dim>2 ? n_child_dofs_1d : 1); ++k)
            for (unsigned int j=0; j<(dim>1 ? n_child_dofs_1d : 1); ++j)
              {
                unsigned int shift = 9*degree_to_3[k] + 3*degree_to_3[j];
                for (unsigned int i=0; i<n_child_dofs_1d; ++i, ++m)
                  weights_on_refined[level][comp*Utilities::fixed_power<dim>(3)+shift+degree_to_3[i]][v] = Number(1.)/
                      this->ghosted_level_vector[level].local_element(level_dof_indices[level][n_child_cell_dofs*c+m]);
              }
        }
    }

  evaluation_data.resize(3*n_child_cell_dofs);
}



template <int dim, typename Number>
void MGTransferMatrixFree<dim,Number>
::prolongate (const unsigned int                           to_level,
              parallel::distributed::Vector<Number>       &dst,
              const parallel::distributed::Vector<Number> &src) const
{
  Assert ((to_level >= 1) && (to_level<=level_dof_indices.size()),
          ExcIndexRange (to_level, 1, level_dof_indices.size()+1));

  AssertDimension(this->ghosted_level_vector[to_level].local_size(),
                  dst.local_size());
  AssertDimension(this->ghosted_level_vector[to_level-1].local_size(),
                  src.local_size());

  this->ghosted_level_vector[to_level-1] = src;
  this->ghosted_level_vector[to_level-1].update_ghost_values();
  this->ghosted_level_vector[to_level] = 0.;

  // the implementation in do_prolongate_add is templated in the degree of the
  // element (for efficiency reasons), so we need to find the appropriate
  // kernel here...
  if (fe_degree == 0)
    do_prolongate_add<0>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 1)
    do_prolongate_add<1>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 2)
    do_prolongate_add<2>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 3)
    do_prolongate_add<3>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 4)
    do_prolongate_add<4>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 5)
    do_prolongate_add<5>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 6)
    do_prolongate_add<6>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 7)
    do_prolongate_add<7>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 8)
    do_prolongate_add<8>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 9)
    do_prolongate_add<9>(to_level, this->ghosted_level_vector[to_level],
                         this->ghosted_level_vector[to_level-1]);
  else if (fe_degree == 10)
    do_prolongate_add<10>(to_level, this->ghosted_level_vector[to_level],
                          this->ghosted_level_vector[to_level-1]);
  else
    AssertThrow(false, ExcNotImplemented("Only degrees 0 up to 10 implemented."));

  this->ghosted_level_vector[to_level].compress(VectorOperation::add);
  dst = this->ghosted_level_vector[to_level];
}



template <int dim, typename Number>
void MGTransferMatrixFree<dim,Number>
::restrict_and_add (const unsigned int                           from_level,
                    parallel::distributed::Vector<Number>       &dst,
                    const parallel::distributed::Vector<Number> &src) const
{
  Assert ((from_level >= 1) && (from_level<=level_dof_indices.size()),
          ExcIndexRange (from_level, 1, level_dof_indices.size()+1));

  AssertDimension(this->ghosted_level_vector[from_level].local_size(),
                  src.local_size());
  AssertDimension(this->ghosted_level_vector[from_level-1].local_size(),
                  dst.local_size());

  this->ghosted_level_vector[from_level] = src;
  this->ghosted_level_vector[from_level].update_ghost_values();
  this->ghosted_level_vector[from_level-1] = 0.;

  if (fe_degree == 0)
    do_restrict_add<0>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 1)
    do_restrict_add<1>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 2)
    do_restrict_add<2>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 3)
    do_restrict_add<3>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 4)
    do_restrict_add<4>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 5)
    do_restrict_add<5>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 6)
    do_restrict_add<6>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 7)
    do_restrict_add<7>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 8)
    do_restrict_add<8>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 9)
    do_restrict_add<9>(from_level, this->ghosted_level_vector[from_level-1],
                       this->ghosted_level_vector[from_level]);
  else if (fe_degree == 10)
    do_restrict_add<10>(from_level, this->ghosted_level_vector[from_level-1],
                        this->ghosted_level_vector[from_level]);
  else
    AssertThrow(false, ExcNotImplemented("Only degrees 0 up to 10 implemented."));

  this->ghosted_level_vector[from_level-1].compress(VectorOperation::add);
  dst += this->ghosted_level_vector[from_level-1];
}



namespace
{
  template <int dim, typename Eval, typename Number, bool prolongate>
  void
  perform_tensorized_op(const Eval &evaluator,
                        const unsigned int n_child_cell_dofs,
                        const unsigned int n_components,
                        AlignedVector<VectorizedArray<Number> > &evaluation_data)
  {
    AssertDimension(n_components * Eval::n_q_points, n_child_cell_dofs);
    VectorizedArray<Number> *t0 = &evaluation_data[0];
    VectorizedArray<Number> *t1 = &evaluation_data[n_child_cell_dofs];
    VectorizedArray<Number> *t2 = &evaluation_data[2*n_child_cell_dofs];

    for (unsigned int c=0; c<n_components; ++c)
      {
        // for the prolongate case, we go from dofs (living on the parent cell) to
        // quads (living on all children) in the FEEvaluation terminology
        if (dim == 1)
          evaluator.template values<0,prolongate,false>(t0, t2);
        else if (dim == 2)
          {
            evaluator.template values<0,prolongate,false>(t0, t1);
            evaluator.template values<1,prolongate,false>(t1, t2);
          }
        else if (dim == 3)
          {
            evaluator.template values<0,prolongate,false>(t0, t2);
            evaluator.template values<1,prolongate,false>(t2, t1);
            evaluator.template values<2,prolongate,false>(t1, t2);
          }
        else
          Assert(false, ExcNotImplemented());
        if (prolongate)
          {
            t0 += Eval::dofs_per_cell;
            t2 += Eval::n_q_points;
          }
        else
          {
            t0 += Eval::n_q_points;
            t2 += Eval::dofs_per_cell;
          }
      }
  }

  template <int dim, int degree, typename Number>
  void weight_dofs_on_child (const VectorizedArray<Number> *weights,
                             const unsigned int n_components,
                             VectorizedArray<Number> *data)
  {
    Assert(degree > 0, ExcNotImplemented());
    const int loop_length = 2*degree+1;
    unsigned int degree_to_3 [loop_length];
    degree_to_3[0] = 0;
    for (int i=1; i<loop_length-1; ++i)
      degree_to_3[i] = 1;
    degree_to_3[loop_length-1] = 2;
    for (unsigned int c=0; c<n_components; ++c)
      for (int k=0; k<(dim>2 ? loop_length : 1); ++k)
        for (int j=0; j<(dim>1 ? loop_length : 1); ++j)
          {
            const unsigned int shift = 9*degree_to_3[k] + 3*degree_to_3[j];
            data[0] *= weights[shift];
            // loop bound as int avoids compiler warnings in case loop_length
            // == 1 (polynomial degree 0)
            for (int i=1; i<loop_length-1; ++i)
              data[i] *= weights[shift+1];
            data[loop_length-1] *= weights[shift+2];
            data += loop_length;
          }
  }
}



template <int dim, typename Number>
template <int degree>
void MGTransferMatrixFree<dim,Number>
::do_prolongate_add (const unsigned int                           to_level,
                     parallel::distributed::Vector<Number>       &dst,
                     const parallel::distributed::Vector<Number> &src) const
{
  const unsigned int vec_size = VectorizedArray<Number>::n_array_elements;
  const unsigned int n_child_dofs_1d = 2*(fe_degree+1) - element_is_continuous;
  const unsigned int n_scalar_cell_dofs = Utilities::fixed_power<dim>(n_child_dofs_1d);
  const unsigned int three_to_dim = Utilities::fixed_int_power<3,dim>::value;

  for (unsigned int cell=0; cell < n_owned_level_cells[to_level-1];
       cell += vec_size)
    {
      const unsigned int n_chunks = cell+vec_size > n_owned_level_cells[to_level-1] ?
                                    n_owned_level_cells[to_level-1] - cell : vec_size;

      // read from source vector
      for (unsigned int v=0; v<n_chunks; ++v)
        {
          const unsigned int shift = compute_shift_within_children<dim>
                                     (parent_child_connect[to_level-1][cell+v].second,
                                      degree+1-element_is_continuous, degree);
          const unsigned int *indices = &level_dof_indices[to_level-1][parent_child_connect[to_level-1][cell+v].first*n_child_cell_dofs+shift];
          for (unsigned int c=0, m=0; c<n_components; ++c)
            {
              for (unsigned int k=0; k<(dim>2 ? (degree+1) : 1); ++k)
                for (unsigned int j=0; j<(dim>1 ? (degree+1) : 1); ++j)
                  for (unsigned int i=0; i<(degree+1); ++i, ++m)
                    evaluation_data[m][v] =
                      src.local_element(indices[c*n_scalar_cell_dofs +
                                                k*n_child_dofs_1d*n_child_dofs_1d+
                                                j*n_child_dofs_1d+i]);

              // apply Dirichlet boundary conditions on parent cell
              for (std::vector<unsigned short>::const_iterator i=dirichlet_indices[to_level-1][cell+v].begin(); i!=dirichlet_indices[to_level-1][cell+v].end(); ++i)
                evaluation_data[*i][v] = 0.;
            }
        }

      // perform tensorized operation
      Assert(shape_info.element_type ==
             internal::MatrixFreeFunctions::tensor_symmetric, ExcNotImplemented());
      if (element_is_continuous)
        {
          AssertDimension(shape_info.shape_val_evenodd.size(),
                          (degree+1)*(degree+1));
          typedef internal::EvaluatorTensorProduct<internal::evaluate_evenodd,dim,degree,2*degree+1,VectorizedArray<Number> > Evaluator;
          Evaluator evaluator(shape_info.shape_val_evenodd,
                              shape_info.shape_val_evenodd,
                              shape_info.shape_val_evenodd);
          perform_tensorized_op<dim,Evaluator,Number,true>(evaluator,
                                                           n_child_cell_dofs,
                                                           n_components,
                                                           evaluation_data);
          weight_dofs_on_child<dim,degree,Number>(&weights_on_refined[to_level][(cell/vec_size)*three_to_dim],
                                                  n_components,
                                                  &evaluation_data[2*n_child_cell_dofs]);
        }
      else
        {
          AssertDimension(shape_info.shape_val_evenodd.size(),
                          (degree+1)*(degree+1));
          typedef internal::EvaluatorTensorProduct<internal::evaluate_evenodd,dim,degree,2*degree+2,VectorizedArray<Number> > Evaluator;
          Evaluator evaluator(shape_info.shape_val_evenodd,
                              shape_info.shape_val_evenodd,
                              shape_info.shape_val_evenodd);
          perform_tensorized_op<dim,Evaluator,Number,true>(evaluator,
                                                           n_child_cell_dofs,
                                                           n_components,
                                                           evaluation_data);
        }

      // write into dst vector
      const unsigned int *indices = &level_dof_indices[to_level][cell*
                                                                 n_child_cell_dofs];
      for (unsigned int v=0; v<n_chunks; ++v)
        {
          for (unsigned int i=0; i<n_child_cell_dofs; ++i)
            dst.local_element(indices[i]) += evaluation_data[2*n_child_cell_dofs+i][v];
          indices += n_child_cell_dofs;
        }
    }
}



template <int dim, typename Number>
template <int degree>
void MGTransferMatrixFree<dim,Number>
::do_restrict_add (const unsigned int                           from_level,
                   parallel::distributed::Vector<Number>       &dst,
                   const parallel::distributed::Vector<Number> &src) const
{
  const unsigned int vec_size = VectorizedArray<Number>::n_array_elements;
  const unsigned int n_child_dofs_1d = 2*(fe_degree+1) - element_is_continuous;
  const unsigned int n_scalar_cell_dofs = Utilities::fixed_power<dim>(n_child_dofs_1d);
  const unsigned int three_to_dim = Utilities::fixed_int_power<3,dim>::value;

  for (unsigned int cell=0; cell < n_owned_level_cells[from_level-1];
       cell += vec_size)
    {
      const unsigned int n_chunks = cell+vec_size > n_owned_level_cells[from_level-1] ?
                                    n_owned_level_cells[from_level-1] - cell : vec_size;

      // read from source vector
      {
        const unsigned int *indices = &level_dof_indices[from_level][cell*
                                      n_child_cell_dofs];
        for (unsigned int v=0; v<n_chunks; ++v)
          {
            for (unsigned int i=0; i<n_child_cell_dofs; ++i)
              evaluation_data[i][v] = src.local_element(indices[i]);
            indices += n_child_cell_dofs;
          }
      }

      // perform tensorized operation
      Assert(shape_info.element_type ==
             internal::MatrixFreeFunctions::tensor_symmetric, ExcNotImplemented());
      if (element_is_continuous)
        {
          AssertDimension(shape_info.shape_val_evenodd.size(),
                          (degree+1)*(degree+1));
          typedef internal::EvaluatorTensorProduct<internal::evaluate_evenodd,dim,degree,2*degree+1,VectorizedArray<Number> > Evaluator;
          Evaluator evaluator(shape_info.shape_val_evenodd,
                              shape_info.shape_val_evenodd,
                              shape_info.shape_val_evenodd);
          weight_dofs_on_child<dim,degree,Number>(&weights_on_refined[from_level][(cell/vec_size)*three_to_dim],
                                                  n_components,
                                                  &evaluation_data[0]);
          perform_tensorized_op<dim,Evaluator,Number,false>(evaluator,
                                                            n_child_cell_dofs,
                                                            n_components,
                                                            evaluation_data);
        }
      else
        {
          AssertDimension(shape_info.shape_val_evenodd.size(),
                          (degree+1)*(degree+1));
          typedef internal::EvaluatorTensorProduct<internal::evaluate_evenodd,dim,degree,2*degree+2,VectorizedArray<Number> > Evaluator;
          Evaluator evaluator(shape_info.shape_val_evenodd,
                              shape_info.shape_val_evenodd,
                              shape_info.shape_val_evenodd);
          perform_tensorized_op<dim,Evaluator,Number,false>(evaluator,
                                                            n_child_cell_dofs,
                                                            n_components,
                                                            evaluation_data);
        }

      // write into dst vector
      for (unsigned int v=0; v<n_chunks; ++v)
        {
          const unsigned int shift = compute_shift_within_children<dim>
                                     (parent_child_connect[from_level-1][cell+v].second,
                                      degree+1-element_is_continuous, degree);
          AssertIndexRange(parent_child_connect[from_level-1][cell+v].first*
                           n_child_cell_dofs+n_child_cell_dofs-1,
                           level_dof_indices[from_level-1].size());
          const unsigned int *indices = &level_dof_indices[from_level-1][parent_child_connect[from_level-1][cell+v].first*n_child_cell_dofs+shift];
          for (unsigned int c=0, m=0; c<n_components; ++c)
            {
              // apply Dirichlet boundary conditions on parent cell
              for (std::vector<unsigned short>::const_iterator i=dirichlet_indices[from_level-1][cell+v].begin(); i!=dirichlet_indices[from_level-1][cell+v].end(); ++i)
                evaluation_data[2*n_child_cell_dofs+(*i)][v] = 0.;

              for (unsigned int k=0; k<(dim>2 ? (degree+1) : 1); ++k)
                for (unsigned int j=0; j<(dim>1 ? (degree+1) : 1); ++j)
                  for (unsigned int i=0; i<(degree+1); ++i, ++m)
                    dst.local_element(indices[c*n_scalar_cell_dofs +
                                              k*n_child_dofs_1d*n_child_dofs_1d+
                                              j*n_child_dofs_1d+i])
                    += evaluation_data[2*n_child_cell_dofs+m][v];
            }
        }
    }
}



template <int dim, typename Number>
std::size_t
MGTransferMatrixFree<dim,Number>::memory_consumption() const
{
  std::size_t memory = MGLevelGlobalTransfer<parallel::distributed::Vector<Number> >::memory_consumption();
  memory += MemoryConsumption::memory_consumption(level_dof_indices);
  memory += MemoryConsumption::memory_consumption(parent_child_connect);
  memory += MemoryConsumption::memory_consumption(n_owned_level_cells);
  memory += shape_info.memory_consumption();
  memory += MemoryConsumption::memory_consumption(evaluation_data);
  memory += MemoryConsumption::memory_consumption(weights_on_refined);
  memory += MemoryConsumption::memory_consumption(dirichlet_indices);
  return memory;
}


// explicit instantiation
#include "mg_transfer_matrix_free.inst"


DEAL_II_NAMESPACE_CLOSE
