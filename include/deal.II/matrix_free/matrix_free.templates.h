// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2018 by the deal.II authors
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

#ifndef dealii_matrix_free_templates_h
#define dealii_matrix_free_templates_h


#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_poly.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/matrix_free/dof_info.templates.h>
#include <deal.II/matrix_free/face_info.h>
#include <deal.II/matrix_free/face_setup_internal.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/shape_info.templates.h>

#ifdef DEAL_II_WITH_THREADS
#  include <deal.II/base/parallel.h>

#  include <tbb/concurrent_unordered_map.h>
#endif

#include <fstream>


DEAL_II_NAMESPACE_OPEN


// --------------------- MatrixFree -----------------------------------

template <int dim, typename Number>
MatrixFree<dim, Number>::MatrixFree() :
  Subscriptor(),
  indices_are_initialized(false),
  mapping_is_initialized(false)
{}



template <int dim, typename Number>
MatrixFree<dim, Number>::MatrixFree(const MatrixFree<dim, Number> &other) :
  Subscriptor()
{
  copy_from(other);
}



template <int dim, typename Number>
std::pair<unsigned int, unsigned int>
MatrixFree<dim, Number>::create_cell_subrange_hp_by_index(
  const std::pair<unsigned int, unsigned int> &range,
  const unsigned int                           fe_index,
  const unsigned int                           vector_component) const
{
  AssertIndexRange(fe_index, dof_info[vector_component].max_fe_index);
  const std::vector<unsigned int> &fe_indices =
    dof_info[vector_component].cell_active_fe_index;
  if (fe_indices.empty() == true)
    return range;
  else
    {
      // the range over which we are searching must be ordered, otherwise we
      // got a range that spans over too many cells
#ifdef DEBUG
      for (unsigned int i = range.first + 1; i < range.second; ++i)
        Assert(
          fe_indices[i] >= fe_indices[i - 1],
          ExcMessage(
            "Cell range must be over sorted range of fe indices in hp case!"));
      AssertIndexRange(range.first, fe_indices.size() + 1);
      AssertIndexRange(range.second, fe_indices.size() + 1);
#endif
      std::pair<unsigned int, unsigned int> return_range;
      return_range.first = std::lower_bound(&fe_indices[0] + range.first,
                                            &fe_indices[0] + range.second,
                                            fe_index) -
                           &fe_indices[0];
      return_range.second =
        std::lower_bound(&fe_indices[0] + return_range.first,
                         &fe_indices[0] + range.second,
                         fe_index + 1) -
        &fe_indices[0];
      Assert(return_range.first >= range.first &&
               return_range.second <= range.second,
             ExcInternalError());
      return return_range;
    }
}



template <int dim, typename Number>
void
MatrixFree<dim, Number>::renumber_dofs(
  std::vector<types::global_dof_index> &renumbering,
  const unsigned int                    vector_component)
{
  AssertIndexRange(vector_component, dof_info.size());
  dof_info[vector_component].compute_dof_renumbering(renumbering);
}



template <int dim, typename Number>
const DoFHandler<dim> &
MatrixFree<dim, Number>::get_dof_handler(const unsigned int dof_index) const
{
  AssertIndexRange(dof_index, n_components());
  if (dof_handlers.active_dof_handler == DoFHandlers::usual)
    {
      AssertDimension(dof_handlers.dof_handler.size(),
                      dof_handlers.n_dof_handlers);
      return *dof_handlers.dof_handler[dof_index];
    }
  else
    {
      AssertThrow(false, ExcNotImplemented());
      // put pseudo return argument to avoid compiler error, but trigger a
      // segfault in case this is only run in optimized mode
      return *dof_handlers.dof_handler[numbers::invalid_unsigned_int];
    }
}



template <int dim, typename Number>
typename DoFHandler<dim>::cell_iterator
MatrixFree<dim, Number>::get_cell_iterator(const unsigned int macro_cell_number,
                                           const unsigned int vector_number,
                                           const unsigned int dof_index) const
{
  const unsigned int vectorization_length =
    VectorizedArray<Number>::n_array_elements;
  AssertIndexRange(dof_index, dof_handlers.n_dof_handlers);
  AssertIndexRange(macro_cell_number, task_info.cell_partition_data.back());
  AssertIndexRange(vector_number, n_components_filled(macro_cell_number));

  const DoFHandler<dim> *dofh = nullptr;
  if (dof_handlers.active_dof_handler == DoFHandlers::usual)
    {
      AssertDimension(dof_handlers.dof_handler.size(),
                      dof_handlers.n_dof_handlers);
      dofh = dof_handlers.dof_handler[dof_index];
    }
  else
    {
      Assert(false,
             ExcMessage("Cannot return DoFHandler<dim>::cell_iterator "
                        "for underlying DoFHandler!"));
    }

  std::pair<unsigned int, unsigned int> index =
    cell_level_index[macro_cell_number * vectorization_length + vector_number];
  return typename DoFHandler<dim>::cell_iterator(
    &dofh->get_triangulation(), index.first, index.second, dofh);
}



template <int dim, typename Number>
typename hp::DoFHandler<dim>::active_cell_iterator
MatrixFree<dim, Number>::get_hp_cell_iterator(
  const unsigned int macro_cell_number,
  const unsigned int vector_number,
  const unsigned int dof_index) const
{
  constexpr unsigned int vectorization_length =
    VectorizedArray<Number>::n_array_elements;
  AssertIndexRange(dof_index, dof_handlers.n_dof_handlers);
  AssertIndexRange(macro_cell_number, task_info.cell_partition_data.back());
  AssertIndexRange(vector_number, n_components_filled(macro_cell_number));

  Assert(dof_handlers.active_dof_handler == DoFHandlers::hp,
         ExcNotImplemented());
  const hp::DoFHandler<dim> *dofh = dof_handlers.hp_dof_handler[dof_index];
  std::pair<unsigned int, unsigned int> index =
    cell_level_index[macro_cell_number * vectorization_length + vector_number];
  return typename hp::DoFHandler<dim>::cell_iterator(
    &dofh->get_triangulation(), index.first, index.second, dofh);
}



template <int dim, typename Number>
void
MatrixFree<dim, Number>::copy_from(const MatrixFree<dim, Number> &v)
{
  clear();
  dof_handlers              = v.dof_handlers;
  dof_info                  = v.dof_info;
  constraint_pool_data      = v.constraint_pool_data;
  constraint_pool_row_index = v.constraint_pool_row_index;
  mapping_info              = v.mapping_info;
  shape_info                = v.shape_info;
  cell_level_index          = v.cell_level_index;
  task_info                 = v.task_info;
  indices_are_initialized   = v.indices_are_initialized;
  mapping_is_initialized    = v.mapping_is_initialized;
}



template <int dim, typename Number>
void
MatrixFree<dim, Number>::internal_reinit(
  const Mapping<dim> &                                    mapping,
  const std::vector<const DoFHandler<dim> *> &            dof_handler,
  const std::vector<const ConstraintMatrix *> &           constraint,
  const std::vector<IndexSet> &                           locally_owned_set,
  const std::vector<hp::QCollection<1>> &                 quad,
  const typename MatrixFree<dim, Number>::AdditionalData &additional_data)
{
  // Reads out the FE information and stores the shape function values,
  // gradients and Hessians for quadrature points.
  {
    unsigned int n_fe = 0;
    for (unsigned int no = 0; no < dof_handler.size(); ++no)
      n_fe += dof_handler[no]->get_fe().n_base_elements();
    const unsigned int n_quad = quad.size();
    shape_info.reinit(TableIndices<4>(n_fe, n_quad, 1, 1));
    for (unsigned int no = 0, c = 0; no < dof_handler.size(); no++)
      for (unsigned int b = 0; b < dof_handler[no]->get_fe().n_base_elements();
           ++b, ++c)
        for (unsigned int nq = 0; nq < n_quad; nq++)
          {
            AssertDimension(quad[nq].size(), 1);
            shape_info(c, nq, 0, 0)
              .reinit(quad[nq][0], dof_handler[no]->get_fe(), b);
          }
  }

  if (additional_data.initialize_indices == true)
    {
      clear();
      Assert(dof_handler.size() > 0, ExcMessage("No DoFHandler is given."));
      AssertDimension(dof_handler.size(), constraint.size());
      AssertDimension(dof_handler.size(), locally_owned_set.size());

      // set variables that are independent of FE
      if (Utilities::MPI::job_supports_mpi() == true)
        {
          const parallel::Triangulation<dim> *dist_tria =
            dynamic_cast<const parallel::Triangulation<dim> *>(
              &(dof_handler[0]->get_triangulation()));
          task_info.communicator = dist_tria != nullptr ?
                                     dist_tria->get_communicator() :
                                     MPI_COMM_SELF;
          task_info.my_pid =
            Utilities::MPI::this_mpi_process(task_info.communicator);
          task_info.n_procs =
            Utilities::MPI::n_mpi_processes(task_info.communicator);
        }
      else
        {
          task_info.communicator = MPI_COMM_SELF;
          task_info.my_pid       = 0;
          task_info.n_procs      = 1;
        }

      initialize_dof_handlers(dof_handler, additional_data);
      for (unsigned int no = 0; no < dof_handler.size(); ++no)
        {
          dof_info[no].store_plain_indices =
            additional_data.store_plain_indices;
          dof_info[no].global_base_element_offset =
            no > 0 ? dof_info[no - 1].global_base_element_offset +
                       dof_handler[no - 1]->get_fe().n_base_elements() :
                     0;
        }

        // initialize the basic multithreading information that needs to be
        // passed to the DoFInfo structure
#ifdef DEAL_II_WITH_THREADS
      if (additional_data.tasks_parallel_scheme != AdditionalData::none &&
          MultithreadInfo::n_threads() > 1)
        {
          task_info.scheme =
            internal::MatrixFreeFunctions::TaskInfo::TasksParallelScheme(
              static_cast<int>(additional_data.tasks_parallel_scheme));
          task_info.block_size = additional_data.tasks_block_size;
        }
      else
#endif
        task_info.scheme = internal::MatrixFreeFunctions::TaskInfo::none;

      // set dof_indices together with constraint_indicator and
      // constraint_pool_data. It also reorders the way cells are gone through
      // (to separate cells with overlap to other processors from others
      // without).
      initialize_indices(constraint, locally_owned_set, additional_data);
    }

  // initialize bare structures
  else if (dof_info.size() != dof_handler.size())
    {
      initialize_dof_handlers(dof_handler, additional_data);
      std::vector<unsigned int>  dummy;
      std::vector<unsigned char> dummy2;
      task_info.collect_boundary_cells(
        cell_level_index.size(),
        cell_level_index.size(),
        VectorizedArray<Number>::n_array_elements,
        dummy);
      task_info.create_blocks_serial(
        dummy, dummy, 1, dummy, false, dummy, dummy2);
      for (unsigned int i = 0; i < dof_info.size(); ++i)
        {
          dof_info[i].dimension = dim;
          dof_info[i].n_base_elements =
            dof_handler[i]->get_fe().n_base_elements();
          dof_info[i].n_components.resize(dof_info[i].n_base_elements);
          dof_info[i].start_components.resize(dof_info[i].n_base_elements + 1);
          for (unsigned int c = 0; c < dof_info[i].n_base_elements; ++c)
            {
              dof_info[i].n_components[c] =
                dof_handler[i]->get_fe().element_multiplicity(c);
              for (unsigned int l = 0; l < dof_info[i].n_components[c]; ++l)
                dof_info[i].component_to_base_index.push_back(c);
              dof_info[i].start_components[c + 1] =
                dof_info[i].start_components[c] + dof_info[i].n_components[c];
            }
          dof_info[i].dofs_per_cell.push_back(
            dof_handler[i]->get_fe().dofs_per_cell);

          // if indices are not initialized, the cell_level_index might not be
          // divisible by the vectorization length. But it must be for
          // mapping_info...
          while (cell_level_index.size() %
                   VectorizedArray<Number>::n_array_elements !=
                 0)
            cell_level_index.push_back(cell_level_index.back());
        }
    }

  // Evaluates transformations from unit to real cell, Jacobian determinants,
  // quadrature points in real space, based on the ordering of the cells
  // determined in @p extract_local_to_global_indices. The algorithm assumes
  // that the active FE index for the transformations is given the active FE
  // index in the zeroth DoFHandler. TODO: how do things look like in the more
  // general case?
  if (additional_data.initialize_mapping == true)
    {
      std::vector<unsigned int> dummy;
      mapping_info.initialize(
        dof_handler[0]->get_triangulation(),
        cell_level_index,
        face_info,
        dummy,
        mapping,
        quad,
        additional_data.mapping_update_flags,
        additional_data.mapping_update_flags_boundary_faces,
        additional_data.mapping_update_flags_inner_faces,
        additional_data.mapping_update_flags_faces_by_cells);

      mapping_is_initialized = true;
    }
}



template <int dim, typename Number>
void
MatrixFree<dim, Number>::internal_reinit(
  const Mapping<dim> &                                    mapping,
  const std::vector<const hp::DoFHandler<dim> *> &        dof_handler,
  const std::vector<const ConstraintMatrix *> &           constraint,
  const std::vector<IndexSet> &                           locally_owned_set,
  const std::vector<hp::QCollection<1>> &                 quad,
  const typename MatrixFree<dim, Number>::AdditionalData &additional_data)
{
  // Reads out the FE information and stores the shape function values,
  // gradients and Hessians for quadrature points.
  {
    unsigned int n_components = 0;
    for (unsigned int no = 0; no < dof_handler.size(); ++no)
      n_components += dof_handler[no]->get_fe()[0].n_base_elements();
    const unsigned int n_quad             = quad.size();
    unsigned int       n_fe_in_collection = 0;
    for (unsigned int i = 0; i < n_components; ++i)
      n_fe_in_collection = std::max(n_fe_in_collection,
                                    dof_handler[i]->get_fe_collection().size());
    unsigned int n_quad_in_collection = 0;
    for (unsigned int q = 0; q < n_quad; ++q)
      n_quad_in_collection = std::max(n_quad_in_collection, quad[q].size());
    shape_info.reinit(TableIndices<4>(
      n_components, n_quad, n_fe_in_collection, n_quad_in_collection));
    for (unsigned int no = 0, c = 0; no < dof_handler.size(); no++)
      for (unsigned int b = 0; b < dof_handler[no]->get_fe(0).n_base_elements();
           ++b, ++c)
        for (unsigned int fe_no = 0;
             fe_no < dof_handler[no]->get_fe_collection().size();
             ++fe_no)
          for (unsigned int nq = 0; nq < n_quad; nq++)
            for (unsigned int q_no = 0; q_no < quad[nq].size(); ++q_no)
              shape_info(c, nq, fe_no, q_no)
                .reinit(quad[nq][q_no], dof_handler[no]->get_fe(fe_no), b);
  }

  if (additional_data.initialize_indices == true)
    {
      clear();
      Assert(dof_handler.size() > 0, ExcMessage("No DoFHandler is given."));
      AssertDimension(dof_handler.size(), constraint.size());
      AssertDimension(dof_handler.size(), locally_owned_set.size());

      // set variables that are independent of FE
      if (Utilities::MPI::job_supports_mpi() == true)
        {
          const parallel::Triangulation<dim> *dist_tria =
            dynamic_cast<const parallel::Triangulation<dim> *>(
              &(dof_handler[0]->get_triangulation()));
          task_info.communicator = dist_tria != nullptr ?
                                     dist_tria->get_communicator() :
                                     MPI_COMM_SELF;
          task_info.my_pid =
            Utilities::MPI::this_mpi_process(task_info.communicator);
          task_info.n_procs =
            Utilities::MPI::n_mpi_processes(task_info.communicator);
        }
      else
        {
          task_info.communicator = MPI_COMM_SELF;
          task_info.my_pid       = 0;
          task_info.n_procs      = 1;
        }

      initialize_dof_handlers(dof_handler, additional_data);
      for (unsigned int no = 0; no < dof_handler.size(); ++no)
        {
          dof_info[no].store_plain_indices =
            additional_data.store_plain_indices;
          dof_info[no].global_base_element_offset =
            no > 0 ? dof_info[no - 1].global_base_element_offset +
                       dof_handler[no - 1]->get_fe()[0].n_base_elements() :
                     0;
        }

        // initialize the basic multithreading information that needs to be
        // passed to the DoFInfo structure
#ifdef DEAL_II_WITH_THREADS
      if (additional_data.tasks_parallel_scheme != AdditionalData::none &&
          MultithreadInfo::n_threads() > 1)
        {
          task_info.scheme =
            internal::MatrixFreeFunctions::TaskInfo::TasksParallelScheme(
              static_cast<int>(additional_data.tasks_parallel_scheme));
          task_info.block_size = additional_data.tasks_block_size;
        }
      else
#endif
        task_info.scheme = internal::MatrixFreeFunctions::TaskInfo::none;

      // set dof_indices together with constraint_indicator and
      // constraint_pool_data. It also reorders the way cells are gone through
      // (to separate cells with overlap to other processors from others
      // without).
      initialize_indices(constraint, locally_owned_set, additional_data);
    }

  // initialize bare structures
  else if (dof_info.size() != dof_handler.size())
    {
      initialize_dof_handlers(dof_handler, additional_data);
      std::vector<unsigned int>  dummy;
      std::vector<unsigned char> dummy2;
      task_info.collect_boundary_cells(
        cell_level_index.size(),
        cell_level_index.size(),
        VectorizedArray<Number>::n_array_elements,
        dummy);
      task_info.create_blocks_serial(
        dummy, dummy, 1, dummy, false, dummy, dummy2);
      for (unsigned int i = 0; i < dof_info.size(); ++i)
        {
          Assert(dof_handler[i]->get_fe_collection().size() == 1,
                 ExcNotImplemented());
          dof_info[i].dimension = dim;
          dof_info[i].n_base_elements =
            dof_handler[i]->get_fe(0).n_base_elements();
          dof_info[i].n_components.resize(dof_info[i].n_base_elements);
          dof_info[i].start_components.resize(dof_info[i].n_base_elements + 1);
          for (unsigned int c = 0; c < dof_info[i].n_base_elements; ++c)
            {
              dof_info[i].n_components[c] =
                dof_handler[i]->get_fe(0).element_multiplicity(c);
              for (unsigned int l = 0; l < dof_info[i].n_components[c]; ++l)
                dof_info[i].component_to_base_index.push_back(c);
              dof_info[i].start_components[c + 1] =
                dof_info[i].start_components[c] + dof_info[i].n_components[c];
            }
          dof_info[i].dofs_per_cell.push_back(
            dof_handler[i]->get_fe(0).dofs_per_cell);

          // if indices are not initialized, the cell_level_index might not be
          // divisible by the vectorization length. But it must be for
          // mapping_info...
          while (cell_level_index.size() %
                   VectorizedArray<Number>::n_array_elements !=
                 0)
            cell_level_index.push_back(cell_level_index.back());
        }
    }

  // Evaluates transformations from unit to real cell, Jacobian determinants,
  // quadrature points in real space, based on the ordering of the cells
  // determined in @p extract_local_to_global_indices.
  if (additional_data.initialize_mapping == true)
    {
      mapping_info.initialize(
        dof_handler[0]->get_triangulation(),
        cell_level_index,
        face_info,
        dof_info[0].cell_active_fe_index,
        mapping,
        quad,
        additional_data.mapping_update_flags,
        additional_data.mapping_update_flags_boundary_faces,
        additional_data.mapping_update_flags_inner_faces,
        additional_data.mapping_update_flags_faces_by_cells);

      mapping_is_initialized = true;
    }
}


template <int dim, typename Number>
template <int spacedim>
bool
MatrixFree<dim, Number>::is_supported(const FiniteElement<dim, spacedim> &fe)
{
  if (dim != spacedim)
    return false;

  // first check for degree, number of base_elemnt and number of its components
  if (fe.degree == 0 || fe.n_base_elements() != 1)
    return false;

  const FiniteElement<dim, spacedim> *fe_ptr = &(fe.base_element(0));
  if (fe_ptr->n_components() != 1)
    return false;

  // then check of the base element is supported
  if (dynamic_cast<const FE_Poly<TensorProductPolynomials<dim>, dim, spacedim>
                     *>(fe_ptr) != nullptr)
    return true;
  if (dynamic_cast<const FE_Poly<
        TensorProductPolynomials<dim, Polynomials::PiecewisePolynomial<double>>,
        dim,
        spacedim> *>(fe_ptr) != nullptr)
    return true;
  if (dynamic_cast<const FE_DGP<dim, spacedim> *>(fe_ptr) != nullptr)
    return true;
  if (dynamic_cast<const FE_Q_DG0<dim, spacedim> *>(fe_ptr) != nullptr)
    return true;

  // if the base element is not in the above list it is not supported
  return false;
}



namespace internal
{
  namespace MatrixFreeFunctions
  {
    // steps through all children and adds the active cells recursively
    template <typename InIterator>
    void
    resolve_cell(const InIterator &                                  cell,
                 std::vector<std::pair<unsigned int, unsigned int>> &cell_its,
                 const unsigned int subdomain_id)
    {
      if (cell->has_children())
        for (unsigned int child = 0; child < cell->n_children(); ++child)
          resolve_cell(cell->child(child), cell_its, subdomain_id);
      else if (subdomain_id == numbers::invalid_subdomain_id ||
               cell->subdomain_id() == subdomain_id)
        {
          Assert(cell->active(), ExcInternalError());
          cell_its.emplace_back(cell->level(), cell->index());
        }
    }
  } // namespace MatrixFreeFunctions
} // namespace internal



template <int dim, typename Number>
void
MatrixFree<dim, Number>::initialize_dof_handlers(
  const std::vector<const DoFHandler<dim> *> &dof_handler,
  const AdditionalData &                      additional_data)
{
  cell_level_index.clear();
  dof_handlers.active_dof_handler = DoFHandlers::usual;
  dof_handlers.level              = additional_data.level_mg_handler;
  dof_handlers.n_dof_handlers     = dof_handler.size();
  dof_handlers.dof_handler.resize(dof_handlers.n_dof_handlers);
  for (unsigned int no = 0; no < dof_handlers.n_dof_handlers; ++no)
    dof_handlers.dof_handler[no] = dof_handler[no];

  dof_info.resize(dof_handlers.n_dof_handlers);
  for (unsigned int no = 0; no < dof_handlers.n_dof_handlers; ++no)
    dof_info[no].vectorization_length =
      VectorizedArray<Number>::n_array_elements;

  // Go through cells on zeroth level and then successively step down into
  // children. This gives a z-ordering of the cells, which is beneficial when
  // setting up neighboring relations between cells for thread parallelization
  const unsigned int n_mpi_procs = task_info.n_procs;
  const unsigned int my_pid      = task_info.my_pid;

  const Triangulation<dim> &tria =
    dof_handlers.dof_handler[0]->get_triangulation();
  const unsigned int level = additional_data.level_mg_handler;
  if (level == numbers::invalid_unsigned_int)
    {
      if (n_mpi_procs == 1)
        cell_level_index.reserve(tria.n_active_cells());
      typename Triangulation<dim>::cell_iterator cell     = tria.begin(0),
                                                 end_cell = tria.end(0);
      // For serial Triangulations always take all cells
      const unsigned int subdomain_id =
        (dynamic_cast<const parallel::Triangulation<dim> *>(
           &dof_handler[0]->get_triangulation()) != nullptr) ?
          my_pid :
          numbers::invalid_subdomain_id;
      for (; cell != end_cell; ++cell)
        internal::MatrixFreeFunctions::resolve_cell(
          cell, cell_level_index, subdomain_id);

      Assert(n_mpi_procs > 1 ||
               cell_level_index.size() == tria.n_active_cells(),
             ExcInternalError());
    }
  else
    {
      AssertIndexRange(level, tria.n_global_levels());
      if (level < tria.n_levels())
        {
          cell_level_index.reserve(tria.n_cells(level));
          typename Triangulation<dim>::cell_iterator cell = tria.begin(level),
                                                     end_cell = tria.end(level);
          for (; cell != end_cell; ++cell)
            if (cell->level_subdomain_id() == my_pid)
              cell_level_index.emplace_back(cell->level(), cell->index());
        }
    }

  // All these are cells local to this processor. Therefore, set
  // cell_level_index_end_local to the size of cell_level_index.
  cell_level_index_end_local = cell_level_index.size();
}



template <int dim, typename Number>
void
MatrixFree<dim, Number>::initialize_dof_handlers(
  const std::vector<const hp::DoFHandler<dim> *> &dof_handler,
  const AdditionalData &                          additional_data)
{
  cell_level_index.clear();
  dof_handlers.active_dof_handler = DoFHandlers::hp;
  dof_handlers.level              = additional_data.level_mg_handler;
  Assert(dof_handlers.level == numbers::invalid_unsigned_int,
         ExcNotImplemented());
  dof_handlers.n_dof_handlers = dof_handler.size();
  dof_handlers.hp_dof_handler.resize(dof_handlers.n_dof_handlers);
  for (unsigned int no = 0; no < dof_handlers.n_dof_handlers; ++no)
    dof_handlers.hp_dof_handler[no] = dof_handler[no];

  dof_info.resize(dof_handlers.n_dof_handlers);
  for (unsigned int no = 0; no < dof_handlers.n_dof_handlers; ++no)
    dof_info[no].vectorization_length =
      VectorizedArray<Number>::n_array_elements;

  // go through cells on zeroth level and then successively step down into
  // children. This gives a z-ordering of the cells, which is beneficial when
  // setting up neighboring relations between cells for thread parallelization
  const unsigned int n_mpi_procs = task_info.n_procs;
  const unsigned int my_pid      = task_info.my_pid;

  // if we have no level given, use the same as for the standard DoFHandler,
  // otherwise we must loop through the respective level
  const Triangulation<dim> &tria = dof_handler[0]->get_triangulation();

  if (n_mpi_procs == 1)
    {
      cell_level_index.reserve(tria.n_active_cells());
    }
  typename hp::DoFHandler<dim>::cell_iterator cell = dof_handler[0]->begin(0),
                                              end_cell = dof_handler[0]->end(0);
  // For serial Triangulations always take all cells
  const unsigned int subdomain_id =
    (dynamic_cast<const parallel::Triangulation<dim> *>(
       &dof_handler[0]->get_triangulation()) != nullptr) ?
      my_pid :
      numbers::invalid_subdomain_id;
  for (; cell != end_cell; ++cell)
    internal::MatrixFreeFunctions::resolve_cell(
      cell, cell_level_index, subdomain_id);

  Assert(n_mpi_procs > 1 || cell_level_index.size() == tria.n_active_cells(),
         ExcInternalError());

  // All these are cells local to this processor. Therefore, set
  // cell_level_index_end_local to the size of cell_level_index.
  cell_level_index_end_local = cell_level_index.size();
}



template <int dim, typename Number>
void
MatrixFree<dim, Number>::initialize_indices(
  const std::vector<const ConstraintMatrix *> &constraint,
  const std::vector<IndexSet> &                locally_owned_set,
  const AdditionalData &                       additional_data)
{
  // insert possible ghost cells and construct face topology
  const bool do_face_integrals =
    (additional_data.mapping_update_flags_inner_faces |
     additional_data.mapping_update_flags_boundary_faces) != update_default;
  internal::MatrixFreeFunctions::FaceSetup<dim> face_setup;
  if (do_face_integrals)
    face_setup.initialize(dof_handlers.active_dof_handler ==
                              DoFHandlers::usual ?
                            dof_handlers.dof_handler[0]->get_triangulation() :
                            dof_handlers.hp_dof_handler[0]->get_triangulation(),
                          additional_data,
                          cell_level_index);

  const unsigned int n_fe           = dof_handlers.n_dof_handlers;
  const unsigned int n_active_cells = cell_level_index.size();

  AssertDimension(n_active_cells, cell_level_index.size());
  AssertDimension(n_fe, locally_owned_set.size());
  AssertDimension(n_fe, constraint.size());

  std::vector<types::global_dof_index>                local_dof_indices;
  std::vector<std::vector<std::vector<unsigned int>>> lexicographic(n_fe);

  internal::MatrixFreeFunctions::ConstraintValues<double> constraint_values;

  bool cell_categorization_enabled =
    !additional_data.cell_vectorization_category.empty();

  for (unsigned int no = 0; no < n_fe; ++no)
    {
      std::vector<const FiniteElement<dim> *> fes;
      if (dof_handlers.active_dof_handler == DoFHandlers::hp)
        {
          const hp::DoFHandler<dim> *  hpdof = dof_handlers.hp_dof_handler[no];
          const hp::FECollection<dim> &fe    = hpdof->get_fe_collection();
          for (unsigned int f = 0; f < fe.size(); ++f)
            fes.push_back(&fe[f]);

          if (fe.size() > 1)
            dof_info[no].cell_active_fe_index.resize(
              n_active_cells, numbers::invalid_unsigned_int);

          Assert(additional_data.cell_vectorization_category.empty(),
                 ExcNotImplemented());
        }
      else
        {
          const DoFHandler<dim> *dofh = &*dof_handlers.dof_handler[no];
          fes.push_back(&dofh->get_fe());
          if (cell_categorization_enabled == true)
            dof_info[no].cell_active_fe_index.resize(
              n_active_cells, numbers::invalid_unsigned_int);
        }
      lexicographic[no].resize(fes.size());

      dof_info[no].fe_index_conversion.resize(fes.size());
      dof_info[no].max_fe_index = fes.size();

      dof_info[no].component_dof_indices_offset.clear();
      dof_info[no].component_dof_indices_offset.resize(fes.size());
      for (unsigned int fe_index = 0; fe_index < fes.size(); ++fe_index)
        {
          const FiniteElement<dim> &fe = *fes[fe_index];
          // cache number of finite elements and dofs_per_cell
          dof_info[no].dofs_per_cell.push_back(fe.dofs_per_cell);
          dof_info[no].dofs_per_face.push_back(fe.dofs_per_face);
          dof_info[no].dimension       = dim;
          dof_info[no].n_base_elements = fe.n_base_elements();
          dof_info[no].n_components.resize(dof_info[no].n_base_elements);
          dof_info[no].start_components.resize(dof_info[no].n_base_elements +
                                               1);
          dof_info[no].component_to_base_index.clear();
          dof_info[no].component_dof_indices_offset[fe_index].push_back(0);
          dof_info[no].fe_index_conversion[fe_index].clear();
          for (unsigned int c = 0; c < dof_info[no].n_base_elements; ++c)
            {
              dof_info[no].n_components[c] = fe.element_multiplicity(c);
              for (unsigned int l = 0; l < dof_info[no].n_components[c]; ++l)
                {
                  dof_info[no].component_to_base_index.push_back(c);
                  dof_info[no].component_dof_indices_offset[fe_index].push_back(
                    dof_info[no].component_dof_indices_offset[fe_index].back() +
                    fe.base_element(c).dofs_per_cell);
                  dof_info[no].fe_index_conversion[fe_index].push_back(
                    fe.base_element(c).degree);
                }
              dof_info[no].start_components[c + 1] =
                dof_info[no].start_components[c] + dof_info[no].n_components[c];
              lexicographic[no][fe_index].insert(
                lexicographic[no][fe_index].end(),
                shape_info(
                  dof_info[no].global_base_element_offset + c, 0, fe_index, 0)
                  .lexicographic_numbering.begin(),
                shape_info(
                  dof_info[no].global_base_element_offset + c, 0, fe_index, 0)
                  .lexicographic_numbering.end());
            }

          AssertDimension(lexicographic[no][fe_index].size(),
                          dof_info[no].dofs_per_cell[fe_index]);
          AssertDimension(
            dof_info[no].component_dof_indices_offset[fe_index].size() - 1,
            dof_info[no].start_components.back());
          AssertDimension(
            dof_info[no].component_dof_indices_offset[fe_index].back(),
            dof_info[no].dofs_per_cell[fe_index]);
        }

      // set locally owned range for each component
      Assert(locally_owned_set[no].is_contiguous(), ExcNotImplemented());
      dof_info[no].vector_partitioner.reset(new Utilities::MPI::Partitioner(
        locally_owned_set[no], task_info.communicator));

      // initialize the arrays for indices
      const unsigned int n_components_total =
        dof_info[no].start_components.back();
      dof_info[no].row_starts.resize(n_active_cells * n_components_total + 1);
      dof_info[no].row_starts[0].first  = 0;
      dof_info[no].row_starts[0].second = 0;
      dof_info[no].dof_indices.reserve(
        (n_active_cells * dof_info[no].dofs_per_cell[0] * 3) / 2);

      // cache the constrained indices for use in matrix-vector products and
      // the like
      const types::global_dof_index
        start_index = dof_info[no].vector_partitioner->local_range().first,
        end_index   = dof_info[no].vector_partitioner->local_range().second;
      for (types::global_dof_index i = start_index; i < end_index; ++i)
        if (constraint[no]->is_constrained(i) == true)
          dof_info[no].constrained_dofs.push_back(
            static_cast<unsigned int>(i - start_index));
    }

  // extract all the global indices associated with the computation, and form
  // the ghost indices
  std::vector<unsigned int> subdomain_boundary_cells;
  for (unsigned int counter = 0; counter < n_active_cells; ++counter)
    {
      bool cell_at_subdomain_boundary =
        (face_setup.at_processor_boundary.size() > counter &&
         face_setup.at_processor_boundary[counter]) ||
        (additional_data.overlap_communication_computation == false &&
         task_info.n_procs > 1);

      for (unsigned int no = 0; no < n_fe; ++no)
        {
          // read indices from standard DoFHandler in the usual way
          if (dof_handlers.active_dof_handler == DoFHandlers::usual &&
              dof_handlers.level == numbers::invalid_unsigned_int)
            {
              const DoFHandler<dim> *dofh = &*dof_handlers.dof_handler[no];
              typename DoFHandler<dim>::active_cell_iterator cell_it(
                &dofh->get_triangulation(),
                cell_level_index[counter].first,
                cell_level_index[counter].second,
                dofh);
              local_dof_indices.resize(dof_info[no].dofs_per_cell[0]);
              cell_it->get_dof_indices(local_dof_indices);
              dof_info[no].read_dof_indices(local_dof_indices,
                                            lexicographic[no][0],
                                            *constraint[no],
                                            counter,
                                            constraint_values,
                                            cell_at_subdomain_boundary);
              if (cell_categorization_enabled)
                {
                  AssertIndexRange(
                    cell_it->active_cell_index(),
                    additional_data.cell_vectorization_category.size());
                  dof_info[no].cell_active_fe_index[counter] =
                    additional_data.cell_vectorization_category
                      [cell_it->active_cell_index()];
                }
            }
          // we are requested to use a multigrid level
          else if (dof_handlers.active_dof_handler == DoFHandlers::usual &&
                   dof_handlers.level != numbers::invalid_unsigned_int)
            {
              const DoFHandler<dim> *dofh = dof_handlers.dof_handler[no];
              AssertIndexRange(dof_handlers.level,
                               dofh->get_triangulation().n_levels());
              typename DoFHandler<dim>::cell_iterator cell_it(
                &dofh->get_triangulation(),
                cell_level_index[counter].first,
                cell_level_index[counter].second,
                dofh);
              local_dof_indices.resize(dof_info[no].dofs_per_cell[0]);
              cell_it->get_mg_dof_indices(local_dof_indices);
              dof_info[no].read_dof_indices(local_dof_indices,
                                            lexicographic[no][0],
                                            *constraint[no],
                                            counter,
                                            constraint_values,
                                            cell_at_subdomain_boundary);
              if (cell_categorization_enabled)
                {
                  AssertIndexRange(
                    cell_it->active_cell_index(),
                    additional_data.cell_vectorization_category.size());
                  dof_info[no].cell_active_fe_index[counter] =
                    additional_data.cell_vectorization_category
                      [cell_level_index[counter].second];
                }
            }
          // hp case where we need to decode the FE index and similar
          else if (dof_handlers.active_dof_handler == DoFHandlers::hp)
            {
              const hp::DoFHandler<dim> *dofh = dof_handlers.hp_dof_handler[no];
              typename hp::DoFHandler<dim>::active_cell_iterator cell_it(
                &dofh->get_triangulation(),
                cell_level_index[counter].first,
                cell_level_index[counter].second,
                dofh);
              if (dofh->get_fe_collection().size() > 1)
                dof_info[no].cell_active_fe_index[counter] =
                  cell_it->active_fe_index();
              local_dof_indices.resize(cell_it->get_fe().dofs_per_cell);
              cell_it->get_dof_indices(local_dof_indices);
              dof_info[no].read_dof_indices(
                local_dof_indices,
                lexicographic[no][cell_it->active_fe_index()],
                *constraint[no],
                counter,
                constraint_values,
                cell_at_subdomain_boundary);
            }
          else
            {
              Assert(false, ExcNotImplemented());
            }
        }

      // if we found dofs on some FE component that belong to other
      // processors, the cell is added to the boundary cells.
      if (cell_at_subdomain_boundary == true &&
          counter < cell_level_index_end_local)
        subdomain_boundary_cells.push_back(counter);
    }

  const unsigned int vectorization_length =
    VectorizedArray<Number>::n_array_elements;
  task_info.collect_boundary_cells(cell_level_index_end_local,
                                   n_active_cells,
                                   vectorization_length,
                                   subdomain_boundary_cells);

  // Finalize the creation of the ghost indices
  {
    std::vector<unsigned int> cells_with_ghosts(subdomain_boundary_cells);
    for (unsigned int c = cell_level_index_end_local; c < n_active_cells; ++c)
      cells_with_ghosts.push_back(c);
    for (unsigned int no = 0; no < n_fe; ++no)
      {
        if (do_face_integrals &&
            additional_data.level_mg_handler != numbers::invalid_unsigned_int)
          {
            // in case of adaptivity, go through the cells on the next finer
            // level and check whether we need to get read access to some of
            // those entries for the mg flux matrices
            const DoFHandler<dim> &dof_handler = *dof_handlers.dof_handler[no];
            std::vector<types::global_dof_index> dof_indices;
            if (additional_data.level_mg_handler + 1 <
                dof_handler.get_triangulation().n_global_levels())
              for (typename DoFHandler<dim>::cell_iterator cell =
                     dof_handler.begin(additional_data.level_mg_handler + 1);
                   cell !=
                   dof_handler.end(additional_data.level_mg_handler + 1);
                   ++cell)
                if (cell->level_subdomain_id() == task_info.my_pid)
                  for (unsigned int f = 0;
                       f < GeometryInfo<dim>::faces_per_cell;
                       ++f)
                    if ((cell->at_boundary(f) == false ||
                         cell->has_periodic_neighbor(f) == true) &&
                        cell->level() >
                          cell->neighbor_or_periodic_neighbor(f)->level() &&
                        cell->neighbor_or_periodic_neighbor(f)
                            ->level_subdomain_id() != task_info.my_pid)
                      {
                        dof_indices.resize(
                          cell->neighbor_or_periodic_neighbor(f)
                            ->get_fe()
                            .dofs_per_cell);
                        cell->neighbor_or_periodic_neighbor(f)
                          ->get_mg_dof_indices(dof_indices);
                        for (unsigned int i = 0; i < dof_indices.size(); ++i)
                          dof_info[no].ghost_dofs.push_back(dof_indices[i]);
                      }
          }
        dof_info[no].assign_ghosts(cells_with_ghosts);
      }
  }

  std::vector<unsigned int>  renumbering;
  std::vector<unsigned char> irregular_cells;
  if (task_info.scheme == internal::MatrixFreeFunctions::TaskInfo::none)
    {
      const bool strict_categories =
        additional_data.cell_vectorization_categories_strict ||
        dof_handlers.active_dof_handler == DoFHandlers::hp;
      unsigned int dofs_per_cell = 0;
      for (unsigned int no = 0; no < dof_info.size(); ++no)
        dofs_per_cell = std::max(dofs_per_cell, dof_info[no].dofs_per_cell[0]);
      task_info.create_blocks_serial(subdomain_boundary_cells,
                                     face_setup.cells_close_to_boundary,
                                     dofs_per_cell,
                                     dof_info[0].cell_active_fe_index,
                                     strict_categories,
                                     renumbering,
                                     irregular_cells);
    }
  else
    {
      // For strategy with blocking before partitioning: reorganize the indices
      // in order to overlap communication in MPI with computations: Place all
      // cells with ghost indices into one chunk. Also reorder cells so that we
      // can parallelize by threads
      Assert(additional_data.cell_vectorization_category.empty(),
             ExcNotImplemented());
      task_info.initial_setup_blocks_tasks(
        subdomain_boundary_cells, renumbering, irregular_cells);
      task_info.guess_block_size(dof_info[0].dofs_per_cell[0]);

      unsigned int n_macro_cells_before =
        *(task_info.cell_partition_data.end() - 2);
      unsigned int n_ghost_slots =
        *(task_info.cell_partition_data.end() - 1) - n_macro_cells_before;

      unsigned int start_nonboundary = numbers::invalid_unsigned_int;
      if (task_info.scheme ==
            internal::MatrixFreeFunctions::TaskInfo::partition_color ||
          task_info.scheme == internal::MatrixFreeFunctions::TaskInfo::color)
        {
          // set up partitions. if we just use coloring without partitions, do
          // nothing here, assume all cells to belong to the zero partition
          // (that we otherwise use for MPI boundary cells)
          if (task_info.scheme ==
              internal::MatrixFreeFunctions::TaskInfo::color)
            {
              start_nonboundary =
                task_info.n_procs > 1 ?
                  std::min(((task_info.cell_partition_data[2] -
                             task_info.cell_partition_data[1] +
                             task_info.block_size - 1) /
                            task_info.block_size) *
                             task_info.block_size,
                           task_info.cell_partition_data[3]) :
                  0;
            }
          else
            {
              if (task_info.n_procs > 1)
                {
                  task_info.cell_partition_data[1] = 0;
                  task_info.cell_partition_data[2] =
                    task_info.cell_partition_data[3];
                }
              start_nonboundary = task_info.cell_partition_data.back();
            }

          if (dof_handlers.active_dof_handler == DoFHandlers::hp)
            {
              irregular_cells.resize(0);
              irregular_cells.resize(task_info.cell_partition_data.back() +
                                     2 * dof_info[0].max_fe_index);
              std::vector<std::vector<unsigned int>> renumbering_fe_index;
              renumbering_fe_index.resize(dof_info[0].max_fe_index);
              unsigned int counter;
              n_macro_cells_before = 0;
              for (counter = 0;
                   counter < std::min(start_nonboundary * vectorization_length,
                                      task_info.n_active_cells);
                   counter++)
                {
                  AssertIndexRange(counter, renumbering.size());
                  AssertIndexRange(renumbering[counter],
                                   dof_info[0].cell_active_fe_index.size());
                  renumbering_fe_index
                    [dof_info[0].cell_active_fe_index[renumbering[counter]]]
                      .push_back(renumbering[counter]);
                }
              counter = 0;
              for (unsigned int j = 0; j < dof_info[0].max_fe_index; j++)
                {
                  for (unsigned int jj = 0; jj < renumbering_fe_index[j].size();
                       jj++)
                    renumbering[counter++] = renumbering_fe_index[j][jj];
                  irregular_cells[renumbering_fe_index[j].size() /
                                    vectorization_length +
                                  n_macro_cells_before] =
                    renumbering_fe_index[j].size() % vectorization_length;
                  n_macro_cells_before += (renumbering_fe_index[j].size() +
                                           vectorization_length - 1) /
                                          vectorization_length;
                  renumbering_fe_index[j].resize(0);
                }

              for (counter = start_nonboundary * vectorization_length;
                   counter < task_info.n_active_cells;
                   counter++)
                {
                  renumbering_fe_index
                    [dof_info[0].cell_active_fe_index.empty() ?
                       0 :
                       dof_info[0].cell_active_fe_index[renumbering[counter]]]
                      .push_back(renumbering[counter]);
                }
              counter = start_nonboundary * vectorization_length;
              for (unsigned int j = 0; j < dof_info[0].max_fe_index; j++)
                {
                  for (unsigned int jj = 0; jj < renumbering_fe_index[j].size();
                       jj++)
                    renumbering[counter++] = renumbering_fe_index[j][jj];
                  irregular_cells[renumbering_fe_index[j].size() /
                                    vectorization_length +
                                  n_macro_cells_before] =
                    renumbering_fe_index[j].size() % vectorization_length;
                  n_macro_cells_before += (renumbering_fe_index[j].size() +
                                           vectorization_length - 1) /
                                          vectorization_length;
                }
              AssertIndexRange(n_macro_cells_before,
                               task_info.cell_partition_data.back() +
                                 2 * dof_info[0].max_fe_index + 1);
              irregular_cells.resize(n_macro_cells_before + n_ghost_slots);
              *(task_info.cell_partition_data.end() - 2) = n_macro_cells_before;
              *(task_info.cell_partition_data.end() - 1) =
                n_macro_cells_before + n_ghost_slots;
            }
        }

      task_info.n_blocks =
        (n_macro_cells() + task_info.block_size - 1) / task_info.block_size;

      DynamicSparsityPattern connectivity;
      connectivity.reinit(task_info.n_active_cells, task_info.n_active_cells);
      if ((additional_data.mapping_update_flags_inner_faces |
           additional_data.mapping_update_flags_boundary_faces) !=
          update_default)
        make_connectivity_graph_faces(connectivity);
      if (task_info.n_active_cells > 0)
        dof_info[0].make_connectivity_graph(
          task_info, renumbering, connectivity);

      task_info.make_thread_graph(dof_info[0].cell_active_fe_index,
                                  connectivity,
                                  renumbering,
                                  irregular_cells,
                                  dof_handlers.active_dof_handler ==
                                    DoFHandlers::hp);

      Assert(irregular_cells.size() >= task_info.cell_partition_data.back(),
             ExcInternalError());

      irregular_cells.resize(task_info.cell_partition_data.back() +
                             n_ghost_slots);
      if (n_ghost_slots > 0)
        {
          for (unsigned int i = task_info.cell_partition_data.back();
               i < task_info.cell_partition_data.back() + n_ghost_slots - 1;
               ++i)
            irregular_cells[i] = 0;
          irregular_cells.back() =
            task_info.n_ghost_cells % vectorization_length;
        }

      {
        unsigned int n_cells = 0;
        for (unsigned int i = 0; i < task_info.cell_partition_data.back(); ++i)
          n_cells +=
            irregular_cells[i] > 0 ? irregular_cells[i] : vectorization_length;
        AssertDimension(n_cells, task_info.n_active_cells);
        n_cells = 0;
        for (unsigned int i = task_info.cell_partition_data.back();
             i < n_ghost_slots + task_info.cell_partition_data.back();
             ++i)
          n_cells +=
            irregular_cells[i] > 0 ? irregular_cells[i] : vectorization_length;
        AssertDimension(n_cells, task_info.n_ghost_cells);
      }

      task_info.cell_partition_data.push_back(
        task_info.cell_partition_data.back() + n_ghost_slots);
    }

    // Finally perform the renumbering. We also want to group several cells
    // together to a batch of cells for SIMD (vectorized) execution (where the
    // arithmetic operations of several cells will then be done simultaneously).
#ifdef DEBUG
  {
    AssertDimension(renumbering.size(),
                    task_info.n_active_cells + task_info.n_ghost_cells);
    std::vector<unsigned int> sorted_renumbering(renumbering);
    std::sort(sorted_renumbering.begin(), sorted_renumbering.end());
    for (unsigned int i = 0; i < sorted_renumbering.size(); ++i)
      Assert(sorted_renumbering[i] == i, ExcInternalError());
  }
#endif
  {
    std::vector<std::pair<unsigned int, unsigned int>> cell_level_index_old;
    cell_level_index.swap(cell_level_index_old);
    cell_level_index.reserve(task_info.cell_partition_data.back() *
                             vectorization_length);
    unsigned int position_cell = 0;
    for (unsigned int i = 0; i < task_info.cell_partition_data.back(); ++i)
      {
        unsigned int n_comp =
          (irregular_cells[i] > 0) ? irregular_cells[i] : vectorization_length;
        for (unsigned int j = 0; j < n_comp; ++j)
          cell_level_index.push_back(
            cell_level_index_old[renumbering[position_cell + j]]);

        // generate a cell and level index also when we have not filled up
        // vectorization_length cells. This is needed for MappingInfo when the
        // transformation data is initialized. We just set the value to the
        // last valid cell in that case.
        for (unsigned int j = n_comp; j < vectorization_length; ++j)
          cell_level_index.push_back(
            cell_level_index_old[renumbering[position_cell + n_comp - 1]]);
        position_cell += n_comp;
      }
    AssertDimension(position_cell,
                    task_info.n_active_cells + task_info.n_ghost_cells);
    AssertDimension(cell_level_index.size(),
                    task_info.cell_partition_data.back() *
                      vectorization_length);
  }

  // set constraint pool from the std::map and reorder the indices
  typename std::map<std::vector<double>,
                    types::global_dof_index,
                    internal::MatrixFreeFunctions::FPArrayComparator<double>>::
    iterator it  = constraint_values.constraints.begin(),
             end = constraint_values.constraints.end();
  std::vector<const std::vector<double> *> constraints(
    constraint_values.constraints.size());
  unsigned int length = 0;
  for (; it != end; ++it)
    {
      AssertIndexRange(it->second, constraints.size());
      constraints[it->second] = &it->first;
      length += it->first.size();
    }
  constraint_pool_data.clear();
  constraint_pool_data.reserve(length);
  constraint_pool_row_index.reserve(constraint_values.constraints.size() + 1);
  constraint_pool_row_index.resize(1, 0);
  for (unsigned int i = 0; i < constraints.size(); ++i)
    {
      Assert(constraints[i] != nullptr, ExcInternalError());
      constraint_pool_data.insert(constraint_pool_data.end(),
                                  constraints[i]->begin(),
                                  constraints[i]->end());
      constraint_pool_row_index.push_back(constraint_pool_data.size());
    }

  AssertDimension(constraint_pool_data.size(), length);
  for (unsigned int no = 0; no < n_fe; ++no)
    dof_info[no].reorder_cells(
      task_info, renumbering, constraint_pool_row_index, irregular_cells);

  // Finally resort the faces and collect several faces for vectorization
  if ((additional_data.mapping_update_flags_inner_faces |
       additional_data.mapping_update_flags_boundary_faces) != update_default)
    {
      face_setup.generate_faces(
        dof_handlers.active_dof_handler == DoFHandlers::usual ?
          dof_handlers.dof_handler[0]->get_triangulation() :
          dof_handlers.hp_dof_handler[0]->get_triangulation(),
        cell_level_index,
        task_info);
      face_info.faces.clear();

      std::vector<bool> hard_vectorization_boundary(
        task_info.face_partition_data.size(), false);
      if (task_info.scheme == internal::MatrixFreeFunctions::TaskInfo::none &&
          task_info.partition_row_index[2] <
            task_info.face_partition_data.size())
        hard_vectorization_boundary[task_info.partition_row_index[2]] = true;
      else
        for (unsigned int i = 0; i < hard_vectorization_boundary.size(); ++i)
          hard_vectorization_boundary[i] = true;

      internal::MatrixFreeFunctions::collect_faces_vectorization(
        face_setup.inner_faces,
        hard_vectorization_boundary,
        task_info.face_partition_data,
        face_info.faces);
      hard_vectorization_boundary.clear();
      hard_vectorization_boundary.resize(
        task_info.boundary_partition_data.size(),
        task_info.scheme != internal::MatrixFreeFunctions::TaskInfo::none);
      internal::MatrixFreeFunctions::collect_faces_vectorization(
        face_setup.boundary_faces,
        hard_vectorization_boundary,
        task_info.boundary_partition_data,
        face_info.faces);
      hard_vectorization_boundary.clear();
      hard_vectorization_boundary.resize(
        task_info.ghost_face_partition_data.size(), false);
      internal::MatrixFreeFunctions::collect_faces_vectorization(
        face_setup.inner_ghost_faces,
        hard_vectorization_boundary,
        task_info.ghost_face_partition_data,
        face_info.faces);
      hard_vectorization_boundary.clear();
      hard_vectorization_boundary.resize(
        task_info.refinement_edge_face_partition_data.size(), false);
      internal::MatrixFreeFunctions::collect_faces_vectorization(
        face_setup.refinement_edge_faces,
        hard_vectorization_boundary,
        task_info.refinement_edge_face_partition_data,
        face_info.faces);

      cell_level_index.resize(
        cell_level_index.size() +
        VectorizedArray<Number>::n_array_elements *
          (task_info.refinement_edge_face_partition_data[1] -
           task_info.refinement_edge_face_partition_data[0]));

      for (unsigned int no = 0; no < n_fe; ++no)
        dof_info[no].compute_face_index_compression(face_info.faces);

      // build the inverse map back from the faces array to
      // cell_and_face_to_plain_faces
      face_info.cell_and_face_to_plain_faces.reinit(
        TableIndices<3>(task_info.cell_partition_data.back(),
                        GeometryInfo<dim>::faces_per_cell,
                        VectorizedArray<Number>::n_array_elements),
        true);
      face_info.cell_and_face_to_plain_faces.fill(
        numbers::invalid_unsigned_int);
      face_info.cell_and_face_boundary_id.reinit(
        TableIndices<3>(task_info.cell_partition_data.back(),
                        GeometryInfo<dim>::faces_per_cell,
                        VectorizedArray<Number>::n_array_elements),
        true);
      face_info.cell_and_face_boundary_id.fill(numbers::invalid_boundary_id);

      for (unsigned int f = 0; f < task_info.boundary_partition_data.back();
           ++f)
        for (unsigned int v = 0;
             v < VectorizedArray<Number>::n_array_elements &&
             face_info.faces[f].cells_interior[v] !=
               numbers::invalid_unsigned_int;
             ++v)
          {
            TableIndices<3> index(face_info.faces[f].cells_interior[v] /
                                    VectorizedArray<Number>::n_array_elements,
                                  face_info.faces[f].interior_face_no,
                                  face_info.faces[f].cells_interior[v] %
                                    VectorizedArray<Number>::n_array_elements);

            // Assert(cell_and_face_to_plain_faces(index) ==
            // numbers::invalid_unsigned_int,
            //       ExcInternalError("Should only visit each face once"));
            face_info.cell_and_face_to_plain_faces(index) =
              f * VectorizedArray<Number>::n_array_elements + v;
            if (face_info.faces[f].cells_exterior[v] !=
                numbers::invalid_unsigned_int)
              {
                TableIndices<3> index(
                  face_info.faces[f].cells_exterior[v] /
                    VectorizedArray<Number>::n_array_elements,
                  face_info.faces[f].exterior_face_no,
                  face_info.faces[f].cells_exterior[v] %
                    VectorizedArray<Number>::n_array_elements);
                // Assert(cell_and_face_to_plain_faces(index) ==
                // numbers::invalid_unsigned_int,
                //       ExcInternalError("Should only visit each face once"));
                face_info.cell_and_face_to_plain_faces(index) =
                  f * VectorizedArray<Number>::n_array_elements + v;
              }
            else
              face_info.cell_and_face_boundary_id(index) =
                types::boundary_id(face_info.faces[f].exterior_face_no);
          }

      // compute tighter index sets for various sets of face integrals
      for (unsigned int no = 0; no < n_fe; ++no)
        {
          const Utilities::MPI::Partitioner &part =
            *dof_info[no].vector_partitioner;

          // partitioner 0: no face integrals, simply use the indices present
          // on the cells
          unsigned int n_macro_cells_before =
            *(task_info.cell_partition_data.end() - 2);
          std::vector<types::global_dof_index> ghost_indices;
          {
            for (unsigned int cell = 0;
                 cell < VectorizedArray<Number>::n_array_elements *
                          n_macro_cells_before;
                 ++cell)
              if (cell > 0 &&
                  cell_level_index[cell] != cell_level_index[cell - 1])
                {
                  for (unsigned int i =
                         dof_info[no]
                           .row_starts[cell *
                                       dof_info[no].start_components.back()]
                           .first;
                       i < dof_info[no]
                             .row_starts[(cell + 1) *
                                         dof_info[no].start_components.back()]
                             .first;
                       ++i)
                    if (dof_info[no].dof_indices[i] > part.local_size())
                      ghost_indices.push_back(
                        part.local_to_global(dof_info[no].dof_indices[i]));
                  for (unsigned int i =
                         dof_info[no].row_starts_plain_indices[cell];
                       i < dof_info[no].row_starts_plain_indices[cell + 1];
                       ++i)
                    if (dof_info[no].plain_dof_indices[i] > part.local_size())
                      ghost_indices.push_back(part.local_to_global(
                        dof_info[no].plain_dof_indices[i]));
                }
            std::sort(ghost_indices.begin(), ghost_indices.end());
            ghost_indices.erase(
              std::unique(ghost_indices.begin(), ghost_indices.end()),
              ghost_indices.end());
            IndexSet compressed_set(part.size());
            compressed_set.add_indices(ghost_indices.begin(),
                                       ghost_indices.end());
            compressed_set.subtract_set(
              dof_info[no].vector_partitioner->locally_owned_range());
            const bool all_ghosts_equal = Utilities::MPI::min(
              (int)(compressed_set.n_elements() ==
                    dof_info[no]
                      .vector_partitioner->ghost_indices()
                      .n_elements()),
              dof_info[no].vector_partitioner->get_mpi_communicator());
            if (all_ghosts_equal)
              dof_info[no].vector_partitioner_face_variants[0] =
                dof_info[no].vector_partitioner;
            else
              {
                dof_info[no].vector_partitioner_face_variants[0].reset(
                  new Utilities::MPI::Partitioner(part.locally_owned_range(),
                                                  part.get_mpi_communicator()));
                const_cast<Utilities::MPI::Partitioner *>(
                  dof_info[no].vector_partitioner_face_variants[0].get())
                  ->set_ghost_indices(compressed_set, part.ghost_indices());
              }
          }

          // partitioner 1: values on faces
          {
            bool all_nodal = true;
            for (unsigned int c = 0; c < dof_info[no].n_base_elements; ++c)
              if (!shape_info(
                     dof_info[no].global_base_element_offset + c, 0, 0, 0)
                     .nodal_at_cell_boundaries)
                all_nodal = false;
            if (all_nodal == false)
              dof_info[no].vector_partitioner_face_variants[1] =
                dof_info[no].vector_partitioner;
            else
              {
                bool has_noncontiguous_cell = false;
                for (unsigned int f = 0; f < n_inner_face_batches(); ++f)
                  for (unsigned int v = 0;
                       v < VectorizedArray<Number>::n_array_elements &&
                       face_info.faces[f].cells_interior[v] !=
                         numbers::invalid_unsigned_int;
                       ++v)
                    {
                      AssertIndexRange(
                        face_info.faces[f].cells_interior[v],
                        n_macro_cells_before *
                          VectorizedArray<Number>::n_array_elements);
                      if (dof_info[no].index_storage_variants
                              [internal::MatrixFreeFunctions::DoFInfo::
                                 dof_access_face_exterior][f] >=
                            internal::MatrixFreeFunctions::DoFInfo::
                              IndexStorageVariants::contiguous &&
                          dof_info[no].dof_indices_contiguous
                              [internal::MatrixFreeFunctions::DoFInfo::
                                 dof_access_cell]
                              [face_info.faces[f].cells_exterior[v]] >=
                            part.local_size())
                        {
                          const unsigned int p =
                            face_info.faces[f].cells_exterior[v];
                          const unsigned int stride = 1;
                          unsigned int       i      = 0;
                          for (unsigned int e = 0;
                               e < dof_info[no].n_base_elements;
                               ++e)
                            for (unsigned int c = 0;
                                 c < dof_info[no].n_components[e];
                                 ++c)
                              {
                                const internal::MatrixFreeFunctions::ShapeInfo<
                                  VectorizedArray<Number>> &shape =
                                  shape_info(
                                    dof_info[no].global_base_element_offset + e,
                                    0,
                                    0,
                                    0);
                                for (unsigned int j = 0;
                                     j < shape.dofs_per_component_on_face;
                                     ++j)
                                  ghost_indices.push_back(part.local_to_global(
                                    dof_info[no].dof_indices_contiguous
                                      [internal::MatrixFreeFunctions::DoFInfo::
                                         dof_access_cell][p] +
                                    i +
                                    shape.face_to_cell_index_nodal(
                                      face_info.faces[f].exterior_face_no, j) *
                                      stride));
                                i += shape.dofs_per_component_on_cell * stride;
                              }
                          AssertDimension(
                            i, dof_info[no].dofs_per_cell[0] * stride);
                        }
                      else if (dof_info[no].index_storage_variants
                                 [internal::MatrixFreeFunctions::DoFInfo::
                                    dof_access_face_exterior][f] <
                               internal::MatrixFreeFunctions::DoFInfo::
                                 IndexStorageVariants::contiguous)
                        has_noncontiguous_cell = true;
                    }
                has_noncontiguous_cell = Utilities::MPI::min(
                  (int)has_noncontiguous_cell, task_info.communicator);

                std::sort(ghost_indices.begin(), ghost_indices.end());
                ghost_indices.erase(
                  std::unique(ghost_indices.begin(), ghost_indices.end()),
                  ghost_indices.end());
                IndexSet compressed_set(part.size());
                compressed_set.add_indices(ghost_indices.begin(),
                                           ghost_indices.end());
                compressed_set.subtract_set(
                  dof_info[no].vector_partitioner->locally_owned_range());
                const bool all_ghosts_equal = Utilities::MPI::min(
                  (int)(compressed_set.n_elements() ==
                        dof_info[no]
                          .vector_partitioner->ghost_indices()
                          .n_elements()),
                  dof_info[no].vector_partitioner->get_mpi_communicator());
                if (all_ghosts_equal || has_noncontiguous_cell)
                  dof_info[no].vector_partitioner_face_variants[1] =
                    dof_info[no].vector_partitioner;
                else
                  {
                    dof_info[no].vector_partitioner_face_variants[1].reset(
                      new Utilities::MPI::Partitioner(
                        part.locally_owned_range(),
                        part.get_mpi_communicator()));
                    const_cast<Utilities::MPI::Partitioner *>(
                      dof_info[no].vector_partitioner_face_variants[1].get())
                      ->set_ghost_indices(compressed_set, part.ghost_indices());
                  }
              }
          }

          // partitioner 2: values and gradients on faces
          {
            bool all_hermite = true;
            for (unsigned int c = 0; c < dof_info[no].n_base_elements; ++c)
              if (shape_info(
                    dof_info[no].global_base_element_offset + c, 0, 0, 0)
                    .element_type !=
                  internal::MatrixFreeFunctions::tensor_symmetric_hermite)
                all_hermite = false;
            if (all_hermite == false ||
                dof_info[no].vector_partitioner_face_variants[1].get() ==
                  dof_info[no].vector_partitioner.get())
              dof_info[no].vector_partitioner_face_variants[2] =
                dof_info[no].vector_partitioner;
            else
              {
                for (unsigned int f = 0; f < n_inner_face_batches(); ++f)
                  for (unsigned int v = 0;
                       v < VectorizedArray<Number>::n_array_elements &&
                       face_info.faces[f].cells_interior[v] !=
                         numbers::invalid_unsigned_int;
                       ++v)
                    {
                      AssertIndexRange(
                        face_info.faces[f].cells_interior[v],
                        n_macro_cells_before *
                          VectorizedArray<Number>::n_array_elements);
                      if (dof_info[no].index_storage_variants
                              [internal::MatrixFreeFunctions::DoFInfo::
                                 dof_access_face_exterior][f] >=
                            internal::MatrixFreeFunctions::DoFInfo::
                              IndexStorageVariants::contiguous &&
                          dof_info[no].dof_indices_contiguous
                              [internal::MatrixFreeFunctions::DoFInfo::
                                 dof_access_cell]
                              [face_info.faces[f].cells_exterior[v]] >=
                            part.local_size())
                        {
                          const unsigned int p =
                            face_info.faces[f].cells_exterior[v];
                          const unsigned int stride = 1;
                          unsigned int       i      = 0;
                          for (unsigned int e = 0;
                               e < dof_info[no].n_base_elements;
                               ++e)
                            for (unsigned int c = 0;
                                 c < dof_info[no].n_components[e];
                                 ++c)
                              {
                                const internal::MatrixFreeFunctions::ShapeInfo<
                                  VectorizedArray<Number>> &shape =
                                  shape_info(
                                    dof_info[no].global_base_element_offset + e,
                                    0,
                                    0,
                                    0);
                                for (unsigned int j = 0;
                                     j < 2 * shape.dofs_per_component_on_face;
                                     ++j)
                                  ghost_indices.push_back(part.local_to_global(
                                    dof_info[no].dof_indices_contiguous
                                      [internal::MatrixFreeFunctions::DoFInfo::
                                         dof_access_cell][p] +
                                    i +
                                    shape.face_to_cell_index_hermite(
                                      face_info.faces[f].exterior_face_no, j) *
                                      stride));
                                i += shape.dofs_per_component_on_cell * stride;
                              }
                          AssertDimension(
                            i, dof_info[no].dofs_per_cell[0] * stride);
                        }
                    }
                std::sort(ghost_indices.begin(), ghost_indices.end());
                ghost_indices.erase(
                  std::unique(ghost_indices.begin(), ghost_indices.end()),
                  ghost_indices.end());
                IndexSet compressed_set(part.size());
                compressed_set.add_indices(ghost_indices.begin(),
                                           ghost_indices.end());
                compressed_set.subtract_set(
                  dof_info[no].vector_partitioner->locally_owned_range());
                const bool all_ghosts_equal = Utilities::MPI::min(
                  (int)(compressed_set.n_elements() ==
                        dof_info[no]
                          .vector_partitioner->ghost_indices()
                          .n_elements()),
                  dof_info[no].vector_partitioner->get_mpi_communicator());
                if (all_ghosts_equal)
                  dof_info[no].vector_partitioner_face_variants[2] =
                    dof_info[no].vector_partitioner;
                else
                  {
                    dof_info[no].vector_partitioner_face_variants[2].reset(
                      new Utilities::MPI::Partitioner(
                        part.locally_owned_range(),
                        part.get_mpi_communicator()));
                    const_cast<Utilities::MPI::Partitioner *>(
                      dof_info[no].vector_partitioner_face_variants[2].get())
                      ->set_ghost_indices(compressed_set, part.ghost_indices());
                  }
              }
          }
        }
    }

  for (unsigned int no = 0; no < n_fe; ++no)
    dof_info[no].compute_vector_zero_access_pattern(task_info, face_info.faces);

  indices_are_initialized = true;
}



template <int dim, typename Number>
void
MatrixFree<dim, Number>::clear()
{
  dof_info.clear();
  mapping_info.clear();
  cell_level_index.clear();
  task_info.clear();
  dof_handlers.dof_handler.clear();
  dof_handlers.hp_dof_handler.clear();
  face_info.clear();
  indices_are_initialized = false;
  mapping_is_initialized  = false;
}



#ifdef DEAL_II_WITH_THREADS

namespace internal
{
  namespace
  {
    void
    fill_index_subrange(
      const unsigned int begin,
      const unsigned int end,
      const std::vector<std::pair<unsigned int, unsigned int>>
        &                                          cell_level_index,
      tbb::concurrent_unordered_map<std::pair<unsigned int, unsigned int>,
                                    unsigned int> &map)
    {
      if (cell_level_index.empty())
        return;
      unsigned int cell = begin;
      if (cell == 0)
        map.insert(std::make_pair(cell_level_index[cell++], 0U));
      for (; cell < end; ++cell)
        if (cell_level_index[cell] != cell_level_index[cell - 1])
          map.insert(std::make_pair(cell_level_index[cell], cell));
    }

    template <int dim>
    void
    fill_connectivity_subrange(
      const unsigned int                begin,
      const unsigned int                end,
      const dealii::Triangulation<dim> &tria,
      const std::vector<std::pair<unsigned int, unsigned int>>
        &                                                cell_level_index,
      const tbb::concurrent_unordered_map<std::pair<unsigned int, unsigned int>,
                                          unsigned int> &map,
      DynamicSparsityPattern &                           connectivity_direct)
    {
      std::vector<types::global_dof_index> new_indices;
      for (unsigned int cell = begin; cell < end; ++cell)
        {
          new_indices.clear();
          typename dealii::Triangulation<dim>::cell_iterator dcell(
            &tria, cell_level_index[cell].first, cell_level_index[cell].second);
          for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
            {
              // Only inner faces couple different cells
              if (dcell->at_boundary(f) == false &&
                  dcell->neighbor_or_periodic_neighbor(f)
                      ->level_subdomain_id() == dcell->level_subdomain_id())
                {
                  std::pair<unsigned int, unsigned int> level_index(
                    dcell->neighbor_or_periodic_neighbor(f)->level(),
                    dcell->neighbor_or_periodic_neighbor(f)->index());
                  auto it = map.find(level_index);
                  if (it != map.end())
                    {
                      const unsigned int neighbor_cell = it->second;
                      if (neighbor_cell != cell)
                        new_indices.push_back(neighbor_cell);
                    }
                }
            }
          std::sort(new_indices.begin(), new_indices.end());
          connectivity_direct.add_entries(
            cell,
            new_indices.begin(),
            std::unique(new_indices.begin(), new_indices.end()));
        }
    }

    void
    fill_connectivity_indirect_subrange(
      const unsigned int            begin,
      const unsigned int            end,
      const DynamicSparsityPattern &connectivity_direct,
      DynamicSparsityPattern &      connectivity)
    {
      std::vector<types::global_dof_index> new_indices;
      for (unsigned int block = begin; block < end; ++block)
        {
          new_indices.clear();
          for (DynamicSparsityPattern::iterator it =
                 connectivity_direct.begin(block);
               it != connectivity_direct.end(block);
               ++it)
            {
              new_indices.push_back(it->column());
              for (DynamicSparsityPattern::iterator it_neigh =
                     connectivity_direct.begin(it->column());
                   it_neigh != connectivity_direct.end(it->column());
                   ++it_neigh)
                if (it_neigh->column() != block)
                  new_indices.push_back(it_neigh->column());
            }
          std::sort(new_indices.begin(), new_indices.end());
          connectivity.add_entries(
            block,
            new_indices.begin(),
            std::unique(new_indices.begin(), new_indices.end()));
        }
    }
  } // namespace
} // namespace internal

#endif



template <int dim, typename Number>
void
MatrixFree<dim, Number>::make_connectivity_graph_faces(
  DynamicSparsityPattern &connectivity)
{
  (void)connectivity;
#ifdef DEAL_II_WITH_THREADS
  // step 1: build map between the index in the matrix-free context and the
  // one in the triangulation
  tbb::concurrent_unordered_map<std::pair<unsigned int, unsigned int>,
                                unsigned int>
    map;
  parallel::apply_to_subranges(0,
                               cell_level_index.size(),
                               std::bind(&internal::fill_index_subrange,
                                         std::placeholders::_1,
                                         std::placeholders::_2,
                                         std::cref(cell_level_index),
                                         std::ref(map)),
                               50);

  // step 2: Make a list for all blocks with other blocks that write to the
  // cell (due to the faces that are associated to it)
  DynamicSparsityPattern    connectivity_direct(connectivity.n_rows(),
                                             connectivity.n_cols());
  const Triangulation<dim> &tria =
    dof_handlers.active_dof_handler == DoFHandlers::usual ?
      dof_handlers.dof_handler[0]->get_triangulation() :
      dof_handlers.hp_dof_handler[0]->get_triangulation();
  parallel::apply_to_subranges(
    0,
    task_info.n_active_cells,
    std::bind(&internal::fill_connectivity_subrange<dim>,
              std::placeholders::_1,
              std::placeholders::_2,
              std::cref(tria),
              std::cref(cell_level_index),
              std::cref(map),
              std::ref(connectivity_direct)),
    20);
  connectivity_direct.symmetrize();

  // step 3: Include also interaction between neighbors one layer away because
  // faces might be assigned to cells differently
  parallel::apply_to_subranges(
    0,
    task_info.n_active_cells,
    std::bind(&internal::fill_connectivity_indirect_subrange,
              std::placeholders::_1,
              std::placeholders::_2,
              std::cref(connectivity_direct),
              std::ref(connectivity)),
    20);
#endif
}



template <int dim, typename Number>
std::size_t
MatrixFree<dim, Number>::memory_consumption() const
{
  std::size_t memory = MemoryConsumption::memory_consumption(dof_info);
  memory += MemoryConsumption::memory_consumption(cell_level_index);
  memory += MemoryConsumption::memory_consumption(face_info);
  memory += MemoryConsumption::memory_consumption(shape_info);
  memory += MemoryConsumption::memory_consumption(constraint_pool_data);
  memory += MemoryConsumption::memory_consumption(constraint_pool_row_index);
  memory += MemoryConsumption::memory_consumption(task_info);
  memory += sizeof(*this);
  memory += mapping_info.memory_consumption();
  return memory;
}



template <int dim, typename Number>
template <typename StreamType>
void
MatrixFree<dim, Number>::print_memory_consumption(StreamType &out) const
{
  out << "  Memory matrix-free data total: --> ";
  task_info.print_memory_statistics(out, memory_consumption());
  out << "   Memory cell index:                ";
  task_info.print_memory_statistics(
    out, MemoryConsumption::memory_consumption(cell_level_index));
  if (Utilities::MPI::sum(face_info.faces.size(), task_info.communicator) > 0)
    {
      out << "   Memory face indicators:           ";
      task_info.print_memory_statistics(
        out, MemoryConsumption::memory_consumption(face_info.faces));
    }
  for (unsigned int j = 0; j < dof_info.size(); ++j)
    {
      out << "   Memory DoFInfo component " << j << std::endl;
      dof_info[j].print_memory_consumption(out, task_info);
    }

  out << "   Memory mapping info" << std::endl;
  mapping_info.print_memory_consumption(out, task_info);

  out << "   Memory unit cell shape data:      ";
  task_info.print_memory_statistics(
    out, MemoryConsumption::memory_consumption(shape_info));
  if (task_info.scheme != internal::MatrixFreeFunctions::TaskInfo::none)
    {
      out << "   Memory task partitioning info:    ";
      task_info.print_memory_statistics(
        out, MemoryConsumption::memory_consumption(task_info));
    }
}



template <int dim, typename Number>
void
MatrixFree<dim, Number>::print(std::ostream &out) const
{
  // print indices local to global
  for (unsigned int no = 0; no < dof_info.size(); ++no)
    {
      out << "\n-- Index data for component " << no << " --" << std::endl;
      dof_info[no].print(constraint_pool_data, constraint_pool_row_index, out);
      out << std::endl;
    }
}



DEAL_II_NAMESPACE_CLOSE

#endif
