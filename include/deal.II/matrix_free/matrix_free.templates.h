// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2017 by the deal.II authors
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


#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/polynomials_piecewise.h>
#include <deal.II/base/mpi.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/shape_info.templates.h>
#include <deal.II/matrix_free/mapping_info.templates.h>
#include <deal.II/matrix_free/dof_info.templates.h>
#include <deal.II/matrix_free/face_info.h>


DEAL_II_NAMESPACE_OPEN


// --------------------- MatrixFree -----------------------------------

template <int dim, typename Number>
MatrixFree<dim, Number>::MatrixFree()
  :
  Subscriptor(),
  indices_are_initialized (false),
  mapping_is_initialized  (false)
{}



template <int dim, typename Number>
MatrixFree<dim, Number>::MatrixFree(const MatrixFree<dim,Number> &other)
  :
  Subscriptor()
{
  copy_from(other);
}



template <int dim, typename Number>
void MatrixFree<dim,Number>::
copy_from (const MatrixFree<dim,Number> &v)
{
  clear ();
  dof_handlers = v.dof_handlers;
  dof_info = v.dof_info;
  constraint_pool_data = v.constraint_pool_data;
  constraint_pool_row_index = v.constraint_pool_row_index;
  mapping_info = v.mapping_info;
  shape_info = v.shape_info;
  cell_level_index = v.cell_level_index;
  task_info = v.task_info;
  indices_are_initialized = v.indices_are_initialized;
  mapping_is_initialized  = v.mapping_is_initialized;
}



template <int dim, typename Number>
void MatrixFree<dim,Number>::
internal_reinit(const Mapping<dim>                          &mapping,
                const std::vector<const DoFHandler<dim> *>  &dof_handler,
                const std::vector<const ConstraintMatrix *> &constraint,
                const std::vector<IndexSet>                 &locally_owned_set,
                const std::vector<hp::QCollection<1> >      &quad,
                const typename MatrixFree<dim,Number>::AdditionalData additional_data)
{

  // Reads out the FE information and stores the shape function values,
  // gradients and Hessians for quadrature points.
  {
    const unsigned int n_fe   = dof_handler.size();
    const unsigned int n_quad = quad.size();
    shape_info.reinit (TableIndices<4>(n_fe, n_quad, 1, 1));
    for (unsigned int no=0; no<n_fe; no++)
      for (unsigned int nq =0; nq<n_quad; nq++)
        {
          AssertDimension (quad[nq].size(), 1);
          shape_info(no,nq,0,0).reinit(quad[nq][0], dof_handler[no]->get_fe());
        }
  }

  if (additional_data.initialize_indices == true)
    {
      clear();
      Assert (dof_handler.size() > 0, ExcMessage("No DoFHandler is given."));
      AssertDimension (dof_handler.size(), constraint.size());
      AssertDimension (dof_handler.size(), locally_owned_set.size());

      // set variables that are independent of FE
      if (Utilities::MPI::job_supports_mpi() == true)
        {
          const parallel::Triangulation<dim> *dist_tria =
            dynamic_cast<const parallel::Triangulation<dim>*>
            (&(dof_handler[0]->get_triangulation()));
          task_info.communicator = dist_tria != nullptr ?
                                   dist_tria->get_communicator() :
                                   MPI_COMM_SELF;
          task_info.my_pid  =
            Utilities::MPI::this_mpi_process(task_info.communicator);
          task_info.n_procs =
            Utilities::MPI::n_mpi_processes(task_info.communicator);
        }
      else
        {
          task_info.communicator = MPI_COMM_SELF;
          task_info.my_pid = 0;
          task_info.n_procs = 1;
        }

      initialize_dof_handlers (dof_handler, additional_data.level_mg_handler);
      for (unsigned int no=0; no<dof_handler.size(); ++no)
        dof_info[no].store_plain_indices = additional_data.store_plain_indices;

      // initialize the basic multithreading information that needs to be
      // passed to the DoFInfo structure
#ifdef DEAL_II_WITH_THREADS
      if (additional_data.tasks_parallel_scheme != AdditionalData::none &&
          MultithreadInfo::n_threads() > 1)
        {
          task_info.scheme = internal::MatrixFreeFunctions::TaskInfo::TasksParallelScheme(static_cast<int>(additional_data.tasks_parallel_scheme));
          task_info.block_size = additional_data.tasks_block_size;
        }
      else
#endif
        task_info.scheme = internal::MatrixFreeFunctions::TaskInfo::none;

      // set dof_indices together with constraint_indicator and
      // constraint_pool_data. It also reorders the way cells are gone through
      // (to separate cells with overlap to other processors from others
      // without).
      initialize_indices (constraint, locally_owned_set);
    }

  // initialize bare structures
  else if (dof_info.size() != dof_handler.size())
    {
      initialize_dof_handlers(dof_handler, additional_data.level_mg_handler);
      std::vector<unsigned int> dummy;
      std::vector<unsigned char> dummy2;
      task_info.collect_boundary_cells (cell_level_index.size(), cell_level_index.size(),
                                        VectorizedArray<Number>::n_array_elements, dummy);
      task_info.create_blocks_serial(dummy, dummy, 1, dummy, false, dummy, dummy2);
      for (unsigned int i=0; i<dof_info.size(); ++i)
        {
          dof_info[i].dimension    = dim;
          dof_info[i].n_components = dof_handler[i]->get_fe().element_multiplicity(0);
          dof_info[i].dofs_per_cell.push_back(dof_handler[i]->get_fe().dofs_per_cell);
          dof_info[i].row_starts.resize(task_info.cell_partition_data.back()+1);
          dof_info[i].row_starts.back()[2] =
            cell_level_index.size() % VectorizedArray<Number>::n_array_elements;

          // if indices are not initialized, the cell_level_index might not be
          // divisible by the vectorization length. But it must be for
          // mapping_info...
          while (cell_level_index.size() % VectorizedArray<Number>::n_array_elements
                 != 0)
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
      mapping_info.initialize (dof_handler[0]->get_triangulation(), cell_level_index,
                               internal::MatrixFreeFunctions:: FaceInfo
                               <VectorizedArray<Number>::n_array_elements>(),
                               dummy, mapping,
                               quad, additional_data.mapping_update_flags,
                               update_default, update_default, update_default);

      mapping_is_initialized = true;
    }
}



template <int dim, typename Number>
void MatrixFree<dim,Number>::
internal_reinit(const Mapping<dim>                            &mapping,
                const std::vector<const hp::DoFHandler<dim>*> &dof_handler,
                const std::vector<const ConstraintMatrix *>    &constraint,
                const std::vector<IndexSet>                   &locally_owned_set,
                const std::vector<hp::QCollection<1> >        &quad,
                const typename MatrixFree<dim,Number>::AdditionalData additional_data)
{
  // Reads out the FE information and stores the shape function values,
  // gradients and Hessians for quadrature points.
  {
    const unsigned int n_components = dof_handler.size();
    const unsigned int n_quad       = quad.size();
    unsigned int n_fe_in_collection = 0;
    for (unsigned int i=0; i<n_components; ++i)
      n_fe_in_collection = std::max (n_fe_in_collection,
                                     dof_handler[i]->get_fe_collection().size());
    unsigned int n_quad_in_collection = 0;
    for (unsigned int q=0; q<n_quad; ++q)
      n_quad_in_collection = std::max (n_quad_in_collection, quad[q].size());
    shape_info.reinit (TableIndices<4>(n_components, n_quad,
                                       n_fe_in_collection,
                                       n_quad_in_collection));
    for (unsigned int no=0; no<n_components; no++)
      for (unsigned int fe_no=0; fe_no<dof_handler[no]->get_fe_collection().size(); ++fe_no)
        for (unsigned int nq =0; nq<n_quad; nq++)
          for (unsigned int q_no=0; q_no<quad[nq].size(); ++q_no)
            shape_info(no,nq,fe_no,q_no).reinit (quad[nq][q_no],
                                                 dof_handler[no]->get_fe(fe_no));
  }

  if (additional_data.initialize_indices == true)
    {
      clear();
      Assert (dof_handler.size() > 0, ExcMessage("No DoFHandler is given."));
      AssertDimension (dof_handler.size(), constraint.size());
      AssertDimension (dof_handler.size(), locally_owned_set.size());

      // set variables that are independent of FE
      if (Utilities::MPI::job_supports_mpi() == true)
        {
          const parallel::Triangulation<dim> *dist_tria =
            dynamic_cast<const parallel::Triangulation<dim>*>
            (&(dof_handler[0]->get_triangulation()));
          task_info.communicator = dist_tria != nullptr ?
                                   dist_tria->get_communicator() :
                                   MPI_COMM_SELF;
          task_info.my_pid  =
            Utilities::MPI::this_mpi_process(task_info.communicator);
          task_info.n_procs =
            Utilities::MPI::n_mpi_processes(task_info.communicator);
        }
      else
        {
          task_info.communicator = MPI_COMM_SELF;
          task_info.my_pid = 0;
          task_info.n_procs = 1;
        }

      initialize_dof_handlers (dof_handler, additional_data.level_mg_handler);
      for (unsigned int no=0; no<dof_handler.size(); ++no)
        dof_info[no].store_plain_indices = additional_data.store_plain_indices;

      // initialize the basic multithreading information that needs to be
      // passed to the DoFInfo structure
#ifdef DEAL_II_WITH_THREADS
      if (additional_data.tasks_parallel_scheme != AdditionalData::none &&
          MultithreadInfo::n_threads() > 1)
        {
          task_info.scheme = internal::MatrixFreeFunctions::TaskInfo::TasksParallelScheme(static_cast<int>(additional_data.tasks_parallel_scheme));
          task_info.block_size = additional_data.tasks_block_size;
        }
      else
#endif
        task_info.scheme = internal::MatrixFreeFunctions::TaskInfo::none;

      // set dof_indices together with constraint_indicator and
      // constraint_pool_data. It also reorders the way cells are gone through
      // (to separate cells with overlap to other processors from others
      // without).
      initialize_indices (constraint, locally_owned_set);
    }

  // initialize bare structures
  else if (dof_info.size() != dof_handler.size())
    {
      initialize_dof_handlers(dof_handler, additional_data.level_mg_handler);
      std::vector<unsigned int> dummy;
      std::vector<unsigned char> dummy2;
      task_info.collect_boundary_cells (cell_level_index.size(), cell_level_index.size(),
                                        VectorizedArray<Number>::n_array_elements, dummy);
      task_info.create_blocks_serial(dummy, dummy, 1, dummy, false, dummy, dummy2);
      for (unsigned int i=0; i<dof_info.size(); ++i)
        {
          Assert(dof_handler[i]->get_fe_collection().size() == 1, ExcNotImplemented());
          dof_info[i].dimension    = dim;
          dof_info[i].n_components = dof_handler[i]->get_fe(0).element_multiplicity(0);
          dof_info[i].dofs_per_cell.push_back(dof_handler[i]->get_fe(0).dofs_per_cell);
          dof_info[i].row_starts.resize(task_info.cell_partition_data.back()+1);
          dof_info[i].row_starts.back()[2] =
            cell_level_index.size() % VectorizedArray<Number>::n_array_elements;

          // if indices are not initialized, the cell_level_index might not be
          // divisible by the vectorization length. But it must be for
          // mapping_info...
          while (cell_level_index.size() % VectorizedArray<Number>::n_array_elements
                 != 0)
            cell_level_index.push_back(cell_level_index.back());
        }
    }

  // Evaluates transformations from unit to real cell, Jacobian determinants,
  // quadrature points in real space, based on the ordering of the cells
  // determined in @p extract_local_to_global_indices.
  if (additional_data.initialize_mapping == true)
    {
      mapping_info.initialize (dof_handler[0]->get_triangulation(), cell_level_index,
                               internal::MatrixFreeFunctions::FaceInfo
                               <VectorizedArray<Number>::n_array_elements>(),
                               dof_info[0].cell_active_fe_index, mapping,
                               quad, additional_data.mapping_update_flags,
                               update_default, update_default, update_default);

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
  if (fe.degree==0 || fe.n_base_elements()!=1)
    return false;

  const FiniteElement<dim, spacedim> *fe_ptr = &(fe.base_element(0));
  if (fe_ptr->n_components() != 1)
    return false;

  // then check of the base element is supported
  if (dynamic_cast<const FE_Poly<TensorProductPolynomials<dim>,dim,spacedim>*>(fe_ptr)!=nullptr)
    return true;
  if (dynamic_cast<const FE_Poly<TensorProductPolynomials<dim,
      Polynomials::PiecewisePolynomial<double> >,dim,spacedim>*>(fe_ptr)!=nullptr)
    return true;
  if (dynamic_cast<const FE_DGP<dim, spacedim>*>(fe_ptr)!=nullptr)
    return true;
  if (dynamic_cast<const FE_Q_DG0<dim, spacedim>*>(fe_ptr)!=nullptr)
    return true;

  // if the base element is not in the above list it is not supported
  return false;
}



namespace internal
{

  // steps through all children and adds the
  // active cells recursively
  template <typename InIterator>
  void resolve_cell (const InIterator   &cell,
                     std::vector<std::pair<unsigned int,unsigned int> > &cell_its,
                     const unsigned int  subdomain_id)
  {
    if (cell->has_children())
      for (unsigned int child=0; child<cell->n_children(); ++child)
        resolve_cell (cell->child(child), cell_its,
                      subdomain_id);
    else if (subdomain_id == numbers::invalid_subdomain_id
             || cell->subdomain_id() == subdomain_id)
      {
        Assert (cell->active(), ExcInternalError());
        cell_its.emplace_back (cell->level(), cell->index());
      }
  }
}



template <int dim, typename Number>
void MatrixFree<dim,Number>::
initialize_dof_handlers (const std::vector<const DoFHandler<dim>*> &dof_handler,
                         const unsigned int level)
{
  dof_handlers.active_dof_handler = DoFHandlers::usual;
  dof_handlers.level = level;
  dof_handlers.n_dof_handlers = dof_handler.size();
  dof_handlers.dof_handler.resize (dof_handlers.n_dof_handlers);
  for (unsigned int no=0; no<dof_handlers.n_dof_handlers; ++no)
    dof_handlers.dof_handler[no] = dof_handler[no];

  dof_info.resize (dof_handlers.n_dof_handlers);

  // Go through cells on zeroth level and then successively step down into
  // children. This gives a z-ordering of the cells, which is beneficial when
  // setting up neighboring relations between cells for thread parallelization
  const unsigned int n_mpi_procs = task_info.n_procs;
  const unsigned int my_pid = task_info.my_pid;

  const Triangulation<dim> &tria = dof_handlers.dof_handler[0]->get_triangulation();
  if (level == numbers::invalid_unsigned_int)
    {
      if (n_mpi_procs == 1)
        cell_level_index.reserve (tria.n_active_cells());
      typename Triangulation<dim>::cell_iterator cell = tria.begin(0),
                                                 end_cell = tria.end(0);
      // For serial Triangulations always take all cells
      const unsigned int subdomain_id
        = (dynamic_cast<const parallel::Triangulation<dim> *>
           (&dof_handler[0]->get_triangulation())!=nullptr)
          ? my_pid : numbers::invalid_subdomain_id;
      for ( ; cell != end_cell; ++cell)
        internal::resolve_cell (cell, cell_level_index, subdomain_id);

      Assert(n_mpi_procs>1 || cell_level_index.size()==tria.n_active_cells(),
             ExcInternalError());
    }
  else
    {
      AssertIndexRange (level, tria.n_global_levels());
      if (level < tria.n_levels())
        {
          cell_level_index.reserve (tria.n_cells(level));
          typename Triangulation<dim>::cell_iterator cell = tria.begin(level),
                                                     end_cell = tria.end(level);
          for ( ; cell != end_cell; ++cell)
            if (cell->level_subdomain_id() == my_pid)
              cell_level_index.emplace_back (cell->level(), cell->index());
        }
    }

  // All these are cells local to this processor. Therefore, set
  // cell_level_index_end_local to the size of cell_level_index.
  cell_level_index_end_local = cell_level_index.size();
}



template <int dim, typename Number>
void MatrixFree<dim,Number>::
initialize_dof_handlers (const std::vector<const hp::DoFHandler<dim>*> &dof_handler,
                         const unsigned int)
{
  dof_handlers.active_dof_handler = DoFHandlers::hp;
  dof_handlers.n_dof_handlers = dof_handler.size();
  dof_handlers.hp_dof_handler.resize (dof_handlers.n_dof_handlers);
  for (unsigned int no=0; no<dof_handlers.n_dof_handlers; ++no)
    dof_handlers.hp_dof_handler[no] = dof_handler[no];

  dof_info.resize (dof_handlers.n_dof_handlers);

  // go through cells on zeroth level and then successively step down into
  // children. This gives a z-ordering of the cells, which is beneficial when
  // setting up neighboring relations between cells for thread parallelization
  const unsigned int n_mpi_procs = task_info.n_procs;
  const unsigned int my_pid = task_info.my_pid;

  // if we have no level given, use the same as for the standard DoFHandler,
  // otherwise we must loop through the respective level
  const Triangulation<dim> &tria = dof_handler[0]->get_triangulation();

  if (n_mpi_procs == 1)
    {
      cell_level_index.reserve (tria.n_active_cells());
    }
  typename hp::DoFHandler<dim>::cell_iterator cell = dof_handler[0]->begin(0),
                                              end_cell = dof_handler[0]->end(0);
  // For serial Triangulations always take all cells
  const unsigned int subdomain_id
    = (dynamic_cast<const parallel::Triangulation<dim> *>
       (&dof_handler[0]->get_triangulation())!=nullptr)
      ? my_pid : numbers::invalid_subdomain_id;
  for ( ; cell != end_cell; ++cell)
    internal::resolve_cell (cell, cell_level_index,
                            subdomain_id);

  Assert(n_mpi_procs>1 || cell_level_index.size()==tria.n_active_cells(),
         ExcInternalError());

  // All these are cells local to this processor. Therefore, set
  // cell_level_index_end_local to the size of cell_level_index.
  cell_level_index_end_local = cell_level_index.size();
}



template <int dim, typename Number>
void MatrixFree<dim,Number>::initialize_indices
(const std::vector<const ConstraintMatrix *> &constraint,
 const std::vector<IndexSet>                 &locally_owned_set)
{
  const unsigned int n_fe = dof_handlers.n_dof_handlers;
  const unsigned int n_active_cells = cell_level_index.size();

  AssertDimension (n_active_cells, cell_level_index.size());
  AssertDimension (n_fe, locally_owned_set.size());
  AssertDimension (n_fe, constraint.size());

  std::vector<types::global_dof_index> local_dof_indices;

  internal::MatrixFreeFunctions::ConstraintValues<double> constraint_values;

  for (unsigned int no=0; no<n_fe; ++no)
    {
      std::vector<const FiniteElement<dim>*> fes;
      if (dof_handlers.active_dof_handler == DoFHandlers::hp)
        {
          const hp::DoFHandler<dim> *hpdof = dof_handlers.hp_dof_handler[no];
          const hp::FECollection<dim> &fe = hpdof->get_fe_collection();
          for (unsigned int f=0; f<fe.size(); ++f)
            fes.push_back (&fe[f]);

          dof_info[no].max_fe_index = fe.size();
          dof_info[no].fe_index_conversion.resize (fe.size());
          for (unsigned int ind=0; ind<hpdof->get_fe_collection().size(); ++ind)
            dof_info[no].fe_index_conversion[ind] =
              std::pair<unsigned int,unsigned int>(fe[ind].degree,
                                                   fe[ind].dofs_per_cell);
          if (fe.size() > 1)
            dof_info[no].cell_active_fe_index.resize(n_active_cells,
                                                     numbers::invalid_unsigned_int);
        }
      else
        {
          const DoFHandler<dim> *dofh =&*dof_handlers.dof_handler[no];
          fes.push_back (&dofh->get_fe());
          dof_info[no].max_fe_index = 1;
          dof_info[no].fe_index_conversion.resize (1);
          dof_info[no].fe_index_conversion[0] =
            std::pair<unsigned int,unsigned int>(fes.back()->degree,
                                                 fes.back()->dofs_per_cell);
        }

      for (unsigned int fe_index = 0; fe_index<fes.size(); ++fe_index)
        {
          const FiniteElement<dim> &fe = *fes[fe_index];
          Assert (fe.n_base_elements() == 1,
                  ExcMessage ("MatrixFree currently only works for DoFHandler with one base element"));
          const unsigned int n_fe_components = fe.element_multiplicity (0);

          // cache number of finite elements and dofs_per_cell
          dof_info[no].dofs_per_cell.push_back (fe.dofs_per_cell);
          dof_info[no].dofs_per_face.push_back (fe.dofs_per_face);
          dof_info[no].dimension    = dim;
          dof_info[no].n_components = n_fe_components;

          AssertDimension (shape_info(no,0,fe_index,0).lexicographic_numbering.size(),
                           dof_info[no].dofs_per_cell[fe_index]);
        }

      // set locally owned range for each component
      Assert (locally_owned_set[no].is_contiguous(), ExcNotImplemented());
      dof_info[no].vector_partitioner.reset
      (new Utilities::MPI::Partitioner(locally_owned_set[no], task_info.communicator));

      // initialize the arrays for indices
      dof_info[no].row_starts.resize (n_active_cells+1);
      dof_info[no].row_starts[0][0] = 0;
      dof_info[no].row_starts[0][1] = 0;
      dof_info[no].row_starts[0][2] = 0;
      dof_info[no].dof_indices.reserve
      ((n_active_cells*dof_info[no].dofs_per_cell[0]*3)/2);

      // cache the constrained indices for use in matrix-vector products
      {
        const types::global_dof_index
        start_index = dof_info[no].vector_partitioner->local_range().first,
        end_index   = dof_info[no].vector_partitioner->local_range().second;
        for (types::global_dof_index i=start_index; i<end_index; ++i)
          if (constraint[no]->is_constrained(i)==true)
            dof_info[no].constrained_dofs.
            push_back(static_cast<unsigned int>(i-start_index));
      }
    }

  // extract all the global indices associated with the computation, and form
  // the ghost indices
  std::vector<unsigned int> boundary_cells;
  for (unsigned int counter = 0 ; counter < n_active_cells ; ++counter)
    {
      bool cell_at_boundary = false;
      for (unsigned int no=0; no<n_fe; ++no)
        {
          // OK, read indices from standard DoFHandler in the usual way
          if (dof_handlers.active_dof_handler == DoFHandlers::usual &&
              dof_handlers.level == numbers::invalid_unsigned_int)
            {
              const DoFHandler<dim> *dofh = &*dof_handlers.dof_handler[no];
              typename DoFHandler<dim>::active_cell_iterator
              cell_it (&dofh->get_triangulation(),
                       cell_level_index[counter].first,
                       cell_level_index[counter].second,
                       dofh);
              local_dof_indices.resize (dof_info[no].dofs_per_cell[0]);
              cell_it->get_dof_indices(local_dof_indices);
              dof_info[no].read_dof_indices (local_dof_indices,
                                             shape_info(no,0,0,0).lexicographic_numbering,
                                             *constraint[no], counter,
                                             constraint_values,
                                             cell_at_boundary);
            }
          // ok, now we are requested to use a level in a MG DoFHandler
          else if (dof_handlers.active_dof_handler == DoFHandlers::usual &&
                   dof_handlers.level != numbers::invalid_unsigned_int)
            {
              const DoFHandler<dim> *dofh = dof_handlers.dof_handler[no];
              AssertIndexRange (dof_handlers.level, dofh->get_triangulation().n_levels());
              typename DoFHandler<dim>::cell_iterator
              cell_it (&dofh->get_triangulation(),
                       cell_level_index[counter].first,
                       cell_level_index[counter].second,
                       dofh);
              local_dof_indices.resize (dof_info[no].dofs_per_cell[0]);
              cell_it->get_mg_dof_indices(local_dof_indices);
              dof_info[no].read_dof_indices (local_dof_indices,
                                             shape_info(no,0,0,0).lexicographic_numbering,
                                             *constraint[no], counter,
                                             constraint_values,
                                             cell_at_boundary);
            }
          else if (dof_handlers.active_dof_handler == DoFHandlers::hp)
            {
              const hp::DoFHandler<dim> *dofh =
                dof_handlers.hp_dof_handler[no];
              typename hp::DoFHandler<dim>::active_cell_iterator
              cell_it (&dofh->get_triangulation(),
                       cell_level_index[counter].first,
                       cell_level_index[counter].second,
                       dofh);
              if (dofh->get_fe_collection().size() > 1)
                dof_info[no].cell_active_fe_index[counter] =
                  cell_it->active_fe_index();
              local_dof_indices.resize (cell_it->get_fe().dofs_per_cell);
              cell_it->get_dof_indices(local_dof_indices);
              dof_info[no].read_dof_indices (local_dof_indices,
                                             shape_info(no,0,cell_it->active_fe_index(),0).lexicographic_numbering,
                                             *constraint[no], counter,
                                             constraint_values,
                                             cell_at_boundary);
            }
          else
            {
              Assert (false, ExcNotImplemented());
            }
        }

      // if we found dofs on some FE component that belong to other
      // processors, the cell is added to the boundary cells.
      if (cell_at_boundary == true && counter < cell_level_index_end_local)
        boundary_cells.push_back(counter);
    }

  const unsigned int vectorization_length =
    VectorizedArray<Number>::n_array_elements;
  task_info.collect_boundary_cells (cell_level_index_end_local,
                                    n_active_cells, vectorization_length,
                                    boundary_cells);

  // finalize the creation of ghosts
  for (unsigned int no=0; no<n_fe; ++no)
    dof_info[no].assign_ghosts (boundary_cells);

  std::vector<unsigned int> renumbering;
  std::vector<unsigned char> irregular_cells;
  if (task_info.scheme == internal::MatrixFreeFunctions::TaskInfo::none)
    {
      const bool strict_categories = dof_handlers.active_dof_handler == DoFHandlers::hp;
      unsigned int dofs_per_cell = 0;
      for (unsigned int no=0; no<dof_info.size(); ++no)
        dofs_per_cell = std::max(dofs_per_cell, dof_info[no].dofs_per_cell[0]);
      task_info.create_blocks_serial(boundary_cells, std::vector<unsigned int>(),
                                     dofs_per_cell,
                                     dof_info[0].cell_active_fe_index,
                                     strict_categories,
                                     renumbering, irregular_cells);
    }
  else
    {
      // For strategy with blocking before partitioning: reorganize the indices
      // in order to overlap communication in MPI with computations: Place all
      // cells with ghost indices into one chunk. Also reorder cells so that we
      // can parallelize by threads
      task_info.initial_setup_blocks_tasks(boundary_cells, renumbering,
                                           irregular_cells);
      task_info.guess_block_size (dof_info[0].dofs_per_cell[0]);

      unsigned int n_macro_cells_before = *(task_info.cell_partition_data.end()-2);
      unsigned int n_ghost_slots = *(task_info.cell_partition_data.end()-1)-
                                   n_macro_cells_before;

      unsigned int start_nonboundary = numbers::invalid_unsigned_int;

      if (task_info.scheme == internal::MatrixFreeFunctions::TaskInfo::partition_color)
        {
          // set up partitions. if we just use coloring without partitions, do
          // nothing here, assume all cells to belong to the zero partition (that
          // we otherwise use for MPI boundary cells)
          if (task_info.scheme == internal::MatrixFreeFunctions::TaskInfo::color)
            {
              start_nonboundary = task_info.n_procs > 1 ?
                                  std::min(((task_info.cell_partition_data[2]-
                                             task_info.cell_partition_data[1]+task_info.block_size-1)/
                                            task_info.block_size)*task_info.block_size,
                                           task_info.cell_partition_data[3]) : 0;
            }
          else
            {
              if (task_info.n_procs > 1)
                {
                  task_info.cell_partition_data[1] = 0;
                  task_info.cell_partition_data[2] = task_info.cell_partition_data[3];
                }
              start_nonboundary = task_info.cell_partition_data.back();
            }

          if (dof_handlers.active_dof_handler == DoFHandlers::hp)
            {
              irregular_cells.resize (0);
              irregular_cells.resize (task_info.cell_partition_data.back()+
                                      2*dof_info[0].max_fe_index);
              std::vector<std::vector<unsigned int> > renumbering_fe_index;
              renumbering_fe_index.resize(dof_info[0].max_fe_index);
              unsigned int counter;
              n_macro_cells_before = 0;
              for (counter=0; counter<std::min(start_nonboundary*vectorization_length,
                                               task_info.n_active_cells); counter++)
                {
                  AssertIndexRange (counter, renumbering.size());
                  AssertIndexRange (renumbering[counter],
                                    dof_info[0].cell_active_fe_index.size());
                  renumbering_fe_index[dof_info[0].cell_active_fe_index[renumbering[counter]]].
                  push_back(renumbering[counter]);
                }
              counter = 0;
              for (unsigned int j=0; j<dof_info[0].max_fe_index; j++)
                {
                  for (unsigned int jj=0; jj<renumbering_fe_index[j].size(); jj++)
                    renumbering[counter++] = renumbering_fe_index[j][jj];
                  irregular_cells[renumbering_fe_index[j].size()/vectorization_length+
                                  n_macro_cells_before] =
                                    renumbering_fe_index[j].size()%vectorization_length;
                  n_macro_cells_before += (renumbering_fe_index[j].size()+vectorization_length-1)/
                                          vectorization_length;
                  renumbering_fe_index[j].resize(0);
                }

              for (counter=start_nonboundary*vectorization_length;
                   counter<task_info.n_active_cells; counter++)
                {
                  renumbering_fe_index[dof_info[0].cell_active_fe_index.empty() ? 0 :
                                       dof_info[0].cell_active_fe_index[renumbering[counter]]].
                  push_back(renumbering[counter]);
                }
              counter = start_nonboundary * vectorization_length;
              for (unsigned int j=0; j<dof_info[0].max_fe_index; j++)
                {
                  for (unsigned int jj=0; jj<renumbering_fe_index[j].size(); jj++)
                    renumbering[counter++] = renumbering_fe_index[j][jj];
                  irregular_cells[renumbering_fe_index[j].size()/vectorization_length+
                                  n_macro_cells_before] =
                                    renumbering_fe_index[j].size()%vectorization_length;
                  n_macro_cells_before += (renumbering_fe_index[j].size()+vectorization_length-1)/
                                          vectorization_length;
                }
              AssertIndexRange (n_macro_cells_before,
                                task_info.cell_partition_data.back() + 2*dof_info[0].max_fe_index+1);
              irregular_cells.resize (n_macro_cells_before+n_ghost_slots);
              *(task_info.cell_partition_data.end()-2) = n_macro_cells_before;
              *(task_info.cell_partition_data.end()-1) = n_macro_cells_before+n_ghost_slots;
            }
        }

      task_info.n_blocks = (n_macro_cells()+task_info.block_size-1)/
                           task_info.block_size;

      DynamicSparsityPattern connectivity;
      connectivity.reinit(task_info.n_active_cells, task_info.n_active_cells);
      if (task_info.n_active_cells > 0)
        dof_info[0].make_connectivity_graph(task_info, renumbering, connectivity);

      task_info.make_thread_graph(dof_info[0].cell_active_fe_index,
                                  connectivity, renumbering, irregular_cells,
                                  dof_handlers.active_dof_handler == DoFHandlers::hp);

      Assert(irregular_cells.size() >= task_info.cell_partition_data.back(),
             ExcInternalError());

      irregular_cells.resize(task_info.cell_partition_data.back()+n_ghost_slots);
      if (n_ghost_slots > 0)
        {
          for (unsigned int i=task_info.cell_partition_data.back();
               i<task_info.cell_partition_data.back()+n_ghost_slots-1; ++i)
            irregular_cells[i] = 0;
          irregular_cells.back() = task_info.n_ghost_cells%vectorization_length;
        }

      {
        unsigned int n_cells = 0;
        for (unsigned int i=0; i<task_info.cell_partition_data.back(); ++i)
          n_cells += irregular_cells[i] > 0 ? irregular_cells[i] : vectorization_length;
        AssertDimension(n_cells, task_info.n_active_cells);
        n_cells = 0;
        for (unsigned int i=task_info.cell_partition_data.back();
             i<n_ghost_slots+task_info.cell_partition_data.back(); ++i)
          n_cells += irregular_cells[i] > 0 ? irregular_cells[i] : vectorization_length;
        AssertDimension(n_cells, task_info.n_ghost_cells);
      }

      task_info.cell_partition_data
      .push_back(task_info.cell_partition_data.back()+n_ghost_slots);
    }

  // Finally perform the renumbering of the degree of freedom number data. We
  // also want to group several cells together to one "macro-cell" for
  // vectorization (where the arithmetic operations will then be done
  // simultaneously).
#ifdef DEBUG
  {
    std::vector<unsigned int> sorted_renumbering (renumbering);
    std::sort (sorted_renumbering.begin(), sorted_renumbering.end());
    for (unsigned int i=0; i<sorted_renumbering.size(); ++i)
      Assert (sorted_renumbering[i] == i, ExcInternalError());
  }
#endif
  {
    std::vector<std::pair<unsigned int,unsigned int> >
    cell_level_index_old;
    cell_level_index.swap (cell_level_index_old);
    cell_level_index.reserve(task_info.cell_partition_data.back()*vectorization_length);
    unsigned int position_cell=0;
    for (unsigned int i=0; i<task_info.cell_partition_data.back(); ++i)
      {
        unsigned int n_comp = (irregular_cells[i]>0)?
                              irregular_cells[i] : vectorization_length;
        for (unsigned int j=0; j<n_comp; ++j)
          cell_level_index.push_back
          (cell_level_index_old[renumbering[position_cell+j]]);

        // generate a cell and level index also when we have not filled up
        // vectorization_length cells. This is needed for MappingInfo when the
        // transformation data is initialized. We just set the value to the
        // last valid cell in that case.
        for (unsigned int j=n_comp; j<vectorization_length; ++j)
          cell_level_index.push_back
          (cell_level_index_old[renumbering[position_cell+n_comp-1]]);
        position_cell += n_comp;
      }
    AssertDimension (position_cell, task_info.n_active_cells + task_info.n_ghost_cells);
    AssertDimension (cell_level_index.size(),task_info.cell_partition_data.back()*
                     vectorization_length);
  }

  // set constraint pool from the std::map and reorder the indices
  typename std::map<std::vector<double>, types::global_dof_index,
           internal::MatrixFreeFunctions::FPArrayComparator<double> >::iterator
           it = constraint_values.constraints.begin(),
           end = constraint_values.constraints.end();
  std::vector<const std::vector<double>*>
  constraints (constraint_values.constraints.size());
  unsigned int length = 0;
  for ( ; it != end; ++it)
    {
      AssertIndexRange(it->second, constraints.size());
      constraints[it->second] = &it->first;
      length += it->first.size();
    }
  constraint_pool_data.clear();
  constraint_pool_data.reserve (length);
  constraint_pool_row_index.reserve(constraint_values.constraints.size()+1);
  constraint_pool_row_index.resize(1, 0);
  for (unsigned int i=0; i<constraints.size(); ++i)
    {
      Assert(constraints[i] != nullptr, ExcInternalError());
      constraint_pool_data.insert(constraint_pool_data.end(),
                                  constraints[i]->begin(),
                                  constraints[i]->end());
      constraint_pool_row_index.push_back(constraint_pool_data.size());
    }
  AssertDimension(constraint_pool_data.size(), length);
  for (unsigned int no=0; no<n_fe; ++no)
    dof_info[no].reorder_cells(task_info, renumbering,
                               constraint_pool_row_index,
                               irregular_cells, vectorization_length);

  indices_are_initialized = true;
}



template <int dim, typename Number>
void MatrixFree<dim,Number>::clear()
{
  dof_info.clear();
  mapping_info.clear();
  cell_level_index.clear();
  task_info.clear();
  dof_handlers.dof_handler.clear();
  dof_handlers.hp_dof_handler.clear();
  indices_are_initialized = false;
  mapping_is_initialized  = false;
}



template <int dim, typename Number>
std::size_t MatrixFree<dim,Number>::memory_consumption () const
{
  std::size_t memory = MemoryConsumption::memory_consumption (dof_info);
  memory += MemoryConsumption::memory_consumption (cell_level_index);
  memory += MemoryConsumption::memory_consumption (shape_info);
  memory += MemoryConsumption::memory_consumption (constraint_pool_data);
  memory += MemoryConsumption::memory_consumption (constraint_pool_row_index);
  memory += MemoryConsumption::memory_consumption (task_info);
  memory += sizeof(*this);
  memory += mapping_info.memory_consumption();
  return memory;
}


template <int dim, typename Number>
template <typename StreamType>
void MatrixFree<dim,Number>::print_memory_consumption (StreamType &out) const
{
  out << "  Memory cell FE operator total: --> ";
  task_info.print_memory_statistics (out, memory_consumption());
  out << "   Memory cell index:                ";
  task_info.print_memory_statistics
  (out, MemoryConsumption::memory_consumption (cell_level_index));
  for (unsigned int j=0; j<dof_info.size(); ++ j)
    {
      out << "   Memory DoFInfo component "<< j << std::endl;
      dof_info[j].print_memory_consumption(out, task_info);
    }

  out << "   Memory mapping info" << std::endl;
  mapping_info.print_memory_consumption(out, task_info);

  out << "   Memory unit cell shape data:      ";
  task_info.print_memory_statistics
  (out, MemoryConsumption::memory_consumption (shape_info));
  if (task_info.scheme != internal::MatrixFreeFunctions::TaskInfo::none)
    {
      out << "   Memory task partitioning info:    ";
      task_info.print_memory_statistics
      (out, MemoryConsumption::memory_consumption (task_info));
    }
}



template <int dim, typename Number>
void MatrixFree<dim,Number>::print (std::ostream &out) const
{
  // print indices local to global
  for (unsigned int no=0; no<dof_info.size(); ++no)
    {
      out << "\n-- Index data for component " << no << " --" << std::endl;
      dof_info[no].print (constraint_pool_data, constraint_pool_row_index, out);
      out << std::endl;
    }
}



DEAL_II_NAMESPACE_CLOSE

#endif
