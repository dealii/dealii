// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2014 by the deal.II authors
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


DEAL_II_NAMESPACE_OPEN


// --------------------- MatrixFree -----------------------------------

template <int dim, typename Number>
MatrixFree<dim, Number>::MatrixFree()
  :
  indices_are_initialized (false),
  mapping_is_initialized  (false)
{}



template <int dim, typename Number>
MatrixFree<dim,Number>::~MatrixFree()
{}



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
  size_info = v.size_info;
  indices_are_initialized = v.indices_are_initialized;
  mapping_is_initialized  = v.mapping_is_initialized;
}



namespace internal
{
  template <int dim>
  void assert_communicator_equality (const dealii::Triangulation<dim> &tria,
                                     const MPI_Comm                  &comm_mf)
  {
#ifdef DEAL_II_WITH_MPI
    const parallel::distributed::Triangulation<dim> *dist_tria =
      dynamic_cast<const parallel::distributed::Triangulation<dim>*>(&tria);
    if (dist_tria != 0)
      {
        if (Utilities::System::job_supports_mpi())
          {
            int communicators_same = 0;
            MPI_Comm_compare (dist_tria->get_communicator(), comm_mf,
                              &communicators_same);
            Assert (communicators_same == MPI_IDENT ||
                    communicators_same == MPI_CONGRUENT,
                    ExcMessage ("MPI communicator in parallel::distributed::Triangulation "
                                "and matrix free class must be the same!"));
          }
      }
#else
    (void)tria;
    (void)comm_mf;
#endif
  }
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
      internal::assert_communicator_equality (dof_handler[0]->get_tria(),
                                              additional_data.mpi_communicator);
      size_info.communicator = additional_data.mpi_communicator;
      if (Utilities::System::job_supports_mpi() == true)
        {
          size_info.my_pid  =
            Utilities::MPI::this_mpi_process(size_info.communicator);
          size_info.n_procs =
            Utilities::MPI::n_mpi_processes(size_info.communicator);
        }
      else
        {
          size_info.my_pid = 0;
          size_info.n_procs = 1;
        }

      initialize_dof_handlers (dof_handler, additional_data.level_mg_handler);
      for (unsigned int no=0; no<dof_handler.size(); ++no)
        dof_info[no].store_plain_indices = additional_data.store_plain_indices;

      // initialize the basic multithreading information that needs to be
      // passed to the DoFInfo structure
#ifdef DEAL_II_WITH_THREADS
      if (additional_data.tasks_parallel_scheme != AdditionalData::none &&
          multithread_info.n_threads() > 1)
        {
          task_info.use_multithreading = true;
          task_info.block_size = additional_data.tasks_block_size;
          task_info.use_partition_partition =
            (additional_data.tasks_parallel_scheme ==
             AdditionalData::partition_partition ? true : false);
          task_info.use_coloring_only =
            (additional_data.tasks_parallel_scheme ==
             AdditionalData::color ? true : false);
        }
      else
#endif
        task_info.use_multithreading = false;

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
      size_info.make_layout (cell_level_index.size(),
                             VectorizedArray<Number>::n_array_elements,
                             dummy, dummy);
      for (unsigned int i=0; i<dof_info.size(); ++i)
        {
          dof_info[i].dimension    = dim;
          dof_info[i].n_components = dof_handler[i]->get_fe().element_multiplicity(0);
          dof_info[i].dofs_per_cell.push_back(dof_handler[i]->get_fe().dofs_per_cell);
          dof_info[i].row_starts.resize(size_info.n_macro_cells+1);
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
      mapping_info.initialize (dof_handler[0]->get_tria(), cell_level_index,
                               dof_info[0].cell_active_fe_index, mapping, quad,
                               additional_data.mapping_update_flags);

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
                                     dof_handler[i]->get_fe().size());
    unsigned int n_quad_in_collection = 0;
    for (unsigned int q=0; q<n_quad; ++q)
      n_quad_in_collection = std::max (n_quad_in_collection, quad[q].size());
    shape_info.reinit (TableIndices<4>(n_components, n_quad,
                                       n_fe_in_collection,
                                       n_quad_in_collection));
    for (unsigned int no=0; no<n_components; no++)
      for (unsigned int fe_no=0; fe_no<dof_handler[no]->get_fe().size(); ++fe_no)
        for (unsigned int nq =0; nq<n_quad; nq++)
          for (unsigned int q_no=0; q_no<quad[nq].size(); ++q_no)
            shape_info(no,nq,fe_no,q_no).reinit (quad[nq][q_no],
                                                 dof_handler[no]->get_fe()[fe_no]);
  }

  if (additional_data.initialize_indices == true)
    {
      clear();
      Assert (dof_handler.size() > 0, ExcMessage("No DoFHandler is given."));
      AssertDimension (dof_handler.size(), constraint.size());
      AssertDimension (dof_handler.size(), locally_owned_set.size());

      // set variables that are independent of FE
      internal::assert_communicator_equality (dof_handler[0]->get_tria(),
                                              additional_data.mpi_communicator);
      size_info.communicator = additional_data.mpi_communicator;
      if (Utilities::System::job_supports_mpi() == true)
        {
          size_info.my_pid  =
            Utilities::MPI::this_mpi_process(size_info.communicator);
          size_info.n_procs =
            Utilities::MPI::n_mpi_processes(size_info.communicator);
        }
      else
        {
          size_info.my_pid = 0;
          size_info.n_procs = 1;
        }

      initialize_dof_handlers (dof_handler, additional_data.level_mg_handler);
      for (unsigned int no=0; no<dof_handler.size(); ++no)
        dof_info[no].store_plain_indices = additional_data.store_plain_indices;

      // initialize the basic multithreading information that needs to be
      // passed to the DoFInfo structure
#ifdef DEAL_II_WITH_THREADS
      if (additional_data.tasks_parallel_scheme != AdditionalData::none &&
          multithread_info.n_threads() > 1)
        {
          task_info.use_multithreading = true;
          task_info.block_size = additional_data.tasks_block_size;
          task_info.use_partition_partition =
            (additional_data.tasks_parallel_scheme ==
             AdditionalData::partition_partition ? true : false);
          task_info.use_coloring_only =
            (additional_data.tasks_parallel_scheme ==
             AdditionalData::color ? true : false);
        }
      else
#endif
        task_info.use_multithreading = false;

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
      size_info.make_layout (cell_level_index.size(),
                             VectorizedArray<Number>::n_array_elements,
                             dummy, dummy);
      for (unsigned int i=0; i<dof_info.size(); ++i)
        {
          Assert(dof_handler[i]->get_fe().size() == 1, ExcNotImplemented());
          dof_info[i].dimension    = dim;
          dof_info[i].n_components = dof_handler[i]->get_fe()[0].element_multiplicity(0);
          dof_info[i].dofs_per_cell.push_back(dof_handler[i]->get_fe()[0].dofs_per_cell);
          dof_info[i].row_starts.resize(size_info.n_macro_cells+1);
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
      mapping_info.initialize (dof_handler[0]->get_tria(), cell_level_index,
                               dof_info[0].cell_active_fe_index, mapping, quad,
                               additional_data.mapping_update_flags);

      mapping_is_initialized = true;
    }
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
    else if (cell->subdomain_id() == subdomain_id)
      {
        Assert (cell->active(), ExcInternalError());
        cell_its.push_back (std::pair<unsigned int,unsigned int>
                            (cell->level(), cell->index()));
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

  // go through cells on zeroth level and then successively step down into
  // children. This gives a z-ordering of the cells, which is beneficial when
  // setting up neighboring relations between cells for thread parallelization
  const unsigned int n_mpi_procs = size_info.n_procs;
  const unsigned int my_pid = size_info.my_pid;

  const Triangulation<dim> &tria = dof_handlers.dof_handler[0]->get_tria();
  if (level == numbers::invalid_unsigned_int)
    {
      if (n_mpi_procs == 1)
        cell_level_index.reserve (tria.n_active_cells());
      typename Triangulation<dim>::cell_iterator cell = tria.begin(0),
                                                 end_cell = tria.end(0);
      for ( ; cell != end_cell; ++cell)
        internal::resolve_cell (cell, cell_level_index, my_pid);
    }
  else
    {
      AssertIndexRange (level, tria.n_levels());
      cell_level_index.reserve (tria.n_cells(level));
      typename Triangulation<dim>::cell_iterator cell = tria.begin(level),
                                                 end_cell = tria.end(level);
      for ( ; cell != end_cell; ++cell)
        cell_level_index.push_back (std::pair<unsigned int,unsigned int>
                                    (cell->level(), cell->index()));
    }
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
  const unsigned int n_mpi_procs = size_info.n_procs;
  const unsigned int my_pid = size_info.my_pid;

  // if we have no level given, use the same as for the standard DoFHandler,
  // otherwise we must loop through the respective level
  const Triangulation<dim> &tria = dof_handler[0]->get_tria();

  if (n_mpi_procs == 1)
    {
      cell_level_index.reserve (tria.n_active_cells());
    }
  typename hp::DoFHandler<dim>::cell_iterator cell = dof_handler[0]->begin(0),
                                              end_cell = dof_handler[0]->end(0);
  for ( ; cell != end_cell; ++cell)
    internal::resolve_cell (cell, cell_level_index,
                            my_pid);
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
  std::vector<std::vector<std::vector<unsigned int> > > lexicographic_inv(n_fe);

  internal::MatrixFreeFunctions::ConstraintValues<double> constraint_values;
  std::vector<unsigned int> constraint_indices;

  for (unsigned int no=0; no<n_fe; ++no)
    {
      std::vector<const FiniteElement<dim>*> fes;
      if (dof_handlers.active_dof_handler == DoFHandlers::hp)
        {
          const hp::DoFHandler<dim> *hpdof = dof_handlers.hp_dof_handler[no];
          const hp::FECollection<dim> &fe = hpdof->get_fe();
          for (unsigned int f=0; f<fe.size(); ++f)
            fes.push_back (&fe[f]);

          dof_info[no].max_fe_index = fe.size();
          dof_info[no].fe_index_conversion.resize (fe.size());
          for (unsigned int ind=0; ind<hpdof->get_fe().size(); ++ind)
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
      (new Utilities::MPI::Partitioner(locally_owned_set[no], size_info.communicator));

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
              cell_it (&dofh->get_tria(),
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
          // ok, now we are requested to use a level in a MGDoFHandler
          else if (dof_handlers.active_dof_handler == DoFHandlers::usual &&
                   dof_handlers.level != numbers::invalid_unsigned_int)
            {
              const DoFHandler<dim> *dofh = dof_handlers.dof_handler[no];
              AssertIndexRange (dof_handlers.level, dofh->get_tria().n_levels());
              typename DoFHandler<dim>::cell_iterator
              cell_it (&dofh->get_tria(),
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
              cell_it (&dofh->get_tria(),
                       cell_level_index[counter].first,
                       cell_level_index[counter].second,
                       dofh);
              if (dofh->get_fe().size() > 1)
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
      if (cell_at_boundary == true)
        boundary_cells.push_back(counter);
    }

  const unsigned int vectorization_length =
    VectorizedArray<Number>::n_array_elements;
  std::vector<unsigned int> irregular_cells;
  size_info.make_layout (n_active_cells, vectorization_length, boundary_cells,
                         irregular_cells);

  for (unsigned int no=0; no<n_fe; ++no)
    dof_info[no].assign_ghosts (boundary_cells);

  // reorganize the indices in order to overlap communication in MPI with
  // computations: Place all cells with ghost indices into one chunk. Also
  // reorder cells so that we can parallelize by threads
  std::vector<unsigned int> renumbering;
  if (task_info.use_multithreading == true)
    {
      dof_info[0].compute_renumber_parallel (boundary_cells, size_info,
                                             renumbering);
      if (task_info.use_partition_partition == true)
        dof_info[0].make_thread_graph_partition_partition
        (size_info, task_info, renumbering, irregular_cells,
         dof_handlers.active_dof_handler == DoFHandlers::hp);
      else
        dof_info[0].make_thread_graph_partition_color
        (size_info, task_info, renumbering, irregular_cells,
         dof_handlers.active_dof_handler == DoFHandlers::hp);
    }
  else
    {
      // In case, we have an hp-dofhandler, we have to reorder the cell
      // according to the polynomial degree on the cell.
      dof_info[0].compute_renumber_serial (boundary_cells, size_info,
                                           renumbering);
      if (dof_handlers.active_dof_handler == DoFHandlers::hp)
        dof_info[0].compute_renumber_hp_serial (size_info, renumbering,
                                                irregular_cells);
    }

  // Finally perform the renumbering. We also want to group several cells
  // together to one "macro-cell" for vectorization (where the arithmetic
  // operations will then be done simultaneously).
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
    cell_level_index.reserve(size_info.n_macro_cells*vectorization_length);
    unsigned int position_cell=0;
    for (unsigned int i=0; i<size_info.n_macro_cells; ++i)
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
    AssertDimension (position_cell, size_info.n_active_cells);
    AssertDimension (cell_level_index.size(),size_info.n_macro_cells*vectorization_length);
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
      Assert(constraints[i] != 0, ExcInternalError());
      constraint_pool_data.insert(constraint_pool_data.end(),
                                  constraints[i]->begin(),
                                  constraints[i]->end());
      constraint_pool_row_index.push_back(constraint_pool_data.size());
    }
  AssertDimension(constraint_pool_data.size(), length);
  for (unsigned int no=0; no<n_fe; ++no)
    dof_info[no].reorder_cells(size_info, renumbering,
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
  size_info.clear();
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
template <typename STREAM>
void MatrixFree<dim,Number>::print_memory_consumption (STREAM &out) const
{
  out << "  Memory cell FE operator total: --> ";
  size_info.print_memory_statistics (out, memory_consumption());
  out << "   Memory cell index:                ";
  size_info.print_memory_statistics
  (out, MemoryConsumption::memory_consumption (cell_level_index));
  for (unsigned int j=0; j<dof_info.size(); ++ j)
    {
      out << "   Memory DoFInfo component "<< j << std::endl;
      dof_info[j].print_memory_consumption(out, size_info);
    }

  out << "   Memory mapping info" << std::endl;
  mapping_info.print_memory_consumption(out, size_info);

  out << "   Memory unit cell shape data:      ";
  size_info.print_memory_statistics
  (out, MemoryConsumption::memory_consumption (shape_info));
  if (task_info.use_multithreading == true)
    {
      out << "   Memory task partitioning info:    ";
      size_info.print_memory_statistics
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



/*-------------------- Implementation of helper functions ------------------*/

namespace internal
{
  namespace MatrixFreeFunctions
  {

    TaskInfo::TaskInfo ()
    {
      clear();
    }



    void TaskInfo::clear ()
    {
      block_size = 0;
      n_blocks = 0;
      block_size_last = 0;
      position_short_block = 0;
      use_multithreading = false;
      use_partition_partition = false;
      use_coloring_only = false;
      partition_color_blocks_row_index.clear();
      partition_color_blocks_data.clear();
      evens = 0;
      odds = 0;
      n_blocked_workers = 0;
      n_workers = 0;
      partition_evens.clear();
      partition_odds.clear();
      partition_n_blocked_workers.clear();
      partition_n_workers.clear();
    }



    std::size_t
    TaskInfo::memory_consumption () const
    {
      return (sizeof(*this)+
              MemoryConsumption::memory_consumption (partition_color_blocks_row_index) +
              MemoryConsumption::memory_consumption (partition_color_blocks_data)+
              MemoryConsumption::memory_consumption (partition_evens) +
              MemoryConsumption::memory_consumption (partition_odds) +
              MemoryConsumption::memory_consumption (partition_n_blocked_workers) +
              MemoryConsumption::memory_consumption (partition_n_workers));
    }



    SizeInfo::SizeInfo ()
    {
      clear();
    }



    void SizeInfo::clear()
    {
      n_active_cells = 0;
      n_macro_cells  = 0;
      boundary_cells_start = 0;
      boundary_cells_end   = 0;
      vectorization_length = 0;
      locally_owned_cells  = IndexSet();
      ghost_cells = IndexSet();
      communicator = MPI_COMM_SELF;
      my_pid = 0;
      n_procs = 0;
    }



    template <typename STREAM>
    void SizeInfo::print_memory_statistics (STREAM     &out,
                                            std::size_t data_length) const
    {
      Utilities::MPI::MinMaxAvg memory_c;
      if (Utilities::System::job_supports_mpi() == true)
        {
          memory_c = Utilities::MPI::min_max_avg (1e-6*data_length,
                                                  communicator);
        }
      else
        {
          memory_c.sum = 1e-6*data_length;
          memory_c.min = memory_c.sum;
          memory_c.max = memory_c.sum;
          memory_c.avg = memory_c.sum;
          memory_c.min_index = 0;
          memory_c.max_index = 0;
        }
      if (n_procs < 2)
        out << memory_c.min;
      else
        out << memory_c.min << "/" << memory_c.avg << "/" << memory_c.max;
      out << " MB" << std::endl;
    }



    inline
    void SizeInfo::make_layout (const unsigned int n_active_cells_in,
                                const unsigned int vectorization_length_in,
                                std::vector<unsigned int> &boundary_cells,
                                std::vector<unsigned int> &irregular_cells)
    {
      vectorization_length = vectorization_length_in;
      n_active_cells = n_active_cells_in;

      unsigned int n_max_boundary_cells = boundary_cells.size();
      unsigned int n_boundary_cells = n_max_boundary_cells;

      // try to make the number of boundary cells divisible by the number of
      // vectors in vectorization

      /*
      // try to balance the number of cells before and after the boundary part
      // on each processor. probably not worth it!
      #ifdef DEAL_II_WITH_MPI
      MPI_Allreduce (&n_boundary_cells, &n_max_boundary_cells, 1, MPI_UNSIGNED,
                     MPI_MAX, size_info.communicator);
      #endif
      if (n_max_boundary_cells > n_active_cells)
        n_max_boundary_cells = n_active_cells;
      */

      unsigned int fillup_needed =
        (vectorization_length - n_boundary_cells%vectorization_length)%vectorization_length;
      if (fillup_needed > 0 && n_boundary_cells < n_active_cells)
        {
          // fill additional cells into the list of boundary cells to get a
          // balanced number. Go through the indices successively until we
          // found enough indices
          std::vector<unsigned int> new_boundary_cells;
          new_boundary_cells.reserve (n_max_boundary_cells);

          unsigned int next_free_slot = 0, bound_index = 0;
          while (fillup_needed > 0 && bound_index < boundary_cells.size())
            {
              if (next_free_slot < boundary_cells[bound_index])
                {
                  // check if there are enough cells to fill with in the
                  // current slot
                  if (next_free_slot + fillup_needed <= boundary_cells[bound_index])
                    {
                      for (unsigned int j=boundary_cells[bound_index]-fillup_needed;
                           j < boundary_cells[bound_index]; ++j)
                        new_boundary_cells.push_back(j);
                      fillup_needed = 0;
                    }
                  // ok, not enough indices, so just take them all up to the
                  // next boundary cell
                  else
                    {
                      for (unsigned int j=next_free_slot;
                           j<boundary_cells[bound_index]; ++j)
                        new_boundary_cells.push_back(j);
                      fillup_needed -= boundary_cells[bound_index]-next_free_slot;
                    }
                }
              new_boundary_cells.push_back(boundary_cells[bound_index]);
              next_free_slot = boundary_cells[bound_index]+1;
              ++bound_index;
            }
          while (fillup_needed > 0 && (new_boundary_cells.size()==0 ||
                                       new_boundary_cells.back()<n_active_cells-1))
            new_boundary_cells.push_back(new_boundary_cells.back()+1);
          while (bound_index<boundary_cells.size())
            new_boundary_cells.push_back(boundary_cells[bound_index++]);

          boundary_cells.swap(new_boundary_cells);
        }

      // set the number of cells
      std::sort (boundary_cells.begin(), boundary_cells.end());
      n_boundary_cells = boundary_cells.size();

      // check that number of boundary cells is divisible by
      // vectorization_length or that it contains all cells
      Assert (n_boundary_cells % vectorization_length == 0 ||
              n_boundary_cells == n_active_cells, ExcInternalError());
      n_macro_cells = (n_active_cells+vectorization_length-1)/vectorization_length;
      irregular_cells.resize (n_macro_cells);
      if (n_macro_cells*vectorization_length > n_active_cells)
        {
          irregular_cells[n_macro_cells-1] =
            vectorization_length - (n_macro_cells*vectorization_length - n_active_cells);
        }
      if (n_procs > 1)
        {
          const unsigned int n_macro_boundary_cells =
            (n_boundary_cells+vectorization_length-1)/vectorization_length;
          boundary_cells_start = (n_macro_cells-n_macro_boundary_cells)/2;
          boundary_cells_end   = boundary_cells_start + n_macro_boundary_cells;
        }
      else
        boundary_cells_start = boundary_cells_end = n_macro_cells;
    }



    /* ------------------------------------------------------------------ */

    template <typename Number>
    FPArrayComparator<Number>::FPArrayComparator (const Number scaling)
      :
      tolerance (scaling * std::numeric_limits<double>::epsilon() * 1024.)
    {}



    template <typename Number>
    bool
    FPArrayComparator<Number>::operator() (const std::vector<Number> &v1,
                                           const std::vector<Number> &v2) const
    {
      const unsigned int s1 = v1.size(), s2 = v2.size();
      if (s1 < s2)
        return true;
      else if (s1 > s2)
        return false;
      else
        for (unsigned int i=0; i<s1; ++i)
          if (v1[i] < v2[i] - tolerance)
            return true;
          else if (v1[i] > v2[i] + tolerance)
            return false;
      return false;
    }



    template <typename Number>
    template <int dim>
    bool
    FPArrayComparator<Number>::
    operator ()(const Tensor<1,dim,Tensor<1,VectorizedArray<Number>::n_array_elements,Number> > &t1,
                const Tensor<1,dim,Tensor<1,VectorizedArray<Number>::n_array_elements,Number> > &t2) const
    {
      for (unsigned int d=0; d<dim; ++d)
        for (unsigned int k=0; k<VectorizedArray<Number>::n_array_elements; ++k)
          if ((t1)[d][k] < (t2)[d][k] - tolerance)
            return true;
          else if ((t1)[d][k] > (t2)[d][k] + tolerance)
            return false;
      return false;
    }



    template <typename Number>
    template <int dim>
    bool
    FPArrayComparator<Number>::
    operator ()(const Tensor<2,dim,Tensor<1,VectorizedArray<Number>::n_array_elements,Number> > &t1,
                const Tensor<2,dim,Tensor<1,VectorizedArray<Number>::n_array_elements,Number> > &t2) const
    {
      for (unsigned int d=0; d<dim; ++d)
        for (unsigned int e=0; e<dim; ++e)
          for (unsigned int k=0; k<VectorizedArray<Number>::n_array_elements; ++k)
            if ((t1)[d][e][k] < (t2)[d][e][k] - tolerance)
              return true;
            else if ((t1)[d][e][k] > (t2)[d][e][k] + tolerance)
              return false;
      return false;
    }
  }
}


DEAL_II_NAMESPACE_CLOSE
