//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <deal.II/base/utilities.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/tensor_product_polynomials.h>
#include <deal.II/base/mpi.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/multigrid/mg_dof_accessor.h>
#include <deal.II/fe/fe_poly.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/fe_evaluation_data.templates.h>
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
  constraint_pool = v.constraint_pool;
  mapping_info = v.mapping_info;
  fe_evaluation_data = v.fe_evaluation_data;
  cell_level_index = v.cell_level_index;
  task_info = v.task_info;
  size_info = v.size_info;
  indices_are_initialized = v.indices_are_initialized;
  mapping_is_initialized  = v.mapping_is_initialized;
}



namespace internal
{
  template <int dim>
  void assert_communicator_equality (const dealii::Triangulation<dim>&tria,
                                     const MPI_Comm                  &comm_mf)
  {
#if DEAL_II_COMPILER_SUPPORTS_MPI
    const parallel::distributed::Triangulation<dim>* dist_tria =
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
template <typename DH>
void MatrixFree<dim,Number>::
internal_reinit(const Mapping<dim>                         &mapping,
                const std::vector<const DH*>               &dof_handler,
                const std::vector<const ConstraintMatrix*> &constraint,
                const std::vector<IndexSet>                &locally_owned_set,
                const std::vector<hp::QCollection<1> >     &quad,
                const MatrixFree<dim,Number>::AdditionalData additional_data)
{
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

                                // initialize the basic multithreading
                                // information that needs to be passed to the
                                // DoFInfo structure
#if DEAL_II_USE_MT == 1
      if (additional_data.tasks_parallel_scheme != AdditionalData::none)
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

                                // set dof_indices together with
                                // constraint_indicator and
                                // constraint_pool. It also reorders the way
                                // cells are gone through (to separate cells
                                // with overlap to other processors from
                                // others without).
      initialize_indices (constraint, locally_owned_set);
    }

                                // Reads out the FE information and stores the
                                // shape function values, gradients and
                                // Hessians for quadrature points.
  const unsigned int n_fe   = dof_handler.size();
  const unsigned int n_quad = quad.size();
  fe_evaluation_data.reinit (TableIndices<4>(n_fe, n_quad, 1, 1));
  for (unsigned int no=0; no<n_fe; no++)
    {
      const FiniteElement<dim> &fe = dof_handler[no]->get_fe();
      for(unsigned int nq =0;nq<n_quad;nq++)
        {
          AssertDimension (quad[nq].size(), 1);
          fe_evaluation_data(no,nq,0,0).reinit(quad[nq][0], fe.base_element(0));
        }
    }

                                // Evaluates transformations from unit to real
                                // cell, Jacobian determinants, quadrature
                                // points in real space, based on the ordering
                                // of the cells determined in @p
                                // extract_local_to_global_indices.
  if(additional_data.initialize_mapping == true)
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
                const std::vector<const ConstraintMatrix*>    &constraint,
                const std::vector<IndexSet>                   &locally_owned_set,
                const std::vector<hp::QCollection<1> >        &quad,
                const MatrixFree<dim,Number>::AdditionalData additional_data)
{
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

                                // initialize the basic multithreading
                                // information that needs to be passed to the
                                // DoFInfo structure
#if DEAL_II_USE_MT == 1
      if (additional_data.tasks_parallel_scheme != AdditionalData::none)
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

                                // set dof_indices together with
                                // constraint_indicator and
                                // constraint_pool. It also reorders the way
                                // cells are gone through (to separate cells
                                // with overlap to other processors from
                                // others without).
      initialize_indices (constraint, locally_owned_set);
    }

                                // Reads out the FE information and stores the
                                // shape function values, gradients and
                                // Hessians for quadrature points.
  const unsigned int n_components = dof_handler.size();
  const unsigned int n_quad       = quad.size();
  unsigned int n_fe_in_collection = 0;
  for (unsigned int i=0; i<n_components; ++i)
    n_fe_in_collection = std::max (n_fe_in_collection,
                                   dof_handler[i]->get_fe().size());
  unsigned int n_quad_in_collection = 0;
  for (unsigned int q=0; q<n_quad; ++q)
    n_quad_in_collection = std::max (n_quad_in_collection, quad[q].size());
  fe_evaluation_data.reinit (TableIndices<4>(n_components, n_quad,
                                             n_fe_in_collection, 
                                             n_quad_in_collection));
  for (unsigned int no=0; no<n_components; no++)
    for (unsigned int fe_no=0; fe_no<dof_handler[no]->get_fe().size(); ++fe_no)
      {
        const FiniteElement<dim> &fe = dof_handler[no]->get_fe()[fe_no];
        for(unsigned int nq =0; nq<n_quad; nq++)
          for (unsigned int q_no=0; q_no<quad[nq].size(); ++q_no)
            fe_evaluation_data(no,nq,fe_no,q_no).reinit (quad[nq][q_no],
                                                         fe.base_element(0));
      }

                                // Evaluates transformations from unit to real
                                // cell, Jacobian determinants, quadrature
                                // points in real space, based on the ordering
                                // of the cells determined in @p
                                // extract_local_to_global_indices.
  if(additional_data.initialize_mapping == true)
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
                     std::vector<std::pair<unsigned int,unsigned int> >&cell_its,
                     const unsigned int  subdomain_id)
  {
    if (cell->has_children())
      for (unsigned int child=0; child<cell->n_children(); ++child)
        resolve_cell (cell->child(child), cell_its,
                      subdomain_id);
    else
      if (cell->subdomain_id() == subdomain_id)
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
                         const unsigned int)
{
  dof_handlers.active_dof_handler = DoFHandlers::usual;
  dof_handlers.n_dof_handlers = dof_handler.size();
  dof_handlers.dof_handler.resize (dof_handlers.n_dof_handlers);
  for (unsigned int no=0; no<dof_handlers.n_dof_handlers; ++no)
    dof_handlers.dof_handler[no] = dof_handler[no];

  dof_info.resize (dof_handlers.n_dof_handlers);

                                // go through cells on zeroth level and then
                                // successively step down into children. This
                                // gives a z-ordering of the cells, which is
                                // beneficial when setting up neighboring
                                // relations between cells for thread
                                // parallelization
  const unsigned int n_mpi_procs = size_info.n_procs;
  const unsigned int my_pid = size_info.my_pid;

  const Triangulation<dim> &tria = dof_handlers.dof_handler[0]->get_tria();
  {
    if (n_mpi_procs == 1)
      cell_level_index.reserve (tria.n_active_cells());
    typename Triangulation<dim>::cell_iterator cell = tria.begin(0),
      end_cell = tria.end(0);
    for ( ; cell != end_cell; ++cell)
      internal::resolve_cell (cell, cell_level_index, my_pid);
  }
}



template <int dim, typename Number>
void MatrixFree<dim,Number>::
initialize_dof_handlers (const std::vector<const MGDoFHandler<dim>*> &dof_handler,
                         const unsigned int                           level)
{
  dof_handlers.active_dof_handler = DoFHandlers::multigrid;
  dof_handlers.level = level;
  dof_handlers.n_dof_handlers = dof_handler.size();
  dof_handlers.mg_dof_handler.resize (dof_handlers.n_dof_handlers);
  for (unsigned int no=0; no<dof_handlers.n_dof_handlers; ++no)
    dof_handlers.mg_dof_handler[no] = dof_handler[no];

  dof_info.resize (dof_handlers.n_dof_handlers);

                                // go through cells on zeroth level and then
                                // successively step down into children. This
                                // gives a z-ordering of the cells, which is
                                // beneficial when setting up neighboring
                                // relations between cells for thread
                                // parallelization
  const unsigned int n_mpi_procs = size_info.n_procs;
  const unsigned int my_pid = size_info.my_pid;

                                // if we have no level given, use the same as
                                // for the standard DoFHandler, otherwise we
                                // must loop through the respective level
  const Triangulation<dim> &tria = dof_handlers.mg_dof_handler[0]->get_tria();

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

                                // go through cells on zeroth level and then
                                // successively step down into children. This
                                // gives a z-ordering of the cells, which is
                                // beneficial when setting up neighboring
                                // relations between cells for thread
                                // parallelization
  const unsigned int n_mpi_procs = size_info.n_procs;
  const unsigned int my_pid = size_info.my_pid;

                                // if we have no level given, use the same as
                                // for the standard DoFHandler, otherwise we
                                // must loop through the respective level
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
(const std::vector<const ConstraintMatrix*> &constraint,
 const std::vector<IndexSet>                &locally_owned_set)
{
  const unsigned int n_fe = dof_handlers.n_dof_handlers;
  const unsigned int n_active_cells = cell_level_index.size();
  const unsigned int n_mpi_procs = size_info.n_procs;

  AssertDimension (n_active_cells, cell_level_index.size());
  AssertDimension (n_fe, locally_owned_set.size());
  AssertDimension (n_fe, constraint.size());

  std::vector<unsigned int> local_dof_indices;
  std::vector<std::vector<unsigned int> > ghost_dofs(n_fe);
  std::vector<std::vector<std::vector<unsigned int> > > lexicographic_inv(n_fe);

  internal::MatrixFreeFunctions::internal::ConstraintValues<double> constraint_values;
  std::vector<unsigned int> constraint_indices;

  for(unsigned int no=0; no<n_fe; ++no)
    {
      std::vector<const FiniteElement<dim>*> fes;
      if (dof_handlers.active_dof_handler == DoFHandlers::hp)
        {
          const hp::DoFHandler<dim> *hpdof = dof_handlers.hp_dof_handler[no];
          const hp::FECollection<dim> &fe = hpdof->get_fe();
          for (unsigned int f=0; f<fe.size(); ++f)
            fes.push_back (&fe[f]);

          dof_info[no].cell_active_fe_index.resize(n_active_cells,
                                                   numbers::invalid_unsigned_int);
          dof_info[no].max_fe_index = fe.size();
          dof_info[no].fe_index_conversion.resize (fe.size());
          for (unsigned int ind=0; ind<hpdof->get_fe().size(); ++ind)
            dof_info[no].fe_index_conversion[ind] =
              std::pair<unsigned int,unsigned int>(fe[ind].degree,
                                                   fe[ind].dofs_per_cell);
        }
      else
        {
          const DoFHandler<dim> * dofh =
            dof_handlers.active_dof_handler == DoFHandlers::usual ?
            &*dof_handlers.dof_handler[no] : &*dof_handlers.mg_dof_handler[no];
          fes.push_back (&dofh->get_fe());
          dof_info[no].max_fe_index = 1;
        }

      lexicographic_inv[no].resize (fes.size());
      for (unsigned int fe_index = 0; fe_index<fes.size(); ++fe_index)
        {
          const FiniteElement<dim> &fe = *fes[fe_index];
          Assert (fe.n_base_elements() == 1,
                  ExcMessage ("MatrixFree only works for DoFHandler with one base element"));
          const unsigned int n_fe_components = fe.element_multiplicity (0);

          // cache number of finite elements and
          // dofs_per_cell
          dof_info[no].dofs_per_cell.push_back (fe.dofs_per_cell);
          dof_info[no].dofs_per_face.push_back (fe.dofs_per_face);
          dof_info[no].n_components  = n_fe_components;


          // get permutation that gives lexicographic
          // renumbering of the cell dofs
          // renumber (this is necessary for FE_Q, for
          // example, since there the vertex DoFs come
          // first, which is incompatible with the
          // lexicographic ordering necessary to apply
          // tensor products efficiently)
          const FE_Poly<TensorProductPolynomials<dim>,dim,dim> *cast_fe =
            dynamic_cast<const FE_Poly<TensorProductPolynomials<dim>,dim,dim>*>
            (&fe.base_element(0));
          // This class currently only works for
          // elements derived from
          // FE_Poly<TensorProductPolynomials<dim>,dim,dim>.
          // For any other element, the dynamic cast
          // above will fail and give cast_fe == 0.
          Assert (cast_fe != 0, ExcNotImplemented());

          // create a derived finite element that gives
          // us access to the inverse numbering (which
          // we need in order to get a lexicographic
          // ordering of local degrees of freedom)
          const internal::MatrixFreeFunctions::internal::FE_PolyAccess<dim,dim>&fe_acc =
            static_cast<const internal::MatrixFreeFunctions::internal::
            FE_PolyAccess<dim,dim> &>(*cast_fe);
          if (n_fe_components == 1)
            {
              lexicographic_inv[no][fe_index] = fe_acc.get_numbering_inverse();
              AssertDimension (lexicographic_inv[no][fe_index].size(),
                               dof_info[no].dofs_per_cell[fe_index]);
            }
          else
            {
              // ok, we have more than one component
              Assert (n_fe_components > 1, ExcInternalError());
              std::vector<unsigned int> scalar_lex=fe_acc.get_numbering();
              AssertDimension (scalar_lex.size() * n_fe_components,
                               dof_info[no].dofs_per_cell[fe_index]);
              lexicographic_inv[no][fe_index].resize (dof_info[no].dofs_per_cell[fe_index]);
              std::vector<unsigned int> lexicographic (dof_info[no].dofs_per_cell[fe_index]);
              for (unsigned int comp=0; comp<n_fe_components; ++comp)
                for (unsigned int i=0; i<scalar_lex.size(); ++i)
                  lexicographic[fe.component_to_system_index(comp,i)]
                    = scalar_lex.size () * comp + scalar_lex[i];

              // invert numbering
              for (unsigned int i=0; i<lexicographic.size(); ++i)
                lexicographic_inv[no][fe_index][lexicographic[i]] = i;

#ifdef DEBUG
              // check that we got a useful permutation
              lexicographic = lexicographic_inv[no][fe_index];
              std::sort(lexicographic.begin(), lexicographic.end());
              for (unsigned int i=0; i<lexicographic.size(); ++i)
                AssertDimension (lexicographic[i], i);
#endif
            }
          AssertDimension (lexicographic_inv[no][fe_index].size(),
                           dof_info[no].dofs_per_cell[fe_index]);
        }

                                // set locally owned range for each component
      Assert (locally_owned_set[no].is_contiguous(), ExcNotImplemented());
      dof_info[no].vector_partitioner.reset
        (new Utilities::MPI::Partitioner(locally_owned_set[no], size_info.communicator));

                                // initialize the arrays for indices
      dof_info[no].row_starts.resize (n_active_cells+1);
      dof_info[no].row_starts[0] = std_cxx1x::tuple<unsigned int,unsigned int,
                                                    unsigned int> (0,0,0);
      dof_info[no].dof_indices.reserve
        ((n_active_cells*dof_info[no].dofs_per_cell[0]*3)/2);

                                // cache the constrained indices for use in
                                // matrix-vector products
      {
        const unsigned int
          start_index = dof_info[no].vector_partitioner->local_range().first,
          end_index   = dof_info[no].vector_partitioner->local_range().second;
        for (unsigned int i=start_index; i<end_index; ++i)
          if (constraint[no]->is_constrained(i)==true)
            dof_info[no].constrained_dofs.push_back(i-start_index);
      }

      if (n_mpi_procs > 1)
        ghost_dofs[no].reserve (locally_owned_set[no].n_elements()/10+1);
    }

                                // extract all the global indices associated
                                // with the computation, and form the ghost
                                // indices
  std::vector<unsigned int> boundary_cells;
  for(unsigned int counter = 0 ; counter < n_active_cells ; ++counter)
    {
      bool cell_at_boundary = false;
      for (unsigned int no=0; no<n_fe; ++no)
        {
                                // OK, read indices from standard DoFHandler
                                // or active indices in MGDoFHandler. That is
                                // the usual stuff
          if (dof_handlers.active_dof_handler == DoFHandlers::usual ||
              (dof_handlers.active_dof_handler == DoFHandlers::multigrid &&
               dof_handlers.level == numbers::invalid_unsigned_int))
            {
              const DoFHandler<dim> * dofh =
                dof_handlers.active_dof_handler == DoFHandlers::usual ?
                &*dof_handlers.dof_handler[no] : &*dof_handlers.mg_dof_handler[no];
              typename DoFHandler<dim>::active_cell_iterator
                cell_it (&dofh->get_tria(),
                         cell_level_index[counter].first,
                         cell_level_index[counter].second,
                         dofh);
              local_dof_indices.resize (dof_info[no].dofs_per_cell[0]);
              cell_it->get_dof_indices(local_dof_indices);
              dof_info[no].read_dof_indices (local_dof_indices,
                                             lexicographic_inv[no][0],
                                             *constraint[no], counter,
                                             constraint_values,
                                             cell_at_boundary);
            }
                                // ok, now we are requested to use a level in
                                // a MGDoFHandler
          else if (dof_handlers.active_dof_handler == DoFHandlers::multigrid &&
                   dof_handlers.level != numbers::invalid_unsigned_int)
            {
              const MGDoFHandler<dim> * dofh =
                dof_handlers.mg_dof_handler[no];
              AssertIndexRange (dof_handlers.level, dofh->get_tria().n_levels());
              typename MGDoFHandler<dim>::cell_iterator
                cell_it (&dofh->get_tria(),
                         cell_level_index[counter].first,
                         cell_level_index[counter].second,
                         dofh);
              local_dof_indices.resize (dof_info[no].dofs_per_cell[0]);
              cell_it->get_mg_dof_indices(local_dof_indices);
              dof_info[no].read_dof_indices (local_dof_indices,
                                             lexicographic_inv[no][0],
                                             *constraint[no], counter,
                                             constraint_values,
                                             cell_at_boundary);
            }
          else if (dof_handlers.active_dof_handler == DoFHandlers::hp)
            {
              const hp::DoFHandler<dim> * dofh =
                dof_handlers.hp_dof_handler[no];
              typename hp::DoFHandler<dim>::active_cell_iterator
                cell_it (&dofh->get_tria(),
                         cell_level_index[counter].first,
                         cell_level_index[counter].second,
                         dofh);
              dof_info[no].cell_active_fe_index[counter] =
                cell_it->active_fe_index();
              local_dof_indices.resize (cell_it->get_fe().dofs_per_cell);
              cell_it->get_dof_indices(local_dof_indices);
              dof_info[no].read_dof_indices (local_dof_indices,
                                             lexicographic_inv[no][cell_it->active_fe_index()],
                                             *constraint[no], counter,
                                             constraint_values,
                                             cell_at_boundary);
            }
          else
            {
              Assert (false, ExcNotImplemented());
            }
        }

                                // if we found dofs on some FE component that
                                // belong to other processors, the cell is
                                // added to the boundary cells.
      if (cell_at_boundary == true)
        boundary_cells.push_back(counter);
    }

                                // try to make the number of boundary cells
                                // divisible by the number of vectors in
                                // vectorization
  const unsigned int n_vectors = VectorizedArray<Number>::n_array_elements;
  {
    unsigned int n_max_boundary_cells = boundary_cells.size();
    unsigned int n_boundary_cells = n_max_boundary_cells;

    /*
                                // try to balance the number of cells before
                                // and after the boundary part on each
                                // processor. probably not worth it!
#if DEAL_II_COMPILER_SUPPORTS_MPI
    MPI_Allreduce (&n_boundary_cells, &n_max_boundary_cells, 1, MPI_UNSIGNED,
                   MPI_MAX, size_info.communicator);
#endif
    if (n_max_boundary_cells > n_active_cells)
      n_max_boundary_cells = n_active_cells;
    */

    unsigned int fillup_needed =
      (n_vectors - n_boundary_cells%n_vectors)%n_vectors;
    /*
    if (task_info.use_multithreading == true)
      fillup_needed =
        (n_vectors - n_boundary_cells%n_vectors)%n_vectors;
    else
      fillup_needed = (n_max_boundary_cells +
                       (n_vectors - n_max_boundary_cells%n_vectors)%n_vectors -
                             n_boundary_cells);
    */
    if (fillup_needed > 0 && n_boundary_cells < n_active_cells)
      {
                                // fill additional cells into the list of
                                // boundary cells to get a balanced number. Go
                                // through the indices successively until we
                                // found enough indices
        std::vector<unsigned int> new_boundary_cells;
        new_boundary_cells.reserve (n_max_boundary_cells);

        unsigned int next_free_slot = 0, bound_index = 0;
        while (fillup_needed > 0 && bound_index < boundary_cells.size())
          {
            if (next_free_slot < boundary_cells[bound_index])
              {
                                // check if there are enough cells to fill
                                // with in the current slot
                if (next_free_slot + fillup_needed <= boundary_cells[bound_index])
                  {
                    for (unsigned int j=boundary_cells[bound_index]-fillup_needed;
                         j < boundary_cells[bound_index]; ++j)
                      new_boundary_cells.push_back(j);
                    fillup_needed = 0;
                  }
                                // ok, not enough indices, so just take them
                                // all up to the next boundary cell
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
  }

                                // set the number of cells
  const unsigned int n_boundary_cells = boundary_cells.size();
  std::sort (boundary_cells.begin(), boundary_cells.end());
  std::vector<unsigned int> irregular_cells;
  size_info.make_layout (n_active_cells, n_boundary_cells, n_vectors,
                         irregular_cells);

  for (unsigned int no=0; no<n_fe; ++no)
    dof_info[no].assign_ghosts (boundary_cells);

                                // reorganize the indices: we want to put the
                                // boundary cells at the beginning for
                                // multithreading. So just renumber the cell
                                // indices and put them at the beginning of
                                // the list that determines the order of the
                                // cells
  std::vector<unsigned int> renumbering (n_active_cells,
                                         numbers::invalid_unsigned_int);
  {
    std::vector<unsigned int> reverse_numbering (n_active_cells,
                                                 numbers::invalid_unsigned_int);
    unsigned int counter;
    if (task_info.use_multithreading == true)
      {
        for (unsigned int j=0; j<n_boundary_cells; ++j)
          reverse_numbering[boundary_cells[j]] = j;
        counter = n_boundary_cells;
        for (unsigned int j=0; j<n_active_cells; ++j)
          if (reverse_numbering[j] == numbers::invalid_unsigned_int)
            reverse_numbering[j] = counter++;

        size_info.boundary_cells_end   = (size_info.boundary_cells_end -
                                          size_info.boundary_cells_start);
        size_info.boundary_cells_start = 0;
      }
                                //  Otherwise, we put the boundary cells to
                                //  the middle.
    else
      {
        for (unsigned int j=0; j<n_boundary_cells; ++j)
          reverse_numbering[boundary_cells[j]] = j+n_vectors*size_info.boundary_cells_start;
        counter = 0;
        unsigned int j = 0;
        while (counter < n_active_cells &&
               counter < n_vectors * size_info.boundary_cells_start)
          {
            if (reverse_numbering[j] == numbers::invalid_unsigned_int)
              reverse_numbering[j] = counter++;
            j++;
          }
        counter = std::min (n_vectors*size_info.boundary_cells_start+n_boundary_cells,
                            n_active_cells);
        if (counter < n_active_cells)
          {
            for ( ; j<n_active_cells; ++j)
              if (reverse_numbering[j] == numbers::invalid_unsigned_int)
                reverse_numbering[j] = counter++;
          }
      }
    AssertDimension (counter, n_active_cells);
    for (unsigned int j=0; j<n_active_cells; ++j)
      {
        AssertIndexRange (reverse_numbering[j], n_active_cells);
        renumbering[reverse_numbering[j]] = j;
      }
  }

                                // reorder cells so that we can parallelize by
                                // threads
  if (task_info.use_multithreading == true)
    {
      if(task_info.use_partition_partition == true)
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
                                // In case, we have an hp-dofhandler, we have
                                // to reorder the cell according to the
                                // polynomial degree on the cell.
      if (dof_handlers.active_dof_handler == DoFHandlers::hp)
      {
        const unsigned int max_fe_index =
          dof_info[0].max_fe_index;
        irregular_cells.resize (0);
        irregular_cells.resize (size_info.n_macro_cells+3*max_fe_index);
        const std::vector<unsigned int> &cell_active_fe_index =
          dof_info[0].cell_active_fe_index;
        std::vector<std::vector<unsigned int> > renumbering_fe_index;
      renumbering_fe_index.resize(max_fe_index);
      unsigned int counter,n_macro_cells_before = 0;
      const unsigned int
        start_bound = std::min (size_info.n_active_cells,
                                size_info.boundary_cells_start*n_vectors),
        end_bound   = std::min (size_info.n_active_cells,
                                size_info.boundary_cells_end*n_vectors);
      for(counter=0; counter<start_bound; counter++)
        {
          renumbering_fe_index[cell_active_fe_index[renumbering[counter]]].
            push_back(renumbering[counter]);
        }
      counter = 0;
      for (unsigned int j=0;j<max_fe_index;j++)
        {
          for(unsigned int jj=0;jj<renumbering_fe_index[j].size();jj++)
            renumbering[counter++] = renumbering_fe_index[j][jj];
          irregular_cells[renumbering_fe_index[j].size()/n_vectors+
                          n_macro_cells_before] =
            renumbering_fe_index[j].size()%n_vectors;
          n_macro_cells_before += (renumbering_fe_index[j].size()+n_vectors-1)/
            n_vectors;
          renumbering_fe_index[j].resize(0);
        }
      unsigned int new_boundary_start = n_macro_cells_before;
      for(counter = start_bound; counter < end_bound; counter++)
        {
          renumbering_fe_index[cell_active_fe_index[renumbering[counter]]].
            push_back(renumbering[counter]);
        }
      counter = start_bound;
      for (unsigned int j=0;j<max_fe_index;j++)
        {
          for(unsigned int jj=0;jj<renumbering_fe_index[j].size();jj++)
            renumbering[counter++] = renumbering_fe_index[j][jj];
          irregular_cells[renumbering_fe_index[j].size()/n_vectors+
                          n_macro_cells_before] =
            renumbering_fe_index[j].size()%n_vectors;
          n_macro_cells_before += (renumbering_fe_index[j].size()+n_vectors-1)/
            n_vectors;
          renumbering_fe_index[j].resize(0);
        }
      unsigned int new_boundary_end = n_macro_cells_before;
      for(counter=end_bound; counter<n_active_cells; counter++)
        {
          renumbering_fe_index[cell_active_fe_index[renumbering[counter]]].
            push_back(renumbering[counter]);
        }
      counter = end_bound;
      for (unsigned int j=0;j<max_fe_index;j++)
        {
          for(unsigned int jj=0;jj<renumbering_fe_index[j].size();jj++)
            renumbering[counter++] = renumbering_fe_index[j][jj];
          irregular_cells[renumbering_fe_index[j].size()/n_vectors+
                          n_macro_cells_before] =
            renumbering_fe_index[j].size()%n_vectors;
          n_macro_cells_before += (renumbering_fe_index[j].size()+n_vectors-1)/
            n_vectors;
        }
      AssertIndexRange (n_macro_cells_before,
                        size_info.n_macro_cells + 3*max_fe_index+1);
      irregular_cells.resize (n_macro_cells_before);
      size_info.n_macro_cells = n_macro_cells_before;
      size_info.boundary_cells_start = new_boundary_start;
      size_info.boundary_cells_end = new_boundary_end;
      }
    }
                                // Finally perform the renumbering. We also
                                // want to group several cells together to one
                                // "macro-cell" for vectorization (where the
                                // arithmetic operations will then be done
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
    cell_level_index.reserve(size_info.n_macro_cells*n_vectors);
    unsigned int position_cell=0;
    for (unsigned int i=0; i<size_info.n_macro_cells; ++i)
      {
        unsigned int n_comp = (irregular_cells[i]>0)?
          irregular_cells[i] : n_vectors;
        for (unsigned int j=0; j<n_comp; ++j)
          cell_level_index.push_back
            (cell_level_index_old[renumbering[position_cell+j]]);

                                // generate a cell and level index
                                // also when we have not filled up
                                // n_vectors cells. This is needed for
                                // MappingInfo when the transformation
                                // data is initialized. We just set
                                // the value to the last valid cell in
                                // that case.
        for (unsigned int j=n_comp; j<n_vectors; ++j)
          cell_level_index.push_back
            (cell_level_index_old[renumbering[position_cell+n_comp-1]]);
        position_cell += n_comp;
      }
    AssertDimension (position_cell, size_info.n_active_cells);
    AssertDimension (cell_level_index.size(),size_info.n_macro_cells*n_vectors);
  }

                                // set constraint pool and reorder the indices
  constraint_pool.row_index =
    constraint_values.constraint_pool.row_index;
  constraint_pool.data.resize (constraint_values.constraint_pool.data.size());
  std::copy (constraint_values.constraint_pool.data.begin(),
             constraint_values.constraint_pool.data.end(),
             constraint_pool.data.begin());
  for (unsigned int no=0; no<n_fe; ++no)
    {
      dof_info[no].reorder_cells(size_info, renumbering,
                                 constraint_pool.row_index,
                                 irregular_cells, n_vectors);
    }

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
  indices_are_initialized = false;
  mapping_is_initialized  = false;
}



template <int dim, typename Number>
std::size_t MatrixFree<dim,Number>::memory_consumption () const
{
  std::size_t memory = MemoryConsumption::memory_consumption (dof_info);
  memory += MemoryConsumption::memory_consumption (cell_level_index);
  memory += MemoryConsumption::memory_consumption (fe_evaluation_data);
  memory += MemoryConsumption::memory_consumption (constraint_pool);
  memory += MemoryConsumption::memory_consumption (task_info);
  memory += sizeof(this);
  memory += mapping_info.memory_consumption();
  return memory;
}


template <int dim, typename Number>
template <typename STREAM>
void MatrixFree<dim,Number>::print_memory_consumption (STREAM &out) const
{
  out << "  Memory cell FE operator total: --> ";
  size_info.print_mem (out, memory_consumption());
  out << "   Memory cell index:                ";
  size_info.print_mem (out, MemoryConsumption::memory_consumption (cell_level_index));
  for (unsigned int j=0; j<dof_info.size(); ++ j)
    {
      out << "   Memory DoFInfo component "<< j << std::endl;
      dof_info[j].print_memory_consumption(out, size_info);
    }

  out << "   Memory mapping info" << std::endl;
  mapping_info.print_memory_consumption(out, size_info);

  out << "   Memory unit cell shape data:      ";
  size_info.print_mem (out, MemoryConsumption::memory_consumption (fe_evaluation_data));
  if (task_info.use_multithreading == true)
    {
      out << "   Memory task partitioning info:    ";
      size_info.print_mem (out, MemoryConsumption::memory_consumption (task_info));
    }
}



template <int dim, typename Number>
void MatrixFree<dim,Number>::print (std::ostream &out) const
{
                                // print indices local to global
  for (unsigned int no=0; no<dof_info.size(); ++no)
    {
      out << "\n-- Index data for component " << no << " --" << std::endl;
      dof_info[no].print (constraint_pool, out);
      out << std::endl;
    }
}


DEAL_II_NAMESPACE_CLOSE
