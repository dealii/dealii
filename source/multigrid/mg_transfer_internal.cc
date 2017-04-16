// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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


#include <deal.II/distributed/tria.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/matrix_free/shape_info.h>
#include <deal.II/multigrid/mg_transfer_internal.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MGTransfer
  {

    // Internal data structure that is used in the MPI communication in
    // fill_copy_indices().  It represents an entry in the copy_indices* map,
    // that associates a level dof index with a global dof index.
    struct DoFPair
    {
      unsigned int level;
      types::global_dof_index global_dof_index;
      types::global_dof_index level_dof_index;

      DoFPair(const unsigned int level,
              const types::global_dof_index global_dof_index,
              const types::global_dof_index level_dof_index)
        :
        level(level),
        global_dof_index(global_dof_index),
        level_dof_index(level_dof_index)
      {}

      DoFPair()
        :
        level (numbers::invalid_unsigned_int),
        global_dof_index (numbers::invalid_dof_index),
        level_dof_index (numbers::invalid_dof_index)
      {}
    };

    // Internal function for filling the copy indices from global to level
    // indices
    template <int dim, int spacedim>
    void fill_copy_indices(const dealii::DoFHandler<dim,spacedim>                                                  &mg_dof,
                           const MGConstrainedDoFs                                                                 *mg_constrained_dofs,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices_global_mine,
                           std::vector<std::vector<std::pair<types::global_dof_index, types::global_dof_index> > > &copy_indices_level_mine)
    {
      // Now we are filling the variables copy_indices*, which are essentially
      // maps from global to mgdof for each level stored as a std::vector of
      // pairs. We need to split this map on each level depending on the
      // ownership of the global and mgdof, so that we later not access
      // non-local elements in copy_to/from_mg.
      // We keep track in the bitfield dof_touched which global dof has been
      // processed already (on the current level). This is the same as the
      // multigrid running in serial.

      // map cpu_index -> vector of data
      // that will be copied into copy_indices_level_mine
      std::vector<DoFPair> send_data_temp;

      const unsigned int n_levels = mg_dof.get_triangulation().n_global_levels();
      copy_indices.resize(n_levels);
      copy_indices_global_mine.resize(n_levels);
      copy_indices_level_mine.resize(n_levels);
      IndexSet globally_relevant;
      DoFTools::extract_locally_relevant_dofs(mg_dof, globally_relevant);

      const unsigned int dofs_per_cell = mg_dof.get_fe().dofs_per_cell;
      std::vector<types::global_dof_index> global_dof_indices (dofs_per_cell);
      std::vector<types::global_dof_index> level_dof_indices  (dofs_per_cell);

      for (unsigned int level=0; level<n_levels; ++level)
        {
          std::vector<bool> dof_touched(globally_relevant.n_elements(), false);
          copy_indices[level].clear();
          copy_indices_level_mine[level].clear();
          copy_indices_global_mine[level].clear();

          typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator
          level_cell = mg_dof.begin_active(level);
          const typename dealii::DoFHandler<dim,spacedim>::active_cell_iterator
          level_end  = mg_dof.end_active(level);

          for (; level_cell!=level_end; ++level_cell)
            {
              if (mg_dof.get_triangulation().locally_owned_subdomain()!=numbers::invalid_subdomain_id
                  &&  (level_cell->level_subdomain_id()==numbers::artificial_subdomain_id
                       ||  level_cell->subdomain_id()==numbers::artificial_subdomain_id)
                 )
                continue;

              // get the dof numbers of this cell for the global and the
              // level-wise numbering
              level_cell->get_dof_indices (global_dof_indices);
              level_cell->get_mg_dof_indices (level_dof_indices);

              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  // we need to ignore if the DoF is on a refinement edge (hanging node)
                  if (mg_constrained_dofs != nullptr
                      && mg_constrained_dofs->at_refinement_edge(level, level_dof_indices[i]))
                    continue;
                  types::global_dof_index global_idx = globally_relevant.index_within_set(global_dof_indices[i]);
                  //skip if we did this global dof already (on this or a coarser level)
                  if (dof_touched[global_idx])
                    continue;
                  bool global_mine = mg_dof.locally_owned_dofs().is_element(global_dof_indices[i]);
                  bool level_mine = mg_dof.locally_owned_mg_dofs(level).is_element(level_dof_indices[i]);


                  if (global_mine && level_mine)
                    {
                      copy_indices[level].emplace_back(global_dof_indices[i],
                                                       level_dof_indices[i]);
                    }
                  else if (global_mine)
                    {
                      copy_indices_global_mine[level].emplace_back(global_dof_indices[i],
                                                                   level_dof_indices[i]);

                      //send this to the owner of the level_dof:
                      send_data_temp.emplace_back(level,
                                                  global_dof_indices[i],
                                                  level_dof_indices[i]);
                    }
                  else
                    {
                      // somebody will send those to me
                    }

                  dof_touched[global_idx] = true;
                }
            }
        }

      const dealii::parallel::distributed::Triangulation<dim,spacedim> *tria =
        (dynamic_cast<const parallel::distributed::Triangulation<dim,spacedim>*>
         (&mg_dof.get_triangulation()));
      AssertThrow(send_data_temp.size()==0 || tria!=nullptr, ExcMessage("parallel Multigrid only works with a distributed Triangulation!"));

#ifdef DEAL_II_WITH_MPI
      if (tria)
        {
          // TODO: Searching the owner for every single DoF becomes quite
          // inefficient. Please fix this, Timo.
          // The list of neighbors is symmetric (our neighbors have us as a
          // neighbor), so we can use it to send and to know how many messages
          // we will get.
          std::set<types::subdomain_id> neighbors = tria->level_ghost_owners();
          std::map<int, std::vector<DoFPair> > send_data;

          // * find owners of the level dofs and insert into send_data accordingly
          for (typename std::vector<DoFPair>::iterator dofpair=send_data_temp.begin(); dofpair != send_data_temp.end(); ++dofpair)
            {
              std::set<types::subdomain_id>::iterator it;
              for (it = neighbors.begin(); it != neighbors.end(); ++it)
                {
                  if (mg_dof.locally_owned_mg_dofs_per_processor(dofpair->level)[*it].is_element(dofpair->level_dof_index))
                    {
                      send_data[*it].push_back(*dofpair);
                      break;
                    }
                }
              // Is this level dof not owned by any of our neighbors? That would
              // certainly be a bug!
              Assert(it!=neighbors.end(), ExcMessage("could not find DoF owner."));
            }

          // * send
          std::vector<MPI_Request> requests;
          {
            for (std::set<types::subdomain_id>::iterator it = neighbors.begin(); it != neighbors.end(); ++it)
              {
                requests.push_back(MPI_Request());
                unsigned int dest = *it;
                std::vector<DoFPair> &data = send_data[dest];
                // If there is nothing to send, we still need to send a message,
                // because the receiving end will be waitng. In that case we
                // just send an empty message.
                if (data.size())
                  {
                    const int ierr = MPI_Isend(&data[0], data.size()*sizeof(data[0]),
                                               MPI_BYTE, dest, 71, tria->get_communicator(),
                                               &*requests.rbegin());
                    AssertThrowMPI(ierr);
                  }
                else
                  {
                    const int ierr = MPI_Isend(nullptr, 0, MPI_BYTE, dest, 71,
                                               tria->get_communicator(), &*requests.rbegin());
                    AssertThrowMPI(ierr);
                  }
              }
          }

          // * receive
          {
            // We should get one message from each of our neighbors
            std::vector<DoFPair> receive_buffer;
            for (unsigned int counter=0; counter<neighbors.size(); ++counter)
              {
                MPI_Status status;
                int len;
                int ierr = MPI_Probe(MPI_ANY_SOURCE, 71, tria->get_communicator(), &status);
                AssertThrowMPI(ierr);
                ierr = MPI_Get_count(&status, MPI_BYTE, &len);
                AssertThrowMPI(ierr);

                if (len==0)
                  {
                    ierr = MPI_Recv(nullptr, 0, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                                    tria->get_communicator(), &status);
                    AssertThrowMPI(ierr);
                    continue;
                  }

                int count = len / sizeof(DoFPair);
                Assert(static_cast<int>(count * sizeof(DoFPair)) == len, ExcInternalError());
                receive_buffer.resize(count);

                void *ptr = &receive_buffer[0];
                ierr = MPI_Recv(ptr, len, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG,
                                tria->get_communicator(), &status);
                AssertThrowMPI(ierr);

                for (unsigned int i=0; i<receive_buffer.size(); ++i)
                  {
                    copy_indices_level_mine[receive_buffer[i].level].push_back(
                      std::make_pair (receive_buffer[i].global_dof_index, receive_buffer[i].level_dof_index)
                    );
                  }
              }
          }

          // * wait for all MPI_Isend to complete
          if (requests.size() > 0)
            {
              const int ierr = MPI_Waitall(requests.size(), &requests[0], MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
              requests.clear();
            }
#ifdef DEBUG
          // Make sure in debug mode, that everybody sent/received all packages
          // on this level. If a deadlock occurs here, the list of expected
          // senders is not computed correctly.
          const int ierr = MPI_Barrier(tria->get_communicator());
          AssertThrowMPI(ierr);
#endif
        }
#endif

      // Sort the indices. This will produce more reliable debug output for
      // regression tests and likely won't hurt performance even in release
      // mode.
      std::less<std::pair<types::global_dof_index, types::global_dof_index> > compare;
      for (unsigned int level=0; level<copy_indices.size(); ++level)
        std::sort(copy_indices[level].begin(), copy_indices[level].end(), compare);
      for (unsigned int level=0; level<copy_indices_level_mine.size(); ++level)
        std::sort(copy_indices_level_mine[level].begin(), copy_indices_level_mine[level].end(), compare);
      for (unsigned int level=0; level<copy_indices_global_mine.size(); ++level)
        std::sort(copy_indices_global_mine[level].begin(), copy_indices_global_mine[level].end(), compare);
    }



    // initialize the vectors needed for the transfer (and merge with the
    // content in copy_indices_global_mine)
    template <typename Number>
    void
    reinit_ghosted_vector(const IndexSet &locally_owned,
                          std::vector<types::global_dof_index> &ghosted_level_dofs,
                          const MPI_Comm &communicator,
                          LinearAlgebra::distributed::Vector<Number> &ghosted_level_vector,
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
          const std::shared_ptr<const Utilities::MPI::Partitioner> part
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
    inline void
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

    // given the collection of child cells in lexicographic ordering as seen
    // from the parent, compute the first index of the given child
    template <int dim>
    unsigned int
    compute_shift_within_children(const unsigned int child,
                                  const unsigned int fe_shift_1d,
                                  const unsigned int fe_degree)
    {
      // we put the degrees of freedom of all child cells in lexicographic
      // ordering
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

    template <int dim, typename Number>
    void setup_element_info(ElementInfo<Number> &elem_info,const FiniteElement<1> &fe,
                            const dealii::DoFHandler<dim> &mg_dof)
    {
      // currently, we have only FE_Q and FE_DGQ type elements implemented
      elem_info.n_components = mg_dof.get_fe().element_multiplicity(0);
      AssertDimension(Utilities::fixed_power<dim>(fe.dofs_per_cell)*elem_info.n_components,
                      mg_dof.get_fe().dofs_per_cell);
      AssertDimension(fe.degree, mg_dof.get_fe().degree);
      elem_info.fe_degree = fe.degree;
      elem_info.element_is_continuous = fe.dofs_per_vertex > 0;
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

      // step 1.3: create a dummy 1D quadrature formula to extract the
      // lexicographic numbering for the elements
      std::vector<Point<1> > basic_support_points = fe.get_unit_support_points();
      Assert(fe.dofs_per_vertex == 0 || fe.dofs_per_vertex == 1,
             ExcNotImplemented());
      const unsigned int shift = fe.dofs_per_cell - fe.dofs_per_vertex;
      const unsigned int n_child_dofs_1d = (fe.dofs_per_vertex > 0 ?
                                            (2 * fe.dofs_per_cell - 1) :
                                            (2 * fe.dofs_per_cell));

      elem_info.n_child_cell_dofs = elem_info.n_components*Utilities::fixed_power<dim>(n_child_dofs_1d);
      const Quadrature<1> dummy_quadrature(std::vector<Point<1> >(1, Point<1>()));
      internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info;
      shape_info.reinit(dummy_quadrature, mg_dof.get_fe(), 0);
      elem_info.lexicographic_numbering = shape_info.lexicographic_numbering;

      // step 1.4: get the 1d prolongation matrix and combine from both children
      elem_info.prolongation_matrix_1d.resize(fe.dofs_per_cell*n_child_dofs_1d);

      for (unsigned int c=0; c<GeometryInfo<1>::max_children_per_cell; ++c)
        for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
          for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
            elem_info.prolongation_matrix_1d[i*n_child_dofs_1d+j+c*shift] =
              fe.get_prolongation_matrix(c)(renumbering[j],renumbering[i]);
    }


    // Sets up most of the internal data structures of the MGTransferMatrixFree
    // class
    template <int dim, typename Number>
    void setup_transfer(const dealii::DoFHandler<dim>                                       &mg_dof,
                        const MGConstrainedDoFs                                             *mg_constrained_dofs,
                        ElementInfo<Number>                                                 &elem_info,
                        std::vector<std::vector<unsigned int> >                             &level_dof_indices,
                        std::vector<std::vector<std::pair<unsigned int,unsigned int> > >    &parent_child_connect,
                        std::vector<unsigned int>                                           &n_owned_level_cells,
                        std::vector<std::vector<std::vector<unsigned short> > >             &dirichlet_indices,
                        std::vector<std::vector<Number> >                                   &weights_on_refined,
                        std::vector<std::vector<std::pair<unsigned int, unsigned int> > >   &copy_indices_global_mine,
                        MGLevelObject<LinearAlgebra::distributed::Vector<Number> >          &ghosted_level_vector)
    {
      level_dof_indices.clear();
      parent_child_connect.clear();
      n_owned_level_cells.clear();
      dirichlet_indices.clear();
      weights_on_refined.clear();

      // we collect all child DoFs of a mother cell together. For faster
      // tensorized operations, we align the degrees of freedom
      // lexicographically. We distinguish FE_Q elements and FE_DGQ elements

      const dealii::Triangulation<dim> &tria = mg_dof.get_triangulation();

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
      std::shared_ptr<FiniteElement<1> > fe_1d
      (FETools::get_fe_by_name<1,1>(fe_name));
      const FiniteElement<1> &fe = *fe_1d;

      setup_element_info(elem_info,fe,mg_dof);


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
      if (ghosted_level_vector.max_level() != n_levels-1)
        ghosted_level_vector.resize(0, n_levels-1);

      for (unsigned int level=n_levels-1; level > 0; --level)
        {
          unsigned int counter = 0;
          std::vector<types::global_dof_index> global_level_dof_indices;
          std::vector<types::global_dof_index> global_level_dof_indices_remote;
          std::vector<types::global_dof_index> ghosted_level_dofs;
          std::vector<types::global_dof_index> global_level_dof_indices_l0;
          std::vector<types::global_dof_index> ghosted_level_dofs_l0;

          // step 2.1: loop over the cells on the coarse side
          typename dealii::DoFHandler<dim>::cell_iterator cell, endc = mg_dof.end(level-1);
          for (cell = mg_dof.begin(level-1); cell != endc; ++cell)
            {
              // need to look into a cell if it has children and it is locally
              // owned
              if (!cell->has_children())
                continue;

              bool consider_cell = false;
              if (tria.locally_owned_subdomain()==numbers::invalid_subdomain_id
                  || cell->level_subdomain_id()==tria.locally_owned_subdomain()
                 )
                consider_cell = true;

              // due to the particular way we store DoF indices (via children),
              // we also need to add the DoF indices for coarse cells where we
              // own at least one child
              bool cell_is_remote = !consider_cell;
              for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
                if (cell->child(c)->level_subdomain_id()==tria.locally_owned_subdomain())
                  {
                    consider_cell = true;
                    break;
                  }

              if (!consider_cell)
                continue;

              // step 2.2: loop through children and append the dof indices to
              // the appropriate list. We need separate lists for the owned
              // coarse cell case (which will be part of
              // restriction/prolongation between level-1 and level) and the
              // remote case (which needs to store DoF indices for the
              // operations between level and level+1).
              AssertDimension(cell->n_children(),
                              GeometryInfo<dim>::max_children_per_cell);
              std::vector<types::global_dof_index> &next_indices =
                cell_is_remote ? global_level_dof_indices_remote : global_level_dof_indices;
              const std::size_t start_index = next_indices.size();
              next_indices.resize(start_index + elem_info.n_child_cell_dofs,
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
                                         fe.degree, elem_info.lexicographic_numbering,
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
                      // be placed at the end of the list; however, we do not
                      // yet know the exact position in the array, so shift
                      // their parent index by the number of cells so we can set
                      // the correct number after the end of this loop
                      if (cell_is_remote)
                        parent_index = start_index/elem_info.n_child_cell_dofs + tria.n_cells(level);
                      parent_child_connect[level][child_index] =
                        std::make_pair(parent_index, c);
                      AssertIndexRange(mg_dof.get_fe().dofs_per_cell,
                                       static_cast<unsigned short>(-1));

                      // set Dirichlet boundary conditions (as a list of
                      // constrained DoFs) for the child
                      if (mg_constrained_dofs != nullptr)
                        for (unsigned int i=0; i<mg_dof.get_fe().dofs_per_cell; ++i)
                          if (mg_constrained_dofs->is_boundary_index(level,
                                                                     local_dof_indices[elem_info.lexicographic_numbering[i]]))
                            dirichlet_indices[level][child_index].push_back(i);
                    }
                }
              if (!cell_is_remote)
                {
                  AssertIndexRange(static_cast<unsigned int>(cell->index()),
                                   coarse_level_indices[level-1].size());
                  coarse_level_indices[level-1][cell->index()] = counter++;
                }

              // step 2.4: include indices for the coarsest cells. we still
              // insert the indices as if they were from a child in order to use
              // the same code (the coarsest level does not matter much in terms
              // of memory, so we gain in code simplicity)
              if (level == 1 && !cell_is_remote)
                {
                  cell->get_mg_dof_indices(local_dof_indices);

                  const IndexSet &owned_level_dofs_l0 = mg_dof.locally_owned_mg_dofs(0);
                  for (unsigned int i=0; i<local_dof_indices.size(); ++i)
                    if (!owned_level_dofs_l0.is_element(local_dof_indices[i]))
                      ghosted_level_dofs_l0.push_back(local_dof_indices[i]);

                  const std::size_t start_index = global_level_dof_indices_l0.size();
                  global_level_dof_indices_l0.resize(start_index+elem_info.n_child_cell_dofs,
                                                     numbers::invalid_dof_index);
                  add_child_indices<dim>(0, fe.dofs_per_cell - fe.dofs_per_vertex,
                                         fe.degree, elem_info.lexicographic_numbering,
                                         local_dof_indices,
                                         &global_level_dof_indices_l0[start_index]);

                  dirichlet_indices[0].emplace_back();
                  if (mg_constrained_dofs != nullptr)
                    for (unsigned int i=0; i<mg_dof.get_fe().dofs_per_cell; ++i)
                      if (mg_constrained_dofs->is_boundary_index(0, local_dof_indices[elem_info.lexicographic_numbering[i]]))
                        dirichlet_indices[0].back().push_back(i);
                }
            }

          // step 2.5: store information about the current level and prepare the
          // Dirichlet indices and parent-child relationship for the next
          // coarser level
          AssertDimension(counter*elem_info.n_child_cell_dofs, global_level_dof_indices.size());
          n_owned_level_cells[level-1] = counter;
          dirichlet_indices[level-1].resize(counter);
          parent_child_connect[level-1].
          resize(counter, std::make_pair(numbers::invalid_unsigned_int,
                                         numbers::invalid_unsigned_int));

          // step 2.6: put the cells with remotely owned parent to the end of
          // the list (these are needed for the transfer from level to level+1
          // but not for the transfer from level-1 to level).
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
            ptria != nullptr ? ptria->get_communicator() : MPI_COMM_SELF;

          reinit_ghosted_vector (mg_dof.locally_owned_mg_dofs(level),
                                 ghosted_level_dofs, communicator,
                                 ghosted_level_vector[level],
                                 copy_indices_global_mine[level]);

          copy_indices_to_mpi_local_numbers(*ghosted_level_vector[level].get_partitioner(),
                                            global_level_dof_indices,
                                            global_level_dof_indices_remote,
                                            level_dof_indices[level]);
          // step 2.8: Initialize the ghosted vector for level 0
          if (level == 1)
            {
              for (unsigned int i = 0; i<parent_child_connect[0].size(); ++i)
                parent_child_connect[0][i] = std::make_pair(i, 0U);

              reinit_ghosted_vector (mg_dof.locally_owned_mg_dofs(0),
                                     ghosted_level_dofs_l0, communicator,
                                     ghosted_level_vector[0],
                                     copy_indices_global_mine[0]);

              copy_indices_to_mpi_local_numbers(*ghosted_level_vector[0].get_partitioner(),
                                                global_level_dof_indices_l0,
                                                std::vector<types::global_dof_index>(),
                                                level_dof_indices[0]);
            }
        }

      // ---------------------- 3. compute weights to make restriction additive

      const unsigned int n_child_dofs_1d = fe.degree + 1 + fe.dofs_per_cell - fe.dofs_per_vertex;

      // get the valence of the individual components and compute the weights as
      // the inverse of the valence
      weights_on_refined.resize(n_levels-1);
      for (unsigned int level = 1; level<n_levels; ++level)
        {
          ghosted_level_vector[level] = 0;
          for (unsigned int c=0; c<n_owned_level_cells[level-1]; ++c)
            for (unsigned int j=0; j<elem_info.n_child_cell_dofs; ++j)
              ghosted_level_vector[level].local_element(level_dof_indices[level][elem_info.n_child_cell_dofs*c+j]) += Number(1.);
          ghosted_level_vector[level].compress(VectorOperation::add);
          ghosted_level_vector[level].update_ghost_values();

          std::vector<unsigned int> degree_to_3 (n_child_dofs_1d);
          degree_to_3[0] = 0;
          for (unsigned int i=1; i<n_child_dofs_1d-1; ++i)
            degree_to_3[i] = 1;
          degree_to_3.back() = 2;

          // we only store 3^dim weights because all dofs on a line have the
          // same valence, and all dofs on a quad have the same valence.
          weights_on_refined[level-1].resize(n_owned_level_cells[level-1]*Utilities::fixed_power<dim>(3));
          for (unsigned int c=0; c<n_owned_level_cells[level-1]; ++c)
            for (unsigned int k=0, m=0; k<(dim>2 ? n_child_dofs_1d : 1); ++k)
              for (unsigned int j=0; j<(dim>1 ? n_child_dofs_1d : 1); ++j)
                {
                  unsigned int shift = 9*degree_to_3[k] + 3*degree_to_3[j];
                  for (unsigned int i=0; i<n_child_dofs_1d; ++i, ++m)
                    weights_on_refined[level-1][c*Utilities::fixed_power<dim>(3)+shift+degree_to_3[i]] = Number(1.)/
                        ghosted_level_vector[level].local_element(level_dof_indices[level][elem_info.n_child_cell_dofs*c+m]);
                }
        }

    }

  }
}

// Explicit instantiations

#include "mg_transfer_internal.inst"

DEAL_II_NAMESPACE_CLOSE
