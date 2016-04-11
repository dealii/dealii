// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2015 by the deal.II authors
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

#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/types.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/fe/fe.h>
#include <deal.II/hp/dof_handler.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>

#include <deal.II/multigrid/mg_tools.h>

#include <deal.II/distributed/tria.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/king_ordering.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <boost/random.hpp>
#include <boost/random/uniform_int_distribution.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <functional>


DEAL_II_NAMESPACE_OPEN

namespace DoFRenumbering
{
  namespace boost
  {
    namespace boosttypes
    {
      using namespace ::boost;
      using namespace std;

      typedef adjacency_list<vecS, vecS, undirectedS,
              property<vertex_color_t, default_color_type,
              property<vertex_degree_t,int> > > Graph;
      typedef graph_traits<Graph>::vertex_descriptor Vertex;
      typedef graph_traits<Graph>::vertices_size_type size_type;

      typedef std::pair<size_type, size_type> Pair;
    }


    namespace internal
    {
      template <typename DoFHandlerType>
      void create_graph
      (const DoFHandlerType                                                          &dof_handler,
       const bool                                                                     use_constraints,
       boosttypes::Graph                                                             &graph,
       boosttypes::property_map<boosttypes::Graph,boosttypes::vertex_degree_t>::type &graph_degree)
      {
        {
          // create intermediate sparsity pattern
          // (faster than directly submitting
          // indices)
          ConstraintMatrix constraints;
          if (use_constraints)
            DoFTools::make_hanging_node_constraints (dof_handler, constraints);
          constraints.close ();
          DynamicSparsityPattern dsp (dof_handler.n_dofs(),
                                      dof_handler.n_dofs());
          DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints);

          // submit the entries to the boost graph
          for (unsigned int row=0; row<dsp.n_rows(); ++row)
            for (unsigned int col=0; col < dsp.row_length(row); ++col)
              add_edge (row, dsp.column_number (row, col), graph);
        }

        boosttypes::graph_traits<boosttypes::Graph>::vertex_iterator ui, ui_end;

        graph_degree = get(::boost::vertex_degree, graph);
        for (::boost::tie(ui, ui_end) = vertices(graph); ui != ui_end; ++ui)
          graph_degree[*ui] = degree(*ui, graph);
      }
    }


    template <typename DoFHandlerType>
    void
    Cuthill_McKee (DoFHandlerType  &dof_handler,
                   const bool       reversed_numbering,
                   const bool       use_constraints)
    {
      std::vector<types::global_dof_index> renumbering(dof_handler.n_dofs(),
                                                       DoFHandlerType::invalid_dof_index);
      compute_Cuthill_McKee(renumbering, dof_handler, reversed_numbering,
                            use_constraints);

      // actually perform renumbering;
      // this is dimension specific and
      // thus needs an own function
      dof_handler.renumber_dofs (renumbering);
    }


    template <typename DoFHandlerType>
    void
    compute_Cuthill_McKee (std::vector<types::global_dof_index> &new_dof_indices,
                           const DoFHandlerType                 &dof_handler,
                           const bool                            reversed_numbering,
                           const bool                            use_constraints)
    {
      boosttypes::Graph
      graph(dof_handler.n_dofs());
      boosttypes::property_map<boosttypes::Graph,boosttypes::vertex_degree_t>::type
      graph_degree;

      internal::create_graph (dof_handler, use_constraints, graph, graph_degree);

      boosttypes::property_map<boosttypes::Graph, boosttypes::vertex_index_t>::type
      index_map = get(::boost::vertex_index, graph);


      std::vector<boosttypes::Vertex> inv_perm(num_vertices(graph));

      if (reversed_numbering == false)
        ::boost::cuthill_mckee_ordering(graph, inv_perm.rbegin(),
                                        get(::boost::vertex_color, graph),
                                        make_degree_map(graph));
      else
        ::boost::cuthill_mckee_ordering(graph, inv_perm.begin(),
                                        get(::boost::vertex_color, graph),
                                        make_degree_map(graph));

      for (boosttypes::size_type c = 0; c != inv_perm.size(); ++c)
        new_dof_indices[index_map[inv_perm[c]]] = c;

      Assert (std::find (new_dof_indices.begin(), new_dof_indices.end(),
                         DoFHandlerType::invalid_dof_index) == new_dof_indices.end(),
              ExcInternalError());
    }



    template <typename DoFHandlerType>
    void
    king_ordering (DoFHandlerType  &dof_handler,
                   const bool       reversed_numbering,
                   const bool       use_constraints)
    {
      std::vector<types::global_dof_index> renumbering(dof_handler.n_dofs(),
                                                       DoFHandlerType::invalid_dof_index);
      compute_king_ordering(renumbering, dof_handler, reversed_numbering,
                            use_constraints);

      // actually perform renumbering;
      // this is dimension specific and
      // thus needs an own function
      dof_handler.renumber_dofs (renumbering);
    }


    template <typename DoFHandlerType>
    void
    compute_king_ordering (std::vector<types::global_dof_index> &new_dof_indices,
                           const DoFHandlerType                 &dof_handler,
                           const bool                            reversed_numbering,
                           const bool                            use_constraints)
    {
      boosttypes::Graph
      graph(dof_handler.n_dofs());
      boosttypes::property_map<boosttypes::Graph,boosttypes::vertex_degree_t>::type
      graph_degree;

      internal::create_graph (dof_handler, use_constraints, graph, graph_degree);

      boosttypes::property_map<boosttypes::Graph, boosttypes::vertex_index_t>::type
      index_map = get(::boost::vertex_index, graph);


      std::vector<boosttypes::Vertex> inv_perm(num_vertices(graph));

      if (reversed_numbering == false)
        ::boost::king_ordering(graph, inv_perm.rbegin());
      else
        ::boost::king_ordering(graph, inv_perm.begin());

      for (boosttypes::size_type c = 0; c != inv_perm.size(); ++c)
        new_dof_indices[index_map[inv_perm[c]]] = c;

      Assert (std::find (new_dof_indices.begin(), new_dof_indices.end(),
                         DoFHandlerType::invalid_dof_index) == new_dof_indices.end(),
              ExcInternalError());
    }



    template <typename DoFHandlerType>
    void
    minimum_degree (DoFHandlerType  &dof_handler,
                    const bool       reversed_numbering,
                    const bool       use_constraints)
    {
      std::vector<types::global_dof_index> renumbering(dof_handler.n_dofs(),
                                                       DoFHandlerType::invalid_dof_index);
      compute_minimum_degree(renumbering, dof_handler, reversed_numbering,
                             use_constraints);

      // actually perform renumbering;
      // this is dimension specific and
      // thus needs an own function
      dof_handler.renumber_dofs (renumbering);
    }


    template <typename DoFHandlerType>
    void
    compute_minimum_degree (std::vector<types::global_dof_index> &new_dof_indices,
                            const DoFHandlerType                 &dof_handler,
                            const bool                            reversed_numbering,
                            const bool                            use_constraints)
    {
      (void)use_constraints;
      Assert (use_constraints == false, ExcNotImplemented());

      // the following code is pretty
      // much a verbatim copy of the
      // sample code for the
      // minimum_degree_ordering manual
      // page from the BOOST Graph
      // Library
      using namespace ::boost;

      int delta = 0;

      typedef double Type;

      // must be BGL directed graph now
      typedef adjacency_list<vecS, vecS, directedS>  Graph;
      typedef graph_traits<Graph>::vertex_descriptor Vertex;

      int n = dof_handler.n_dofs();

      Graph G(n);

      std::vector<dealii::types::global_dof_index> dofs_on_this_cell;

      typename DoFHandlerType::active_cell_iterator cell = dof_handler.begin_active(),
                                                    endc = dof_handler.end();

      for (; cell!=endc; ++cell)
        {

          const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;

          dofs_on_this_cell.resize (dofs_per_cell);

          cell->get_active_or_mg_dof_indices (dofs_on_this_cell);
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              if (dofs_on_this_cell[i] > dofs_on_this_cell[j])
                {
                  add_edge (dofs_on_this_cell[i], dofs_on_this_cell[j], G);
                  add_edge (dofs_on_this_cell[j], dofs_on_this_cell[i], G);
                }
        }


      typedef std::vector<int> Vector;


      Vector inverse_perm(n, 0);

      Vector perm(n, 0);


      Vector supernode_sizes(n, 1);
      // init has to be 1

      ::boost::property_map<Graph, vertex_index_t>::type
      id = get(vertex_index, G);


      Vector degree(n, 0);


      minimum_degree_ordering
      (G,
       make_iterator_property_map(&degree[0], id, degree[0]),
       &inverse_perm[0],
       &perm[0],
       make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]),
       delta, id);


      for (int i=0; i<n; ++i)
        {
          Assert (std::find (perm.begin(), perm.end(), i)
                  != perm.end(),
                  ExcInternalError());
          Assert (std::find (inverse_perm.begin(), inverse_perm.end(), i)
                  != inverse_perm.end(),
                  ExcInternalError());
          Assert (inverse_perm[perm[i]] == i, ExcInternalError());
        }

      if (reversed_numbering == true)
        std::copy (perm.begin(), perm.end(),
                   new_dof_indices.begin());
      else
        std::copy (inverse_perm.begin(), inverse_perm.end(),
                   new_dof_indices.begin());
    }

  }  // namespace boost



  template <typename DoFHandlerType>
  void
  Cuthill_McKee (DoFHandlerType                             &dof_handler,
                 const bool                                  reversed_numbering,
                 const bool                                  use_constraints,
                 const std::vector<types::global_dof_index> &starting_indices)
  {
    std::vector<types::global_dof_index> renumbering(dof_handler.locally_owned_dofs().n_elements(),
                                                     DoFHandlerType::invalid_dof_index);
    compute_Cuthill_McKee(renumbering, dof_handler, reversed_numbering,
                          use_constraints, starting_indices);

    // actually perform renumbering;
    // this is dimension specific and
    // thus needs an own function
    dof_handler.renumber_dofs (renumbering);
  }



  template <typename DoFHandlerType>
  void
  compute_Cuthill_McKee (std::vector<types::global_dof_index>       &new_indices,
                         const DoFHandlerType                       &dof_handler,
                         const bool                                  reversed_numbering,
                         const bool                                  use_constraints,
                         const std::vector<types::global_dof_index> &starting_indices)
  {
    // see if there is anything to do at all or whether we can skip the work on this processor
    if (dof_handler.locally_owned_dofs().n_elements() == 0)
      {
        Assert (new_indices.size() == 0, ExcInternalError());
        return;
      }

    // make the connection graph. in 2d/3d use an intermediate compressed
    // sparsity pattern since the we don't have very good estimates for
    // max_couplings_between_dofs() in 3d and this then leads to excessive
    // memory consumption
    //
    // note that if constraints are not requested, then the 'constraints'
    // object will be empty and nothing happens
    ConstraintMatrix constraints;
    if (use_constraints)
      {
        IndexSet relevant_dofs;
        DoFTools::extract_locally_relevant_dofs(dof_handler, relevant_dofs);
        constraints.reinit(relevant_dofs);
        DoFTools::make_hanging_node_constraints (dof_handler, constraints);
      }
    constraints.close ();

    const IndexSet locally_owned = dof_handler.locally_owned_dofs();

    // otherwise compute the Cuthill-McKee permutation
    DynamicSparsityPattern dsp (dof_handler.n_dofs(),
                                dof_handler.n_dofs(),
                                locally_owned);
    DoFTools::make_sparsity_pattern (dof_handler, dsp, constraints);

    // constraints are not needed anymore
    constraints.clear ();

    // If the index set is not complete, need to get indices in local index
    // space.
    if (locally_owned.n_elements() != locally_owned.size())
      {
        // Create sparsity pattern from dsp by transferring its indices to
        // processor-local index space and doing Cuthill-McKee there
        DynamicSparsityPattern sparsity(locally_owned.n_elements(),
                                        locally_owned.n_elements());
        std::vector<types::global_dof_index> row_entries;
        for (unsigned int i=0; i<locally_owned.n_elements(); ++i)
          {
            const types::global_dof_index row = locally_owned.nth_index_in_set(i);
            const unsigned int row_length = dsp.row_length(row);
            row_entries.clear();
            for (unsigned int j=0; j<row_length; ++j)
              {
                const unsigned int col = dsp.column_number(row, j);
                if (col != row && locally_owned.is_element(col))
                  row_entries.push_back(locally_owned.index_within_set(col));
              }
            sparsity.add_entries(i, row_entries.begin(), row_entries.end(),
                                 true);
          }

        // translate starting indices from global to local indices
        std::vector<types::global_dof_index> local_starting_indices (starting_indices.size());
        for (unsigned int i=0; i<starting_indices.size(); ++i)
          {
            Assert (locally_owned.is_element (starting_indices[i]),
                    ExcMessage ("You specified global degree of freedom "
                                + Utilities::to_string(starting_indices[i]) +
                                " as a starting index, but this index is not among the "
                                "locally owned ones on this processor."));
            local_starting_indices[i] = locally_owned.index_within_set(starting_indices[i]);
          }

        // then do the renumbering on the locally owned portion
        AssertDimension(new_indices.size(), locally_owned.n_elements());
        SparsityTools::reorder_Cuthill_McKee (sparsity, new_indices,
                                              local_starting_indices);
        if (reversed_numbering)
          new_indices = Utilities::reverse_permutation (new_indices);

        // convert indices back to global index space
        for (std::size_t i=0; i<new_indices.size(); ++i)
          new_indices[i] = locally_owned.nth_index_in_set(new_indices[i]);
      }
    else
      {
        AssertDimension(new_indices.size(), dsp.n_rows());
        SparsityTools::reorder_Cuthill_McKee (dsp, new_indices,
                                              starting_indices);
        if (reversed_numbering)
          new_indices = Utilities::reverse_permutation (new_indices);
      }
  }



  template <typename DoFHandlerType>
  void Cuthill_McKee (DoFHandlerType                             &dof_handler,
                      const unsigned int                          level,
                      const bool                                  reversed_numbering,
                      const std::vector<types::global_dof_index> &starting_indices)
  {
    Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
           ExcNotInitialized());

    // make the connection graph
    DynamicSparsityPattern dsp (dof_handler.n_dofs(level),
                                dof_handler.n_dofs(level));
    MGTools::make_sparsity_pattern (dof_handler, dsp, level);

    std::vector<types::global_dof_index> new_indices(dsp.n_rows());
    SparsityTools::reorder_Cuthill_McKee (dsp, new_indices,
                                          starting_indices);

    if (reversed_numbering)
      new_indices = Utilities::reverse_permutation (new_indices);

    // actually perform renumbering;
    // this is dimension specific and
    // thus needs an own function
    dof_handler.renumber_dofs (level, new_indices);
  }



  template <int dim, int spacedim>
  void
  component_wise (DoFHandler<dim,spacedim>        &dof_handler,
                  const std::vector<unsigned int> &component_order_arg)
  {
    std::vector<types::global_dof_index> renumbering (dof_handler.n_locally_owned_dofs(),
                                                      DoFHandler<dim>::invalid_dof_index);

    typename DoFHandler<dim,spacedim>::active_cell_iterator
    start = dof_handler.begin_active();
    const typename DoFHandler<dim,spacedim>::level_cell_iterator
    end = dof_handler.end();

    const types::global_dof_index result =
      compute_component_wise<dim,spacedim,
      typename DoFHandler<dim,spacedim>::active_cell_iterator,
      typename DoFHandler<dim,spacedim>::level_cell_iterator>
      (renumbering, start, end, component_order_arg, false);
    if (result == 0)
      return;

    // verify that the last numbered
    // degree of freedom is either
    // equal to the number of degrees
    // of freedom in total (the
    // sequential case) or in the
    // distributed case at least
    // makes sense
    Assert ((result == dof_handler.n_locally_owned_dofs())
            ||
            ((dof_handler.n_locally_owned_dofs() < dof_handler.n_dofs())
             &&
             (result <= dof_handler.n_dofs())),
            ExcRenumberingIncomplete());

    dof_handler.renumber_dofs (renumbering);

    // for (unsigned int level=0;level<dof_handler.get_triangulation().n_levels();++level)
    //   if (dof_handler.n_dofs(level) != numbers::invalid_dof_index)
    //  component_wise(dof_handler, level, component_order_arg);
  }



  template <int dim>
  void
  component_wise (hp::DoFHandler<dim>             &dof_handler,
                  const std::vector<unsigned int> &component_order_arg)
  {
//TODO: Merge with previous function
    std::vector<types::global_dof_index> renumbering (dof_handler.n_dofs(),
                                                      hp::DoFHandler<dim>::invalid_dof_index);

    typename hp::DoFHandler<dim>::active_cell_iterator
    start = dof_handler.begin_active();
    const typename hp::DoFHandler<dim>::level_cell_iterator
    end = dof_handler.end();

    const types::global_dof_index result =
      compute_component_wise<hp::DoFHandler<dim>::dimension,hp::DoFHandler<dim>::space_dimension,
      typename hp::DoFHandler<dim>::active_cell_iterator,
      typename hp::DoFHandler<dim>::level_cell_iterator>
      (renumbering, start, end, component_order_arg, false);

    if (result == 0) return;

    Assert (result == dof_handler.n_dofs(),
            ExcRenumberingIncomplete());

    dof_handler.renumber_dofs (renumbering);
  }



  template <typename DoFHandlerType>
  void
  component_wise (DoFHandlerType                  &dof_handler,
                  const unsigned int               level,
                  const std::vector<unsigned int> &component_order_arg)
  {
    Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
           ExcNotInitialized());

    std::vector<types::global_dof_index> renumbering (dof_handler.n_dofs(level),
                                                      DoFHandlerType::invalid_dof_index);

    typename DoFHandlerType::level_cell_iterator start =dof_handler.begin(level);
    typename DoFHandlerType::level_cell_iterator end = dof_handler.end(level);

    const types::global_dof_index result =
      compute_component_wise<DoFHandlerType::dimension, DoFHandlerType::space_dimension,
      typename DoFHandlerType::level_cell_iterator, typename DoFHandlerType::level_cell_iterator>
      (renumbering, start, end, component_order_arg, true);

    if (result == 0) return;

    Assert (result == dof_handler.n_dofs(level),
            ExcRenumberingIncomplete());

    if (renumbering.size()!=0)
      dof_handler.renumber_dofs (level, renumbering);
  }



  template <int dim, int spacedim, class ITERATOR, class ENDITERATOR>
  types::global_dof_index
  compute_component_wise (std::vector<types::global_dof_index> &new_indices,
                          const ITERATOR    &start,
                          const ENDITERATOR &end,
                          const std::vector<unsigned int> &component_order_arg,
                          bool is_level_operation)
  {
    const hp::FECollection<dim,spacedim>
    fe_collection (start->get_dof_handler().get_fe ());

    // do nothing if the FE has only
    // one component
    if (fe_collection.n_components() == 1)
      {
        new_indices.resize(0);
        return 0;
      }

    // Copy last argument into a
    // writable vector.
    std::vector<unsigned int> component_order (component_order_arg);
    // If the last argument was an
    // empty vector, set up things to
    // store components in the order
    // found in the system.
    if (component_order.size() == 0)
      for (unsigned int i=0; i<fe_collection.n_components(); ++i)
        component_order.push_back (i);

    Assert (component_order.size() == fe_collection.n_components(),
            ExcDimensionMismatch(component_order.size(), fe_collection.n_components()));

    for (unsigned int i=0; i<component_order.size(); ++i)
      Assert(component_order[i] < fe_collection.n_components(),
             ExcIndexRange(component_order[i], 0, fe_collection.n_components()));

    // vector to hold the dof indices on
    // the cell we visit at a time
    std::vector<types::global_dof_index> local_dof_indices;

    // prebuilt list to which component
    // a given dof on a cell
    // should go. note that we get into
    // trouble here if the shape
    // function is not primitive, since
    // then there is no single vector
    // component to which it
    // belongs. in this case, assign it
    // to the first vector component to
    // which it belongs
    std::vector<std::vector<unsigned int> > component_list (fe_collection.size());
    for (unsigned int f=0; f<fe_collection.size(); ++f)
      {
        const FiniteElement<dim,spacedim> &fe = fe_collection[f];
        const unsigned int dofs_per_cell = fe.dofs_per_cell;
        component_list[f].resize(dofs_per_cell);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          if (fe.is_primitive(i))
            component_list[f][i]
              = component_order[fe.system_to_component_index(i).first];
          else
            {
              const unsigned int comp
                = fe.get_nonzero_components(i).first_selected_component();

              // then associate this degree
              // of freedom with this
              // component
              component_list[f][i] = component_order[comp];
            }
      }

    // set up a map where for each
    // component the respective degrees
    // of freedom are collected.
    //
    // note that this map is sorted by
    // component but that within each
    // component it is NOT sorted by
    // dof index. note also that some
    // dof indices are entered
    // multiply, so we will have to
    // take care of that
    std::vector<std::vector<types::global_dof_index> >
    component_to_dof_map (fe_collection.n_components());
    for (ITERATOR cell=start; cell!=end; ++cell)
      {
        if (is_level_operation)
          {
            //we are dealing with mg dofs, skip foreign level cells:
            if (!cell->is_locally_owned_on_level())
              continue;
          }
        else
          {
            //we are dealing with active dofs, skip the loop if not locally
            // owned:
            if (!cell->is_locally_owned())
              continue;
          }
        // on each cell: get dof indices
        // and insert them into the global
        // list using their component
        const unsigned int fe_index = cell->active_fe_index();
        const unsigned int dofs_per_cell = fe_collection[fe_index].dofs_per_cell;
        local_dof_indices.resize (dofs_per_cell);
        cell->get_active_or_mg_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          if (start->get_dof_handler().locally_owned_dofs().is_element(local_dof_indices[i]))
            component_to_dof_map[component_list[fe_index][i]].
            push_back (local_dof_indices[i]);
      }

    // now we've got all indices sorted
    // into buckets labeled by their
    // target component number. we've
    // only got to traverse this list
    // and assign the new indices
    //
    // however, we first want to sort
    // the indices entered into the
    // buckets to preserve the order
    // within each component and during
    // this also remove duplicate
    // entries
    //
    // note that we no longer have to
    // care about non-primitive shape
    // functions since the buckets
    // corresponding to the second and
    // following vector components of a
    // non-primitive FE will simply be
    // empty, everything being shoved
    // into the first one. The same
    // holds if several components were
    // joined into a single target.
    for (unsigned int component=0; component<fe_collection.n_components();
         ++component)
      {
        std::sort (component_to_dof_map[component].begin(),
                   component_to_dof_map[component].end());
        component_to_dof_map[component]
        .erase (std::unique (component_to_dof_map[component].begin(),
                             component_to_dof_map[component].end()),
                component_to_dof_map[component].end());
      }

    // calculate the number of locally owned
    // DoFs per bucket
    const unsigned int n_buckets = fe_collection.n_components();
    std::vector<types::global_dof_index> shifts(n_buckets);

    if (const parallel::Triangulation<dim,spacedim> *tria
        = (dynamic_cast<const parallel::Triangulation<dim,spacedim>*>
           (&start->get_dof_handler().get_triangulation())))
      {
#ifdef DEAL_II_WITH_MPI
        std::vector<types::global_dof_index> local_dof_count(n_buckets);

        for (unsigned int c=0; c<n_buckets; ++c)
          local_dof_count[c] = component_to_dof_map[c].size();


        // gather information from all CPUs
        std::vector<types::global_dof_index>
        all_dof_counts(fe_collection.n_components() *
                       Utilities::MPI::n_mpi_processes (tria->get_communicator()));

        MPI_Allgather ( &local_dof_count[0],
                        n_buckets, DEAL_II_DOF_INDEX_MPI_TYPE,
                        &all_dof_counts[0],
                        n_buckets, DEAL_II_DOF_INDEX_MPI_TYPE,
                        tria->get_communicator());

        for (unsigned int i=0; i<n_buckets; ++i)
          Assert (all_dof_counts[n_buckets*tria->locally_owned_subdomain()+i]
                  ==
                  local_dof_count[i],
                  ExcInternalError());

        //calculate shifts
        unsigned int cumulated = 0;
        for (unsigned int c=0; c<n_buckets; ++c)
          {
            shifts[c]=cumulated;
            for (types::subdomain_id i=0; i<tria->locally_owned_subdomain(); ++i)
              shifts[c] += all_dof_counts[c+n_buckets*i];
            for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes (tria->get_communicator()); ++i)
              cumulated += all_dof_counts[c+n_buckets*i];
          }
#else
        (void)tria;
        Assert (false, ExcInternalError());
#endif
      }
    else
      {
        shifts[0] = 0;
        for (unsigned int c=1; c<fe_collection.n_components(); ++c)
          shifts[c] = shifts[c-1] + component_to_dof_map[c-1].size();
      }




    // now concatenate all the
    // components in the order the user
    // desired to see
    types::global_dof_index next_free_index = 0;
    for (unsigned int component=0; component<fe_collection.n_components(); ++component)
      {
        const typename std::vector<types::global_dof_index>::const_iterator
        begin_of_component = component_to_dof_map[component].begin(),
        end_of_component   = component_to_dof_map[component].end();

        next_free_index = shifts[component];

        for (typename std::vector<types::global_dof_index>::const_iterator
             dof_index = begin_of_component;
             dof_index != end_of_component; ++dof_index)
          {
            Assert (start->get_dof_handler().locally_owned_dofs()
                    .index_within_set(*dof_index)
                    <
                    new_indices.size(),
                    ExcInternalError());
            new_indices[start->get_dof_handler().locally_owned_dofs()
                        .index_within_set(*dof_index)]
              = next_free_index++;
          }
      }

    return next_free_index;
  }



  template <int dim, int spacedim>
  void
  block_wise (DoFHandler<dim,spacedim> &dof_handler)
  {
    std::vector<types::global_dof_index> renumbering (dof_handler.n_locally_owned_dofs(),
                                                      DoFHandler<dim>::invalid_dof_index);

    typename DoFHandler<dim,spacedim>::active_cell_iterator
    start = dof_handler.begin_active();
    const typename DoFHandler<dim,spacedim>::level_cell_iterator
    end = dof_handler.end();

    const types::global_dof_index result =
      compute_block_wise<dim, spacedim, typename DoFHandler<dim,spacedim>::active_cell_iterator,
      typename DoFHandler<dim,spacedim>::level_cell_iterator>
      (renumbering, start, end, false);
    if (result == 0)
      return;

    // verify that the last numbered
    // degree of freedom is either
    // equal to the number of degrees
    // of freedom in total (the
    // sequential case) or in the
    // distributed case at least
    // makes sense
    Assert ((result == dof_handler.n_locally_owned_dofs())
            ||
            ((dof_handler.n_locally_owned_dofs() < dof_handler.n_dofs())
             &&
             (result <= dof_handler.n_dofs())),
            ExcRenumberingIncomplete());

    dof_handler.renumber_dofs (renumbering);
  }



  template <int dim, int spacedim>
  void
  block_wise (hp::DoFHandler<dim,spacedim> &dof_handler)
  {
    std::vector<types::global_dof_index> renumbering (dof_handler.n_dofs(),
                                                      hp::DoFHandler<dim,spacedim>::invalid_dof_index);

    typename hp::DoFHandler<dim,spacedim>::active_cell_iterator
    start = dof_handler.begin_active();
    const typename hp::DoFHandler<dim,spacedim>::level_cell_iterator
    end = dof_handler.end();

    const types::global_dof_index result =
      compute_block_wise<dim, spacedim, typename hp::DoFHandler<dim,spacedim>::active_cell_iterator,
      typename hp::DoFHandler<dim,spacedim>::level_cell_iterator>(renumbering,
          start, end, false);

    if (result == 0)
      return;

    Assert (result == dof_handler.n_dofs(),
            ExcRenumberingIncomplete());

    dof_handler.renumber_dofs (renumbering);
  }



  template <int dim, int spacedim>
  void
  block_wise (DoFHandler<dim,spacedim> &dof_handler, const unsigned int level)
  {
    Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
           ExcNotInitialized());

    std::vector<types::global_dof_index> renumbering (dof_handler.n_dofs(level),
                                                      DoFHandler<dim, spacedim>::invalid_dof_index);

    typename DoFHandler<dim, spacedim>::level_cell_iterator
    start =dof_handler.begin(level);
    typename DoFHandler<dim, spacedim>::level_cell_iterator
    end = dof_handler.end(level);

    const types::global_dof_index result =
      compute_block_wise<dim, spacedim, typename DoFHandler<dim, spacedim>::level_cell_iterator,
      typename DoFHandler<dim, spacedim>::level_cell_iterator>(
        renumbering, start, end, true);

    if (result == 0) return;

    Assert (result == dof_handler.n_dofs(level),
            ExcRenumberingIncomplete());

    if (renumbering.size()!=0)
      dof_handler.renumber_dofs (level, renumbering);
  }



  template <int dim, int spacedim, class ITERATOR, class ENDITERATOR>
  types::global_dof_index
  compute_block_wise (std::vector<types::global_dof_index> &new_indices,
                      const ITERATOR    &start,
                      const ENDITERATOR &end,
                      const bool is_level_operation)
  {
    const hp::FECollection<dim,spacedim>
    fe_collection (start->get_dof_handler().get_fe ());

    // do nothing if the FE has only
    // one component
    if (fe_collection.n_blocks() == 1)
      {
        new_indices.resize(0);
        return 0;
      }

    // vector to hold the dof indices on
    // the cell we visit at a time
    std::vector<types::global_dof_index> local_dof_indices;

    // prebuilt list to which block
    // a given dof on a cell
    // should go.
    std::vector<std::vector<types::global_dof_index> > block_list (fe_collection.size());
    for (unsigned int f=0; f<fe_collection.size(); ++f)
      {
        const FiniteElement<dim,spacedim> &fe = fe_collection[f];
        block_list[f].resize(fe.dofs_per_cell);
        for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
          block_list[f][i]
            = fe.system_to_block_index(i).first;
      }

    // set up a map where for each
    // block the respective degrees
    // of freedom are collected.
    //
    // note that this map is sorted by
    // block but that within each
    // block it is NOT sorted by
    // dof index. note also that some
    // dof indices are entered
    // multiply, so we will have to
    // take care of that
    std::vector<std::vector<types::global_dof_index> >
    block_to_dof_map (fe_collection.n_blocks());
    for (ITERATOR cell=start; cell!=end; ++cell)
      {
        if (is_level_operation)
          {
            //we are dealing with mg dofs, skip foreign level cells:
            if (!cell->is_locally_owned_on_level())
              continue;
          }
        else
          {
            //we are dealing with active dofs, skip the loop if not locally
            // owned:
            if (!cell->is_locally_owned())
              continue;
          }

        // on each cell: get dof indices
        // and insert them into the global
        // list using their component
        const unsigned int fe_index = cell->active_fe_index();
        const unsigned int dofs_per_cell =fe_collection[fe_index].dofs_per_cell;
        local_dof_indices.resize (dofs_per_cell);
        cell->get_active_or_mg_dof_indices (local_dof_indices);
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          if (start->get_dof_handler().locally_owned_dofs().is_element(local_dof_indices[i]))
            block_to_dof_map[block_list[fe_index][i]].
            push_back (local_dof_indices[i]);
      }

    // now we've got all indices sorted
    // into buckets labeled by their
    // target block number. we've
    // only got to traverse this list
    // and assign the new indices
    //
    // however, we first want to sort
    // the indices entered into the
    // buckets to preserve the order
    // within each component and during
    // this also remove duplicate
    // entries
    for (unsigned int block=0; block<fe_collection.n_blocks();
         ++block)
      {
        std::sort (block_to_dof_map[block].begin(),
                   block_to_dof_map[block].end());
        block_to_dof_map[block]
        .erase (std::unique (block_to_dof_map[block].begin(),
                             block_to_dof_map[block].end()),
                block_to_dof_map[block].end());
      }

    // calculate the number of locally owned
    // DoFs per bucket
    const unsigned int n_buckets = fe_collection.n_blocks();
    std::vector<types::global_dof_index> shifts(n_buckets);

    if (const parallel::Triangulation<dim,spacedim> *tria
        = (dynamic_cast<const parallel::Triangulation<dim,spacedim>*>
           (&start->get_dof_handler().get_triangulation())))
      {
#ifdef DEAL_II_WITH_MPI
        std::vector<types::global_dof_index> local_dof_count(n_buckets);

        for (unsigned int c=0; c<n_buckets; ++c)
          local_dof_count[c] = block_to_dof_map[c].size();


        // gather information from all CPUs
        std::vector<types::global_dof_index>
        all_dof_counts(fe_collection.n_components() *
                       Utilities::MPI::n_mpi_processes (tria->get_communicator()));

        MPI_Allgather ( &local_dof_count[0],
                        n_buckets, DEAL_II_DOF_INDEX_MPI_TYPE,
                        &all_dof_counts[0],
                        n_buckets, DEAL_II_DOF_INDEX_MPI_TYPE,
                        tria->get_communicator());

        for (unsigned int i=0; i<n_buckets; ++i)
          Assert (all_dof_counts[n_buckets*tria->locally_owned_subdomain()+i]
                  ==
                  local_dof_count[i],
                  ExcInternalError());

        //calculate shifts
        types::global_dof_index cumulated = 0;
        for (unsigned int c=0; c<n_buckets; ++c)
          {
            shifts[c]=cumulated;
            for (types::subdomain_id i=0; i<tria->locally_owned_subdomain(); ++i)
              shifts[c] += all_dof_counts[c+n_buckets*i];
            for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes (tria->get_communicator()); ++i)
              cumulated += all_dof_counts[c+n_buckets*i];
          }
#else
        (void)tria;
        Assert (false, ExcInternalError());
#endif
      }
    else
      {
        shifts[0] = 0;
        for (unsigned int c=1; c<fe_collection.n_blocks(); ++c)
          shifts[c] = shifts[c-1] + block_to_dof_map[c-1].size();
      }




    // now concatenate all the
    // components in the order the user
    // desired to see
    types::global_dof_index next_free_index = 0;
    for (unsigned int block=0; block<fe_collection.n_blocks(); ++block)
      {
        const typename std::vector<types::global_dof_index>::const_iterator
        begin_of_component = block_to_dof_map[block].begin(),
        end_of_component   = block_to_dof_map[block].end();

        next_free_index = shifts[block];

        for (typename std::vector<types::global_dof_index>::const_iterator
             dof_index = begin_of_component;
             dof_index != end_of_component; ++dof_index)
          {
            Assert (start->get_dof_handler().locally_owned_dofs()
                    .index_within_set(*dof_index)
                    <
                    new_indices.size(),
                    ExcInternalError());
            new_indices[start->get_dof_handler().locally_owned_dofs()
                        .index_within_set(*dof_index)]
              = next_free_index++;
          }
      }

    return next_free_index;
  }



  namespace
  {
    // helper function for hierarchical()

// Note that this function only works for active dofs.
    template <int dim, class iterator>
    types::global_dof_index
    compute_hierarchical_recursive (
      types::global_dof_index next_free,
      std::vector<types::global_dof_index> &new_indices,
      const iterator &cell,
      const IndexSet &locally_owned)
    {
      if (cell->has_children())
        {
          //recursion
          for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell; ++c)
            next_free = compute_hierarchical_recursive<dim> (
                          next_free,
                          new_indices,
                          cell->child (c),
                          locally_owned);
        }
      else
        {
          if (cell->is_locally_owned())
            {
              const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
              std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
              cell->get_dof_indices (local_dof_indices);

              for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                  if (locally_owned.is_element (local_dof_indices[i]))
                    {
                      // this is a locally owned DoF, assign new number if not assigned a number yet
                      unsigned int idx = locally_owned.index_within_set (local_dof_indices[i]);
                      if (new_indices[idx] == DoFHandler<dim>::invalid_dof_index)
                        {
                          new_indices[idx] = locally_owned.nth_index_in_set (next_free);
                          next_free++;
                        }
                    }
                }
            }
        }
      return next_free;
    }
  }



  template <int dim>
  void
  hierarchical (DoFHandler<dim> &dof_handler)
  {
    std::vector<types::global_dof_index> renumbering (dof_handler.n_locally_owned_dofs(),
                                                      DoFHandler<dim>::invalid_dof_index);

    typename DoFHandler<dim>::level_cell_iterator cell;

    types::global_dof_index next_free = 0;
    const IndexSet locally_owned = dof_handler.locally_owned_dofs();

    const parallel::distributed::Triangulation<dim> *tria
      = dynamic_cast<const parallel::distributed::Triangulation<dim>*>
        (&dof_handler.get_triangulation());

    if (tria)
      {
#ifdef DEAL_II_WITH_P4EST
        // this is a distributed Triangulation. We need to traverse the coarse
        // cells in the order p4est does
        for (unsigned int c = 0; c < tria->n_cells (0); ++c)
          {
            const unsigned int coarse_cell_index =
              tria->get_p4est_tree_to_coarse_cell_permutation() [c];

            const typename DoFHandler<dim>::level_cell_iterator
            this_cell (tria, 0, coarse_cell_index, &dof_handler);

            next_free = compute_hierarchical_recursive<dim> (next_free,
                                                             renumbering,
                                                             this_cell,
                                                             locally_owned);
          }
#else
        Assert (false, ExcNotImplemented());
#endif
      }
    else
      {
        //this is not a distributed Triangulation. Traverse coarse cells in the
        //normal order
        for (cell = dof_handler.begin (0); cell != dof_handler.end (0); ++cell)
          next_free = compute_hierarchical_recursive<dim> (next_free,
                                                           renumbering,
                                                           cell,
                                                           locally_owned);
      }

    // verify that the last numbered
    // degree of freedom is either
    // equal to the number of degrees
    // of freedom in total (the
    // sequential case) or in the
    // distributed case at least
    // makes sense
    Assert ((next_free == dof_handler.n_locally_owned_dofs())
            ||
            ((dof_handler.n_locally_owned_dofs() < dof_handler.n_dofs())
             &&
             (next_free <= dof_handler.n_dofs())),
            ExcRenumberingIncomplete());

    // make sure that all local DoFs got new numbers assigned
    Assert (std::find (renumbering.begin(), renumbering.end(),
                       numbers::invalid_dof_index)
            == renumbering.end(),
            ExcInternalError());

    dof_handler.renumber_dofs(renumbering);
  }



  template <typename DoFHandlerType>
  void
  sort_selected_dofs_back (DoFHandlerType          &dof_handler,
                           const std::vector<bool> &selected_dofs)
  {
    std::vector<types::global_dof_index> renumbering(dof_handler.n_dofs(),
                                                     DoFHandlerType::invalid_dof_index);
    compute_sort_selected_dofs_back(renumbering, dof_handler, selected_dofs);

    dof_handler.renumber_dofs(renumbering);
  }



  template <typename DoFHandlerType>
  void
  sort_selected_dofs_back (DoFHandlerType          &dof_handler,
                           const std::vector<bool> &selected_dofs,
                           const unsigned int       level)
  {
    Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
           ExcNotInitialized());

    std::vector<types::global_dof_index> renumbering(dof_handler.n_dofs(level),
                                                     DoFHandlerType::invalid_dof_index);
    compute_sort_selected_dofs_back(renumbering, dof_handler, selected_dofs, level);

    dof_handler.renumber_dofs(level, renumbering);
  }



  template <typename DoFHandlerType>
  void
  compute_sort_selected_dofs_back (std::vector<types::global_dof_index> &new_indices,
                                   const DoFHandlerType                 &dof_handler,
                                   const std::vector<bool>              &selected_dofs)
  {
    const types::global_dof_index n_dofs = dof_handler.n_dofs();
    Assert (selected_dofs.size() == n_dofs,
            ExcDimensionMismatch (selected_dofs.size(), n_dofs));

    // re-sort the dofs according to
    // their selection state
    Assert (new_indices.size() == n_dofs,
            ExcDimensionMismatch(new_indices.size(), n_dofs));

    const types::global_dof_index   n_selected_dofs = std::count (selected_dofs.begin(),
                                                      selected_dofs.end(),
                                                      false);

    types::global_dof_index next_unselected = 0;
    types::global_dof_index next_selected   = n_selected_dofs;
    for (types::global_dof_index i=0; i<n_dofs; ++i)
      if (selected_dofs[i] == false)
        {
          new_indices[i] = next_unselected;
          ++next_unselected;
        }
      else
        {
          new_indices[i] = next_selected;
          ++next_selected;
        };
    Assert (next_unselected == n_selected_dofs, ExcInternalError());
    Assert (next_selected == n_dofs, ExcInternalError());
  }



  template <typename DoFHandlerType>
  void
  compute_sort_selected_dofs_back (std::vector<types::global_dof_index> &new_indices,
                                   const DoFHandlerType                 &dof_handler,
                                   const std::vector<bool>              &selected_dofs,
                                   const unsigned int                    level)
  {
    Assert(dof_handler.n_dofs(level) != numbers::invalid_dof_index,
           ExcNotInitialized());

    const unsigned int n_dofs = dof_handler.n_dofs(level);
    Assert (selected_dofs.size() == n_dofs,
            ExcDimensionMismatch (selected_dofs.size(), n_dofs));

    // re-sort the dofs according to
    // their selection state
    Assert (new_indices.size() == n_dofs,
            ExcDimensionMismatch(new_indices.size(), n_dofs));

    const unsigned int   n_selected_dofs = std::count (selected_dofs.begin(),
                                                       selected_dofs.end(),
                                                       false);

    unsigned int next_unselected = 0;
    unsigned int next_selected   = n_selected_dofs;
    for (unsigned int i=0; i<n_dofs; ++i)
      if (selected_dofs[i] == false)
        {
          new_indices[i] = next_unselected;
          ++next_unselected;
        }
      else
        {
          new_indices[i] = next_selected;
          ++next_selected;
        };
    Assert (next_unselected == n_selected_dofs, ExcInternalError());
    Assert (next_selected == n_dofs, ExcInternalError());
  }



  template <typename DoFHandlerType>
  void
  cell_wise (DoFHandlerType &dof,
             const std::vector<typename DoFHandlerType::active_cell_iterator> &cells)
  {
    std::vector<types::global_dof_index> renumbering(dof.n_dofs());
    std::vector<types::global_dof_index> reverse(dof.n_dofs());
    compute_cell_wise(renumbering, reverse, dof, cells);

    dof.renumber_dofs(renumbering);
  }


  template <typename DoFHandlerType>
  void
  compute_cell_wise
  (std::vector<types::global_dof_index>                                      &new_indices,
   std::vector<types::global_dof_index>                                      &reverse,
   const DoFHandlerType                                                      &dof,
   const typename std::vector<typename DoFHandlerType::active_cell_iterator> &cells)
  {
    Assert(cells.size() == dof.get_triangulation().n_active_cells(),
           ExcDimensionMismatch(cells.size(),
                                dof.get_triangulation().n_active_cells()));

    types::global_dof_index n_global_dofs = dof.n_dofs();

    // Actually, we compute the
    // inverse of the reordering
    // vector, called reverse here.
    // Later, irs inverse is computed
    // into new_indices, which is the
    // return argument.

    Assert(new_indices.size() == n_global_dofs,
           ExcDimensionMismatch(new_indices.size(), n_global_dofs));
    Assert(reverse.size() == n_global_dofs,
           ExcDimensionMismatch(reverse.size(), n_global_dofs));

    // For continuous elements, we must
    // make sure, that each dof is
    // reordered only once.
    std::vector<bool> already_sorted(n_global_dofs, false);
    std::vector<types::global_dof_index> cell_dofs;

    unsigned int global_index = 0;

    typename std::vector<typename DoFHandlerType::active_cell_iterator>::const_iterator cell;

    for (cell = cells.begin(); cell != cells.end(); ++cell)
      {
        // Determine the number of dofs
        // on this cell and reinit the
        // vector storing these
        // numbers.
        unsigned int n_cell_dofs = (*cell)->get_fe().n_dofs_per_cell();
        cell_dofs.resize(n_cell_dofs);

        (*cell)->get_active_or_mg_dof_indices(cell_dofs);

        // Sort here to make sure that
        // degrees of freedom inside a
        // single cell are in the same
        // order after renumbering.
        std::sort(cell_dofs.begin(), cell_dofs.end());

        for (unsigned int i=0; i<n_cell_dofs; ++i)
          {
            if (!already_sorted[cell_dofs[i]])
              {
                already_sorted[cell_dofs[i]] = true;
                reverse[global_index++] = cell_dofs[i];
              }
          }
      }
    Assert(global_index == n_global_dofs, ExcRenumberingIncomplete());

    for (types::global_dof_index i=0; i<reverse.size(); ++i)
      new_indices[reverse[i]] = i;
  }



  template <typename DoFHandlerType>
  void cell_wise
  (DoFHandlerType                                                           &dof,
   const unsigned int                                                        level,
   const typename std::vector<typename DoFHandlerType::level_cell_iterator> &cells)
  {
    Assert(dof.n_dofs(level) != numbers::invalid_dof_index,
           ExcNotInitialized());

    std::vector<types::global_dof_index> renumbering(dof.n_dofs(level));
    std::vector<types::global_dof_index> reverse(dof.n_dofs(level));

    compute_cell_wise(renumbering, reverse, dof, level, cells);
    dof.renumber_dofs(level, renumbering);
  }



  template <typename DoFHandlerType>
  void compute_cell_wise
  (std::vector<types::global_dof_index>                                     &new_order,
   std::vector<types::global_dof_index>                                     &reverse,
   const DoFHandlerType                                                     &dof,
   const unsigned int                                                        level,
   const typename std::vector<typename DoFHandlerType::level_cell_iterator> &cells)
  {
    Assert(cells.size() == dof.get_triangulation().n_cells(level),
           ExcDimensionMismatch(cells.size(),
                                dof.get_triangulation().n_cells(level)));
    Assert (new_order.size() == dof.n_dofs(level),
            ExcDimensionMismatch(new_order.size(), dof.n_dofs(level)));
    Assert (reverse.size() == dof.n_dofs(level),
            ExcDimensionMismatch(reverse.size(), dof.n_dofs(level)));

    unsigned int n_global_dofs = dof.n_dofs(level);
    unsigned int n_cell_dofs = dof.get_fe().n_dofs_per_cell();

    std::vector<bool> already_sorted(n_global_dofs, false);
    std::vector<types::global_dof_index> cell_dofs(n_cell_dofs);

    unsigned int global_index = 0;

    typename std::vector<typename DoFHandlerType::level_cell_iterator>::const_iterator cell;

    for (cell = cells.begin(); cell != cells.end(); ++cell)
      {
        Assert ((*cell)->level() == (int) level, ExcInternalError());

        (*cell)->get_active_or_mg_dof_indices(cell_dofs);
        std::sort(cell_dofs.begin(), cell_dofs.end());

        for (unsigned int i=0; i<n_cell_dofs; ++i)
          {
            if (!already_sorted[cell_dofs[i]])
              {
                already_sorted[cell_dofs[i]] = true;
                reverse[global_index++] = cell_dofs[i];
              }
          }
      }
    Assert(global_index == n_global_dofs, ExcRenumberingIncomplete());

    for (types::global_dof_index i=0; i<new_order.size(); ++i)
      new_order[reverse[i]] = i;
  }







  template <typename DoFHandlerType>
  void
  compute_downstream
  (std::vector<types::global_dof_index>         &new_indices,
   std::vector<types::global_dof_index>         &reverse,
   const DoFHandlerType                         &dof,
   const Point<DoFHandlerType::space_dimension> &direction,
   const bool                                    dof_wise_renumbering)
  {
    if (dof_wise_renumbering == false)
      {
        std::vector<typename DoFHandlerType::active_cell_iterator> ordered_cells;
        ordered_cells.reserve(dof.get_triangulation().n_active_cells());
        const CompareDownstream<typename DoFHandlerType::active_cell_iterator,
              DoFHandlerType::space_dimension> comparator(direction);

        typename DoFHandlerType::active_cell_iterator p = dof.begin_active();
        typename DoFHandlerType::active_cell_iterator end = dof.end();

        while (p!=end)
          {
            ordered_cells.push_back(p);
            ++p;
          }
        std::sort (ordered_cells.begin(), ordered_cells.end(), comparator);

        compute_cell_wise(new_indices, reverse, dof, ordered_cells);
      }
    else
      {
        // similar code as for
        // DoFTools::map_dofs_to_support_points, but
        // need to do this for general DoFHandlerType classes and
        // want to be able to sort the result
        // (otherwise, could use something like
        // DoFTools::map_support_points_to_dofs)
        const unsigned int n_dofs = dof.n_dofs();
        std::vector<std::pair<Point<DoFHandlerType::space_dimension>,unsigned int> > support_point_list
        (n_dofs);

        const hp::FECollection<DoFHandlerType::dimension> fe_collection (dof.get_fe ());
        Assert (fe_collection[0].has_support_points(),
                typename FiniteElement<DoFHandlerType::dimension>::ExcFEHasNoSupportPoints());
        hp::QCollection<DoFHandlerType::dimension> quadrature_collection;
        for (unsigned int comp=0; comp<fe_collection.size(); ++comp)
          {
            Assert (fe_collection[comp].has_support_points(),
                    typename FiniteElement<DoFHandlerType::dimension>::ExcFEHasNoSupportPoints());
            quadrature_collection.push_back
            (Quadrature<DoFHandlerType::dimension> (fe_collection[comp].
                                                    get_unit_support_points()));
          }
        hp::FEValues<DoFHandlerType::dimension,DoFHandlerType::space_dimension>
        hp_fe_values (fe_collection, quadrature_collection,
                      update_quadrature_points);

        std::vector<bool> already_touched (n_dofs, false);

        std::vector<types::global_dof_index> local_dof_indices;
        typename DoFHandlerType::active_cell_iterator begin = dof.begin_active();
        typename DoFHandlerType::active_cell_iterator end = dof.end();
        for ( ; begin != end; ++begin)
          {
            const unsigned int dofs_per_cell = begin->get_fe().dofs_per_cell;
            local_dof_indices.resize (dofs_per_cell);
            hp_fe_values.reinit (begin);
            const FEValues<DoFHandlerType::dimension> &fe_values =
              hp_fe_values.get_present_fe_values ();
            begin->get_active_or_mg_dof_indices(local_dof_indices);
            const std::vector<Point<DoFHandlerType::space_dimension> > &points
              = fe_values.get_quadrature_points ();
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              if (!already_touched[local_dof_indices[i]])
                {
                  support_point_list[local_dof_indices[i]].first = points[i];
                  support_point_list[local_dof_indices[i]].second =
                    local_dof_indices[i];
                  already_touched[local_dof_indices[i]] = true;
                }
          }

        ComparePointwiseDownstream<DoFHandlerType::space_dimension> comparator (direction);
        std::sort (support_point_list.begin(), support_point_list.end(),
                   comparator);
        for (types::global_dof_index i=0; i<n_dofs; ++i)
          new_indices[support_point_list[i].second] = i;
      }
  }



  template <typename DoFHandlerType>
  void downstream (DoFHandlerType                               &dof,
                   const unsigned int                            level,
                   const Point<DoFHandlerType::space_dimension> &direction,
                   const bool                                    dof_wise_renumbering)
  {
    std::vector<types::global_dof_index> renumbering(dof.n_dofs(level));
    std::vector<types::global_dof_index> reverse(dof.n_dofs(level));
    compute_downstream(renumbering, reverse, dof, level, direction,
                       dof_wise_renumbering);

    dof.renumber_dofs(level, renumbering);
  }



  template <typename DoFHandlerType>
  void
  compute_downstream
  (std::vector<types::global_dof_index>         &new_indices,
   std::vector<types::global_dof_index>         &reverse,
   const DoFHandlerType                         &dof,
   const unsigned int                            level,
   const Point<DoFHandlerType::space_dimension> &direction,
   const bool                                    dof_wise_renumbering)
  {
    if (dof_wise_renumbering == false)
      {
        std::vector<typename DoFHandlerType::level_cell_iterator> ordered_cells;
        ordered_cells.reserve (dof.get_triangulation().n_cells(level));
        const CompareDownstream<typename DoFHandlerType::level_cell_iterator,
              DoFHandlerType::space_dimension> comparator(direction);

        typename DoFHandlerType::level_cell_iterator p = dof.begin(level);
        typename DoFHandlerType::level_cell_iterator end = dof.end(level);

        while (p!=end)
          {
            ordered_cells.push_back(p);
            ++p;
          }
        std::sort (ordered_cells.begin(), ordered_cells.end(), comparator);

        compute_cell_wise(new_indices, reverse, dof, level, ordered_cells);
      }
    else
      {
        Assert (dof.get_fe().has_support_points(),
                typename FiniteElement<DoFHandlerType::dimension>::ExcFEHasNoSupportPoints());
        const unsigned int n_dofs = dof.n_dofs(level);
        std::vector<std::pair<Point<DoFHandlerType::space_dimension>,unsigned int> > support_point_list
        (n_dofs);

        Quadrature<DoFHandlerType::dimension>   q_dummy(dof.get_fe().get_unit_support_points());
        FEValues<DoFHandlerType::dimension,DoFHandlerType::space_dimension> fe_values (dof.get_fe(), q_dummy,
            update_quadrature_points);

        std::vector<bool> already_touched (dof.n_dofs(), false);

        const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
        std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
        typename DoFHandlerType::level_cell_iterator begin = dof.begin(level);
        typename DoFHandlerType::level_cell_iterator end = dof.end(level);
        for ( ; begin != end; ++begin)
          {
            const typename Triangulation<DoFHandlerType::dimension,
                  DoFHandlerType::space_dimension>::cell_iterator &begin_tria = begin;
            begin->get_active_or_mg_dof_indices(local_dof_indices);
            fe_values.reinit (begin_tria);
            const std::vector<Point<DoFHandlerType::space_dimension> > &points
              = fe_values.get_quadrature_points ();
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              if (!already_touched[local_dof_indices[i]])
                {
                  support_point_list[local_dof_indices[i]].first = points[i];
                  support_point_list[local_dof_indices[i]].second =
                    local_dof_indices[i];
                  already_touched[local_dof_indices[i]] = true;
                }
          }

        ComparePointwiseDownstream<DoFHandlerType::space_dimension> comparator (direction);
        std::sort (support_point_list.begin(), support_point_list.end(),
                   comparator);
        for (types::global_dof_index i=0; i<n_dofs; ++i)
          new_indices[support_point_list[i].second] = i;
      }
  }



  /**
   * Provide comparator for DoFCellAccessors
   */
  namespace internal
  {
    template <int dim>
    struct ClockCells
    {
      /**
       * Center of rotation.
       */
      const Point<dim> &center;
      /**
       * Revert sorting order.
       */
      bool counter;

      /**
       * Constructor.
       */
      ClockCells (const Point<dim> &center, bool counter) :
        center(center),
        counter(counter)
      {}

      /**
       * Comparison operator
       */
      template <class DHCellIterator>
      bool operator () (const DHCellIterator &c1,
                        const DHCellIterator &c2) const
      {
        // dispatch to
        // dimension-dependent functions
        return compare (c1, c2, dealii::internal::int2type<dim>());
      }

    private:
      /**
       * Comparison operator for dim>=2
       */
      template <class DHCellIterator, int xdim>
      bool compare (const DHCellIterator &c1,
                    const DHCellIterator &c2,
                    dealii::internal::int2type<xdim>) const
      {
        const Tensor<1,dim> v1 = c1->center() - center;
        const Tensor<1,dim> v2 = c2->center() - center;
        const double s1 = std::atan2(v1[0], v1[1]);
        const double s2 = std::atan2(v2[0], v2[1]);
        return ( counter ? (s1>s2) : (s2>s1));
      }


      /**
       * Comparison operator for dim==1
       * where this function makes no sense
       */
      template <class DHCellIterator>
      bool compare (const DHCellIterator &,
                    const DHCellIterator &,
                    dealii::internal::int2type<1>) const
      {
        Assert (dim >= 2,
                ExcMessage ("This operation only makes sense for dim>=2."));
        return false;
      }

    };
  }



  template <typename DoFHandlerType>
  void
  clockwise_dg (
    DoFHandlerType &dof,
    const Point<DoFHandlerType::space_dimension> &center,
    const bool counter)
  {
    std::vector<types::global_dof_index> renumbering(dof.n_dofs());
    compute_clockwise_dg(renumbering, dof, center, counter);

    dof.renumber_dofs(renumbering);
  }



  template <typename DoFHandlerType>
  void
  compute_clockwise_dg
  (std::vector<types::global_dof_index>         &new_indices,
   const DoFHandlerType                         &dof,
   const Point<DoFHandlerType::space_dimension> &center,
   const bool                                    counter)
  {
    std::vector<typename DoFHandlerType::active_cell_iterator> ordered_cells;
    ordered_cells.reserve (dof.get_triangulation().n_active_cells());
    internal::ClockCells<DoFHandlerType::space_dimension> comparator(center, counter);

    typename DoFHandlerType::active_cell_iterator p = dof.begin_active();
    typename DoFHandlerType::active_cell_iterator end = dof.end();

    while (p!=end)
      {
        ordered_cells.push_back(p);
        ++p;
      }
    std::sort (ordered_cells.begin(), ordered_cells.end(), comparator);

    std::vector<types::global_dof_index> reverse(new_indices.size());
    compute_cell_wise(new_indices, reverse, dof, ordered_cells);
  }



  template <typename DoFHandlerType>
  void clockwise_dg (DoFHandlerType                               &dof,
                     const unsigned int                            level,
                     const Point<DoFHandlerType::space_dimension> &center,
                     const bool                                    counter)
  {
    std::vector<typename DoFHandlerType::level_cell_iterator> ordered_cells;
    ordered_cells.reserve(dof.get_triangulation().n_active_cells());
    internal::ClockCells<DoFHandlerType::space_dimension> comparator(center, counter);

    typename DoFHandlerType::level_cell_iterator p = dof.begin(level);
    typename DoFHandlerType::level_cell_iterator end = dof.end(level);

    while (p!=end)
      {
        ordered_cells.push_back(p);
        ++p;
      }
    std::sort (ordered_cells.begin(), ordered_cells.end(), comparator);

    cell_wise(dof, level, ordered_cells);
  }



  template <typename DoFHandlerType>
  void
  random (DoFHandlerType &dof_handler)
  {
    std::vector<types::global_dof_index> renumbering(dof_handler.n_dofs(),
                                                     DoFHandlerType::invalid_dof_index);
    compute_random(renumbering, dof_handler);

    dof_handler.renumber_dofs(renumbering);
  }



  template <typename DoFHandlerType>
  void
  compute_random (
    std::vector<types::global_dof_index> &new_indices,
    const DoFHandlerType      &dof_handler)
  {
    const types::global_dof_index n_dofs = dof_handler.n_dofs();
    Assert(new_indices.size() == n_dofs,
           ExcDimensionMismatch(new_indices.size(), n_dofs));

    for (unsigned int i=0; i<n_dofs; ++i)
      new_indices[i] = i;

    // shuffle the elements; the following is essentially the
    // std::random_shuffle algorithm but uses a predictable
    // random number generator
    ::boost::mt19937 random_number_generator;
    for (unsigned int i=1; i<n_dofs; ++i)
      {
        // get a random number between 0 and i (inclusive)
        const unsigned int j
          = ::boost::random::uniform_int_distribution<>(0, i)(random_number_generator);

        // if possible, swap the elements
        if (i != j)
          std::swap (new_indices[i], new_indices[j]);
      }
  }



  template <typename DoFHandlerType>
  void
  subdomain_wise (DoFHandlerType &dof_handler)
  {
    std::vector<types::global_dof_index> renumbering(dof_handler.n_dofs(),
                                                     DoFHandlerType::invalid_dof_index);
    compute_subdomain_wise(renumbering, dof_handler);

    dof_handler.renumber_dofs(renumbering);
  }



  template <typename DoFHandlerType>
  void
  compute_subdomain_wise (std::vector<types::global_dof_index> &new_dof_indices,
                          const DoFHandlerType      &dof_handler)
  {
    const types::global_dof_index n_dofs = dof_handler.n_dofs();
    Assert (new_dof_indices.size() == n_dofs,
            ExcDimensionMismatch (new_dof_indices.size(), n_dofs));

    // first get the association of each dof
    // with a subdomain and determine the total
    // number of subdomain ids used
    std::vector<types::subdomain_id> subdomain_association (n_dofs);
    DoFTools::get_subdomain_association (dof_handler,
                                         subdomain_association);
    const unsigned int n_subdomains
      = *std::max_element (subdomain_association.begin(),
                           subdomain_association.end()) + 1;

    // then renumber the subdomains by first
    // looking at those belonging to subdomain
    // 0, then those of subdomain 1, etc. note
    // that the algorithm is stable, i.e. if
    // two dofs i,j have i<j and belong to the
    // same subdomain, then they will be in
    // this order also after reordering
    std::fill (new_dof_indices.begin(), new_dof_indices.end(),
               numbers::invalid_dof_index);
    types::global_dof_index next_free_index = 0;
    for (types::subdomain_id subdomain=0; subdomain<n_subdomains; ++subdomain)
      for (types::global_dof_index i=0; i<n_dofs; ++i)
        if (subdomain_association[i] == subdomain)
          {
            Assert (new_dof_indices[i] == numbers::invalid_dof_index,
                    ExcInternalError());
            new_dof_indices[i] = next_free_index;
            ++next_free_index;
          }

    // we should have numbered all dofs
    Assert (next_free_index == n_dofs, ExcInternalError());
    Assert (std::find (new_dof_indices.begin(), new_dof_indices.end(),
                       numbers::invalid_dof_index)
            == new_dof_indices.end(),
            ExcInternalError());
  }

} // namespace DoFRenumbering



/*-------------- Explicit Instantiations -------------------------------*/
#include "dof_renumbering.inst"


DEAL_II_NAMESPACE_CLOSE
