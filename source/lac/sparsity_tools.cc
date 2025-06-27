// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/exceptions.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <algorithm>
#include <functional>
#include <memory>
#include <set>

#ifdef DEAL_II_WITH_MPI
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/utilities.h>

#  include <deal.II/lac/block_sparsity_pattern.h>
#  include <deal.II/lac/dynamic_sparsity_pattern.h>
#endif

#ifdef DEAL_II_WITH_METIS
extern "C"
{
#  include <metis.h>
}
#endif

#ifdef DEAL_II_TRILINOS_WITH_ZOLTAN
#  include <zoltan_cpp.h>
#endif

#include <string>


DEAL_II_NAMESPACE_OPEN

namespace SparsityTools
{
  namespace
  {
    void
    partition_metis(const SparsityPattern           &sparsity_pattern,
                    const std::vector<unsigned int> &cell_weights,
                    const unsigned int               n_partitions,
                    std::vector<unsigned int>       &partition_indices)
    {
      // Make sure that METIS is actually
      // installed and detected
#ifndef DEAL_II_WITH_METIS
      (void)sparsity_pattern;
      (void)cell_weights;
      (void)n_partitions;
      (void)partition_indices;
      AssertThrow(false, ExcMETISNotInstalled());
#else

      // Generate the data structures for METIS. Note that this is particularly
      // simple, since METIS wants exactly our compressed row storage format.
      // We only have to set up a few auxiliary arrays and convert from our
      // unsigned cell weights to signed ones.
      idx_t n = static_cast<signed int>(sparsity_pattern.n_rows());

      idx_t ncon = 1; // number of balancing constraints (should be >0)

      // We can not partition n items into more than n parts. METIS will
      // generate non-sensical output (everything is owned by a single process)
      // and complain with a message (but won't return an error code!):
      // ***Cannot bisect a graph with 0 vertices!
      // ***You are trying to partition a graph into too many parts!
      idx_t nparts =
        std::min(n,
                 static_cast<idx_t>(
                   n_partitions)); // number of subdomains to create

      // use default options for METIS
      idx_t options[METIS_NOPTIONS];
      METIS_SetDefaultOptions(options);

      // one more nuisance: we have to copy our own data to arrays that store
      // signed integers :-(
      std::vector<idx_t> int_rowstart(1);
      int_rowstart.reserve(sparsity_pattern.n_rows() + 1);
      std::vector<idx_t> int_colnums;
      int_colnums.reserve(sparsity_pattern.n_nonzero_elements());
      for (SparsityPattern::size_type row = 0; row < sparsity_pattern.n_rows();
           ++row)
        {
          for (SparsityPattern::iterator col = sparsity_pattern.begin(row);
               col < sparsity_pattern.end(row);
               ++col)
            int_colnums.push_back(col->column());
          int_rowstart.push_back(int_colnums.size());
        }

      std::vector<idx_t> int_partition_indices(sparsity_pattern.n_rows());

      // Set up cell weighting option
      std::vector<idx_t> int_cell_weights;
      if (cell_weights.size() > 0)
        {
          Assert(cell_weights.size() == sparsity_pattern.n_rows(),
                 ExcDimensionMismatch(cell_weights.size(),
                                      sparsity_pattern.n_rows()));
          int_cell_weights.resize(cell_weights.size());
          std::copy(cell_weights.begin(),
                    cell_weights.end(),
                    int_cell_weights.begin());
        }
      // Set a pointer to the optional cell weighting information.
      // METIS expects a null pointer if there are no weights to be considered.
      idx_t *const p_int_cell_weights =
        (cell_weights.size() > 0 ? int_cell_weights.data() : nullptr);


      // Make use of METIS' error code.
      int ierr;

      // Select which type of partitioning to create

      // Use recursive if the number of partitions is less than or equal to 8
      idx_t dummy; // output: # of edges cut by the resulting partition
      if (nparts <= 8)
        ierr = METIS_PartGraphRecursive(&n,
                                        &ncon,
                                        int_rowstart.data(),
                                        int_colnums.data(),
                                        p_int_cell_weights,
                                        nullptr,
                                        nullptr,
                                        &nparts,
                                        nullptr,
                                        nullptr,
                                        options,
                                        &dummy,
                                        int_partition_indices.data());

      // Otherwise use kway
      else
        ierr = METIS_PartGraphKway(&n,
                                   &ncon,
                                   int_rowstart.data(),
                                   int_colnums.data(),
                                   p_int_cell_weights,
                                   nullptr,
                                   nullptr,
                                   &nparts,
                                   nullptr,
                                   nullptr,
                                   options,
                                   &dummy,
                                   int_partition_indices.data());

      // If metis returns normally, an error code METIS_OK=1 is returned from
      // the above functions (see metish.h)
      AssertThrow(ierr == 1, ExcMETISError(ierr));

      // now copy back generated indices into the output array
      std::copy(int_partition_indices.begin(),
                int_partition_indices.end(),
                partition_indices.begin());
#endif
    }


// Query functions unused if zoltan is not installed
#ifdef DEAL_II_TRILINOS_WITH_ZOLTAN
    // Query functions for partition_zoltan
    int
    get_number_of_objects(void *data, int *ierr)
    {
      SparsityPattern *graph = reinterpret_cast<SparsityPattern *>(data);

      *ierr = ZOLTAN_OK;

      return graph->n_rows();
    }


    void
    get_object_list(void *data,
                    int /*sizeGID*/,
                    int /*sizeLID*/,
                    ZOLTAN_ID_PTR globalID,
                    ZOLTAN_ID_PTR localID,
                    int /*wgt_dim*/,
                    float * /*obj_wgts*/,
                    int *ierr)
    {
      SparsityPattern *graph = reinterpret_cast<SparsityPattern *>(data);
      *ierr                  = ZOLTAN_OK;

      Assert(globalID != nullptr, ExcInternalError());
      Assert(localID != nullptr, ExcInternalError());

      // set global degrees of freedom
      auto n_dofs = graph->n_rows();

      for (unsigned int i = 0; i < n_dofs; ++i)
        {
          globalID[i] = i;
          localID[i]  = i; // Same as global ids.
        }
    }


    void
    get_num_edges_list(void *data,
                       int /*sizeGID*/,
                       int /*sizeLID*/,
                       int           num_obj,
                       ZOLTAN_ID_PTR globalID,
                       ZOLTAN_ID_PTR /*localID*/,
                       int *numEdges,
                       int *ierr)
    {
      SparsityPattern *graph = reinterpret_cast<SparsityPattern *>(data);

      *ierr = ZOLTAN_OK;

      Assert(numEdges != nullptr, ExcInternalError());

      for (int i = 0; i < num_obj; ++i)
        {
          if (graph->exists(i, i)) // Check if diagonal element is present
            numEdges[i] = graph->row_length(globalID[i]) - 1;
          else
            numEdges[i] = graph->row_length(globalID[i]);
        }
    }



    void
    get_edge_list(void *data,
                  int /*sizeGID*/,
                  int /*sizeLID*/,
                  int num_obj,
                  ZOLTAN_ID_PTR /*globalID*/,
                  ZOLTAN_ID_PTR /*localID*/,
                  int * /*num_edges*/,
                  ZOLTAN_ID_PTR nborGID,
                  int          *nborProc,
                  int /*wgt_dim*/,
                  float * /*ewgts*/,
                  int *ierr)
    {
      SparsityPattern *graph = reinterpret_cast<SparsityPattern *>(data);
      *ierr                  = ZOLTAN_OK;

      ZOLTAN_ID_PTR nextNborGID  = nborGID;
      int          *nextNborProc = nborProc;

      // Loop through rows corresponding to indices in globalID implicitly
      for (SparsityPattern::size_type i = 0;
           i < static_cast<SparsityPattern::size_type>(num_obj);
           ++i)
        {
          // Loop through each column to find neighbours
          for (SparsityPattern::iterator col = graph->begin(i);
               col < graph->end(i);
               ++col)
            // Ignore diagonal entries. Not needed for partitioning.
            if (i != col->column())
              {
                Assert(nextNborGID != nullptr, ExcInternalError());
                Assert(nextNborProc != nullptr, ExcInternalError());

                *nextNborGID++  = col->column();
                *nextNborProc++ = 0; // All the vertices on processor 0
              }
        }
    }
#endif


    void
    partition_zoltan(const SparsityPattern           &sparsity_pattern,
                     const std::vector<unsigned int> &cell_weights,
                     const unsigned int               n_partitions,
                     std::vector<unsigned int>       &partition_indices)
    {
      // Make sure that ZOLTAN is actually
      // installed and detected
#ifndef DEAL_II_TRILINOS_WITH_ZOLTAN
      (void)sparsity_pattern;
      (void)cell_weights;
      (void)n_partitions;
      (void)partition_indices;
      AssertThrow(false, ExcZOLTANNotInstalled());
#else

      Assert(
        cell_weights.empty(),
        ExcMessage(
          "The cell weighting functionality for Zoltan has not yet been implemented."));

      // MPI environment must have been initialized by this point.
      std::unique_ptr<Zoltan> zz = std::make_unique<Zoltan>(MPI_COMM_SELF);

      // General parameters
      // DEBUG_LEVEL call must precede the call to LB_METHOD
      zz->Set_Param("DEBUG_LEVEL", "0"); // set level of debug info
      zz->Set_Param(
        "LB_METHOD",
        "GRAPH"); // graph based partition method (LB-load balancing)
      zz->Set_Param("NUM_LOCAL_PARTS",
                    std::to_string(n_partitions)); // set number of partitions

      // The PHG partitioner is a hypergraph partitioner that Zoltan could use
      // for graph partitioning.
      // If number of vertices in hyperedge divided by total vertices in
      // hypergraph exceeds PHG_EDGE_SIZE_THRESHOLD,
      // then the hyperedge will be omitted as such (dense) edges will likely
      // incur high communication costs regardless of the partition.
      // PHG_EDGE_SIZE_THRESHOLD value is raised to 0.5 from the default
      // value of 0.25 so that the PHG partitioner doesn't throw warning saying
      // "PHG_EDGE_SIZE_THRESHOLD is low ..." after removing all dense edges.
      // For instance, in two dimensions if the triangulation being partitioned
      // is two quadrilaterals sharing an edge and if PHG_EDGE_SIZE_THRESHOLD
      // value is set to 0.25, PHG will remove all the edges throwing the
      // above warning.
      zz->Set_Param("PHG_EDGE_SIZE_THRESHOLD", "0.5");

      // Need a non-const object equal to sparsity_pattern
      SparsityPattern graph;
      graph.copy_from(sparsity_pattern);

      // Set query functions
      zz->Set_Num_Obj_Fn(get_number_of_objects, &graph);
      zz->Set_Obj_List_Fn(get_object_list, &graph);
      zz->Set_Num_Edges_Multi_Fn(get_num_edges_list, &graph);
      zz->Set_Edge_List_Multi_Fn(get_edge_list, &graph);

      // Variables needed by partition function
      int           changes           = 0;
      int           num_gid_entries   = 1;
      int           num_lid_entries   = 1;
      int           num_import        = 0;
      ZOLTAN_ID_PTR import_global_ids = nullptr;
      ZOLTAN_ID_PTR import_local_ids  = nullptr;
      int          *import_procs      = nullptr;
      int          *import_to_part    = nullptr;
      int           num_export        = 0;
      ZOLTAN_ID_PTR export_global_ids = nullptr;
      ZOLTAN_ID_PTR export_local_ids  = nullptr;
      int          *export_procs      = nullptr;
      int          *export_to_part    = nullptr;

      // call partitioner
      const int rc = zz->LB_Partition(changes,
                                      num_gid_entries,
                                      num_lid_entries,
                                      num_import,
                                      import_global_ids,
                                      import_local_ids,
                                      import_procs,
                                      import_to_part,
                                      num_export,
                                      export_global_ids,
                                      export_local_ids,
                                      export_procs,
                                      export_to_part);

      // check for error code in partitioner
      Assert(rc == ZOLTAN_OK, ExcInternalError());

      // By default, all indices belong to part 0. After zoltan partition
      // some are migrated to different part ID, which is stored in
      // export_to_part array.
      std::fill(partition_indices.begin(), partition_indices.end(), 0);

      // copy from export_to_part to partition_indices, whose part_ids != 0.
      Assert(export_to_part != nullptr, ExcInternalError());
      for (int i = 0; i < num_export; ++i)
        partition_indices[export_local_ids[i]] = export_to_part[i];
#endif
    }
  } // namespace


  void
  partition(const SparsityPattern     &sparsity_pattern,
            const unsigned int         n_partitions,
            std::vector<unsigned int> &partition_indices,
            const Partitioner          partitioner)
  {
    std::vector<unsigned int> cell_weights;

    // Call the other more general function
    partition(sparsity_pattern,
              cell_weights,
              n_partitions,
              partition_indices,
              partitioner);
  }


  void
  partition(const SparsityPattern           &sparsity_pattern,
            const std::vector<unsigned int> &cell_weights,
            const unsigned int               n_partitions,
            std::vector<unsigned int>       &partition_indices,
            const Partitioner                partitioner)
  {
    Assert(sparsity_pattern.n_rows() == sparsity_pattern.n_cols(),
           ExcNotQuadratic());
    Assert(sparsity_pattern.is_compressed(),
           SparsityPattern::ExcNotCompressed());

    Assert(n_partitions > 0, ExcInvalidNumberOfPartitions(n_partitions));
    Assert(partition_indices.size() == sparsity_pattern.n_rows(),
           ExcInvalidArraySize(partition_indices.size(),
                               sparsity_pattern.n_rows()));

    // check for an easy return
    if (n_partitions == 1 || (sparsity_pattern.n_rows() == 1))
      {
        std::fill_n(partition_indices.begin(), partition_indices.size(), 0U);
        return;
      }

    if (partitioner == Partitioner::metis)
      partition_metis(sparsity_pattern,
                      cell_weights,
                      n_partitions,
                      partition_indices);
    else if (partitioner == Partitioner::zoltan)
      partition_zoltan(sparsity_pattern,
                       cell_weights,
                       n_partitions,
                       partition_indices);
    else
      AssertThrow(false, ExcInternalError());
  }


  unsigned int
  color_sparsity_pattern(const SparsityPattern     &sparsity_pattern,
                         std::vector<unsigned int> &color_indices)
  {
    // Make sure that ZOLTAN is actually
    // installed and detected
#ifndef DEAL_II_TRILINOS_WITH_ZOLTAN
    (void)sparsity_pattern;
    (void)color_indices;
    AssertThrow(false, ExcZOLTANNotInstalled());
    return 0;
#else
    // coloring algorithm is run in serial by each processor.
    std::unique_ptr<Zoltan> zz = std::make_unique<Zoltan>(MPI_COMM_SELF);

    // Coloring parameters
    // DEBUG_LEVEL must precede all other calls
    zz->Set_Param("DEBUG_LEVEL", "0");               // level of debug info
    zz->Set_Param("COLORING_PROBLEM", "DISTANCE-1"); // Standard coloring
    zz->Set_Param("NUM_GID_ENTRIES", "1"); // 1 entry represents global ID
    zz->Set_Param("NUM_LID_ENTRIES", "1"); // 1 entry represents local ID
    zz->Set_Param("OBJ_WEIGHT_DIM", "0");  // object weights not used
    zz->Set_Param("RECOLORING_NUM_OF_ITERATIONS", "0");

    // Zoltan::Color function requires a non-const SparsityPattern object
    SparsityPattern graph;
    graph.copy_from(sparsity_pattern);

    // Set query functions required by coloring function
    zz->Set_Num_Obj_Fn(get_number_of_objects, &graph);
    zz->Set_Obj_List_Fn(get_object_list, &graph);
    zz->Set_Num_Edges_Multi_Fn(get_num_edges_list, &graph);
    zz->Set_Edge_List_Multi_Fn(get_edge_list, &graph);

    // Variables needed by coloring function
    int       num_gid_entries = 1;
    const int num_objects     = graph.n_rows();

    // Preallocate input variables. Element type fixed by ZOLTAN.
    std::vector<ZOLTAN_ID_TYPE> global_ids(num_objects);
    std::vector<int>            color_exp(num_objects);

    // Set ids for which coloring needs to be done
    for (int i = 0; i < num_objects; ++i)
      global_ids[i] = i;

    // Call ZOLTAN coloring algorithm
    int rc = zz->Color(num_gid_entries,
                       num_objects,
                       global_ids.data(),
                       color_exp.data());
    // Check for error code
    Assert(rc == ZOLTAN_OK, ExcInternalError());

    // Allocate and assign color indices
    color_indices.resize(num_objects);
    Assert(color_exp.size() == color_indices.size(),
           ExcDimensionMismatch(color_exp.size(), color_indices.size()));

    std::copy(color_exp.begin(), color_exp.end(), color_indices.begin());

    unsigned int n_colors =
      *(std::max_element(color_indices.begin(), color_indices.end()));
    return n_colors;
#endif
  }


  namespace internal
  {
    /**
     * Given a connectivity graph and a list of indices (where
     * invalid_size_type indicates that a node has not been numbered yet),
     * pick a valid starting index among the as-yet unnumbered one.
     */
    DynamicSparsityPattern::size_type
    find_unnumbered_starting_index(
      const DynamicSparsityPattern                         &sparsity,
      const std::vector<DynamicSparsityPattern::size_type> &new_indices)
    {
      DynamicSparsityPattern::size_type starting_point =
        numbers::invalid_size_type;
      DynamicSparsityPattern::size_type min_coordination = sparsity.n_rows();
      for (DynamicSparsityPattern::size_type row = 0; row < sparsity.n_rows();
           ++row)
        // look over all as-yet unnumbered indices
        if (new_indices[row] == numbers::invalid_size_type)
          {
            if (sparsity.row_length(row) < min_coordination)
              {
                min_coordination = sparsity.row_length(row);
                starting_point   = row;
              }
          }

      // now we still have to care for the case that no unnumbered dof has a
      // coordination number less than sparsity.n_rows(). this rather exotic
      // case only happens if we only have one cell, as far as I can see,
      // but there may be others as well.
      //
      // if that should be the case, we can chose an arbitrary dof as
      // starting point, e.g. the first unnumbered one
      if (starting_point == numbers::invalid_size_type)
        {
          for (DynamicSparsityPattern::size_type i = 0; i < new_indices.size();
               ++i)
            if (new_indices[i] == numbers::invalid_size_type)
              {
                starting_point = i;
                break;
              }

          Assert(starting_point != numbers::invalid_size_type,
                 ExcInternalError());
        }

      return starting_point;
    }
  } // namespace internal



  void
  reorder_Cuthill_McKee(
    const DynamicSparsityPattern                         &sparsity,
    std::vector<DynamicSparsityPattern::size_type>       &new_indices,
    const std::vector<DynamicSparsityPattern::size_type> &starting_indices)
  {
    Assert(sparsity.n_rows() == sparsity.n_cols(),
           ExcDimensionMismatch(sparsity.n_rows(), sparsity.n_cols()));
    Assert(sparsity.n_rows() == new_indices.size(),
           ExcDimensionMismatch(sparsity.n_rows(), new_indices.size()));
    Assert(starting_indices.size() <= sparsity.n_rows(),
           ExcMessage(
             "You can't specify more starting indices than there are rows"));
    Assert(sparsity.row_index_set().size() == 0 ||
             sparsity.row_index_set().size() == sparsity.n_rows(),
           ExcMessage(
             "Only valid for sparsity patterns which store all rows."));
    for (const auto starting_index : starting_indices)
      {
        (void)starting_index;
        Assert(starting_index < sparsity.n_rows(),
               ExcMessage("Invalid starting index: All starting indices need "
                          "to be between zero and the number of rows in the "
                          "sparsity pattern."));
      }

    // store the indices of the dofs renumbered in the last round. Default to
    // starting points
    std::vector<DynamicSparsityPattern::size_type> last_round_dofs(
      starting_indices);

    // initialize the new_indices array with invalid values
    std::fill(new_indices.begin(),
              new_indices.end(),
              numbers::invalid_size_type);

    // if no starting indices were given: find dof with lowest coordination
    // number
    if (last_round_dofs.empty())
      last_round_dofs.push_back(
        internal::find_unnumbered_starting_index(sparsity, new_indices));

    // store next free dof index
    DynamicSparsityPattern::size_type next_free_number = 0;

    // enumerate the first round dofs
    for (const auto &last_round_dof : last_round_dofs)
      new_indices[last_round_dof] = next_free_number++;

    // store the indices of the dofs to be renumbered in the next round
    std::vector<DynamicSparsityPattern::size_type> next_round_dofs;

    // store for each coordination number the dofs with these coordination
    // number
    std::vector<std::pair<unsigned int, DynamicSparsityPattern::size_type>>
      dofs_by_coordination;

    // now do as many steps as needed to renumber all dofs
    while (true)
      {
        next_round_dofs.clear();

        // find all neighbors of the dofs numbered in the last round
        for (const auto dof : last_round_dofs)
          {
            const unsigned int row_length = sparsity.row_length(dof);
            for (unsigned int i = 0; i < row_length; ++i)
              {
                // skip dofs which are already numbered
                const auto column = sparsity.column_number(dof, i);
                if (new_indices[column] == numbers::invalid_size_type)
                  {
                    next_round_dofs.push_back(column);

                    // assign a dummy value to 'new_indices' to avoid adding
                    // the same index again; those will get the right number
                    // at the end of the outer 'while' loop
                    new_indices[column] = 0;
                  }
              }
          }

        // check whether there are any new dofs in the list. if there are
        // none, then we have completely numbered the current component of the
        // graph. check if there are as yet unnumbered components of the graph
        // that we would then have to do next
        if (next_round_dofs.empty())
          {
            if (std::find(new_indices.begin(),
                          new_indices.end(),
                          numbers::invalid_size_type) == new_indices.end())
              // no unnumbered indices, so we can leave now
              break;

            // otherwise find a valid starting point for the next component of
            // the graph and continue with numbering that one. we only do so
            // if no starting indices were provided by the user (see the
            // documentation of this function) so produce an error if we got
            // here and starting indices were given
            Assert(starting_indices.empty(),
                   ExcMessage("The input graph appears to have more than one "
                              "component, but as stated in the documentation "
                              "we only want to reorder such graphs if no "
                              "starting indices are given. The function was "
                              "called with starting indices, however."));

            next_round_dofs.push_back(
              internal::find_unnumbered_starting_index(sparsity, new_indices));
          }


        // find coordination number for each of these dofs
        dofs_by_coordination.clear();
        for (const types::global_dof_index next_round_dof : next_round_dofs)
          dofs_by_coordination.emplace_back(sparsity.row_length(next_round_dof),
                                            next_round_dof);
        std::sort(dofs_by_coordination.begin(), dofs_by_coordination.end());

        // assign new DoF numbers to the elements of the present front:
        for (const auto &i : dofs_by_coordination)
          new_indices[i.second] = next_free_number++;

        // after that: use this round's dofs for the next round
        last_round_dofs.swap(next_round_dofs);
      }

    // test for all indices numbered. this mostly tests whether the
    // front-marching-algorithm (which Cuthill-McKee actually is) has reached
    // all points.
    Assert((std::find(new_indices.begin(),
                      new_indices.end(),
                      numbers::invalid_size_type) == new_indices.end()) &&
             (next_free_number == sparsity.n_rows()),
           ExcInternalError());
  }



  namespace internal
  {
    void
    reorder_hierarchical(
      const DynamicSparsityPattern                   &connectivity,
      std::vector<DynamicSparsityPattern::size_type> &renumbering)
    {
      AssertDimension(connectivity.n_rows(), connectivity.n_cols());
      AssertDimension(connectivity.n_rows(), renumbering.size());
      Assert(connectivity.row_index_set().size() == 0 ||
               connectivity.row_index_set().size() == connectivity.n_rows(),
             ExcMessage(
               "Only valid for sparsity patterns which store all rows."));

      // The algorithm below works by partitioning the rows in the
      // connectivity graph, called nodes, into groups. The groups are defined
      // as those nodes in immediate neighborhood of some pivot node, which we
      // choose by minimal adjacency below.

      // We define two types of node categories for nodes not yet classified,
      // one consisting of all nodes we've not seen at all, and one for nodes
      // identified as neighbors (variable current_neighbors below) but not
      // yet grouped. We use this classification in combination with an
      // unsorted vector, which is much faster than keeping a sorted data
      // structure (e.g. std::set)
      constexpr types::global_dof_index unseen_node =
        numbers::invalid_dof_index;
      constexpr types::global_dof_index    available_node = unseen_node - 1;
      const types::global_dof_index        n_nodes = connectivity.n_rows();
      std::vector<types::global_dof_index> touched_nodes(n_nodes, unseen_node);

      std::vector<unsigned int>            row_lengths(n_nodes);
      std::vector<types::global_dof_index> current_neighbors;
      std::vector<types::global_dof_index> group_starts(1);
      std::vector<types::global_dof_index> group_indices;
      group_indices.reserve(n_nodes);

      // First collect the number of neighbors for each node. We use this
      // field to find the next node with the minimum number of non-touched
      // neighbors in the field n_remaining_neighbors, so we will count down
      // on this field. We also cache the row lengths because we need this
      // data frequently and getting it from the sparsity pattern is more
      // expensive.
      for (types::global_dof_index row = 0; row < n_nodes; ++row)
        {
          row_lengths[row] = connectivity.row_length(row);
          Assert(row_lengths[row] > 0, ExcInternalError());
        }
      std::vector<unsigned int> n_remaining_neighbors(row_lengths);

      // This outer loop is typically traversed only once, unless the global
      // graph is not connected
      while (true)
        {
          // Find node with the minimal number of neighbors (typically a
          // corner node when based on FEM meshes). If no node is left, we are
          // done. Together with the outer while loop, this loop can possibly
          // be of quadratic complexity in the number of disconnected
          // partitions, i.e. up to n_nodes in the worst case,
          // but that is not the usual use case of this loop and thus not
          // optimized for.
          {
            unsigned int candidate_valence = numbers::invalid_unsigned_int;
            types::global_dof_index candidate_index =
              numbers::invalid_dof_index;
            for (types::global_dof_index i = 0; i < n_nodes; ++i)
              if (touched_nodes[i] == unseen_node)
                if (row_lengths[i] < candidate_valence)
                  {
                    candidate_index   = i;
                    candidate_valence = n_remaining_neighbors[i];
                    if (candidate_valence <= 1)
                      break;
                  }
            if (candidate_index == numbers::invalid_dof_index)
              break;

            Assert(candidate_valence > 0, ExcInternalError());

            current_neighbors              = {candidate_index};
            touched_nodes[candidate_index] = available_node;
          }

          while (true)
            {
              // Find node with minimum number of untouched neighbors among
              // the next set of possible neighbors (= valence), and among the
              // set of nodes with the minimal number of neighbors, choose the
              // one with the largest number of touched neighbors (i.e., the
              // largest row length).
              //
              // This loop is typically the most expensive part for large
              // graphs and thus only run once. We also do some cleanup, i.e.,
              // the indices added to a group in the previous round need to be
              // removed at this point.
              unsigned int candidate_valence = numbers::invalid_unsigned_int;
              types::global_dof_index candidate_index =
                numbers::invalid_dof_index;
              unsigned int       candidate_row_length = 0;
              const unsigned int loop_length = current_neighbors.size();
              unsigned int       write_index = 0;
              for (unsigned int i = 0; i < loop_length; ++i)
                {
                  const types::global_dof_index node = current_neighbors[i];
                  Assert(touched_nodes[node] != unseen_node,
                         ExcInternalError());
                  if (touched_nodes[node] == available_node)
                    {
                      current_neighbors[write_index] = node;
                      ++write_index;
                      if (n_remaining_neighbors[node] < candidate_valence ||
                          (n_remaining_neighbors[node] == candidate_valence &&
                           (row_lengths[node] > candidate_row_length ||
                            (row_lengths[node] == candidate_row_length &&
                             node < candidate_index))))
                        {
                          candidate_index      = node;
                          candidate_valence    = n_remaining_neighbors[node];
                          candidate_row_length = row_lengths[node];
                        }
                    }
                }
              current_neighbors.resize(write_index);

              if constexpr (running_in_debug_mode())
                {
                  for (const types::global_dof_index node : current_neighbors)
                    Assert(touched_nodes[node] == available_node,
                           ExcInternalError());
                }

              // No more neighbors left -> terminate loop
              if (current_neighbors.empty())
                break;

              // Add the pivot and all direct neighbors of the pivot node not
              // yet touched to the list of new entries.
              group_indices.push_back(candidate_index);
              touched_nodes[candidate_index] = group_starts.size() - 1;
              const auto end_it = connectivity.end(candidate_index);
              for (auto it = connectivity.begin(candidate_index); it != end_it;
                   ++it)
                if (touched_nodes[it->column()] >= available_node)
                  {
                    group_indices.push_back(it->column());
                    touched_nodes[it->column()] = group_starts.size() - 1;
                  }
              group_starts.push_back(group_indices.size());

              // Add all neighbors of the current list not yet seen to the set
              // of possible next nodes. The added node is grouped and thus no
              // longer a valid neighbor (here we assume symmetry of the
              // connectivity). It will be removed from the list of neighbors
              // by the code further up in the next iteration of the
              // surrounding loop.
              for (types::global_dof_index index =
                     group_starts[group_starts.size() - 2];
                   index < group_starts.back();
                   ++index)
                {
                  auto       it      = connectivity.begin(group_indices[index]);
                  const auto end_row = connectivity.end(group_indices[index]);
                  for (; it != end_row; ++it)
                    {
                      if (touched_nodes[it->column()] == unseen_node)
                        {
                          current_neighbors.push_back(it->column());
                          touched_nodes[it->column()] = available_node;
                        }
                      n_remaining_neighbors[it->column()]--;
                    }
                }
            }
        }

      // Sanity check: for all nodes, there should not be any neighbors left
      for (types::global_dof_index row = 0; row < n_nodes; ++row)
        Assert(n_remaining_neighbors[row] == 0, ExcInternalError());

      // If the number of groups is smaller than the number of nodes, we
      // continue by recursively calling this method
      const unsigned int n_groups = group_starts.size() - 1;
      if (n_groups < n_nodes)
        {
          // Form the connectivity of the groups
          DynamicSparsityPattern connectivity_next(n_groups, n_groups);
          for (types::global_dof_index row = 0; row < n_groups; ++row)
            for (types::global_dof_index index = group_starts[row];
                 index < group_starts[row + 1];
                 ++index)
              {
                auto       it     = connectivity.begin(group_indices[index]);
                const auto end_it = connectivity.end(group_indices[index]);
                for (; it != end_it; ++it)
                  connectivity_next.add(row, touched_nodes[it->column()]);
              }

          // Recursively call the reordering
          std::vector<types::global_dof_index> renumbering_next(n_groups);
          reorder_hierarchical(connectivity_next, renumbering_next);

          // Renumber the indices group by group according to the incoming
          // ordering for the groups
          for (types::global_dof_index row = 0, c = 0; row < n_groups; ++row)
            for (types::global_dof_index index =
                   group_starts[renumbering_next[row]];
                 index < group_starts[renumbering_next[row] + 1];
                 ++index, ++c)
              renumbering[c] = group_indices[index];
        }
      else
        {
          // All groups should have size one and no more recursion is possible,
          // so use the numbering of the groups
          unsigned int c = 0;
          for (const types::global_dof_index i : group_indices)
            renumbering[c++] = i;
        }
    }
  } // namespace internal

  void
  reorder_hierarchical(
    const DynamicSparsityPattern                   &connectivity,
    std::vector<DynamicSparsityPattern::size_type> &renumbering)
  {
    // the internal renumbering keeps the numbering the wrong way around (but
    // we cannot invert the numbering inside that method because it is used
    // recursively), so invert it here
    internal::reorder_hierarchical(connectivity, renumbering);
    renumbering = Utilities::invert_permutation(renumbering);
  }



#ifdef DEAL_II_WITH_MPI

  void
  gather_sparsity_pattern(DynamicSparsityPattern &dsp,
                          const IndexSet         &locally_owned_rows,
                          const MPI_Comm          mpi_comm,
                          const IndexSet         &locally_relevant_rows)
  {
    using map_vec_t =
      std::map<unsigned int, std::vector<DynamicSparsityPattern::size_type>>;

    // 1. limit rows to non owned:
    IndexSet requested_rows(locally_relevant_rows);
    requested_rows.subtract_set(locally_owned_rows);

    std::vector<unsigned int> index_owner =
      Utilities::MPI::compute_index_owner(locally_owned_rows,
                                          requested_rows,
                                          mpi_comm);

    // 2. go through requested_rows, figure out the owner and add the row to
    // request
    map_vec_t rows_data;
    for (DynamicSparsityPattern::size_type i = 0;
         i < requested_rows.n_elements();
         ++i)
      {
        const DynamicSparsityPattern::size_type row =
          requested_rows.nth_index_in_set(i);

        rows_data[index_owner[i]].push_back(row);
      }

    // 3. get what others ask us to send
    const auto rows_data_received =
      Utilities::MPI::some_to_some(mpi_comm, rows_data);

    // 4. now prepare data to be sent in the same format as in
    // distribute_sparsity_pattern() below, i.e.
    // rX,num_rX,cols_rX
    map_vec_t send_data;
    for (const auto &data : rows_data_received)
      {
        for (const auto &row : data.second)
          {
            const auto rlen = dsp.row_length(row);

            // skip empty lines
            if (rlen == 0)
              continue;

            // save entries
            send_data[data.first].push_back(row);  // row index
            send_data[data.first].push_back(rlen); // number of entries
            for (DynamicSparsityPattern::size_type c = 0; c < rlen; ++c)
              send_data[data.first].push_back(
                dsp.column_number(row, c)); // columns
          }                                 // loop over rows
      }                                     // loop over received data

    // 5. communicate rows
    const auto received_data =
      Utilities::MPI::some_to_some(mpi_comm, send_data);

    // 6. add result to our sparsity
    for (const auto &data : received_data)
      {
        const auto &recv_buf = data.second;
        auto        ptr      = recv_buf.begin();
        const auto  end      = recv_buf.end();
        while (ptr != end)
          {
            const auto row = *(ptr++);
            Assert(ptr != end, ExcInternalError());

            const auto n_entries = *(ptr++);
            Assert(n_entries > 0, ExcInternalError());
            Assert(ptr != end, ExcInternalError());

            // make sure we clear whatever was previously stored
            // in these rows. Otherwise we can't guarantee that the
            // data is consistent across MPI communicator.
            dsp.clear_row(row);

            Assert(ptr + (n_entries - 1) != end, ExcInternalError());
            dsp.add_entries(row, ptr, ptr + n_entries, true);
            ptr += n_entries;
          }
        Assert(ptr == end, ExcInternalError());
      }
  }



  void
  distribute_sparsity_pattern(
    DynamicSparsityPattern                               &dsp,
    const std::vector<DynamicSparsityPattern::size_type> &rows_per_cpu,
    const MPI_Comm                                        mpi_comm,
    const IndexSet                                       &myrange)
  {
    const unsigned int myid = Utilities::MPI::this_mpi_process(mpi_comm);
    std::vector<DynamicSparsityPattern::size_type> start_index(
      rows_per_cpu.size() + 1);
    start_index[0] = 0;
    for (DynamicSparsityPattern::size_type i = 0; i < rows_per_cpu.size(); ++i)
      start_index[i + 1] = start_index[i] + rows_per_cpu[i];

    IndexSet owned(start_index.back());
    owned.add_range(start_index[myid], start_index[myid] + rows_per_cpu[myid]);

    distribute_sparsity_pattern(dsp, owned, mpi_comm, myrange);
  }



  void
  distribute_sparsity_pattern(DynamicSparsityPattern &dsp,
                              const IndexSet         &locally_owned_rows,
                              const MPI_Comm          mpi_comm,
                              const IndexSet         &locally_relevant_rows)
  {
    AssertThrow(
      dsp.row_index_set() == locally_relevant_rows,
      ExcMessage(
        "The DynamicSparsityPattern must be initialized with an IndexSet that contains locally relevant indices."));

    IndexSet requested_rows(locally_relevant_rows);
    requested_rows.subtract_set(locally_owned_rows);

    std::vector<unsigned int> index_owner =
      Utilities::MPI::compute_index_owner(locally_owned_rows,
                                          requested_rows,
                                          mpi_comm);

    using map_vec_t =
      std::map<unsigned int, std::vector<DynamicSparsityPattern::size_type>>;

    map_vec_t send_data;

    for (DynamicSparsityPattern::size_type i = 0;
         i < requested_rows.n_elements();
         ++i)
      {
        const DynamicSparsityPattern::size_type row =
          requested_rows.nth_index_in_set(i);

        const auto rlen = dsp.row_length(row);

        // skip empty lines
        if (rlen == 0)
          continue;

        // save entries
        send_data[index_owner[i]].push_back(row);  // row index
        send_data[index_owner[i]].push_back(rlen); // number of entries
        for (DynamicSparsityPattern::size_type c = 0; c < rlen; ++c)
          {
            // columns
            const auto column = dsp.column_number(row, c);
            send_data[index_owner[i]].push_back(column);
          }
      }

    const auto receive_data = Utilities::MPI::some_to_some(mpi_comm, send_data);

    // add what we received
    for (const auto &data : receive_data)
      {
        const auto &recv_buf = data.second;
        auto        ptr      = recv_buf.begin();
        const auto  end      = recv_buf.end();
        while (ptr != end)
          {
            const auto row = *(ptr++);
            Assert(ptr != end, ExcInternalError());
            const auto n_entries = *(ptr++);

            Assert(ptr + (n_entries - 1) != end, ExcInternalError());
            dsp.add_entries(row, ptr, ptr + n_entries, true);
            ptr += n_entries;
          }
        Assert(ptr == end, ExcInternalError());
      }
  }



  void
  distribute_sparsity_pattern(BlockDynamicSparsityPattern &dsp,
                              const std::vector<IndexSet> &owned_set_per_cpu,
                              const MPI_Comm               mpi_comm,
                              const IndexSet              &myrange)
  {
    const unsigned int myid = Utilities::MPI::this_mpi_process(mpi_comm);
    distribute_sparsity_pattern(dsp,
                                owned_set_per_cpu[myid],
                                mpi_comm,
                                myrange);
  }



  void
  distribute_sparsity_pattern(BlockDynamicSparsityPattern &dsp,
                              const IndexSet              &locally_owned_rows,
                              const MPI_Comm               mpi_comm,
                              const IndexSet &locally_relevant_rows)
  {
    using map_vec_t =
      std::map<BlockDynamicSparsityPattern::size_type,
               std::vector<BlockDynamicSparsityPattern::size_type>>;
    map_vec_t send_data;

    IndexSet requested_rows(locally_relevant_rows);
    requested_rows.subtract_set(locally_owned_rows);

    std::vector<unsigned int> index_owner =
      Utilities::MPI::compute_index_owner(locally_owned_rows,
                                          requested_rows,
                                          mpi_comm);

    for (DynamicSparsityPattern::size_type i = 0;
         i < requested_rows.n_elements();
         ++i)
      {
        const DynamicSparsityPattern::size_type row =
          requested_rows.nth_index_in_set(i);

        BlockDynamicSparsityPattern::size_type rlen = dsp.row_length(row);

        // skip empty lines
        if (rlen == 0)
          continue;

        // save entries
        std::vector<BlockDynamicSparsityPattern::size_type> &dst =
          send_data[index_owner[i]];

        dst.push_back(rlen); // number of entries
        dst.push_back(row);  // row index
        for (BlockDynamicSparsityPattern::size_type c = 0; c < rlen; ++c)
          {
            // columns
            BlockDynamicSparsityPattern::size_type column =
              dsp.column_number(row, c);
            dst.push_back(column);
          }
      }

    unsigned int num_receive = 0;
    {
      std::vector<unsigned int> send_to;
      send_to.reserve(send_data.size());
      for (const auto &sparsity_line : send_data)
        send_to.push_back(sparsity_line.first);

      num_receive =
        Utilities::MPI::compute_n_point_to_point_communications(mpi_comm,
                                                                send_to);
    }

    std::vector<MPI_Request> requests(send_data.size());


    // send data

    static Utilities::MPI::CollectiveMutex      mutex;
    Utilities::MPI::CollectiveMutex::ScopedLock lock(mutex, mpi_comm);

    const int mpi_tag = Utilities::MPI::internal::Tags::
      sparsity_tools_distribute_sparsity_pattern;

    {
      unsigned int idx = 0;
      for (const auto &sparsity_line : send_data)
        {
          const int ierr = MPI_Isend(
            sparsity_line.second.data(),
            sparsity_line.second.size(),
            Utilities::MPI::mpi_type_id_for_type<types::global_dof_index>,
            sparsity_line.first,
            mpi_tag,
            mpi_comm,
            &requests[idx++]);
          AssertThrowMPI(ierr);
        }
    }

    {
      // receive
      std::vector<BlockDynamicSparsityPattern::size_type> recv_buf;
      for (unsigned int index = 0; index < num_receive; ++index)
        {
          MPI_Status status;
          int ierr = MPI_Probe(MPI_ANY_SOURCE, mpi_tag, mpi_comm, &status);
          AssertThrowMPI(ierr);

          int len;
          ierr = MPI_Get_count(
            &status,
            Utilities::MPI::mpi_type_id_for_type<types::global_dof_index>,
            &len);
          AssertThrowMPI(ierr);

          recv_buf.resize(len);
          ierr = MPI_Recv(
            recv_buf.data(),
            len,
            Utilities::MPI::mpi_type_id_for_type<types::global_dof_index>,
            status.MPI_SOURCE,
            status.MPI_TAG,
            mpi_comm,
            &status);
          AssertThrowMPI(ierr);

          std::vector<BlockDynamicSparsityPattern::size_type>::const_iterator
            ptr = recv_buf.begin();
          std::vector<BlockDynamicSparsityPattern::size_type>::const_iterator
            end = recv_buf.end();
          while (ptr != end)
            {
              BlockDynamicSparsityPattern::size_type num = *(ptr++);
              Assert(ptr != end, ExcInternalError());
              BlockDynamicSparsityPattern::size_type row = *(ptr++);
              for (unsigned int c = 0; c < num; ++c)
                {
                  Assert(ptr != end, ExcInternalError());
                  dsp.add(row, *ptr);
                  ++ptr;
                }
            }
          Assert(ptr == end, ExcInternalError());
        }
    }

    // complete all sends, so that we can safely destroy the buffers.
    if (requests.size() > 0)
      {
        const int ierr =
          MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
        AssertThrowMPI(ierr);
      }
  }
#endif
} // namespace SparsityTools

DEAL_II_NAMESPACE_CLOSE
