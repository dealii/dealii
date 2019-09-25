// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#ifndef dealii_face_setup_internal_h
#define dealii_face_setup_internal_h

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/matrix_free/face_info.h>
#include <deal.II/matrix_free/task_info.h>

#include <fstream>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * A struct that is used to represent a collection of faces of a process
     * with one of its neighbor within the setup done in struct FaceInfo.
     *
     * @author Katharina Kormann, Martin Kronbichler, 2018
     */
    struct FaceIdentifier
    {
      FaceIdentifier()
        : n_hanging_faces_smaller_subdomain(0)
        , n_hanging_faces_larger_subdomain(0)
      {}

      std::vector<std::pair<CellId, CellId>> shared_faces;
      unsigned int                           n_hanging_faces_smaller_subdomain;
      unsigned int                           n_hanging_faces_larger_subdomain;
    };



    /**
     * A struct that extracts the faces relevant to a given set of cells,
     * including the assignment of which of the two neighboring processors at
     * a subdomain boundary with MPI should do the integration (from both
     * sides). This data structure is used for the setup of the connectivity
     * between faces and cells and for identification of the dof indices to be
     * used for face integrals.
     *
     * @author Katharina Kormann, Martin Kronbichler, 2018
     */
    template <int dim>
    struct FaceSetup
    {
      FaceSetup();

      /**
       * Perform the initial detection of faces before reading the indices on
       * the cells. This does not add the faces yet but only decides on
       * whether some of the faces should be considered for processing
       * locally.
       */
      template <typename MFAddData>
      void
      initialize(
        const dealii::Triangulation<dim> &                  triangulation,
        const MFAddData &                                   additional_data,
        std::vector<std::pair<unsigned int, unsigned int>> &cell_levels);

      /**
       * Upon completion of the dof indices, this function extracts the
       * information relevant for FaceToCellTopology and categorizes the faces
       * into interior faces, boundary faces, and ghost faces (not processed
       * locally but adjacent to some of the cells present locally).
       */
      void
      generate_faces(
        const dealii::Triangulation<dim> &                        triangulation,
        const std::vector<std::pair<unsigned int, unsigned int>> &cell_levels,
        TaskInfo &                                                task_info);

      /**
       * Fills the information about the cell, the face number, and numbers
       * within the plain array representation in MatrixFree into
       * FaceToCellTopology (without vectorization, which is something applied
       * later).
       */
      FaceToCellTopology<1>
      create_face(
        const unsigned int                                        face_no,
        const typename dealii::Triangulation<dim>::cell_iterator &cell,
        const unsigned int number_cell_interior,
        const typename dealii::Triangulation<dim>::cell_iterator &neighbor,
        const unsigned int number_cell_exterior);

      bool use_active_cells;

      /**
       * A type that categorizes faces in the first initialize() function such
       * that we can later get their correct value in generate_faces().
       */
      enum class FaceCategory : char
      {
        locally_active_at_boundary,
        locally_active_done_here,
        locally_active_done_elsewhere,
        ghosted,
        multigrid_refinement_edge
      };

      std::vector<FaceCategory>          face_is_owned;
      std::vector<bool>                  at_processor_boundary;
      std::vector<unsigned int>          cells_close_to_boundary;
      std::vector<FaceToCellTopology<1>> inner_faces;
      std::vector<FaceToCellTopology<1>> boundary_faces;
      std::vector<FaceToCellTopology<1>> inner_ghost_faces;
      std::vector<FaceToCellTopology<1>> refinement_edge_faces;
    };



    /**
     * Actually form the batches for vectorized execution of face integrals.
     */
    template <int vectorization_width>
    void
    collect_faces_vectorization(
      const std::vector<FaceToCellTopology<1>> &faces_in,
      const std::vector<bool> &                 hard_vectorization_boundary,
      std::vector<unsigned int> &               face_partition_data,
      std::vector<FaceToCellTopology<vectorization_width>> &faces_out);



    /* -------------------------------------------------------------------- */

#ifndef DOXYGEN

    template <int dim>
    FaceSetup<dim>::FaceSetup()
      : use_active_cells(true)
    {}



    template <int dim>
    template <typename MFAddData>
    void
    FaceSetup<dim>::initialize(
      const dealii::Triangulation<dim> &                  triangulation,
      const MFAddData &                                   additional_data,
      std::vector<std::pair<unsigned int, unsigned int>> &cell_levels)
    {
      use_active_cells =
        additional_data.level_mg_handler == numbers::invalid_unsigned_int;

#  ifdef DEBUG
      // safety check
      if (use_active_cells)
        for (const auto &cell_level : cell_levels)
          {
            typename dealii::Triangulation<dim>::cell_iterator dcell(
              &triangulation, cell_level.first, cell_level.second);
            Assert(dcell->active(), ExcInternalError());
          }
#  endif

      // step 1: add ghost cells for those cells that we identify as
      // interesting

      at_processor_boundary.resize(cell_levels.size(), false);
      cells_close_to_boundary.clear();
      face_is_owned.resize(dim > 1 ? triangulation.n_raw_faces() :
                                     triangulation.n_vertices(),
                           FaceCategory::locally_active_done_elsewhere);

      // go through the mesh and divide the faces on the processor
      // boundaries as evenly as possible between the processors
      std::map<types::subdomain_id, FaceIdentifier>
        inner_faces_at_proc_boundary;
      if (triangulation.locally_owned_subdomain() !=
          numbers::invalid_subdomain_id)
        {
          const types::subdomain_id my_domain =
            triangulation.locally_owned_subdomain();
          for (unsigned int i = 0; i < cell_levels.size(); ++i)
            {
              if (i > 0 && cell_levels[i] == cell_levels[i - 1])
                continue;
              typename dealii::Triangulation<dim>::cell_iterator dcell(
                &triangulation, cell_levels[i].first, cell_levels[i].second);
              for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell;
                   ++f)
                {
                  if (dcell->at_boundary(f) && !dcell->has_periodic_neighbor(f))
                    continue;
                  typename dealii::Triangulation<dim>::cell_iterator neighbor =
                    dcell->neighbor_or_periodic_neighbor(f);

                  // faces at hanging nodes are always treated by the processor
                  // who owns the element on the fine side. but we need to count
                  // the number of inner faces in order to balance the remaining
                  // faces properly
                  const CellId id_mine = dcell->id();
                  if (use_active_cells && neighbor->has_children())
                    for (unsigned int c = 0; c < dcell->face(f)->n_children();
                         ++c)
                      {
                        typename dealii::Triangulation<dim>::cell_iterator
                          neighbor_c =
                            dcell->at_boundary(f) ?
                              dcell->periodic_neighbor_child_on_subface(f, c) :
                              dcell->neighbor_child_on_subface(f, c);
                        const types::subdomain_id neigh_domain =
                          neighbor_c->subdomain_id();
                        if (my_domain < neigh_domain)
                          inner_faces_at_proc_boundary[neigh_domain]
                            .n_hanging_faces_larger_subdomain++;
                        else if (my_domain > neigh_domain)
                          inner_faces_at_proc_boundary[neigh_domain]
                            .n_hanging_faces_smaller_subdomain++;
                      }
                  else
                    {
                      const types::subdomain_id neigh_domain =
                        use_active_cells ? neighbor->subdomain_id() :
                                           neighbor->level_subdomain_id();
                      if (neighbor->level() < dcell->level() &&
                          use_active_cells)
                        {
                          if (my_domain < neigh_domain)
                            inner_faces_at_proc_boundary[neigh_domain]
                              .n_hanging_faces_smaller_subdomain++;
                          else if (my_domain > neigh_domain)
                            inner_faces_at_proc_boundary[neigh_domain]
                              .n_hanging_faces_larger_subdomain++;
                        }
                      else if (neighbor->level() == dcell->level() &&
                               my_domain != neigh_domain)
                        {
                          // always list the cell whose owner has the lower
                          // subdomain id first. this applies to both processors
                          // involved, so both processors will generate the same
                          // list that we will later order
                          const CellId id_neigh = neighbor->id();
                          if (my_domain < neigh_domain)
                            inner_faces_at_proc_boundary[neigh_domain]
                              .shared_faces.emplace_back(id_mine, id_neigh);
                          else
                            inner_faces_at_proc_boundary[neigh_domain]
                              .shared_faces.emplace_back(id_neigh, id_mine);
                        }
                    }
                }
            }

          // sort the cell ids related to each neighboring processor. This
          // algorithm is symmetric so every processor combination should
          // arrive here and no deadlock should be possible
          for (auto &inner_face : inner_faces_at_proc_boundary)
            {
              Assert(inner_face.first != my_domain,
                     ExcInternalError("Should not send info to myself"));
              std::sort(inner_face.second.shared_faces.begin(),
                        inner_face.second.shared_faces.end());
              inner_face.second.shared_faces.erase(
                std::unique(inner_face.second.shared_faces.begin(),
                            inner_face.second.shared_faces.end()),
                inner_face.second.shared_faces.end());

              // safety check: both involved processors should see the same list
              // because the pattern of ghosting is symmetric. We test this by
              // looking at the length of the lists of faces
#  if defined(DEAL_II_WITH_MPI) && defined(DEBUG)
              MPI_Comm comm = MPI_COMM_SELF;
              if (const parallel::TriangulationBase<dim> *ptria =
                    dynamic_cast<const parallel::TriangulationBase<dim> *>(
                      &triangulation))
                comm = ptria->get_communicator();

              MPI_Status   status;
              unsigned int mysize    = inner_face.second.shared_faces.size();
              unsigned int othersize = numbers::invalid_unsigned_int;
              MPI_Sendrecv(&mysize,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           600 + my_domain,
                           &othersize,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           600 + inner_face.first,
                           comm,
                           &status);
              AssertDimension(mysize, othersize);
              mysize = inner_face.second.n_hanging_faces_smaller_subdomain;
              MPI_Sendrecv(&mysize,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           700 + my_domain,
                           &othersize,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           700 + inner_face.first,
                           comm,
                           &status);
              AssertDimension(mysize, othersize);
              mysize = inner_face.second.n_hanging_faces_larger_subdomain;
              MPI_Sendrecv(&mysize,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           800 + my_domain,
                           &othersize,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           800 + inner_face.first,
                           comm,
                           &status);
              AssertDimension(mysize, othersize);
#  endif

              // Arrange the face "ownership" such that cells that are access
              // by more than one face (think of a cell in a corner) get
              // ghosted. This arrangement has the advantage that we need to
              // send less data because the same data is used twice. The
              // strategy applied here is to ensure the same order of face
              // pairs on both processors that share some faces, and make the
              // same decision on both sides.

              // Create a vector with cell ids sorted over the processor with
              // the larger rank. In the code below we need to be able to
              // identify the same cell once for the processor with higher
              // rank and once for the processor with the lower rank. The
              // format for the processor with the higher rank is already
              // contained in `shared_faces`, whereas we need a copy that we
              // sort differently for the other way around.
              std::vector<std::tuple<CellId, CellId, unsigned int>> other_range(
                inner_face.second.shared_faces.size());
              for (unsigned int i = 0; i < other_range.size(); ++i)
                other_range[i] =
                  std::make_tuple(inner_face.second.shared_faces[i].second,
                                  inner_face.second.shared_faces[i].first,
                                  i);
              std::sort(other_range.begin(), other_range.end());

              // the vector 'assignment' sets whether a particular cell
              // appears more often and acts as a pre-selection of the rank. A
              // value of 1 means that the process with the higher rank gets
              // those faces, a value -1 means that the process with the lower
              // rank gets it, whereas a value 0 means that the decision can
              // be made in an arbitrary way.
              unsigned int      n_faces_lower_proc = 0, n_faces_higher_proc = 0;
              std::vector<char> assignment(other_range.size(), 0);
              if (inner_face.second.shared_faces.size() > 0)
                {
                  // identify faces that go to the processor with the higher
                  // rank
                  unsigned int count = 0;
                  for (unsigned int i = 1;
                       i < inner_face.second.shared_faces.size();
                       ++i)
                    if (inner_face.second.shared_faces[i].first ==
                        inner_face.second.shared_faces[i - 1 - count].first)
                      ++count;
                    else
                      {
                        AssertThrow(count < 2 * dim, ExcInternalError());
                        if (count > 0)
                          {
                            for (unsigned int k = 0; k <= count; ++k)
                              assignment[i - 1 - k] = 1;
                            n_faces_higher_proc += count + 1;
                          }
                        count = 0;
                      }

                  // identify faces that definitely go to the processor with
                  // the lower rank - this must use the sorting of CellId
                  // variables from the processor with the higher rank, i.e.,
                  // other_range rather than `shared_faces`.
                  count = 0;
                  for (unsigned int i = 1; i < other_range.size(); ++i)
                    if (std::get<0>(other_range[i]) ==
                        std::get<0>(other_range[i - 1 - count]))
                      ++count;
                    else
                      {
                        AssertThrow(count < 2 * dim, ExcInternalError());
                        if (count > 0)
                          {
                            for (unsigned int k = 0; k <= count; ++k)
                              {
                                Assert(inner_face.second
                                           .shared_faces[std::get<2>(
                                             other_range[i - 1])]
                                           .second ==
                                         inner_face.second
                                           .shared_faces[std::get<2>(
                                             other_range[i - 1 - k])]
                                           .second,
                                       ExcInternalError());
                                // only assign to -1 if higher rank was not
                                // yet set
                                if (assignment[std::get<2>(
                                      other_range[i - 1 - k])] == 0)
                                  {
                                    assignment[std::get<2>(
                                      other_range[i - 1 - k])] = -1;
                                    ++n_faces_lower_proc;
                                  }
                              }
                          }
                        count = 0;
                      }
                }


              // divide the faces evenly between the two processors. the
              // processor with small rank takes the first half, the processor
              // with larger rank the second half. Adjust for the hanging
              // faces that always get assigned to one side, and the faces we
              // have already assigned due to the criterion above
              n_faces_lower_proc +=
                inner_face.second.n_hanging_faces_smaller_subdomain;
              n_faces_higher_proc +=
                inner_face.second.n_hanging_faces_larger_subdomain;
              const unsigned int n_total_faces_at_proc_boundary =
                (inner_face.second.shared_faces.size() +
                 inner_face.second.n_hanging_faces_smaller_subdomain +
                 inner_face.second.n_hanging_faces_larger_subdomain);
              unsigned int split_index = n_total_faces_at_proc_boundary / 2;
              if (split_index < n_faces_lower_proc)
                split_index = 0;
              else if (split_index <
                       n_total_faces_at_proc_boundary - n_faces_higher_proc)
                split_index -= n_faces_lower_proc;
              else
                split_index = n_total_faces_at_proc_boundary -
                              n_faces_higher_proc - n_faces_lower_proc;

                // make sure the splitting is consistent between both sides
#  if defined(DEAL_II_WITH_MPI) && defined(DEBUG)
              MPI_Sendrecv(&split_index,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           900 + my_domain,
                           &othersize,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           900 + inner_face.first,
                           comm,
                           &status);
              AssertDimension(split_index, othersize);
              MPI_Sendrecv(&n_faces_lower_proc,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           1000 + my_domain,
                           &othersize,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           1000 + inner_face.first,
                           comm,
                           &status);
              AssertDimension(n_faces_lower_proc, othersize);
              MPI_Sendrecv(&n_faces_higher_proc,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           1100 + my_domain,
                           &othersize,
                           1,
                           MPI_UNSIGNED,
                           inner_face.first,
                           1100 + inner_face.first,
                           comm,
                           &status);
              AssertDimension(n_faces_higher_proc, othersize);
#  endif

              // collect the faces on both sides
              std::vector<std::pair<CellId, CellId>> owned_faces_lower,
                owned_faces_higher;
              for (unsigned int i = 0; i < assignment.size(); ++i)
                if (assignment[i] < 0)
                  owned_faces_lower.push_back(
                    inner_face.second.shared_faces[i]);
                else if (assignment[i] > 0)
                  owned_faces_higher.push_back(
                    inner_face.second.shared_faces[i]);
              AssertIndexRange(split_index,
                               inner_face.second.shared_faces.size() + 1 -
                                 owned_faces_lower.size() -
                                 owned_faces_higher.size());

              unsigned int i = 0, c = 0;
              for (; i < assignment.size() && c < split_index; ++i)
                if (assignment[i] == 0)
                  {
                    owned_faces_lower.push_back(
                      inner_face.second.shared_faces[i]);
                    ++c;
                  }
              for (; i < assignment.size(); ++i)
                if (assignment[i] == 0)
                  {
                    owned_faces_higher.push_back(
                      inner_face.second.shared_faces[i]);
                  }

#  ifdef DEBUG
              // check consistency of faces on both sides
              std::vector<std::pair<CellId, CellId>> check_faces;
              check_faces.insert(check_faces.end(),
                                 owned_faces_lower.begin(),
                                 owned_faces_lower.end());
              check_faces.insert(check_faces.end(),
                                 owned_faces_higher.begin(),
                                 owned_faces_higher.end());
              std::sort(check_faces.begin(), check_faces.end());
              AssertDimension(check_faces.size(),
                              inner_face.second.shared_faces.size());
              for (unsigned int i = 0; i < check_faces.size(); ++i)
                Assert(check_faces[i] == inner_face.second.shared_faces[i],
                       ExcInternalError());
#  endif

              // now only set half of the faces as the ones to keep
              if (my_domain < inner_face.first)
                inner_face.second.shared_faces.swap(owned_faces_lower);
              else
                inner_face.second.shared_faces.swap(owned_faces_higher);

              std::sort(inner_face.second.shared_faces.begin(),
                        inner_face.second.shared_faces.end());
            }
        }

      // fill in the additional cells that we need access to via ghosting to
      // cell_levels
      std::set<std::pair<unsigned int, unsigned int>> ghost_cells;
      for (unsigned int i = 0; i < cell_levels.size(); ++i)
        {
          typename dealii::Triangulation<dim>::cell_iterator dcell(
            &triangulation, cell_levels[i].first, cell_levels[i].second);
          if (use_active_cells)
            Assert(dcell->active(), ExcNotImplemented());
          for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
            {
              if (dcell->at_boundary(f) && !dcell->has_periodic_neighbor(f))
                face_is_owned[dcell->face(f)->index()] =
                  FaceCategory::locally_active_at_boundary;

              // treat boundaries of cells of different refinement level
              // inside the domain in case of multigrid separately
              else if ((dcell->at_boundary(f) == false ||
                        dcell->has_periodic_neighbor(f)) &&
                       additional_data.level_mg_handler !=
                         numbers::invalid_unsigned_int &&
                       dcell->neighbor_or_periodic_neighbor(f)->level() <
                         dcell->level())
                {
                  face_is_owned[dcell->face(f)->index()] =
                    FaceCategory::multigrid_refinement_edge;
                }
              else
                {
                  typename dealii::Triangulation<dim>::cell_iterator neighbor =
                    dcell->neighbor_or_periodic_neighbor(f);

                  // neighbor is refined -> face will be treated by neighbor
                  if (use_active_cells && neighbor->has_children() &&
                      additional_data.hold_all_faces_to_owned_cells == false)
                    continue;

                  bool add_to_ghost = false;
                  const types::subdomain_id
                    id1 = use_active_cells ? dcell->subdomain_id() :
                                             dcell->level_subdomain_id(),
                    id2 = use_active_cells ?
                            (neighbor->has_children() ?
                               dcell->neighbor_child_on_subface(f, 0)
                                 ->subdomain_id() :
                               neighbor->subdomain_id()) :
                            neighbor->level_subdomain_id();

                  // Check whether the current face should be processed
                  // locally (instead of being processed from the other
                  // side). We process a face locally when we are more refined
                  // (in the active cell case) or when the face is listed in
                  // the `shared_faces` data structure that we built above.
                  if ((id1 == id2 && (use_active_cells == false ||
                                      neighbor->has_children() == false)) ||
                      dcell->level() > neighbor->level() ||
                      std::binary_search(
                        inner_faces_at_proc_boundary[id2].shared_faces.begin(),
                        inner_faces_at_proc_boundary[id2].shared_faces.end(),
                        std::make_pair(id1 < id2 ? dcell->id() : neighbor->id(),
                                       id1 < id2 ? neighbor->id() :
                                                   dcell->id())))
                    {
                      face_is_owned[dcell->face(f)->index()] =
                        FaceCategory::locally_active_done_here;
                      if (dcell->level() == neighbor->level())
                        face_is_owned
                          [neighbor
                             ->face(dcell->has_periodic_neighbor(f) ?
                                      dcell->periodic_neighbor_face_no(f) :
                                      dcell->neighbor_face_no(f))
                             ->index()] =
                            FaceCategory::locally_active_done_here;

                      // If neighbor is a ghost element (i.e.
                      // dcell->subdomain_id !
                      // dcell->neighbor(f)->subdomain_id()), we need to add its
                      // index into cell level list.
                      if (use_active_cells)
                        add_to_ghost =
                          (dcell->subdomain_id() != neighbor->subdomain_id());
                      else
                        add_to_ghost = (dcell->level_subdomain_id() !=
                                        neighbor->level_subdomain_id());
                    }
                  else if (additional_data.hold_all_faces_to_owned_cells ==
                           false)
                    {
                      // mark the cell to be close to the boundary
                      cells_close_to_boundary.emplace_back(i);
                    }
                  else
                    {
                      // add all cells to ghost layer...
                      face_is_owned[dcell->face(f)->index()] =
                        FaceCategory::ghosted;
                      if (use_active_cells)
                        {
                          if (neighbor->has_children())
                            for (unsigned int s = 0;
                                 s < dcell->face(f)->n_children();
                                 ++s)
                              if (dcell->at_boundary(f))
                                {
                                  if (dcell
                                        ->periodic_neighbor_child_on_subface(f,
                                                                             s)
                                        ->subdomain_id() !=
                                      dcell->subdomain_id())
                                    add_to_ghost = true;
                                }
                              else
                                {
                                  if (dcell->neighbor_child_on_subface(f, s)
                                        ->subdomain_id() !=
                                      dcell->subdomain_id())
                                    add_to_ghost = true;
                                }
                          else
                            add_to_ghost = (dcell->subdomain_id() !=
                                            neighbor->subdomain_id());
                        }
                      else
                        add_to_ghost = (dcell->level_subdomain_id() !=
                                        neighbor->level_subdomain_id());
                    }

                  if (add_to_ghost)
                    {
                      if (use_active_cells && neighbor->has_children())
                        for (unsigned int s = 0;
                             s < dcell->face(f)->n_children();
                             ++s)
                          {
                            typename dealii::Triangulation<dim>::cell_iterator
                              neighbor_child =
                                dcell->at_boundary(f) ?
                                  dcell->periodic_neighbor_child_on_subface(f,
                                                                            s) :
                                  dcell->neighbor_child_on_subface(f, s);
                            if (neighbor_child->subdomain_id() !=
                                dcell->subdomain_id())
                              ghost_cells.insert(
                                std::pair<unsigned int, unsigned int>(
                                  neighbor_child->level(),
                                  neighbor_child->index()));
                          }
                      else
                        ghost_cells.insert(
                          std::pair<unsigned int, unsigned int>(
                            neighbor->level(), neighbor->index()));
                      at_processor_boundary[i] = true;
                    }
                }
            }
        }

      // step 2: append the ghost cells at the end of the locally owned
      // cells
      for (const auto &ghost_cell : ghost_cells)
        cell_levels.push_back(ghost_cell);

      // step 3: clean up the cells close to the boundary
      std::sort(cells_close_to_boundary.begin(), cells_close_to_boundary.end());
      cells_close_to_boundary.erase(std::unique(cells_close_to_boundary.begin(),
                                                cells_close_to_boundary.end()),
                                    cells_close_to_boundary.end());
      std::vector<unsigned int> final_cells;
      final_cells.reserve(cells_close_to_boundary.size());
      for (unsigned int i = 0; i < cells_close_to_boundary.size(); ++i)
        if (at_processor_boundary[cells_close_to_boundary[i]] == false)
          final_cells.push_back(cells_close_to_boundary[i]);
      cells_close_to_boundary = std::move(final_cells);
    }



    template <int dim>
    void
    FaceSetup<dim>::generate_faces(
      const dealii::Triangulation<dim> &                        triangulation,
      const std::vector<std::pair<unsigned int, unsigned int>> &cell_levels,
      TaskInfo &                                                task_info)
    {
      // step 1: create the inverse map between cell iterators and the
      // cell_level_index field
      std::map<std::pair<unsigned int, unsigned int>, unsigned int>
        map_to_vectorized;
      for (unsigned int cell = 0; cell < cell_levels.size(); ++cell)
        if (cell == 0 || cell_levels[cell] != cell_levels[cell - 1])
          {
            typename dealii::Triangulation<dim>::cell_iterator dcell(
              &triangulation,
              cell_levels[cell].first,
              cell_levels[cell].second);
            std::pair<unsigned int, unsigned int> level_index(dcell->level(),
                                                              dcell->index());
            map_to_vectorized[level_index] = cell;
          }

      // step 2: fill the information about inner faces and boundary faces
      const unsigned int vectorization_length = task_info.vectorization_length;
      task_info.face_partition_data.resize(
        task_info.cell_partition_data.size() - 1, 0);
      task_info.boundary_partition_data.resize(
        task_info.cell_partition_data.size() - 1, 0);
      std::vector<unsigned char> face_visited(face_is_owned.size(), 0);
      for (unsigned int partition = 0;
           partition < task_info.cell_partition_data.size() - 2;
           ++partition)
        {
          unsigned int boundary_counter = 0;
          unsigned int inner_counter    = 0;
          for (unsigned int cell = task_info.cell_partition_data[partition] *
                                   vectorization_length;
               cell < task_info.cell_partition_data[partition + 1] *
                        vectorization_length;
               ++cell)
            if (cell == 0 || cell_levels[cell] != cell_levels[cell - 1])
              {
                typename dealii::Triangulation<dim>::cell_iterator dcell(
                  &triangulation,
                  cell_levels[cell].first,
                  cell_levels[cell].second);
                for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell;
                     ++f)
                  {
                    // boundary face
                    if (face_is_owned[dcell->face(f)->index()] ==
                        FaceCategory::locally_active_at_boundary)
                      {
                        Assert(dcell->at_boundary(f), ExcInternalError());
                        ++boundary_counter;
                        FaceToCellTopology<1> info;
                        info.cells_interior[0] = cell;
                        info.cells_exterior[0] = numbers::invalid_unsigned_int;
                        info.interior_face_no  = f;
                        info.exterior_face_no  = dcell->face(f)->boundary_id();
                        info.subface_index =
                          GeometryInfo<dim>::max_children_per_cell;
                        info.face_orientation = 0;
                        boundary_faces.push_back(info);

                        face_visited[dcell->face(f)->index()]++;
                      }
                    // interior face, including faces over periodic boundaries
                    else
                      {
                        typename dealii::Triangulation<dim>::cell_iterator
                          neighbor = dcell->neighbor_or_periodic_neighbor(f);
                        if (use_active_cells && neighbor->has_children())
                          {
                            for (unsigned int c = 0;
                                 c < dcell->face(f)->n_children();
                                 ++c)
                              {
                                typename dealii::Triangulation<
                                  dim>::cell_iterator neighbor_c =
                                  dcell->at_boundary(f) ?
                                    dcell->periodic_neighbor_child_on_subface(
                                      f, c) :
                                    dcell->neighbor_child_on_subface(f, c);
                                const types::subdomain_id neigh_domain =
                                  neighbor_c->subdomain_id();
                                const unsigned int neighbor_face_no =
                                  dcell->has_periodic_neighbor(f) ?
                                    dcell->periodic_neighbor_face_no(f) :
                                    dcell->neighbor_face_no(f);
                                if (neigh_domain != dcell->subdomain_id() ||
                                    face_visited
                                        [dcell->face(f)->child(c)->index()] ==
                                      1)
                                  {
                                    std::pair<unsigned int, unsigned int>
                                      level_index(neighbor_c->level(),
                                                  neighbor_c->index());
                                    if (face_is_owned
                                          [dcell->face(f)->child(c)->index()] ==
                                        FaceCategory::locally_active_done_here)
                                      {
                                        ++inner_counter;
                                        inner_faces.push_back(create_face(
                                          neighbor_face_no,
                                          neighbor_c,
                                          map_to_vectorized[level_index],
                                          dcell,
                                          cell));
                                      }
                                    else if (face_is_owned[dcell->face(f)
                                                             ->child(c)
                                                             ->index()] ==
                                             FaceCategory::ghosted)
                                      {
                                        inner_ghost_faces.push_back(create_face(
                                          neighbor_face_no,
                                          neighbor_c,
                                          map_to_vectorized[level_index],
                                          dcell,
                                          cell));
                                      }
                                    else
                                      Assert(
                                        face_is_owned[dcell->face(f)
                                                        ->index()] ==
                                            FaceCategory::
                                              locally_active_done_elsewhere ||
                                          face_is_owned[dcell->face(f)
                                                          ->index()] ==
                                            FaceCategory::ghosted,
                                        ExcInternalError());
                                  }
                                else
                                  {
                                    face_visited
                                      [dcell->face(f)->child(c)->index()] = 1;
                                  }
                              }
                          }
                        else
                          {
                            const types::subdomain_id my_domain =
                              use_active_cells ? dcell->subdomain_id() :
                                                 dcell->level_subdomain_id();
                            const types::subdomain_id neigh_domain =
                              use_active_cells ? neighbor->subdomain_id() :
                                                 neighbor->level_subdomain_id();
                            if (neigh_domain != my_domain ||
                                face_visited[dcell->face(f)->index()] == 1)
                              {
                                std::pair<unsigned int, unsigned int>
                                  level_index(neighbor->level(),
                                              neighbor->index());
                                if (face_is_owned[dcell->face(f)->index()] ==
                                    FaceCategory::locally_active_done_here)
                                  {
                                    ++inner_counter;
                                    inner_faces.push_back(create_face(
                                      f,
                                      dcell,
                                      cell,
                                      neighbor,
                                      map_to_vectorized[level_index]));
                                  }
                                else if (face_is_owned[dcell->face(f)
                                                         ->index()] ==
                                         FaceCategory::ghosted)
                                  {
                                    inner_ghost_faces.push_back(create_face(
                                      f,
                                      dcell,
                                      cell,
                                      neighbor,
                                      map_to_vectorized[level_index]));
                                  }
                              }
                            else
                              {
                                face_visited[dcell->face(f)->index()] = 1;
                                if (dcell->has_periodic_neighbor(f))
                                  face_visited
                                    [neighbor
                                       ->face(
                                         dcell->periodic_neighbor_face_no(f))
                                       ->index()] = 1;
                              }
                            if (face_is_owned[dcell->face(f)->index()] ==
                                FaceCategory::multigrid_refinement_edge)
                              {
                                refinement_edge_faces.push_back(
                                  create_face(f,
                                              dcell,
                                              cell,
                                              neighbor,
                                              refinement_edge_faces.size()));
                              }
                          }
                      }
                  }
              }
          task_info.face_partition_data[partition + 1] =
            task_info.face_partition_data[partition] + inner_counter;
          task_info.boundary_partition_data[partition + 1] =
            task_info.boundary_partition_data[partition] + boundary_counter;
        }
      task_info.ghost_face_partition_data.resize(2);
      task_info.ghost_face_partition_data[0] = 0;
      task_info.ghost_face_partition_data[1] = inner_ghost_faces.size();
      task_info.refinement_edge_face_partition_data.resize(2);
      task_info.refinement_edge_face_partition_data[0] = 0;
      task_info.refinement_edge_face_partition_data[1] =
        refinement_edge_faces.size();
    }



    template <int dim>
    FaceToCellTopology<1>
    FaceSetup<dim>::create_face(
      const unsigned int                                        face_no,
      const typename dealii::Triangulation<dim>::cell_iterator &cell,
      const unsigned int number_cell_interior,
      const typename dealii::Triangulation<dim>::cell_iterator &neighbor,
      const unsigned int number_cell_exterior)
    {
      FaceToCellTopology<1> info;
      info.cells_interior[0] = number_cell_interior;
      info.cells_exterior[0] = number_cell_exterior;
      info.interior_face_no  = face_no;
      if (cell->has_periodic_neighbor(face_no))
        info.exterior_face_no = cell->periodic_neighbor_face_no(face_no);
      else
        info.exterior_face_no = cell->neighbor_face_no(face_no);

      info.subface_index = GeometryInfo<dim>::max_children_per_cell;
      Assert(neighbor->level() <= cell->level(), ExcInternalError());
      if (cell->level() > neighbor->level())
        {
          if (cell->has_periodic_neighbor(face_no))
            info.subface_index =
              cell->periodic_neighbor_of_coarser_periodic_neighbor(face_no)
                .second;
          else
            info.subface_index =
              cell->neighbor_of_coarser_neighbor(face_no).second;
        }

      info.face_orientation = 0;
      const unsigned int left_face_orientation =
        !cell->face_orientation(face_no) + 2 * cell->face_flip(face_no) +
        4 * cell->face_rotation(face_no);
      const unsigned int right_face_orientation =
        !neighbor->face_orientation(info.exterior_face_no) +
        2 * neighbor->face_flip(info.exterior_face_no) +
        4 * neighbor->face_rotation(info.exterior_face_no);
      if (left_face_orientation != 0)
        {
          info.face_orientation = 8 + left_face_orientation;
          Assert(right_face_orientation == 0,
                 ExcMessage(
                   "Face seems to be wrongly oriented from both sides"));
        }
      else
        info.face_orientation = right_face_orientation;
      return info;
    }



    /**
     * This simple comparison for collect_faces_vectorization() identifies
     * faces of the same type, i.e., where all of the interior and exterior
     * face number, subface index and orientation are the same. This is used
     * to batch similar faces together for vectorization.
     */
    bool
    compare_faces_for_vectorization(const FaceToCellTopology<1> &face1,
                                    const FaceToCellTopology<1> &face2)
    {
      if (face1.interior_face_no != face2.interior_face_no)
        return false;
      if (face1.exterior_face_no != face2.exterior_face_no)
        return false;
      if (face1.subface_index != face2.subface_index)
        return false;
      if (face1.face_orientation != face2.face_orientation)
        return false;
      return true;
    }



    /**
     * This comparator is used within collect_faces_vectorization() to create
     * a sorting of FaceToCellTopology objects based on their
     * identifiers. This is used to obtain a good data locality when
     * processing the face integrals.
     */
    template <int length>
    struct FaceComparator
    {
      bool
      operator()(const FaceToCellTopology<length> &face1,
                 const FaceToCellTopology<length> &face2) const
      {
        for (unsigned int i = 0; i < length; ++i)
          if (face1.cells_interior[i] < face2.cells_interior[i])
            return true;
          else if (face1.cells_interior[i] > face2.cells_interior[i])
            return false;
        for (unsigned int i = 0; i < length; ++i)
          if (face1.cells_exterior[i] < face2.cells_exterior[i])
            return true;
          else if (face1.cells_exterior[i] > face2.cells_exterior[i])
            return false;
        if (face1.interior_face_no < face2.interior_face_no)
          return true;
        else if (face1.interior_face_no > face2.interior_face_no)
          return false;
        if (face1.exterior_face_no < face2.exterior_face_no)
          return true;
        else if (face1.exterior_face_no > face2.exterior_face_no)
          return false;

        // we do not need to check for subface_index and orientation because
        // those cannot be different if when all the other values are the
        // same.
        AssertDimension(face1.subface_index, face2.subface_index);
        AssertDimension(face1.face_orientation, face2.face_orientation);

        return false;
      }
    };



    template <int vectorization_width>
    void
    collect_faces_vectorization(
      const std::vector<FaceToCellTopology<1>> &faces_in,
      const std::vector<bool> &                 hard_vectorization_boundary,
      std::vector<unsigned int> &               face_partition_data,
      std::vector<FaceToCellTopology<vectorization_width>> &faces_out)
    {
      FaceToCellTopology<vectorization_width> macro_face;
      std::vector<std::vector<unsigned int>>  faces_type;

      unsigned int face_start = face_partition_data[0],
                   face_end   = face_partition_data[0];

      face_partition_data[0] = faces_out.size();
      for (unsigned int partition = 0;
           partition < face_partition_data.size() - 1;
           ++partition)
        {
          std::vector<std::vector<unsigned int>> new_faces_type;

          // start with the end point for the last partition
          face_start = face_end;
          face_end   = face_partition_data[partition + 1];

          // set the partitioner to the new vectorized lengths
          face_partition_data[partition + 1] = face_partition_data[partition];

          // loop over the faces in the current partition and reorder according
          // to the face type
          for (unsigned int face = face_start; face < face_end; ++face)
            {
              for (auto &face_type : faces_type)
                {
                  // Compare current face with first face of type type
                  if (compare_faces_for_vectorization(faces_in[face],
                                                      faces_in[face_type[0]]))
                    {
                      face_type.push_back(face);
                      goto face_found;
                    }
                }
              faces_type.emplace_back(1, face);
            face_found:
              {}
            }

          // insert new faces in sorted list to get good data locality
          std::set<FaceToCellTopology<vectorization_width>,
                   FaceComparator<vectorization_width>>
            new_faces;
          for (const auto &face_type : faces_type)
            {
              macro_face.interior_face_no =
                faces_in[face_type[0]].interior_face_no;
              macro_face.exterior_face_no =
                faces_in[face_type[0]].exterior_face_no;
              macro_face.subface_index = faces_in[face_type[0]].subface_index;
              macro_face.face_orientation =
                faces_in[face_type[0]].face_orientation;
              unsigned int               no_faces = face_type.size();
              std::vector<unsigned char> touched(no_faces, 0);

              // do two passes through the data. The first is to identify
              // similar faces within the same index range as the cells which
              // will allow for vectorized read operations, the second picks up
              // all the rest
              unsigned int n_vectorized = 0;
              for (unsigned int f = 0; f < no_faces; ++f)
                if (faces_in[face_type[f]].cells_interior[0] %
                      vectorization_width ==
                    0)
                  {
                    bool is_contiguous = true;
                    if (f + vectorization_width > no_faces)
                      is_contiguous = false;
                    else
                      for (unsigned int v = 1; v < vectorization_width; ++v)
                        if (faces_in[face_type[f + v]].cells_interior[0] !=
                            faces_in[face_type[f]].cells_interior[0] + v)
                          is_contiguous = false;
                    if (is_contiguous)
                      {
                        AssertIndexRange(f,
                                         face_type.size() -
                                           vectorization_width + 1);
                        for (unsigned int v = 0; v < vectorization_width; ++v)
                          {
                            macro_face.cells_interior[v] =
                              faces_in[face_type[f + v]].cells_interior[0];
                            macro_face.cells_exterior[v] =
                              faces_in[face_type[f + v]].cells_exterior[0];
                            touched[f + v] = 1;
                          }
                        new_faces.insert(macro_face);
                        f += vectorization_width - 1;
                        n_vectorized += vectorization_width;
                      }
                  }

              std::vector<unsigned int> untouched;
              untouched.reserve(no_faces - n_vectorized);
              for (unsigned int f = 0; f < no_faces; ++f)
                if (touched[f] == 0)
                  untouched.push_back(f);
              unsigned int v = 0;
              for (const auto f : untouched)
                {
                  macro_face.cells_interior[v] =
                    faces_in[face_type[f]].cells_interior[0];
                  macro_face.cells_exterior[v] =
                    faces_in[face_type[f]].cells_exterior[0];
                  ++v;
                  if (v == vectorization_width)
                    {
                      new_faces.insert(macro_face);
                      v = 0;
                    }
                }
              if (v > 0 && v < vectorization_width)
                {
                  // must add non-filled face
                  if (hard_vectorization_boundary[partition + 1] ||
                      partition == face_partition_data.size() - 2)
                    {
                      for (; v < vectorization_width; ++v)
                        {
                          // Dummy cell, not used
                          macro_face.cells_interior[v] =
                            numbers::invalid_unsigned_int;
                          macro_face.cells_exterior[v] =
                            numbers::invalid_unsigned_int;
                        }
                      new_faces.insert(macro_face);
                    }
                  else
                    {
                      // postpone to the next partition
                      std::vector<unsigned int> untreated(v);
                      for (unsigned int f = 0; f < v; ++f)
                        untreated[f] = face_type[*(untouched.end() - 1 - f)];
                      new_faces_type.push_back(untreated);
                    }
                }
            }

          // insert sorted list to vector of faces
          for (auto it = new_faces.begin(); it != new_faces.end(); ++it)
            faces_out.push_back(*it);
          face_partition_data[partition + 1] += new_faces.size();

          // set the faces that were left over to faces_type for the next round
          faces_type = std::move(new_faces_type);
        }

#  ifdef DEBUG
      // final safety checks
      for (const auto &face_type : faces_type)
        AssertDimension(face_type.size(), 0U);

      AssertDimension(faces_out.size(), face_partition_data.back());
      unsigned int nfaces = 0;
      for (unsigned int i = face_partition_data[0];
           i < face_partition_data.back();
           ++i)
        for (unsigned int v = 0; v < vectorization_width; ++v)
          nfaces +=
            (faces_out[i].cells_interior[v] != numbers::invalid_unsigned_int);
      AssertDimension(nfaces, faces_in.size());

      std::vector<std::pair<unsigned int, unsigned int>> in_faces, out_faces;
      for (const auto &face_in : faces_in)
        in_faces.emplace_back(face_in.cells_interior[0],
                              face_in.cells_exterior[0]);
      for (unsigned int i = face_partition_data[0];
           i < face_partition_data.back();
           ++i)
        for (unsigned int v = 0;
             v < vectorization_width &&
             faces_out[i].cells_interior[v] != numbers::invalid_unsigned_int;
             ++v)
          out_faces.emplace_back(faces_out[i].cells_interior[v],
                                 faces_out[i].cells_exterior[v]);
      std::sort(in_faces.begin(), in_faces.end());
      std::sort(out_faces.begin(), out_faces.end());
      AssertDimension(in_faces.size(), out_faces.size());
      for (unsigned int i = 0; i < in_faces.size(); ++i)
        {
          AssertDimension(in_faces[i].first, out_faces[i].first);
          AssertDimension(in_faces[i].second, out_faces[i].second);
        }
#  endif
    }

#endif // ifndef DOXYGEN

  } // namespace MatrixFreeFunctions
} // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
