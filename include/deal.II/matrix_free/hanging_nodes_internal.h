// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2022 by the deal.II authors
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

#ifndef dealii_hanging_nodes_internal_h
#define dealii_hanging_nodes_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/ndarray.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/hp/fe_collection.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * Here is the system for how we store constraint types in a binary mask.
     * This is not a complete contradiction-free system, i.e., there are
     * invalid states. You can use internal::MatrixFreeFunctions::check()
     * to check if the mask is in a valid state.
     *
     * If the mask is zero, there are no constraints. Then, there are three
     * different fields with one bit per dimension. The first field determines
     * the subcell, or the position of an element along each direction. The
     * second field determines if there is a constrained face with that
     * direction as normal. The last field determines if there is a
     * constrained edge in that direction (only valid in 3D).
     */
    enum class ConstraintKinds : std::uint16_t
    {
      // default: unconstrained cell
      unconstrained = 0,

      // subcell
      subcell_x = 1 << 0,
      subcell_y = 1 << 1,
      subcell_z = 1 << 2,

      // face is constrained
      face_x = 1 << 3,
      face_y = 1 << 4,
      face_z = 1 << 5,

      // edge is constrained
      edge_x = 1 << 6,
      edge_y = 1 << 7,
      edge_z = 1 << 8
    };



    /**
     * Type of the 8-bit representation of the refinement configuration.
     */
    using compressed_constraint_kind = std::uint8_t;



    /**
     * Value of compressed_constraint_kind for the unconstrained case.
     */
    constexpr compressed_constraint_kind
      unconstrained_compressed_constraint_kind = 0;



    /**
     * Check if the combinations of the bits in @p kind_in are valid.
     */
    inline bool
    check(const ConstraintKinds kind_in, const unsigned int dim)
    {
      const std::uint16_t kind    = static_cast<std::uint16_t>(kind_in);
      const std::uint16_t subcell = (kind >> 0) & 7;
      const std::uint16_t face    = (kind >> 3) & 7;
      const std::uint16_t edge    = (kind >> 6) & 7;

      if ((kind >> 9) > 0)
        return false;

      if (dim == 2)
        {
          if (edge > 0)
            return false; // in 2D there are no edge constraints

          if (subcell == 0 && face == 0)
            return true; // no constraints
          else if (0 < face)
            return true; // at least one face is constrained
        }
      else if (dim == 3)
        {
          if (subcell == 0 && face == 0 && edge == 0)
            return true; // no constraints
          else if (0 < face && edge == 0)
            return true; // at least one face is constrained
          else if (0 == face && 0 < edge)
            return true; // at least one edge is constrained
          else if ((face == edge) && (face == 1 || face == 2 || face == 4))
            return true; // one face and its orthogonal edge is constrained
        }

      return false;
    }



    /**
     * Compress a 9-bit representation of the refinement configuration
     * to a 8-bit representation.
     */
    inline compressed_constraint_kind
    compress(const ConstraintKinds kind_in, const unsigned int dim)
    {
      Assert(check(kind_in, dim), ExcInternalError());

      if (dim == 2)
        return static_cast<compressed_constraint_kind>(kind_in);

      if (kind_in == ConstraintKinds::unconstrained)
        return unconstrained_compressed_constraint_kind;

      const std::uint16_t kind    = static_cast<std::uint16_t>(kind_in);
      const std::uint16_t subcell = (kind >> 0) & 7;
      const std::uint16_t face    = (kind >> 3) & 7;
      const std::uint16_t edge    = (kind >> 6) & 7;

      return subcell + ((face > 0) << 3) + ((edge > 0) << 4) +
             (std::max(face, edge) << 5);
    }



    /**
     * Decompress a 8-bit representation of the refinement configuration
     * to a 9-bit representation.
     */
    inline ConstraintKinds
    decompress(const compressed_constraint_kind kind_in, const unsigned int dim)
    {
      if (dim == 2)
        return static_cast<ConstraintKinds>(kind_in);

      if (kind_in == unconstrained_compressed_constraint_kind)
        return ConstraintKinds::unconstrained;

      const std::uint16_t subcell = (kind_in >> 0) & 7;
      const std::uint16_t flag_0  = (kind_in >> 3) & 3;
      const std::uint16_t flag_1  = (kind_in >> 5) & 7;

      const auto result = static_cast<ConstraintKinds>(
        subcell + (((flag_0 & 0b01) ? flag_1 : 0) << 3) +
        (((flag_0 & 0b10) ? flag_1 : 0) << 6));

      Assert(check(result, dim), ExcInternalError());

      return result;
    }



    /**
     * Return the memory consumption in bytes of this enum class.
     */
    inline std::size_t
    memory_consumption(const ConstraintKinds &)
    {
      return sizeof(ConstraintKinds);
    }



    /**
     * Global operator which returns an object in which all bits are set which
     * are either set in the first or the second argument. This operator exists
     * since if it did not then the result of the bit-or <tt>operator |</tt>
     * would be an integer which would in turn trigger a compiler warning when
     * we tried to assign it to an object of type UpdateFlags.
     */
    DEAL_II_CUDA_HOST_DEV inline ConstraintKinds
    operator|(const ConstraintKinds f1, const ConstraintKinds f2)
    {
      return static_cast<ConstraintKinds>(static_cast<std::uint16_t>(f1) |
                                          static_cast<std::uint16_t>(f2));
    }



    /**
     * Global operator which sets the bits from the second argument also in the
     * first one.
     */
    DEAL_II_CUDA_HOST_DEV inline ConstraintKinds &
    operator|=(ConstraintKinds &f1, const ConstraintKinds f2)
    {
      f1 = f1 | f2;
      return f1;
    }



    /**
     * Global operator which checks inequality.
     */
    DEAL_II_CUDA_HOST_DEV inline bool
    operator!=(const ConstraintKinds f1, const ConstraintKinds f2)
    {
      return static_cast<std::uint16_t>(f1) != static_cast<std::uint16_t>(f2);
    }



    /**
     * Global operator which checks if the first argument is less than the
     * second.
     */
    DEAL_II_CUDA_HOST_DEV inline bool
    operator<(const ConstraintKinds f1, const ConstraintKinds f2)
    {
      return static_cast<std::uint16_t>(f1) < static_cast<std::uint16_t>(f2);
    }



    /**
     * Global operator which performs a binary and for the provided arguments.
     */
    DEAL_II_CUDA_HOST_DEV inline ConstraintKinds
    operator&(const ConstraintKinds f1, const ConstraintKinds f2)
    {
      return static_cast<ConstraintKinds>(static_cast<std::uint16_t>(f1) &
                                          static_cast<std::uint16_t>(f2));
    }



    /**
     * This class creates the mask used in the treatment of hanging nodes in
     * CUDAWrappers::MatrixFree.
     * The implementation of this class is explained in Section 3 of
     * @cite ljungkvist2017matrix and in Section 3.4 of
     * @cite kronbichler2019multigrid.
     */
    template <int dim>
    class HangingNodes
    {
    public:
      /**
       * Constructor.
       */
      HangingNodes(const Triangulation<dim> &triangualtion);

      /**
       * Compute the value of the constraint mask for a given cell.
       */
      template <typename CellIterator>
      bool
      setup_constraints(
        const CellIterator &                                      cell,
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
        const std::vector<std::vector<unsigned int>> &lexicographic_mapping,
        std::vector<types::global_dof_index> &        dof_indices,
        const ArrayView<ConstraintKinds> &            mask) const;

      /**
       * Compute the supported components of all entries of the given
       * hp::FECollection object.
       */
      static std::vector<std::vector<bool>>
      compute_supported_components(const dealii::hp::FECollection<dim> &fe);

      /**
       * Determine the refinement configuration of the given cell.
       */
      template <typename CellIterator>
      ConstraintKinds
      compute_refinement_configuration(const CellIterator &cell) const;

      /**
       * Update the DoF indices of a given cell for the given refinement
       * configuration and for the given components.
       */
      template <typename CellIterator>
      void
      update_dof_indices(
        const CellIterator &                                      cell,
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
        const std::vector<std::vector<unsigned int>> &lexicographic_mapping,
        const std::vector<std::vector<bool>> &        component_mask,
        const ConstraintKinds &                       refinement_configuration,
        std::vector<types::global_dof_index> &        dof_indices) const;

    private:
      /**
       * Set up line-to-cell mapping for edge constraints in 3D.
       */
      void
      setup_line_to_cell(const Triangulation<dim> &triangulation);

      void
      rotate_subface_index(int times, unsigned int &subface_index) const;

      void
      rotate_face(int                                   times,
                  unsigned int                          n_dofs_1d,
                  std::vector<types::global_dof_index> &dofs) const;

      unsigned int
      line_dof_idx(int          local_line,
                   unsigned int dof,
                   unsigned int n_dofs_1d) const;

      void
      transpose_face(const unsigned int                    fe_degree,
                     std::vector<types::global_dof_index> &dofs) const;

      void
      transpose_subface_index(unsigned int &subface) const;

      std::vector<std::vector<
        std::pair<typename Triangulation<dim>::cell_iterator, unsigned int>>>
        line_to_cells;

      const dealii::ndarray<unsigned int, 3, 2, 2> local_lines = {
        {{{{{7, 3}}, {{6, 2}}}},
         {{{{5, 1}}, {{4, 0}}}},
         {{{{11, 9}}, {{10, 8}}}}}};
    };



    template <int dim>
    inline HangingNodes<dim>::HangingNodes(
      const Triangulation<dim> &triangulation)
    {
      // Set up line-to-cell mapping for edge constraints (only if dim = 3 and
      // for pure hex meshes)
      if (triangulation.all_reference_cells_are_hyper_cube())
        setup_line_to_cell(triangulation);
    }



    template <int dim>
    inline void
    HangingNodes<dim>::setup_line_to_cell(
      const Triangulation<dim> &triangulation)
    {
      (void)triangulation;
    }



    template <>
    inline void
    HangingNodes<3>::setup_line_to_cell(const Triangulation<3> &triangulation)
    {
      // Check if we there are no hanging nodes on the current MPI process,
      // which we do by checking if the second finest level holds no active
      // non-artificial cell
      if (triangulation.n_levels() <= 1 ||
          std::none_of(triangulation.begin_active(triangulation.n_levels() - 2),
                       triangulation.end_active(triangulation.n_levels() - 2),
                       [](const CellAccessor<3, 3> &cell) {
                         return !cell.is_artificial();
                       }))
        return;

      const unsigned int n_raw_lines = triangulation.n_raw_lines();
      this->line_to_cells.resize(n_raw_lines);

      // In 3D, we can have DoFs on only an edge being constrained (e.g. in a
      // cartesian 2x2x2 grid, where only the upper left 2 cells are refined).
      // This sets up a helper data structure in the form of a mapping from
      // edges (i.e. lines) to neighboring cells.

      // Mapping from an edge to which children that share that edge.
      const unsigned int line_to_children[12][2] = {{0, 2},
                                                    {1, 3},
                                                    {0, 1},
                                                    {2, 3},
                                                    {4, 6},
                                                    {5, 7},
                                                    {4, 5},
                                                    {6, 7},
                                                    {0, 4},
                                                    {1, 5},
                                                    {2, 6},
                                                    {3, 7}};

      std::vector<std::vector<
        std::pair<typename Triangulation<3>::cell_iterator, unsigned int>>>
        line_to_inactive_cells(n_raw_lines);

      // First add active and inactive cells to their lines:
      for (const auto &cell : triangulation.cell_iterators())
        {
          for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_cell;
               ++line)
            {
              const unsigned int line_idx = cell->line(line)->index();
              if (cell->is_active())
                line_to_cells[line_idx].push_back(std::make_pair(cell, line));
              else
                line_to_inactive_cells[line_idx].push_back(
                  std::make_pair(cell, line));
            }
        }

      // Now, we can access edge-neighboring active cells on same level to also
      // access of an edge to the edges "children". These are found from looking
      // at the corresponding edge of children of inactive edge neighbors.
      for (unsigned int line_idx = 0; line_idx < n_raw_lines; ++line_idx)
        {
          if ((line_to_cells[line_idx].size() > 0) &&
              line_to_inactive_cells[line_idx].size() > 0)
            {
              // We now have cells to add (active ones) and edges to which they
              // should be added (inactive cells).
              const auto &inactive_cell =
                line_to_inactive_cells[line_idx][0].first;
              const unsigned int neighbor_line =
                line_to_inactive_cells[line_idx][0].second;

              for (unsigned int c = 0; c < 2; ++c)
                {
                  const auto &child =
                    inactive_cell->child(line_to_children[neighbor_line][c]);
                  const unsigned int child_line_idx =
                    child->line(neighbor_line)->index();

                  // Now add all active cells
                  for (const auto &cl : line_to_cells[line_idx])
                    line_to_cells[child_line_idx].push_back(cl);
                }
            }
        }
    }



    template <int dim>
    inline std::vector<std::vector<bool>>
    HangingNodes<dim>::compute_supported_components(
      const dealii::hp::FECollection<dim> &fe_collection)
    {
      std::vector<std::vector<bool>> supported_components(
        fe_collection.size(),
        std::vector<bool>(fe_collection.n_components(), false));

      for (unsigned int i = 0; i < fe_collection.size(); ++i)
        {
          for (unsigned int base_element_index = 0, comp = 0;
               base_element_index < fe_collection[i].n_base_elements();
               ++base_element_index)
            for (unsigned int c = 0;
                 c < fe_collection[i].element_multiplicity(base_element_index);
                 ++c, ++comp)
              if (dim == 1 || dynamic_cast<const FE_Q<dim> *>(
                                &fe_collection[i].base_element(
                                  base_element_index)) == nullptr)
                supported_components[i][comp] = false;
              else
                supported_components[i][comp] = true;
        }

      return supported_components;
    }



    template <int dim>
    template <typename CellIterator>
    inline ConstraintKinds
    HangingNodes<dim>::compute_refinement_configuration(
      const CellIterator &cell) const
    {
      // TODO: for simplex or mixed meshes: nothing to do
      if ((dim == 3 && line_to_cells.size() == 0) ||
          (cell->reference_cell().is_hyper_cube() == false))
        return ConstraintKinds::unconstrained;

      if (cell->level() == 0)
        return ConstraintKinds::unconstrained;

      const std::uint16_t subcell =
        cell->parent()->child_iterator_to_index(cell);
      const std::uint16_t subcell_x = (subcell >> 0) & 1;
      const std::uint16_t subcell_y = (subcell >> 1) & 1;
      const std::uint16_t subcell_z = (subcell >> 2) & 1;

      std::uint16_t face = 0;
      std::uint16_t edge = 0;

      for (unsigned int direction = 0; direction < dim; ++direction)
        {
          const auto side    = (subcell >> direction) & 1U;
          const auto face_no = direction * 2 + side;

          // ignore if at boundary
          if (cell->at_boundary(face_no))
            continue;

          const auto &neighbor = cell->neighbor(face_no);

          // ignore neighbors that are artificial or have the same level or
          // have children
          if (neighbor->has_children() || neighbor->is_artificial() ||
              neighbor->level() == cell->level())
            continue;

          face |= 1 << direction;
        }

      if (dim == 3)
        for (unsigned int direction = 0; direction < dim; ++direction)
          if (face == 0 || face == (1 << direction))
            {
              const unsigned int line_no =
                direction == 0 ?
                  (local_lines[0][subcell_y == 0][subcell_z == 0]) :
                  (direction == 1 ?
                     (local_lines[1][subcell_x == 0][subcell_z == 0]) :
                     (local_lines[2][subcell_x == 0][subcell_y == 0]));

              const unsigned int line_index = cell->line(line_no)->index();

              const auto edge_neighbor =
                std::find_if(line_to_cells[line_index].begin(),
                             line_to_cells[line_index].end(),
                             [&cell](const auto &edge_neighbor) {
                               return edge_neighbor.first->is_artificial() ==
                                        false &&
                                      edge_neighbor.first->level() <
                                        cell->level();
                             });

              if (edge_neighbor == line_to_cells[line_index].end())
                continue;

              edge |= 1 << direction;
            }

      if ((face == 0) && (edge == 0))
        return ConstraintKinds::unconstrained;

      const std::uint16_t inverted_subcell = (subcell ^ (dim == 2 ? 3 : 7));

      const auto refinement_configuration = static_cast<ConstraintKinds>(
        inverted_subcell + (face << 3) + (edge << 6));
      Assert(check(refinement_configuration, dim), ExcInternalError());
      return refinement_configuration;
    }



    template <int dim>
    template <typename CellIterator>
    inline void
    HangingNodes<dim>::update_dof_indices(
      const CellIterator &                                      cell,
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
      const std::vector<std::vector<unsigned int>> &lexicographic_mapping,
      const std::vector<std::vector<bool>> &        supported_components,
      const ConstraintKinds &                       refinement_configuration,
      std::vector<types::global_dof_index> &        dof_indices) const
    {
      if (std::find(supported_components[cell->active_fe_index()].begin(),
                    supported_components[cell->active_fe_index()].end(),
                    true) ==
          supported_components[cell->active_fe_index()].end())
        return;

      const auto &fe = cell->get_fe();

      AssertDimension(fe.n_unique_faces(), 1);

      std::vector<std::vector<unsigned int>>
        component_to_system_index_face_array(fe.n_components());

      for (unsigned int i = 0; i < fe.n_dofs_per_face(0); ++i)
        component_to_system_index_face_array
          [fe.face_system_to_component_index(i, /*face_no=*/0).first]
            .push_back(i);

      std::vector<unsigned int> idx_offset = {0};

      for (unsigned int base_element_index = 0;
           base_element_index < cell->get_fe().n_base_elements();
           ++base_element_index)
        for (unsigned int c = 0;
             c < cell->get_fe().element_multiplicity(base_element_index);
             ++c)
          idx_offset.push_back(
            idx_offset.back() +
            cell->get_fe().base_element(base_element_index).n_dofs_per_cell());

      std::vector<types::global_dof_index> neighbor_dofs_all(idx_offset.back());
      std::vector<types::global_dof_index> neighbor_dofs_all_temp(
        idx_offset.back());
      std::vector<types::global_dof_index> neighbor_dofs_face(
        fe.n_dofs_per_face(/*face_no=*/0));


      const auto get_face_idx = [](const auto n_dofs_1d,
                                   const auto face_no,
                                   const auto i,
                                   const auto j) -> unsigned int {
        const auto direction = face_no / 2;
        const auto side      = face_no % 2;
        const auto offset    = (side == 1) ? (n_dofs_1d - 1) : 0;

        if (dim == 2)
          return (direction == 0) ? (n_dofs_1d * i + offset) :
                                    (n_dofs_1d * offset + i);
        else if (dim == 3)
          switch (direction)
            {
              case 0:
                return n_dofs_1d * n_dofs_1d * i + n_dofs_1d * j + offset;
              case 1:
                return n_dofs_1d * n_dofs_1d * j + n_dofs_1d * offset + i;
              case 2:
                return n_dofs_1d * n_dofs_1d * offset + n_dofs_1d * i + j;
              default:
                Assert(false, ExcNotImplemented());
            }

        Assert(false, ExcNotImplemented());

        return 0;
      };

      const std::uint16_t kind =
        static_cast<std::uint16_t>(refinement_configuration);
      const std::uint16_t subcell   = (kind >> 0) & 7;
      const std::uint16_t subcell_x = (subcell >> 0) & 1;
      const std::uint16_t subcell_y = (subcell >> 1) & 1;
      const std::uint16_t subcell_z = (subcell >> 2) & 1;
      const std::uint16_t face      = (kind >> 3) & 7;
      const std::uint16_t edge      = (kind >> 6) & 7;

      for (unsigned int direction = 0; direction < dim; ++direction)
        if ((face >> direction) & 1U)
          {
            const auto side    = ((subcell >> direction) & 1U) == 0;
            const auto face_no = direction * 2 + side;

            // read DoFs of parent of face, ...
            cell->neighbor(face_no)
              ->face(cell->neighbor_face_no(face_no))
              ->get_dof_indices(neighbor_dofs_face,
                                cell->neighbor(face_no)->active_fe_index());

            // ... convert the global DoFs to serial ones, and ...
            if (partitioner)
              for (auto &index : neighbor_dofs_face)
                index = partitioner->global_to_local(index);

            for (unsigned int base_element_index = 0, comp = 0;
                 base_element_index < cell->get_fe().n_base_elements();
                 ++base_element_index)
              for (unsigned int c = 0;
                   c < cell->get_fe().element_multiplicity(base_element_index);
                   ++c, ++comp)
                {
                  if (supported_components[cell->active_fe_index()][comp] ==
                      false)
                    continue;

                  const unsigned int n_dofs_1d =
                    cell->get_fe()
                      .base_element(base_element_index)
                      .tensor_degree() +
                    1;
                  const unsigned int dofs_per_face =
                    Utilities::pow(n_dofs_1d, dim - 1);
                  std::vector<types::global_dof_index> neighbor_dofs(
                    dofs_per_face);
                  const auto lex_face_mapping =
                    FETools::lexicographic_to_hierarchic_numbering<dim - 1>(
                      n_dofs_1d - 1);

                  // ... extract the DoFs of the current component
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    neighbor_dofs[i] = neighbor_dofs_face
                      [component_to_system_index_face_array[comp][i]];

                  // fix DoFs depending on orientation, flip, and rotation
                  if (dim == 2)
                    {
                      // TODO: for mixed meshes we need to take care of
                      // orientation here
                      Assert(cell->face_orientation(face_no),
                             ExcNotImplemented());
                    }
                  else if (dim == 3)
                    {
                      int rotate = 0;                   // TODO
                      if (cell->face_rotation(face_no)) //
                        rotate -= 1;                    //
                      if (cell->face_flip(face_no))     //
                        rotate -= 2;                    //

                      rotate_face(rotate, n_dofs_1d, neighbor_dofs);

                      if (cell->face_orientation(face_no) == false)
                        transpose_face(n_dofs_1d - 1, neighbor_dofs);
                    }
                  else
                    {
                      Assert(false, ExcNotImplemented());
                    }

                  // update DoF map
                  for (unsigned int i = 0, k = 0; i < n_dofs_1d; ++i)
                    for (unsigned int j = 0; j < (dim == 2 ? 1 : n_dofs_1d);
                         ++j, ++k)
                      dof_indices[get_face_idx(n_dofs_1d, face_no, i, j) +
                                  idx_offset[comp]] =
                        neighbor_dofs[lex_face_mapping[k]];
                }
          }

      if (dim == 3)
        for (unsigned int direction = 0; direction < dim; ++direction)
          if ((edge >> direction) & 1U)
            {
              const unsigned int line_no =
                direction == 0 ?
                  (local_lines[0][subcell_y][subcell_z]) :
                  (direction == 1 ? (local_lines[1][subcell_x][subcell_z]) :
                                    (local_lines[2][subcell_x][subcell_y]));

              const unsigned int line_index = cell->line(line_no)->index();

              const auto edge_neighbor =
                std::find_if(line_to_cells[line_index].begin(),
                             line_to_cells[line_index].end(),
                             [&cell](const auto &edge_neighbor) {
                               return edge_neighbor.first->is_artificial() ==
                                        false &&
                                      edge_neighbor.first->level() <
                                        cell->level();
                             });

              if (edge_neighbor == line_to_cells[line_index].end())
                continue;

              const auto neighbor_cell       = edge_neighbor->first;
              const auto local_line_neighbor = edge_neighbor->second;

              DoFCellAccessor<dim, dim, false>(
                &neighbor_cell->get_triangulation(),
                neighbor_cell->level(),
                neighbor_cell->index(),
                &cell->get_dof_handler())
                .get_dof_indices(neighbor_dofs_all);

              if (partitioner)
                for (auto &index : neighbor_dofs_all)
                  index = partitioner->global_to_local(index);

              for (unsigned int i = 0; i < neighbor_dofs_all_temp.size(); ++i)
                neighbor_dofs_all_temp[i] = neighbor_dofs_all
                  [lexicographic_mapping[cell->active_fe_index()][i]];

              const bool flipped =
                cell->line_orientation(line_no) !=
                neighbor_cell->line_orientation(local_line_neighbor);

              for (unsigned int base_element_index = 0, comp = 0;
                   base_element_index < cell->get_fe().n_base_elements();
                   ++base_element_index)
                for (unsigned int c = 0;
                     c <
                     cell->get_fe().element_multiplicity(base_element_index);
                     ++c, ++comp)
                  {
                    if (supported_components[cell->active_fe_index()][comp] ==
                        false)
                      continue;

                    const unsigned int n_dofs_1d =
                      cell->get_fe()
                        .base_element(base_element_index)
                        .tensor_degree() +
                      1;

                    for (unsigned int i = 0; i < n_dofs_1d; ++i)
                      dof_indices[line_dof_idx(line_no, i, n_dofs_1d) +
                                  idx_offset[comp]] = neighbor_dofs_all_temp
                        [line_dof_idx(local_line_neighbor,
                                      flipped ? (n_dofs_1d - 1 - i) : i,
                                      n_dofs_1d) +
                         idx_offset[comp]];
                  }
            }
    }



    template <int dim>
    template <typename CellIterator>
    inline bool
    HangingNodes<dim>::setup_constraints(
      const CellIterator &                                      cell,
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
      const std::vector<std::vector<unsigned int>> &lexicographic_mapping,
      std::vector<types::global_dof_index> &        dof_indices,
      const ArrayView<ConstraintKinds> &            masks) const
    {
      // 1) check if finite elements support fast hanging-node algorithm
      const auto supported_components = compute_supported_components(
        cell->get_dof_handler().get_fe_collection());

      if ([](const auto &supported_components) {
            return std::none_of(supported_components.begin(),
                                supported_components.end(),
                                [](const auto &a) {
                                  return *std::max_element(a.begin(), a.end());
                                });
          }(supported_components))
        return false;

      // 2) determine the refinement configuration of the cell
      const auto refinement_configuration =
        compute_refinement_configuration(cell);

      if (refinement_configuration == ConstraintKinds::unconstrained)
        return false;

      // 3) update DoF indices of cell for specified components
      update_dof_indices(cell,
                         partitioner,
                         lexicographic_mapping,
                         supported_components,
                         refinement_configuration,
                         dof_indices);

      // 4)  TODO: copy refinement configuration to all components
      for (unsigned int c = 0; c < supported_components[0].size(); ++c)
        if (supported_components[cell->active_fe_index()][c])
          masks[c] = refinement_configuration;

      return true;
    }



    template <int dim>
    inline void
    HangingNodes<dim>::rotate_subface_index(int           times,
                                            unsigned int &subface_index) const
    {
      const unsigned int rot_mapping[4] = {2, 0, 3, 1};

      times = times % 4;
      times = times < 0 ? times + 4 : times;
      for (int t = 0; t < times; ++t)
        subface_index = rot_mapping[subface_index];
    }



    template <int dim>
    inline void
    HangingNodes<dim>::rotate_face(
      int                                   times,
      unsigned int                          n_dofs_1d,
      std::vector<types::global_dof_index> &dofs) const
    {
      const unsigned int rot_mapping[4] = {2, 0, 3, 1};

      times = times % 4;
      times = times < 0 ? times + 4 : times;

      std::vector<types::global_dof_index> copy(dofs.size());
      for (int t = 0; t < times; ++t)
        {
          std::swap(copy, dofs);

          // Vertices
          for (unsigned int i = 0; i < 4; ++i)
            dofs[rot_mapping[i]] = copy[i];

          // Edges
          const unsigned int n_int  = n_dofs_1d - 2;
          unsigned int       offset = 4;
          for (unsigned int i = 0; i < n_int; ++i)
            {
              // Left edge
              dofs[offset + i] = copy[offset + 2 * n_int + (n_int - 1 - i)];
              // Right edge
              dofs[offset + n_int + i] =
                copy[offset + 3 * n_int + (n_int - 1 - i)];
              // Bottom edge
              dofs[offset + 2 * n_int + i] = copy[offset + n_int + i];
              // Top edge
              dofs[offset + 3 * n_int + i] = copy[offset + i];
            }

          // Interior points
          offset += 4 * n_int;

          for (unsigned int i = 0; i < n_int; ++i)
            for (unsigned int j = 0; j < n_int; ++j)
              dofs[offset + i * n_int + j] =
                copy[offset + j * n_int + (n_int - 1 - i)];
        }
    }



    template <int dim>
    inline unsigned int
    HangingNodes<dim>::line_dof_idx(int          local_line,
                                    unsigned int dof,
                                    unsigned int n_dofs_1d) const
    {
      unsigned int x, y, z;

      const unsigned int fe_degree = n_dofs_1d - 1;

      if (local_line < 8)
        {
          x = (local_line % 4 == 0) ? 0 :
              (local_line % 4 == 1) ? fe_degree :
                                      dof;
          y = (local_line % 4 == 2) ? 0 :
              (local_line % 4 == 3) ? fe_degree :
                                      dof;
          z = (local_line / 4) * fe_degree;
        }
      else
        {
          x = ((local_line - 8) % 2) * fe_degree;
          y = ((local_line - 8) / 2) * fe_degree;
          z = dof;
        }

      return n_dofs_1d * n_dofs_1d * z + n_dofs_1d * y + x;
    }



    template <int dim>
    inline void
    HangingNodes<dim>::transpose_face(
      const unsigned int                    fe_degree,
      std::vector<types::global_dof_index> &dofs) const
    {
      const std::vector<types::global_dof_index> copy(dofs);

      // Vertices
      dofs[1] = copy[2];
      dofs[2] = copy[1];

      // Edges
      const unsigned int n_int  = fe_degree - 1;
      unsigned int       offset = 4;
      for (unsigned int i = 0; i < n_int; ++i)
        {
          // Right edge
          dofs[offset + i] = copy[offset + 2 * n_int + i];
          // Left edge
          dofs[offset + n_int + i] = copy[offset + 3 * n_int + i];
          // Bottom edge
          dofs[offset + 2 * n_int + i] = copy[offset + i];
          // Top edge
          dofs[offset + 3 * n_int + i] = copy[offset + n_int + i];
        }

      // Interior
      offset += 4 * n_int;
      for (unsigned int i = 0; i < n_int; ++i)
        for (unsigned int j = 0; j < n_int; ++j)
          dofs[offset + i * n_int + j] = copy[offset + j * n_int + i];
    }



    template <int dim>
    void
    HangingNodes<dim>::transpose_subface_index(unsigned int &subface) const
    {
      if (subface == 1)
        subface = 2;
      else if (subface == 2)
        subface = 1;
    }
  } // namespace MatrixFreeFunctions
} // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
