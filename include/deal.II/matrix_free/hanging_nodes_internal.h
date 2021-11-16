// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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

#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

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
      type_x = 1 << 0,
      type_y = 1 << 1,
      type_z = 1 << 2,

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
     * Check if the combinations of the bits in @p kind_in are valid.
     */
    inline bool
    check(const ConstraintKinds &kind_in, const unsigned int dim)
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
        const std::vector<unsigned int> &     lexicographic_mapping,
        std::vector<types::global_dof_index> &dof_indices,
        const ArrayView<ConstraintKinds> &    mask) const;

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
    template <typename CellIterator>
    inline bool
    HangingNodes<dim>::setup_constraints(
      const CellIterator &                                      cell,
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
      const std::vector<unsigned int> &     lexicographic_mapping,
      std::vector<types::global_dof_index> &dof_indices,
      const ArrayView<ConstraintKinds> &    masks) const
    {
      bool cell_has_hanging_node_constraints = false;

      // for simplex or mixed meshes: nothing to do
      if (dim == 3 && line_to_cells.size() == 0)
        return cell_has_hanging_node_constraints;

      const auto &fe = cell->get_fe();

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

      for (unsigned int base_element_index = 0, comp = 0;
           base_element_index < cell->get_fe().n_base_elements();
           ++base_element_index)
        for (unsigned int c = 0;
             c < cell->get_fe().element_multiplicity(base_element_index);
             ++c, ++comp)
          {
            auto &mask = masks[comp];
            mask       = ConstraintKinds::unconstrained;

            const auto &fe_base =
              cell->get_fe().base_element(base_element_index);

            if (dim == 1 ||
                dynamic_cast<const FE_Q<dim> *>(&fe_base) == nullptr)
              continue;

            const unsigned int fe_degree = fe_base.tensor_degree();
            const unsigned int n_dofs_1d = fe_degree + 1;
            const unsigned int dofs_per_face =
              Utilities::fixed_power<dim - 1>(n_dofs_1d);

            std::vector<types::global_dof_index> neighbor_dofs_all(
              idx_offset.back());
            std::vector<types::global_dof_index> neighbor_dofs_all_temp(
              idx_offset.back());

            std::vector<types::global_dof_index> neighbor_dofs(dofs_per_face);

            const auto lex_face_mapping =
              FETools::lexicographic_to_hierarchic_numbering<dim - 1>(
                fe_degree);

            for (const unsigned int face : GeometryInfo<dim>::face_indices())
              {
                if ((!cell->at_boundary(face)) &&
                    (cell->neighbor(face)->has_children() == false))
                  {
                    const auto &neighbor = cell->neighbor(face);

                    if (neighbor->is_artificial())
                      continue;

                    // Neighbor is coarser than us, i.e., face is constrained
                    if (neighbor->level() < cell->level())
                      {
                        const unsigned int neighbor_face =
                          cell->neighbor_face_no(face);

                        // Find position of face on neighbor
                        unsigned int subface = 0;
                        for (;
                             subface < GeometryInfo<dim>::max_children_per_face;
                             ++subface)
                          if (neighbor->neighbor_child_on_subface(neighbor_face,
                                                                  subface) ==
                              cell)
                            break;

                        // Get indices to read
                        DoFAccessor<dim - 1, dim, dim, false>(
                          &neighbor->face(neighbor_face)->get_triangulation(),
                          neighbor->face(neighbor_face)->level(),
                          neighbor->face(neighbor_face)->index(),
                          &cell->get_dof_handler())
                          .get_dof_indices(neighbor_dofs_all);

                        for (unsigned int i = 0; i < dofs_per_face; ++i)
                          neighbor_dofs[i] = neighbor_dofs_all
                            [component_to_system_index_face_array[comp][i]];

                        // If the vector is distributed, we need to transform
                        // the global indices to local ones.
                        if (partitioner)
                          for (auto &index : neighbor_dofs)
                            index = partitioner->global_to_local(index);

                        if (dim == 2)
                          {
                            if (face < 2)
                              {
                                mask |= ConstraintKinds::face_x;
                                if (face == 0)
                                  mask |= ConstraintKinds::type_x;
                                if (subface == 0)
                                  mask |= ConstraintKinds::type_y;
                              }
                            else
                              {
                                mask |= ConstraintKinds::face_y;
                                if (face == 2)
                                  mask |= ConstraintKinds::type_y;
                                if (subface == 0)
                                  mask |= ConstraintKinds::type_x;
                              }

                            // Reorder neighbor_dofs and copy into faceth face
                            // of dof_indices

                            // Offset if upper/right face
                            unsigned int offset =
                              (face % 2 == 1) ? fe_degree : 0;

                            for (unsigned int i = 0; i < n_dofs_1d; ++i)
                              {
                                unsigned int idx = 0;
                                // If X-line, i.e., if y = 0 or y = fe_degree
                                if (face > 1)
                                  idx = n_dofs_1d * offset + i;
                                // If Y-line, i.e., if x = 0 or x = fe_degree
                                else
                                  idx = n_dofs_1d * i + offset;

                                dof_indices[idx + idx_offset[comp]] =
                                  neighbor_dofs[lex_face_mapping[i]];
                              }
                          }
                        else if (dim == 3)
                          {
                            const bool transpose =
                              !(cell->face_orientation(face));

                            int rotate = 0;

                            if (cell->face_rotation(face))
                              rotate -= 1;
                            if (cell->face_flip(face))
                              rotate -= 2;

                            rotate_face(rotate, n_dofs_1d, neighbor_dofs);
                            rotate_subface_index(rotate, subface);

                            if (transpose)
                              {
                                transpose_face(fe_degree, neighbor_dofs);
                                transpose_subface_index(subface);
                              }

                            // YZ-plane
                            if (face < 2)
                              {
                                mask |= ConstraintKinds::face_x;
                                if (face == 0)
                                  mask |= ConstraintKinds::type_x;
                                if (subface % 2 == 0)
                                  mask |= ConstraintKinds::type_y;
                                if (subface / 2 == 0)
                                  mask |= ConstraintKinds::type_z;
                              }
                            // XZ-plane
                            else if (face < 4)
                              {
                                mask |= ConstraintKinds::face_y;
                                if (face == 2)
                                  mask |= ConstraintKinds::type_y;
                                if (subface % 2 == 0)
                                  mask |= ConstraintKinds::type_z;
                                if (subface / 2 == 0)
                                  mask |= ConstraintKinds::type_x;
                              }
                            // XY-plane
                            else
                              {
                                mask |= ConstraintKinds::face_z;
                                if (face == 4)
                                  mask |= ConstraintKinds::type_z;
                                if (subface % 2 == 0)
                                  mask |= ConstraintKinds::type_x;
                                if (subface / 2 == 0)
                                  mask |= ConstraintKinds::type_y;
                              }

                            // Offset if upper/right/back face
                            unsigned int offset =
                              (face % 2 == 1) ? fe_degree : 0;

                            for (unsigned int i = 0; i < n_dofs_1d; ++i)
                              {
                                for (unsigned int j = 0; j < n_dofs_1d; ++j)
                                  {
                                    unsigned int idx = 0;
                                    // If YZ-plane, i.e., if x = 0 or x =
                                    // fe_degree, and orientation standard
                                    if (face < 2)
                                      idx = n_dofs_1d * n_dofs_1d * i +
                                            n_dofs_1d * j + offset;
                                    // If XZ-plane, i.e., if y = 0 or y =
                                    // fe_degree, and orientation standard
                                    else if (face < 4)
                                      idx = n_dofs_1d * n_dofs_1d * j +
                                            n_dofs_1d * offset + i;
                                    // If XY-plane, i.e., if z = 0 or z =
                                    // fe_degree, and orientation standard
                                    else
                                      idx = n_dofs_1d * n_dofs_1d * offset +
                                            n_dofs_1d * i + j;

                                    dof_indices[idx + idx_offset[comp]] =
                                      neighbor_dofs
                                        [lex_face_mapping[n_dofs_1d * i + j]];
                                  }
                              }
                          }
                        else
                          ExcNotImplemented();
                      }
                  }
              }

            // In 3D we can have a situation where only DoFs on an edge are
            // constrained. Append these here.
            if (dim == 3)
              {
                // For each line on cell, which faces does it belong to, what is
                // the edge mask, what is the types of the faces it belong to,
                // and what is the type along the edge.
                const ConstraintKinds line_to_edge[12][4] = {
                  {ConstraintKinds::face_x | ConstraintKinds::face_z,
                   ConstraintKinds::edge_y,
                   ConstraintKinds::type_x | ConstraintKinds::type_z,
                   ConstraintKinds::type_y},
                  {ConstraintKinds::face_x | ConstraintKinds::face_z,
                   ConstraintKinds::edge_y,
                   ConstraintKinds::type_z,
                   ConstraintKinds::type_y},
                  {ConstraintKinds::face_y | ConstraintKinds::face_z,
                   ConstraintKinds::edge_x,
                   ConstraintKinds::type_y | ConstraintKinds::type_z,
                   ConstraintKinds::type_x},
                  {ConstraintKinds::face_y | ConstraintKinds::face_z,
                   ConstraintKinds::edge_x,
                   ConstraintKinds::type_z,
                   ConstraintKinds::type_x},
                  {ConstraintKinds::face_x | ConstraintKinds::face_z,
                   ConstraintKinds::edge_y,
                   ConstraintKinds::type_x,
                   ConstraintKinds::type_y},
                  {ConstraintKinds::face_x | ConstraintKinds::face_z,
                   ConstraintKinds::edge_y,
                   ConstraintKinds::unconstrained,
                   ConstraintKinds::type_y},
                  {ConstraintKinds::face_y | ConstraintKinds::face_z,
                   ConstraintKinds::edge_x,
                   ConstraintKinds::type_y,
                   ConstraintKinds::type_x},
                  {ConstraintKinds::face_y | ConstraintKinds::face_z,
                   ConstraintKinds::edge_x,
                   ConstraintKinds::unconstrained,
                   ConstraintKinds::type_x},
                  {ConstraintKinds::face_x | ConstraintKinds::face_y,
                   ConstraintKinds::edge_z,
                   ConstraintKinds::type_x | ConstraintKinds::type_y,
                   ConstraintKinds::type_z},
                  {ConstraintKinds::face_x | ConstraintKinds::face_y,
                   ConstraintKinds::edge_z,
                   ConstraintKinds::type_y,
                   ConstraintKinds::type_z},
                  {ConstraintKinds::face_x | ConstraintKinds::face_y,
                   ConstraintKinds::edge_z,
                   ConstraintKinds::type_x,
                   ConstraintKinds::type_z},
                  {ConstraintKinds::face_x | ConstraintKinds::face_y,
                   ConstraintKinds::edge_z,
                   ConstraintKinds::unconstrained,
                   ConstraintKinds::type_z}};

                for (unsigned int local_line = 0;
                     local_line < GeometryInfo<dim>::lines_per_cell;
                     ++local_line)
                  {
                    // If we don't already have a constraint for as part of a
                    // face
                    if ((mask & line_to_edge[local_line][0]) ==
                        ConstraintKinds::unconstrained)
                      {
                        // For each cell which share that edge
                        const unsigned int line =
                          cell->line(local_line)->index();
                        for (const auto &edge_neighbor : line_to_cells[line])
                          {
                            // If one of them is coarser than us
                            const auto neighbor_cell = edge_neighbor.first;

                            if (neighbor_cell->is_artificial())
                              continue;

                            if (neighbor_cell->level() < cell->level())
                              {
                                const unsigned int local_line_neighbor =
                                  edge_neighbor.second;
                                mask |= line_to_edge[local_line][1] |
                                        line_to_edge[local_line][2];

                                bool flipped = false;
                                if (cell->line(local_line)->vertex_index(0) ==
                                    neighbor_cell->line(local_line_neighbor)
                                      ->vertex_index(0))
                                  {
                                    // Assuming line directions match axes
                                    // directions, we have an unflipped edge of
                                    // first type
                                    mask |= line_to_edge[local_line][3];
                                  }
                                else if (cell->line(local_line)
                                           ->vertex_index(1) ==
                                         neighbor_cell
                                           ->line(local_line_neighbor)
                                           ->vertex_index(1))
                                  {
                                    // We have an unflipped edge of second type
                                  }
                                else if (cell->line(local_line)
                                           ->vertex_index(1) ==
                                         neighbor_cell
                                           ->line(local_line_neighbor)
                                           ->vertex_index(0))
                                  {
                                    // We have a flipped edge of second type
                                    flipped = true;
                                  }
                                else if (cell->line(local_line)
                                           ->vertex_index(0) ==
                                         neighbor_cell
                                           ->line(local_line_neighbor)
                                           ->vertex_index(1))
                                  {
                                    // We have a flipped edge of first type
                                    mask |= line_to_edge[local_line][3];
                                    flipped = true;
                                  }
                                else
                                  ExcInternalError();

                                // Copy the unconstrained values
                                DoFCellAccessor<dim, dim, false>(
                                  &neighbor_cell->get_triangulation(),
                                  neighbor_cell->level(),
                                  neighbor_cell->index(),
                                  &cell->get_dof_handler())
                                  .get_dof_indices(neighbor_dofs_all);
                                // If the vector is distributed, we need to
                                // transform the global indices to local ones.
                                if (partitioner)
                                  for (auto &index : neighbor_dofs_all)
                                    index = partitioner->global_to_local(index);

                                for (unsigned int i = 0;
                                     i < neighbor_dofs_all_temp.size();
                                     ++i)
                                  neighbor_dofs_all_temp[i] =
                                    neighbor_dofs_all[lexicographic_mapping[i]];

                                for (unsigned int i = 0; i < n_dofs_1d; ++i)
                                  {
                                    // Get local dof index along line
                                    const unsigned int idx =
                                      line_dof_idx(local_line, i, n_dofs_1d);

                                    dof_indices[idx + idx_offset[comp]] =
                                      neighbor_dofs_all_temp
                                        [line_dof_idx(local_line_neighbor,
                                                      flipped ? fe_degree - i :
                                                                i,
                                                      n_dofs_1d) +
                                         idx_offset[comp]];
                                  }

                                // Stop looping over edge neighbors
                                break;
                              }
                          }
                      }
                  }
              }
            Assert(check(mask, dim), ExcInternalError());

            cell_has_hanging_node_constraints |=
              mask != ConstraintKinds::unconstrained;
          }
      return cell_has_hanging_node_constraints;
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
