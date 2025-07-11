// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_hanging_nodes_internal_h
#define dealii_hanging_nodes_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/ndarray.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_iso_q1.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/hp/fe_collection.h>

#include <boost/container/small_vector.hpp>


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
     * constrained edge in that direction (only valid in 3d).
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
            return false; // in 2d there are no edge constraints

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
    DEAL_II_HOST_DEVICE inline ConstraintKinds
    operator|(const ConstraintKinds f1, const ConstraintKinds f2)
    {
      return static_cast<ConstraintKinds>(static_cast<std::uint16_t>(f1) |
                                          static_cast<std::uint16_t>(f2));
    }



    /**
     * Global operator which sets the bits from the second argument also in the
     * first one.
     */
    DEAL_II_HOST_DEVICE inline ConstraintKinds &
    operator|=(ConstraintKinds &f1, const ConstraintKinds f2)
    {
      f1 = f1 | f2;
      return f1;
    }



    /**
     * Global operator which checks inequality.
     */
    DEAL_II_HOST_DEVICE inline bool
    operator!=(const ConstraintKinds f1, const ConstraintKinds f2)
    {
      return static_cast<std::uint16_t>(f1) != static_cast<std::uint16_t>(f2);
    }



    /**
     * Global operator which checks if the first argument is less than the
     * second.
     */
    DEAL_II_HOST_DEVICE inline bool
    operator<(const ConstraintKinds f1, const ConstraintKinds f2)
    {
      return static_cast<std::uint16_t>(f1) < static_cast<std::uint16_t>(f2);
    }



    /**
     * Global operator which performs a binary and for the provided arguments.
     */
    DEAL_II_HOST_DEVICE inline ConstraintKinds
    operator&(const ConstraintKinds f1, const ConstraintKinds f2)
    {
      return static_cast<ConstraintKinds>(static_cast<std::uint16_t>(f1) &
                                          static_cast<std::uint16_t>(f2));
    }



    /**
     * This class creates the mask used in the treatment of hanging nodes in
     * Portable::MatrixFree.
     * The implementation of this class is explained in detail in
     * @cite munch2022hn.
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
       * Return the memory consumption of the allocated memory in this class.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Compute the value of the constraint mask for a given cell.
       */
      template <typename CellIterator>
      bool
      setup_constraints(
        const CellIterator                                       &cell,
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
        const std::vector<std::vector<unsigned int>> &lexicographic_mapping,
        std::vector<types::global_dof_index>         &dof_indices,
        const ArrayView<ConstraintKinds>             &mask) const;

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
        const CellIterator                                       &cell,
        const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
        const std::vector<std::vector<unsigned int>> &lexicographic_mapping,
        const std::vector<std::vector<bool>>         &component_mask,
        const ConstraintKinds                        &refinement_configuration,
        std::vector<types::global_dof_index>         &dof_indices) const;

    private:
      /**
       * Set up line-to-cell mapping for edge constraints in 3d.
       */
      void
      setup_line_to_cell(const Triangulation<dim> &triangulation);

      void
      rotate_subface_index(int times, unsigned int &subface_index) const;

      void
      orient_face(const types::geometric_orientation    combined_orientation,
                  const unsigned int                    n_dofs_1d,
                  std::vector<types::global_dof_index> &dofs) const;

      unsigned int
      line_dof_idx(int          local_line,
                   unsigned int dof,
                   unsigned int n_dofs_1d) const;

      void
      transpose_subface_index(unsigned int &subface) const;

      std::vector<
        boost::container::small_vector<std::array<unsigned int, 3>, 6>>
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
    inline std::size_t
    HangingNodes<dim>::memory_consumption() const
    {
      std::size_t size = 0;
      for (const auto &a : line_to_cells)
        size +=
          (a.capacity() > 6 ? a.capacity() : 0) * sizeof(a[0]) + sizeof(a);
      return size;
    }



    template <int dim>
    inline void
    HangingNodes<dim>::setup_line_to_cell(
      const Triangulation<dim> & /*triangulation*/)
    {}



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

      // In 3d, we can have DoFs on only an edge being constrained (e.g. in a
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

      std::vector<
        boost::container::small_vector<std::array<unsigned int, 3>, 6>>
        line_to_inactive_cells(n_raw_lines);

      // First add active and inactive cells to their lines:
      for (const auto &cell : triangulation.cell_iterators())
        {
          const unsigned int cell_level = cell->level();
          const unsigned int cell_index = cell->index();
          for (unsigned int line = 0; line < GeometryInfo<3>::lines_per_cell;
               ++line)
            {
              const unsigned int line_idx = cell->line_index(line);
              if (cell->is_active())
                line_to_cells[line_idx].push_back(
                  {{cell_level, cell_index, line}});
              else
                line_to_inactive_cells[line_idx].push_back(
                  {{cell_level, cell_index, line}});
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
              const Triangulation<3>::cell_iterator inactive_cell(
                &triangulation,
                line_to_inactive_cells[line_idx][0][0],
                line_to_inactive_cells[line_idx][0][1]);
              const unsigned int neighbor_line =
                line_to_inactive_cells[line_idx][0][2];

              for (unsigned int c = 0; c < 2; ++c)
                {
                  const auto &child =
                    inactive_cell->child(line_to_children[neighbor_line][c]);
                  const unsigned int child_line_idx =
                    child->line_index(neighbor_line);

                  Assert(child->is_active(), ExcInternalError());

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
              if (dim == 1 ||
                  (dynamic_cast<const FE_Q<dim> *>(
                     &fe_collection[i].base_element(base_element_index)) ==
                     nullptr &&
                   dynamic_cast<const FE_Q_iso_Q1<dim> *>(
                     &fe_collection[i].base_element(base_element_index)) ==
                     nullptr))
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
      if ((dim == 3 && line_to_cells.empty()) ||
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

          // Ignore if the neighbors are FE_Nothing
          if (neighbor->get_fe().n_dofs_per_cell() == 0)
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

              const unsigned int line_index = cell->line_index(line_no);

              const auto edge_neighbor =
                std::find_if(line_to_cells[line_index].begin(),
                             line_to_cells[line_index].end(),
                             [&cell](const auto &edge_neighbor) {
                               DoFCellAccessor<dim, dim, false> dof_cell(
                                 &cell->get_triangulation(),
                                 edge_neighbor[0],
                                 edge_neighbor[1],
                                 &cell->get_dof_handler());
                               return dof_cell.is_artificial() == false &&
                                      dof_cell.level() < cell->level() &&
                                      dof_cell.get_fe().n_dofs_per_cell() > 0;
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
      const CellIterator                                       &cell,
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
      const std::vector<std::vector<unsigned int>> &lexicographic_mapping,
      const std::vector<std::vector<bool>>         &supported_components,
      const ConstraintKinds                        &refinement_configuration,
      std::vector<types::global_dof_index>         &dof_indices) const
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
           base_element_index < fe.n_base_elements();
           ++base_element_index)
        for (unsigned int c = 0;
             c < fe.element_multiplicity(base_element_index);
             ++c)
          idx_offset.push_back(
            idx_offset.back() +
            fe.base_element(base_element_index).n_dofs_per_cell());

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
                DEAL_II_NOT_IMPLEMENTED();
            }

        DEAL_II_NOT_IMPLEMENTED();

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
                 base_element_index < fe.n_base_elements();
                 ++base_element_index)
              for (unsigned int c = 0;
                   c < fe.element_multiplicity(base_element_index);
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

                  // fix DoFs depending on orientation, rotation, and flip
                  if (dim == 2)
                    {
                      // TODO: this needs to be implemented for simplices but
                      // all-quad meshes are OK
                      Assert(cell->combined_face_orientation(face_no) ==
                               numbers::default_geometric_orientation,
                             ExcNotImplemented());
                    }
                  else if (dim == 3)
                    {
                      orient_face(cell->combined_face_orientation(face_no),
                                  n_dofs_1d,
                                  neighbor_dofs);
                    }
                  else
                    {
                      DEAL_II_NOT_IMPLEMENTED();
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
                             [&cell](const auto &edge_array) {
                               const typename Triangulation<dim>::cell_iterator
                                 edge_neighbor(&cell->get_triangulation(),
                                               edge_array[0],
                                               edge_array[1]);
                               return edge_neighbor->is_artificial() == false &&
                                      edge_neighbor->level() < cell->level();
                             });

              if (edge_neighbor == line_to_cells[line_index].end())
                continue;

              const DoFCellAccessor<dim, dim, false> neighbor_cell(
                &cell->get_triangulation(),
                (*edge_neighbor)[0],
                (*edge_neighbor)[1],
                &cell->get_dof_handler());
              const auto local_line_neighbor = (*edge_neighbor)[2];

              neighbor_cell.get_dof_indices(neighbor_dofs_all);

              if (partitioner)
                for (auto &index : neighbor_dofs_all)
                  index = partitioner->global_to_local(index);

              for (unsigned int i = 0; i < neighbor_dofs_all_temp.size(); ++i)
                neighbor_dofs_all_temp[i] = neighbor_dofs_all
                  [lexicographic_mapping[cell->active_fe_index()][i]];

              const bool flipped =
                cell->line_orientation(line_no) !=
                neighbor_cell.line_orientation(local_line_neighbor);

              for (unsigned int base_element_index = 0, comp = 0;
                   base_element_index < fe.n_base_elements();
                   ++base_element_index)
                for (unsigned int c = 0;
                     c < fe.element_multiplicity(base_element_index);
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
      const CellIterator                                       &cell,
      const std::shared_ptr<const Utilities::MPI::Partitioner> &partitioner,
      const std::vector<std::vector<unsigned int>> &lexicographic_mapping,
      std::vector<types::global_dof_index>         &dof_indices,
      const ArrayView<ConstraintKinds>             &masks) const
    {
      // 1) check if finite elements support fast hanging-node algorithm
      const auto supported_components = compute_supported_components(
        cell->get_dof_handler().get_fe_collection());

      if (std::none_of(supported_components.begin(),
                       supported_components.end(),
                       [](const auto &a) {
                         return *std::max_element(a.begin(), a.end());
                       }))
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
    HangingNodes<dim>::orient_face(
      const types::geometric_orientation    combined_orientation,
      const unsigned int                    n_dofs_1d,
      std::vector<types::global_dof_index> &dofs) const
    {
      const auto [orientation, rotation, flip] =
        ::dealii::internal::split_face_orientation(combined_orientation);
      const int n_rotations =
        rotation || flip ? 4 - int(rotation) - 2 * int(flip) : 0;

      const unsigned int rot_mapping[4] = {2, 0, 3, 1};
      Assert(n_dofs_1d > 1, ExcInternalError());
      // 'per line' has the same meaning here as in FiniteElementData, i.e., the
      // number of dofs assigned to a line (not including vertices).
      const unsigned int dofs_per_line = n_dofs_1d - 2;

      // rotate:
      std::vector<types::global_dof_index> copy(dofs.size());
      for (int t = 0; t < n_rotations; ++t)
        {
          std::swap(copy, dofs);

          // Vertices
          for (unsigned int i = 0; i < 4; ++i)
            dofs[rot_mapping[i]] = copy[i];

          // Edges
          unsigned int offset = 4;
          for (unsigned int i = 0; i < dofs_per_line; ++i)
            {
              // Left edge
              dofs[offset + i] =
                copy[offset + 2 * dofs_per_line + (dofs_per_line - 1 - i)];
              // Right edge
              dofs[offset + dofs_per_line + i] =
                copy[offset + 3 * dofs_per_line + (dofs_per_line - 1 - i)];
              // Bottom edge
              dofs[offset + 2 * dofs_per_line + i] =
                copy[offset + dofs_per_line + i];
              // Top edge
              dofs[offset + 3 * dofs_per_line + i] = copy[offset + i];
            }

          // Interior points
          offset += 4 * dofs_per_line;

          for (unsigned int i = 0; i < dofs_per_line; ++i)
            for (unsigned int j = 0; j < dofs_per_line; ++j)
              dofs[offset + i * dofs_per_line + j] =
                copy[offset + j * dofs_per_line + (dofs_per_line - 1 - i)];
        }

      // transpose (note that we are using the standard geometric orientation
      // here so orientation = true is the default):
      if (!orientation)
        {
          copy = dofs;

          // Vertices
          dofs[1] = copy[2];
          dofs[2] = copy[1];

          // Edges
          unsigned int offset = 4;
          for (unsigned int i = 0; i < dofs_per_line; ++i)
            {
              // Right edge
              dofs[offset + i] = copy[offset + 2 * dofs_per_line + i];
              // Left edge
              dofs[offset + dofs_per_line + i] =
                copy[offset + 3 * dofs_per_line + i];
              // Bottom edge
              dofs[offset + 2 * dofs_per_line + i] = copy[offset + i];
              // Top edge
              dofs[offset + 3 * dofs_per_line + i] =
                copy[offset + dofs_per_line + i];
            }

          // Interior
          offset += 4 * dofs_per_line;
          for (unsigned int i = 0; i < dofs_per_line; ++i)
            for (unsigned int j = 0; j < dofs_per_line; ++j)
              dofs[offset + i * dofs_per_line + j] =
                copy[offset + j * dofs_per_line + i];
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
