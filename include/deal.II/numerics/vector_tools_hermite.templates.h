// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2023 by the deal.II authors
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

#ifndef vector_tools_hermite_template_h
#define vector_tools_hermite_template_h

#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_face.h>
#include <deal.II/fe/fe_hermite.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/sparse_direct.h> //?
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools_boundary.h>
#include <deal.II/numerics/vector_tools_hermite.h>
#include <deal.II/numerics/vector_tools_project.h>
#include <deal.II/numerics/vector_tools_rhs.h>

#include <string>


DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  // Internal namespace for implement Hermite boundary projection methods
  namespace internal
  {
    /*
     * Hermite elements position degrees of freedom corresponding to several
     * orders of derivatives as well as shape values along the boundary, which
     * means generally not all boundary d.o.f.s should be constrained. This
     * function determines which degrees of freedom on the boundary should be
     * constrained given the current boundary conditions, and flags the other
     * boundary d.o.f.s with a marker not to constrain.
     */
    template <int dim, int spacedim = dim, typename Number = double>
    void
    get_constrained_hermite_boundary_dofs(
      const DoFHandler<dim, spacedim> &dof_handler,
      const std::map<types::boundary_id, const Function<spacedim, Number> *>
                                           &boundary_functions,
      const unsigned int                    position,
      std::vector<types::global_dof_index> &dof_to_boundary,
      types::global_dof_index              &no_constrained_dofs)
    {
      AssertDimension(dof_handler.n_dofs(), dof_to_boundary.size());

      std::fill_n(dof_to_boundary.begin(),
                  dof_handler.n_dofs(),
                  numbers::unconstrained_boundary_dof_index);

      const unsigned int degree = dof_handler.get_fe().degree;
      const unsigned int regularity =
        dynamic_cast<const FE_Hermite<dim> &>(dof_handler.get_fe())
          .get_regularity();

      const unsigned int dofs_per_face = dof_handler.get_fe().n_dofs_per_face();
      const unsigned int constrained_dofs_per_face =
        dofs_per_face / (regularity + 1);
      AssertDimension(dofs_per_face,
                      (regularity + 1) * Utilities::pow(degree + 1, dim - 1));

      const unsigned int dofs_per_cell = Utilities::pow(degree + 1, dim);
      std::vector<types::global_dof_index> dofs_on_cell(dofs_per_cell);

      /*
       * Create a table of values corresponding to normal derivatives of the
       * chosen order, so that only those are constrained by the boundary method
       */
      Table<2, unsigned int> constrained_to_local_indices =
        dynamic_cast<const FE_Hermite<dim> &>(dof_handler.get_fe())
          .get_dofs_corresponding_to_outward_normal_derivatives(position);

      std::set<types::boundary_id> selected_boundary_components;

      for (auto i = boundary_functions.cbegin(); i != boundary_functions.cend();
           ++i)
        selected_boundary_components.insert(i->first);

      Assert(selected_boundary_components.find(
               numbers::internal_face_boundary_id) ==
               selected_boundary_components.end(),
             DoFTools::ExcInvalidBoundaryIndicator());

      /*
       * Index variable showing the position of the current constrained d.o.f
       * in the list of constrained d.o.f.s for the whole domain, and to keep
       * track of the current number of constrained d.o.f.s that have been
       * added to the list.
       */
      unsigned int global_constrained_dof_counter = 0;

      for (const auto &cell : dof_handler.active_cell_iterators())
        for (const unsigned int f : cell->face_indices())
          {
            AssertDimension(cell->get_fe().n_dofs_per_face(f), dofs_per_face);

            // Check if face is on selected boundary section
            if (selected_boundary_components.find(
                  cell->face(f)->boundary_id()) !=
                selected_boundary_components.end())
              {
                cell->get_dof_indices(dofs_on_cell);

                for (unsigned int i = 0; i < constrained_dofs_per_face; ++i)
                  {
                    const types::global_dof_index index =
                      dofs_on_cell[constrained_to_local_indices(f, i)];
                    Assert(index < dof_to_boundary.size(),
                           ExcDimensionMismatch(index, dof_to_boundary.size()));

                    if (dof_to_boundary[index] == numbers::invalid_dof_index)
                      dof_to_boundary[index] = global_constrained_dof_counter++;
                  }
              }
          }

      Assert((no_constrained_dofs !=
              dof_handler.n_boundary_dofs(boundary_functions)) ||
               (regularity == 0),
             ExcInternalError());
      no_constrained_dofs = global_constrained_dof_counter;
    }



    template <int dim, int spacedim = dim, typename Number = double>
    void
    do_hermite_direct_projection(
      const MappingCartesian<dim>     &mapping_h,
      const DoFHandler<dim, spacedim> &dof_handler,
      const std::map<types::boundary_id, const Function<spacedim, Number> *>
                                                &boundary_functions,
      const Quadrature<dim - 1>                 &quadrature,
      const unsigned int                         position,
      std::map<types::global_dof_index, Number> &boundary_values,
      std::vector<unsigned int>                  component_mapping = {})
    {
      // Return immediately if no constraint functions are provided
      if (boundary_functions.size() == 0)
        return;

      // For dim=1, the problem simplifies to interpolation at the boundaries
      if (dim == 1)
        {
          Assert(component_mapping.size() == 0, ExcNotImplemented());

          for (const auto &cell : dof_handler.active_cell_iterators())
            for (const unsigned int direction :
                 GeometryInfo<dim>::face_indices())
              if (cell->at_boundary(direction) &&
                  (boundary_functions.find(
                     cell->face(direction)->boundary_id()) !=
                   boundary_functions.end()))
                {
                  const Function<spacedim, Number> &current_function =
                    *boundary_functions
                       .find(cell->face(direction)->boundary_id())
                       ->second;

                  const FiniteElement<dim, spacedim> &fe_herm =
                    dof_handler.get_fe();

                  AssertDimension(fe_herm.n_components(),
                                  current_function.n_components);

                  Vector<Number> boundary_value(fe_herm.n_components());

                  if (boundary_value.size() == 1)
                    boundary_value(0) =
                      current_function.value(cell->vertex(direction));
                  else
                    current_function.vector_value(cell->vertex(direction),
                                                  boundary_value);

                  boundary_values[cell->vertex_dof_index(
                    direction, position, cell->active_fe_index())] =
                    boundary_value(
                      fe_herm.face_system_to_component_index(position).first);
                }

          return;
        }

      // dim=2 or higher, needs actual projection
      if (component_mapping.size() == 0)
        {
          AssertDimension(dof_handler.get_fe().n_components(),
                          boundary_functions.begin()->second->n_components);

          component_mapping.resize(dof_handler.get_fe().n_components());

          for (unsigned int i = 0; i < component_mapping.size(); ++i)
            component_mapping[i] = i;
        }
      else
        AssertDimension(component_mapping.size(),
                        dof_handler.get_fe().n_components());

      std::vector<types::global_dof_index> dof_to_boundary_mapping(
        dof_handler.n_dofs(), numbers::invalid_dof_index);
      types::global_dof_index next_boundary_index;

      get_constrained_hermite_boundary_dofs(dof_handler,
                                            boundary_functions,
                                            position,
                                            dof_to_boundary_mapping,
                                            next_boundary_index);

      // Check to see if the boundaries are homogeneous
      bool zero_indicator = true;
      for (const auto &func_pair : boundary_functions)
        {
          zero_indicator = zero_indicator &
                           (dynamic_cast<const Functions::ZeroFunction<dim> *>(
                              func_pair.second) != nullptr);
          if (!zero_indicator)
            break;
        }

      if (zero_indicator)
        {
          for (unsigned int i = 0; i < dof_to_boundary_mapping.size(); ++i)
            if (dof_to_boundary_mapping[i] != numbers::invalid_dof_index)
              boundary_values[i] = 0;

          return;
        }

      SparsityPattern sparsity;

      {
        DynamicSparsityPattern dsp(next_boundary_index, next_boundary_index);
        DoFTools::make_boundary_sparsity_pattern(dof_handler,
                                                 boundary_functions,
                                                 dof_to_boundary_mapping,
                                                 dsp);

        sparsity.copy_from(dsp);
      }

      // Assert mesh is not partially refined
      // TODO: Add functionality on partially refined meshes
      int level = -1;

      for (const auto &cell : dof_handler.active_cell_iterators())
        for (auto f : cell->face_indices())
          if (cell->at_boundary(f))
            {
              if (level == -1)
                level = cell->level();
              else
                Assert(level == cell->level(),
                       ExcMessage("The mesh you use in projecting "
                                  "boundary values has hanging nodes"
                                  " at the boundary. This would require "
                                  "dealing with hanging node constraints"
                                  " when solving the linear system on "
                                  "the boundary, but this is not "
                                  "currently implemented."));
            }

      // make mass matrix and right hand side
      SparseMatrix<Number> mass_matrix(sparsity);
      Vector<Number>       rhs(sparsity.n_rows());
      Vector<Number>       boundary_projection(rhs.size());

      MatrixCreator::create_boundary_mass_matrix(
        mapping_h,
        dof_handler,
        quadrature,
        mass_matrix,
        boundary_functions,
        rhs,
        dof_to_boundary_mapping,
        static_cast<const Function<spacedim, Number> *>(nullptr),
        component_mapping);

      {
        /*
         * Current implementation uses a direct solver to better account for
         * the mass matrix being potentially sufficiently ill-conditioned to
         * prevent iterations from converging.
         * TODO: Find a good preconditioner to allow more efficient iterative
         * methods
         */
        SparseDirectUMFPACK mass_inv;
        mass_inv.initialize(mass_matrix);
        mass_inv.vmult(boundary_projection, rhs);
      }

      // fill in boundary values
      for (unsigned int i = 0; i < dof_to_boundary_mapping.size(); ++i)
        if (dof_to_boundary_mapping[i] != numbers::invalid_dof_index)
          {
            AssertIsFinite(boundary_projection(dof_to_boundary_mapping[i]));

            /*
             * This is where the boundary solution, which only considers the
             * DOFs that are actually constrained, is applied to the boundary
             * DOFs in the global solution. Remember: i is the global dof
             * number, dof_to_boundary_mapping[i] is the number on the boundary
             * and thus in the solution vector.
             */
            boundary_values[i] =
              boundary_projection(dof_to_boundary_mapping[i]);
          }
      return;
    }
  } // namespace internal



  template <int dim, int spacedim, typename Number>
  void
  project_hermite_boundary_values(
    const Mapping<dim, spacedim>    &mapping_h,
    const DoFHandler<dim, spacedim> &dof_handler,
    const std::map<types::boundary_id, const Function<spacedim, Number> *>
                                              &boundary_functions,
    const Quadrature<dim - 1>                 &quadrature,
    const unsigned int                         boundary_norm_deriv_order,
    std::map<types::global_dof_index, Number> &boundary_values,
    std::vector<unsigned int>                  component_mapping)
  {
    /*
     * This version implements projected values directly,
     * so it's necessary to check that this is possible
     */

    Assert((2 * boundary_norm_deriv_order < dof_handler.get_fe().degree),
           ExcMessage("To constrain r order derivatives on the boundary "
                      "requires a Hermite FE object of at least order "
                      "2r+1."));
    Assert((dynamic_cast<const MappingCartesian<dim> *>(&mapping_h) != nullptr),
           ExcMessage("project_hermite_boundary_values should only be"
                      " used with a mapping compatible with Hermite "
                      "finite element methods."));
    Assert((dynamic_cast<const FE_Hermite<dim> *>(&dof_handler.get_fe()) !=
            nullptr),
           ExcMessage("project_hermite_boundary_values should only"
                      " be used with a DOF handler linked to a "
                      "FE_Hermite object."));

    internal::do_hermite_direct_projection(
      *static_cast<const MappingCartesian<dim> *>(&mapping_h),
      dof_handler,
      boundary_functions,
      quadrature,
      boundary_norm_deriv_order,
      boundary_values,
      component_mapping);
  }



  template <int dim, int spacedim, typename Number>
  void
  project_hermite_boundary_values(
    const Mapping<dim, spacedim>    &mapping_h,
    const DoFHandler<dim, spacedim> &dof_handler,
    const std::map<types::boundary_id, const Function<spacedim, Number> *>
                                              &boundary_functions,
    const Quadrature<dim - 1>                 &quadrature,
    std::map<types::global_dof_index, Number> &boundary_values,
    std::vector<unsigned int>                  component_mapping)
  {
    Assert((dynamic_cast<const MappingCartesian<dim> *>(&mapping_h) != nullptr),
           ExcMessage("project_hermite_boundary_values should only be"
                      " used with a mapping compatible with Hermite "
                      "finite element methods."));
    Assert((dynamic_cast<const FE_Hermite<dim> *>(&dof_handler.get_fe()) !=
            nullptr),
           ExcMessage("project_hermite_boundary_values should only"
                      " be used with a DOF handler linked to a "
                      "FE_Hermite object."));

    internal::do_hermite_direct_projection(
      *static_cast<const MappingCartesian<dim> *>(&mapping_h),
      dof_handler,
      boundary_functions,
      quadrature,
      0,
      boundary_values,
      component_mapping);
  }



  template <int dim, typename VectorType, int spacedim>
  void
  project_hermite(
    const Mapping<dim, spacedim>                              &mapping_h,
    const DoFHandler<dim, spacedim>                           &dof,
    const AffineConstraints<typename VectorType::value_type>  &constraints,
    const Quadrature<dim>                                     &quadrature,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const bool                 enforce_zero_boundary,
    const Quadrature<dim - 1> &q_boundary,
    const bool                 project_to_boundary_first)
  {
    using number = typename VectorType::value_type;

    Assert((dynamic_cast<const MappingCartesian<dim> *>(&mapping_h) != nullptr),
           ExcMessage(
             "The function project_hermite() can only be used with a mapping "
             "compatible with Hermite finite element methods."));

    Assert((dynamic_cast<const FE_Hermite<dim, spacedim> *>(&dof.get_fe()) !=
            nullptr),
           ExcMessage("The project_hermite function should only be used"
                      " with a DoF handler associated with an FE_Hermite"
                      " object with the same dim and spacedim as the "
                      "mapping."));
    Assert(dof.get_fe(0).n_components() == function.n_components,
           ExcDimensionMismatch(dof.get_fe(0).n_components(),
                                function.n_components));
    Assert(vec.size() == dof.n_dofs(),
           ExcDimensionMismatch(vec.size(), dof.n_dofs()));

    // make up boundary values
    std::map<types::global_dof_index, number> boundary_values;
    const std::vector<types::boundary_id>     active_boundary_ids =
      dof.get_triangulation().get_boundary_ids();

    if (enforce_zero_boundary)
      {
        std::map<types::boundary_id, const Function<spacedim, number> *>
          boundary;

        for (const auto it : active_boundary_ids)
          boundary.emplace(std::make_pair(it, nullptr));

        const std::map<types::boundary_id, const Function<spacedim, number> *>
          new_boundary =
            static_cast<const std::map<types::boundary_id,
                                       const Function<spacedim, number> *>>(
              boundary);

        std::vector<types::global_dof_index> dof_to_boundary(
          dof.n_dofs(), numbers::invalid_dof_index);
        types::global_dof_index end_boundary_dof = 0;

        internal::get_constrained_hermite_boundary_dofs(
          dof, new_boundary, 0, dof_to_boundary, end_boundary_dof);

        if (end_boundary_dof != 0)
          for (types::global_dof_index i = 0; i < dof.n_dofs(); ++i)
            if (dof_to_boundary[i] != numbers::invalid_dof_index)
              boundary_values.emplace(std::make_pair(i, 0));
      }
    else if (project_to_boundary_first)
      {
        std::map<types::boundary_id, const Function<spacedim, number> *>
          boundary_function;

        for (const auto it : active_boundary_ids)
          boundary_function.emplace(std::make_pair(it, &function));

        const std::map<types::boundary_id, const Function<spacedim, number> *>
          new_boundary_2 =
            static_cast<const std::map<types::boundary_id,
                                       const Function<spacedim, number> *>>(
              boundary_function);

        /* project_hermite_boundary_values<dim, spacedim, number>(
          mapping_h,
          dof,
          new_boundary_2,
          q_boundary,
          boundary_values); */
        internal::do_hermite_direct_projection(
          *static_cast<const MappingCartesian<dim> *>(&mapping_h),
          dof,
          new_boundary_2,
          q_boundary,
          0,
          boundary_values);
      }

    // check if constraints are compatible (see below)
    bool constraints_are_compatible = true;

    for (const auto &value : boundary_values)
      if (constraints.is_constrained(value.first))
        if ((constraints.get_constraint_entries(value.first)->size() > 0) &&
            (constraints.get_inhomogeneity(value.first) != value.second))
          constraints_are_compatible = false;

    // set up mass matrix and right hand side
    Vector<number>  vec_result(dof.n_dofs());
    SparsityPattern sparsity;

    {
      DynamicSparsityPattern dsp(dof.n_dofs(), dof.n_dofs());
      DoFTools::make_sparsity_pattern(dof,
                                      dsp,
                                      constraints,
                                      !constraints_are_compatible);

      sparsity.copy_from(dsp);
    }

    SparseMatrix<number> mass_matrix(sparsity);
    Vector<number>       tmp(mass_matrix.n());

    /*
     * If the constraints object does not conflict with the given boundary
     * values (i.e., it either does not contain boundary values or it contains
     * the same as boundary_values), we can let it call
     * distribute_local_to_global straight away, otherwise we need to first
     * interpolate the boundary values and then condense the matrix and vector
     */
    if (constraints_are_compatible)
      {
        const Function<spacedim, number> *dummy = nullptr;

        MatrixCreator::create_mass_matrix(mapping_h,
                                          dof,
                                          quadrature,
                                          mass_matrix,
                                          function,
                                          tmp,
                                          dummy,
                                          constraints);

        if (boundary_values.size() > 0)
          MatrixTools::apply_boundary_values(
            boundary_values, mass_matrix, vec_result, tmp, true);
      }
    else
      {
        // create mass matrix and rhs at once, which is faster.
        MatrixCreator::create_mass_matrix(
          mapping_h, dof, quadrature, mass_matrix, function, tmp);
        MatrixTools::apply_boundary_values(
          boundary_values, mass_matrix, vec_result, tmp, true);

        constraints.condense(mass_matrix, tmp);
      }

    /*
     * Current implementation uses a direct solver to better deal
     * with potentially ill-conditioned mass matrices.
     * TODO: Find a good preconditioner for the Hermite mass matrix
     */
    SparseDirectUMFPACK mass_inv;
    mass_inv.initialize(mass_matrix);
    mass_inv.vmult(vec_result, tmp);

    constraints.distribute(vec_result);

    // copy vec_result into vec. we can't use vec itself above, since
    // it may be of another type than Vector<double> and that wouldn't
    // necessarily go together with the matrix and other functions
    for (unsigned int i = 0; i < vec.size(); ++i)
      ::dealii::internal::ElementAccess<VectorType>::set(vec_result(i), i, vec);
  }
} // namespace VectorTools



DEAL_II_NAMESPACE_CLOSE

#endif
