// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors

#include <deal.II/base/mpi.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer_matrix_free.h>

#include <algorithm>
#include <cmath>
#include <vector>

using namespace dealii;

namespace RotatedPeriodicityTest
{
  using Vector = LinearAlgebra::distributed::Vector<double>;

  struct Errors
  {
    double constraint   = 0.;
    double trace        = 0.;
    double prolongation = 0.;
    double adjoint      = 0.;
  };

  template <int dim>
  double
  exact_value(const Point<dim>  &point,
              const unsigned int component,
              const unsigned int components,
              const unsigned int angle)
  {
    if (components == 1 || component == dim)
      return 2.;
    if (angle == 0)
      return component == 0 ? 0. : 1. + point[component];
    if constexpr (dim == 3)
      if (component == 2)
        return 1. + point[2];
    return component == 0 ? -point[1] : point[0];
  }

  template <int dim>
  void
  fill_exact(const DoFHandler<dim> &dofs,
             const unsigned int     level,
             const unsigned int     angle,
             Vector                &values)
  {
    MappingQ1<dim> mapping;
    const auto    &fe     = dofs.get_fe();
    const auto    &points = fe.get_unit_support_points();
    const auto    &owned  = dofs.locally_owned_mg_dofs(level);
    std::vector<types::global_dof_index> indices(fe.n_dofs_per_cell());
    for (const auto &cell : dofs.mg_cell_iterators_on_level(level))
      if (!cell->is_artificial_on_level())
        {
          cell->get_mg_dof_indices(indices);
          for (unsigned int i = 0; i < indices.size(); ++i)
            if (owned.is_element(indices[i]))
              values[indices[i]] =
                exact_value(mapping.transform_unit_to_real_cell(cell,
                                                                points[i]),
                            fe.system_to_component_index(i).first,
                            fe.n_components(),
                            angle);
        }
    values.compress(VectorOperation::insert);
    values.update_ghost_values();
  }

  template <typename TriangulationType>
  double
  check_trace(const DoFHandler<TriangulationType::dimension> &dofs,
              const TriangulationType                        &tria,
              const unsigned int                              level,
              const unsigned int                              angle,
              const AffineConstraints<double>                &constraints,
              const Vector                                   &layout)
  {
    constexpr unsigned int dim = TriangulationType::dimension;
    const auto            &fe  = dofs.get_fe();
    Vector                 values;
    values.reinit(layout);
    for (const auto i : dofs.locally_owned_mg_dofs(level))
      values[i] = std::sin(0.17 * (i + 1));
    values.compress(VectorOperation::insert);
    values.update_ghost_values();
    constraints.distribute(values);
    values.update_ghost_values();
    const Vector         &read_values = values;
    const QGauss<dim - 1> quadrature(fe.degree + 1);
    FEFaceValues<dim>     first_values(fe,
                                   quadrature,
                                   update_values | update_quadrature_points);
    FEFaceValues<dim>     second_values(fe,
                                    quadrature,
                                    update_values | update_quadrature_points);
    const double          cosine = std::cos(angle * numbers::PI / 180.);
    const double          sine   = std::sin(angle * numbers::PI / 180.);
    double                error  = 0.;
    for (const auto &[first, second] : tria.get_periodic_face_map())
      {
        if (first.first->is_artificial_on_level() ||
            second.first.first->is_artificial_on_level() ||
            first.first->level() != static_cast<int>(level) ||
            second.first.first->level() != static_cast<int>(level) ||
            first.first->face(first.second)->boundary_id() !=
              (angle == 0 ? 0 : 2) ||
            second.first.first->face(second.first.second)->boundary_id() !=
              (angle == 0 ? 1 : 3))
          continue;
        const auto cell1 = first.first->as_dof_handler_level_iterator(dofs);
        const auto cell2 =
          second.first.first->as_dof_handler_level_iterator(dofs);
        first_values.reinit(first.first, first.second);
        second_values.reinit(second.first.first, second.first.second);
        std::vector<types::global_dof_index> indices1(fe.n_dofs_per_cell()),
          indices2(fe.n_dofs_per_cell());
        cell1->get_mg_dof_indices(indices1);
        cell2->get_mg_dof_indices(indices2);
        std::vector<std::vector<double>> v1(
          quadrature.size(), std::vector<double>(fe.n_components()));
        auto v2 = v1;
        for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
          for (unsigned int q = 0; q < quadrature.size(); ++q)
            {
              const auto component = fe.system_to_component_index(i).first;
              v1[q][component] +=
                first_values.shape_value(i, q) * read_values[indices1[i]];
              v2[q][component] +=
                second_values.shape_value(i, q) * read_values[indices2[i]];
            }
        for (unsigned int q = 0; q < quadrature.size(); ++q)
          {
            const auto point  = first_values.quadrature_point(q);
            Point<dim> mapped = point;
            if (angle == 0)
              mapped[0] += 1.;
            else
              {
                mapped[0] = cosine * point[0] - sine * point[1];
                mapped[1] = sine * point[0] + cosine * point[1];
              }
            unsigned int other = numbers::invalid_unsigned_int;
            for (unsigned int j = 0; j < quadrature.size(); ++j)
              if (mapped.distance(second_values.quadrature_point(j)) < 1.e-12)
                other = j;
            AssertThrow(other != numbers::invalid_unsigned_int,
                        ExcMessage("Periodic quadrature points do not match."));
            for (unsigned int component = 0; component < fe.n_components();
                 ++component)
              {
                double expected = v1[q][component];
                if (angle != 0 && fe.n_components() > 1 && component < 2)
                  expected = component == 0 ?
                               cosine * v1[q][0] - sine * v1[q][1] :
                               sine * v1[q][0] + cosine * v1[q][1];
                error =
                  std::max(error, std::abs(v2[other][component] - expected));
              }
          }
      }
    return error;
  }

  template <typename TriangulationType>
  Errors
  check_hierarchy(TriangulationType &tria,
                  const MPI_Comm     communicator,
                  const unsigned int angle,
                  const unsigned int degree,
                  const unsigned int components)
  {
    constexpr unsigned int dim = TriangulationType::dimension;
    FESystem<dim>          fe(FE_Q<dim>(degree), components);
    DoFHandler<dim>        dofs(tria);
    dofs.distribute_dofs(fe);
    dofs.distribute_mg_dofs();
    MGConstrainedDoFs mg;
    mg.initialize(dofs, MGLevelObject<IndexSet>(), false);
    const unsigned int                       levels = tria.n_global_levels();
    MGLevelObject<AffineConstraints<double>> constraints(0, levels - 1);
    std::vector<Vector>                      exact(levels);
    FullMatrix<double>                       rotation(dim);
    rotation             = IdentityMatrix(dim);
    const double radians = angle * numbers::PI / 180.;
    rotation[0][0] = rotation[1][1] = std::cos(radians);
    rotation[0][1]                  = -std::sin(radians);
    rotation[1][0]                  = std::sin(radians);
    const types::boundary_id left   = angle == 0 ? 0 : 2;
    const types::boundary_id right  = angle == 0 ? 1 : 3;
    Errors                   errors;
    for (unsigned int level = 0; level < levels; ++level)
      {
        const auto relevant =
          DoFTools::extract_locally_relevant_level_dofs(dofs, level);
        const auto &owned = dofs.locally_owned_mg_dofs(level);
        constraints[level].reinit(owned, relevant);
        for (const auto &[first, second] : tria.get_periodic_face_map())
          {
            if (first.first->is_artificial_on_level() ||
                second.first.first->is_artificial_on_level() ||
                first.first->level() != static_cast<int>(level) ||
                second.first.first->level() != static_cast<int>(level))
              continue;
            const auto b1 = first.first->face(first.second)->boundary_id();
            const auto b2 =
              second.first.first->face(second.first.second)->boundary_id();
            if (b1 != left || b2 != right)
              continue;
            const auto face1 =
              first.first->as_dof_handler_level_iterator(dofs)->face(
                first.second);
            const auto face2 =
              second.first.first->as_dof_handler_level_iterator(dofs)->face(
                second.first.second);
            if (angle == 0 || components == 1)
              DoFTools::make_periodicity_constraints_on_level(
                face1,
                face2,
                level,
                constraints[level],
                ComponentMask(),
                second.second);
            else
              DoFTools::make_periodicity_constraints_on_level(
                face1,
                face2,
                level,
                constraints[level],
                ComponentMask(),
                second.second,
                rotation,
                {0});
          }
        constraints[level].close();
        mg.add_user_constraints(level, constraints[level]);
        AssertThrow(Utilities::MPI::sum(constraints[level].n_constraints(),
                                        communicator) > 0,
                    ExcMessage("The level has no periodic constraints."));
        exact[level].reinit(owned, relevant, communicator);
        fill_exact(dofs, level, angle, exact[level]);
        // Test an arbitrary constrained vector on the actual paired faces.
        // A preserved exact field alone would not detect missing constraints.
        errors.trace = std::max(
          errors.trace,
          check_trace(
            dofs, tria, level, angle, constraints[level], exact[level]));
        for (const auto &line : constraints[level].get_lines())
          if (owned.is_element(line.index))
            {
              double residual = exact[level][line.index] - line.inhomogeneity;
              for (const auto &entry : line.entries)
                residual -= entry.second * exact[level][entry.first];
              errors.constraint =
                std::max(errors.constraint, std::abs(residual));
            }
      }

    MGTransferMatrixFree<dim, double> transfer(mg);
    transfer.build(dofs);
    for (unsigned int level = 1; level < levels; ++level)
      {
        Vector prolonged;
        prolonged.reinit(exact[level]);
        transfer.prolongate(level, prolonged, exact[level - 1]);
        constraints[level].distribute(prolonged);
        prolonged -= exact[level];
        errors.prolongation =
          std::max(errors.prolongation, prolonged.linfty_norm());

        // On the independent DoFs, restriction must be the transpose of
        // prolongation. This checks the reverse transfer without reusing its
        // implementation as a reference.
        Vector coarse, fine, pc, rf;
        coarse.reinit(exact[level - 1]);
        fine.reinit(exact[level]);
        pc.reinit(fine);
        rf.reinit(coarse);
        const auto &coarse_owned = dofs.locally_owned_mg_dofs(level - 1);
        const auto &fine_owned   = dofs.locally_owned_mg_dofs(level);
        for (const auto i : coarse_owned)
          coarse[i] = constraints[level - 1].is_constrained(i) ?
                        0. :
                        std::sin(0.3 * (i + 1));
        for (const auto i : fine_owned)
          fine[i] =
            constraints[level].is_constrained(i) ? 0. : std::cos(0.2 * (i + 1));
        transfer.prolongate(level, pc, coarse);
        transfer.restrict_and_add(level, rf, fine);
        const double lhs = pc * fine;
        const double rhs = coarse * rf;
        errors.adjoint =
          std::max(errors.adjoint,
                   std::abs(lhs - rhs) /
                     std::max({1., std::abs(lhs), std::abs(rhs)}));
      }
    errors.constraint = Utilities::MPI::max(errors.constraint, communicator);
    errors.trace      = Utilities::MPI::max(errors.trace, communicator);
    errors.prolongation =
      Utilities::MPI::max(errors.prolongation, communicator);
    errors.adjoint = Utilities::MPI::max(errors.adjoint, communicator);
    return errors;
  }

  template <typename TriangulationType>
  std::vector<Errors>
  check_mesh(TriangulationType &tria,
             const MPI_Comm     communicator,
             const unsigned int angle,
             const unsigned int degree,
             const unsigned int components)
  {
    constexpr unsigned int dim = TriangulationType::dimension;
    if (angle == 0)
      GridGenerator::hyper_cube(tria, 0., 1., true);
    else
      {
        Point<dim> lower, upper;
        lower[0] = 1.;
        upper[0] = 2.;
        upper[1] = angle * numbers::PI / 180.;
        if constexpr (dim == 3)
          upper[2] = 1.;
        std::vector<unsigned int> repetitions(dim, 1);
        repetitions[1] = std::max(2U, angle / 22U);
        GridGenerator::subdivided_hyper_rectangle(
          tria, repetitions, lower, upper, true);
        GridTools::transform(
          [](const Point<dim> &point) {
            Point<dim> mapped = point;
            mapped[0]         = point[0] * std::cos(point[1]);
            mapped[1]         = point[0] * std::sin(point[1]);
            return mapped;
          },
          tria);
      }

    // Straight-sided nested cells make Cartesian linear fields exactly
    // representable on all levels; curved geometry is tested in ASPECT.
    tria.reset_all_manifolds();
    for (const auto &cell : tria.active_cell_iterators())
      cell->set_all_manifold_ids(numbers::flat_manifold_id);
    FullMatrix<double> rotation;
    if (angle != 0)
      {
        rotation.reinit(dim, dim);
        rotation             = IdentityMatrix(dim);
        const double radians = angle * numbers::PI / 180.;
        rotation[0][0] = rotation[1][1] = std::cos(radians);
        rotation[0][1]                  = -std::sin(radians);
        rotation[1][0]                  = std::sin(radians);
      }
    std::vector<
      GridTools::PeriodicFacePair<typename TriangulationType::cell_iterator>>
      pairs;
    GridTools::collect_periodic_faces(tria,
                                      angle == 0 ? 0 : 2,
                                      angle == 0 ? 1 : 3,
                                      angle == 0 ? 0 : 1,
                                      pairs,
                                      Tensor<1, dim>(),
                                      rotation);
    tria.add_periodicity(pairs);
    tria.refine_global(1);
    std::vector<Errors> results;
    for (unsigned int stage = 0; stage < 3; ++stage)
      {
        if (stage > 0)
          {
            for (const auto &cell : tria.active_cell_iterators())
              if (cell->is_locally_owned())
                {
                  const auto   center = cell->center();
                  const double radius = std::hypot(center[0], center[1]);
                  if (stage == 1 &&
                      (angle == 0 ? center[1] > 0.5 : radius > 1.5))
                    cell->set_refine_flag();
                  else if (stage == 2 && cell->level() > 1)
                    cell->set_coarsen_flag();
                }
            tria.execute_coarsening_and_refinement();
          }
        AssertThrow(
          tria.n_global_levels() == (stage == 1 ? 3 : 2),
          ExcMessage(
            "The refine/coarsen cycle did not change the mesh as intended."));
        const auto error =
          check_hierarchy(tria, communicator, angle, degree, components);
        AssertThrow(error.constraint < 1.e-12,
                    ExcMessage("Incorrect analytic periodic relation."));
        AssertThrow(error.trace < 1.e-11,
                    ExcMessage(
                      "The constrained traces do not satisfy periodicity."));
        AssertThrow(error.prolongation < 1.e-11,
                    ExcMessage("Linear field not preserved by transfer."));
        AssertThrow(error.adjoint < 1.e-11,
                    ExcMessage(
                      "Restriction is not the transpose of prolongation."));
        results.push_back(error);
      }
    return results;
  }
} // namespace RotatedPeriodicityTest
