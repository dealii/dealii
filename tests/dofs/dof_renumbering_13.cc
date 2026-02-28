// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values_extractors.h>
#include <deal.II/fe/fe_values_views.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim = dim>
void
print_local_order(const FiniteElement<dim, spacedim> &fe)
{
  const unsigned int ndofs = fe.n_dofs_per_cell();
  const bool         nodal = fe.has_support_points();
  const auto        &unit_pts =
    nodal ? fe.get_unit_support_points() : std::vector<Point<dim>>{};
  deallog << "FE: " << fe.get_name() << ", n_dofs_per_cell=" << ndofs
          << std::endl;

  for (unsigned int i = 0; i < ndofs; ++i)
    {
      const auto comp_info =
        fe.system_to_component_index(i); // {component, base_index}
      const unsigned int comp = comp_info.first;
      const unsigned int base = comp_info.second;
      deallog << "FESystem local dof index i=" << i
              << " | system component=" << comp;
      //        << " | component shape function index=" << base;
      // if (nodal && i < unit_pts.size())
      //  deallog << " | unit_point=" << unit_pts[i];
      deallog << std::endl;
    }
}



template <int dim>
void
print_dofs(const DoFHandler<dim> &dof)
{
  std::vector<types::global_dof_index> v(dof.get_fe().dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
       cell != dof.end();
       ++cell)
    {
      deallog << "Cell " << cell << " -- ";
      cell->get_dof_indices(v);

      for (unsigned int i = 0; i < v.size(); ++i)
        {
          deallog << v[i] << ' ';
        }
      deallog << std::endl;
    }
}



template <int dim>
void
print_dofs(const DoFHandler<dim> &dof, unsigned int level)
{
  std::vector<types::global_dof_index> v(dof.get_fe().dofs_per_cell);
  for (typename DoFHandler<dim>::cell_iterator cell = dof.begin(level);
       cell != dof.end(level);
       ++cell)
    {
      deallog << "Cell " << cell << " -- ";
      cell->get_mg_dof_indices(v);
      for (unsigned int i = 0; i < v.size(); ++i)
        {
          deallog << v[i] << ' ';
        }
      deallog << std::endl;
    }
}

template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
check_dofs_bounds(const DoFHandler<dim> &dof_handler, unsigned int level)
{
  // Primitive implementation to keep testing straightforward. Maybe implement
  // with ComponentMask.
  std::vector<std::vector<types::global_dof_index>> component_dofs(
    dof_handler.get_fe().n_components());

  std::vector<types::global_dof_index> v(dof_handler.get_fe().dofs_per_cell);
  for (typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin(level);
       cell != dof_handler.end(level);
       ++cell)
    {
      cell->get_mg_dof_indices(v);
      for (unsigned int i = 0; i < v.size(); ++i)
        {
          const auto comp_info =
            dof_handler.get_fe().system_to_component_index(i);
          const unsigned int comp = comp_info.first;
          component_dofs[comp].push_back(v[i]);
        }
    }
  std::vector<std::pair<unsigned int, unsigned int>> component_minmax_dofs(
    dof_handler.get_fe().n_components());

  deallog.push("DoF bounds");
  for (unsigned int i = 0; i < dof_handler.get_fe().n_components(); ++i)
    {
      auto asd =
        *std::min_element(component_dofs[i].begin(), component_dofs[i].end());
      component_minmax_dofs[i] = std::make_pair(
        *std::min_element(component_dofs[i].begin(), component_dofs[i].end()),
        *std::max_element(component_dofs[i].begin(), component_dofs[i].end()));
      deallog << "Component " << i << ": [" << component_minmax_dofs[i].first
              << " " << component_minmax_dofs[i].second << "]" << std::endl;
    }

  bool overlap = false;
  for (unsigned int i = 1; i < dof_handler.get_fe().n_components(); ++i)
    {
      const auto &a = component_minmax_dofs[i - 1];
      const auto &b = component_minmax_dofs[i];

      // Two ranges [a.first, a.second] and [b.first, b.second] overlap if:
      // a.first <= b.second && b.first <= a.second
      if (a.first <= b.second && b.first <= a.second)
        overlap = true;
    }

  deallog << "Non-overlapping: " << std::boolalpha << !overlap << std::endl;
  deallog.pop();

  return component_minmax_dofs;
}


template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
check_dofs_bounds(const DoFHandler<dim> &dof_handler)
{
  // Primitive implementation to keep testing straightforward. Maybe implement
  // with ComponentMask.
  std::vector<std::vector<types::global_dof_index>> component_dofs(
    dof_handler.get_fe().n_components());

  std::vector<types::global_dof_index> v(dof_handler.get_fe().dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell)
    {
      cell->get_dof_indices(v);
      for (unsigned int i = 0; i < v.size(); ++i)
        {
          const auto comp_info =
            dof_handler.get_fe().system_to_component_index(i);
          const unsigned int comp = comp_info.first;
          component_dofs[comp].push_back(v[i]);
        }
    }

  std::vector<std::pair<unsigned int, unsigned int>> component_minmax_dofs(
    dof_handler.get_fe().n_components());

  deallog.push("DoF bounds");
  for (unsigned int i = 0; i < dof_handler.get_fe().n_components(); ++i)
    {
      auto asd =
        *std::min_element(component_dofs[i].begin(), component_dofs[i].end());
      component_minmax_dofs[i] = std::make_pair(
        *std::min_element(component_dofs[i].begin(), component_dofs[i].end()),
        *std::max_element(component_dofs[i].begin(), component_dofs[i].end()));
      deallog << "Component " << i << ": [" << component_minmax_dofs[i].first
              << " " << component_minmax_dofs[i].second << "]" << std::endl;
    }

  bool overlap = false;
  for (unsigned int i = 1; i < dof_handler.get_fe().n_components(); ++i)
    {
      const auto &a = component_minmax_dofs[i - 1];
      const auto &b = component_minmax_dofs[i];

      // Two ranges [a.first, a.second] and [b.first, b.second] overlap if:
      // a.first <= b.second && b.first <= a.second
      if (a.first <= b.second && b.first <= a.second)
        overlap = true;
    }

  deallog << "Overlap: " << std::boolalpha << overlap << std::endl;
  deallog.pop();

  return component_minmax_dofs;
}

// Helper function that returns the components `order` given dof_ranges.
std::vector<unsigned int>
dofs_component_order(
  std::vector<std::pair<unsigned int, unsigned int>> component_dof_ranges)
{
  // Build order (indices of components sorted by min, then max).
  std::vector<unsigned int> order(component_dof_ranges.size());

  // 0, 1, 2, ...
  std::iota(order.begin(), order.end(), 0);

  // Sort by comparing range's min with other range's min.
  std::sort(order.begin(), order.end(), [&](unsigned int a, unsigned int b) {
    const auto &ra = component_dof_ranges[a];
    const auto &rb = component_dof_ranges[b];
    return ra.first < rb.first;
  });

  return order;
}


template <int dim>
void
check_renumbering(
  DoFHandler<dim>                                     &dof_handler,
  const std::vector<FEValuesExtractors::AnyExtractor> &extractors,
  const std::vector<unsigned int>                     &requested_order)
{
  // Check global ordering
  deallog.push("Initial  ");
  print_dofs(dof_handler);
  auto dof_bounds = check_dofs_bounds(dof_handler);
  deallog.pop();

  DoFRenumbering::component_wise(dof_handler, extractors);
  deallog.push("Reordered");
  print_dofs(dof_handler);
  auto dof_bounds_ordered = check_dofs_bounds(dof_handler);
  auto resultant_order    = dofs_component_order(dof_bounds_ordered);
  deallog.push("Component DoF Order");
  deallog << "Expected: " << requested_order << std::endl;
  deallog << "Obtained: " << resultant_order << std::endl;
  deallog << "Correct: " << std::boolalpha << " "
          << bool(resultant_order == requested_order) << std::endl;
  deallog.pop();
  deallog.pop();

  deallog.push("Multigrid");
  // Check level ordering
  for (unsigned int level = 0;
       level < dof_handler.get_triangulation().n_levels();
       ++level)
    {
      deallog.push("level " + std::to_string(level));
      deallog.push("Initial  ");
      print_dofs(dof_handler, level);
      deallog.pop();

      DoFRenumbering::component_wise(dof_handler, level, extractors);

      deallog.push("Reordered");
      print_dofs(dof_handler, level);
      auto dof_bounds_ordered = check_dofs_bounds(dof_handler, level);
      auto resultant_order    = dofs_component_order(dof_bounds_ordered);
      deallog.push("Component DoF Order");
      deallog << "Expected: " << requested_order << std::endl;
      deallog << "Obtained: " << resultant_order << std::endl;
      deallog << "Correct: " << std::boolalpha << " "
              << bool(resultant_order == requested_order) << std::endl;
      deallog.pop();
      deallog.pop();
      deallog.pop();
    }
  deallog.pop();
}


template <int dim>
void
check()
{
  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  if (dim == 2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1, 1);

  tr.reset_all_manifolds();
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  if (dim == 1)
    tr.refine_global(2);

  // Setup extractors
  auto requested_order = std::vector<unsigned int>{1, 0, 2};
  std::vector<FEValuesExtractors::AnyExtractor> extractors = {
    FEValuesExtractors::Scalar(
      1), // system component 1 global dofs will have the lowest indices
    FEValuesExtractors::Scalar(
      0), // system component 0 global dofs will have second lowest indices
    FEValuesExtractors::Scalar(
      2) // system component 2 global dofs will have highest indices
  };

  DoFHandler<dim> dof_handler(tr);

  FESystem<dim> e1(FE_Q<dim>(2), 1, FE_DGQ<dim>(1), 1, FE_Q<dim>(2), 1);
  print_local_order(e1);
  dof_handler.distribute_dofs(e1);
  dof_handler.distribute_mg_dofs();
  check_renumbering(dof_handler, extractors, requested_order);
  dof_handler.clear();

  FESystem<dim> e2(FE_DGP<dim>(2), 2, FE_DGQ<dim>(1), 1);
  print_local_order(e2);
  dof_handler.distribute_dofs(e2);
  dof_handler.distribute_mg_dofs();
  check_renumbering(dof_handler, extractors, requested_order);
  dof_handler.clear();
}



template <int dim, int spacedim = dim>
void
check_all_supported()
{
  const auto scalar_components = 1;

  const auto vector_components =
    FEValuesViews::Vector<dim, spacedim>::value_type::n_independent_components;

  const auto tensor_components = FEValuesViews::Tensor<2, dim, spacedim>::
    value_type::n_independent_components;

  const auto symtensor_components = FEValuesViews::
    SymmetricTensor<2, dim, spacedim>::value_type::n_independent_components;


  int total_number_of_components = scalar_components + vector_components +
                                   tensor_components + symtensor_components;



  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  if (dim == 2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1, 1);

  tr.reset_all_manifolds();
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  if (dim == 1)
    tr.refine_global(2);

  DoFHandler<dim> dof_handler(tr);
  FESystem<dim>   e_all(FE_Q<dim>(1), total_number_of_components);

  // Setup extractors
  //! 0
  auto pressure = FEValuesExtractors::Scalar(0);

  // 1, [2]
  auto velocity = FEValuesExtractors::Vector(0 + scalar_components);

  // 3, [4, 5, 6]
  auto stress = FEValuesExtractors::Tensor<2>(0 + 1 + vector_components);

  // 7, [8, 9]
  auto stress_sym = FEValuesExtractors::SymmetricTensor<2>(
    0 + 1 + vector_components + tensor_components); //

  std::vector<FEValuesExtractors::AnyExtractor> extractors = {stress,
                                                              pressure,
                                                              stress_sym,
                                                              velocity};

  std::vector<unsigned int> requested_order;
  if (dim == 1)
    {
      requested_order = std::vector<unsigned int>{
        2, // stress
        0, // pressure
        3, // stress_sym
        1  // velocity
      };
    }
  if (dim == 2)
    {
      requested_order = std::vector<unsigned int>{
        // stress
        3,
        4,
        5,
        6,

        // pressure
        0,

        // stress_sym
        7,
        8,
        9,

        // velocity
        1,
        2
        // dofs of components 3, 4, 5, 6 will have lowest indices, staggered.
        // dofs of components 1, 2 will have highest indices, staggered.
      };
    }
  else if (dim == 3)
    {
      requested_order = std::vector<unsigned int>{
        // stress
        4,
        5,
        6,
        7,
        8,
        9,
        10,
        11,
        12,

        // pressure
        0,

        // stress_sym
        13,
        14,
        15,
        16,
        17,
        18,

        // velocity
        1,
        2,
        3

        // dofs of component 4, 5, 6, 7, 8, 9, 10, 11, 12 will have lowest
        // indices, staggered dofs of component 1, 2, 3 will have highest
        // indices
      };
    }

  print_local_order(e_all);
  dof_handler.distribute_dofs(e_all);
  dof_handler.distribute_mg_dofs();
  check_renumbering(dof_handler, extractors, requested_order);
  dof_handler.clear();
}

template <int dim, int spacedim = dim>
void
check_invalid_input()
{
  const auto scalar_components = 1;

  const auto vector_components =
    FEValuesViews::Vector<dim, spacedim>::value_type::n_independent_components;

  const auto tensor_components = FEValuesViews::Tensor<2, dim, spacedim>::
    value_type::n_independent_components;

  const auto symtensor_components = FEValuesViews::
    SymmetricTensor<2, dim, spacedim>::value_type::n_independent_components;


  int total_number_of_components = scalar_components + vector_components +
                                   tensor_components + symtensor_components;

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  if (dim == 2)
    GridGenerator::hyper_ball(tr, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(tr, -1, 1);

  tr.reset_all_manifolds();
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();
  if (dim == 1)
    tr.refine_global(2);

  DoFHandler<dim> dof_handler(tr);
  FESystem<dim>   e_all(FE_Q<dim>(1), total_number_of_components);
  dof_handler.distribute_dofs(e_all);

  // Setup extractors
  //! 0
  auto pressure = FEValuesExtractors::Scalar(0);

  // 1, [2]
  auto velocity = FEValuesExtractors::Vector(0 + scalar_components);

  // 3, [4, 5, 6,]
  auto stress = FEValuesExtractors::Tensor<2>(0 + 1 + vector_components);

  // 7, [8, 9]
  auto stress_sym = FEValuesExtractors::SymmetricTensor<2>(
    0 + 1 + vector_components + tensor_components);


  // Invalid input case 1: Missing an extractor
  try
    {
      DoFRenumbering::component_wise(dof_handler,
                                     {
                                       stress, pressure, stress_sym,
                                       // velocity
                                     });
    }
  catch (const ExcMessage &e)
    {
      deallog << "Missing extractors. " << std::endl;
      e.print_info(deallog.get_file_stream());
    }
  catch (...)
    {
      Assert(false, ExcInternalError());
    }

  // Invalid input case 2: duplicate extractor
  try
    {
      DoFRenumbering::component_wise(
        dof_handler, {stress, pressure, stress_sym, stress_sym});
    }
  catch (const ExcMessage &e)
    {
      deallog << "Duplicate extractors. " << std::endl;
      e.print_info(deallog.get_file_stream());
    }
  catch (...)
    {
      Assert(false, ExcInternalError());
    }

  // Invalid input case 3: Overlapping extractor
  try
    {
      DoFRenumbering::component_wise(dof_handler,
                                     {FEValuesExtractors::Tensor<2>(
                                        0 + 1 + vector_components - 1),
                                      pressure,
                                      stress_sym,
                                      velocity});
    }
  catch (const ExcMessage &e)
    {
      deallog << "Overlapping extractors. " << std::endl;
      e.print_info(deallog.get_file_stream());
      // deallog << e.get_exc_name() << std::endl;
    }
  catch (...)
    {
      Assert(false, ExcInternalError());
    }

  // Invalid input case 4: Out of bounds component.
  try
    {
      DoFRenumbering::component_wise(dof_handler,
                                     {stress,
                                      pressure,
                                      stress_sym,
                                      velocity,
                                      FEValuesExtractors::Scalar(1234)});
    }
  catch (const ExcIndexRange &e)
    {
      deallog << "Out of bounds component. " << std::endl;
      e.print_info(deallog.get_file_stream());
    }
  catch (...)
    {
      Assert(false, ExcInternalError());
    }
}

int
main()
{
  initlog();
  deallog << std::setprecision(2) << std::fixed;

  // Test case 1 - Using only `Scalar` and `Vector` extractors
  deallog.push("Basic");
  deallog.push("1d");
  check<1>();
  deallog.pop();
  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();

  // Test case 2 - Using all supported AnyExtractors
  deallog.push("All");
  deallog.push("1d");
  check_all_supported<1>();
  deallog.pop();
  deallog.push("2d");
  check_all_supported<2>();
  deallog.pop();
  deallog.push("3d");
  check_all_supported<3>();
  deallog.pop();
  deallog.pop();

  // Test case 3 - Invalid inputs
  deallog.push("Invalid input");
  deallog.push("2d");
  check_invalid_input<2>();
  deallog.pop();
  deallog.pop();
}
