// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/sparse_ilu.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/loop.h>
#include <deal.II/meshworker/output.h>

#include <deal.II/multigrid/mg_coarse.h>
#include <deal.II/multigrid/mg_matrix.h>
#include <deal.II/multigrid/mg_smoother.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_component.h>
#include <deal.II/multigrid/multigrid.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <functional>
#include <sstream>

#include "../tests.h"



template <int dim, typename number, int spacedim>
void
reinit_vector(const dealii::DoFHandler<dim, spacedim> &mg_dof,
              MGLevelObject<dealii::Vector<number>>   &v)
{
  for (unsigned int level = v.min_level(); level <= v.max_level(); ++level)
    {
      unsigned int n = mg_dof.n_dofs(level);
      v[level].reinit(n);
    }
}

template <int dim>
void
initialize(const DoFHandler<dim> &dof, Vector<double> &u)
{
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
       cell != dof.end();
       ++cell)
    {
      cell->get_dof_indices(dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        {
          const unsigned int comp =
            dof.get_fe().system_to_component_index(i).first;
          u(dof_indices[i]) = comp + 1;
        }
    }
}



template <int dim>
void
initialize(const DoFHandler<dim> &dof, MGLevelObject<Vector<double>> &u)
{
  unsigned int       counter       = 0;
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index>    dof_indices(dofs_per_cell);
  typename DoFHandler<dim>::cell_iterator cell = dof.begin(0);
  cell->get_mg_dof_indices(dof_indices);
  for (unsigned int i = 0; i < dofs_per_cell; ++i)
    u[0](dof_indices[i]) = ++counter;
}


template <int dim>
void
diff(Vector<double>        &diff,
     const DoFHandler<dim> &dof,
     const Vector<double>  &v,
     const unsigned int     level)
{
  diff.reinit(v);
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> mg_dof_indices(dofs_per_cell);
  const unsigned int                   n_comp = dof.get_fe().n_components();
  for (typename DoFHandler<dim>::cell_iterator cell = dof.begin(level);
       cell != dof.end(level);
       ++cell)
    {
      cell->get_mg_dof_indices(mg_dof_indices);
      for (unsigned int i = 0; i < dofs_per_cell / n_comp; ++i)
        {
          const unsigned int ni = n_comp * i;
          const unsigned int i0 = mg_dof_indices[ni];
          const unsigned int i1 = mg_dof_indices[ni + 1];
          const unsigned int i2 = mg_dof_indices[ni + 2];
          diff(i0)              = 2 * v(i0) - v(i1);
          diff(i1)              = 3 * v(i1) - 2 * v(i2);
          diff(i2)              = v(i2) - 3 * v(i0);
        }
    }
}

template <int dim>
void
diff(Vector<double>        &diff,
     const DoFHandler<dim> &dof_1,
     const DoFHandler<dim> &dof_2,
     const Vector<double>  &u,
     const Vector<double>  &v,
     const unsigned int     level)
{
  diff.reinit(u);
  const unsigned int        dofs_per_cell = dof_1.get_fe().dofs_per_cell;
  std::vector<unsigned int> dof_indices_1(dofs_per_cell);
  std::vector<unsigned int> dof_indices_2(dofs_per_cell);
  for (typename DoFHandler<dim>::cell_iterator cell_1 = dof_1.begin(level),
                                               cell_2 = dof_2.begin(level);
       cell_1 != dof_1.end(level);
       ++cell_1, ++cell_2)
    {
      cell_1->get_mg_dof_indices(dof_indices_1);
      cell_2->get_mg_dof_indices(dof_indices_2);
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        diff(dof_indices_1[i]) = u(dof_indices_1[i]) - v(dof_indices_2[i]);
    }
}

template <int dim>
void
print(const DoFHandler<dim>          &dof,
      std::vector<std::vector<bool>> &interface_dofs)
{
  const unsigned int dofs_per_cell = dof.get_fe().dofs_per_cell;
  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);
  for (unsigned int l = 0; l < dof.get_triangulation().n_levels(); ++l)
    {
      deallog << std::endl;
      deallog << "Level " << l << std::endl;
      for (typename DoFHandler<dim>::cell_iterator cell = dof.begin(l);
           cell != dof.end(l);
           ++cell)
        {
          cell->get_mg_dof_indices(dof_indices);
          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            deallog << ' ' << interface_dofs[l][dof_indices[i]];
        }
    }
}

template <int dim>
class OutputCreator : public EnableObserverPointer
{
public:
  void
  cell(MeshWorker::DoFInfo<dim> &dinfo, MeshWorker::IntegrationInfo<dim> &info);
};

template <int dim>
void
OutputCreator<dim>::cell(MeshWorker::DoFInfo<dim>         &dinfo,
                         MeshWorker::IntegrationInfo<dim> &info)
{
  const FEValuesBase<dim>                &fe = info.fe_values();
  const std::vector<std::vector<double>> &uh = info.values[0];

  const unsigned int square_root =
    static_cast<unsigned int>(std::pow(fe.n_quadrature_points, 1. / dim) + .5);
  for (unsigned int k1 = 0; k1 < square_root; ++k1)
    {
      for (unsigned int k2 = 0; k2 < square_root; ++k2)
        {
          for (unsigned int d = 0; d < dim; ++d)
            dinfo.quadrature_value(k1 * square_root + k2, d) =
              fe.quadrature_point(k1 * square_root + k2)[d];
          for (unsigned int i = 0; i < uh.size(); ++i)
            dinfo.quadrature_value(k1 * square_root + k2, i + dim) =
              uh[i][k1 * square_root + k2];
        }
    }
}


template <int dim>
class LaplaceProblem
{
public:
  using CellInfo = MeshWorker::IntegrationInfo<dim>;

  LaplaceProblem(const unsigned int deg);
  void
  run();

private:
  void
  setup_system();
  void
  test();
  void
  output_gpl(const DoFHandler<dim> &dof, MGLevelObject<Vector<double>> &v);
  void
  refine_local();

  Triangulation<dim>  triangulation;
  const MappingQ<dim> mapping;
  FESystem<dim>       fe;
  DoFHandler<dim>     mg_dof_handler;
  DoFHandler<dim>     mg_dof_handler_renumbered;

  const unsigned int                             degree;
  std::vector<std::set<types::global_dof_index>> boundary_indices,
    boundary_indices_renumbered;
};


template <int dim>
LaplaceProblem<dim>::LaplaceProblem(const unsigned int deg)
  : triangulation(Triangulation<dim>::limit_level_difference_at_vertices)
  , mapping(1)
  , fe(FE_Q<dim>(deg), 3)
  , mg_dof_handler(triangulation)
  , mg_dof_handler_renumbered(triangulation)
  , degree(deg)
{}


template <int dim>
void
LaplaceProblem<dim>::setup_system()
{
  mg_dof_handler.distribute_dofs(fe);
  mg_dof_handler.distribute_mg_dofs();
  mg_dof_handler_renumbered.distribute_dofs(fe);
  mg_dof_handler_renumbered.distribute_mg_dofs();

  const unsigned int nlevels = triangulation.n_levels();

  DoFHandler<dim> &dof = mg_dof_handler_renumbered;
  DoFRenumbering::component_wise(dof);
  for (unsigned int level = 0; level < nlevels; ++level)
    DoFRenumbering::component_wise(mg_dof_handler_renumbered, level);

  boundary_indices.resize(nlevels);
  boundary_indices_renumbered.resize(nlevels);

  for (unsigned int l = 0; l < nlevels; ++l)
    {
      boundary_indices_renumbered[l].clear();
      boundary_indices[l].clear();
    }

  deallog << "Number of degrees of freedom: " << mg_dof_handler.n_dofs();

  for (unsigned int l = 0; l < triangulation.n_levels(); ++l)
    deallog << "   " << 'L' << l << ": " << mg_dof_handler.n_dofs(l);
  deallog << std::endl;
}

template <int dim>
void
LaplaceProblem<dim>::output_gpl(const DoFHandler<dim>         &dof,
                                MGLevelObject<Vector<double>> &v)
{
  MeshWorker::IntegrationInfoBox<dim> info_box;
  const unsigned int n_gauss_points = dof.get_fe().tensor_degree();
  QTrapezoid<1>      trapez;
  QIterated<dim>     quadrature(trapez, n_gauss_points);
  info_box.cell_quadrature = quadrature;
  AnyData data;
  data.add<MGLevelObject<Vector<double>> *>(&v, "mg_vector");
  info_box.cell_selector.add("mg_vector");
  info_box.initialize_update_flags();
  UpdateFlags update_flags =
    update_quadrature_points | update_values | update_gradients;
  info_box.add_update_flags(update_flags, true, true, true, true);
  info_box.initialize(fe, mapping, data, v);

  MeshWorker::DoFInfo<dim> dof_info(dof);

  MeshWorker::Assembler::GnuplotPatch assembler;
  assembler.initialize(n_gauss_points + 1, dof.get_fe().n_components());

  for (unsigned int l = 0; l < triangulation.n_levels(); ++l)
    {
      OutputCreator<dim> matrix_integrator;
      MeshWorker::loop<dim,
                       dim,
                       MeshWorker::DoFInfo<dim>,
                       MeshWorker::IntegrationInfoBox<dim>>(
        dof.begin_mg(l),
        dof.end_mg(l),
        dof_info,
        info_box,
        std::bind(&OutputCreator<dim>::cell,
                  &matrix_integrator,
                  std::placeholders::_1,
                  std::placeholders::_2),
        nullptr,
        nullptr,
        assembler);
    }
}



template <int dim>
void
LaplaceProblem<dim>::test()
{
  MGConstrainedDoFs mg_constrained_dofs;
  mg_constrained_dofs.initialize(mg_dof_handler);
  mg_constrained_dofs.make_zero_boundary_constraints(mg_dof_handler, {0});
  MGConstrainedDoFs mg_constrained_dofs_renumbered;
  mg_constrained_dofs_renumbered.initialize(mg_dof_handler_renumbered);
  mg_constrained_dofs_renumbered.make_zero_boundary_constraints(
    mg_dof_handler_renumbered, {0});

  MGTransferPrebuilt<Vector<double>> mg_transfer(mg_constrained_dofs);
  mg_transfer.build(mg_dof_handler);
  MGTransferPrebuilt<Vector<double>> mg_transfer_renumbered(
    mg_constrained_dofs_renumbered);
  mg_transfer_renumbered.build(mg_dof_handler_renumbered);

  Vector<double> test;
  test.reinit(mg_dof_handler.n_dofs());

  MGLevelObject<Vector<double>> v(0, triangulation.n_levels() - 1);
  MGLevelObject<Vector<double>> u(0, triangulation.n_levels() - 1);
  MGLevelObject<Vector<double>> d(0, triangulation.n_levels() - 1);

  initialize(mg_dof_handler, test);
  mg_transfer.copy_to_mg(mg_dof_handler, v, test);

  initialize(mg_dof_handler_renumbered, test);
  mg_transfer_renumbered.copy_to_mg(mg_dof_handler_renumbered, u, test);
  for (unsigned int l = 0; l < triangulation.n_levels(); ++l)
    {
      diff(d[l], mg_dof_handler_renumbered, u[l], l);
      deallog << l << ' ' << u[l].l2_norm() << '\t' << v[l].l2_norm() << '\t'
              << d[l].l2_norm() << std::endl;
      for (unsigned int i = 0; i < d[l].size(); ++i)
        if (d[l](i) != 0)
          deallog << i << ' ' << d[l](i) << std::endl;
    }
  output_gpl(mg_dof_handler_renumbered, d);
}



template <int dim>
void
LaplaceProblem<dim>::refine_local()
{
  bool cell_refined = false;
  for (typename Triangulation<dim>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    {
      for (const unsigned int vertex : GeometryInfo<dim>::vertex_indices())
        {
          const Point<dim> p = cell->vertex(vertex);
          const Point<dim> origin =
            (dim == 2 ? Point<dim>(0, 0) : Point<dim>(0, 0, 0));
          const double dist = p.distance(origin);
          if (dist < 0.25 / numbers::PI)
            {
              cell->set_refine_flag();
              cell_refined = true;
              break;
            }
        }
    }
  // Wenn nichts verfeinert wurde bisher, global verfeinern!
  if (!cell_refined)
    for (typename Triangulation<dim>::active_cell_iterator cell =
           triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      cell->set_refine_flag();


  triangulation.execute_coarsening_and_refinement();
}

template <int dim>
void
LaplaceProblem<dim>::run()
{
  for (unsigned int cycle = 0; cycle < 6; ++cycle)
    {
      deallog << "Cycle " << cycle << std::endl;

      if (cycle == 0)
        {
          GridGenerator::hyper_cube(triangulation, -1, 1);
          triangulation.refine_global(1);
        }
      // triangulation.refine_global (1);
      refine_local();
      setup_system();
      test();
    };
}

int
main()
{
  initlog();

  LaplaceProblem<2> laplace_problem_2d(1);
  laplace_problem_2d.run();
}
