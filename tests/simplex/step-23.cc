// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Step-23 with tetrahedron mesh. Following incompatible modifications had to be
// made:
// - Change FE_Q and QGauss to FE_SimplexP and QGaussSimplex.
// - Explicit use of MappingFE instead of the default mapping.
// - Grid generation by subdivided_hyper_cube_with_simplices instead of
// hyper_cube, because global refinement is not allowed.

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"

// simplex
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_simplex_p_bubbles.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_generator.h>

// #define HEX

const unsigned int degree = 1;

namespace Step23
{

  template <int dim>
  class WaveEquation
  {
  public:
    WaveEquation();
    void
    run();

  private:
    void
    setup_system();
    void
    solve_u();
    void
    solve_v();
    void
    output_results() const;

    Triangulation<dim> triangulation;
#ifdef HEX
    MappingQ<dim, dim> mapping;
    FE_Q<dim>          fe;
    QGauss<dim>        quadrature;
#else
    MappingFE<dim, dim> mapping;
    FE_SimplexP<dim>    fe;
    QGaussSimplex<dim>  quadrature;
#endif
    DoFHandler<dim> dof_handler;

    AffineConstraints<double> constraints;

    SparsityPattern      sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> matrix_u;
    SparseMatrix<double> matrix_v;

    Vector<double> solution_u, solution_v;
    Vector<double> old_solution_u, old_solution_v;
    Vector<double> system_rhs;

    double       time_step;
    double       time;
    unsigned int timestep_number;
    const double theta;
  };

  template <int dim>
  class InitialValuesU : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int component = 0) const override
    {
      (void)component;
      Assert(component == 0, ExcIndexRange(component, 0, 1));
      return 0;
    }
  };



  template <int dim>
  class InitialValuesV : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int component = 0) const override
    {
      (void)component;
      Assert(component == 0, ExcIndexRange(component, 0, 1));
      return 0;
    }
  };



  template <int dim>
  class RightHandSide : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> & /*p*/,
          const unsigned int component = 0) const override
    {
      (void)component;
      Assert(component == 0, ExcIndexRange(component, 0, 1));
      return 0;
    }
  };



  template <int dim>
  class BoundaryValuesU : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      (void)component;
      Assert(component == 0, ExcIndexRange(component, 0, 1));

      if ((this->get_time() <= 0.5) && (p[0] < 0) && (p[1] < 1. / 3) &&
          (p[1] > -1. / 3))
        return std::sin(this->get_time() * 4 * numbers::PI);
      else
        return 0;
    }
  };



  template <int dim>
  class BoundaryValuesV : public Function<dim>
  {
  public:
    virtual double
    value(const Point<dim> &p, const unsigned int component = 0) const override
    {
      (void)component;
      Assert(component == 0, ExcIndexRange(component, 0, 1));

      if ((this->get_time() <= 0.5) && (p[0] < 0) && (p[1] < 1. / 3) &&
          (p[1] > -1. / 3))
        return (std::cos(this->get_time() * 4 * numbers::PI) * 4 * numbers::PI);
      else
        return 0;
    }
  };

  template <int dim>
  WaveEquation<dim>::WaveEquation()
#ifdef HEX
    : mapping(1)
#else
    : mapping(FE_SimplexP<dim>(1))
#endif
    , fe(1)
    , quadrature(fe.degree + 1)
    , dof_handler(triangulation)
    , time_step(1. / 64)
    , time(time_step)
    , timestep_number(1)
    , theta(0.5)
  {}

  template <int dim>
  void
  WaveEquation<dim>::setup_system()
  {
    const unsigned int n_subdivisions = (1 << 3);
#ifdef HEX
    GridGenerator::subdivided_hyper_cube(triangulation, n_subdivisions, -1, 1);
#else
    GridGenerator::subdivided_hyper_cube_with_simplices(triangulation,
                                                        n_subdivisions,
                                                        -1,
                                                        1);
#endif

    deallog << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;

    dof_handler.distribute_dofs(fe);

    deallog << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl
            << std::endl;

    DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, dsp);
    sparsity_pattern.copy_from(dsp);

    mass_matrix.reinit(sparsity_pattern);
    laplace_matrix.reinit(sparsity_pattern);
    matrix_u.reinit(sparsity_pattern);
    matrix_v.reinit(sparsity_pattern);

    MatrixCreator::create_mass_matrix(mapping,
                                      dof_handler,
                                      quadrature,
                                      mass_matrix);
    MatrixCreator::create_laplace_matrix(mapping,
                                         dof_handler,
                                         quadrature,
                                         laplace_matrix);
    solution_u.reinit(dof_handler.n_dofs());
    solution_v.reinit(dof_handler.n_dofs());
    old_solution_u.reinit(dof_handler.n_dofs());
    old_solution_v.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.close();
  }


  template <int dim>
  void
  WaveEquation<dim>::solve_u()
  {
    SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm());
    SolverCG<Vector<double>> cg(solver_control);

    cg.solve(matrix_u, solution_u, system_rhs, PreconditionIdentity());

    deallog << "   u-equation: " << solver_control.last_step()
            << " CG iterations." << std::endl;
  }



  template <int dim>
  void
  WaveEquation<dim>::solve_v()
  {
    SolverControl            solver_control(1000, 1e-8 * system_rhs.l2_norm());
    SolverCG<Vector<double>> cg(solver_control);

    cg.solve(matrix_v, solution_v, system_rhs, PreconditionIdentity());

    deallog << "   v-equation: " << solver_control.last_step()
            << " CG iterations." << std::endl;
  }


  template <int dim>
  void
  WaveEquation<dim>::output_results() const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution_u, "U");
    data_out.add_data_vector(solution_v, "V");

    data_out.build_patches(mapping);

    const std::string filename =
      "solution-" + Utilities::int_to_string(timestep_number, 3) + ".vtu";
    DataOutBase::VtkFlags vtk_flags;
    vtk_flags.compression_level = DataOutBase::CompressionLevel::best_speed;
    data_out.set_flags(vtk_flags);
    std::ofstream output(filename);
    data_out.write_vtu(output);
  }


  template <int dim>
  void
  WaveEquation<dim>::run()
  {
    setup_system();

    VectorTools::project(mapping,
                         dof_handler,
                         constraints,
                         quadrature,
                         InitialValuesU<dim>(),
                         old_solution_u);
    VectorTools::project(mapping,
                         dof_handler,
                         constraints,
                         quadrature,
                         InitialValuesV<dim>(),
                         old_solution_v);

    Vector<double> tmp(solution_u.size());
    Vector<double> forcing_terms(solution_u.size());

    for (; time <= 5; time += time_step, ++timestep_number)
      {
        deallog << "Time step " << timestep_number << " at t=" << time
                << std::endl;

        mass_matrix.vmult(system_rhs, old_solution_u);

        mass_matrix.vmult(tmp, old_solution_v);
        system_rhs.add(time_step, tmp);

        laplace_matrix.vmult(tmp, old_solution_u);
        system_rhs.add(-theta * (1 - theta) * time_step * time_step, tmp);

        RightHandSide<dim> rhs_function;
        rhs_function.set_time(time);
        VectorTools::create_right_hand_side(
          mapping, dof_handler, quadrature, rhs_function, tmp);
        forcing_terms = tmp;
        forcing_terms *= theta * time_step;

        rhs_function.set_time(time - time_step);
        VectorTools::create_right_hand_side(
          mapping, dof_handler, quadrature, rhs_function, tmp);

        forcing_terms.add((1 - theta) * time_step, tmp);

        system_rhs.add(theta * time_step, forcing_terms);

        {
          BoundaryValuesU<dim> boundary_values_u_function;
          boundary_values_u_function.set_time(time);

          std::map<types::global_dof_index, double> boundary_values;
          VectorTools::interpolate_boundary_values(mapping,
                                                   dof_handler,
                                                   0,
                                                   boundary_values_u_function,
                                                   boundary_values);

          matrix_u.copy_from(mass_matrix);
          matrix_u.add(theta * theta * time_step * time_step, laplace_matrix);
          MatrixTools::apply_boundary_values(boundary_values,
                                             matrix_u,
                                             solution_u,
                                             system_rhs);
        }
        solve_u();

        laplace_matrix.vmult(system_rhs, solution_u);
        system_rhs *= -theta * time_step;

        mass_matrix.vmult(tmp, old_solution_v);
        system_rhs += tmp;

        laplace_matrix.vmult(tmp, old_solution_u);
        system_rhs.add(-time_step * (1 - theta), tmp);

        system_rhs += forcing_terms;

        {
          BoundaryValuesV<dim> boundary_values_v_function;
          boundary_values_v_function.set_time(time);

          std::map<types::global_dof_index, double> boundary_values;
          VectorTools::interpolate_boundary_values(mapping,
                                                   dof_handler,
                                                   0,
                                                   boundary_values_v_function,
                                                   boundary_values);
          matrix_v.copy_from(mass_matrix);
          MatrixTools::apply_boundary_values(boundary_values,
                                             matrix_v,
                                             solution_v,
                                             system_rhs);
        }
        solve_v();

        output_results();

        deallog << "   Total energy: "
                << (mass_matrix.matrix_norm_square(solution_v) +
                    laplace_matrix.matrix_norm_square(solution_u)) /
                     2
                << std::endl;

        old_solution_u = solution_u;
        old_solution_v = solution_v;
      }
  }
} // namespace Step23


int
main(int argc, char **argv)
{
  initlog();

  deallog.depth_file(1);

  try
    {
      using namespace Step23;

      WaveEquation<2> wave_equation_solver;
      wave_equation_solver.run();
    }
  catch (const std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
