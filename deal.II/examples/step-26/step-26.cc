/* Author: Wolfgang Bangerth, Texas A&M University, 2008 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2013 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>

#include <fstream>
#include <iostream>


namespace Step26
{
  using namespace dealii;

  template<int dim>
  class HeatEquation
  {
  public:
    HeatEquation();
    void run();

  private:
    void setup_system();
    void solve_u();
    void output_results() const;

    Triangulation<dim> triangulation;
    FE_Q<dim> fe;
    DoFHandler<dim> dof_handler;

    ConstraintMatrix constraints;

    SparsityPattern sparsity_pattern;
    SparseMatrix<double> mass_matrix;
    SparseMatrix<double> laplace_matrix;
    SparseMatrix<double> matrix_u;

    Vector<double> solution_u;
    Vector<double> old_solution_u;
    Vector<double> system_rhs;

    double time, time_step;
    unsigned int timestep_number;
    const double theta;
  };

  //-------------------------------------

  template<int dim>
  class RightHandSide: public Function<dim>
  {
  public:
    RightHandSide() :
    Function<dim>(),
    period (0.2)
    {}

    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;

  private:
    const double period;
  };

  template<int dim>
  double RightHandSide<dim>::value(const Point<dim> &p,
                                   const unsigned int component) const
  {
    Assert (component == 0, ExcInternalError());
    Assert (dim == 2, ExcNotImplemented());

    return std::sin(p[0]*numbers::PI) * std::sin(p[1]*numbers::PI);


    const double time = this->get_time();
    const double point_within_period = (time/period - std::floor(time/period));

    if ((point_within_period >= 0.0) && (point_within_period <= 0.2))
      {
	if ((p[0] > 0.5) && (p[1] > -0.5))
	  return 1;
	else
	  return 0;
      }
    else if ((point_within_period >= 0.5) && (point_within_period <= 0.7))
      {
	if ((p[0] > -0.5) && (p[1] > 0.5))
	  return 1;
	else
	  return 0;
      }
    else
      return 0;
  }

  template<int dim>
  class BoundaryValuesU: public Function<dim>
  {
  public:
    BoundaryValuesU() :
      Function<dim>()
    {
    }
    virtual ~BoundaryValuesU()
    {
    }
    virtual double value(const Point<dim> &p,
                         const unsigned int component = 0) const;
  };

  template<int dim>
  double BoundaryValuesU<dim>::value(const Point<dim> &/*p*/,
                                     const unsigned int component) const
  {
    Assert(component == 0, ExcInternalError());
    return 0; // Zero-Dirichlet Boundary
  }

  template<int dim>
  HeatEquation<dim>::HeatEquation() :
    fe(1), dof_handler(triangulation), time_step(1. / 500), theta(0.5)
  {
  }

  template<int dim>
  void HeatEquation<dim>::setup_system()
  {
    //    GridGenerator::hyper_L (triangulation);
    GridGenerator::hyper_cube (triangulation, 0, 1);
    triangulation.refine_global (5);

    std::cout << "Number of active cells: " << triangulation.n_active_cells()
              << std::endl;

    dof_handler.distribute_dofs(fe);

    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl << std::endl;

    sparsity_pattern.reinit(dof_handler.n_dofs(), dof_handler.n_dofs(),
                            dof_handler.max_couplings_between_dofs());
    DoFTools::make_sparsity_pattern(dof_handler, sparsity_pattern);
    sparsity_pattern.compress();

    mass_matrix.reinit(sparsity_pattern);
    laplace_matrix.reinit(sparsity_pattern);
    matrix_u.reinit(sparsity_pattern);

    MatrixCreator::create_mass_matrix(dof_handler, QGauss<dim>(3), mass_matrix);
    MatrixCreator::create_laplace_matrix(dof_handler, QGauss<dim>(3),
                                         laplace_matrix);

    solution_u.reinit(dof_handler.n_dofs());
    old_solution_u.reinit(dof_handler.n_dofs());
    system_rhs.reinit(dof_handler.n_dofs());

    constraints.close();
  }

  template<int dim>
  void HeatEquation<dim>::solve_u()
  {
    SolverControl solver_control(1000, 1e-8 * system_rhs.l2_norm());
    SolverCG<> cg(solver_control);

    cg.solve(matrix_u, solution_u, system_rhs, PreconditionIdentity());

    std::cout << "   u-equation: " << solver_control.last_step()
              << " CG iterations." << std::endl;
  }

  template<int dim>
  void HeatEquation<dim>::output_results() const
  {
    DataOut<dim> data_out;

    data_out.attach_dof_handler(dof_handler);
    data_out.add_data_vector(solution_u, "U");

    data_out.build_patches();

    const std::string filename = "solution-"
                                 + Utilities::int_to_string(timestep_number, 3) + ".vtk";
    std::ofstream output(filename.c_str());
    data_out.write_vtk(output);

    std::cout << "    max= " << time << ' ' << solution_u.linfty_norm() << std::endl;
  }

  template<int dim>
  void HeatEquation<dim>::run()
  {
    setup_system();

    VectorTools::interpolate(dof_handler, ZeroFunction<dim>(), solution_u);

    timestep_number = 0;
    output_results();

    VectorTools::interpolate(dof_handler, ZeroFunction<dim>(),
                             old_solution_u);

    Vector<double> tmp(solution_u.size());
    Vector<double> forcing_terms(solution_u.size());

    for (timestep_number = 1, time = time_step; time <= 0.5;
         time += time_step, ++timestep_number)
      {
        std::cout << "Time step " << timestep_number << " at t=" << time
                  << std::endl;

        mass_matrix.vmult(system_rhs, old_solution_u);

        laplace_matrix.vmult(tmp, old_solution_u);
        system_rhs.add(-(1 - theta) * time_step, tmp);

        RightHandSide<dim> rhs_function;
        rhs_function.set_time(time);
        VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(2),
                                            rhs_function, tmp);
        forcing_terms = tmp;
        forcing_terms *= theta;

        rhs_function.set_time(time - time_step);
        VectorTools::create_right_hand_side(dof_handler, QGauss<dim>(2),
                                            rhs_function, tmp);

        forcing_terms.add(1 - theta, tmp);

        system_rhs.add(time_step, forcing_terms);

        {
          BoundaryValuesU<dim> boundary_values_u_function;
          boundary_values_u_function.set_time(time);

          std::map<unsigned int, double> boundary_values;
          VectorTools::interpolate_boundary_values(dof_handler, 0,
                                                   boundary_values_u_function, boundary_values);

          matrix_u.copy_from(mass_matrix);
          matrix_u.add(theta * time_step, laplace_matrix); // I omit here a time_step*theta
          MatrixTools::apply_boundary_values(boundary_values, matrix_u,
                                             solution_u, system_rhs);
        }
        solve_u();

        output_results();

        old_solution_u = solution_u;
      }
  }
}

int main()
{
  try
    {
      using namespace dealii;
      using namespace Step26;

      deallog.depth_console(0);

      HeatEquation<2> heat_equation_solver;
      heat_equation_solver.run();

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl << exc.what()
                << std::endl << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl << "Aborting!"
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
