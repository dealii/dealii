/* ---------------------------------------------------------------------
 * $Id$
 *
 * Copyright (C) 2012 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Authors: Joerg Frohne, Texas A&M University and
 *                        University of Siegen, 2012, 2013
 *          Wolfgang Bangerth, Texas A&M University, 2012, 2013
 *          Timo Heister, Texas A&M University, 2013
 */

// @sect3{Include files}
// We are using the the same
// include files as in step-41:
#include <deal.II/grid/tria.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_tools.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/fe_field_function.h>
#include <deal.II/distributed/solution_transfer.h>
#include <deal.II/base/timer.h>
#include <fstream>
#include <iostream>
#include <list>
#include <time.h>

#include <sys/stat.h>

#include <deal.II/base/logstream.h>

namespace Step42
{
  using namespace dealii;

// @sect3{The <code>Input</code> class template}

// This class has the the only purpose
// to read in data from a picture file
// that has to be stored in pbm ascii
// format. This data will be bilinear
// interpolated and provides in this way
// a function which describes an obstacle.
//
// The data which we read in by the
// function read_obstacle () from the file
// "obstacle_file.pbm" will be stored
// in a double std::vector named
// obstacle_data.
// This vector composes the base
// to calculate a piecewise bilinear
// function as a polynomial interpolation.
// This will be done by obstacle_function ().
//
// In the function run () of the class
// PlasticityContactProblem we create
// an object of the class Input which will
// be used in the class Obstacle to
// supply the obstacle function in
// update_solution_and_constraints () of
// the class PlasticityContactProblem.

  template <int dim>
    class Input
    {
      public:
        Input (const std::string &name)
            :
                mpi_communicator(MPI_COMM_WORLD),
                pcout(std::cout,
                    (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
                obstacle_data(0),
                hx(0),
                hy(0),
                nx(0),
                ny(0)
        {
          read_obstacle(name);
        }

        double
        hv (
            int i, int j);

        double
        obstacle_function (const double x,
        		            const double y);

        void
        read_obstacle (const std::string name);

      private:
        MPI_Comm mpi_communicator;
        ConditionalOStream pcout;
        std::vector<double> obstacle_data;
        double hx, hy;
        int nx, ny;
    };

// This function is used in obstacle_function ()
// to provide the proper value of the obstacle.
  template <int dim>
    double
    Input<dim>::hv (const int i,
    		        const int j)
    {
      assert(i >= 0 && i < nx);
      assert(j >= 0 && j < ny);
      return obstacle_data[nx * (ny - 1 - j) + i]; // i indiziert x-werte, j indiziert y-werte
    }

// obstacle_function () calculates the bilinear interpolated
// value in the point (x,y).
  template <int dim>
    double
    Input<dim>::obstacle_function (const double x,
    		                       const double y)
    {
      int ix = (int) (x / hx);
      int iy = (int) (y / hy);

      if (ix < 0)
        ix = 0;

      if (iy < 0)
        iy = 0;

      if (ix >= nx - 1)
        ix = nx - 2;

      if (iy >= ny - 1)
        iy = ny - 2;

      double val = 0.0;
        {
          FullMatrix<double> H(4, 4);
          Vector<double> X(4);
          Vector<double> b(4);

          double xx = 0.0;
          double yy = 0.0;

          xx = ix * hx;
          yy = iy * hy;
          H(0, 0) = xx;
          H(0, 1) = yy;
          H(0, 2) = xx * yy;
          H(0, 3) = 1.0;
          b(0) = hv(ix, iy);

          xx = (ix + 1) * hx;
          yy = iy * hy;
          H(1, 0) = xx;
          H(1, 1) = yy;
          H(1, 2) = xx * yy;
          H(1, 3) = 1.0;
          b(1) = hv(ix + 1, iy);

          xx = (ix + 1) * hx;
          yy = (iy + 1) * hy;
          H(2, 0) = xx;
          H(2, 1) = yy;
          H(2, 2) = xx * yy;
          H(2, 3) = 1.0;
          b(2) = hv(ix + 1, iy + 1);

          xx = ix * hx;
          yy = (iy + 1) * hy;
          H(3, 0) = xx;
          H(3, 1) = yy;
          H(3, 2) = xx * yy;
          H(3, 3) = 1.0;
          b(3) = hv(ix, iy + 1);

          H.gauss_jordan();
          H.vmult(X, b);

          val = X(0) * x + X(1) * y + X(2) * x * y + X(3);
        }

      return val;
    }

// As mentioned above this function reads in the
// obstacle data and stores them in the std::vector
// obstacle_data. It will be used only in run ().
  template <int dim>
    void
    Input<dim>::read_obstacle (const std::string name)
    {
      std::ifstream f(name.c_str());

      std::string temp;
      f >> temp >> nx >> ny;
      assert(nx > 0 && ny > 0);

      for (int k = 0; k < nx * ny; k++)
        {
          double val;
          f >> val;
          obstacle_data.push_back(val);
        }

      hx = 1.0 / (nx - 1);
      hy = 1.0 / (ny - 1);

      pcout << "Resolution of the scanned obstacle picture: " << nx << " x "
          << ny << std::endl;
    }

// @sect3{The <code>ConstitutiveLaw</code> class template}

// This class provides an interface
// for a constitutive law. In this
// example we are using an elastoplastic
// material behavior with linear,
// isotropic hardening.
// For gamma = 0 we obtain perfect elastoplastic
// behavior.
  template <int dim>
    class ConstitutiveLaw
    {
      public:
        ConstitutiveLaw (
            double _E, double _nu, double _sigma_0, double _gamma,
            MPI_Comm _mpi_communicator, ConditionalOStream _pcout);

        void
        plast_linear_hardening (
            SymmetricTensor<4, dim> &stress_strain_tensor,
            const SymmetricTensor<2, dim> &strain_tensor,
            unsigned int &elast_points, unsigned int &plast_points,
            double &yield);
        void
        linearized_plast_linear_hardening (
            SymmetricTensor<4, dim> &stress_strain_tensor_linearized,
            SymmetricTensor<4, dim> &stress_strain_tensor,
            const SymmetricTensor<2, dim> &strain_tensor);
        inline SymmetricTensor<2, dim>
        get_strain (
            const FEValues<dim> &fe_values, const unsigned int shape_func,
            const unsigned int q_point) const;
        void
        set_sigma_0 (
            double sigma_hlp)
        {
          sigma_0 = sigma_hlp;
        }

      private:
        SymmetricTensor<4, dim> stress_strain_tensor_mu;
        SymmetricTensor<4, dim> stress_strain_tensor_kappa;
        double E;
        double nu;
        double sigma_0;
        double gamma;
        double mu;
        double kappa;
        MPI_Comm mpi_communicator;
        ConditionalOStream pcout;
    };

// The constructor of the ConstitutiveLaw class sets the
// required material parameter for our deformable body:
// E -> elastic modulus
// nu -> Passion's number
// sigma_0 -> yield stress
// gamma -> hardening parameter.
// Also it supplies the stress strain tensor of forth order
// of the volumetric and deviator part. For further details
// see the documentation above.
  template <int dim>
    ConstitutiveLaw<dim>::ConstitutiveLaw (
        double _E, double _nu, double _sigma_0, double _gamma,
        MPI_Comm _mpi_communicator, ConditionalOStream _pcout)
        :
            E(_E),
            nu(_nu),
            sigma_0(_sigma_0),
            gamma(_gamma),
            mpi_communicator(_mpi_communicator),
            pcout(_pcout)
    {
      mu = E / (2 * (1 + nu));
      kappa = E / (3 * (1 - 2 * nu));
      stress_strain_tensor_kappa = kappa
          * outer_product(unit_symmetric_tensor<dim>(),
              unit_symmetric_tensor<dim>());
      stress_strain_tensor_mu = 2 * mu
          * (identity_tensor<dim>()
              - outer_product(unit_symmetric_tensor<dim>(),
                  unit_symmetric_tensor<dim>()) / 3.0);
    }

// @sect3{ConstitutiveLaw::ConstitutiveLaw}

// Calculates the strain $\varepsilon(\varphi)=\dfrac{1}{2}\left(\nabla\varphi + \nabla\varphi^T$
// for the shape functions $\varphi$.
  template <int dim>
    inline SymmetricTensor<2, dim>
    ConstitutiveLaw<dim>::get_strain (
        const FEValues<dim> &fe_values, const unsigned int shape_func,
        const unsigned int q_point) const
    {
      const FEValuesExtractors::Vector displacement(0);
      SymmetricTensor<2, dim> tmp;

      tmp = fe_values[displacement].symmetric_gradient(shape_func, q_point);

      return tmp;
    }

// @sect3{ConstitutiveLaw::plast_linear_hardening}

// This is the implemented constitutive law. It projects the
// deviator part of the stresses in a quadrature point back to
// the yield stress plus the linear isotropic hardening.
// Also we sum up the elastic and the plastic quadrature
// points. We need this function to calculate the nonlinear
// residual in
// PlasticityContactProblem::residual_nl_system(TrilinosWrappers::MPI::Vector &u).
  template <int dim>
    void
    ConstitutiveLaw<dim>::plast_linear_hardening (
        SymmetricTensor<4, dim> &stress_strain_tensor,
        const SymmetricTensor<2, dim> &strain_tensor,
        unsigned int &elast_points, unsigned int &plast_points, double &yield)
    {
      if (dim == 3)
        {
          SymmetricTensor<2, dim> stress_tensor;
          stress_tensor = (stress_strain_tensor_kappa + stress_strain_tensor_mu)
              * strain_tensor;

          SymmetricTensor<2, dim> deviator_stress_tensor = deviator(
              stress_tensor);

          double deviator_stress_tensor_norm = deviator_stress_tensor.norm();

          yield = 0;
          stress_strain_tensor = stress_strain_tensor_mu;
          double beta = 1.0;
          if (deviator_stress_tensor_norm > sigma_0)
            {
              beta = sigma_0 / deviator_stress_tensor_norm;
              stress_strain_tensor *= (gamma + (1 - gamma) * beta);
              yield = 1;
              plast_points += 1;
            }
          else
            elast_points += 1;

          stress_strain_tensor += stress_strain_tensor_kappa;
        }
    }

// @sect3{ConstitutiveLaw::linearized_plast_linear_hardening}

// This function returns the linearized stress strain tensor
// in the solution $u^{i-1}$ of the previous Newton $i-1$ step.
// The parameter strain_tensor $\varepsilon(u^{i-1})$ is calculated
// by $u^{i-1}$. It contains the derivative of the nonlinear
// constitutive law. As the result this function returns
// the stress_strain_tensor of the nonlinear problem as well as
// the stress_strain_tensor_linearized of the linearized problem.
// See
// PlasticityContactProblem::assemble_nl_system(TrilinosWrappers::MPI::Vector &u)
// where this function is used.
  template <int dim>
    void
    ConstitutiveLaw<dim>::linearized_plast_linear_hardening (
        SymmetricTensor<4, dim> &stress_strain_tensor_linearized,
        SymmetricTensor<4, dim> &stress_strain_tensor,
        const SymmetricTensor<2, dim> &strain_tensor)
    {
      if (dim == 3)
        {
          SymmetricTensor<2, dim> stress_tensor;
          stress_tensor = (stress_strain_tensor_kappa + stress_strain_tensor_mu)
              * strain_tensor;

          SymmetricTensor<2, dim> deviator_stress_tensor = deviator(
              stress_tensor);

          double deviator_stress_tensor_norm = deviator_stress_tensor.norm();

          stress_strain_tensor = stress_strain_tensor_mu;
          stress_strain_tensor_linearized = stress_strain_tensor_mu;
          double beta = 1.0;
          if (deviator_stress_tensor_norm > sigma_0)
            {
              beta = sigma_0 / deviator_stress_tensor_norm;
              stress_strain_tensor *= (gamma + (1 - gamma) * beta);
              stress_strain_tensor_linearized *= (gamma + (1 - gamma) * beta);
              deviator_stress_tensor /= deviator_stress_tensor_norm;
              stress_strain_tensor_linearized -= (1 - gamma) * beta * 2 * mu
                  * outer_product(deviator_stress_tensor,
                      deviator_stress_tensor);
            }

          stress_strain_tensor += stress_strain_tensor_kappa;
          stress_strain_tensor_linearized += stress_strain_tensor_kappa;
        }
    }

  namespace EquationData
  {
// It possible to apply an additional body force
// but in here it is set to zero.
    template <int dim>
      class RightHandSide : public Function<dim>
      {
        public:
          RightHandSide ()
              :
                  Function<dim>(dim)
          {
          }

          virtual double
          value (
              const Point<dim> &p, const unsigned int component = 0) const;

          virtual void
          vector_value (
              const Point<dim> &p, Vector<double> &values) const;
      };

    template <int dim>
      double
      RightHandSide<dim>::value (
          const Point<dim> &p, const unsigned int component) const
      {
        double return_value = 0.0;

        if (component == 0)
          return_value = 0.0;
        if (component == 1)
          return_value = 0.0;
        if (component == 2)
          return_value = 0.0;

        return return_value;
      }

    template <int dim>
      void
      RightHandSide<dim>::vector_value (
          const Point<dim> &p, Vector<double> &values) const
      {
        for (unsigned int c = 0; c < this->n_components; ++c)
          values(c) = RightHandSide<dim>::value(p, c);
      }

// This function class is used to describe the prescribed displacements
// at the boundary. But again we set this to zero.
    template <int dim>
      class BoundaryValues : public Function<dim>
      {
        public:
          BoundaryValues ()
              :
                  Function<dim>(dim)
          {
          }
          ;

          virtual double
          value (
              const Point<dim> &p, const unsigned int component = 0) const;

          virtual void
          vector_value (
              const Point<dim> &p, Vector<double> &values) const;
      };

    template <int dim>
      double
      BoundaryValues<dim>::value (
          const Point<dim> &p, const unsigned int component) const
      {
        double return_value = 0;

        if (component == 0)
          return_value = 0.0;
        if (component == 1)
          return_value = 0.0;
        if (component == 2)
          return_value = 0.0;

        return return_value;
      }

    template <int dim>
      void
      BoundaryValues<dim>::vector_value (
          const Point<dim> &p, Vector<double> &values) const
      {
        for (unsigned int c = 0; c < this->n_components; ++c)
          values(c) = BoundaryValues<dim>::value(p, c);
      }

// This function is obviously implemented to
// define the obstacle that penetrates our deformable
// body. You can choose between two ways to define
// your obstacle: to read it from a file or to use
// a function (here a ball).
// z_max_domain is the z value of the surface of the work piece  
    template <int dim>
      class Obstacle : public Function<dim>
      {
        public:
          Obstacle (
              std_cxx1x::shared_ptr<Input<dim> > const &_input,
              bool _use_read_obstacle, double z_max_domain)
              :
                  Function<dim>(dim),
                  input_obstacle_copy(_input),
                  use_read_obstacle(_use_read_obstacle),
                  z_max_domain(z_max_domain)
          {
          }

          virtual double
          value (
              const Point<dim> &p, const unsigned int component = 0) const;

          virtual void
          vector_value (
              const Point<dim> &p, Vector<double> &values) const;

        private:
          std_cxx1x::shared_ptr<Input<dim> > const &input_obstacle_copy;
          bool use_read_obstacle;
          double z_max_domain;
      };

    template <int dim>
      double
      Obstacle<dim>::value (
          const Point<dim> &p, const unsigned int component) const
      {
        if (component == 0)
          return p(0);
        if (component == 1)
          return p(1);

        //component==2:
        if (use_read_obstacle)
          {
            if (p(0) >= 0.0 && p(0) <= 1.0 && p(1) >= 0.0 && p(1) <= 1.0)
              return z_max_domain + 0.999
                  - input_obstacle_copy->obstacle_function(p(0), p(1));
            else
              return 10000.0;
          }
        else
          {
            //sphere:
            return -std::sqrt(
                0.36 - (p(0) - 0.5) * (p(0) - 0.5)
                    - (p(1) - 0.5) * (p(1) - 0.5)) + z_max_domain + 0.59;
          }
      }

    template <int dim>
      void
      Obstacle<dim>::vector_value (
          const Point<dim> &p, Vector<double> &values) const
      {
        for (unsigned int c = 0; c < this->n_components; ++c)
          values(c) = Obstacle<dim>::value(p, c);
      }
  }

// @sect3{The <code>PlasticityContactProblem</code> class template}

// This class supplies all function
// and variables needed to describe
// the nonlinear contact problem. It is
// close to step-41 but with some additional
// features like: handling hanging nodes,
// a newton method, using Trilinos and p4est
// for parallel distributed computing.
// To deal with hanging nodes makes
// life a bit more complicated since
// we need an other ConstraintMatrix now.
// We create a newton method for the
// active set method for the contact
// situation and to handle the nonlinear
// operator for the constitutive law.

  template <int dim>
    class PlasticityContactProblem
    {
      public:
        PlasticityContactProblem (
            const ParameterHandler &prm);
        void
        run ();

        static void
        declare (
            ParameterHandler &prm);

      private:
        void
        make_grid ();
        void
        setup_system ();
        void
        assemble_nl_system (
            TrilinosWrappers::MPI::Vector &u);
        void
        residual_nl_system (
            TrilinosWrappers::MPI::Vector &u);
        void
        assemble_mass_matrix_diagonal (
            TrilinosWrappers::SparseMatrix &mass_matrix);
        void
        update_solution_and_constraints ();
        void
        dirichlet_constraints ();
        void
        solve ();
        void
        solve_newton ();
        void
        refine_grid ();
        void
        move_mesh (
            const TrilinosWrappers::MPI::Vector &_complete_displacement) const;
        void
        output_results (
            const std::string &title);
        void
        output_contact_force (
            const unsigned int cycle);

        double to_refine_factor;
        double to_coarsen_factor;
        unsigned int cycle;

        MPI_Comm mpi_communicator;

        parallel::distributed::Triangulation<dim> triangulation;

        FE_Q<dim> u;
        FESystem<dim> fe;
        DoFHandler<dim> dof_handler;

        // We are using the SolutionTransfer class to interpolate the
        // solution on the new refined mesh. It appears in th refine_grid()
        // and the run() function.
        std_cxx1x::shared_ptr<
            parallel::distributed::SolutionTransfer<dim,
                TrilinosWrappers::MPI::Vector> > soltrans;

        IndexSet locally_owned_dofs;
        IndexSet locally_relevant_dofs;

        unsigned int number_iterations;

        ConstraintMatrix constraints;
        ConstraintMatrix constraints_hanging_nodes;
        ConstraintMatrix constraints_dirichlet_hanging_nodes;

        TrilinosWrappers::SparseMatrix system_matrix_newton;

        TrilinosWrappers::MPI::Vector solution;
        TrilinosWrappers::MPI::Vector system_rhs_newton;
        TrilinosWrappers::MPI::Vector system_rhs_lambda;
        TrilinosWrappers::MPI::Vector resid_vector;
        TrilinosWrappers::MPI::Vector diag_mass_matrix_vector;
        Vector<float> cell_constitution;
        IndexSet active_set;

        ConditionalOStream pcout;

        TrilinosWrappers::PreconditionAMG::AdditionalData additional_data;
        TrilinosWrappers::PreconditionAMG preconditioner_u;

        std_cxx1x::shared_ptr<Input<dim> > input_obstacle;
        std_cxx1x::shared_ptr<ConstitutiveLaw<dim> > plast_lin_hard;

        double sigma_0; // Yield stress
        double gamma; // Parameter for the linear isotropic hardening
        double e_modul; // E-Modul
        double nu; // Poisson ratio

        TimerOutput computing_timer;

        unsigned int degree;
        unsigned int n_initial_refinements;
        struct RefinementStrategy
        {
            enum value
            {
              refine_global,
              refine_percentage,
              refine_fix_dofs
            };
        };
        typename RefinementStrategy::value refinement_strategy;
        unsigned int n_cycles;
        std::string obstacle_filename;
        std::string output_dir;
        bool transfer_solution;
        std::string base_mesh;
    };

// @sect3{Implementation of the <code>PlasticityContactProblem</code> class}

// Next for the implementation of the class
// template that makes use of the functions
// above. As before, we will write everything

  template <int dim>
    PlasticityContactProblem<dim>::PlasticityContactProblem (
        const ParameterHandler &prm)
        :
            mpi_communicator(MPI_COMM_WORLD),
            triangulation(mpi_communicator),
            u(QGaussLobatto<1>(prm.get_integer("polynomial degree") + 1)),
            fe(u, dim),
            dof_handler(triangulation),
            pcout(std::cout,
                (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
            sigma_0(400.0),
            gamma(0.01),
            e_modul(2.0e+5),
            nu(0.3),
            computing_timer(MPI_COMM_WORLD, pcout, TimerOutput::never,
                TimerOutput::wall_times)
    {
      // double _E, double _nu, double _sigma_0, double _gamma
      plast_lin_hard.reset(
          new ConstitutiveLaw<dim>(e_modul, nu, sigma_0, gamma,
              mpi_communicator, pcout));

      degree = prm.get_integer("polynomial degree");
      n_initial_refinements = prm.get_integer("number of initial refinements");
      std::string strat = prm.get("refinement strategy");
      if (strat == "global")
        refinement_strategy = RefinementStrategy::refine_global;
      else if (strat == "percentage")
        refinement_strategy = RefinementStrategy::refine_percentage;
      else
        throw ExcNotImplemented();

      n_cycles = prm.get_integer("number of cycles");
      obstacle_filename = prm.get("obstacle filename");
      output_dir = prm.get("output directory");
      if (output_dir != "" && *(output_dir.rbegin()) != '/')
        output_dir += "/";
      mkdir(output_dir.c_str(), 0777);

      transfer_solution = prm.get_bool("transfer solution");
      base_mesh = prm.get("base mesh");

      pcout << "    Using output directory '" << output_dir << "'" << std::endl;
      pcout << "    FE degree " << degree << std::endl;
      pcout << "    Obstacle '" << obstacle_filename << "'" << std::endl;
      pcout << "    transfer solution "
          << (transfer_solution ? "true" : "false") << std::endl;
    }

// @sect4{PlasticityContactProblem::declare}

  template <int dim>
    void
    PlasticityContactProblem<dim>::declare (
        ParameterHandler &prm)
    {
      prm.declare_entry("polynomial degree", "1", Patterns::Integer(),
          "polynomial degree of the FE_Q finite element space, typically 1 or 2");
      prm.declare_entry("number of initial refinements", "2",
          Patterns::Integer(),
          "number of initial global refinements before the first computation");
      prm.declare_entry("refinement strategy", "percentage",
          Patterns::Selection("global|percentage|fix dofs"),
          "refinement strategy for each cycle:\n"
              " global: one global refinement\n"
              "percentage: fixed percentage gets refined using kelly\n"
              " fix dofs: tries to achieve 2^initial_refinement*300 dofs after cycle 1 (only use 2 cycles!). Changes the coarse mesh!");
      prm.declare_entry("number of cycles", "5", Patterns::Integer(),
          "number of adaptive cycles to run");
      prm.declare_entry("obstacle filename", "", Patterns::Anything(),
          "obstacle file to read, use 'obstacle_file.pbm' or leave empty to use a sphere");
      prm.declare_entry("output directory", "", Patterns::Anything(),
          "directory to put output files (graphical output and benchmark statistics), leave empty to put into current directory");
      prm.declare_entry("transfer solution", "false", Patterns::Bool(),
          "decide if the solution should be used as a starting guess for the finer mesh, use 0 otherwise.");
      prm.declare_entry("base mesh", "box",
          Patterns::Selection("box|half sphere"),
          "select the shape of the work piece: 'box' or 'half sphere'");

    }

  Point<3>
  rotate_half_sphere (
      const Point<3> &in)
  {
    return Point<3>(in(2), in(1), -in(0));
  }

// @sect4{PlasticityContactProblem::make_grid}

  template <int dim>
    void
    PlasticityContactProblem<dim>::make_grid ()
    {

      if (base_mesh == "half sphere")
        {
          Point<dim> center(0, 0, 0);
          double radius = 0.8;
          GridGenerator::half_hyper_ball(triangulation, center, radius);
          GridTools::transform(&rotate_half_sphere, triangulation);
          Point<dim> shift(0.5, 0.5, 0.5);
          GridTools::shift(shift, triangulation);
          static HyperBallBoundary<dim> boundary_description(
              Point<dim>(0.5, 0.5, 0.5), radius);
          triangulation.set_boundary(0, boundary_description);

          triangulation.refine_global(n_initial_refinements);

          to_refine_factor = 0.3;
          to_coarsen_factor = 0.03;
          return;
        }

      Point<dim> p1(0, 0, 0);
      Point<dim> p2(1.0, 1.0, 1.0);

      GridGenerator::hyper_rectangle(triangulation, p1, p2);
      to_refine_factor = 0.3;
      to_coarsen_factor = 0.03;

      Triangulation<3>::active_cell_iterator cell =
          triangulation.begin_active(), endc = triangulation.end();

      /* boundary_indicators:
           _______
         /  1    /|
        /______ / |
       |       | 8|
       |   8   | /
       |_______|/
           6

       The boundary indicators of the sides of the cube are 8.
       The boundary indicator of the bottom is indicated with 6
       and the top with 1.
       */

      for (; cell != endc; ++cell)
        for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
            ++face)
          {
            if (cell->face(face)->center()[2] == p2(2))
              cell->face(face)->set_boundary_indicator(1);
            if (cell->face(face)->center()[0] == p1(0)
                || cell->face(face)->center()[0] == p2(0)
                || cell->face(face)->center()[1] == p1(1)
                || cell->face(face)->center()[1] == p2(1))
              cell->face(face)->set_boundary_indicator(8);
            if (cell->face(face)->center()[2] == p1(2))
              cell->face(face)->set_boundary_indicator(6);
          }

      triangulation.refine_global(n_initial_refinements);
    }

  template <int dim>
    void
    PlasticityContactProblem<dim>::setup_system ()
    {
      // setup dofs
        {
          TimerOutput::Scope t(computing_timer, "Setup: distribute DoFs");
          dof_handler.distribute_dofs(fe);

          locally_owned_dofs = dof_handler.locally_owned_dofs();
          locally_relevant_dofs.clear();
          DoFTools::extract_locally_relevant_dofs(dof_handler,
              locally_relevant_dofs);
        }

      // setup hanging nodes and Dirichlet constraints
        {
          TimerOutput::Scope t(computing_timer, "Setup: constraints");
          constraints_hanging_nodes.reinit(locally_relevant_dofs);
          DoFTools::make_hanging_node_constraints(dof_handler,
              constraints_hanging_nodes);
          constraints_hanging_nodes.close();

          pcout << "   Number of active cells: "
              << triangulation.n_global_active_cells() << std::endl
              << "   Number of degrees of freedom: " << dof_handler.n_dofs()
              << std::endl;

          dirichlet_constraints();
        }

      // Initialization for matrices and vectors
        {
          TimerOutput::Scope t(computing_timer, "Setup: vectors");
          solution.reinit(locally_relevant_dofs, mpi_communicator);
          system_rhs_newton.reinit(locally_owned_dofs, mpi_communicator);
          system_rhs_lambda.reinit(system_rhs_newton);
          resid_vector.reinit(system_rhs_newton);
          diag_mass_matrix_vector.reinit(system_rhs_newton);
          cell_constitution.reinit(triangulation.n_active_cells());
          active_set.clear();
          active_set.set_size(locally_relevant_dofs.size());
        }

      // setup sparsity pattern
        {
          TimerOutput::Scope t(computing_timer, "Setup: matrix");
          TrilinosWrappers::SparsityPattern sp(locally_owned_dofs,
              mpi_communicator);

          DoFTools::make_sparsity_pattern(dof_handler, sp,
              constraints_dirichlet_hanging_nodes, false,
              Utilities::MPI::this_mpi_process(mpi_communicator));

          sp.compress();

          system_matrix_newton.reinit(sp);

          // we are going to reuse the system
          // matrix for assembling the diagonal
          // of the mass matrix so that we do not
          // need to allocate two sparse matrices
          // at the same time:
          TrilinosWrappers::SparseMatrix & mass_matrix = system_matrix_newton;
          assemble_mass_matrix_diagonal(mass_matrix);
          const unsigned int start = (system_rhs_newton.local_range().first),
              end = (system_rhs_newton.local_range().second);
          for (unsigned int j = start; j < end; j++)
            diag_mass_matrix_vector(j) = mass_matrix.diag_element(j);

          number_iterations = 0;

          diag_mass_matrix_vector.compress(VectorOperation::insert);

          // remove the mass matrix entries from the matrix:
          mass_matrix = 0;
        }
    }

  template <int dim>
    void
    PlasticityContactProblem<dim>::assemble_nl_system (
        TrilinosWrappers::MPI::Vector &u)
    {
      TimerOutput::Scope t(computing_timer, "Assembling");

      QGauss<dim> quadrature_formula(fe.degree + 1);
      QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

      FEValues<dim> fe_values(fe, quadrature_formula,
          UpdateFlags(
              update_values | update_gradients | update_q_points
                  | update_JxW_values));

      FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
          update_values | update_quadrature_points | update_JxW_values);

      const unsigned int dofs_per_cell = fe.dofs_per_cell;
      const unsigned int n_q_points = quadrature_formula.size();
      const unsigned int n_face_q_points = face_quadrature_formula.size();

      const EquationData::RightHandSide<dim> right_hand_side;
      std::vector<Vector<double> > right_hand_side_values(n_q_points,
          Vector<double>(dim));
      std::vector<Vector<double> > right_hand_side_values_face(n_face_q_points,
          Vector<double>(dim));

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Vector<double> cell_rhs(dofs_per_cell);

      std::vector<unsigned int> local_dof_indices(dofs_per_cell);

      typename DoFHandler<dim>::active_cell_iterator cell =
          dof_handler.begin_active(), endc = dof_handler.end();

      const FEValuesExtractors::Vector displacement(0);

      const double kappa = 1.0;
      for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            cell_matrix = 0;
            cell_rhs = 0;

            right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                right_hand_side_values);

            std::vector<SymmetricTensor<2, dim> > strain_tensor(n_q_points);
            fe_values[displacement].get_function_symmetric_gradients(u,
                strain_tensor);

            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
              {
                SymmetricTensor<4, dim> stress_strain_tensor_linearized;
                SymmetricTensor<4, dim> stress_strain_tensor;
                SymmetricTensor<2, dim> stress_tensor;

                plast_lin_hard->linearized_plast_linear_hardening(
                    stress_strain_tensor_linearized, stress_strain_tensor,
                    strain_tensor[q_point]);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    stress_tensor = stress_strain_tensor_linearized
                        * plast_lin_hard->get_strain(fe_values, i, q_point);

                    for (unsigned int j = 0; j < dofs_per_cell; ++j)
                      {
                        cell_matrix(i, j) += (stress_tensor
                            * plast_lin_hard->get_strain(fe_values, j, q_point)
                            * fe_values.JxW(q_point));
                      }

                    // the linearized part a(v^i;v^i,v) of the rhs
                    cell_rhs(i) += (stress_tensor * strain_tensor[q_point]
                        * fe_values.JxW(q_point));

                    // the residual part a(v^i;v) of the rhs
                    cell_rhs(i) -= (strain_tensor[q_point]
                        * stress_strain_tensor
                        * plast_lin_hard->get_strain(fe_values, i, q_point)
                        * fe_values.JxW(q_point));

                    // the residual part F(v) of the rhs
                    Tensor<1, dim> rhs_values;
                    rhs_values = 0;
                    cell_rhs(i) += (fe_values[displacement].value(i, q_point)
                        * rhs_values * fe_values.JxW(q_point));
                  }
              }

            for (unsigned int face = 0;
                face < GeometryInfo<dim>::faces_per_cell; ++face)
              {
                if (cell->face(face)->at_boundary()
                    && cell->face(face)->boundary_indicator() == 1)
                  {
                    fe_values_face.reinit(cell, face);

                    right_hand_side.vector_value_list(
                        fe_values_face.get_quadrature_points(),
                        right_hand_side_values_face);

                    for (unsigned int q_point = 0; q_point < n_face_q_points;
                        ++q_point)
                      {
                        Tensor<1, dim> rhs_values;
                        rhs_values[2] = right_hand_side_values[q_point][2];
                        for (unsigned int i = 0; i < dofs_per_cell; ++i)
                          cell_rhs(i) += (fe_values_face[displacement].value(i,
                              q_point) * rhs_values
                              * fe_values_face.JxW(q_point));
                      }
                  }
              }

            cell->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(cell_matrix, cell_rhs,
                local_dof_indices, system_matrix_newton, system_rhs_newton,
                true);

          };

      system_matrix_newton.compress(VectorOperation::add);
      system_rhs_newton.compress(VectorOperation::add);
    }

  template <int dim>
    void
    PlasticityContactProblem<dim>::residual_nl_system (
        TrilinosWrappers::MPI::Vector &u)
    {
      QGauss<dim> quadrature_formula(fe.degree + 1);
      QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

      FEValues<dim> fe_values(fe, quadrature_formula,
          UpdateFlags(
              update_values | update_gradients | update_q_points
                  | update_JxW_values));

      FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
          update_values | update_quadrature_points | update_JxW_values);

      const unsigned int dofs_per_cell = fe.dofs_per_cell;
      const unsigned int n_q_points = quadrature_formula.size();
      const unsigned int n_face_q_points = face_quadrature_formula.size();

      const EquationData::RightHandSide<dim> right_hand_side;
      std::vector<Vector<double> > right_hand_side_values(n_q_points,
          Vector<double>(dim));
      std::vector<Vector<double> > right_hand_side_values_face(n_face_q_points,
          Vector<double>(dim));

      Vector<double> cell_rhs(dofs_per_cell);

      std::vector<unsigned int> local_dof_indices(dofs_per_cell);

      const FEValuesExtractors::Vector displacement(0);

      typename DoFHandler<dim>::active_cell_iterator cell =
          dof_handler.begin_active(), endc = dof_handler.end();

      unsigned int elast_points = 0;
      unsigned int plast_points = 0;
      double yield = 0;
      unsigned int cell_number = 0;
      cell_constitution = 0;

      for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
          {
            fe_values.reinit(cell);
            cell_rhs = 0;

            right_hand_side.vector_value_list(fe_values.get_quadrature_points(),
                right_hand_side_values);

            std::vector<SymmetricTensor<2, dim> > strain_tensor(n_q_points);
            fe_values[displacement].get_function_symmetric_gradients(u,
                strain_tensor);

            for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
              {
                SymmetricTensor<4, dim> stress_strain_tensor;
                SymmetricTensor<2, dim> stress_tensor;

                plast_lin_hard->plast_linear_hardening(stress_strain_tensor,
                    strain_tensor[q_point], elast_points, plast_points, yield);

                cell_constitution(cell_number) += yield;
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    cell_rhs(i) -= (strain_tensor[q_point]
                        * stress_strain_tensor
                        * //(stress_tensor) *
                        plast_lin_hard->get_strain(fe_values, i, q_point)
                        * fe_values.JxW(q_point));

                    Tensor<1, dim> rhs_values;
                    rhs_values = 0;
                    cell_rhs(i) += ((fe_values[displacement].value(i, q_point)
                        * rhs_values) * fe_values.JxW(q_point));
                  };
              };

            for (unsigned int face = 0;
                face < GeometryInfo<dim>::faces_per_cell; ++face)
              {
                if (cell->face(face)->at_boundary()
                    && cell->face(face)->boundary_indicator() == 1)
                  {
                    fe_values_face.reinit(cell, face);

                    right_hand_side.vector_value_list(
                        fe_values_face.get_quadrature_points(),
                        right_hand_side_values_face);

                    for (unsigned int q_point = 0; q_point < n_face_q_points;
                        ++q_point)
                      {
                        Tensor<1, dim> rhs_values;
                        rhs_values[2] = right_hand_side_values[q_point][2];
                        for (unsigned int i = 0; i < dofs_per_cell; ++i)
                          cell_rhs(i) += (fe_values_face[displacement].value(i,
                              q_point) * rhs_values
                              * fe_values_face.JxW(q_point));
                      }
                  }
              }

            cell->get_dof_indices(local_dof_indices);
            constraints_dirichlet_hanging_nodes.distribute_local_to_global(
                cell_rhs, local_dof_indices, system_rhs_newton);

            for (unsigned int i = 0; i < dofs_per_cell; i++)
              system_rhs_lambda(local_dof_indices[i]) += cell_rhs(i);

            cell_number += 1;
          }
        else
          {
            cell_constitution(cell_number) = 0;
            cell_number += 1;
          };

      cell_constitution /= n_q_points;
      cell_constitution.compress(VectorOperation::add);
      system_rhs_newton.compress(VectorOperation::add);
      system_rhs_lambda.compress(VectorOperation::add);

//  	constraints_hanging_nodes.condense(system_rhs_lambda);

      unsigned int sum_elast_points = Utilities::MPI::sum(elast_points,
          mpi_communicator);
      unsigned int sum_plast_points = Utilities::MPI::sum(plast_points,
          mpi_communicator);
      pcout << "      Number of elastic quadrature points: " << sum_elast_points
          << " and plastic quadrature points: " << sum_plast_points
          << std::endl;
    }

  template <int dim>
    void
    PlasticityContactProblem<dim>::assemble_mass_matrix_diagonal (
        TrilinosWrappers::SparseMatrix &mass_matrix)
    {
      QGaussLobatto<dim - 1> face_quadrature_formula(fe.degree + 1);

      FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
          update_values | update_quadrature_points | update_JxW_values);

      const unsigned int dofs_per_cell = fe.dofs_per_cell;
      const unsigned int n_face_q_points = face_quadrature_formula.size();

      FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
      Tensor<1, dim, double> ones(dim);
      for (unsigned i = 0; i < dim; i++)
        ones[i] = 1.0;

      std::vector<unsigned int> local_dof_indices(dofs_per_cell);

      const FEValuesExtractors::Vector displacement(0);

      typename DoFHandler<dim>::active_cell_iterator cell =
          dof_handler.begin_active(), endc = dof_handler.end();

      for (; cell != endc; ++cell)
        if (cell->is_locally_owned())
          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
              ++face)
            if (cell->face(face)->at_boundary()
                && cell->face(face)->boundary_indicator() == 1)
              {
                fe_values_face.reinit(cell, face);
                cell_matrix = 0;

                for (unsigned int q_point = 0; q_point < n_face_q_points;
                    ++q_point)
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    cell_matrix(i, i) += (fe_values_face[displacement].value(i,
                        q_point) * ones * fe_values_face.JxW(q_point));

                cell->get_dof_indices(local_dof_indices);

// 					constraints_dirichlet_hanging_nodes.distribute_local_to_global(
// 							cell_matrix, local_dof_indices, mass_matrix);

                for (unsigned int i = 0; i < dofs_per_cell; i++)
                  mass_matrix.add(local_dof_indices[i], local_dof_indices[i],
                      cell_matrix(i, i));
              }
      mass_matrix.compress(VectorOperation::add);
    }

// @sect4{PlasticityContactProblem::update_solution_and_constraints}

// Projection and updating of the active set
// for the dofs which penetrates the obstacle.
  template <int dim>
    void
    PlasticityContactProblem<dim>::update_solution_and_constraints ()
    {
      const EquationData::Obstacle<dim> obstacle(input_obstacle,
          (obstacle_filename != ""), (base_mesh == "box" ? 1.0 : 0.5));
      std::vector<bool> vertex_touched(dof_handler.n_dofs(), false);

      typename DoFHandler<dim>::active_cell_iterator cell =
          dof_handler.begin_active(), endc = dof_handler.end();

      TrilinosWrappers::MPI::Vector distributed_solution(system_rhs_newton);
      distributed_solution = solution;
      TrilinosWrappers::MPI::Vector lambda(solution);
      lambda = resid_vector;
      TrilinosWrappers::MPI::Vector diag_mass_matrix_vector_relevant(solution);
      diag_mass_matrix_vector_relevant = diag_mass_matrix_vector;

      constraints.reinit(locally_relevant_dofs);
      active_set.clear();
      IndexSet active_set_locally_owned;
      active_set_locally_owned.set_size(locally_owned_dofs.size());
      const double c = 100.0 * e_modul;

      Quadrature<dim - 1> face_quadrature(fe.get_unit_face_support_points());
      FEFaceValues<dim> fe_values_face(fe, face_quadrature,
          update_quadrature_points);

      const unsigned int dofs_per_face = fe.dofs_per_face;
      const unsigned int n_face_q_points = face_quadrature.size();

      // pcout<< "dofs_per_face = " << dofs_per_face
      //      << "n_face_q_points = " << n_face_q_points
      //      <<std::endl;
      unsigned int counter_hanging_nodes = 0;
      for (; cell != endc; ++cell)
        if (!cell->is_artificial())
          for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
              ++face)
            if (cell->face(face)->at_boundary()
                && cell->face(face)->boundary_indicator() == 1)
              {
                fe_values_face.reinit(cell, face);
                std::vector<unsigned int> dof_indices(dofs_per_face);
                cell->face(face)->get_dof_indices(dof_indices);

                for (unsigned int q_point = 0; q_point < n_face_q_points;
                    ++q_point)
                  {
                    unsigned int component = fe.face_system_to_component_index(
                        q_point).first;

                    if (component == 2)
                      {
                        unsigned int index_z = dof_indices[q_point];

                        if (vertex_touched[index_z] == false)
                          vertex_touched[index_z] = true;
                        else
                          continue;

                        // the local row where
                        Point<dim> point(
                            fe_values_face.quadrature_point(q_point));

                        double obstacle_value = obstacle.value(point, 2);
                        double solution_index_z = solution(index_z);
                        double gap = obstacle_value - point(2);

                        if (lambda(index_z)
                            / diag_mass_matrix_vector_relevant(index_z)
                            + c * (solution_index_z - gap) > 0
                            && !(constraints_hanging_nodes.is_constrained(
                                index_z)))
                          {
                            constraints.add_line(index_z);
                            constraints.set_inhomogeneity(index_z, gap);
                            distributed_solution(index_z) = gap;

                            if (locally_owned_dofs.is_element(index_z))
                              {
                                active_set_locally_owned.add_index(index_z);
                                if (locally_relevant_dofs.is_element(index_z))
                                  active_set.add_index(index_z);
                              }

                          }
                        else if (lambda(index_z)
                            / diag_mass_matrix_vector_relevant(index_z)
                            + c * (solution_index_z - gap) > 0
                            && constraints_hanging_nodes.is_constrained(
                                index_z))
                          {
                            if (locally_owned_dofs.is_element(index_z))
                              {
                                counter_hanging_nodes += 1;

// 							std::cout << "index_z = " << index_z
// 								<< ", lambda = " << lambda (index_z)
//  								<< ", solution_index_z - gap = " << solution_index_z - gap
// 								<< ", diag_mass_matrix_vector_relevant = " << diag_mass_matrix_vector_relevant (index_z)
//  								<< ", x = " << point(0)
//  								<< ", y = " << point(1)
//  								<< std::endl;
                              }
                          }
                      }
                  }
              }
      distributed_solution.compress(VectorOperation::insert);

      unsigned int sum_contact_constraints = Utilities::MPI::sum(
          active_set_locally_owned.n_elements(), mpi_communicator);
      pcout << "         Size of active set: " << sum_contact_constraints
          << std::endl;
      unsigned int sum_contact_hanging_nodes = Utilities::MPI::sum(
          counter_hanging_nodes, mpi_communicator);
      pcout << "         Number of hanging nodes in contact: "
          << sum_contact_hanging_nodes << std::endl;

      solution = distributed_solution;

      constraints.close();

      //  	constraints_dirichlet_hanging_nodes.print (std::cout);

      constraints.merge(constraints_dirichlet_hanging_nodes);

      //constraints.print (std::cout);
    }

// @sect4{PlasticityContactProblem::dirichlet_constraints}

// This function defines the new ConstraintMatrix
// constraints_dirichlet_hanging_nodes. It contains
// the Dirichlet boundary values as well as the
// hanging nodes constraints.
  template <int dim>
    void
    PlasticityContactProblem<dim>::dirichlet_constraints ()
    {
      /* boundary_indicators:
       _______
       /  1    /|
       /______ / |
       8|       | 8|
       |   8   | /
       |_______|/
       6
       */

      constraints_dirichlet_hanging_nodes.reinit(locally_relevant_dofs);
      constraints_dirichlet_hanging_nodes.merge(constraints_hanging_nodes);

      // interpolate all components of the solution
      VectorTools::interpolate_boundary_values(dof_handler,
          base_mesh == "box" ? 6 : 0, EquationData::BoundaryValues<dim>(),
          constraints_dirichlet_hanging_nodes, ComponentMask());

      // interpolate x- and y-components of the
      // solution (this is a bit mask, so apply
      // operator| )
      FEValuesExtractors::Scalar x_displacement(0);
      FEValuesExtractors::Scalar y_displacement(1);
      VectorTools::interpolate_boundary_values(dof_handler, 8,
          EquationData::BoundaryValues<dim>(),
          constraints_dirichlet_hanging_nodes,
          (fe.component_mask(x_displacement) | fe.component_mask(y_displacement)));
      constraints_dirichlet_hanging_nodes.close();
    }

// @sect4{PlasticityContactProblem::solve}

// In addition to step-41 we have
// to deal with the hanging node
// constraints. Again we also consider
// the locally_owned_dofs only by
// creating the vector distributed_solution.
//
// For the hanging nodes we have to apply
// the set_zero function to system_rhs_newton.
// This is necessary if a hanging node value x_0
// has one neighbor which is in contact with
// value x_0 and one neighbor which is not with
// value x_1. This leads to an inhomogeneity
// constraint with value x_1/2 = gap/2 in the
// ConstraintMatrix.
// So the corresponding entries in the
// ride-hang-side are non-zero with a
// meaningless value. These values have to
// to set to zero.

// The rest of the function is similar to
// step-41 except that we use a FGMRES-solver
// instead of CG. For a very small hardening
// value gamma the linear system becomes
// almost semi definite but still symmetric.
  template <int dim>
    void
    PlasticityContactProblem<dim>::solve ()
    {
      TimerOutput::Scope t(computing_timer, "Solve");

      TrilinosWrappers::MPI::Vector distributed_solution(system_rhs_newton);
      distributed_solution = solution;

      constraints_hanging_nodes.set_zero(distributed_solution);
      constraints_hanging_nodes.set_zero(system_rhs_newton);
      distributed_solution.compress(VectorOperation::insert);
      system_rhs_newton.compress(VectorOperation::insert);

        {
          TimerOutput::Scope t(computing_timer, "Solve: setup preconditioner");
          preconditioner_u.initialize(system_matrix_newton, additional_data);
        }

        {
          TimerOutput::Scope t(computing_timer, "Solve: iterate");

          PrimitiveVectorMemory<TrilinosWrappers::MPI::Vector> mem;
          TrilinosWrappers::MPI::Vector tmp(system_rhs_newton);
          // 1e-4 seems to be the fasted option altogether, but to get more
          // reproducible parallel benchmark results, we use a small residual:
          double relative_accuracy = 1e-8;
          if (output_dir.compare("its/") == 0)
            relative_accuracy = 1e-4;

          const double solver_tolerance = relative_accuracy
              * system_matrix_newton.residual(tmp, distributed_solution,
                  system_rhs_newton);

          SolverControl solver_control(system_matrix_newton.m(),
              solver_tolerance);
          SolverBicgstab<TrilinosWrappers::MPI::Vector> solver(solver_control,
              mem/*,
               SolverFGMRES<TrilinosWrappers::MPI::Vector>::
               AdditionalData(30, true)*/);
          solver.solve(system_matrix_newton, distributed_solution,
              system_rhs_newton, preconditioner_u);

          pcout << "         Error: " << solver_control.initial_value()
              << " -> " << solver_control.last_value() << " in "
              << solver_control.last_step() << " Bicgstab iterations."
              << std::endl;

          number_iterations += solver_control.last_step();
        }

      constraints.distribute(distributed_solution);

      solution = distributed_solution;
    }

// @sect4{PlasticityContactProblem::solve_newton}

// In this function the damped Newton method is implemented.
// That means two nested loops: the outer loop for the newton
// iteration and the inner loop for the damping steps which
// will be used only if necessary. To obtain a good and reasonable
// starting value we solve an elastic problem in very first step (j=1).
  template <int dim>
    void
    PlasticityContactProblem<dim>::solve_newton ()
    {
      TimerOutput::Scope t(computing_timer, "solve newton setup");

      double resid = 0;
      double resid_old = 100000;
      TrilinosWrappers::MPI::Vector old_solution(system_rhs_newton);
      TrilinosWrappers::MPI::Vector res(system_rhs_newton);
      TrilinosWrappers::MPI::Vector tmp_vector(system_rhs_newton);

      std::vector < std::vector<bool> > constant_modes;
      DoFTools::extract_constant_modes(dof_handler, ComponentMask(),
          constant_modes);

      double sigma_hlp = sigma_0;

      additional_data.constant_modes = constant_modes;
      additional_data.elliptic = true;
      additional_data.n_cycles = 1;
      additional_data.w_cycle = false;
      additional_data.output_details = false;
      additional_data.smoother_sweeps = 2;
      additional_data.aggregation_threshold = 1e-2;

      IndexSet active_set_old(active_set);

      t.stop(); // stop newton setup timer

      unsigned int j = 1;
      unsigned int number_assemble_system = 0;
      for (; j <= 100; j++)
        {
          if (transfer_solution)
            {
              if (transfer_solution && j == 1 && cycle == 0)
                plast_lin_hard->set_sigma_0(1e+10);
              else if (transfer_solution && (j == 2 || cycle > 0))
                plast_lin_hard->set_sigma_0(sigma_hlp);
            }
          else
            {
              if (j == 1)
                plast_lin_hard->set_sigma_0(1e+10);
              else
                plast_lin_hard->set_sigma_0(sigma_hlp);
            }

          pcout << " " << std::endl;
          pcout << "   Newton iteration " << j << std::endl;
          pcout << "      Updating active set..." << std::endl;

            {
              TimerOutput::Scope t(computing_timer, "update active set");
              update_solution_and_constraints();
            }

          pcout << "      Assembling system... " << std::endl;
          system_matrix_newton = 0;
          system_rhs_newton = 0;
          assemble_nl_system(solution); //compute Newton-Matrix

          number_assemble_system += 1;

          pcout << "      Solving system... " << std::endl;
          solve();

          TrilinosWrappers::MPI::Vector distributed_solution(system_rhs_newton);
          distributed_solution = solution;

          // We handle a highly nonlinear problem so we have to damp
          // the Newtons method. We refer that we iterate the new solution
          // in each Newton step and not only the solution update.
          // Since the solution set is a convex set and not a space we
          // compute for the damping a linear combination of the
          // previous and the current solution to guarantee that the
          // damped solution is in our solution set again.
          // At most we apply 10 damping steps.
          bool damped = false;
          tmp_vector = old_solution;
          double a = 0;
          for (unsigned int i = 0; (i < 5) && (!damped); i++)
            {
              a = std::pow(0.5, static_cast<double>(i));
              old_solution = tmp_vector;
              old_solution.sadd(1 - a, a, distributed_solution);
              old_solution.compress(VectorOperation::add);

              TimerOutput::Scope t(computing_timer, "Residual and lambda");

              system_rhs_newton = 0;
              system_rhs_lambda = 0;

              solution = old_solution;
              residual_nl_system(solution);
              res = system_rhs_newton;

              const unsigned int start_res = (res.local_range().first),
                  end_res = (res.local_range().second);
              for (unsigned int n = start_res; n < end_res; ++n)
                if (constraints.is_inhomogeneously_constrained(n))
                  res(n) = 0;

              res.compress(VectorOperation::insert);

              resid = res.l2_norm();

              if (resid < resid_old)
                damped = true;

              pcout << "      Residual of the non-contact part of the system: "
                  << resid << std::endl
                  << "         with a damping parameter alpha = " << a
                  << std::endl;

              // The previous iteration of step 0 is the solution of an elastic problem.
              // So a linear combination of a plastic and an elastic solution makes no sense
              // since the elastic solution is not in the convex set of the plastic solution.
              if (!transfer_solution && j == 2)
                break;
              if (transfer_solution && j == 2 && cycle == 0)
                break;
            }

          resid_old = resid;

          resid_vector = system_rhs_lambda;
          resid_vector.compress(VectorOperation::insert);

          int is_my_set_changed = (active_set == active_set_old) ? 0 : 1;
          int num_changed = Utilities::MPI::sum(is_my_set_changed,
              MPI_COMM_WORLD);
          if (num_changed == 0)
            {
              pcout << "      Active set did not change!" << std::endl;
              if (output_dir.compare("its/") != 0 && resid < 1e-7)
                break;
              else if (output_dir.compare("its/") == 0 && resid < 1e-10)
                break;
            }
          active_set_old = active_set;
        }

      pcout << "" << std::endl << "      Number of assembled systems = "
          << number_assemble_system << std::endl
          << "      Number of Solver-Iterations = " << number_iterations
          << std::endl;
    }

// @sect3{The <code>refine_grid</code> function}

  template <int dim>
    void
    PlasticityContactProblem<dim>::refine_grid ()
    {
      if (refinement_strategy == RefinementStrategy::refine_global)
        {
          triangulation.refine_global(1);
        }
      else
        {
          Vector<float> estimated_error_per_cell(
              triangulation.n_active_cells());
          KellyErrorEstimator<dim>::estimate(dof_handler,
              QGauss<dim - 1>(fe.degree + 2), typename FunctionMap<dim>::type(),
              solution, estimated_error_per_cell);

          parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
              triangulation, estimated_error_per_cell, 0.3, 0.03);

          triangulation.prepare_coarsening_and_refinement();
          if (transfer_solution)
            soltrans->prepare_for_coarsening_and_refinement(solution);

          triangulation.execute_coarsening_and_refinement();
        }
    }

// @sect3{The <code>move_mesh</code> function}

  template <int dim>
    void
    PlasticityContactProblem<dim>::move_mesh (
        const TrilinosWrappers::MPI::Vector &_complete_displacement) const
    {
      std::vector<bool> vertex_touched(triangulation.n_vertices(), false);

      for (typename DoFHandler<dim>::active_cell_iterator cell =
          dof_handler.begin_active(); cell != dof_handler.end(); ++cell)
        if (cell->is_locally_owned())
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
              ++v)
            {
              if (vertex_touched[cell->vertex_index(v)] == false)
                {
                  vertex_touched[cell->vertex_index(v)] = true;

                  Point<dim> vertex_displacement;
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      if (_complete_displacement(cell->vertex_dof_index(v, d))
                          != 0)
                        vertex_displacement[d] = _complete_displacement(
                            cell->vertex_dof_index(v, d));
                    }

                  cell->vertex(v) += vertex_displacement;
                }
            }
    }

// @sect4{PlasticityContactProblem::output_results}

  template <int dim>
    void
    PlasticityContactProblem<dim>::output_results (
        const std::string &title)
    {
      move_mesh(solution);

      // Calculation of the contact forces
      TrilinosWrappers::MPI::Vector lambda(solution);
      TrilinosWrappers::MPI::Vector distributed_lambda(system_rhs_newton);
      const unsigned int start_res = (resid_vector.local_range().first),
          end_res = (resid_vector.local_range().second);
      for (unsigned int n = start_res; n < end_res; ++n)
        if (constraints.is_inhomogeneously_constrained(n))
          distributed_lambda(n) = resid_vector(n) / diag_mass_matrix_vector(n);
      distributed_lambda.compress(VectorOperation::insert);
      constraints_hanging_nodes.distribute(distributed_lambda);
      lambda = distributed_lambda;
      TrilinosWrappers::MPI::Vector resid_vector_relevant(solution);
      TrilinosWrappers::MPI::Vector distributed_resid_vector(resid_vector);
      constraints_hanging_nodes.distribute(distributed_resid_vector);
      resid_vector_relevant = distributed_resid_vector;

      DataOut<dim> data_out;

      data_out.attach_dof_handler(dof_handler);

      const std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation(
          dim, DataComponentInterpretation::component_is_part_of_vector);
      data_out.add_data_vector(solution,
          std::vector < std::string > (dim, "Displacement"),
          DataOut<dim>::type_dof_data, data_component_interpretation);
      data_out.add_data_vector(lambda,
          std::vector < std::string > (dim, "ContactForce"),
          DataOut<dim>::type_dof_data, data_component_interpretation);
      data_out.add_data_vector(active_set,
          std::vector < std::string > (dim, "ActiveSet"),
          DataOut<dim>::type_dof_data, data_component_interpretation);
      data_out.add_data_vector(resid_vector_relevant,
          std::vector < std::string > (dim, "Residual"),
          DataOut<dim>::type_dof_data, data_component_interpretation);

      Vector<float> subdomain(triangulation.n_active_cells());
      for (unsigned int i = 0; i < subdomain.size(); ++i)
        subdomain(i) = triangulation.locally_owned_subdomain();
      data_out.add_data_vector(subdomain, "subdomain");

      data_out.add_data_vector(cell_constitution, "CellConstitution");

      data_out.build_patches();

      const std::string filename =
          (output_dir + title + "-"
              + Utilities::int_to_string(
                  triangulation.locally_owned_subdomain(), 4));

      std::ofstream output_vtu((filename + ".vtu").c_str());
      data_out.write_vtu(output_vtu);
      pcout << output_dir + title << ".pvtu" << std::endl;

      if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
        {
          std::vector < std::string > filenames;
          for (unsigned int i = 0;
              i < Utilities::MPI::n_mpi_processes(mpi_communicator); ++i)
            filenames.push_back(
                title + "-" + Utilities::int_to_string(i, 4) + ".vtu");

          std::ofstream master_output((output_dir + title + ".pvtu").c_str());
          data_out.write_pvtu_record(master_output, filenames);
        }

      TrilinosWrappers::MPI::Vector tmp(solution);
      tmp *= -1;
      move_mesh(tmp);
    }

// @sect4{PlasticityContactProblem::output_contact_force}

// This function provides the contact force by calculating
// an integral over the contact pressure in z-directions
// over the contact area. For this purpose we set the contact
// pressure lambda to 0 for all inactive dofs. For all
// active dofs we lambda contains the quotient of the nonlinear
// residual (resid_vector) and corresponding diagonal entry
// of the mass matrix (diag_mass_matrix_vector). Because it is
// not unlikely that hanging nodes shows up in the contact area
// it is important to apply contraints_hanging_nodes.distribute
// to the distributed_lambda vector.
// To calculate the contact pressure in a certain point in the
// contact area, we apply the Functions::FEFieldFunction.
// In parallel this is a little tricky because we have to find the
// process with the right cell which contains this point. If
// a processor does not own the cell with the point we have to
// catch these cases.
  template <int dim>
    void
    PlasticityContactProblem<dim>::output_contact_force (
        const unsigned int cycle)
    {
      Functions::FEFieldFunction<dim, DoFHandler<dim>,
          TrilinosWrappers::MPI::Vector> solution_function(dof_handler,
          solution);
      std::cout.precision(10);

      Vector<double> solution_p1(dim);
      std::vector<Tensor<1, dim> > solution_gradient_p1(dim);

      // Here we calculate the contact pressure as a vector lambda.
      // If a dof is element of the active set lambda contains the
      // nonlinear residual this dof divided by the according entry
      // of the mass matrix. In all other dofs lambda will be set to
      // zero.
      TrilinosWrappers::MPI::Vector lambda(solution);
      TrilinosWrappers::MPI::Vector distributed_lambda(system_rhs_newton);
      const unsigned int start_res = (resid_vector.local_range().first),
          end_res = (resid_vector.local_range().second);
      for (unsigned int n = start_res; n < end_res; ++n)
        if (constraints.is_inhomogeneously_constrained(n))
          distributed_lambda(n) = resid_vector(n) / diag_mass_matrix_vector(n);
        else
          distributed_lambda(n) = 0;
      distributed_lambda.compress(VectorOperation::insert);
      constraints_hanging_nodes.distribute(distributed_lambda);
      lambda = distributed_lambda;
      Functions::FEFieldFunction<dim, DoFHandler<dim>,
          TrilinosWrappers::MPI::Vector> lambda_function(dof_handler, lambda);

      // Here we try to find the MPI-process which owns the cell
      // with the point_of_interest. If it is the wrong MPI-process
      // we catch this case and set point_found to false.
      const Point<dim> point_of_interest(0.49, 0.5001, 1.0);
      Vector<double> contact_pressure_in_point(dim);
      bool point_found = true;

      MPI_Barrier(MPI_COMM_WORLD);
      try
        {
          lambda_function.vector_value(point_of_interest,
              contact_pressure_in_point);
        }
      catch (const typename Functions::FEFieldFunction<dim, DoFHandler<dim>,
          TrilinosWrappers::MPI::Vector>::ExcPointNotAvailableHere &)
        {
          point_found = false;
        }

      if (point_found == true)
        {
          std::cout << "PoI contact pressure: " << contact_pressure_in_point(2)
              << std::endl;
        }

      // To obtain the contact force we have to compute an integral of the contact pressure
      // in z-direction over the whole contact area. To be accurate enough we use the
      // Gaussian quadrature rule with fe.degree + 1.
      double contact_force = 0.0;
        {
          QGauss<dim - 1> face_quadrature_formula(fe.degree + 1);

          FEFaceValues<dim> fe_values_face(fe, face_quadrature_formula,
              update_values | update_quadrature_points | update_JxW_values);

          const unsigned int n_face_q_points = face_quadrature_formula.size();

          const FEValuesExtractors::Vector displacement(0);

          typename DoFHandler<dim>::active_cell_iterator cell =
              dof_handler.begin_active(), endc = dof_handler.end();
          for (; cell != endc; ++cell)
            if (cell->is_locally_owned())
              for (unsigned int face = 0;
                  face < GeometryInfo<dim>::faces_per_cell; ++face)
                if (cell->face(face)->at_boundary()
                    && cell->face(face)->boundary_indicator() == 1)
                  {
                    fe_values_face.reinit(cell, face);

                    std::vector<Tensor<1, dim> > lambda_values(n_face_q_points);
                    fe_values_face[displacement].get_function_values(lambda,
                        lambda_values);

                    for (unsigned int q_point = 0; q_point < n_face_q_points;
                        ++q_point)
                      {
                        contact_force += lambda_values[q_point][2]
                            * fe_values_face.JxW(q_point);
                      }
                  }
          contact_force = Utilities::MPI::sum(contact_force, MPI_COMM_WORLD);
          pcout << "Contact force = " << contact_force << std::endl;
        }
      MPI_Barrier(MPI_COMM_WORLD);
    }

// @sect4{PlasticityContactProblem::run}

  template <int dim>
    void
    PlasticityContactProblem<dim>::run ()
    {

      if (obstacle_filename != "")
        {
          pcout << "Read the obstacle from '" << obstacle_filename << "' ... "
              << std::flush;
          input_obstacle.reset(new Input<dim>(obstacle_filename.c_str()));
          pcout << "done." << std::endl;
        }

      computing_timer.reset();
      for (cycle = 0; cycle < n_cycles; ++cycle)
        {
            {
              TimerOutput::Scope t(computing_timer, "Setup");

              pcout << std::endl;
              pcout << "Cycle " << cycle << ':' << std::endl;

              if (cycle == 0)
                {
                  make_grid();
                }
              else
                {
                  TimerOutput::Scope t(computing_timer, "Setup: refine mesh");
                  if (transfer_solution)
                    soltrans.reset(
                        new parallel::distributed::SolutionTransfer<dim,
                            TrilinosWrappers::MPI::Vector>(dof_handler));
                  refine_grid();
                }

              setup_system();

              if (transfer_solution && cycle > 0)
                {
                  TrilinosWrappers::MPI::Vector distributed_solution(
                      system_rhs_newton);
                  distributed_solution = solution;
                  soltrans->interpolate(distributed_solution);
                  solution = distributed_solution;
                  residual_nl_system(solution);
                  resid_vector = system_rhs_lambda;
                  resid_vector.compress(VectorOperation::insert);
                }

            }

          solve_newton();

          if (true) //Utilities::MPI::n_mpi_processes(mpi_communicator) <= 64)
            {
              pcout << "      Writing graphical output... " << std::flush;

              TimerOutput::Scope t(computing_timer, "Graphical output");

              std::ostringstream filename_solution;
              filename_solution << "solution-";
              filename_solution << Utilities::int_to_string(cycle, 2);
              output_results(filename_solution.str());
            }

          computing_timer.print_summary();
          computing_timer.reset();

          Utilities::System::MemoryStats stats;
          Utilities::System::get_memory_stats(stats);
          pcout << "VMPEAK, Resident in kB: " << stats.VmSize << " "
              << stats.VmRSS << std::endl;

          if (base_mesh == "box")
            output_contact_force(cycle);
        }
    }
}

// @sect3{The <code>main</code> function}

int
main (
    int argc, char *argv[])
{
  using namespace dealii;
  using namespace Step42;

  deallog.depth_console(0);
  ParameterHandler prm;
  PlasticityContactProblem<3>::declare(prm);
  if (argc != 2)
    {
      prm.print_parameters(std::cout, ParameterHandler::Text);
      return 0;
    }

  prm.read_input(argv[1]);
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
    {
      PlasticityContactProblem<3> problem(prm);
      problem.run();
    }

  return 0;
}
