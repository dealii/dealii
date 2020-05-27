/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2020 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

/*
 * Authors: TODO, 2020
 */


// We start by including all the necessary deal.II header files and some C++
// related ones. They have been discussed in detail in previous tutorial
// programs, so you need only refer to past tutorials for details.
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/quadrature_point_data.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>

#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_eulerian.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <deal.II/base/config.h>
// #if DEAL_II_VERSION_MAJOR >= 9 && (defined(DEAL_II_WITH_ADOLC) || defined(DEAL_II_WITH_TRILINOS))
#include <deal.II/differentiation/ad.h>
#define ENABLE_AD_FORMULATION
// #endif

// These must be included below the AD headers so that
// their math functions are available for use in the
// definition of tensors and kinematic quantities
#include <deal.II/physics/elasticity/kinematics.h>
#include <deal.II/physics/elasticity/standard_tensors.h>

#include <iostream>
#include <fstream>


// We then stick everything that relates to this tutorial program into a
// namespace of its own, and import all the deal.II function and class names
// into it:
namespace Cook_Membrane
{
  using namespace dealii;

// @sect3{Run-time parameters}
//
// There are several parameters that can be set in the code so we set up a
// ParameterHandler object to read in the choices at run-time.
  namespace Parameters
  {
// @sect4{Assembly method}

// Here we specify whether automatic differentiation is to be used to assemble
// the linear system, and if so then what order of differentiation is to be
// employed.
    struct AssemblyMethod
    {
      unsigned int automatic_differentiation_order;
      std::string automatic_differentiation_library;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };


    void AssemblyMethod::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Assembly method");
      {
        prm.declare_entry("Automatic differentiation order", "0",
                          Patterns::Integer(0,2),
                          "The automatic differentiation order to be used in the assembly of the linear system.\n"
                          "# Order = 0: Both the residual and linearisation are computed manually.\n"
                          "# Order = 1: The residual is computed manually but the linearisation is performed using AD.\n"
                          "# Order = 2: Both the residual and linearisation are computed using AD.");

        prm.declare_entry("Automatic differentiation library", "None",
                          Patterns::Selection("None|ADOL-C|Sacado"),
                          "The automatic differentiation library to be used in the assembly of the linear system.");
      }
      prm.leave_subsection();
    }

    void AssemblyMethod::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Assembly method");
      {
        automatic_differentiation_order = prm.get_integer("Automatic differentiation order");
        automatic_differentiation_library = prm.get("Automatic differentiation library");
      }
      prm.leave_subsection();
    }

// @sect4{Finite Element system}

// Here we specify the polynomial order used to approximate the solution.
// The quadrature order should be adjusted accordingly.
    struct FESystem
    {
      unsigned int poly_degree;
      unsigned int quad_order;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };


    void FESystem::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Finite element system");
      {
        prm.declare_entry("Polynomial degree", "2",
                          Patterns::Integer(0),
                          "Displacement system polynomial order");

        prm.declare_entry("Quadrature order", "3",
                          Patterns::Integer(0),
                          "Gauss quadrature order");
      }
      prm.leave_subsection();
    }

    void FESystem::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Finite element system");
      {
        poly_degree = prm.get_integer("Polynomial degree");
        quad_order = prm.get_integer("Quadrature order");
      }
      prm.leave_subsection();
    }

// @sect4{Geometry}

// Make adjustments to the problem geometry and its discretisation.
    struct Geometry
    {
      unsigned int elements_per_edge;
      double       scale;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    void Geometry::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry");
      {
        prm.declare_entry("Elements per edge", "32",
                          Patterns::Integer(0),
                          "Number of elements per long edge of the beam");

        prm.declare_entry("Grid scale", "1e-3",
                          Patterns::Double(0.0),
                          "Global grid scaling factor");
      }
      prm.leave_subsection();
    }

    void Geometry::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry");
      {
        elements_per_edge = prm.get_integer("Elements per edge");
        scale = prm.get_double("Grid scale");
      }
      prm.leave_subsection();
    }

// @sect4{Materials}

// We also need the shear modulus $ \mu $ and Poisson ration $ \nu $ for the
// neo-Hookean material.
    struct Materials
    {
      double nu;
      double mu;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    void Materials::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material properties");
      {
        prm.declare_entry("Poisson's ratio", "0.3",
                          Patterns::Double(-1.0,0.5),
                          "Poisson's ratio");

        prm.declare_entry("Shear modulus", "0.4225e6",
                          Patterns::Double(),
                          "Shear modulus");
      }
      prm.leave_subsection();
    }

    void Materials::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Material properties");
      {
        nu = prm.get_double("Poisson's ratio");
        mu = prm.get_double("Shear modulus");
      }
      prm.leave_subsection();
    }

// @sect4{Linear solver}

// Next, we choose both solver and preconditioner settings.  The use of an
// effective preconditioner is critical to ensure convergence when a large
// nonlinear motion occurs within a Newton increment.
    struct LinearSolver
    {
      std::string type_lin;
      double      tol_lin;
      std::string preconditioner_type;
      double      preconditioner_relaxation;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    void LinearSolver::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Linear solver");
      {
        prm.declare_entry("Solver type", "CG",
                          Patterns::Selection("CG|Direct"),
                          "Type of solver used to solve the linear system");

        prm.declare_entry("Residual", "1e-6",
                          Patterns::Double(0.0),
                          "Linear solver residual (scaled by residual norm)");

        prm.declare_entry("Preconditioner type", "ssor",
                          Patterns::Selection("jacobi|ssor"),
                          "Type of preconditioner");

        prm.declare_entry("Preconditioner relaxation", "0.65",
                          Patterns::Double(0.0),
                          "Preconditioner relaxation value");
      }
      prm.leave_subsection();
    }

    void LinearSolver::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Linear solver");
      {
        type_lin = prm.get("Solver type");
        tol_lin = prm.get_double("Residual");
        preconditioner_type = prm.get("Preconditioner type");
        preconditioner_relaxation = prm.get_double("Preconditioner relaxation");
      }
      prm.leave_subsection();
    }

// @sect4{Nonlinear solver}

// A Newton-Raphson scheme is used to solve the nonlinear system of governing
// equations.  We now define the tolerances and the maximum number of
// iterations for the Newton-Raphson nonlinear solver.
    struct NonlinearSolver
    {
      unsigned int max_iterations_NR;
      double       tol_f;
      double       tol_u;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    void NonlinearSolver::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Nonlinear solver");
      {
        prm.declare_entry("Max iterations Newton-Raphson", "10",
                          Patterns::Integer(0),
                          "Number of Newton-Raphson iterations allowed");

        prm.declare_entry("Tolerance force", "1.0e-9",
                          Patterns::Double(0.0),
                          "Force residual tolerance");

        prm.declare_entry("Tolerance displacement", "1.0e-6",
                          Patterns::Double(0.0),
                          "Displacement error tolerance");
      }
      prm.leave_subsection();
    }

    void NonlinearSolver::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Nonlinear solver");
      {
        max_iterations_NR = prm.get_integer("Max iterations Newton-Raphson");
        tol_f = prm.get_double("Tolerance force");
        tol_u = prm.get_double("Tolerance displacement");
      }
      prm.leave_subsection();
    }

// @sect4{Time}

// Set the timestep size $ \varDelta t $ and the simulation end-time.
    struct Time
    {
      double delta_t;
      double end_time;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    void Time::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Time");
      {
        prm.declare_entry("End time", "1",
                          Patterns::Double(),
                          "End time");

        prm.declare_entry("Time step size", "0.1",
                          Patterns::Double(),
                          "Time step size");
      }
      prm.leave_subsection();
    }

    void Time::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Time");
      {
        end_time = prm.get_double("End time");
        delta_t = prm.get_double("Time step size");
      }
      prm.leave_subsection();
    }

// @sect4{All parameters}

// Finally we consolidate all of the above structures into a single container
// that holds all of our run-time selections.
    struct AllParameters :
      public AssemblyMethod,
      public FESystem,
      public Geometry,
      public Materials,
      public LinearSolver,
      public NonlinearSolver,
      public Time

    {
      AllParameters(const std::string &input_file);

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    AllParameters::AllParameters(const std::string &input_file)
    {
      ParameterHandler prm;
      declare_parameters(prm);
      prm.parse_input(input_file);
      parse_parameters(prm);
    }

    void AllParameters::declare_parameters(ParameterHandler &prm)
    {
      AssemblyMethod::declare_parameters(prm);
      FESystem::declare_parameters(prm);
      Geometry::declare_parameters(prm);
      Materials::declare_parameters(prm);
      LinearSolver::declare_parameters(prm);
      NonlinearSolver::declare_parameters(prm);
      Time::declare_parameters(prm);
    }

    void AllParameters::parse_parameters(ParameterHandler &prm)
    {
      AssemblyMethod::parse_parameters(prm);
      FESystem::parse_parameters(prm);
      Geometry::parse_parameters(prm);
      Materials::parse_parameters(prm);
      LinearSolver::parse_parameters(prm);
      NonlinearSolver::parse_parameters(prm);
      Time::parse_parameters(prm);
    }
  }


// @sect3{Time class}

// A simple class to store time data. Its functioning is transparent so no
// discussion is necessary. For simplicity we assume a constant time step
// size.
  class Time
  {
  public:
    Time (const double time_end,
          const double delta_t)
      :
      timestep(0),
      time_current(0.0),
      time_end(time_end),
      delta_t(delta_t)
    {}

    virtual ~Time()
    {}

    double current() const
    {
      return time_current;
    }
    double end() const
    {
      return time_end;
    }
    double get_delta_t() const
    {
      return delta_t;
    }
    unsigned int get_timestep() const
    {
      return timestep;
    }
    void increment()
    {
      time_current += delta_t;
      ++timestep;
    }

  private:
    unsigned int timestep;
    double       time_current;
    const double time_end;
    const double delta_t;
  };

// @sect3{Compressible neo-Hookean material within a one-field formulation}

// As discussed in the literature and step-44, Neo-Hookean materials are a type
// of hyperelastic materials.  The entire domain is assumed to be composed of a
// compressible neo-Hookean material.  This class defines the behaviour of
// this material within a one-field formulation.  Compressible neo-Hookean
// materials can be described by a strain-energy function (SEF) $ \Psi =
// \Psi_{\text{iso}}(\overline{\mathbf{b}}) + \Psi_{\text{vol}}(J)
// $.
//
// The isochoric response is given by $
// \Psi_{\text{iso}}(\overline{\mathbf{b}}) = c_{1} [\overline{I}_{1} - 3] $
// where $ c_{1} = \frac{\mu}{2} $ and $\overline{I}_{1}$ is the first
// invariant of the left- or right-isochoric Cauchy-Green deformation tensors.
// That is $\overline{I}_1 :=\textrm{tr}(\overline{\mathbf{b}})$.  In this
// example the SEF that governs the volumetric response is defined as $
// \Psi_{\text{vol}}(J) = \kappa \frac{1}{4} [ J^2 - 1
// - 2\textrm{ln}\; J ]$,  where $\kappa:= \lambda + 2/3 \mu$ is
// the <a href="http://en.wikipedia.org/wiki/Bulk_modulus">bulk modulus</a>
// and $\lambda$ is <a
// href="http://en.wikipedia.org/wiki/Lam%C3%A9_parameters">Lame's first
// parameter</a>.
//
// The following class will be used to characterize the material we work with,
// and provides a central point that one would need to modify if one were to
// implement a different material model. For it to work, we will store one
// object of this type per quadrature point, and in each of these objects
// store the current state (characterized by the values or measures  of the
// displacement field) so that we can compute the elastic coefficients
// linearized around the current state.
  template <int dim,typename NumberType>
  class Material_Compressible_Neo_Hook_One_Field
  {
  public:
    Material_Compressible_Neo_Hook_One_Field(const double mu,
                                             const double nu)
      :
      kappa((2.0 * mu * (1.0 + nu)) / (3.0 * (1.0 - 2.0 * nu))),
      c_1(mu / 2.0)
    {
      Assert(kappa > 0, ExcInternalError());
    }

    ~Material_Compressible_Neo_Hook_One_Field()
    {}

    // The first function is the total energy
    // $\Psi = \Psi_{\textrm{iso}} + \Psi_{\textrm{vol}}$.
    NumberType
    get_Psi(const NumberType                        &det_F,
            const SymmetricTensor<2,dim,NumberType> &b_bar) const
    {
      return get_Psi_vol(det_F) + get_Psi_iso(b_bar);
    }

    // The second function determines the Kirchhoff stress $\boldsymbol{\tau}
    // = \boldsymbol{\tau}_{\textrm{iso}} + \boldsymbol{\tau}_{\textrm{vol}}$
    SymmetricTensor<2,dim,NumberType>
    get_tau(const NumberType                        &det_F,
            const SymmetricTensor<2,dim,NumberType> &b_bar)
    {
      // See Holzapfel p231 eq6.98 onwards
      return get_tau_vol(det_F) + get_tau_iso(b_bar);
    }

    // The fourth-order elasticity tensor in the spatial setting
    // $\mathfrak{c}$ is calculated from the SEF $\Psi$ as $ J
    // \mathfrak{c}_{ijkl} = F_{iA} F_{jB} \mathfrak{C}_{ABCD} F_{kC} F_{lD}$
    // where $ \mathfrak{C} = 4 \frac{\partial^2 \Psi(\mathbf{C})}{\partial
    // \mathbf{C} \partial \mathbf{C}}$
    SymmetricTensor<4,dim,NumberType>
    get_Jc(const NumberType                        &det_F,
           const SymmetricTensor<2,dim,NumberType> &b_bar) const
    {
      return get_Jc_vol(det_F) + get_Jc_iso(b_bar);
    }

  private:
    // Define constitutive model parameters $\kappa$ (bulk modulus) and the
    // neo-Hookean model parameter $c_1$:
    const double kappa;
    const double c_1;

    // Value of the volumetric free energy
    NumberType
    get_Psi_vol(const NumberType &det_F) const
    {
      return (kappa / 4.0) * (det_F*det_F - 1.0 - 2.0*std::log(det_F));
    }

    // Value of the isochoric free energy
    NumberType
    get_Psi_iso(const SymmetricTensor<2,dim,NumberType> &b_bar) const
    {
      return c_1 * (trace(b_bar) - dim);
    }

    // Derivative of the volumetric free energy with respect to
    // $J$ return $\frac{\partial
    // \Psi_{\text{vol}}(J)}{\partial J}$
    NumberType
    get_dPsi_vol_dJ(const NumberType &det_F) const
    {
      return (kappa / 2.0) * (det_F - 1.0 / det_F);
    }

    // The following functions are used internally in determining the result
    // of some of the public functions above. The first one determines the
    // volumetric Kirchhoff stress $\boldsymbol{\tau}_{\textrm{vol}}$.
    // Note the difference in its definition when compared to step-44.
    SymmetricTensor<2,dim,NumberType>
    get_tau_vol(const NumberType &det_F) const
    {
      return NumberType(get_dPsi_vol_dJ(det_F) * det_F) * Physics::Elasticity::StandardTensors<dim>::I;
    }

    // Next, determine the isochoric Kirchhoff stress
    // $\boldsymbol{\tau}_{\textrm{iso}} =
    // \mathcal{P}:\overline{\boldsymbol{\tau}}$:
    SymmetricTensor<2,dim,NumberType>
    get_tau_iso(const SymmetricTensor<2,dim,NumberType> &b_bar) const
    {
      return Physics::Elasticity::StandardTensors<dim>::dev_P * get_tau_bar(b_bar);
    }

    // Then, determine the fictitious Kirchhoff stress
    // $\overline{\boldsymbol{\tau}}$:
    SymmetricTensor<2,dim,NumberType>
    get_tau_bar(const SymmetricTensor<2,dim,NumberType> &b_bar) const
    {
      return 2.0 * c_1 * b_bar;
    }

    // Second derivative of the volumetric free energy wrt $J$. We
    // need the following computation explicitly in the tangent so we make it
    // public.  We calculate $\frac{\partial^2
    // \Psi_{\textrm{vol}}(J)}{\partial J \partial
    // J}$
    NumberType
    get_d2Psi_vol_dJ2(const NumberType &det_F) const
    {
      return ( (kappa / 2.0) * (1.0 + 1.0 / (det_F * det_F)));
    }

    // Calculate the volumetric part of the tangent $J
    // \mathfrak{c}_\textrm{vol}$. Again, note the difference in its
    // definition when compared to step-44. The extra terms result from two
    // quantities in $\boldsymbol{\tau}_{\textrm{vol}}$ being dependent on
    // $\boldsymbol{F}$.
    SymmetricTensor<4,dim,NumberType>
    get_Jc_vol(const NumberType &det_F) const
    {
      // See Holzapfel p265
      return det_F
             * ( (get_dPsi_vol_dJ(det_F) + det_F * get_d2Psi_vol_dJ2(det_F))*Physics::Elasticity::StandardTensors<dim>::IxI
                 - (2.0 * get_dPsi_vol_dJ(det_F))*Physics::Elasticity::StandardTensors<dim>::S );
    }

    // Calculate the isochoric part of the tangent $J
    // \mathfrak{c}_\textrm{iso}$:
    SymmetricTensor<4,dim,NumberType>
    get_Jc_iso(const SymmetricTensor<2,dim,NumberType> &b_bar) const
    {
      const SymmetricTensor<2, dim> tau_bar = get_tau_bar(b_bar);
      const SymmetricTensor<2, dim> tau_iso = get_tau_iso(b_bar);
      const SymmetricTensor<4, dim> tau_iso_x_I
        = outer_product(tau_iso,
                        Physics::Elasticity::StandardTensors<dim>::I);
      const SymmetricTensor<4, dim> I_x_tau_iso
        = outer_product(Physics::Elasticity::StandardTensors<dim>::I,
                        tau_iso);
      const SymmetricTensor<4, dim> c_bar = get_c_bar();

      return (2.0 / dim) * trace(tau_bar)
             * Physics::Elasticity::StandardTensors<dim>::dev_P
             - (2.0 / dim) * (tau_iso_x_I + I_x_tau_iso)
             + Physics::Elasticity::StandardTensors<dim>::dev_P * c_bar
             * Physics::Elasticity::StandardTensors<dim>::dev_P;
    }

    // Calculate the fictitious elasticity tensor $\overline{\mathfrak{c}}$.
    // For the material model chosen this is simply zero:
    SymmetricTensor<4,dim,double>
    get_c_bar() const
    {
      return SymmetricTensor<4, dim>();
    }
  };

// @sect3{Quadrature point history}

// As seen in step-18, the <code> PointHistory </code> class offers a method
// for storing data at the quadrature points.  Here each quadrature point
// holds a pointer to a material description.  Thus, different material models
// can be used in different regions of the domain.  Among other data, we
// choose to store the Kirchhoff stress $\boldsymbol{\tau}$ and the tangent
// $J\mathfrak{c}$ for the quadrature points.
  template <int dim,typename NumberType>
  class PointHistory
  {
  public:
    PointHistory()
    {}

    virtual ~PointHistory()
    {}

    // The first function is used to create a material object and to
    // initialize all tensors correctly: The second one updates the stored
    // values and stresses based on the current deformation measure
    // $\textrm{Grad}\mathbf{u}_{\textrm{n}}$.
    void setup_lqp (const Parameters::AllParameters &parameters)
    {
      material.reset(new Material_Compressible_Neo_Hook_One_Field<dim,NumberType>(parameters.mu,
                     parameters.nu));
    }

    // We offer an interface to retrieve certain data.
    // This is the strain energy:
    NumberType
    get_Psi(const NumberType                        &det_F,
            const SymmetricTensor<2,dim,NumberType> &b_bar) const
    {
      return material->get_Psi(det_F,b_bar);
    }

    // Here are the kinetic variables. These are used in the material and
    // global tangent matrix and residual assembly operations:
    // First is the Kirchhoff stress:
    SymmetricTensor<2,dim,NumberType>
    get_tau(const NumberType                        &det_F,
            const SymmetricTensor<2,dim,NumberType> &b_bar) const
    {
      return material->get_tau(det_F,b_bar);
    }

    // And the tangent:
    SymmetricTensor<4,dim,NumberType>
    get_Jc(const NumberType                        &det_F,
           const SymmetricTensor<2,dim,NumberType> &b_bar) const
    {
      return material->get_Jc(det_F,b_bar);
    }

    // In terms of member functions, this class stores for the quadrature
    // point it represents a copy of a material type in case different
    // materials are used in different regions of the domain, as well as the
    // inverse of the deformation gradient...
  private:
    std::shared_ptr< Material_Compressible_Neo_Hook_One_Field<dim,NumberType> > material;
  };


// @sect3{Quasi-static compressible finite-strain solid}

  // Forward declarations for classes that will
  // perform assembly of the linear system.
  template <int dim,typename NumberType>
  struct Assembler_Base;
  template <int dim,typename NumberType>
  struct Assembler;

// The Solid class is the central class in that it represents the problem at
// hand. It follows the usual scheme in that all it really has is a
// constructor, destructor and a <code>run()</code> function that dispatches
// all the work to private functions of this class:
  template <int dim,typename NumberType>
  class Solid
  {
  public:
    Solid(const Parameters::AllParameters &parameters);

    virtual
    ~Solid();

    void
    run();

  private:

    // We start the collection of member functions with one that builds the
    // grid:
    void
    make_grid();

    // Set up the finite element system to be solved:
    void
    system_setup();

    // Several functions to assemble the system and right hand side matrices
    // using multithreading. Each of them comes as a wrapper function, one
    // that is executed to do the work in the WorkStream model on one cell,
    // and one that copies the work done on this one cell into the global
    // object that represents it:
    void
    assemble_system(const BlockVector<double> &solution_delta);

    // We use a separate data structure to perform the assembly. It needs access
    // to some low-level data, so we simply befriend the class instead of
    // creating a complex interface to provide access as necessary.
    friend struct Assembler_Base<dim,NumberType>;
    friend struct Assembler<dim,NumberType>;

    // Apply Dirichlet boundary conditions on the displacement field
    void
    make_constraints(const int &it_nr);

    // Create and update the quadrature points. Here, no data needs to be
    // copied into a global object, so the copy_local_to_global function is
    // empty:
    void
    setup_qph();

    // Solve for the displacement using a Newton-Raphson method. We break this
    // function into the nonlinear loop and the function that solves the
    // linearized Newton-Raphson step:
    void
    solve_nonlinear_timestep(BlockVector<double> &solution_delta);

    std::pair<unsigned int, double>
    solve_linear_system(BlockVector<double> &newton_update);

    // Solution retrieval as well as post-processing and writing data to file:
    BlockVector<double>
    get_total_solution(const BlockVector<double> &solution_delta) const;

    void
    output_results() const;

    // Finally, some member variables that describe the current state: A
    // collection of the parameters used to describe the problem setup...
    const Parameters::AllParameters &parameters;

    // ...the volume of the reference and current configurations...
    double                           vol_reference;
    double                           vol_current;

    // ...and description of the geometry on which the problem is solved:
    Triangulation<dim>               triangulation;

    // Also, keep track of the current time and the time spent evaluating
    // certain functions
    Time                             time;
    TimerOutput                      timer;

    // A storage object for quadrature point information. As opposed to
    // step-18, deal.II's native quadrature point data manager is employed here.
    CellDataStorage<typename Triangulation<dim>::cell_iterator,
                    PointHistory<dim,NumberType> > quadrature_point_history;

    // A description of the finite-element system including the displacement
    // polynomial degree, the degree-of-freedom handler, number of DoFs per
    // cell and the extractor objects used to retrieve information from the
    // solution vectors:
    const unsigned int               degree;
    const FESystem<dim>              fe;
    DoFHandler<dim>                  dof_handler_ref;
    const unsigned int               dofs_per_cell;
    const FEValuesExtractors::Vector u_fe;

    // Description of how the block-system is arranged. There is just 1 block,
    // that contains a vector DOF $\mathbf{u}$.
    // There are two reasons that we retain the block system in this problem.
    // The first is pure laziness to perform further modifications to the
    // code from which this work originated. The second is that a block system
    // would typically necessary when extending this code to multiphysics
    // problems.
    static const unsigned int        n_blocks = 1;
    static const unsigned int        n_components = dim;
    static const unsigned int        first_u_component = 0;

    enum
    {
      u_dof = 0
    };

    std::vector<types::global_dof_index>  dofs_per_block;

    // Rules for Gauss-quadrature on both the cell and faces. The number of
    // quadrature points on both cells and faces is recorded.
    const QGauss<dim>                qf_cell;
    const QGauss<dim - 1>            qf_face;
    const unsigned int               n_q_points;
    const unsigned int               n_q_points_f;

    // Objects that store the converged solution and right-hand side vectors,
    // as well as the tangent matrix. There is a AffineConstraints<double> object used
    // to keep track of constraints.  We make use of a sparsity pattern
    // designed for a block system.
    AffineConstraints<double>                 constraints;
    BlockSparsityPattern             sparsity_pattern;
    BlockSparseMatrix<double>        tangent_matrix;
    BlockVector<double>              system_rhs;
    BlockVector<double>              solution_n;

    // Then define a number of variables to store norms and update norms and
    // normalisation factors.
    struct Errors
    {
      Errors()
        :
        norm(1.0), u(1.0)
      {}

      void reset()
      {
        norm = 1.0;
        u = 1.0;
      }
      void normalise(const Errors &rhs)
      {
        if (rhs.norm != 0.0)
          norm /= rhs.norm;
        if (rhs.u != 0.0)
          u /= rhs.u;
      }

      double norm, u;
    };

    Errors error_residual, error_residual_0, error_residual_norm, error_update,
           error_update_0, error_update_norm;

    // Methods to calculate error measures
    void
    get_error_residual(Errors &error_residual);

    void
    get_error_update(const BlockVector<double> &newton_update,
                     Errors &error_update);

    // Print information to screen in a pleasing way...
    static
    void
    print_conv_header();

    void
    print_conv_footer();

    void
    print_vertical_tip_displacement();
  };

// @sect3{Implementation of the <code>Solid</code> class}

// @sect4{Public interface}

// We initialise the Solid class using data extracted from the parameter file.
  template <int dim,typename NumberType>
  Solid<dim,NumberType>::Solid(const Parameters::AllParameters &parameters)
    :
    parameters(parameters),
    vol_reference (0.0),
    vol_current (0.0),
    triangulation(Triangulation<dim>::maximum_smoothing),
    time(parameters.end_time, parameters.delta_t),
    timer(std::cout,
          TimerOutput::summary,
          TimerOutput::wall_times),
    degree(parameters.poly_degree),
    // The Finite Element System is composed of dim continuous displacement
    // DOFs.
    fe(FE_Q<dim>(parameters.poly_degree), dim), // displacement
    dof_handler_ref(triangulation),
    dofs_per_cell (fe.dofs_per_cell),
    u_fe(first_u_component),
    dofs_per_block(n_blocks),
    qf_cell(parameters.quad_order),
    qf_face(parameters.quad_order),
    n_q_points (qf_cell.size()),
    n_q_points_f (qf_face.size())
  {

  }

// The class destructor simply clears the data held by the DOFHandler
  template <int dim,typename NumberType>
  Solid<dim,NumberType>::~Solid()
  {
    dof_handler_ref.clear();
  }


// In solving the quasi-static problem, the time becomes a loading parameter,
// i.e. we increasing the loading linearly with time, making the two concepts
// interchangeable. We choose to increment time linearly using a constant time
// step size.
//
// We start the function with preprocessing, and then output the initial grid
// before starting the simulation proper with the first time (and loading)
// increment.
//
  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::run()
  {
    make_grid();
    system_setup();
    output_results();
    time.increment();

    // We then declare the incremental solution update $\varDelta
    // \mathbf{\Xi}:= \{\varDelta \mathbf{u}\}$ and start the loop over the
    // time domain.
    //
    // At the beginning, we reset the solution update for this time step...
    BlockVector<double> solution_delta(dofs_per_block);
    while (time.current() <= time.end())
      {
        solution_delta = 0.0;

        // ...solve the current time step and update total solution vector
        // $\mathbf{\Xi}_{\textrm{n}} = \mathbf{\Xi}_{\textrm{n-1}} +
        // \varDelta \mathbf{\Xi}$...
        solve_nonlinear_timestep(solution_delta);
        solution_n += solution_delta;

        // ...and plot the results before moving on happily to the next time
        // step:
        output_results();
        time.increment();
      }

    // Lastly, we print the vertical tip displacement of the Cook cantilever
    // after the full load is applied
    print_vertical_tip_displacement();
  }


// @sect3{Private interface}

// @sect4{Solid::make_grid}

// On to the first of the private member functions. Here we create the
// triangulation of the domain, for which we choose a scaled an anisotripically
// discretised rectangle which is subsequently transformed into the correct
// of the Cook cantilever. Each relevant boundary face is then given a boundary
// ID number.
//
// We then determine the volume of the reference configuration and print it
// for comparison.

  template <int dim>
  Point<dim> grid_y_transform (const Point<dim> &pt_in)
  {
    const double &x = pt_in[0];
    const double &y = pt_in[1];

    const double y_upper = 44.0 + (16.0/48.0)*x; // Line defining upper edge of beam
    const double y_lower =  0.0 + (44.0/48.0)*x; // Line defining lower edge of beam
    const double theta = y/44.0; // Fraction of height along left side of beam
    const double y_transform = (1-theta)*y_lower + theta*y_upper; // Final transformation

    Point<dim> pt_out = pt_in;
    pt_out[1] = y_transform;

    return pt_out;
  }

  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::make_grid()
  {
    // Divide the beam, but only along the x- and y-coordinate directions
    std::vector< unsigned int > repetitions(dim, parameters.elements_per_edge);
    // Only allow one element through the thickness
    // (modelling a plane strain condition)
    if (dim == 3)
      repetitions[dim-1] = 1;

    const Point<dim> bottom_left = (dim == 3 ? Point<dim>(0.0, 0.0, -0.5) : Point<dim>(0.0, 0.0));
    const Point<dim> top_right = (dim == 3 ? Point<dim>(48.0, 44.0, 0.5) : Point<dim>(48.0, 44.0));

    GridGenerator::subdivided_hyper_rectangle(triangulation,
                                              repetitions,
                                              bottom_left,
                                              top_right);

    // Since we wish to apply a Neumann BC to the right-hand surface, we
    // must find the cell faces in this part of the domain and mark them with
    // a distinct boundary ID number.  The faces we are looking for are on the
    // +x surface and will get boundary ID 11.
    // Dirichlet boundaries exist on the left-hand face of the beam (this fixed
    // boundary will get ID 1) and on the +Z and -Z faces (which correspond to
    // ID 2 and we will use to impose the plane strain condition)
    const double tol_boundary = 1e-6;
    typename Triangulation<dim>::active_cell_iterator cell =
      triangulation.begin_active(), endc = triangulation.end();
    for (; cell != endc; ++cell)
      for (unsigned int face = 0;
           face < GeometryInfo<dim>::faces_per_cell; ++face)
        if (cell->face(face)->at_boundary() == true)
          {
            if (std::abs(cell->face(face)->center()[0] - 0.0) < tol_boundary)
              cell->face(face)->set_boundary_id(1); // -X faces
            else if (std::abs(cell->face(face)->center()[0] - 48.0) < tol_boundary)
              cell->face(face)->set_boundary_id(11); // +X faces
            else if (std::abs(std::abs(cell->face(face)->center()[0]) - 0.5) < tol_boundary)
              cell->face(face)->set_boundary_id(2); // +Z and -Z faces
          }

    // Transform the hyper-rectangle into the beam shape
    GridTools::transform(&grid_y_transform<dim>, triangulation);

    GridTools::scale(parameters.scale, triangulation);

    vol_reference = GridTools::volume(triangulation);
    vol_current = vol_reference;
    std::cout << "Grid:\n\t Reference volume: " << vol_reference << std::endl;
  }


// @sect4{Solid::system_setup}

// Next we describe how the FE system is setup.  We first determine the number
// of components per block. Since the displacement is a vector component, the
// first dim components belong to it.
  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::system_setup()
  {
    timer.enter_subsection("Setup system");

    std::vector<unsigned int> block_component(n_components, u_dof); // Displacement

    // The DOF handler is then initialised and we renumber the grid in an
    // efficient manner. We also record the number of DOFs per block.
    dof_handler_ref.distribute_dofs(fe);
    DoFRenumbering::Cuthill_McKee(dof_handler_ref);
    DoFRenumbering::component_wise(dof_handler_ref, block_component);
    DoFTools::count_dofs_per_block(dof_handler_ref, dofs_per_block,
                                   block_component);

    std::cout << "Triangulation:"
              << "\n\t Number of active cells: " << triangulation.n_active_cells()
              << "\n\t Number of degrees of freedom: " << dof_handler_ref.n_dofs()
              << std::endl;

    // Setup the sparsity pattern and tangent matrix
    tangent_matrix.clear();
    {
      const types::global_dof_index n_dofs_u = dofs_per_block[u_dof];

      BlockDynamicSparsityPattern csp(n_blocks, n_blocks);

      csp.block(u_dof, u_dof).reinit(n_dofs_u, n_dofs_u);
      csp.collect_sizes();

      // Naturally, for a one-field vector-valued problem, all of the
      // components of the system are coupled.
      Table<2, DoFTools::Coupling> coupling(n_components, n_components);
      for (unsigned int ii = 0; ii < n_components; ++ii)
        for (unsigned int jj = 0; jj < n_components; ++jj)
          coupling[ii][jj] = DoFTools::always;
      DoFTools::make_sparsity_pattern(dof_handler_ref,
                                      coupling,
                                      csp,
                                      constraints,
                                      false);
      sparsity_pattern.copy_from(csp);
    }

    tangent_matrix.reinit(sparsity_pattern);

    // We then set up storage vectors
    system_rhs.reinit(dofs_per_block);
    system_rhs.collect_sizes();

    solution_n.reinit(dofs_per_block);
    solution_n.collect_sizes();

    // ...and finally set up the quadrature
    // point history:
    setup_qph();

    timer.leave_subsection();
  }


// @sect4{Solid::setup_qph}
// The method used to store quadrature information is already described in
// step-18 and step-44. Here we implement a similar setup for a SMP machine.
//
// Firstly the actual QPH data objects are created. This must be done only
// once the grid is refined to its finest level.
  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::setup_qph()
  {
    std::cout << "    Setting up quadrature point data..." << std::endl;

    quadrature_point_history.initialize(triangulation.begin_active(),
                                        triangulation.end(),
                                        n_q_points);

    // Next we setup the initial quadrature point data. Note that when
    // the quadrature point data is retrieved, it is returned as a vector
    // of smart pointers.
    for (typename Triangulation<dim>::active_cell_iterator cell =
           triangulation.begin_active(); cell != triangulation.end(); ++cell)
      {
        const std::vector<std::shared_ptr<PointHistory<dim,NumberType> > > lqph =
          quadrature_point_history.get_data(cell);
        Assert(lqph.size() == n_q_points, ExcInternalError());

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          lqph[q_point]->setup_lqp(parameters);
      }
  }


// @sect4{Solid::solve_nonlinear_timestep}

// The next function is the driver method for the Newton-Raphson scheme. At
// its top we create a new vector to store the current Newton update step,
// reset the error storage objects and print solver header.
  template <int dim,typename NumberType>
  void
  Solid<dim,NumberType>::solve_nonlinear_timestep(BlockVector<double> &solution_delta)
  {
    std::cout << std::endl << "Timestep " << time.get_timestep() << " @ "
              << time.current() << "s" << std::endl;

    BlockVector<double> newton_update(dofs_per_block);

    error_residual.reset();
    error_residual_0.reset();
    error_residual_norm.reset();
    error_update.reset();
    error_update_0.reset();
    error_update_norm.reset();

    print_conv_header();

    // We now perform a number of Newton iterations to iteratively solve the
    // nonlinear problem.  Since the problem is fully nonlinear and we are
    // using a full Newton method, the data stored in the tangent matrix and
    // right-hand side vector is not reusable and must be cleared at each
    // Newton step.  We then initially build the right-hand side vector to
    // check for convergence (and store this value in the first iteration).
    // The unconstrained DOFs of the rhs vector hold the out-of-balance
    // forces. The building is done before assembling the system matrix as the
    // latter is an expensive operation and we can potentially avoid an extra
    // assembly process by not assembling the tangent matrix when convergence
    // is attained.
    unsigned int newton_iteration = 0;
    for (; newton_iteration < parameters.max_iterations_NR;
         ++newton_iteration)
      {
        std::cout << " " << std::setw(2) << newton_iteration << " " << std::flush;

        // If we have decided that we want to continue with the iteration, we
        // assemble the tangent, make and impose the Dirichlet constraints,
        // and do the solve of the linearized system:
        make_constraints(newton_iteration);
        assemble_system(solution_delta);

        get_error_residual(error_residual);

        if (newton_iteration == 0)
          error_residual_0 = error_residual;

        // We can now determine the normalised residual error and check for
        // solution convergence:
        error_residual_norm = error_residual;
        error_residual_norm.normalise(error_residual_0);

        if (newton_iteration > 0 && error_update_norm.u <= parameters.tol_u
            && error_residual_norm.u <= parameters.tol_f)
          {
            std::cout << " CONVERGED! " << std::endl;
            print_conv_footer();

            break;
          }

        const std::pair<unsigned int, double>
        lin_solver_output = solve_linear_system(newton_update);

        get_error_update(newton_update, error_update);
        if (newton_iteration == 0)
          error_update_0 = error_update;

        // We can now determine the normalised Newton update error, and
        // perform the actual update of the solution increment for the current
        // time step, update all quadrature point information pertaining to
        // this new displacement and stress state and continue iterating:
        error_update_norm = error_update;
        error_update_norm.normalise(error_update_0);

        solution_delta += newton_update;

        std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(7)
                  << std::scientific << lin_solver_output.first << "  "
                  << lin_solver_output.second << "  " << error_residual_norm.norm
                  << "  " << error_residual_norm.u << "  "
                  << "  " << error_update_norm.norm << "  " << error_update_norm.u
                  << "  " << std::endl;
      }

    // At the end, if it turns out that we have in fact done more iterations
    // than the parameter file allowed, we raise an exception that can be
    // caught in the main() function. The call <code>AssertThrow(condition,
    // exc_object)</code> is in essence equivalent to <code>if (!cond) throw
    // exc_object;</code> but the former form fills certain fields in the
    // exception object that identify the location (filename and line number)
    // where the exception was raised to make it simpler to identify where the
    // problem happened.
    AssertThrow (newton_iteration <= parameters.max_iterations_NR,
                 ExcMessage("No convergence in nonlinear solver!"));
  }


// @sect4{Solid::print_conv_header, Solid::print_conv_footer and Solid::print_vertical_tip_displacement}

// This program prints out data in a nice table that is updated
// on a per-iteration basis. The next two functions set up the table
// header and footer:
  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::print_conv_header()
  {
    static const unsigned int l_width = 87;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;

    std::cout << "    SOLVER STEP    "
              << " |  LIN_IT   LIN_RES    RES_NORM    "
              << " RES_U     NU_NORM     "
              << " NU_U " << std::endl;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;
  }



  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::print_conv_footer()
  {
    static const unsigned int l_width = 87;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;

    std::cout << "Relative errors:" << std::endl
              << "Displacement:\t" << error_update.u / error_update_0.u << std::endl
              << "Force: \t\t" << error_residual.u / error_residual_0.u << std::endl
              << "v / V_0:\t" << vol_current << " / " << vol_reference
              << std::endl;
  }

// At the end we also output the result that can be compared to that found in
// the literature, namely the displacement at the upper right corner of the
// beam.
  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::print_vertical_tip_displacement()
  {
    static const unsigned int l_width = 87;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;

    Point<dim> soln_pt (48.0*parameters.scale,60.0*parameters.scale);
    if (dim == 3)
      soln_pt[2] = 0.5*parameters.scale;
    double vertical_tip_displacement = 0.0;
    double vertical_tip_displacement_check = 0.0;

    typename DoFHandler<dim>::active_cell_iterator cell =
      dof_handler_ref.begin_active(), endc = dof_handler_ref.end();
    for (; cell != endc; ++cell)
      {
        // if (cell->point_inside(soln_pt) == true)
        for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
          if (cell->vertex(v).distance(soln_pt) < 1e-6)
            {
              // Extract y-component of solution at the given point
              // This point is coindicent with a vertex, so we can
              // extract it directly as we're using FE_Q finite elements
              // that have support at the vertices
              vertical_tip_displacement = solution_n(cell->vertex_dof_index(v,u_dof+1));

              // Sanity check using alternate method to extract the solution
              // at the given point. To do this, we must create an FEValues instance
              // to help us extract the solution value at the desired point
              const MappingQ<dim> mapping (parameters.poly_degree);
              const Point<dim> qp_unit = mapping.transform_real_to_unit_cell(cell,soln_pt);
              const Quadrature<dim> soln_qrule (qp_unit);
              AssertThrow(soln_qrule.size() == 1, ExcInternalError());
              FEValues<dim> fe_values_soln (fe, soln_qrule, update_values);
              fe_values_soln.reinit(cell);

              // Extract y-component of solution at given point
              std::vector< Tensor<1,dim> > soln_values (soln_qrule.size());
              fe_values_soln[u_fe].get_function_values(solution_n,
                                                       soln_values);
              vertical_tip_displacement_check = soln_values[0][u_dof+1];

              break;
            }
      }
    AssertThrow(vertical_tip_displacement > 0.0, ExcMessage("Found no cell with point inside!"))

    std::cout << "Vertical tip displacement: " << vertical_tip_displacement
              << "\t Check: " << vertical_tip_displacement_check
              << std::endl;
  }


// @sect4{Solid::get_error_residual}

// Determine the true residual error for the problem.  That is, determine the
// error in the residual for the unconstrained degrees of freedom.  Note that to
// do so, we need to ignore constrained DOFs by setting the residual in these
// vector components to zero.
  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::get_error_residual(Errors &error_residual)
  {
    BlockVector<double> error_res(dofs_per_block);

    for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
      if (!constraints.is_constrained(i))
        error_res(i) = system_rhs(i);

    error_residual.norm = error_res.l2_norm();
    error_residual.u = error_res.block(u_dof).l2_norm();
  }


// @sect4{Solid::get_error_udpate}

// Determine the true Newton update error for the problem
  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::get_error_update(const BlockVector<double> &newton_update,
                                               Errors &error_update)
  {
    BlockVector<double> error_ud(dofs_per_block);
    for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
      if (!constraints.is_constrained(i))
        error_ud(i) = newton_update(i);

    error_update.norm = error_ud.l2_norm();
    error_update.u = error_ud.block(u_dof).l2_norm();
  }



// @sect4{Solid::get_total_solution}

// This function provides the total solution, which is valid at any Newton step.
// This is required as, to reduce computational error, the total solution is
// only updated at the end of the timestep.
  template <int dim,typename NumberType>
  BlockVector<double>
  Solid<dim,NumberType>::get_total_solution(const BlockVector<double> &solution_delta) const
  {
    BlockVector<double> solution_total(solution_n);
    solution_total += solution_delta;
    return solution_total;
  }


// @sect4{Solid::assemble_system}

  template <int dim,typename NumberType>
  struct Assembler_Base
  {
    virtual ~Assembler_Base() {}

    // Here we deal with the tangent matrix assembly structures. The
    // PerTaskData object stores local contributions.
    struct PerTaskData_ASM
    {
      const Solid<dim,NumberType>          *solid;
      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;

      PerTaskData_ASM(const Solid<dim,NumberType> *solid)
        :
        solid (solid),
        cell_matrix(solid->dofs_per_cell, solid->dofs_per_cell),
        cell_rhs(solid->dofs_per_cell),
        local_dof_indices(solid->dofs_per_cell)
      {}

      void reset()
      {
        cell_matrix = 0.0;
        cell_rhs = 0.0;
      }
    };

    // On the other hand, the ScratchData object stores the larger objects such as
    // the shape-function values array (<code>Nx</code>) and a shape function
    // gradient and symmetric gradient vector which we will use during the
    // assembly.
    struct ScratchData_ASM
    {
      const BlockVector<double>               &solution_total;
      std::vector<Tensor<2, dim,NumberType> >  solution_grads_u_total;

      FEValues<dim>                fe_values_ref;
      FEFaceValues<dim>            fe_face_values_ref;

      std::vector<std::vector<Tensor<2, dim,NumberType> > >         grad_Nx;
      std::vector<std::vector<SymmetricTensor<2,dim,NumberType> > > symm_grad_Nx;

      ScratchData_ASM(const FiniteElement<dim> &fe_cell,
                      const QGauss<dim> &qf_cell,
                      const UpdateFlags uf_cell,
                      const QGauss<dim-1> & qf_face,
                      const UpdateFlags uf_face,
                      const BlockVector<double> &solution_total)
        :
        solution_total(solution_total),
        solution_grads_u_total(qf_cell.size()),
        fe_values_ref(fe_cell, qf_cell, uf_cell),
        fe_face_values_ref(fe_cell, qf_face, uf_face),
        grad_Nx(qf_cell.size(),
                std::vector<Tensor<2,dim,NumberType> >(fe_cell.dofs_per_cell)),
        symm_grad_Nx(qf_cell.size(),
                     std::vector<SymmetricTensor<2,dim,NumberType> >
                     (fe_cell.dofs_per_cell))
      {}

      ScratchData_ASM(const ScratchData_ASM &rhs)
        :
        solution_total (rhs.solution_total),
        solution_grads_u_total(rhs.solution_grads_u_total),
        fe_values_ref(rhs.fe_values_ref.get_fe(),
                      rhs.fe_values_ref.get_quadrature(),
                      rhs.fe_values_ref.get_update_flags()),
        fe_face_values_ref(rhs.fe_face_values_ref.get_fe(),
                           rhs.fe_face_values_ref.get_quadrature(),
                           rhs.fe_face_values_ref.get_update_flags()),
        grad_Nx(rhs.grad_Nx),
        symm_grad_Nx(rhs.symm_grad_Nx)
      {}

      void reset()
      {
        const unsigned int n_q_points = fe_values_ref.get_quadrature().size();
        const unsigned int n_dofs_per_cell = fe_values_ref.dofs_per_cell;
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            Assert( grad_Nx[q_point].size() == n_dofs_per_cell,
                    ExcInternalError());
            Assert( symm_grad_Nx[q_point].size() == n_dofs_per_cell,
                    ExcInternalError());

            solution_grads_u_total[q_point] = Tensor<2,dim,NumberType>();
            for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
              {
                grad_Nx[q_point][k] = Tensor<2,dim,NumberType>();
                symm_grad_Nx[q_point][k] = SymmetricTensor<2,dim,NumberType>();
              }
          }
      }

    };

    // Of course, we still have to define how we assemble the tangent matrix
    // contribution for a single cell.
    void
    assemble_system_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                             ScratchData_ASM &scratch,
                             PerTaskData_ASM &data)
    {
      // Due to the C++ specialization rules, we need one more
      // level of indirection in order to define the assembly
      // routine for all different number. The next function call
      // is specialized for each NumberType, but to prevent having
      // to specialize the whole class along with it we have inlined
      // the definition of the other functions that are common to
      // all implementations.
      assemble_system_tangent_residual_one_cell(cell, scratch, data);
      assemble_neumann_contribution_one_cell(cell, scratch, data);
    }

    // This function adds the local contribution to the system matrix.
    void
    copy_local_to_global_ASM(const PerTaskData_ASM &data)
    {
      const AffineConstraints<double> &constraints = data.solid->constraints;
      BlockSparseMatrix<double> &tangent_matrix = const_cast<Solid<dim,NumberType> *>(data.solid)->tangent_matrix;
      BlockVector<double> &system_rhs =  const_cast<Solid<dim,NumberType> *>(data.solid)->system_rhs;

      constraints.distribute_local_to_global(
        data.cell_matrix, data.cell_rhs,
        data.local_dof_indices,
        tangent_matrix, system_rhs);
    }

  protected:

    // This function needs to exist in the base class for
    // Workstream to work with a reference to the base class.
    virtual void
    assemble_system_tangent_residual_one_cell(const typename DoFHandler<dim>::active_cell_iterator &/*cell*/,
                                              ScratchData_ASM &/*scratch*/,
                                              PerTaskData_ASM &/*data*/)
    {
      AssertThrow(false, ExcPureFunctionCalled());
    }

    void
    assemble_neumann_contribution_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                           ScratchData_ASM &scratch,
                                           PerTaskData_ASM &data)
    {
      // Aliases for data referenced from the Solid class
      const unsigned int &n_q_points_f = data.solid->n_q_points_f;
      const unsigned int &dofs_per_cell = data.solid->dofs_per_cell;
      const Parameters::AllParameters &parameters = data.solid->parameters;
      const Time &time = data.solid->time;
      const FESystem<dim> &fe = data.solid->fe;
      const unsigned int &u_dof = data.solid->u_dof;

      // Next we assemble the Neumann contribution. We first check to see it the
      // cell face exists on a boundary on which a traction is applied and add
      // the contribution if this is the case.
      for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
           ++face)
        if (cell->face(face)->at_boundary() == true
            && cell->face(face)->boundary_id() == 11)
          {
            scratch.fe_face_values_ref.reinit(cell, face);

            for (unsigned int f_q_point = 0; f_q_point < n_q_points_f;
                 ++f_q_point)
              {
                // We specify the traction in reference configuration.
                // For this problem, a defined total vertical force is applied
                // in the reference configuration.
                // The direction of the applied traction is assumed not to
                // evolve with the deformation of the domain.

                // Note that the contributions to the right hand side vector we
                // compute here only exist in the displacement components of the
                // vector.
                const double time_ramp = (time.current() / time.end());
                const double magnitude  = (1.0/(16.0*parameters.scale*1.0*parameters.scale))*time_ramp; // (Total force) / (RHS surface area)
                Tensor<1,dim> dir;
                dir[1] = 1.0;
                const Tensor<1, dim> traction  = magnitude*dir;

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  {
                    const unsigned int i_group =
                      fe.system_to_base_index(i).first.first;

                    if (i_group == u_dof)
                      {
                        const unsigned int component_i =
                          fe.system_to_component_index(i).first;
                        const double Ni =
                          scratch.fe_face_values_ref.shape_value(i,
                                                                 f_q_point);
                        const double JxW = scratch.fe_face_values_ref.JxW(
                                             f_q_point);

                        data.cell_rhs(i) += (Ni * traction[component_i])
                                            * JxW;
                      }
                  }
              }
          }
    }

  };

  template <int dim>
  struct Assembler<dim,double> : Assembler_Base<dim,double>
  {
    typedef double NumberType;
    using typename Assembler_Base<dim,NumberType>::ScratchData_ASM;
    using typename Assembler_Base<dim,NumberType>::PerTaskData_ASM;

    virtual ~Assembler() {}

    virtual void
    assemble_system_tangent_residual_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                              ScratchData_ASM &scratch,
                                              PerTaskData_ASM &data)
    {
      // Aliases for data referenced from the Solid class
      const unsigned int &n_q_points = data.solid->n_q_points;
      const unsigned int &dofs_per_cell = data.solid->dofs_per_cell;
      const FESystem<dim> &fe = data.solid->fe;
      const unsigned int &u_dof = data.solid->u_dof;
      const FEValuesExtractors::Vector &u_fe = data.solid->u_fe;

      data.reset();
      scratch.reset();
      scratch.fe_values_ref.reinit(cell);
      cell->get_dof_indices(data.local_dof_indices);

      const std::vector<std::shared_ptr<const PointHistory<dim,NumberType> > > lqph =
        const_cast<const Solid<dim,NumberType> *>(data.solid)->quadrature_point_history.get_data(cell);
      Assert(lqph.size() == n_q_points, ExcInternalError());

      // We first need to find the solution gradients at quadrature points
      // inside the current cell and then we update each local QP using the
      // displacement gradient:
      scratch.fe_values_ref[u_fe].get_function_gradients(scratch.solution_total,
                                                         scratch.solution_grads_u_total);

      // Now we build the local cell stiffness matrix. Since the global and
      // local system matrices are symmetric, we can exploit this property by
      // building only the lower half of the local matrix and copying the values
      // to the upper half.
      //
      // In doing so, we first extract some configuration dependent variables
      // from our QPH history objects for the current quadrature point.
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          const Tensor<2,dim,NumberType> &grad_u = scratch.solution_grads_u_total[q_point];
          const Tensor<2,dim,NumberType> F = Physics::Elasticity::Kinematics::F(grad_u);
          const NumberType               det_F = determinant(F);
          const Tensor<2,dim,NumberType> F_bar = Physics::Elasticity::Kinematics::F_iso(F);
          const SymmetricTensor<2,dim,NumberType> b_bar = Physics::Elasticity::Kinematics::b(F_bar);
          const Tensor<2,dim,NumberType> F_inv = invert(F);
          Assert(det_F > NumberType(0.0), ExcInternalError());

          for (unsigned int k = 0; k < dofs_per_cell; ++k)
            {
              const unsigned int k_group = fe.system_to_base_index(k).first.first;

              if (k_group == u_dof)
                {
                  scratch.grad_Nx[q_point][k] = scratch.fe_values_ref[u_fe].gradient(k, q_point) * F_inv;
                  scratch.symm_grad_Nx[q_point][k] = symmetrize(scratch.grad_Nx[q_point][k]);
                }
              else
                Assert(k_group <= u_dof, ExcInternalError());
            }

          const SymmetricTensor<2,dim,NumberType> tau = lqph[q_point]->get_tau(det_F,b_bar);
          const SymmetricTensor<4,dim,NumberType> Jc  = lqph[q_point]->get_Jc(det_F,b_bar);
          const Tensor<2,dim,NumberType> tau_ns (tau);

          // Next we define some aliases to make the assembly process easier to
          // follow
          const std::vector<SymmetricTensor<2, dim> > &symm_grad_Nx = scratch.symm_grad_Nx[q_point];
          const std::vector<Tensor<2, dim> > &grad_Nx = scratch.grad_Nx[q_point];
          const double JxW = scratch.fe_values_ref.JxW(q_point);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
            {
              const unsigned int component_i = fe.system_to_component_index(i).first;
              const unsigned int i_group     = fe.system_to_base_index(i).first.first;

              if (i_group == u_dof)
                data.cell_rhs(i) -= (symm_grad_Nx[i] * tau) * JxW;
              else
                Assert(i_group <= u_dof, ExcInternalError());

              for (unsigned int j = 0; j <= i; ++j)
                {
                  const unsigned int component_j = fe.system_to_component_index(j).first;
                  const unsigned int j_group     = fe.system_to_base_index(j).first.first;

                  // This is the $\mathsf{\mathbf{k}}_{\mathbf{u} \mathbf{u}}$
                  // contribution. It comprises a material contribution, and a
                  // geometrical stress contribution which is only added along
                  // the local matrix diagonals:
                  if ((i_group == j_group) && (i_group == u_dof))
                    {
                      data.cell_matrix(i, j) += symm_grad_Nx[i] * Jc // The material contribution:
                                                * symm_grad_Nx[j] * JxW;
                      if (component_i == component_j) // geometrical stress contribution
                        data.cell_matrix(i, j) += grad_Nx[i][component_i] * tau_ns
                                                  * grad_Nx[j][component_j] * JxW;
                    }
                  else
                    Assert((i_group <= u_dof) && (j_group <= u_dof),
                           ExcInternalError());
                }
            }
        }


      // Finally, we need to copy the lower half of the local matrix into the
      // upper half:
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
          data.cell_matrix(i, j) = data.cell_matrix(j, i);
    }

  };

#ifdef ENABLE_AD_FORMULATION

  namespace AD = dealii::Differentiation::AD;

  template <int dim>
  struct Assembler<dim,Sacado::Fad::DFad<double> > : Assembler_Base<dim,Sacado::Fad::DFad<double> >
  {
    using ADHelper = AD::ResidualLinearization<AD::NumberTypes::sacado_dfad,double>;
    using ADNumberType = typename ADHelper::ad_type;

    using typename Assembler_Base<dim,ADNumberType>::ScratchData_ASM;
    using typename Assembler_Base<dim,ADNumberType>::PerTaskData_ASM;

    virtual ~Assembler() {}

    virtual void
    assemble_system_tangent_residual_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                              ScratchData_ASM &scratch,
                                              PerTaskData_ASM &data)
    {
      const FEValuesExtractors::Vector &u_fe = data.solid->u_fe;
      const unsigned int &n_q_points = data.solid->n_q_points;
      const unsigned int &dofs_per_cell = data.solid->dofs_per_cell;
      const FESystem<dim> &fe = data.solid->fe;
      const unsigned int &u_dof = data.solid->u_dof;

      data.reset();
      scratch.reset();
      scratch.fe_values_ref.reinit(cell);
      cell->get_dof_indices(data.local_dof_indices);
      const unsigned int n_independent_variables =
          data.local_dof_indices.size();
      const unsigned int n_dependent_variables = dofs_per_cell;

      const std::vector<std::shared_ptr<const PointHistory<dim,ADNumberType> > > lqph =
          data.solid->quadrature_point_history.get_data(cell);
      Assert(lqph.size() == n_q_points, ExcInternalError());

      // Create and initialize an instance of the helper class.
      ADHelper ad_helper(n_independent_variables, n_dependent_variables);
      // Initialize the local data structures for assembly.
      // This is also taken care of by the ADHelper, so this step could
      // be skipped.
      Assert(data.cell_rhs.size() == n_independent_variables, ExcDimensionMismatch(data.cell_rhs.size(),n_independent_variables));
      Assert(data.cell_matrix.m() == n_independent_variables, ExcDimensionMismatch(data.cell_matrix.m(),n_independent_variables));
      Assert(data.cell_matrix.n() == n_dependent_variables, ExcDimensionMismatch(data.cell_matrix.n(),n_dependent_variables));

      // An optional call to set the amount of memory to be allocated to
      // storing taped data.
      // If using a taped AD number then we would likely want to increase
      // the buffer size from the default values as the expression for each
      // residual component will likely be lengthy, and therefore memory
      // intensive.
      ad_helper.set_tape_buffer_sizes();
      // If using a taped AD number, then at this point we would initiate
      // taping of the expression for the energy for this FE type and
      // material combination:
      // Select a tape number to record to
      const typename AD::Types<ADNumberType>::tape_index tape_index = 1;
      // Indicate that we are about to start tracing the operations for
      // function evaluation on the tape. If this tape has already been
      // used (i.e. the operations are already recorded) then we
      // (optionally) load the tape and reuse this data.
      const bool is_recording
        = ad_helper.start_recording_operations(tape_index);

      // The steps that follow in the recording phase are required for
      // tapeless methods as well.
      if (is_recording == true)
      {
        // This is the "recording" phase of the operations.
        // First, we set the values for all DoFs.
        ad_helper.register_dof_values(scratch.solution_total, data.local_dof_indices);
        // Then we get the complete set of degree of freedom values as
        // represented by auto-differentiable numbers. The operations
        // performed with these variables are tracked by the AD library
        // from this point until stop_recording_operations() is called.
        const std::vector<ADNumberType> &dof_values_ad
          = ad_helper.get_sensitive_dof_values();

        // Then we do some problem specific tasks, the first being to
        // compute all values, gradients, etc. based on sensitive AD DoF
        // values. Here we are fetching the displacement gradients at each
        // quadrature point.
        std::vector<Tensor<2, dim, ADNumberType>> Grad_u(
          n_q_points, Tensor<2, dim, ADNumberType>());
        scratch.fe_values_ref[u_fe].get_function_gradients_from_local_dof_values(
          dof_values_ad, Grad_u);

        // This variable stores the cell residual vector contributions.
        // IMPORTANT: Note that each entry is hand-initialized with a value
        // of zero. This is a highly recommended practise, as some AD
        // numbers appear not to safely initialize their internal data
        // structures.
        std::vector<ADNumberType> residual_ad (
          n_dependent_variables, ADNumberType(0.0));
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          // Calculate the deformation gradient at this quadrature point
          const Tensor<2,dim,ADNumberType> F = Physics::Elasticity::Kinematics::F(Grad_u[q_point]);
          Assert(numbers::value_is_greater_than(determinant(F), 0.0),
                 ExcMessage("Negative determinant of the deformation "
                            "gradient detected!"));
          const ADNumberType               det_F = determinant(F);
          const Tensor<2,dim,ADNumberType> F_bar = Physics::Elasticity::Kinematics::F_iso(F);
          const SymmetricTensor<2,dim,ADNumberType> b_bar = Physics::Elasticity::Kinematics::b(F_bar);
          Assert(det_F > ADNumberType(0.0), ExcInternalError());
          const Tensor<2,dim,ADNumberType> F_inv = invert(F);

          for (unsigned int k = 0; k < dofs_per_cell; ++k)
          {
            const unsigned int k_group = fe.system_to_base_index(k).first.first;

            if (k_group == u_dof)
            {
              scratch.grad_Nx[q_point][k] = scratch.fe_values_ref[u_fe].gradient(k, q_point) * F_inv;
              scratch.symm_grad_Nx[q_point][k] = symmetrize(scratch.grad_Nx[q_point][k]);
            }
            else
              Assert(k_group <= u_dof, ExcInternalError());
          }

          const SymmetricTensor<2,dim,ADNumberType> tau = lqph[q_point]->get_tau(det_F,b_bar);

          // Next we define some position-dependent aliases, again to
          // make the assembly process easier to follow.
          const std::vector<SymmetricTensor<2, dim,ADNumberType> > &symm_grad_Nx = scratch.symm_grad_Nx[q_point];
          const double JxW = scratch.fe_values_ref.JxW(q_point);

          for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const unsigned int i_group     = fe.system_to_base_index(i).first.first;

            if (i_group == u_dof)
              residual_ad[i] += (symm_grad_Nx[i] * tau) * JxW;
            else
              Assert(i_group <= u_dof, ExcInternalError());
          }
        }

        // Register the definition of the cell residual
        ad_helper.register_residual_vector(residual_ad);
        // Indicate that we have completed tracing the operations onto
        // the tape.
        ad_helper.stop_recording_operations(false); // write_tapes_to_file
      }
      else
      {
        // This is the "tape reuse" phase of the operations.
        // Here we will leverage the already traced operations that reside
        // on a tape, and simply re-evaluate the tape at a different point
        // to get the function values and their derivatives.
        // Load the existing tape to be reused
        ad_helper.activate_recorded_tape(tape_index);
        // Set the new values of the independent variables where the
        // recorded dependent functions are to be evaluated (and
        // differentiated around).
        ad_helper.set_dof_values(scratch.solution_total, data.local_dof_indices);
      }

      // Compute the residual values and their Jacobian at the
      // evaluation point
      ad_helper.compute_residual(data.cell_rhs);
      data.cell_rhs *= -1.0; // RHS = - residual
      ad_helper.compute_linearization(data.cell_matrix);
    }

  };


  template <int dim>
  struct Assembler<dim,Sacado::Rad::ADvar<Sacado::Fad::DFad<double> > > : Assembler_Base<dim,Sacado::Rad::ADvar<Sacado::Fad::DFad<double> > >
  {
    using ADHelper = AD::EnergyFunctional<AD::NumberTypes::sacado_rad_dfad,double>;
    using ADNumberType = typename ADHelper::ad_type;

    using typename Assembler_Base<dim,ADNumberType>::ScratchData_ASM;
    using typename Assembler_Base<dim,ADNumberType>::PerTaskData_ASM;

    virtual ~Assembler() {}

    virtual void
    assemble_system_tangent_residual_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                              ScratchData_ASM &scratch,
                                              PerTaskData_ASM &data)
    {
      const FEValuesExtractors::Vector &u_fe = data.solid->u_fe;
      const unsigned int &n_q_points = data.solid->n_q_points;

      data.reset();
      scratch.reset();
      scratch.fe_values_ref.reinit(cell);
      cell->get_dof_indices(data.local_dof_indices);
      const unsigned int n_independent_variables =
          data.local_dof_indices.size();

      const std::vector<std::shared_ptr<const PointHistory<dim,ADNumberType> > > lqph =
          data.solid->quadrature_point_history.get_data(cell);
      Assert(lqph.size() == n_q_points, ExcInternalError());

      // Create and initialize an instance of the helper class.
      ADHelper ad_helper(n_independent_variables);
      // Initialize the local data structures for assembly.
      // This is also taken care of by the ADHelper, so this step could
      // be skipped.
      Assert(data.cell_rhs.size() == n_independent_variables, ExcDimensionMismatch(data.cell_rhs.size(),n_independent_variables));
      Assert(data.cell_matrix.m() == n_independent_variables, ExcDimensionMismatch(data.cell_matrix.m(),n_independent_variables));
      Assert(data.cell_matrix.n() == n_independent_variables, ExcDimensionMismatch(data.cell_matrix.n(),n_independent_variables));

      // An optional call to set the amount of memory to be allocated to
      // storing taped data.
      // If using a taped AD number then we would likely want to increase
      // the buffer size from the default values as the expression for each
      // residual component will likely be lengthy, and therefore memory
      // intensive.
      ad_helper.set_tape_buffer_sizes();
      // If using a taped AD number, then at this point we would initiate
      // taping of the expression for the energy for this FE type and
      // material combination:
      // Select a tape number to record to
      const typename AD::Types<ADNumberType>::tape_index tape_index = 1;
      // Indicate that we are about to start tracing the operations for
      // function evaluation on the tape. If this tape has already been
      // used (i.e. the operations are already recorded) then we
      // (optionally) load the tape and reuse this data.
      const bool is_recording
        = ad_helper.start_recording_operations(tape_index);

      // The steps that follow in the recording phase are required for
      // tapeless methods as well.
      if (is_recording == true)
      {
        // This is the "recording" phase of the operations.
        // First, we set the values for all DoFs.
        ad_helper.register_dof_values(scratch.solution_total, data.local_dof_indices);
        // Then we get the complete set of degree of freedom values as
        // represented by auto-differentiable numbers. The operations
        // performed with these variables are tracked by the AD library
        // from this point until stop_recording_operations() is called.
        const std::vector<ADNumberType> &dof_values_ad
          = ad_helper.get_sensitive_dof_values();

        // Then we do some problem specific tasks, the first being to
        // compute all values, gradients, etc. based on sensitive AD DoF
        // values. Here we are fetching the displacement gradients at each
        // quadrature point.
        std::vector<Tensor<2, dim, ADNumberType>> Grad_u(
          n_q_points, Tensor<2, dim, ADNumberType>());
        scratch.fe_values_ref[u_fe].get_function_gradients_from_local_dof_values(
          dof_values_ad, Grad_u);

        // This variable stores the cell total energy.
        // IMPORTANT: Note that it is hand-initialized with a value of
        // zero. This is a highly recommended practise, as some AD numbers
        // appear not to safely initialize their internal data structures.
        ADNumberType energy_ad = ADNumberType(0.0);
        // Compute the cell total energy = (internal + external) energies
        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          // Calculate the deformation gradient at this quadrature point
          const Tensor<2,dim,ADNumberType> F = Physics::Elasticity::Kinematics::F(Grad_u[q_point]);
          Assert(numbers::value_is_greater_than(determinant(F), 0.0),
                 ExcMessage("Negative determinant of the deformation "
                            "gradient detected!"));
          const ADNumberType               det_F = determinant(F);
          const Tensor<2,dim,ADNumberType> F_bar = Physics::Elasticity::Kinematics::F_iso(F);
          const SymmetricTensor<2,dim,ADNumberType> b_bar = Physics::Elasticity::Kinematics::b(F_bar);
          Assert(det_F > ADNumberType(0.0), ExcInternalError());

          // Next we define some position-dependent aliases, again to
          // make the assembly process easier to follow.
          const double JxW = scratch.fe_values_ref.JxW(q_point);

          const ADNumberType Psi = lqph[q_point]->get_Psi(det_F,b_bar);

          // We extract the configuration-dependent material energy
          // from our QPH history objects for the current quadrature point
          // and integrate its contribution to increment the total
          // cell energy.
          energy_ad += Psi * JxW;
        }

        // Register the definition of the total cell energy
        ad_helper.register_energy_functional(energy_ad);
        // Indicate that we have completed tracing the operations onto
        // the tape.
        ad_helper.stop_recording_operations(false); // write_tapes_to_file
      }
      else
      {
        // This is the "tape reuse" phase of the operations.
        // Here we will leverage the already traced operations that reside
        // on a tape, and simply re-evaluate the tape at a different point
        // to get the function values and their derivatives.
        // Load the existing tape to be reused
        ad_helper.activate_recorded_tape(tape_index);
        // Set the new values of the independent variables where the
        // recorded dependent functions are to be evaluated (and
        // differentiated around).
        ad_helper.set_dof_values(scratch.solution_total, data.local_dof_indices);
      }

      // Compute the residual values and their Jacobian at the
      // evaluation point
      ad_helper.compute_residual(data.cell_rhs);
      data.cell_rhs *= -1.0; // RHS = - residual
      ad_helper.compute_linearization(data.cell_matrix);
    }

  };


#endif


// Since we use TBB for assembly, we simply setup a copy of the
// data structures required for the process and pass them, along
// with the memory addresses of the assembly functions to the
// WorkStream object for processing. Note that we must ensure that
// the matrix is reset before any assembly operations can occur.
  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::assemble_system(const BlockVector<double> &solution_delta)
  {
    timer.enter_subsection("Assemble linear system");
    std::cout << " ASM " << std::flush;

    tangent_matrix = 0.0;
    system_rhs = 0.0;

    const UpdateFlags uf_cell(update_gradients |
                              update_JxW_values);
    const UpdateFlags uf_face(update_values |
                              update_JxW_values);

    const BlockVector<double> solution_total(get_total_solution(solution_delta));
    typename Assembler_Base<dim,NumberType>::PerTaskData_ASM per_task_data(this);
    typename Assembler_Base<dim,NumberType>::ScratchData_ASM scratch_data(fe, qf_cell, uf_cell, qf_face, uf_face, solution_total);
    Assembler<dim,NumberType> assembler;

    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    static_cast<Assembler_Base<dim,NumberType>&>(assembler),
                    &Assembler_Base<dim,NumberType>::assemble_system_one_cell,
                    &Assembler_Base<dim,NumberType>::copy_local_to_global_ASM,
                    scratch_data,
                    per_task_data);

    timer.leave_subsection();
  }


// @sect4{Solid::make_constraints}
// The constraints for this problem are simple to describe.
// However, since we are dealing with an iterative Newton method,
// it should be noted that any displacement constraints should only
// be specified at the zeroth iteration and subsequently no
// additional contributions are to be made since the constraints
// are already exactly satisfied.
  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::make_constraints(const int &it_nr)
  {
    std::cout << " CST " << std::flush;

    // Since the constraints are different at different Newton iterations, we
    // need to clear the constraints matrix and completely rebuild
    // it. However, after the first iteration, the constraints remain the same
    // and we can simply skip the rebuilding step if we do not clear it.
    if (it_nr > 1)
      return;
    constraints.clear();
    const bool apply_dirichlet_bc = (it_nr == 0);

    // The boundary conditions for the indentation problem are as follows: On
    // the -x, -y and -z faces (ID's 0,2,4) we set up a symmetry condition to
    // allow only planar movement while the +x and +y faces (ID's 1,3) are
    // traction free. In this contrived problem, part of the +z face (ID 5) is
    // set to have no motion in the x- and y-component. Finally, as described
    // earlier, the other part of the +z face has an the applied pressure but
    // is also constrained in the x- and y-directions.
    //
    // In the following, we will have to tell the function interpolation
    // boundary values which components of the solution vector should be
    // constrained (i.e., whether it's the x-, y-, z-displacements or
    // combinations thereof). This is done using ComponentMask objects (see
    // @ref GlossComponentMask) which we can get from the finite element if we
    // provide it with an extractor object for the component we wish to
    // select. To this end we first set up such extractor objects and later
    // use it when generating the relevant component masks:

    // Fixed left hand side of the beam
    {
      const int boundary_id = 1;

      if (apply_dirichlet_bc == true)
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 fe.component_mask(u_fe));
      else
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 fe.component_mask(u_fe));
    }

    // Zero Z-displacement through thickness direction
    // This corresponds to a plane strain condition being imposed on the beam
    if (dim == 3)
      {
        const int boundary_id = 2;
        const FEValuesExtractors::Scalar z_displacement(2);

        if (apply_dirichlet_bc == true)
          VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                   boundary_id,
                                                   ZeroFunction<dim>(n_components),
                                                   constraints,
                                                   fe.component_mask(z_displacement));
        else
          VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                   boundary_id,
                                                   ZeroFunction<dim>(n_components),
                                                   constraints,
                                                   fe.component_mask(z_displacement));
      }

    constraints.close();
  }

// @sect4{Solid::solve_linear_system}
// As the system is composed of a single block, defining a solution scheme
// for the linear problem is straight-forward.
  template <int dim,typename NumberType>
  std::pair<unsigned int, double>
  Solid<dim,NumberType>::solve_linear_system(BlockVector<double> &newton_update)
  {
    BlockVector<double> A(dofs_per_block);
    BlockVector<double> B(dofs_per_block);

    unsigned int lin_it = 0;
    double lin_res = 0.0;

    // We solve for the incremental displacement $d\mathbf{u}$.
    {
      timer.enter_subsection("Linear solver");
      std::cout << " SLV " << std::flush;
      if (parameters.type_lin == "CG")
        {
          const int solver_its = tangent_matrix.block(u_dof, u_dof).m();
          const double tol_sol = parameters.tol_lin
                                 * system_rhs.block(u_dof).l2_norm();

          SolverControl solver_control(solver_its, tol_sol);

          GrowingVectorMemory<Vector<double> > GVM;
          SolverCG<Vector<double> > solver_CG(solver_control, GVM);

          // We've chosen by default a SSOR preconditioner as it appears to
          // provide the fastest solver convergence characteristics for this
          // problem on a single-thread machine.  However, for multicore
          // computing, the Jacobi preconditioner which is multithreaded may
          // converge quicker for larger linear systems.
          PreconditionSelector<SparseMatrix<double>, Vector<double> >
          preconditioner (parameters.preconditioner_type,
                          parameters.preconditioner_relaxation);
          preconditioner.use_matrix(tangent_matrix.block(u_dof, u_dof));

          solver_CG.solve(tangent_matrix.block(u_dof, u_dof),
                          newton_update.block(u_dof),
                          system_rhs.block(u_dof),
                          preconditioner);

          lin_it = solver_control.last_step();
          lin_res = solver_control.last_value();
        }
      else if (parameters.type_lin == "Direct")
        {
          // Otherwise if the problem is small
          // enough, a direct solver can be
          // utilised.
          SparseDirectUMFPACK A_direct;
          A_direct.initialize(tangent_matrix.block(u_dof, u_dof));
          A_direct.vmult(newton_update.block(u_dof), system_rhs.block(u_dof));

          lin_it = 1;
          lin_res = 0.0;
        }
      else
        Assert (false, ExcMessage("Linear solver type not implemented"));

      timer.leave_subsection();
    }

    // Now that we have the displacement update, distribute the constraints
    // back to the Newton update:
    constraints.distribute(newton_update);

    return std::make_pair(lin_it, lin_res);
  }

// @sect4{Solid::output_results}
// Here we present how the results are written to file to be viewed
// using ParaView or Visit. The method is similar to that shown in the
// tutorials so will not be discussed in detail.
  template <int dim,typename NumberType>
  void Solid<dim,NumberType>::output_results() const
  {
    DataOut<dim> data_out;
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(dim,
                                  DataComponentInterpretation::component_is_part_of_vector);

    std::vector<std::string> solution_name(dim, "displacement");

    data_out.attach_dof_handler(dof_handler_ref);
    data_out.add_data_vector(solution_n,
                             solution_name,
                             DataOut<dim>::type_dof_data,
                             data_component_interpretation);

    // Since we are dealing with a large deformation problem, it would be nice
    // to display the result on a displaced grid!  The MappingQEulerian class
    // linked with the DataOut class provides an interface through which this
    // can be achieved without physically moving the grid points in the
    // Triangulation object ourselves.  We first need to copy the solution to
    // a temporary vector and then create the Eulerian mapping. We also
    // specify the polynomial degree to the DataOut object in order to produce
    // a more refined output data set when higher order polynomials are used.
    Vector<double> soln(solution_n.size());
    for (unsigned int i = 0; i < soln.size(); ++i)
      soln(i) = solution_n(i);
    MappingQEulerian<dim> q_mapping(degree, dof_handler_ref, soln);
    data_out.build_patches(q_mapping, degree);

    std::ostringstream filename;
    filename << "solution-" << time.get_timestep() << ".vtk";

    std::ofstream output(filename.str().c_str());
    data_out.write_vtk(output);
  }

}


// @sect3{Main function}
// Lastly we provide the main driver function which appears
// no different to the other tutorials.
int main (int argc, char *argv[])
{
  using namespace dealii;
  using namespace Cook_Membrane;

  const unsigned int dim = 2;
  try
    {
      deallog.depth_console(0);
      Parameters::AllParameters parameters("parameters.prm");
      if (parameters.automatic_differentiation_order == 0)
        {
          std::cout << "Assembly method: Residual and linearisation are computed manually." << std::endl;

          // Allow multi-threading
          Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,
                                                              dealii::numbers::invalid_unsigned_int);

          typedef double NumberType;
          Solid<dim,NumberType> solid_3d(parameters);
          solid_3d.run();
        }
#ifdef ENABLE_AD_FORMULATION
      else if (parameters.automatic_differentiation_order == 1)
        {
          std::cout << "Assembly method: Residual computed manually; linearisation performed using AD." << std::endl;

#ifdef DEAL_II_WITH_TRILINOS
          if (parameters.automatic_differentiation_library == "Sacado")
          {
            // Allow multi-threading
            Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,
                                                                dealii::numbers::invalid_unsigned_int);

            typedef Sacado::Fad::DFad<double> NumberType;
            Solid<dim,NumberType> solid_3d(parameters);
            solid_3d.run();
          }
          else
#endif
#ifdef DEAL_II_WITH_ADOLC
          if (parameters.automatic_differentiation_library == "ADOL-C")
          {
            AssertThrow(false, ExcNotImplemented());
            // ADOL-C is not thread-safe, so disable threading.
            // Parallisation using MPI would be possible though.
//            Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
//
//            typedef Sacado::Fad::DFad<double> NumberType;
//            Solid<dim,NumberType> solid_3d(parameters);
//            solid_3d.run();
          }
          else
#endif
          {
            AssertThrow(false, ExcMessage("The selected auto-differentiable library is not supported."));
          }
        }
      else if (parameters.automatic_differentiation_order == 2)
        {
          std::cout << "Assembly method: Residual and linearisation computed using AD." << std::endl;

#ifdef DEAL_II_WITH_TRILINOS

          if (parameters.automatic_differentiation_library == "Sacado")
          {
            // Sacado Rad-Fad is not thread-safe, so disable threading.
            // Parallisation using MPI would be possible though.
            Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,
                                                                1);

            typedef Sacado::Rad::ADvar< Sacado::Fad::DFad<double> > NumberType;
            Solid<dim,NumberType> solid_3d(parameters);
            solid_3d.run();
          }
          else
#endif
#ifdef DEAL_II_WITH_ADOLC
          if (parameters.automatic_differentiation_library == "ADOL-C")
          {
            AssertThrow(false, ExcNotImplemented());
//            // ADOL-C is not thread-safe, so disable threading.
//            // Parallisation using MPI would be possible though.
//            Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,
//                                                                1);
//
//            typedef Sacado::Rad::ADvar< Sacado::Fad::DFad<double> > NumberType;
//            Solid<dim,NumberType> solid_3d(parameters);
//            solid_3d.run();
          }
          else
#endif
          {
            AssertThrow(false, ExcMessage("The selected auto-differentiable library is not supported."));
          }
        }
#endif
      else
        {
          AssertThrow(false,
                      ExcMessage("The selected assembly method is not supported. "
                                 "You need deal.II 9.0 and Trilinos with the Sacado package "
                                 "to enable assembly using automatic differentiation."));
        }
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
