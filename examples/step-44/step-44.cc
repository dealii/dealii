/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2010 - 2013 by the deal.II authors and
 *                              & Jean-Paul Pelteret and Andrew McBride
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
 * Authors: Jean-Paul Pelteret, University of Cape Town,
 *          Andrew McBride, University of Erlangen-Nuremberg, 2010
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

#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition_selector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_direct.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <iostream>
#include <fstream>


// We then stick everything that relates to this tutorial program into a
// namespace of its own, and import all the deal.II function and class names
// into it:
namespace Step44
{
  using namespace dealii;

// @sect3{Run-time parameters}
//
// There are several parameters that can be set in the code so we set up a
// ParameterHandler object to read in the choices at run-time.
  namespace Parameters
  {
// @sect4{Finite Element system}

// As mentioned in the introduction, a different order interpolation should be
// used for the displacement $\mathbf{u}$ than for the pressure
// $\widetilde{p}$ and the dilatation $\widetilde{J}$.  Choosing
// $\widetilde{p}$ and $\widetilde{J}$ as discontinuous (constant) functions
// at the element level leads to the mean-dilatation method. The discontinuous
// approximation allows $\widetilde{p}$ and $\widetilde{J}$ to be condensed
// out and a classical displacement based method is recovered.  Here we
// specify the polynomial order used to approximate the solution.  The
// quadrature order should be adjusted accordingly.
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

// Make adjustments to the problem geometry and the applied load.  Since the
// problem modelled here is quite specific, the load scale can be altered to
// specific values to compare with the results given in the literature.
    struct Geometry
    {
      unsigned int global_refinement;
      double       scale;
      double       p_p0;

      static void
      declare_parameters(ParameterHandler &prm);

      void
      parse_parameters(ParameterHandler &prm);
    };

    void Geometry::declare_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry");
      {
        prm.declare_entry("Global refinement", "2",
                          Patterns::Integer(0),
                          "Global refinement level");

        prm.declare_entry("Grid scale", "1e-3",
                          Patterns::Double(0.0),
                          "Global grid scaling factor");

        prm.declare_entry("Pressure ratio p/p0", "100",
                          Patterns::Selection("20|40|60|80|100"),
                          "Ratio of applied pressure to reference pressure");
      }
      prm.leave_subsection();
    }

    void Geometry::parse_parameters(ParameterHandler &prm)
    {
      prm.enter_subsection("Geometry");
      {
        global_refinement = prm.get_integer("Global refinement");
        scale = prm.get_double("Grid scale");
        p_p0 = prm.get_double("Pressure ratio p/p0");
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
        prm.declare_entry("Poisson's ratio", "0.4999",
                          Patterns::Double(-1.0,0.5),
                          "Poisson's ratio");

        prm.declare_entry("Shear modulus", "80.194e6",
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
      double      max_iterations_lin;
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

        prm.declare_entry("Max iteration multiplier", "1",
                          Patterns::Double(0.0),
                          "Linear solver iterations (multiples of the system matrix size)");

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
        max_iterations_lin = prm.get_double("Max iteration multiplier");
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
    struct AllParameters : public FESystem,
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
      prm.read_input(input_file);
      parse_parameters(prm);
    }

    void AllParameters::declare_parameters(ParameterHandler &prm)
    {
      FESystem::declare_parameters(prm);
      Geometry::declare_parameters(prm);
      Materials::declare_parameters(prm);
      LinearSolver::declare_parameters(prm);
      NonlinearSolver::declare_parameters(prm);
      Time::declare_parameters(prm);
    }

    void AllParameters::parse_parameters(ParameterHandler &prm)
    {
      FESystem::parse_parameters(prm);
      Geometry::parse_parameters(prm);
      Materials::parse_parameters(prm);
      LinearSolver::parse_parameters(prm);
      NonlinearSolver::parse_parameters(prm);
      Time::parse_parameters(prm);
    }
  }

// @sect3{Some standard tensors}

// Now we define some frequently used second and fourth-order tensors:
  template <int dim>
  class StandardTensors
  {
  public:

    // $\mathbf{I}$
    static const SymmetricTensor<2, dim> I;
    // $\mathbf{I} \otimes \mathbf{I}$
    static const SymmetricTensor<4, dim> IxI;
    // $\mathcal{S}$, note that as we only use this fourth-order unit tensor
    // to operate on symmetric second-order tensors.  To maintain notation
    // consistent with Holzapfel (2001) we name the tensor $\mathcal{I}$
    static const SymmetricTensor<4, dim> II;
    // Fourth-order deviatoric tensor such that
    // $\textrm{dev} \{ \bullet \} = \{ \bullet \} -
    //  [1/\textrm{dim}][ \{ \bullet\} :\mathbf{I}]\mathbf{I}$
    static const SymmetricTensor<4, dim> dev_P;
  };

  template <int dim>
  const SymmetricTensor<2, dim>
  StandardTensors<dim>::I = unit_symmetric_tensor<dim>();

  template <int dim>
  const SymmetricTensor<4, dim>
  StandardTensors<dim>::IxI = outer_product(I, I);

  template <int dim>
  const SymmetricTensor<4, dim>
  StandardTensors<dim>::II = identity_tensor<dim>();

  template <int dim>
  const SymmetricTensor<4, dim>
  StandardTensors<dim>::dev_P = deviator_tensor<dim>();

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

// @sect3{Compressible neo-Hookean material within a three-field formulation}

// As discussed in the Introduction, Neo-Hookean materials are a type of
// hyperelastic materials.  The entire domain is assumed to be composed of a
// compressible neo-Hookean material.  This class defines the behaviour of
// this material within a three-field formulation.  Compressible neo-Hookean
// materials can be described by a strain-energy function (SEF) $ \Psi =
// \Psi_{\text{iso}}(\overline{\mathbf{b}}) + \Psi_{\text{vol}}(\widetilde{J})
// $.
//
// The isochoric response is given by $
// \Psi_{\text{iso}}(\overline{\mathbf{b}}) = c_{1} [\overline{I}_{1} - 3] $
// where $ c_{1} = \frac{\mu}{2} $ and $\overline{I}_{1}$ is the first
// invariant of the left- or right-isochoric Cauchy-Green deformation tensors.
// That is $\overline{I}_1 :=\textrm{tr}(\overline{\mathbf{b}})$.  In this
// example the SEF that governs the volumetric response is defined as $
// \Psi_{\text{vol}}(\widetilde{J}) = \kappa \frac{1}{4} [ \widetilde{J}^2 - 1
// - 2\textrm{ln}\; \widetilde{J} ]$,  where $\kappa:= \lambda + 2/3 \mu$ is
// the <a href="http://en.wikipedia.org/wiki/Bulk_modulus">bulk modulus</a>
// and $\lambda$ is <a
// href="http://en.wikipedia.org/wiki/Lam%C3%A9_parameters">Lame's first
// parameter</a>.
//
// The following class will be used to characterize the material we work with,
// and provides a central point that one would need to modify if one were to
// implement a different material model. For it to work, we will store one
// object of this type per quadrature point, and in each of these objects
// store the current state (characterized by the values or measures  of the three fields)
// so that we can compute the elastic coefficients linearized around the
// current state.
  template <int dim>
  class Material_Compressible_Neo_Hook_Three_Field
  {
  public:
    Material_Compressible_Neo_Hook_Three_Field(const double mu,
                                               const double nu)
      :
      kappa((2.0 * mu * (1.0 + nu)) / (3.0 * (1.0 - 2.0 * nu))),
      c_1(mu / 2.0),
      det_F(1.0),
      p_tilde(0.0),
      J_tilde(1.0),
      b_bar(StandardTensors<dim>::I)
    {
      Assert(kappa > 0, ExcInternalError());
    }

    ~Material_Compressible_Neo_Hook_Three_Field()
    {}

    // We update the material model with various deformation dependent data
    // based on $F$ and the pressure $\widetilde{p}$ and dilatation
    // $\widetilde{J}$, and at the end of the function include a physical
    // check for internal consistency:
    void update_material_data(const Tensor<2, dim> &F,
                              const double p_tilde_in,
                              const double J_tilde_in)
    {
      det_F = determinant(F);
      b_bar = std::pow(det_F, -2.0 / 3.0) * symmetrize(F * transpose(F));
      p_tilde = p_tilde_in;
      J_tilde = J_tilde_in;

      Assert(det_F > 0, ExcInternalError());
    }

    // The second function determines the Kirchhoff stress $\boldsymbol{\tau}
    // = \boldsymbol{\tau}_{\textrm{iso}} + \boldsymbol{\tau}_{\textrm{vol}}$
    SymmetricTensor<2, dim> get_tau()
    {
      return get_tau_iso() + get_tau_vol();
    }

    // The fourth-order elasticity tensor in the spatial setting
    // $\mathfrak{c}$ is calculated from the SEF $\Psi$ as $ J
    // \mathfrak{c}_{ijkl} = F_{iA} F_{jB} \mathfrak{C}_{ABCD} F_{kC} F_{lD}$
    // where $ \mathfrak{C} = 4 \frac{\partial^2 \Psi(\mathbf{C})}{\partial
    // \mathbf{C} \partial \mathbf{C}}$
    SymmetricTensor<4, dim> get_Jc() const
    {
      return get_Jc_vol() + get_Jc_iso();
    }

    // Derivative of the volumetric free energy with respect to
    // $\widetilde{J}$ return $\frac{\partial
    // \Psi_{\text{vol}}(\widetilde{J})}{\partial \widetilde{J}}$
    double get_dPsi_vol_dJ() const
    {
      return (kappa / 2.0) * (J_tilde - 1.0 / J_tilde);
    }

    // Second derivative of the volumetric free energy wrt $\widetilde{J}$. We
    // need the following computation explicitly in the tangent so we make it
    // public.  We calculate $\frac{\partial^2
    // \Psi_{\textrm{vol}}(\widetilde{J})}{\partial \widetilde{J} \partial
    // \widetilde{J}}$
    double get_d2Psi_vol_dJ2() const
    {
      return ( (kappa / 2.0) * (1.0 + 1.0 / (J_tilde * J_tilde)));
    }

    // The next few functions return various data that we choose to store with
    // the material:
    double get_det_F() const
    {
      return det_F;
    }

    double get_p_tilde() const
    {
      return p_tilde;
    }

    double get_J_tilde() const
    {
      return J_tilde;
    }

  protected:
    // Define constitutive model parameters $\kappa$ (bulk modulus) and the
    // neo-Hookean model parameter $c_1$:
    const double kappa;
    const double c_1;

    // Model specific data that is convenient to store with the material:
    double det_F;
    double p_tilde;
    double J_tilde;
    SymmetricTensor<2, dim> b_bar;

    // The following functions are used internally in determining the result
    // of some of the public functions above. The first one determines the
    // volumetric Kirchhoff stress $\boldsymbol{\tau}_{\textrm{vol}}$:
    SymmetricTensor<2, dim> get_tau_vol() const
    {
      return p_tilde * det_F * StandardTensors<dim>::I;
    }

    // Next, determine the isochoric Kirchhoff stress
    // $\boldsymbol{\tau}_{\textrm{iso}} =
    // \mathcal{P}:\overline{\boldsymbol{\tau}}$:
    SymmetricTensor<2, dim> get_tau_iso() const
    {
      return StandardTensors<dim>::dev_P * get_tau_bar();
    }

    // Then, determine the fictitious Kirchhoff stress
    // $\overline{\boldsymbol{\tau}}$:
    SymmetricTensor<2, dim> get_tau_bar() const
    {
      return 2.0 * c_1 * b_bar;
    }

    // Calculate the volumetric part of the tangent $J
    // \mathfrak{c}_\textrm{vol}$:
    SymmetricTensor<4, dim> get_Jc_vol() const
    {

      return p_tilde * det_F
             * ( StandardTensors<dim>::IxI
                 - (2.0 * StandardTensors<dim>::II) );
    }

    // Calculate the isochoric part of the tangent $J
    // \mathfrak{c}_\textrm{iso}$:
    SymmetricTensor<4, dim> get_Jc_iso() const
    {
      const SymmetricTensor<2, dim> tau_bar = get_tau_bar();
      const SymmetricTensor<2, dim> tau_iso = get_tau_iso();
      const SymmetricTensor<4, dim> tau_iso_x_I
        = outer_product(tau_iso,
                        StandardTensors<dim>::I);
      const SymmetricTensor<4, dim> I_x_tau_iso
        = outer_product(StandardTensors<dim>::I,
                        tau_iso);
      const SymmetricTensor<4, dim> c_bar = get_c_bar();

      return (2.0 / 3.0) * trace(tau_bar)
             * StandardTensors<dim>::dev_P
             - (2.0 / 3.0) * (tau_iso_x_I + I_x_tau_iso)
             + StandardTensors<dim>::dev_P * c_bar
             * StandardTensors<dim>::dev_P;
    }

    // Calculate the fictitious elasticity tensor $\overline{\mathfrak{c}}$.
    // For the material model chosen this is simply zero:
    SymmetricTensor<4, dim> get_c_bar() const
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
  template <int dim>
  class PointHistory
  {
  public:
    PointHistory()
      :
      material(NULL),
      F_inv(StandardTensors<dim>::I),
      tau(SymmetricTensor<2, dim>()),
      d2Psi_vol_dJ2(0.0),
      dPsi_vol_dJ(0.0),
      Jc(SymmetricTensor<4, dim>())
    {}

    virtual ~PointHistory()
    {
      delete material;
      material = NULL;
    }

    // The first function is used to create a material object and to
    // initialize all tensors correctly: The second one updates the stored
    // values and stresses based on the current deformation measure
    // $\textrm{Grad}\mathbf{u}_{\textrm{n}}$, pressure $\widetilde{p}$ and
    // dilation $\widetilde{J}$ field values.
    void setup_lqp (const Parameters::AllParameters &parameters)
    {
      material = new Material_Compressible_Neo_Hook_Three_Field<dim>(parameters.mu,
          parameters.nu);
      update_values(Tensor<2, dim>(), 0.0, 1.0);
    }

    // To this end, we calculate the deformation gradient $\mathbf{F}$ from
    // the displacement gradient $\textrm{Grad}\ \mathbf{u}$, i.e.
    // $\mathbf{F}(\mathbf{u}) = \mathbf{I} + \textrm{Grad}\ \mathbf{u}$ and
    // then let the material model associated with this quadrature point
    // update itself. When computing the deformation gradient, we have to take
    // care with which data types we compare the sum $\mathbf{I} +
    // \textrm{Grad}\ \mathbf{u}$: Since $I$ has data type SymmetricTensor,
    // just writing <code>I + Grad_u_n</code> would convert the second
    // argument to a symmetric tensor, perform the sum, and then cast the
    // result to a Tensor (i.e., the type of a possibly nonsymmetric
    // tensor). However, since <code>Grad_u_n</code> is nonsymmetric in
    // general, the conversion to SymmetricTensor will fail. We can avoid this
    // back and forth by converting $I$ to Tensor first, and then performing
    // the addition as between nonsymmetric tensors:
    void update_values (const Tensor<2, dim> &Grad_u_n,
                        const double p_tilde,
                        const double J_tilde)
    {
      const Tensor<2, dim> F
        = (Tensor<2, dim>(StandardTensors<dim>::I) +
           Grad_u_n);
      material->update_material_data(F, p_tilde, J_tilde);

      // The material has been updated so we now calculate the Kirchhoff
      // stress $\mathbf{\tau}$, the tangent $J\mathfrak{c}$ and the first and
      // second derivatives of the volumetric free energy.
      //
      // We also store the inverse of the deformation gradient since we
      // frequently use it:
      F_inv = invert(F);
      tau = material->get_tau();
      Jc = material->get_Jc();
      dPsi_vol_dJ = material->get_dPsi_vol_dJ();
      d2Psi_vol_dJ2 = material->get_d2Psi_vol_dJ2();

    }

    // We offer an interface to retrieve certain data.  Here are the kinematic
    // variables:
    double get_J_tilde() const
    {
      return material->get_J_tilde();
    }

    double get_det_F() const
    {
      return material->get_det_F();
    }

    const Tensor<2, dim> &get_F_inv() const
    {
      return F_inv;
    }

    // ...and the kinetic variables.  These are used in the material and
    // global tangent matrix and residual assembly operations:
    double get_p_tilde() const
    {
      return material->get_p_tilde();
    }

    const SymmetricTensor<2, dim> &get_tau() const
    {
      return tau;
    }

    double get_dPsi_vol_dJ() const
    {
      return dPsi_vol_dJ;
    }

    double get_d2Psi_vol_dJ2() const
    {
      return d2Psi_vol_dJ2;
    }

    // And finally the tangent:
    const SymmetricTensor<4, dim> &get_Jc() const
    {
      return Jc;
    }

    // In terms of member functions, this class stores for the quadrature
    // point it represents a copy of a material type in case different
    // materials are used in different regions of the domain, as well as the
    // inverse of the deformation gradient...
  private:
    Material_Compressible_Neo_Hook_Three_Field<dim> *material;

    Tensor<2, dim> F_inv;

    // ... and stress-type variables along with the tangent $J\mathfrak{c}$:
    SymmetricTensor<2, dim> tau;
    double                  d2Psi_vol_dJ2;
    double                  dPsi_vol_dJ;

    SymmetricTensor<4, dim> Jc;
  };


// @sect3{Quasi-static quasi-incompressible finite-strain solid}

// The Solid class is the central class in that it represents the problem at
// hand. It follows the usual scheme in that all it really has is a
// constructor, destructor and a <code>run()</code> function that dispatches
// all the work to private functions of this class:
  template <int dim>
  class Solid
  {
  public:
    Solid(const std::string &input_file);

    virtual
    ~Solid();

    void
    run();

  private:

    // In the private section of this class, we first forward declare a number
    // of objects that are used in parallelizing work using the WorkStream
    // object (see the @ref threads module for more information on this).
    //
    // We declare such structures for the computation of tangent (stiffness)
    // matrix, right hand side, static condensation, and for updating
    // quadrature points:
    struct PerTaskData_K;
    struct ScratchData_K;

    struct PerTaskData_RHS;
    struct ScratchData_RHS;

    struct PerTaskData_SC;
    struct ScratchData_SC;

    struct PerTaskData_UQPH;
    struct ScratchData_UQPH;

    // We start the collection of member functions with one that builds the
    // grid:
    void
    make_grid();

    // Set up the finite element system to be solved:
    void
    system_setup();

    void
    determine_component_extractors();

    // Several functions to assemble the system and right hand side matrices
    // using multithreading. Each of them comes as a wrapper function, one
    // that is executed to do the work in the WorkStream model on one cell,
    // and one that copies the work done on this one cell into the global
    // object that represents it:
    void
    assemble_system_tangent();

    void
    assemble_system_tangent_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                     ScratchData_K &scratch,
                                     PerTaskData_K &data);

    void
    copy_local_to_global_K(const PerTaskData_K &data);

    void
    assemble_system_rhs();

    void
    assemble_system_rhs_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                 ScratchData_RHS &scratch,
                                 PerTaskData_RHS &data);

    void
    copy_local_to_global_rhs(const PerTaskData_RHS &data);

    void
    assemble_sc();

    void
    assemble_sc_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                         ScratchData_SC &scratch,
                         PerTaskData_SC &data);

    void
    copy_local_to_global_sc(const PerTaskData_SC &data);

    // Apply Dirichlet boundary conditions on the displacement field
    void
    make_constraints(const int &it_nr);

    // Create and update the quadrature points. Here, no data needs to be
    // copied into a global object, so the copy_local_to_global function is
    // empty:
    void
    setup_qph();

    void
    update_qph_incremental(const BlockVector<double> &solution_delta);

    void
    update_qph_incremental_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                    ScratchData_UQPH &scratch,
                                    PerTaskData_UQPH &data);

    void
    copy_local_to_global_UQPH(const PerTaskData_UQPH &/*data*/)
    {}

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
    Parameters::AllParameters        parameters;

    // ...the volume of the reference and current configurations...
    double                           vol_reference;
    double                           vol_current;

    // ...and description of the geometry on which the problem is solved:
    Triangulation<dim>               triangulation;

    // Also, keep track of the current time and the time spent evaluating
    // certain functions
    Time                             time;
    TimerOutput                      timer;

    // A storage object for quadrature point information.  See step-18 for
    // more on this:
    std::vector<PointHistory<dim> >  quadrature_point_history;

    // A description of the finite-element system including the displacement
    // polynomial degree, the degree-of-freedom handler, number of DoFs per
    // cell and the extractor objects used to retrieve information from the
    // solution vectors:
    const unsigned int               degree;
    const FESystem<dim>              fe;
    DoFHandler<dim>                  dof_handler_ref;
    const unsigned int               dofs_per_cell;
    const FEValuesExtractors::Vector u_fe;
    const FEValuesExtractors::Scalar p_fe;
    const FEValuesExtractors::Scalar J_fe;

    // Description of how the block-system is arranged. There are 3 blocks,
    // the first contains a vector DOF $\mathbf{u}$ while the other two
    // describe scalar DOFs, $\widetilde{p}$ and $\widetilde{J}$.
    static const unsigned int        n_blocks = 3;
    static const unsigned int        n_components = dim + 2;
    static const unsigned int        first_u_component = 0;
    static const unsigned int        p_component = dim;
    static const unsigned int        J_component = dim + 1;

    enum
    {
      u_dof = 0,
      p_dof = 1,
      J_dof = 2
    };

    std::vector<types::global_dof_index>  dofs_per_block;
    std::vector<types::global_dof_index>        element_indices_u;
    std::vector<types::global_dof_index>        element_indices_p;
    std::vector<types::global_dof_index>        element_indices_J;

    // Rules for Gauss-quadrature on both the cell and faces. The number of
    // quadrature points on both cells and faces is recorded.
    const QGauss<dim>                qf_cell;
    const QGauss<dim - 1>            qf_face;
    const unsigned int               n_q_points;
    const unsigned int               n_q_points_f;

    // Objects that store the converged solution and right-hand side vectors,
    // as well as the tangent matrix. There is a ConstraintMatrix object used
    // to keep track of constraints.  We make use of a sparsity pattern
    // designed for a block system.
    ConstraintMatrix                 constraints;
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
        norm(1.0), u(1.0), p(1.0), J(1.0)
      {}

      void reset()
      {
        norm = 1.0;
        u = 1.0;
        p = 1.0;
        J = 1.0;
      }
      void normalise(const Errors &rhs)
      {
        if (rhs.norm != 0.0)
          norm /= rhs.norm;
        if (rhs.u != 0.0)
          u /= rhs.u;
        if (rhs.p != 0.0)
          p /= rhs.p;
        if (rhs.J != 0.0)
          J /= rhs.J;
      }

      double norm, u, p, J;
    };

    Errors error_residual, error_residual_0, error_residual_norm, error_update,
           error_update_0, error_update_norm;

    // Methods to calculate error measures
    void
    get_error_residual(Errors &error_residual);

    void
    get_error_update(const BlockVector<double> &newton_update,
                     Errors &error_update);

    std::pair<double, double>
    get_error_dilation();

    // Print information to screen in a pleasing way...
    static
    void
    print_conv_header();

    void
    print_conv_footer();
  };

// @sect3{Implementation of the <code>Solid</code> class}

// @sect4{Public interface}

// We initialise the Solid class using data extracted from the parameter file.
  template <int dim>
  Solid<dim>::Solid(const std::string &input_file)
    :
    parameters(input_file),
    triangulation(Triangulation<dim>::maximum_smoothing),
    time(parameters.end_time, parameters.delta_t),
    timer(std::cout,
          TimerOutput::summary,
          TimerOutput::wall_times),
    degree(parameters.poly_degree),
    // The Finite Element System is composed of dim continuous displacement
    // DOFs, and discontinuous pressure and dilatation DOFs. In an attempt to
    // satisfy the Babuska-Brezzi or LBB stability conditions (see Hughes
    // (2000)), we setup a $Q_n \times DGPM_{n-1} \times DGPM_{n-1}$
    // system. $Q_2 \times DGPM_1 \times DGPM_1$ elements satisfy this
    // condition, while $Q_1 \times DGPM_0 \times DGPM_0$ elements do
    // not. However, it has been shown that the latter demonstrate good
    // convergence characteristics nonetheless.
    fe(FE_Q<dim>(parameters.poly_degree), dim, // displacement
       FE_DGPMonomial<dim>(parameters.poly_degree - 1), 1, // pressure
       FE_DGPMonomial<dim>(parameters.poly_degree - 1), 1), // dilatation
    dof_handler_ref(triangulation),
    dofs_per_cell (fe.dofs_per_cell),
    u_fe(first_u_component),
    p_fe(p_component),
    J_fe(J_component),
    dofs_per_block(n_blocks),
    qf_cell(parameters.quad_order),
    qf_face(parameters.quad_order),
    n_q_points (qf_cell.size()),
    n_q_points_f (qf_face.size())
  {
    determine_component_extractors();
  }

// The class destructor simply clears the data held by the DOFHandler
  template <int dim>
  Solid<dim>::~Solid()
  {
    dof_handler_ref.clear();
  }


// In solving the quasi-static problem, the time becomes a loading parameter,
// i.e. we increasing the loading linearly with time, making the two concepts
// interchangeable. We choose to increment time linearly using a constant time
// step size.
//
// We start the function with preprocessing, setting the initial dilatation
// values, and then output the initial grid before starting the simulation
//  proper with the first time (and loading)
// increment.
//
// Care must be taken (or at least some thought given) when imposing the
// constraint $\widetilde{J}=1$ on the initial solution field. The constraint
// corresponds to the determinant of the deformation gradient in the undeformed
// configuration, which is the identity tensor.
// We use FE_DGPMonomial bases to interpolate the dilatation field, thus we can't
// simply set the corresponding dof to unity as they correspond to the
// monomial coefficients. Thus we use the VectorTools::project function to do
// the work for us. The VectorTools::project function requires an argument
// indicating the hanging node constraints. We have none in this program
// So we have to create a constraint object. In its original state, constraint
// objects are unsorted, and have to be sorted (using the ConstraintMatrix::close function)
// before they can be used. Have a look at step-21 for more information.
// We only need to enforce the initial condition on the dilatation.
// In order to do this, we make use of a ComponentSelectFunction which acts
// as a mask and sets the J_component of n_components to 1. This is exactly what
// we want. Have a look at its usage in step-20 for more information.
  template <int dim>
  void Solid<dim>::run()
  {
    make_grid();
    system_setup();
    {
      ConstraintMatrix constraints;
      constraints.close();

      const ComponentSelectFunction<dim>
      J_mask (J_component, n_components);

      VectorTools::project (dof_handler_ref,
                            constraints,
                            QGauss<dim>(degree+2),
                            J_mask,
                            solution_n);
    }
    output_results();
    time.increment();

    // We then declare the incremental solution update $\varDelta
    // \mathbf{\Xi}:= \{\varDelta \mathbf{u},\varDelta \widetilde{p},
    // \varDelta \widetilde{J} \}$ and start the loop over the time domain.
    //
    // At the beginning, we reset the solution update for this time step...
    BlockVector<double> solution_delta(dofs_per_block);
    while (time.current() < time.end())
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
  }


// @sect3{Private interface}

// @sect4{Threading-building-blocks structures}

// The first group of private member functions is related to parallization.
// We use the Threading Building Blocks library (TBB) to perform as many
// computationally intensive distributed tasks as possible. In particular, we
// assemble the tangent matrix and right hand side vector, the static
// condensation contributions, and update data stored at the quadrature points
// using TBB. Our main tool for this is the WorkStream class (see the @ref
// threads module for more information).

// Firstly we deal with the tangent matrix assembly structures.  The
// PerTaskData object stores local contributions.
  template <int dim>
  struct Solid<dim>::PerTaskData_K
  {
    FullMatrix<double>        cell_matrix;
    std::vector<types::global_dof_index> local_dof_indices;

    PerTaskData_K(const unsigned int dofs_per_cell)
      :
      cell_matrix(dofs_per_cell, dofs_per_cell),
      local_dof_indices(dofs_per_cell)
    {}

    void reset()
    {
      cell_matrix = 0.0;
    }
  };


// On the other hand, the ScratchData object stores the larger objects such as
// the shape-function values array (<code>Nx</code>) and a shape function
// gradient and symmetric gradient vector which we will use during the
// assembly.
  template <int dim>
  struct Solid<dim>::ScratchData_K
  {
    FEValues<dim> fe_values_ref;

    std::vector<std::vector<double> >                   Nx;
    std::vector<std::vector<Tensor<2, dim> > >          grad_Nx;
    std::vector<std::vector<SymmetricTensor<2, dim> > > symm_grad_Nx;

    ScratchData_K(const FiniteElement<dim> &fe_cell,
                  const QGauss<dim> &qf_cell,
                  const UpdateFlags uf_cell)
      :
      fe_values_ref(fe_cell, qf_cell, uf_cell),
      Nx(qf_cell.size(),
         std::vector<double>(fe_cell.dofs_per_cell)),
      grad_Nx(qf_cell.size(),
              std::vector<Tensor<2, dim> >(fe_cell.dofs_per_cell)),
      symm_grad_Nx(qf_cell.size(),
                   std::vector<SymmetricTensor<2, dim> >
                   (fe_cell.dofs_per_cell))
    {}

    ScratchData_K(const ScratchData_K &rhs)
      :
      fe_values_ref(rhs.fe_values_ref.get_fe(),
                    rhs.fe_values_ref.get_quadrature(),
                    rhs.fe_values_ref.get_update_flags()),
      Nx(rhs.Nx),
      grad_Nx(rhs.grad_Nx),
      symm_grad_Nx(rhs.symm_grad_Nx)
    {}

    void reset()
    {
      const unsigned int n_q_points = Nx.size();
      const unsigned int n_dofs_per_cell = Nx[0].size();
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          Assert( Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
          Assert( grad_Nx[q_point].size() == n_dofs_per_cell,
                  ExcInternalError());
          Assert( symm_grad_Nx[q_point].size() == n_dofs_per_cell,
                  ExcInternalError());
          for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
            {
              Nx[q_point][k] = 0.0;
              grad_Nx[q_point][k] = 0.0;
              symm_grad_Nx[q_point][k] = 0.0;
            }
        }
    }

  };

// Next, the same approach is used for the right-hand side assembly.  The
// PerTaskData object again stores local contributions and the ScratchData
// object the shape function object and precomputed values vector:
  template <int dim>
  struct Solid<dim>::PerTaskData_RHS
  {
    Vector<double>            cell_rhs;
    std::vector<types::global_dof_index> local_dof_indices;

    PerTaskData_RHS(const unsigned int dofs_per_cell)
      :
      cell_rhs(dofs_per_cell),
      local_dof_indices(dofs_per_cell)
    {}

    void reset()
    {
      cell_rhs = 0.0;
    }
  };


  template <int dim>
  struct Solid<dim>::ScratchData_RHS
  {
    FEValues<dim>     fe_values_ref;
    FEFaceValues<dim> fe_face_values_ref;

    std::vector<std::vector<double> >                   Nx;
    std::vector<std::vector<SymmetricTensor<2, dim> > > symm_grad_Nx;

    ScratchData_RHS(const FiniteElement<dim> &fe_cell,
                    const QGauss<dim> &qf_cell, const UpdateFlags uf_cell,
                    const QGauss<dim - 1> & qf_face, const UpdateFlags uf_face)
      :
      fe_values_ref(fe_cell, qf_cell, uf_cell),
      fe_face_values_ref(fe_cell, qf_face, uf_face),
      Nx(qf_cell.size(),
         std::vector<double>(fe_cell.dofs_per_cell)),
      symm_grad_Nx(qf_cell.size(),
                   std::vector<SymmetricTensor<2, dim> >
                   (fe_cell.dofs_per_cell))
    {}

    ScratchData_RHS(const ScratchData_RHS &rhs)
      :
      fe_values_ref(rhs.fe_values_ref.get_fe(),
                    rhs.fe_values_ref.get_quadrature(),
                    rhs.fe_values_ref.get_update_flags()),
      fe_face_values_ref(rhs.fe_face_values_ref.get_fe(),
                         rhs.fe_face_values_ref.get_quadrature(),
                         rhs.fe_face_values_ref.get_update_flags()),
      Nx(rhs.Nx),
      symm_grad_Nx(rhs.symm_grad_Nx)
    {}

    void reset()
    {
      const unsigned int n_q_points      = Nx.size();
      const unsigned int n_dofs_per_cell = Nx[0].size();
      for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
        {
          Assert( Nx[q_point].size() == n_dofs_per_cell, ExcInternalError());
          Assert( symm_grad_Nx[q_point].size() == n_dofs_per_cell,
                  ExcInternalError());
          for (unsigned int k = 0; k < n_dofs_per_cell; ++k)
            {
              Nx[q_point][k] = 0.0;
              symm_grad_Nx[q_point][k] = 0.0;
            }
        }
    }

  };

// Then we define structures to assemble the statically condensed tangent
// matrix. Recall that we wish to solve for a displacement-based formulation.
// We do the condensation at the element level as the $\widetilde{p}$ and
// $\widetilde{J}$ fields are element-wise discontinuous.  As these operations
// are matrix-based, we need to setup a number of matrices to store the local
// contributions from a number of the tangent matrix sub-blocks.  We place
// these in the PerTaskData struct.
//
// We choose not to reset any data in the <code>reset()</code> function as the
// matrix extraction and replacement tools will take care of this
  template <int dim>
  struct Solid<dim>::PerTaskData_SC
  {
    FullMatrix<double>        cell_matrix;
    std::vector<types::global_dof_index> local_dof_indices;

    FullMatrix<double>        k_orig;
    FullMatrix<double>        k_pu;
    FullMatrix<double>        k_pJ;
    FullMatrix<double>        k_JJ;
    FullMatrix<double>        k_pJ_inv;
    FullMatrix<double>        k_bbar;
    FullMatrix<double>        A;
    FullMatrix<double>        B;
    FullMatrix<double>        C;

    PerTaskData_SC(const unsigned int dofs_per_cell,
                   const unsigned int n_u,
                   const unsigned int n_p,
                   const unsigned int n_J)
      :
      cell_matrix(dofs_per_cell, dofs_per_cell),
      local_dof_indices(dofs_per_cell),
      k_orig(dofs_per_cell, dofs_per_cell),
      k_pu(n_p, n_u),
      k_pJ(n_p, n_J),
      k_JJ(n_J, n_J),
      k_pJ_inv(n_p, n_J),
      k_bbar(n_u, n_u),
      A(n_J,n_u),
      B(n_J, n_u),
      C(n_p, n_u)
    {}

    void reset()
    {}
  };


// The ScratchData object for the operations we wish to perform here is empty
// since we need no temporary data, but it still needs to be defined for the
// current implementation of TBB in deal.II.  So we create a dummy struct for
// this purpose.
  template <int dim>
  struct Solid<dim>::ScratchData_SC
  {
    void reset()
    {}
  };


// And finally we define the structures to assist with updating the quadrature
// point information. Similar to the SC assembly process, we do not need the
// PerTaskData object (since there is nothing to store here) but must define
// one nonetheless. Note that this is because for the operation that we have
// here -- updating the data on quadrature points -- the operation is purely
// local: the things we do on every cell get consumed on every cell, without
// any global aggregation operation as is usually the case when using the
// WorkStream class. The fact that we still have to define a per-task data
// structure points to the fact that the WorkStream class may be ill-suited to
// this operation (we could, in principle simply create a new task using
// Threads::new_task for each cell) but there is not much harm done to doing
// it this way anyway.
// Furthermore, should there be different material models associated with a
// quadrature point, requiring varying levels of computational expense, then
// the method used here could be advantageous.
  template <int dim>
  struct Solid<dim>::PerTaskData_UQPH
  {
    void reset()
    {}
  };


// The ScratchData object will be used to store an alias for the solution
// vector so that we don't have to copy this large data structure. We then
// define a number of vectors to extract the solution values and gradients at
// the quadrature points.
  template <int dim>
  struct Solid<dim>::ScratchData_UQPH
  {
    const BlockVector<double>   &solution_total;

    std::vector<Tensor<2, dim> > solution_grads_u_total;
    std::vector<double>          solution_values_p_total;
    std::vector<double>          solution_values_J_total;

    FEValues<dim>                fe_values_ref;

    ScratchData_UQPH(const FiniteElement<dim> &fe_cell,
                     const QGauss<dim> &qf_cell,
                     const UpdateFlags uf_cell,
                     const BlockVector<double> &solution_total)
      :
      solution_total(solution_total),
      solution_grads_u_total(qf_cell.size()),
      solution_values_p_total(qf_cell.size()),
      solution_values_J_total(qf_cell.size()),
      fe_values_ref(fe_cell, qf_cell, uf_cell)
    {}

    ScratchData_UQPH(const ScratchData_UQPH &rhs)
      :
      solution_total(rhs.solution_total),
      solution_grads_u_total(rhs.solution_grads_u_total),
      solution_values_p_total(rhs.solution_values_p_total),
      solution_values_J_total(rhs.solution_values_J_total),
      fe_values_ref(rhs.fe_values_ref.get_fe(),
                    rhs.fe_values_ref.get_quadrature(),
                    rhs.fe_values_ref.get_update_flags())
    {}

    void reset()
    {
      const unsigned int n_q_points = solution_grads_u_total.size();
      for (unsigned int q = 0; q < n_q_points; ++q)
        {
          solution_grads_u_total[q] = 0.0;
          solution_values_p_total[q] = 0.0;
          solution_values_J_total[q] = 0.0;
        }
    }
  };


// @sect4{Solid::make_grid}

// On to the first of the private member functions. Here we create the
// triangulation of the domain, for which we choose the scaled cube with each
// face given a boundary ID number.  The grid must be refined at least once
// for the indentation problem.
//
// We then determine the volume of the reference configuration and print it
// for comparison:
  template <int dim>
  void Solid<dim>::make_grid()
  {
    GridGenerator::hyper_rectangle(triangulation,
                                   Point<dim>(0.0, 0.0, 0.0),
                                   Point<dim>(1.0, 1.0, 1.0),
                                   true);
    GridTools::scale(parameters.scale, triangulation);
    triangulation.refine_global(std::max (1U, parameters.global_refinement));

    vol_reference = GridTools::volume(triangulation);
    vol_current = vol_reference;
    std::cout << "Grid:\n\t Reference volume: " << vol_reference << std::endl;

    // Since we wish to apply a Neumann BC to a patch on the top surface, we
    // must find the cell faces in this part of the domain and mark them with
    // a distinct boundary ID number.  The faces we are looking for are on the
    // +y surface and will get boundary ID 6 (zero through five are already
    // used when creating the six faces of the cube domain):
    typename Triangulation<dim>::active_cell_iterator cell =
      triangulation.begin_active(), endc = triangulation.end();
    for (; cell != endc; ++cell)
      for (unsigned int face = 0;
           face < GeometryInfo<dim>::faces_per_cell; ++face)
        if (cell->face(face)->at_boundary() == true
            &&
            cell->face(face)->center()[2] == 1.0 * parameters.scale)
          if (cell->face(face)->center()[0] < 0.5 * parameters.scale
              &&
              cell->face(face)->center()[1] < 0.5 * parameters.scale)
            cell->face(face)->set_boundary_indicator(6);
  }


// @sect4{Solid::system_setup}

// Next we describe how the FE system is setup.  We first determine the number
// of components per block. Since the displacement is a vector component, the
// first dim components belong to it, while the next two describe scalar
// pressure and dilatation DOFs.
  template <int dim>
  void Solid<dim>::system_setup()
  {
    timer.enter_subsection("Setup system");

    std::vector<unsigned int> block_component(n_components, u_dof); // Displacement
    block_component[p_component] = p_dof; // Pressure
    block_component[J_component] = J_dof; // Dilatation

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
      const types::global_dof_index n_dofs_p = dofs_per_block[p_dof];
      const types::global_dof_index n_dofs_J = dofs_per_block[J_dof];

      BlockCompressedSimpleSparsityPattern csp(n_blocks, n_blocks);

      csp.block(u_dof, u_dof).reinit(n_dofs_u, n_dofs_u);
      csp.block(u_dof, p_dof).reinit(n_dofs_u, n_dofs_p);
      csp.block(u_dof, J_dof).reinit(n_dofs_u, n_dofs_J);

      csp.block(p_dof, u_dof).reinit(n_dofs_p, n_dofs_u);
      csp.block(p_dof, p_dof).reinit(n_dofs_p, n_dofs_p);
      csp.block(p_dof, J_dof).reinit(n_dofs_p, n_dofs_J);

      csp.block(J_dof, u_dof).reinit(n_dofs_J, n_dofs_u);
      csp.block(J_dof, p_dof).reinit(n_dofs_J, n_dofs_p);
      csp.block(J_dof, J_dof).reinit(n_dofs_J, n_dofs_J);
      csp.collect_sizes();

      // The global system matrix initially has the following structure
      // @f{align*}
      // \underbrace{\begin{bmatrix}
      //   \mathbf{\mathsf{K}}_{uu}  & \mathbf{\mathsf{K}}_{u\widetilde{p}} & \mathbf{0}
      //   \\ \mathbf{\mathsf{K}}_{\widetilde{p}u} & \mathbf{0} & \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}
      //   \\ \mathbf{0} & \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}} & \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
      // \end{bmatrix}}_{\mathbf{\mathsf{K}}(\mathbf{\Xi}_{\textrm{i}})}
      //      \underbrace{\begin{bmatrix}
      //          d \mathbf{\mathsf{u}}
      //      \\  d \widetilde{\mathbf{\mathsf{p}}}
      //      \\  d \widetilde{\mathbf{\mathsf{J}}}
      //      \end{bmatrix}}_{d \mathbf{\Xi}}
      // =
      // \underbrace{\begin{bmatrix}
      //  \mathbf{\mathsf{F}}_{u}(\mathbf{u}_{\textrm{i}})
      //  \\ \mathbf{\mathsf{F}}_{\widetilde{p}}(\widetilde{p}_{\textrm{i}})
      //  \\ \mathbf{\mathsf{F}}_{\widetilde{J}}(\widetilde{J}_{\textrm{i}})
      //\end{bmatrix}}_{ \mathbf{\mathsf{F}}(\mathbf{\Xi}_{\textrm{i}}) } \, .
      // @f}
      // We optimise the sparsity pattern to reflect this structure
      // and prevent unnecessary data creation for the right-diagonal
      // block components.
      Table<2, DoFTools::Coupling> coupling(n_components, n_components);
      for (unsigned int ii = 0; ii < n_components; ++ii)
        for (unsigned int jj = 0; jj < n_components; ++jj)
          if (((ii < p_component) && (jj == J_component))
              || ((ii == J_component) && (jj < p_component))
              || ((ii == p_component) && (jj == p_component)))
            coupling[ii][jj] = DoFTools::none;
          else
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


// @sect4{Solid::determine_component_extractors}
// Next we compute some information from the FE system that describes which local
// element DOFs are attached to which block component.  This is used later to
// extract sub-blocks from the global matrix.
//
// In essence, all we need is for the FESystem object to indicate to which
// block component a DOF on the reference cell is attached to.  Currently, the
// interpolation fields are setup such that 0 indicates a displacement DOF, 1
// a pressure DOF and 2 a dilatation DOF.
  template <int dim>
  void
  Solid<dim>::determine_component_extractors()
  {
    element_indices_u.clear();
    element_indices_p.clear();
    element_indices_J.clear();

    for (unsigned int k = 0; k < fe.dofs_per_cell; ++k)
      {
        const unsigned int k_group = fe.system_to_base_index(k).first.first;
        if (k_group == u_dof)
          element_indices_u.push_back(k);
        else if (k_group == p_dof)
          element_indices_p.push_back(k);
        else if (k_group == J_dof)
          element_indices_J.push_back(k);
        else
          {
            Assert(k_group <= J_dof, ExcInternalError());
          }
      }
  }

// @sect4{Solid::setup_qph}
// The method used to store quadrature information is already described in
// step-18. Here we implement a similar setup for a SMP machine.
//
// Firstly the actual QPH data objects are created. This must be done only
// once the grid is refined to its finest level.
  template <int dim>
  void Solid<dim>::setup_qph()
  {
    std::cout << "    Setting up quadrature point data..." << std::endl;

    {
      triangulation.clear_user_data();
      {
        std::vector<PointHistory<dim> > tmp;
        tmp.swap(quadrature_point_history);
      }

      quadrature_point_history
      .resize(triangulation.n_active_cells() * n_q_points);

      unsigned int history_index = 0;
      for (typename Triangulation<dim>::active_cell_iterator cell =
             triangulation.begin_active(); cell != triangulation.end();
           ++cell)
        {
          cell->set_user_pointer(&quadrature_point_history[history_index]);
          history_index += n_q_points;
        }

      Assert(history_index == quadrature_point_history.size(),
             ExcInternalError());
    }

    // Next we setup the initial quadrature
    // point data:
    for (typename Triangulation<dim>::active_cell_iterator cell =
           triangulation.begin_active(); cell != triangulation.end(); ++cell)
      {
        PointHistory<dim> *lqph =
          reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

        Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
        Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          lqph[q_point].setup_lqp(parameters);
      }
  }

// @sect4{Solid::update_qph_incremental}
// As the update of QP information occurs frequently and involves a number of
// expensive operations, we define a multithreaded approach to distributing
// the task across a number of CPU cores.
//
// To start this, we first we need to obtain the total solution as it stands
// at this Newton increment and then create the initial copy of the scratch and
// copy data objects:
  template <int dim>
  void Solid<dim>::update_qph_incremental(const BlockVector<double> &solution_delta)
  {
    timer.enter_subsection("Update QPH data");
    std::cout << " UQPH " << std::flush;

    const BlockVector<double> solution_total(get_total_solution(solution_delta));

    const UpdateFlags uf_UQPH(update_values | update_gradients);
    PerTaskData_UQPH per_task_data_UQPH;
    ScratchData_UQPH scratch_data_UQPH(fe, qf_cell, uf_UQPH, solution_total);

    // We then pass them and the one-cell update function to the WorkStream to
    // be processed:
    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    *this,
                    &Solid::update_qph_incremental_one_cell,
                    &Solid::copy_local_to_global_UQPH,
                    scratch_data_UQPH,
                    per_task_data_UQPH);

    timer.leave_subsection();
  }


// Now we describe how we extract data from the solution vector and pass it
// along to each QP storage object for processing.
  template <int dim>
  void
  Solid<dim>::update_qph_incremental_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                              ScratchData_UQPH &scratch,
                                              PerTaskData_UQPH &/*data*/)
  {
    PointHistory<dim> *lqph =
      reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
    Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

    Assert(scratch.solution_grads_u_total.size() == n_q_points,
           ExcInternalError());
    Assert(scratch.solution_values_p_total.size() == n_q_points,
           ExcInternalError());
    Assert(scratch.solution_values_J_total.size() == n_q_points,
           ExcInternalError());

    scratch.reset();

    // We first need to find the values and gradients at quadrature points
    // inside the current cell and then we update each local QP using the
    // displacement gradient and total pressure and dilatation solution
    // values:
    scratch.fe_values_ref.reinit(cell);
    scratch.fe_values_ref[u_fe].get_function_gradients(scratch.solution_total,
                                                       scratch.solution_grads_u_total);
    scratch.fe_values_ref[p_fe].get_function_values(scratch.solution_total,
                                                    scratch.solution_values_p_total);
    scratch.fe_values_ref[J_fe].get_function_values(scratch.solution_total,
                                                    scratch.solution_values_J_total);

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      lqph[q_point].update_values(scratch.solution_grads_u_total[q_point],
                                  scratch.solution_values_p_total[q_point],
                                  scratch.solution_values_J_total[q_point]);
  }


// @sect4{Solid::solve_nonlinear_timestep}

// The next function is the driver method for the Newton-Raphson scheme. At
// its top we create a new vector to store the current Newton update step,
// reset the error storage objects and print solver header.
  template <int dim>
  void
  Solid<dim>::solve_nonlinear_timestep(BlockVector<double> &solution_delta)
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

        tangent_matrix = 0.0;
        system_rhs = 0.0;

        assemble_system_rhs();
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

        // If we have decided that we want to continue with the iteration, we
        // assemble the tangent, make and impose the Dirichlet constraints,
        // and do the solve of the linearised system:
        assemble_system_tangent();
        make_constraints(newton_iteration);
        constraints.condense(tangent_matrix, system_rhs);

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
        update_qph_incremental(solution_delta);

        std::cout << " | " << std::fixed << std::setprecision(3) << std::setw(7)
                  << std::scientific << lin_solver_output.first << "  "
                  << lin_solver_output.second << "  " << error_residual_norm.norm
                  << "  " << error_residual_norm.u << "  "
                  << error_residual_norm.p << "  " << error_residual_norm.J
                  << "  " << error_update_norm.norm << "  " << error_update_norm.u
                  << "  " << error_update_norm.p << "  " << error_update_norm.J
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


// @sect4{Solid::print_conv_header and Solid::print_conv_footer}

// This program prints out data in a nice table that is updated
// on a per-iteration basis. The next two functions set up the table
// header and footer:
  template <int dim>
  void Solid<dim>::print_conv_header()
  {
    static const unsigned int l_width = 155;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;

    std::cout << "                 SOLVER STEP                  "
              << " |  LIN_IT   LIN_RES    RES_NORM    "
              << " RES_U     RES_P      RES_J     NU_NORM     "
              << " NU_U       NU_P       NU_J " << std::endl;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;
  }



  template <int dim>
  void Solid<dim>::print_conv_footer()
  {
    static const unsigned int l_width = 155;

    for (unsigned int i = 0; i < l_width; ++i)
      std::cout << "_";
    std::cout << std::endl;

    const std::pair <double,double> error_dil = get_error_dilation();

    std::cout << "Relative errors:" << std::endl
              << "Displacement:\t" << error_update.u / error_update_0.u << std::endl
              << "Force: \t\t" << error_residual.u / error_residual_0.u << std::endl
              << "Dilatation:\t" << error_dil.first << std::endl
              << "v / V_0:\t" << vol_current << " / " << vol_reference
              << " = " << error_dil.second << std::endl;
  }


// @sect4{Solid::get_error_dilation}

// Calculate how well the dilatation $\widetilde{J}$ agrees with $J :=
// \textrm{det}\ \mathbf{F}$ from the $L^2$ error $ \bigl[ \int_{\Omega_0} {[ J
// - \widetilde{J}]}^{2}\textrm{d}V \bigr]^{1/2}$.
// We also return the ratio of the current volume of the
// domain to the reference volume. This is of interest for incompressible
// media where we want to check how well the isochoric constraint has been
// enforced.
  template <int dim>
  std::pair<double, double>
  Solid<dim>::get_error_dilation()
  {
    double dil_L2_error = 0.0;
    vol_current = 0.0;

    FEValues<dim> fe_values_ref(fe, qf_cell, update_JxW_values);

    for (typename Triangulation<dim>::active_cell_iterator
         cell = triangulation.begin_active();
         cell != triangulation.end(); ++cell)
      {
        fe_values_ref.reinit(cell);

        PointHistory<dim> *lqph =
          reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

        Assert(lqph >= &quadrature_point_history.front(), ExcInternalError());
        Assert(lqph <= &quadrature_point_history.back(), ExcInternalError());

        for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
          {
            const double det_F_qp = lqph[q_point].get_det_F();
            const double J_tilde_qp = lqph[q_point].get_J_tilde();
            const double the_error_qp_squared = std::pow((det_F_qp - J_tilde_qp),
                                                         2);
            const double JxW = fe_values_ref.JxW(q_point);

            dil_L2_error += the_error_qp_squared * JxW;
            vol_current += det_F_qp * JxW;
          }
        Assert(vol_current > 0, ExcInternalError());
      }

    std::pair<double, double> error_dil;
    error_dil.first = std::sqrt(dil_L2_error);
    error_dil.second = vol_current / vol_reference;

    return error_dil;
  }


// @sect4{Solid::get_error_residual}

// Determine the true residual error for the problem.  That is, determine the
// error in the residual for the unconstrained degrees of freedom.  Note that to
// do so, we need to ignore constrained DOFs by setting the residual in these
// vector components to zero.
  template <int dim>
  void Solid<dim>::get_error_residual(Errors &error_residual)
  {
    BlockVector<double> error_res(dofs_per_block);

    for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
      if (!constraints.is_constrained(i))
        error_res(i) = system_rhs(i);

    error_residual.norm = error_res.l2_norm();
    error_residual.u = error_res.block(u_dof).l2_norm();
    error_residual.p = error_res.block(p_dof).l2_norm();
    error_residual.J = error_res.block(J_dof).l2_norm();
  }


// @sect4{Solid::get_error_udpate}

// Determine the true Newton update error for the problem
  template <int dim>
  void Solid<dim>::get_error_update(const BlockVector<double> &newton_update,
                                    Errors &error_update)
  {
    BlockVector<double> error_ud(dofs_per_block);
    for (unsigned int i = 0; i < dof_handler_ref.n_dofs(); ++i)
      if (!constraints.is_constrained(i))
        error_ud(i) = newton_update(i);

    error_update.norm = error_ud.l2_norm();
    error_update.u = error_ud.block(u_dof).l2_norm();
    error_update.p = error_ud.block(p_dof).l2_norm();
    error_update.J = error_ud.block(J_dof).l2_norm();
  }



// @sect4{Solid::get_total_solution}

// This function provides the total solution, which is valid at any Newton step.
// This is required as, to reduce computational error, the total solution is
// only updated at the end of the timestep.
  template <int dim>
  BlockVector<double>
  Solid<dim>::get_total_solution(const BlockVector<double> &solution_delta) const
  {
    BlockVector<double> solution_total(solution_n);
    solution_total += solution_delta;
    return solution_total;
  }


// @sect4{Solid::assemble_system_tangent}

// Since we use TBB for assembly, we simply setup a copy of the
// data structures required for the process and pass them, along
// with the memory addresses of the assembly functions to the
// WorkStream object for processing. Note that we must ensure that
// the matrix is reset before any assembly operations can occur.
  template <int dim>
  void Solid<dim>::assemble_system_tangent()
  {
    timer.enter_subsection("Assemble tangent matrix");
    std::cout << " ASM_K " << std::flush;

    tangent_matrix = 0.0;

    const UpdateFlags uf_cell(update_values    |
                              update_gradients |
                              update_JxW_values);

    PerTaskData_K per_task_data(dofs_per_cell);
    ScratchData_K scratch_data(fe, qf_cell, uf_cell);

    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    *this,
                    &Solid::assemble_system_tangent_one_cell,
                    &Solid::copy_local_to_global_K,
                    scratch_data,
                    per_task_data);

    timer.leave_subsection();
  }

// This function adds the local contribution to the system matrix.
// Note that we choose not to use the constraint matrix to do the
// job for us because the tangent matrix and residual processes have
// been split up into two separate functions.
  template <int dim>
  void Solid<dim>::copy_local_to_global_K(const PerTaskData_K &data)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        tangent_matrix.add(data.local_dof_indices[i],
                           data.local_dof_indices[j],
                           data.cell_matrix(i, j));
  }

// Of course, we still have to define how we assemble the tangent matrix
// contribution for a single cell.  We first need to reset and initialise some
// of the scratch data structures and retrieve some basic information
// regarding the DOF numbering on this cell.  We can precalculate the cell
// shape function values and gradients. Note that the shape function gradients
// are defined with regard to the current configuration.  That is
// $\textrm{grad}\ \boldsymbol{\varphi} = \textrm{Grad}\ \boldsymbol{\varphi}
// \ \mathbf{F}^{-1}$.
  template <int dim>
  void
  Solid<dim>::assemble_system_tangent_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                               ScratchData_K &scratch,
                                               PerTaskData_K &data)
  {
    data.reset();
    scratch.reset();
    scratch.fe_values_ref.reinit(cell);
    cell->get_dof_indices(data.local_dof_indices);
    PointHistory<dim> *lqph =
      reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        const Tensor<2, dim> F_inv = lqph[q_point].get_F_inv();
        for (unsigned int k = 0; k < dofs_per_cell; ++k)
          {
            const unsigned int k_group = fe.system_to_base_index(k).first.first;

            if (k_group == u_dof)
              {
                scratch.grad_Nx[q_point][k] = scratch.fe_values_ref[u_fe].gradient(k, q_point)
                                              * F_inv;
                scratch.symm_grad_Nx[q_point][k] = symmetrize(scratch.grad_Nx[q_point][k]);
              }
            else if (k_group == p_dof)
              scratch.Nx[q_point][k] = scratch.fe_values_ref[p_fe].value(k,
                                                                         q_point);
            else if (k_group == J_dof)
              scratch.Nx[q_point][k] = scratch.fe_values_ref[J_fe].value(k,
                                                                         q_point);
            else
              Assert(k_group <= J_dof, ExcInternalError());
          }
      }

    // Now we build the local cell stiffness matrix. Since the global and
    // local system matrices are symmetric, we can exploit this property by
    // building only the lower half of the local matrix and copying the values
    // to the upper half.  So we only assemble half of the
    // $\mathsf{\mathbf{k}}_{uu}$, $\mathsf{\mathbf{k}}_{\widetilde{p}
    // \widetilde{p}} = \mathbf{0}$, $\mathsf{\mathbf{k}}_{\widetilde{J}
    // \widetilde{J}}$ blocks, while the whole
    // $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$,
    // $\mathsf{\mathbf{k}}_{\mathbf{u} \widetilde{J}} = \mathbf{0}$,
    // $\mathsf{\mathbf{k}}_{\mathbf{u} \widetilde{p}}$ blocks are built.
    //
    // In doing so, we first extract some configuration dependent variables
    // from our QPH history objects for the current quadrature point.
    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        const Tensor<2, dim> tau         = lqph[q_point].get_tau();
        const SymmetricTensor<4, dim> Jc = lqph[q_point].get_Jc();
        const double d2Psi_vol_dJ2       = lqph[q_point].get_d2Psi_vol_dJ2();
        const double det_F               = lqph[q_point].get_det_F();

        // Next we define some aliases to make the assembly process easier to
        // follow
        const std::vector<double>
        &N = scratch.Nx[q_point];
        const std::vector<SymmetricTensor<2, dim> >
        &symm_grad_Nx = scratch.symm_grad_Nx[q_point];
        const std::vector<Tensor<2, dim> >
        &grad_Nx = scratch.grad_Nx[q_point];
        const double JxW = scratch.fe_values_ref.JxW(q_point);

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const unsigned int component_i = fe.system_to_component_index(i).first;
            const unsigned int i_group     = fe.system_to_base_index(i).first.first;

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
                      data.cell_matrix(i, j) += grad_Nx[i][component_i] * tau
                                                * grad_Nx[j][component_j] * JxW;
                  }
                // Next is the $\mathsf{\mathbf{k}}_{ \widetilde{p} \mathbf{u}}$ contribution
                else if ((i_group == p_dof) && (j_group == u_dof))
                  {
                    data.cell_matrix(i, j) += N[i] * det_F
                                              * (symm_grad_Nx[j]
                                                 * StandardTensors<dim>::I)
                                              * JxW;
                  }
                // and lastly the $\mathsf{\mathbf{k}}_{ \widetilde{J} \widetilde{p}}$
                // and $\mathsf{\mathbf{k}}_{ \widetilde{J} \widetilde{J}}$
                // contributions:
                else if ((i_group == J_dof) && (j_group == p_dof))
                  data.cell_matrix(i, j) -= N[i] * N[j] * JxW;
                else if ((i_group == j_group) && (i_group == J_dof))
                  data.cell_matrix(i, j) += N[i] * d2Psi_vol_dJ2 * N[j] * JxW;
                else
                  Assert((i_group <= J_dof) && (j_group <= J_dof),
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

// @sect4{Solid::assemble_system_rhs}
// The assembly of the right-hand side process is similar to the
// tangent matrix, so we will not describe it in too much detail.
// Note that since we are describing a problem with Neumann BCs,
// we will need the face normals and so must specify this in the
// update flags.
  template <int dim>
  void Solid<dim>::assemble_system_rhs()
  {
    timer.enter_subsection("Assemble system right-hand side");
    std::cout << " ASM_R " << std::flush;

    system_rhs = 0.0;

    const UpdateFlags uf_cell(update_values |
                              update_gradients |
                              update_JxW_values);
    const UpdateFlags uf_face(update_values |
                              update_normal_vectors |
                              update_JxW_values);

    PerTaskData_RHS per_task_data(dofs_per_cell);
    ScratchData_RHS scratch_data(fe, qf_cell, uf_cell, qf_face, uf_face);

    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    *this,
                    &Solid::assemble_system_rhs_one_cell,
                    &Solid::copy_local_to_global_rhs,
                    scratch_data,
                    per_task_data);

    timer.leave_subsection();
  }



  template <int dim>
  void Solid<dim>::copy_local_to_global_rhs(const PerTaskData_RHS &data)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      system_rhs(data.local_dof_indices[i]) += data.cell_rhs(i);
  }



  template <int dim>
  void
  Solid<dim>::assemble_system_rhs_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                           ScratchData_RHS &scratch,
                                           PerTaskData_RHS &data)
  {
    data.reset();
    scratch.reset();
    scratch.fe_values_ref.reinit(cell);
    cell->get_dof_indices(data.local_dof_indices);
    PointHistory<dim> *lqph =
      reinterpret_cast<PointHistory<dim>*>(cell->user_pointer());

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        const Tensor<2, dim> F_inv = lqph[q_point].get_F_inv();

        for (unsigned int k = 0; k < dofs_per_cell; ++k)
          {
            const unsigned int k_group = fe.system_to_base_index(k).first.first;

            if (k_group == u_dof)
              scratch.symm_grad_Nx[q_point][k]
                = symmetrize(scratch.fe_values_ref[u_fe].gradient(k, q_point)
                             * F_inv);
            else if (k_group == p_dof)
              scratch.Nx[q_point][k] = scratch.fe_values_ref[p_fe].value(k,
                                                                         q_point);
            else if (k_group == J_dof)
              scratch.Nx[q_point][k] = scratch.fe_values_ref[J_fe].value(k,
                                                                         q_point);
            else
              Assert(k_group <= J_dof, ExcInternalError());
          }
      }

    for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
      {
        const SymmetricTensor<2, dim> tau = lqph[q_point].get_tau();
        const double det_F = lqph[q_point].get_det_F();
        const double J_tilde = lqph[q_point].get_J_tilde();
        const double p_tilde = lqph[q_point].get_p_tilde();
        const double dPsi_vol_dJ = lqph[q_point].get_dPsi_vol_dJ();

        const std::vector<double>
        &N = scratch.Nx[q_point];
        const std::vector<SymmetricTensor<2, dim> >
        &symm_grad_Nx = scratch.symm_grad_Nx[q_point];
        const double JxW = scratch.fe_values_ref.JxW(q_point);

        // We first compute the contributions
        // from the internal forces.  Note, by
        // definition of the rhs as the negative
        // of the residual, these contributions
        // are subtracted.
        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          {
            const unsigned int i_group = fe.system_to_base_index(i).first.first;

            if (i_group == u_dof)
              data.cell_rhs(i) -= (symm_grad_Nx[i] * tau) * JxW;
            else if (i_group == p_dof)
              data.cell_rhs(i) -= N[i] * (det_F - J_tilde) * JxW;
            else if (i_group == J_dof)
              data.cell_rhs(i) -= N[i] * (dPsi_vol_dJ - p_tilde) * JxW;
            else
              Assert(i_group <= J_dof, ExcInternalError());
          }
      }

    // Next we assemble the Neumann contribution. We first check to see it the
    // cell face exists on a boundary on which a traction is applied and add
    // the contribution if this is the case.
    for (unsigned int face = 0; face < GeometryInfo<dim>::faces_per_cell;
         ++face)
      if (cell->face(face)->at_boundary() == true
          && cell->face(face)->boundary_indicator() == 6)
        {
          scratch.fe_face_values_ref.reinit(cell, face);

          for (unsigned int f_q_point = 0; f_q_point < n_q_points_f;
               ++f_q_point)
            {
              const Tensor<1, dim> &N =
                scratch.fe_face_values_ref.normal_vector(f_q_point);

              // Using the face normal at this quadrature point we specify the
              // traction in reference configuration. For this problem, a
              // defined pressure is applied in the reference configuration.
              // The direction of the applied traction is assumed not to
              // evolve with the deformation of the domain. The traction is
              // defined using the first Piola-Kirchhoff stress is simply
              // $\mathbf{t} = \mathbf{P}\mathbf{N} = [p_0 \mathbf{I}]
              // \mathbf{N} = p_0 \mathbf{N}$ We use the time variable to
              // linearly ramp up the pressure load.
              //
              // Note that the contributions to the right hand side vector we
              // compute here only exist in the displacement components of the
              // vector.
              static const double  p0        = -4.0
                                               /
                                               (parameters.scale * parameters.scale);
              const double         time_ramp = (time.current() / time.end());
              const double         pressure  = p0 * parameters.p_p0 * time_ramp;
              const Tensor<1, dim> traction  = pressure * N;

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

// @sect4{Solid::make_constraints}
// The constraints for this problem are simple to describe.
// However, since we are dealing with an iterative Newton method,
// it should be noted that any displacement constraints should only
// be specified at the zeroth iteration and subsequently no
// additional contributions are to be made since the constraints
// are already exactly satisfied.
  template <int dim>
  void Solid<dim>::make_constraints(const int &it_nr)
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
    const FEValuesExtractors::Scalar x_displacement(0);
    const FEValuesExtractors::Scalar y_displacement(1);
    const FEValuesExtractors::Scalar z_displacement(2);

    {
      const int boundary_id = 0;

      if (apply_dirichlet_bc == true)
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 fe.component_mask(x_displacement));
      else
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 fe.component_mask(x_displacement));
    }
    {
      const int boundary_id = 2;

      if (apply_dirichlet_bc == true)
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 fe.component_mask(y_displacement));
      else
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 fe.component_mask(y_displacement));
    }
    {
      const int boundary_id = 4;

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
    {
      const int boundary_id = 5;

      if (apply_dirichlet_bc == true)
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 (fe.component_mask(x_displacement)
                                                  |
                                                  fe.component_mask(y_displacement)));
      else
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 (fe.component_mask(x_displacement)
                                                  |
                                                  fe.component_mask(y_displacement)));
    }
    {
      const int boundary_id = 6;

      if (apply_dirichlet_bc == true)
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 (fe.component_mask(x_displacement)
                                                  |
                                                  fe.component_mask(y_displacement)));
      else
        VectorTools::interpolate_boundary_values(dof_handler_ref,
                                                 boundary_id,
                                                 ZeroFunction<dim>(n_components),
                                                 constraints,
                                                 (fe.component_mask(x_displacement)
                                                  |
                                                  fe.component_mask(y_displacement)));
    }

    constraints.close();
  }

// @sect4{Solid::solve_linear_system}
// Solving the entire block system is a bit problematic as there are no
// contributions to the $\mathsf{\mathbf{k}}_{ \widetilde{J} \widetilde{J}}$
// block, rendering it noninvertible.
// Since the pressure and dilatation variables DOFs are discontinuous, we can
// condense them out to form a smaller displacement-only system which
// we will then solve and subsequently post-process to retrieve the
// pressure and dilatation solutions.
//
// At the top, we allocate two temporary vectors to help with the static
// condensation, and variables to store the number of linear solver iterations
// and the (hopefully converged) residual.
//
// For the following, recall that
// @f{align*}
//  \mathbf{\mathsf{K}}_{\textrm{store}}
//:=
//  \begin{bmatrix}
//      \mathbf{\mathsf{K}}_{\textrm{con}}      &       \mathbf{\mathsf{K}}_{u\widetilde{p}}    & \mathbf{0}
//  \\  \mathbf{\mathsf{K}}_{\widetilde{p}u}    &       \mathbf{0}      &       \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1}
//  \\  \mathbf{0}      &       \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}                & \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
//  \end{bmatrix} \, .
// @f}
// and
//  @f{align*}
//              d \widetilde{\mathbf{\mathsf{p}}}
//              & = \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \bigl[
//                       \mathbf{\mathsf{F}}_{\widetilde{J}}
//                       - \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}} d \widetilde{\mathbf{\mathsf{J}}} \bigr]
//              \\ d \widetilde{\mathbf{\mathsf{J}}}
//              & = \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1} \bigl[
//                      \mathbf{\mathsf{F}}_{\widetilde{p}}
//                      - \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
//                      \bigr]
//               \\ \Rightarrow d \widetilde{\mathbf{\mathsf{p}}}
//              &=  \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{F}}_{\widetilde{J}}
//              - \underbrace{\bigl[\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
//              \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1}\bigr]}_{\overline{\mathbf{\mathsf{K}}}}\bigl[ \mathbf{\mathsf{F}}_{\widetilde{p}}
//              - \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}} \bigr]
//  @f}
//  and thus
//  @f[
//              \underbrace{\bigl[ \mathbf{\mathsf{K}}_{uu} + \overline{\overline{\mathbf{\mathsf{K}}}}~ \bigr]
//              }_{\mathbf{\mathsf{K}}_{\textrm{con}}} d \mathbf{\mathsf{u}}
//              =
//          \underbrace{
//              \Bigl[
//              \mathbf{\mathsf{F}}_{u}
//                      - \mathbf{\mathsf{K}}_{u\widetilde{p}} \bigl[ \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{F}}_{\widetilde{J}}
//                      - \overline{\mathbf{\mathsf{K}}}\mathbf{\mathsf{F}}_{\widetilde{p}} \bigr]
//              \Bigr]}_{\mathbf{\mathsf{F}}_{\textrm{con}}}
//  @f]
//  where
//  @f[
//              \overline{\overline{\mathbf{\mathsf{K}}}} :=
//                      \mathbf{\mathsf{K}}_{u\widetilde{p}} \overline{\mathbf{\mathsf{K}}} \mathbf{\mathsf{K}}_{\widetilde{p}u} \, .
//  @f]
  template <int dim>
  std::pair<unsigned int, double>
  Solid<dim>::solve_linear_system(BlockVector<double> &newton_update)
  {
    BlockVector<double> A(dofs_per_block);
    BlockVector<double> B(dofs_per_block);

    unsigned int lin_it = 0;
    double lin_res = 0.0;

    // In the first step of this function, we solve for the incremental
    // displacement $d\mathbf{u}$.  To this end, we perform static
    // condensation to make
    //    $\mathbf{\mathsf{K}}_{\textrm{con}}
    //    = \bigl[ \mathbf{\mathsf{K}}_{uu} + \overline{\overline{\mathbf{\mathsf{K}}}}~ \bigr]$
    // and put
    // $\mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}$
    // in the original $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$ block.
    // That is, we make $\mathbf{\mathsf{K}}_{\textrm{store}}$.
    {
      assemble_sc();

      //              $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //      =
      //              \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //              \mathsf{\mathbf{F}}_{\widetilde{p}}
      //              $
      tangent_matrix.block(p_dof, J_dof).vmult(A.block(J_dof),
                                               system_rhs.block(p_dof));
      //      $
      //      \mathsf{\mathbf{B}}_{\widetilde{J}}
      //      =
      //      \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //      \mathsf{\mathbf{F}}_{\widetilde{p}}
      //      $
      tangent_matrix.block(J_dof, J_dof).vmult(B.block(J_dof),
                                               A.block(J_dof));
      //      $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //      =
      //      \mathsf{\mathbf{F}}_{\widetilde{J}}
      //      -
      //      \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //      \mathsf{\mathbf{F}}_{\widetilde{p}}
      //      $
      A.block(J_dof).equ(1.0, system_rhs.block(J_dof), -1.0, B.block(J_dof));
      //      $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //      =
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
      //      [
      //      \mathsf{\mathbf{F}}_{\widetilde{J}}
      //      -
      //      \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //      \mathsf{\mathbf{F}}_{\widetilde{p}}
      //      ]
      //      $
      tangent_matrix.block(p_dof, J_dof).Tvmult(A.block(p_dof),
                                                A.block(J_dof));
      //      $
      //      \mathsf{\mathbf{A}}_{\mathbf{u}}
      //      =
      //      \mathsf{\mathbf{K}}_{\mathbf{u} \widetilde{p}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
      //      [
      //      \mathsf{\mathbf{F}}_{\widetilde{J}}
      //      -
      //      \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //      \mathsf{\mathbf{F}}_{\widetilde{p}}
      //      ]
      //      $
      tangent_matrix.block(u_dof, p_dof).vmult(A.block(u_dof),
                                               A.block(p_dof));
      //      $
      //      \mathsf{\mathbf{F}}_{\text{con}}
      //      =
      //      \mathsf{\mathbf{F}}_{\mathbf{u}}
      //      -
      //      \mathsf{\mathbf{K}}_{\mathbf{u} \widetilde{p}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{J} \widetilde{p}}
      //      [
      //      \mathsf{\mathbf{F}}_{\widetilde{J}}
      //      -
      //      \mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{J}}
      //      \mathsf{\mathbf{K}}^{-1}_{\widetilde{p} \widetilde{J}}
      //      \mathsf{\mathbf{K}}_{\widetilde{p}}
      //      ]
      //      $
      system_rhs.block(u_dof) -= A.block(u_dof);

      timer.enter_subsection("Linear solver");
      std::cout << " SLV " << std::flush;
      if (parameters.type_lin == "CG")
        {
          const int solver_its = tangent_matrix.block(u_dof, u_dof).m()
                                 * parameters.max_iterations_lin;
          const double tol_sol = parameters.tol_lin
                                 * system_rhs.block(u_dof).l2_norm();

          SolverControl solver_control(solver_its, tol_sol);

          GrowingVectorMemory<Vector<double> > GVM;
          SolverCG<Vector<double> > solver_CG(solver_control, GVM);

          // We've chosen by default a SSOR preconditioner as it appears to
          // provide the fastest solver convergence characteristics for this
          // problem on a single-thread machine.  However, this might not be
          // true for different problem sizes.
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

    timer.enter_subsection("Linear solver postprocessing");
    std::cout << " PP " << std::flush;

    // The next step after solving the displacement
    // problem is to post-process to get the
    // dilatation solution from the
    // substitution:
    //    $
    //     d \widetilde{\mathbf{\mathsf{J}}}
    //      = \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}^{-1} \bigl[
    //       \mathbf{\mathsf{F}}_{\widetilde{p}}
    //     - \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
    //      \bigr]
    //    $
    {
      //      $
      //      \mathbf{\mathsf{A}}_{\widetilde{p}}
      //      =
      //      \mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
      //      $
      tangent_matrix.block(p_dof, u_dof).vmult(A.block(p_dof),
                                               newton_update.block(u_dof));
      //      $
      //      \mathbf{\mathsf{A}}_{\widetilde{p}}
      //      =
      //      -\mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
      //      $
      A.block(p_dof) *= -1.0;
      //      $
      //      \mathbf{\mathsf{A}}_{\widetilde{p}}
      //      =
      //      \mathbf{\mathsf{F}}_{\widetilde{p}}
      //      -\mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
      //      $
      A.block(p_dof) += system_rhs.block(p_dof);
      //      $
      //      d\mathbf{\mathsf{\widetilde{J}}}
      //      =
      //      \mathbf{\mathsf{K}}^{-1}_{\widetilde{p}\widetilde{J}}
      //      [
      //      \mathbf{\mathsf{F}}_{\widetilde{p}}
      //      -\mathbf{\mathsf{K}}_{\widetilde{p}u} d \mathbf{\mathsf{u}}
      //      ]
      //      $
      tangent_matrix.block(p_dof, J_dof).vmult(newton_update.block(J_dof),
                                               A.block(p_dof));
    }

    // we insure here that any Dirichlet constraints
    // are distributed on the updated solution:
    constraints.distribute(newton_update);

    // Finally we solve for the pressure
    // update with the substitution:
    //    $
    //    d \widetilde{\mathbf{\mathsf{p}}}
    //     =
    //    \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
    //    \bigl[
    //     \mathbf{\mathsf{F}}_{\widetilde{J}}
    //      - \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
    //    d \widetilde{\mathbf{\mathsf{J}}}
    //    \bigr]
    //    $
    {
      //      $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //       =
      //      \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
      //      d \widetilde{\mathbf{\mathsf{J}}}
      //      $
      tangent_matrix.block(J_dof, J_dof).vmult(A.block(J_dof),
                                               newton_update.block(J_dof));
      //      $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //       =
      //      -\mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
      //      d \widetilde{\mathbf{\mathsf{J}}}
      //      $
      A.block(J_dof) *= -1.0;
      //      $
      //      \mathsf{\mathbf{A}}_{\widetilde{J}}
      //       =
      //      \mathsf{\mathbf{F}}_{\widetilde{J}}
      //      -
      //      \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
      //      d \widetilde{\mathbf{\mathsf{J}}}
      //      $
      A.block(J_dof) += system_rhs.block(J_dof);
      // and finally....
      //    $
      //    d \widetilde{\mathbf{\mathsf{p}}}
      //     =
      //    \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{p}}^{-1}
      //    \bigl[
      //     \mathbf{\mathsf{F}}_{\widetilde{J}}
      //      - \mathbf{\mathsf{K}}_{\widetilde{J}\widetilde{J}}
      //    d \widetilde{\mathbf{\mathsf{J}}}
      //    \bigr]
      //    $
      tangent_matrix.block(p_dof, J_dof).Tvmult(newton_update.block(p_dof),
                                                A.block(J_dof));
    }

    // We are now at the end, so we distribute all
    // constrained dofs back to the Newton
    // update:
    constraints.distribute(newton_update);

    timer.leave_subsection();

    return std::make_pair(lin_it, lin_res);
  }

// @sect4{Solid::assemble_system_SC}

// The static condensation process could be performed at a global level but we
// need the inverse of one of the blocks. However, since the pressure and
// dilatation variables are discontinuous, the static condensation (SC)
// operation can be done on a per-cell basis and we can produce the inverse of
// the block-diagonal $ \mathbf{\mathsf{K}}_{\widetilde{p}\widetilde{J}}$
  // block by inverting the local blocks. We can again
// use TBB to do this since each operation will be independent of one another.
//
// Using the TBB via the WorkStream class, we assemble the contributions to form
//  $
//  \mathbf{\mathsf{K}}_{\textrm{con}}
//  = \bigl[ \mathbf{\mathsf{K}}_{uu} + \overline{\overline{\mathbf{\mathsf{K}}}}~ \bigr]
//  $
// from each element's contributions. These
// contributions are then added to the global stiffness matrix. Given this
// description, the following two functions should be clear:
  template <int dim>
  void Solid<dim>::assemble_sc()
  {
    timer.enter_subsection("Perform static condensation");
    std::cout << " ASM_SC " << std::flush;

    PerTaskData_SC per_task_data(dofs_per_cell, element_indices_u.size(),
                                 element_indices_p.size(),
                                 element_indices_J.size());
    ScratchData_SC scratch_data;

    WorkStream::run(dof_handler_ref.begin_active(),
                    dof_handler_ref.end(),
                    *this,
                    &Solid::assemble_sc_one_cell,
                    &Solid::copy_local_to_global_sc,
                    scratch_data,
                    per_task_data);

    timer.leave_subsection();
  }


  template <int dim>
  void Solid<dim>::copy_local_to_global_sc(const PerTaskData_SC &data)
  {
    for (unsigned int i = 0; i < dofs_per_cell; ++i)
      for (unsigned int j = 0; j < dofs_per_cell; ++j)
        tangent_matrix.add(data.local_dof_indices[i],
                           data.local_dof_indices[j],
                           data.cell_matrix(i, j));
  }


// Now we describe the static condensation process.  As per usual, we must
// first find out which global numbers the degrees of freedom on this cell
// have and reset some data structures:
  template <int dim>
  void
  Solid<dim>::assemble_sc_one_cell(const typename DoFHandler<dim>::active_cell_iterator &cell,
                                   ScratchData_SC &scratch,
                                   PerTaskData_SC &data)
  {
    data.reset();
    scratch.reset();
    cell->get_dof_indices(data.local_dof_indices);

    // We now extract the contribution of the dofs associated with the current
    // cell to the global stiffness matrix.  The discontinuous nature of the
    // $\widetilde{p}$ and $\widetilde{J}$ interpolations mean that their is
    // no coupling of the local contributions at the global level. This is not
    // the case with the u dof.  In other words,
    // $\mathsf{\mathbf{k}}_{\widetilde{J} \widetilde{p}}$,
    // $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{p}}$ and
    // $\mathsf{\mathbf{k}}_{\widetilde{J} \widetilde{p}}$, when extracted
    // from the global stiffness matrix are the element contributions.  This
    // is not the case for $\mathsf{\mathbf{k}}_{\mathbf{u} \mathbf{u}}$
    //
    // Note: A lower-case symbol is used to denote element stiffness matrices.

    // Currently the matrix corresponding to
    // the dof associated with the current element
    // (denoted somewhat loosely as $\mathsf{\mathbf{k}}$)
    // is of the form:
    // @f{align*}
    //    \begin{bmatrix}
    //       \mathbf{\mathsf{k}}_{uu}  &  \mathbf{\mathsf{k}}_{u\widetilde{p}}    & \mathbf{0}
    //    \\ \mathbf{\mathsf{k}}_{\widetilde{p}u} & \mathbf{0}  &  \mathbf{\mathsf{k}}_{\widetilde{p}\widetilde{J}}
    //    \\ \mathbf{0}  &  \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{p}}  & \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{J}}
    //    \end{bmatrix}
    // @f}
    //
    // We now need to modify it such that it appear as
    // @f{align*}
    //    \begin{bmatrix}
    //       \mathbf{\mathsf{k}}_{\textrm{con}}   & \mathbf{\mathsf{k}}_{u\widetilde{p}}    & \mathbf{0}
    //    \\ \mathbf{\mathsf{k}}_{\widetilde{p}u} & \mathbf{0} & \mathbf{\mathsf{k}}_{\widetilde{p}\widetilde{J}}^{-1}
    //    \\ \mathbf{0} & \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{p}} & \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{J}}
    //    \end{bmatrix}
    // @f}
    // with $\mathbf{\mathsf{k}}_{\textrm{con}} = \bigl[ \mathbf{\mathsf{k}}_{uu} +\overline{\overline{\mathbf{\mathsf{k}}}}~ \bigr]$
    // where
    // $               \overline{\overline{\mathbf{\mathsf{k}}}} :=
    // \mathbf{\mathsf{k}}_{u\widetilde{p}} \overline{\mathbf{\mathsf{k}}} \mathbf{\mathsf{k}}_{\widetilde{p}u}
    // $
    // and
    // $
    //    \overline{\mathbf{\mathsf{k}}} =
    //     \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{p}}^{-1} \mathbf{\mathsf{k}}_{\widetilde{J}\widetilde{J}}
    //    \mathbf{\mathsf{k}}_{\widetilde{p}\widetilde{J}}^{-1}
    // $.
    //
    // At this point, we need to take note of
    // the fact that global data already exists
    // in the $\mathsf{\mathbf{K}}_{uu}$,
    // $\mathsf{\mathbf{K}}_{\widetilde{p} \widetilde{J}}$
    // and
    //  $\mathsf{\mathbf{K}}_{\widetilde{J} \widetilde{p}}$
    // sub-blocks.  So if we are to modify them, we must account for the data
    // that is already there (i.e. simply add to it or remove it if
    // necessary).  Since the copy_local_to_global operation is a "+="
    // operation, we need to take this into account
    //
    // For the $\mathsf{\mathbf{K}}_{uu}$ block in particular, this means that
    // contributions have been added from the surrounding cells, so we need to
    // be careful when we manipulate this block.  We can't just erase the
    // sub-blocks.
    //
    // This is the strategy we will employ to get the sub-blocks we want:
    //
    // - $ {\mathbf{\mathsf{k}}}_{\textrm{store}}$:
    // Since we don't have access to $\mathsf{\mathbf{k}}_{uu}$,
    // but we know its contribution is added to
    // the global $\mathsf{\mathbf{K}}_{uu}$ matrix, we just want
    // to add the element wise
    // static-condensation $\overline{\overline{\mathbf{\mathsf{k}}}}$.
    //
    // - $\mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}$:
    //                      Similarly, $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$ exists in
    //          the subblock. Since the copy
    //          operation is a += operation, we
    //          need to subtract the existing
    //          $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$
    //                      submatrix in addition to
    //          "adding" that which we wish to
    //          replace it with.
    //
    // - $\mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}$:
    //              Since the global matrix
    //          is symmetric, this block is the
    //          same as the one above and we
    //          can simply use
    //              $\mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}$
    //          as a substitute for this one.
    //
    // We first extract element data from the
    // system matrix. So first we get the
    // entire subblock for the cell, then
    // extract $\mathsf{\mathbf{k}}$
    // for the dofs associated with
    // the current element
    data.k_orig.extract_submatrix_from(tangent_matrix,
                                       data.local_dof_indices,
                                       data.local_dof_indices);
    // and next the local matrices for
    // $\mathsf{\mathbf{k}}_{ \widetilde{p} \mathbf{u}}$
    // $\mathsf{\mathbf{k}}_{ \widetilde{p} \widetilde{J}}$
    // and
    // $\mathsf{\mathbf{k}}_{ \widetilde{J} \widetilde{J}}$:
    data.k_pu.extract_submatrix_from(data.k_orig,
                                     element_indices_p,
                                     element_indices_u);
    data.k_pJ.extract_submatrix_from(data.k_orig,
                                     element_indices_p,
                                     element_indices_J);
    data.k_JJ.extract_submatrix_from(data.k_orig,
                                     element_indices_J,
                                     element_indices_J);

    // To get the inverse of $\mathsf{\mathbf{k}}_{\widetilde{p}
    // \widetilde{J}}$, we invert it directly.  This operation is relatively
    // inexpensive since $\mathsf{\mathbf{k}}_{\widetilde{p} \widetilde{J}}$
    // since block-diagonal.
    data.k_pJ_inv.invert(data.k_pJ);

    // Now we can make condensation terms to
    // add to the $\mathsf{\mathbf{k}}_{\mathbf{u} \mathbf{u}}$
    // block and put them in
    // the cell local matrix
    //    $
    //    \mathsf{\mathbf{A}}
    //    =
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
    //    \mathsf{\mathbf{k}}_{\widetilde{p} \mathbf{u}}
    //    $:
    data.k_pJ_inv.mmult(data.A, data.k_pu);
    //      $
    //      \mathsf{\mathbf{B}}
    //      =
    //      \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
    //      \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
    //      \mathsf{\mathbf{k}}_{\widetilde{p} \mathbf{u}}
    //      $
    data.k_JJ.mmult(data.B, data.A);
    //    $
    //    \mathsf{\mathbf{C}}
    //    =
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
    //    \mathsf{\mathbf{k}}_{\widetilde{p} \mathbf{u}}
    //    $
    data.k_pJ_inv.Tmmult(data.C, data.B);
    //    $
    //    \overline{\overline{\mathsf{\mathbf{k}}}}
    //    =
    //    \mathsf{\mathbf{k}}_{\mathbf{u} \widetilde{p}}
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{p}}
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{J} \widetilde{J}}
    //    \mathsf{\mathbf{k}}^{-1}_{\widetilde{p} \widetilde{J}}
    //    \mathsf{\mathbf{k}}_{\widetilde{p} \mathbf{u}}
    //    $
    data.k_pu.Tmmult(data.k_bbar, data.C);
    data.k_bbar.scatter_matrix_to(element_indices_u,
                                  element_indices_u,
                                  data.cell_matrix);

    // Next we place
    // $\mathsf{\mathbf{k}}^{-1}_{ \widetilde{p} \widetilde{J}}$
    // in the
    // $\mathsf{\mathbf{k}}_{ \widetilde{p} \widetilde{J}}$
    // block for post-processing.  Note again
    // that we need to remove the
    // contribution that already exists there.
    data.k_pJ_inv.add(-1.0, data.k_pJ);
    data.k_pJ_inv.scatter_matrix_to(element_indices_p,
                                    element_indices_J,
                                    data.cell_matrix);
  }

// @sect4{Solid::output_results}
// Here we present how the results are written to file to be viewed
// using ParaView or Visit. The method is similar to that shown in previous
// tutorials so will not be discussed in detail.
  template <int dim>
  void Solid<dim>::output_results() const
  {
    DataOut<dim> data_out;
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(dim,
                                  DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
    data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

    std::vector<std::string> solution_name(dim, "displacement");
    solution_name.push_back("pressure");
    solution_name.push_back("dilatation");

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
    MappingQEulerian<dim> q_mapping(degree, soln, dof_handler_ref);
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
  using namespace Step44;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, dealii::numbers::invalid_unsigned_int);

  try
    {
      deallog.depth_console(0);

      Solid<3> solid_3d("parameters.prm");
      solid_3d.run();
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
