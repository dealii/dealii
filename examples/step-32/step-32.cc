/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2008 - 2014 by the deal.II authors
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
 * Authors: Martin Kronbichler, Uppsala University,
 *          Wolfgang Bangerth, Texas A&M University,
 *          Timo Heister, University of Goettingen, 2008-2011
 */


// @sect3{Include files}

// The first task as usual is to include the functionality of these well-known
// deal.II library files and some C++ header files.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/work_stream.h>
#include <deal.II/base/timer.h>
#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_bicgstab.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
#include <deal.II/numerics/solution_transfer.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <limits>
#include <locale>
#include <string>

// This is the only include file that is new: It introduces the
// parallel::distributed::SolutionTransfer equivalent of the
// dealii::SolutionTransfer class to take a solution from on mesh to the next
// one upon mesh refinement, but in the case of parallel distributed
// triangulations:
#include <deal.II/distributed/solution_transfer.h>

// The following classes are used in parallel distributed computations and
// have all already been introduced in step-40:
#include <deal.II/base/index_set.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/grid_refinement.h>


// The next step is like in all previous tutorial programs: We put everything
// into a namespace of its own and then import the deal.II classes and
// functions into it:
namespace Step32
{
  using namespace dealii;

  // @sect3{Equation data}

  // In the following namespace, we define the various pieces of equation data
  // that describe the problem. This corresponds to the various aspects of
  // making the problem at least slightly realistic and that were exhaustively
  // discussed in the description of the testcase in the introduction.
  //
  // We start with a few coefficients that have constant values (the comment
  // after the value indicates its physical units):
  namespace EquationData
  {
    const double eta                   = 1e21;    /* Pa s       */
    const double kappa                 = 1e-6;    /* m^2 / s    */
    const double reference_density     = 3300;    /* kg / m^3   */
    const double reference_temperature = 293;     /* K          */
    const double expansion_coefficient = 2e-5;    /* 1/K        */
    const double specific_heat         = 1250;    /* J / K / kg */
    const double radiogenic_heating    = 7.4e-12; /* W / kg     */


    const double R0      = 6371000.-2890000.;     /* m          */
    const double R1      = 6371000.-  35000.;     /* m          */

    const double T0      = 4000+273;              /* K          */
    const double T1      =  700+273;              /* K          */


    // The next set of definitions are for functions that encode the density
    // as a function of temperature, the gravity vector, and the initial
    // values for the temperature. Again, all of these (along with the values
    // they compute) are discussed in the introduction:
    double density (const double temperature)
    {
      return (reference_density *
              (1 - expansion_coefficient * (temperature -
                                            reference_temperature)));
    }


    template <int dim>
    Tensor<1,dim> gravity_vector (const Point<dim> &p)
    {
      const double r = p.norm();
      return -(1.245e-6 * r + 7.714e13/r/r) * p / r;
    }



    template <int dim>
    class TemperatureInitialValues : public Function<dim>
    {
    public:
      TemperatureInitialValues () : Function<dim>(1) {}

      virtual double value (const Point<dim>   &p,
                            const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
                                 Vector<double>   &value) const;
    };



    template <int dim>
    double
    TemperatureInitialValues<dim>::value (const Point<dim>  &p,
                                          const unsigned int) const
    {
      const double r = p.norm();
      const double h = R1-R0;

      const double s = (r-R0)/h;
      const double q = (dim==3)?std::max(0.0,cos(numbers::PI*abs(p(2)/R1))):1.0;
      const double phi   = std::atan2(p(0),p(1));
      const double tau = s
                         +
                         0.2 * s * (1-s) * std::sin(6*phi) * q;

      return T0*(1.0-tau) + T1*tau;
    }


    template <int dim>
    void
    TemperatureInitialValues<dim>::vector_value (const Point<dim> &p,
                                                 Vector<double>   &values) const
    {
      for (unsigned int c=0; c<this->n_components; ++c)
        values(c) = TemperatureInitialValues<dim>::value (p, c);
    }


    // As mentioned in the introduction we need to rescale the pressure to
    // avoid the relative ill-conditioning of the momentum and mass
    // conservation equations. The scaling factor is $\frac{\eta}{L}$ where
    // $L$ was a typical length scale. By experimenting it turns out that a
    // good length scale is the diameter of plumes, which is around 10 km:
    const double pressure_scaling = eta / 10000;

    // The final number in this namespace is a constant that denotes the
    // number of seconds per (average, tropical) year. We use this only when
    // generating screen output: internally, all computations of this program
    // happen in SI units (kilogram, meter, seconds) but writing geological
    // times in seconds yields numbers that one can't relate to reality, and
    // so we convert to years using the factor defined here:
    const double year_in_seconds  = 60*60*24*365.2425;

  }



  // @sect3{Preconditioning the Stokes system}

  // This namespace implements the preconditioner. As discussed in the
  // introduction, this preconditioner differs in a number of key portions
  // from the one used in step-31. Specifically, it is a right preconditioner,
  // implementing the matrix
  // @f{align*}
  //   \left(\begin{array}{cc}A^{-1} & B^T
  //                        \\0 & S^{-1}
  // \end{array}\right)
  // @f}
  // where the two inverse matrix operations
  // are approximated by linear solvers or, if the right flag is given to the
  // constructor of this class, by a single AMG V-cycle for the velocity
  // block. The three code blocks of the <code>vmult</code> function implement
  // the multiplications with the three blocks of this preconditioner matrix
  // and should be self explanatory if you have read through step-31 or the
  // discussion of composing solvers in step-20.
  namespace LinearSolvers
  {
    template <class PreconditionerA, class PreconditionerMp>
    class BlockSchurPreconditioner : public Subscriptor
    {
    public:
      BlockSchurPreconditioner (const TrilinosWrappers::BlockSparseMatrix  &S,
                                const TrilinosWrappers::BlockSparseMatrix  &Spre,
                                const PreconditionerMp                     &Mppreconditioner,
                                const PreconditionerA                      &Apreconditioner,
                                const bool                                  do_solve_A)
        :
        stokes_matrix     (&S),
        stokes_preconditioner_matrix     (&Spre),
        mp_preconditioner (Mppreconditioner),
        a_preconditioner  (Apreconditioner),
        do_solve_A        (do_solve_A)
      {}

      void vmult (TrilinosWrappers::MPI::BlockVector       &dst,
                  const TrilinosWrappers::MPI::BlockVector &src) const
      {
        TrilinosWrappers::MPI::Vector utmp(src.block(0));

        {
          SolverControl solver_control(5000, 1e-6 * src.block(1).l2_norm());

          SolverCG<TrilinosWrappers::MPI::Vector> solver(solver_control);

          solver.solve(stokes_preconditioner_matrix->block(1,1),
                       dst.block(1), src.block(1),
                       mp_preconditioner);

          dst.block(1) *= -1.0;
        }

        {
          stokes_matrix->block(0,1).vmult(utmp, dst.block(1));
          utmp*=-1.0;
          utmp.add(src.block(0));
        }

        if (do_solve_A == true)
          {
            SolverControl solver_control(5000, utmp.l2_norm()*1e-2);
            TrilinosWrappers::SolverCG solver(solver_control);
            solver.solve(stokes_matrix->block(0,0), dst.block(0), utmp,
                         a_preconditioner);
          }
        else
          a_preconditioner.vmult (dst.block(0), utmp);
      }

    private:
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> stokes_matrix;
      const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> stokes_preconditioner_matrix;
      const PreconditionerMp &mp_preconditioner;
      const PreconditionerA  &a_preconditioner;
      const bool do_solve_A;
    };
  }



  // @sect3{Definition of assembly data structures}
  //
  // As described in the introduction, we will use the WorkStream mechanism
  // discussed in the @ref threads module to parallelize operations among the
  // processors of a single machine. The WorkStream class requires that data
  // is passed around in two kinds of data structures, one for scratch data
  // and one to pass data from the assembly function to the function that
  // copies local contributions into global objects.
  //
  // The following namespace (and the two sub-namespaces) contains a
  // collection of data structures that serve this purpose, one pair for each
  // of the four operations discussed in the introduction that we will want to
  // parallelize. Each assembly routine gets two sets of data: a Scratch array
  // that collects all the classes and arrays that are used for the
  // calculation of the cell contribution, and a CopyData array that keeps
  // local matrices and vectors which will be written into the global
  // matrix. Whereas CopyData is a container for the final data that is
  // written into the global matrices and vector (and, thus, absolutely
  // necessary), the Scratch arrays are merely there for performance reasons
  // &mdash; it would be much more expensive to set up a FEValues object on
  // each cell, than creating it only once and updating some derivative data.
  //
  // Step-31 had four assembly routines: One for the preconditioner matrix of
  // the Stokes system, one for the Stokes matrix and right hand side, one for
  // the temperature matrices and one for the right hand side of the
  // temperature equation. We here organize the scratch arrays and CopyData
  // objects for each of those four assembly components using a
  // <code>struct</code> environment (since we consider these as temporary
  // objects we pass around, rather than classes that implement functionality
  // of their own, though this is a more subjective point of view to
  // distinguish between <code>struct</code>s and <code>class</code>es).
  //
  // Regarding the Scratch objects, each struct is equipped with a constructor
  // that creates an FEValues object for a @ref FiniteElement "finite
  // element", a @ref Quadrature "quadrature formula", the @ref Mapping
  // "mapping" that describes the interpolation of curved boundaries, and some
  // @ref UpdateFlags "update flags".  Moreover, we manually implement a copy
  // constructor (since the FEValues class is not copyable by itself), and
  // provide some additional vector fields that are used to hold intermediate
  // data during the computation of local contributions.
  //
  // Let us start with the scratch arrays and, specifically, the one used for
  // assembly of the Stokes preconditioner:
  namespace Assembly
  {
    namespace Scratch
    {
      template <int dim>
      struct StokesPreconditioner
      {
        StokesPreconditioner (const FiniteElement<dim> &stokes_fe,
                              const Quadrature<dim>    &stokes_quadrature,
                              const Mapping<dim>       &mapping,
                              const UpdateFlags         update_flags);

        StokesPreconditioner (const StokesPreconditioner &data);


        FEValues<dim>               stokes_fe_values;

        std::vector<Tensor<2,dim> > grad_phi_u;
        std::vector<double>         phi_p;
      };

      template <int dim>
      StokesPreconditioner<dim>::
      StokesPreconditioner (const FiniteElement<dim> &stokes_fe,
                            const Quadrature<dim>    &stokes_quadrature,
                            const Mapping<dim>       &mapping,
                            const UpdateFlags         update_flags)
        :
        stokes_fe_values (mapping, stokes_fe, stokes_quadrature,
                          update_flags),
        grad_phi_u (stokes_fe.dofs_per_cell),
        phi_p (stokes_fe.dofs_per_cell)
      {}



      template <int dim>
      StokesPreconditioner<dim>::
      StokesPreconditioner (const StokesPreconditioner &scratch)
        :
        stokes_fe_values (scratch.stokes_fe_values.get_mapping(),
                          scratch.stokes_fe_values.get_fe(),
                          scratch.stokes_fe_values.get_quadrature(),
                          scratch.stokes_fe_values.get_update_flags()),
        grad_phi_u (scratch.grad_phi_u),
        phi_p (scratch.phi_p)
      {}



      // The next one is the scratch object used for the assembly of the full
      // Stokes system. Observe that we derive the StokesSystem scratch class
      // from the StokesPreconditioner class above. We do this because all the
      // objects that are necessary for the assembly of the preconditioner are
      // also needed for the actual matrix system and right hand side, plus
      // some extra data. This makes the program more compact. Note also that
      // the assembly of the Stokes system and the temperature right hand side
      // further down requires data from temperature and velocity,
      // respectively, so we actually need two FEValues objects for those two
      // cases.
      template <int dim>
      struct StokesSystem : public StokesPreconditioner<dim>
      {
        StokesSystem (const FiniteElement<dim> &stokes_fe,
                      const Mapping<dim>       &mapping,
                      const Quadrature<dim>    &stokes_quadrature,
                      const UpdateFlags         stokes_update_flags,
                      const FiniteElement<dim> &temperature_fe,
                      const UpdateFlags         temperature_update_flags);

        StokesSystem (const StokesSystem<dim> &data);


        FEValues<dim>                        temperature_fe_values;

        std::vector<Tensor<1,dim> >          phi_u;
        std::vector<SymmetricTensor<2,dim> > grads_phi_u;
        std::vector<double>                  div_phi_u;

        std::vector<double>                  old_temperature_values;
      };


      template <int dim>
      StokesSystem<dim>::
      StokesSystem (const FiniteElement<dim> &stokes_fe,
                    const Mapping<dim>       &mapping,
                    const Quadrature<dim>    &stokes_quadrature,
                    const UpdateFlags         stokes_update_flags,
                    const FiniteElement<dim> &temperature_fe,
                    const UpdateFlags         temperature_update_flags)
        :
        StokesPreconditioner<dim> (stokes_fe, stokes_quadrature,
                                   mapping,
                                   stokes_update_flags),
        temperature_fe_values (mapping, temperature_fe, stokes_quadrature,
                               temperature_update_flags),
        phi_u (stokes_fe.dofs_per_cell),
        grads_phi_u (stokes_fe.dofs_per_cell),
        div_phi_u (stokes_fe.dofs_per_cell),
        old_temperature_values (stokes_quadrature.size())
      {}


      template <int dim>
      StokesSystem<dim>::
      StokesSystem (const StokesSystem<dim> &scratch)
        :
        StokesPreconditioner<dim> (scratch),
        temperature_fe_values (scratch.temperature_fe_values.get_mapping(),
                               scratch.temperature_fe_values.get_fe(),
                               scratch.temperature_fe_values.get_quadrature(),
                               scratch.temperature_fe_values.get_update_flags()),
        phi_u (scratch.phi_u),
        grads_phi_u (scratch.grads_phi_u),
        div_phi_u (scratch.div_phi_u),
        old_temperature_values (scratch.old_temperature_values)
      {}


      // After defining the objects used in the assembly of the Stokes system,
      // we do the same for the assembly of the matrices necessary for the
      // temperature system. The general structure is very similar:
      template <int dim>
      struct TemperatureMatrix
      {
        TemperatureMatrix (const FiniteElement<dim> &temperature_fe,
                           const Mapping<dim>       &mapping,
                           const Quadrature<dim>    &temperature_quadrature);

        TemperatureMatrix (const TemperatureMatrix &data);


        FEValues<dim>               temperature_fe_values;

        std::vector<double>         phi_T;
        std::vector<Tensor<1,dim> > grad_phi_T;
      };


      template <int dim>
      TemperatureMatrix<dim>::
      TemperatureMatrix (const FiniteElement<dim> &temperature_fe,
                         const Mapping<dim>       &mapping,
                         const Quadrature<dim>    &temperature_quadrature)
        :
        temperature_fe_values (mapping,
                               temperature_fe, temperature_quadrature,
                               update_values    | update_gradients |
                               update_JxW_values),
        phi_T (temperature_fe.dofs_per_cell),
        grad_phi_T (temperature_fe.dofs_per_cell)
      {}


      template <int dim>
      TemperatureMatrix<dim>::
      TemperatureMatrix (const TemperatureMatrix &scratch)
        :
        temperature_fe_values (scratch.temperature_fe_values.get_mapping(),
                               scratch.temperature_fe_values.get_fe(),
                               scratch.temperature_fe_values.get_quadrature(),
                               scratch.temperature_fe_values.get_update_flags()),
        phi_T (scratch.phi_T),
        grad_phi_T (scratch.grad_phi_T)
      {}


      // The final scratch object is used in the assembly of the right hand
      // side of the temperature system. This object is significantly larger
      // than the ones above because a lot more quantities enter the
      // computation of the right hand side of the temperature equation. In
      // particular, the temperature values and gradients of the previous two
      // time steps need to be evaluated at the quadrature points, as well as
      // the velocities and the strain rates (i.e. the symmetric gradients of
      // the velocity) that enter the right hand side as friction heating
      // terms. Despite the number of terms, the following should be rather
      // self explanatory:
      template <int dim>
      struct TemperatureRHS
      {
        TemperatureRHS (const FiniteElement<dim> &temperature_fe,
                        const FiniteElement<dim> &stokes_fe,
                        const Mapping<dim>       &mapping,
                        const Quadrature<dim>    &quadrature);

        TemperatureRHS (const TemperatureRHS &data);


        FEValues<dim>                        temperature_fe_values;
        FEValues<dim>                        stokes_fe_values;

        std::vector<double>                  phi_T;
        std::vector<Tensor<1,dim> >          grad_phi_T;

        std::vector<Tensor<1,dim> >          old_velocity_values;
        std::vector<Tensor<1,dim> >          old_old_velocity_values;

        std::vector<SymmetricTensor<2,dim> > old_strain_rates;
        std::vector<SymmetricTensor<2,dim> > old_old_strain_rates;

        std::vector<double>                  old_temperature_values;
        std::vector<double>                  old_old_temperature_values;
        std::vector<Tensor<1,dim> >          old_temperature_grads;
        std::vector<Tensor<1,dim> >          old_old_temperature_grads;
        std::vector<double>                  old_temperature_laplacians;
        std::vector<double>                  old_old_temperature_laplacians;
      };


      template <int dim>
      TemperatureRHS<dim>::
      TemperatureRHS (const FiniteElement<dim> &temperature_fe,
                      const FiniteElement<dim> &stokes_fe,
                      const Mapping<dim>       &mapping,
                      const Quadrature<dim>    &quadrature)
        :
        temperature_fe_values (mapping,
                               temperature_fe, quadrature,
                               update_values    |
                               update_gradients |
                               update_hessians  |
                               update_quadrature_points |
                               update_JxW_values),
        stokes_fe_values (mapping,
                          stokes_fe, quadrature,
                          update_values | update_gradients),
        phi_T (temperature_fe.dofs_per_cell),
        grad_phi_T (temperature_fe.dofs_per_cell),

        old_velocity_values (quadrature.size()),
        old_old_velocity_values (quadrature.size()),
        old_strain_rates (quadrature.size()),
        old_old_strain_rates (quadrature.size()),

        old_temperature_values (quadrature.size()),
        old_old_temperature_values(quadrature.size()),
        old_temperature_grads(quadrature.size()),
        old_old_temperature_grads(quadrature.size()),
        old_temperature_laplacians(quadrature.size()),
        old_old_temperature_laplacians(quadrature.size())
      {}


      template <int dim>
      TemperatureRHS<dim>::
      TemperatureRHS (const TemperatureRHS &scratch)
        :
        temperature_fe_values (scratch.temperature_fe_values.get_mapping(),
                               scratch.temperature_fe_values.get_fe(),
                               scratch.temperature_fe_values.get_quadrature(),
                               scratch.temperature_fe_values.get_update_flags()),
        stokes_fe_values (scratch.stokes_fe_values.get_mapping(),
                          scratch.stokes_fe_values.get_fe(),
                          scratch.stokes_fe_values.get_quadrature(),
                          scratch.stokes_fe_values.get_update_flags()),
        phi_T (scratch.phi_T),
        grad_phi_T (scratch.grad_phi_T),

        old_velocity_values (scratch.old_velocity_values),
        old_old_velocity_values (scratch.old_old_velocity_values),
        old_strain_rates (scratch.old_strain_rates),
        old_old_strain_rates (scratch.old_old_strain_rates),

        old_temperature_values (scratch.old_temperature_values),
        old_old_temperature_values (scratch.old_old_temperature_values),
        old_temperature_grads (scratch.old_temperature_grads),
        old_old_temperature_grads (scratch.old_old_temperature_grads),
        old_temperature_laplacians (scratch.old_temperature_laplacians),
        old_old_temperature_laplacians (scratch.old_old_temperature_laplacians)
      {}
    }


    // The CopyData objects are even simpler than the Scratch objects as all
    // they have to do is to store the results of local computations until
    // they can be copied into the global matrix or vector objects. These
    // structures therefore only need to provide a constructor, a copy
    // operation, and some arrays for local matrix, local vectors and the
    // relation between local and global degrees of freedom (a.k.a.
    // <code>local_dof_indices</code>). Again, we have one such structure for
    // each of the four operations we will parallelize using the WorkStream
    // class:
    namespace CopyData
    {
      template <int dim>
      struct StokesPreconditioner
      {
        StokesPreconditioner (const FiniteElement<dim> &stokes_fe);
        StokesPreconditioner (const StokesPreconditioner &data);

        FullMatrix<double>          local_matrix;
        std::vector<types::global_dof_index> local_dof_indices;
      };

      template <int dim>
      StokesPreconditioner<dim>::
      StokesPreconditioner (const FiniteElement<dim> &stokes_fe)
        :
        local_matrix (stokes_fe.dofs_per_cell,
                      stokes_fe.dofs_per_cell),
        local_dof_indices (stokes_fe.dofs_per_cell)
      {}

      template <int dim>
      StokesPreconditioner<dim>::
      StokesPreconditioner (const StokesPreconditioner &data)
        :
        local_matrix (data.local_matrix),
        local_dof_indices (data.local_dof_indices)
      {}



      template <int dim>
      struct StokesSystem : public StokesPreconditioner<dim>
      {
        StokesSystem (const FiniteElement<dim> &stokes_fe);
        StokesSystem (const StokesSystem<dim> &data);

        Vector<double> local_rhs;
      };

      template <int dim>
      StokesSystem<dim>::
      StokesSystem (const FiniteElement<dim> &stokes_fe)
        :
        StokesPreconditioner<dim> (stokes_fe),
        local_rhs (stokes_fe.dofs_per_cell)
      {}

      template <int dim>
      StokesSystem<dim>::
      StokesSystem (const StokesSystem<dim> &data)
        :
        StokesPreconditioner<dim> (data),
        local_rhs (data.local_rhs)
      {}



      template <int dim>
      struct TemperatureMatrix
      {
        TemperatureMatrix (const FiniteElement<dim> &temperature_fe);
        TemperatureMatrix (const TemperatureMatrix &data);

        FullMatrix<double>          local_mass_matrix;
        FullMatrix<double>          local_stiffness_matrix;
        std::vector<types::global_dof_index>   local_dof_indices;
      };

      template <int dim>
      TemperatureMatrix<dim>::
      TemperatureMatrix (const FiniteElement<dim> &temperature_fe)
        :
        local_mass_matrix (temperature_fe.dofs_per_cell,
                           temperature_fe.dofs_per_cell),
        local_stiffness_matrix (temperature_fe.dofs_per_cell,
                                temperature_fe.dofs_per_cell),
        local_dof_indices (temperature_fe.dofs_per_cell)
      {}

      template <int dim>
      TemperatureMatrix<dim>::
      TemperatureMatrix (const TemperatureMatrix &data)
        :
        local_mass_matrix (data.local_mass_matrix),
        local_stiffness_matrix (data.local_stiffness_matrix),
        local_dof_indices (data.local_dof_indices)
      {}



      template <int dim>
      struct TemperatureRHS
      {
        TemperatureRHS (const FiniteElement<dim> &temperature_fe);
        TemperatureRHS (const TemperatureRHS &data);

        Vector<double>              local_rhs;
        std::vector<types::global_dof_index> local_dof_indices;
        FullMatrix<double>          matrix_for_bc;
      };

      template <int dim>
      TemperatureRHS<dim>::
      TemperatureRHS (const FiniteElement<dim> &temperature_fe)
        :
        local_rhs (temperature_fe.dofs_per_cell),
        local_dof_indices (temperature_fe.dofs_per_cell),
        matrix_for_bc (temperature_fe.dofs_per_cell,
                       temperature_fe.dofs_per_cell)
      {}

      template <int dim>
      TemperatureRHS<dim>::
      TemperatureRHS (const TemperatureRHS &data)
        :
        local_rhs (data.local_rhs),
        local_dof_indices (data.local_dof_indices),
        matrix_for_bc (data.matrix_for_bc)
      {}
    }
  }



  // @sect3{The <code>BoussinesqFlowProblem</code> class template}
  //
  // This is the declaration of the main class. It is very similar to step-31
  // but there are a number differences we will comment on below.
  //
  // The top of the class is essentially the same as in step-31, listing the
  // public methods and a set of private functions that do the heavy
  // lifting. Compared to step-31 there are only two additions to this
  // section: the function <code>get_cfl_number()</code> that computes the
  // maximum CFL number over all cells which we then compute the global time
  // step from, and the function <code>get_entropy_variation()</code> that is
  // used in the computation of the entropy stabilization. It is akin to the
  // <code>get_extrapolated_temperature_range()</code> we have used in step-31
  // for this purpose, but works on the entropy instead of the temperature
  // instead.
  template <int dim>
  class BoussinesqFlowProblem
  {
  public:
    struct Parameters;
    BoussinesqFlowProblem (Parameters &parameters);
    void run ();

  private:
    void setup_dofs ();
    void assemble_stokes_preconditioner ();
    void build_stokes_preconditioner ();
    void assemble_stokes_system ();
    void assemble_temperature_matrix ();
    void assemble_temperature_system (const double maximal_velocity);
    void project_temperature_field ();
    double get_maximal_velocity () const;
    double get_cfl_number () const;
    double get_entropy_variation (const double average_temperature) const;
    std::pair<double,double> get_extrapolated_temperature_range () const;
    void solve ();
    void output_results ();
    void refine_mesh (const unsigned int max_grid_level);

    double
    compute_viscosity(const std::vector<double>          &old_temperature,
                      const std::vector<double>          &old_old_temperature,
                      const std::vector<Tensor<1,dim> >  &old_temperature_grads,
                      const std::vector<Tensor<1,dim> >  &old_old_temperature_grads,
                      const std::vector<double>          &old_temperature_laplacians,
                      const std::vector<double>          &old_old_temperature_laplacians,
                      const std::vector<Tensor<1,dim> >  &old_velocity_values,
                      const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
                      const std::vector<SymmetricTensor<2,dim> >  &old_strain_rates,
                      const std::vector<SymmetricTensor<2,dim> >  &old_old_strain_rates,
                      const double                        global_u_infty,
                      const double                        global_T_variation,
                      const double                        average_temperature,
                      const double                        global_entropy_variation,
                      const double                        cell_diameter) const;

  public:

    // The first significant new component is the definition of a struct for
    // the parameters according to the discussion in the introduction. This
    // structure is initialized by reading from a parameter file during
    // construction of this object.
    struct Parameters
    {
      Parameters (const std::string &parameter_filename);

      static void declare_parameters (ParameterHandler &prm);
      void parse_parameters (ParameterHandler &prm);

      double       end_time;

      unsigned int initial_global_refinement;
      unsigned int initial_adaptive_refinement;

      bool         generate_graphical_output;
      unsigned int graphical_output_interval;

      unsigned int adaptive_refinement_interval;

      double       stabilization_alpha;
      double       stabilization_c_R;
      double       stabilization_beta;

      unsigned int stokes_velocity_degree;
      bool         use_locally_conservative_discretization;

      unsigned int temperature_degree;
    };

  private:
    Parameters                               &parameters;

    // The <code>pcout</code> (for <i>%parallel <code>std::cout</code></i>)
    // object is used to simplify writing output: each MPI process can use
    // this to generate output as usual, but since each of these processes
    // will (hopefully) produce the same output it will just be replicated
    // many times over; with the ConditionalOStream class, only the output
    // generated by one MPI process will actually be printed to screen,
    // whereas the output by all the other threads will simply be forgotten.
    ConditionalOStream                        pcout;

    // The following member variables will then again be similar to those in
    // step-31 (and to other tutorial programs). As mentioned in the
    // introduction, we fully distribute computations, so we will have to use
    // the parallel::distributed::Triangulation class (see step-40) but the
    // remainder of these variables is rather standard with two exceptions:
    //
    // - The <code>mapping</code> variable is used to denote a higher-order
    // polynomial mapping. As mentioned in the introduction, we use this
    // mapping when forming integrals through quadrature for all cells that
    // are adjacent to either the inner or outer boundaries of our domain
    // where the boundary is curved.
    //
    // - In a bit of naming confusion, you will notice below that some of the
    // variables from namespace TrilinosWrappers are taken from namespace
    // TrilinosWrappers::MPI (such as the right hand side vectors) whereas
    // others are not (such as the various matrices). For the matrices, we
    // happen to use the same class names for %parallel and sequential data
    // structures, i.e., all matrices will actually be considered %parallel
    // below. On the other hand, for vectors, only those from namespace
    // TrilinosWrappers::MPI are actually distributed. In particular, we will
    // frequently have to query velocities and temperatures at arbitrary
    // quadrature points; consequently, rather than importing ghost
    // information of a vector whenever we need access to degrees of freedom
    // that are relevant locally but owned by another processor, we solve
    // linear systems in %parallel but then immediately initialize a vector
    // including ghost entries of the solution for further processing. The
    // various <code>*_solution</code> vectors are therefore filled
    // immediately after solving their respective linear system in %parallel
    // and will always contain values for all @ref GlossLocallyRelevantDof
    // "locally relevant degrees of freedom"; the fully distributed vectors
    // that we obtain from the solution process and that only ever contain the
    // @ref GlossLocallyOwnedDof "locally owned degrees of freedom" are
    // destroyed immediately after the solution process and after we have
    // copied the relevant values into the member variable vectors.
    parallel::distributed::Triangulation<dim> triangulation;
    double                                    global_Omega_diameter;

    const MappingQ<dim>                       mapping;

    const FESystem<dim>                       stokes_fe;
    DoFHandler<dim>                           stokes_dof_handler;
    ConstraintMatrix                          stokes_constraints;

    TrilinosWrappers::BlockSparseMatrix       stokes_matrix;
    TrilinosWrappers::BlockSparseMatrix       stokes_preconditioner_matrix;

    TrilinosWrappers::MPI::BlockVector        stokes_solution;
    TrilinosWrappers::MPI::BlockVector        old_stokes_solution;
    TrilinosWrappers::MPI::BlockVector        stokes_rhs;


    FE_Q<dim>                                 temperature_fe;
    DoFHandler<dim>                           temperature_dof_handler;
    ConstraintMatrix                          temperature_constraints;

    TrilinosWrappers::SparseMatrix            temperature_mass_matrix;
    TrilinosWrappers::SparseMatrix            temperature_stiffness_matrix;
    TrilinosWrappers::SparseMatrix            temperature_matrix;

    TrilinosWrappers::MPI::Vector             temperature_solution;
    TrilinosWrappers::MPI::Vector             old_temperature_solution;
    TrilinosWrappers::MPI::Vector             old_old_temperature_solution;
    TrilinosWrappers::MPI::Vector             temperature_rhs;


    double                                    time_step;
    double                                    old_time_step;
    unsigned int                              timestep_number;

    std_cxx11::shared_ptr<TrilinosWrappers::PreconditionAMG>    Amg_preconditioner;
    std_cxx11::shared_ptr<TrilinosWrappers::PreconditionJacobi> Mp_preconditioner;
    std_cxx11::shared_ptr<TrilinosWrappers::PreconditionJacobi> T_preconditioner;

    bool                                      rebuild_stokes_matrix;
    bool                                      rebuild_stokes_preconditioner;
    bool                                      rebuild_temperature_matrices;
    bool                                      rebuild_temperature_preconditioner;

    // The next member variable, <code>computing_timer</code> is used to
    // conveniently account for compute time spent in certain "sections" of
    // the code that are repeatedly entered. For example, we will enter (and
    // leave) sections for Stokes matrix assembly and would like to accumulate
    // the run time spent in this section over all time steps. Every so many
    // time steps as well as at the end of the program (through the destructor
    // of the TimerOutput class) we will then produce a nice summary of the
    // times spent in the different sections into which we categorize the
    // run-time of this program.
    TimerOutput                               computing_timer;

    // After these member variables we have a number of auxiliary functions
    // that have been broken out of the ones listed above. Specifically, there
    // are first three functions that we call from <code>setup_dofs</code> and
    // then the ones that do the assembling of linear systems:
    void setup_stokes_matrix (const std::vector<IndexSet> &stokes_partitioning,
                              const std::vector<IndexSet> &stokes_relevant_partitioning);
    void setup_stokes_preconditioner (const std::vector<IndexSet> &stokes_partitioning,
                                      const std::vector<IndexSet> &stokes_relevant_partitioning);
    void setup_temperature_matrices (const IndexSet &temperature_partitioning,
                                     const IndexSet &temperature_relevant_partitioning);


    // Following the @ref MTWorkStream "task-based parallelization" paradigm,
    // we split all the assembly routines into two parts: a first part that
    // can do all the calculations on a certain cell without taking care of
    // other threads, and a second part (which is writing the local data into
    // the global matrices and vectors) which can be entered by only one
    // thread at a time. In order to implement that, we provide functions for
    // each of those two steps for all the four assembly routines that we use
    // in this program. The following eight functions do exactly this:
    void
    local_assemble_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                          Assembly::Scratch::StokesPreconditioner<dim> &scratch,
                                          Assembly::CopyData::StokesPreconditioner<dim> &data);

    void
    copy_local_to_global_stokes_preconditioner (const Assembly::CopyData::StokesPreconditioner<dim> &data);


    void
    local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                  Assembly::Scratch::StokesSystem<dim>  &scratch,
                                  Assembly::CopyData::StokesSystem<dim> &data);

    void
    copy_local_to_global_stokes_system (const Assembly::CopyData::StokesSystem<dim> &data);


    void
    local_assemble_temperature_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                       Assembly::Scratch::TemperatureMatrix<dim>  &scratch,
                                       Assembly::CopyData::TemperatureMatrix<dim> &data);

    void
    copy_local_to_global_temperature_matrix (const Assembly::CopyData::TemperatureMatrix<dim> &data);



    void
    local_assemble_temperature_rhs (const std::pair<double,double> global_T_range,
                                    const double                   global_max_velocity,
                                    const double                   global_entropy_variation,
                                    const typename DoFHandler<dim>::active_cell_iterator &cell,
                                    Assembly::Scratch::TemperatureRHS<dim> &scratch,
                                    Assembly::CopyData::TemperatureRHS<dim> &data);

    void
    copy_local_to_global_temperature_rhs (const Assembly::CopyData::TemperatureRHS<dim> &data);

    // Finally, we forward declare a member class that we will define later on
    // and that will be used to compute a number of quantities from our
    // solution vectors that we'd like to put into the output files for
    // visualization.
    class Postprocessor;
  };


  // @sect3{BoussinesqFlowProblem class implementation}

  // @sect4{BoussinesqFlowProblem::Parameters}
  //
  // Here comes the definition of the parameters for the Stokes problem. We
  // allow to set the end time for the simulation, the level of refinements
  // (both global and adaptive, which in the sum specify what maximum level
  // the cells are allowed to have), and the interval between refinements in
  // the time stepping.
  //
  // Then, we let the user specify constants for the stabilization parameters
  // (as discussed in the introduction), the polynomial degree for the Stokes
  // velocity space, whether to use the locally conservative discretization
  // based on FE_DGP elements for the pressure or not (FE_Q elements for
  // pressure), and the polynomial degree for the temperature interpolation.
  //
  // The constructor checks for a valid input file (if not, a file with
  // default parameters for the quantities is written), and eventually parses
  // the parameters.
  template <int dim>
  BoussinesqFlowProblem<dim>::Parameters::Parameters (const std::string &parameter_filename)
    :
    end_time (1e8),
    initial_global_refinement (2),
    initial_adaptive_refinement (2),
    adaptive_refinement_interval (10),
    stabilization_alpha (2),
    stabilization_c_R (0.11),
    stabilization_beta (0.078),
    stokes_velocity_degree (2),
    use_locally_conservative_discretization (true),
    temperature_degree (2)
  {
    ParameterHandler prm;
    BoussinesqFlowProblem<dim>::Parameters::declare_parameters (prm);

    std::ifstream parameter_file (parameter_filename.c_str());

    if (!parameter_file)
      {
        parameter_file.close ();

        std::ostringstream message;
        message << "Input parameter file <"
                << parameter_filename << "> not found. Creating a"
                << std::endl
                << "template file of the same name."
                << std::endl;

        std::ofstream parameter_out (parameter_filename.c_str());
        prm.print_parameters (parameter_out,
                              ParameterHandler::Text);

        AssertThrow (false, ExcMessage (message.str().c_str()));
      }

    const bool success = prm.read_input (parameter_file);
    AssertThrow (success, ExcMessage ("Invalid input parameter file."));

    parse_parameters (prm);
  }



  // Next we have a function that declares the parameters that we expect in
  // the input file, together with their data types, default values and a
  // description:
  template <int dim>
  void
  BoussinesqFlowProblem<dim>::Parameters::
  declare_parameters (ParameterHandler &prm)
  {
    prm.declare_entry ("End time", "1e8",
                       Patterns::Double (0),
                       "The end time of the simulation in years.");
    prm.declare_entry ("Initial global refinement", "2",
                       Patterns::Integer (0),
                       "The number of global refinement steps performed on "
                       "the initial coarse mesh, before the problem is first "
                       "solved there.");
    prm.declare_entry ("Initial adaptive refinement", "2",
                       Patterns::Integer (0),
                       "The number of adaptive refinement steps performed after "
                       "initial global refinement.");
    prm.declare_entry ("Time steps between mesh refinement", "10",
                       Patterns::Integer (1),
                       "The number of time steps after which the mesh is to be "
                       "adapted based on computed error indicators.");
    prm.declare_entry ("Generate graphical output", "false",
                       Patterns::Bool (),
                       "Whether graphical output is to be generated or not. "
                       "You may not want to get graphical output if the number "
                       "of processors is large.");
    prm.declare_entry ("Time steps between graphical output", "50",
                       Patterns::Integer (1),
                       "The number of time steps between each generation of "
                       "graphical output files.");

    prm.enter_subsection ("Stabilization parameters");
    {
      prm.declare_entry ("alpha", "2",
                         Patterns::Double (1, 2),
                         "The exponent in the entropy viscosity stabilization.");
      prm.declare_entry ("c_R", "0.11",
                         Patterns::Double (0),
                         "The c_R factor in the entropy viscosity "
                         "stabilization.");
      prm.declare_entry ("beta", "0.078",
                         Patterns::Double (0),
                         "The beta factor in the artificial viscosity "
                         "stabilization. An appropriate value for 2d is 0.052 "
                         "and 0.078 for 3d.");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Discretization");
    {
      prm.declare_entry ("Stokes velocity polynomial degree", "2",
                         Patterns::Integer (1),
                         "The polynomial degree to use for the velocity variables "
                         "in the Stokes system.");
      prm.declare_entry ("Temperature polynomial degree", "2",
                         Patterns::Integer (1),
                         "The polynomial degree to use for the temperature variable.");
      prm.declare_entry ("Use locally conservative discretization", "true",
                         Patterns::Bool (),
                         "Whether to use a Stokes discretization that is locally "
                         "conservative at the expense of a larger number of degrees "
                         "of freedom, or to go with a cheaper discretization "
                         "that does not locally conserve mass (although it is "
                         "globally conservative.");
    }
    prm.leave_subsection ();
  }



  // And then we need a function that reads the contents of the
  // ParameterHandler object we get by reading the input file and puts the
  // results into variables that store the values of the parameters we have
  // previously declared:
  template <int dim>
  void
  BoussinesqFlowProblem<dim>::Parameters::
  parse_parameters (ParameterHandler &prm)
  {
    end_time                    = prm.get_double ("End time");
    initial_global_refinement   = prm.get_integer ("Initial global refinement");
    initial_adaptive_refinement = prm.get_integer ("Initial adaptive refinement");

    adaptive_refinement_interval= prm.get_integer ("Time steps between mesh refinement");

    generate_graphical_output   = prm.get_bool ("Generate graphical output");
    graphical_output_interval   = prm.get_integer ("Time steps between graphical output");

    prm.enter_subsection ("Stabilization parameters");
    {
      stabilization_alpha = prm.get_double ("alpha");
      stabilization_c_R   = prm.get_double ("c_R");
      stabilization_beta  = prm.get_double ("beta");
    }
    prm.leave_subsection ();

    prm.enter_subsection ("Discretization");
    {
      stokes_velocity_degree = prm.get_integer ("Stokes velocity polynomial degree");
      temperature_degree     = prm.get_integer ("Temperature polynomial degree");
      use_locally_conservative_discretization
        = prm.get_bool ("Use locally conservative discretization");
    }
    prm.leave_subsection ();
  }




  // @sect4{BoussinesqFlowProblem::BoussinesqFlowProblem}
  //
  // The constructor of the problem is very similar to the constructor in
  // step-31. What is different is the %parallel communication: Trilinos uses
  // a message passing interface (MPI) for data distribution. When entering
  // the BoussinesqFlowProblem class, we have to decide how the parallization
  // is to be done. We choose a rather simple strategy and let all processors
  // that are running the program work together, specified by the communicator
  // <code>MPI_COMM_WORLD</code>. Next, we create the output stream (as we
  // already did in step-18) that only generates output on the first MPI
  // process and is completely forgetful on all others. The implementation of
  // this idea is to check the process number when <code>pcout</code> gets a
  // true argument, and it uses the <code>std::cout</code> stream for
  // output. If we are one processor five, for instance, then we will give a
  // <code>false</code> argument to <code>pcout</code>, which means that the
  // output of that processor will not be printed. With the exception of the
  // mapping object (for which we use polynomials of degree 4) all but the
  // final member variable are exactly the same as in step-31.
  //
  // This final object, the TimerOutput object, is then told to restrict
  // output to the <code>pcout</code> stream (processor 0), and then we
  // specify that we want to get a summary table at the end of the program
  // which shows us wallclock times (as opposed to CPU times). We will
  // manually also request intermediate summaries every so many time steps in
  // the <code>run()</code> function below.
  template <int dim>
  BoussinesqFlowProblem<dim>::BoussinesqFlowProblem (Parameters &parameters_)
    :
    parameters (parameters_),
    pcout (std::cout,
           (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)
            == 0)),

    triangulation (MPI_COMM_WORLD,
                   typename Triangulation<dim>::MeshSmoothing
                   (Triangulation<dim>::smoothing_on_refinement |
                    Triangulation<dim>::smoothing_on_coarsening)),

    mapping (4),

    stokes_fe (FE_Q<dim>(parameters.stokes_velocity_degree),
               dim,
               (parameters.use_locally_conservative_discretization
                ?
                static_cast<const FiniteElement<dim> &>
                (FE_DGP<dim>(parameters.stokes_velocity_degree-1))
                :
                static_cast<const FiniteElement<dim> &>
                (FE_Q<dim>(parameters.stokes_velocity_degree-1))),
               1),

    stokes_dof_handler (triangulation),

    temperature_fe (parameters.temperature_degree),
    temperature_dof_handler (triangulation),

    time_step (0),
    old_time_step (0),
    timestep_number (0),
    rebuild_stokes_matrix (true),
    rebuild_stokes_preconditioner (true),
    rebuild_temperature_matrices (true),
    rebuild_temperature_preconditioner (true),

    computing_timer (pcout,
                     TimerOutput::summary,
                     TimerOutput::wall_times)
  {}



  // @sect4{The BoussinesqFlowProblem helper functions}
  // @sect5{BoussinesqFlowProblem::get_maximal_velocity}

  // Except for two small details, the function to compute the global maximum
  // of the velocity is the same as in step-31. The first detail is actually
  // common to all functions that implement loops over all cells in the
  // triangulation: When operating in %parallel, each processor can only work
  // on a chunk of cells since each processor only has a certain part of the
  // entire triangulation. This chunk of cells that we want to work on is
  // identified via a so-called <code>subdomain_id</code>, as we also did in
  // step-18. All we need to change is hence to perform the cell-related
  // operations only on cells that are owned by the current process (as
  // opposed to ghost or artificial cells), i.e. for which the subdomain id
  // equals the number of the process ID. Since this is a commonly used
  // operation, there is a shortcut for this operation: we can ask whether the
  // cell is owned by the current processor using
  // <code>cell-@>is_locally_owned()</code>.
  //
  // The second difference is the way we calculate the maximum value. Before,
  // we could simply have a <code>double</code> variable that we checked
  // against on each quadrature point for each cell. Now, we have to be a bit
  // more careful since each processor only operates on a subset of
  // cells. What we do is to first let each processor calculate the maximum
  // among its cells, and then do a global communication operation
  // <code>Utilities::MPI::max</code> that computes the maximum value among
  // all the maximum values of the individual processors. MPI provides such a
  // call, but it's even simpler to use the respective function in namespace
  // Utilities::MPI using the MPI communicator object since that will do the
  // right thing even if we work without MPI and on a single machine only. The
  // call to <code>Utilities::MPI::max</code> needs two arguments, namely the
  // local maximum (input) and the MPI communicator, which is MPI_COMM_WORLD
  // in this example.
  template <int dim>
  double BoussinesqFlowProblem<dim>::get_maximal_velocity () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.stokes_velocity_degree);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping, stokes_fe, quadrature_formula, update_values);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    const FEValuesExtractors::Vector velocities (0);

    double max_local_velocity = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = stokes_dof_handler.begin_active(),
    endc = stokes_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[velocities].get_function_values (stokes_solution,
                                                     velocity_values);

          for (unsigned int q=0; q<n_q_points; ++q)
            max_local_velocity = std::max (max_local_velocity,
                                           velocity_values[q].norm());
        }

    return Utilities::MPI::max (max_local_velocity, MPI_COMM_WORLD);
  }


  // @sect5{BoussinesqFlowProblem::get_cfl_number}

  // The next function does something similar, but we now compute the CFL
  // number, i.e., maximal velocity on a cell divided by the cell
  // diameter. This number is necessary to determine the time step size, as we
  // use a semi-explicit time stepping scheme for the temperature equation
  // (see step-31 for a discussion). We compute it in the same way as above:
  // Compute the local maximum over all locally owned cells, then exchange it
  // via MPI to find the global maximum.
  template <int dim>
  double BoussinesqFlowProblem<dim>::get_cfl_number () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.stokes_velocity_degree);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping, stokes_fe, quadrature_formula, update_values);
    std::vector<Tensor<1,dim> > velocity_values(n_q_points);

    const FEValuesExtractors::Vector velocities (0);

    double max_local_cfl = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = stokes_dof_handler.begin_active(),
    endc = stokes_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values[velocities].get_function_values (stokes_solution,
                                                     velocity_values);

          double max_local_velocity = 1e-10;
          for (unsigned int q=0; q<n_q_points; ++q)
            max_local_velocity = std::max (max_local_velocity,
                                           velocity_values[q].norm());
          max_local_cfl = std::max(max_local_cfl,
                                   max_local_velocity / cell->diameter());
        }

    return Utilities::MPI::max (max_local_cfl, MPI_COMM_WORLD);
  }


  // @sect5{BoussinesqFlowProblem::get_entropy_variation}

  // Next comes the computation of the global entropy variation
  // $\|E(T)-\bar{E}(T)\|_\infty$ where the entropy $E$ is defined as
  // discussed in the introduction.  This is needed for the evaluation of the
  // stabilization in the temperature equation as explained in the
  // introduction. The entropy variation is actually only needed if we use
  // $\alpha=2$ as a power in the residual computation. The infinity norm is
  // computed by the maxima over quadrature points, as usual in discrete
  // computations.
  //
  // In order to compute this quantity, we first have to find the
  // space-average $\bar{E}(T)$ and then evaluate the maximum. However, that
  // means that we would need to perform two loops. We can avoid the overhead
  // by noting that $\|E(T)-\bar{E}(T)\|_\infty =
  // \max\big(E_{\textrm{max}}(T)-\bar{E}(T),
  // \bar{E}(T)-E_{\textrm{min}}(T)\big)$, i.e., the maximum out of the
  // deviation from the average entropy in positive and negative
  // directions. The four quantities we need for the latter formula (maximum
  // entropy, minimum entropy, average entropy, area) can all be evaluated in
  // the same loop over all cells, so we choose this simpler variant.
  template <int dim>
  double
  BoussinesqFlowProblem<dim>::get_entropy_variation (const double average_temperature) const
  {
    if (parameters.stabilization_alpha != 2)
      return 1.;

    const QGauss<dim> quadrature_formula (parameters.temperature_degree+1);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (temperature_fe, quadrature_formula,
                             update_values | update_JxW_values);
    std::vector<double> old_temperature_values(n_q_points);
    std::vector<double> old_old_temperature_values(n_q_points);

    // In the two functions above we computed the maximum of numbers that were
    // all non-negative, so we knew that zero was certainly a lower bound. On
    // the other hand, here we need to find the maximum deviation from the
    // average value, i.e., we will need to know the maximal and minimal
    // values of the entropy for which we don't a priori know the sign.
    //
    // To compute it, we can therefore start with the largest and smallest
    // possible values we can store in a double precision number: The minimum
    // is initialized with a bigger and the maximum with a smaller number than
    // any one that is going to appear. We are then guaranteed that these
    // numbers will be overwritten in the loop on the first cell or, if this
    // processor does not own any cells, in the communication step at the
    // latest. The following loop then computes the minimum and maximum local
    // entropy as well as keeps track of the area/volume of the part of the
    // domain we locally own and the integral over the entropy on it:
    double min_entropy = std::numeric_limits<double>::max(),
           max_entropy = -std::numeric_limits<double>::max(),
           area = 0,
           entropy_integrated = 0;

    typename DoFHandler<dim>::active_cell_iterator
    cell = temperature_dof_handler.begin_active(),
    endc = temperature_dof_handler.end();
    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          fe_values.reinit (cell);
          fe_values.get_function_values (old_temperature_solution,
                                         old_temperature_values);
          fe_values.get_function_values (old_old_temperature_solution,
                                         old_old_temperature_values);
          for (unsigned int q=0; q<n_q_points; ++q)
            {
              const double T = (old_temperature_values[q] +
                                old_old_temperature_values[q]) / 2;
              const double entropy = ((T-average_temperature) *
                                      (T-average_temperature));

              min_entropy = std::min (min_entropy, entropy);
              max_entropy = std::max (max_entropy, entropy);
              area += fe_values.JxW(q);
              entropy_integrated += fe_values.JxW(q) * entropy;
            }
        }

    // Now we only need to exchange data between processors: we need to sum
    // the two integrals (<code>area</code>, <code>entropy_integrated</code>),
    // and get the extrema for maximum and minimum. We could do this through
    // four different data exchanges, but we can it with two:
    // Utilities::MPI::sum also exists in a variant that takes an array of
    // values that are all to be summed up. And we can also utilize the
    // Utilities::MPI::max function by realizing that forming the minimum over
    // the minimal entropies equals forming the negative of the maximum over
    // the negative of the minimal entropies; this maximum can then be
    // combined with forming the maximum over the maximal entropies.
    const double local_sums[2]   = { entropy_integrated, area },
                                   local_maxima[2] = { -min_entropy, max_entropy };
    double global_sums[2], global_maxima[2];

    Utilities::MPI::sum (local_sums,   MPI_COMM_WORLD, global_sums);
    Utilities::MPI::max (local_maxima, MPI_COMM_WORLD, global_maxima);

    // Having computed everything this way, we can then compute the average
    // entropy and find the $L^\infty$ norm by taking the larger of the
    // deviation of the maximum or minimum from the average:
    const double average_entropy = global_sums[0] / global_sums[1];
    const double entropy_diff = std::max(global_maxima[1] - average_entropy,
                                         average_entropy - (-global_maxima[0]));
    return entropy_diff;
  }



  // @sect5{BoussinesqFlowProblem::get_extrapolated_temperature_range}

  // The next function computes the minimal and maximal value of the
  // extrapolated temperature over the entire domain. Again, this is only a
  // slightly modified version of the respective function in step-31. As in
  // the function above, we collect local minima and maxima and then compute
  // the global extrema using the same trick as above.
  //
  // As already discussed in step-31, the function needs to distinguish
  // between the first and all following time steps because it uses a higher
  // order temperature extrapolation scheme when at least two previous time
  // steps are available.
  template <int dim>
  std::pair<double,double>
  BoussinesqFlowProblem<dim>::get_extrapolated_temperature_range () const
  {
    const QIterated<dim> quadrature_formula (QTrapez<1>(),
                                             parameters.temperature_degree);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (mapping, temperature_fe, quadrature_formula,
                             update_values);
    std::vector<double> old_temperature_values(n_q_points);
    std::vector<double> old_old_temperature_values(n_q_points);

    double min_local_temperature = std::numeric_limits<double>::max(),
           max_local_temperature = -std::numeric_limits<double>::max();

    if (timestep_number != 0)
      {
        typename DoFHandler<dim>::active_cell_iterator
        cell = temperature_dof_handler.begin_active(),
        endc = temperature_dof_handler.end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values.get_function_values (old_temperature_solution,
                                             old_temperature_values);
              fe_values.get_function_values (old_old_temperature_solution,
                                             old_old_temperature_values);

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double temperature =
                    (1. + time_step/old_time_step) * old_temperature_values[q]-
                    time_step/old_time_step * old_old_temperature_values[q];

                  min_local_temperature = std::min (min_local_temperature,
                                                    temperature);
                  max_local_temperature = std::max (max_local_temperature,
                                                    temperature);
                }
            }
      }
    else
      {
        typename DoFHandler<dim>::active_cell_iterator
        cell = temperature_dof_handler.begin_active(),
        endc = temperature_dof_handler.end();
        for (; cell!=endc; ++cell)
          if (cell->is_locally_owned())
            {
              fe_values.reinit (cell);
              fe_values.get_function_values (old_temperature_solution,
                                             old_temperature_values);

              for (unsigned int q=0; q<n_q_points; ++q)
                {
                  const double temperature = old_temperature_values[q];

                  min_local_temperature = std::min (min_local_temperature,
                                                    temperature);
                  max_local_temperature = std::max (max_local_temperature,
                                                    temperature);
                }
            }
      }

    double local_extrema[2] = { -min_local_temperature,
                                max_local_temperature
                              };
    double global_extrema[2];
    Utilities::MPI::max (local_extrema, MPI_COMM_WORLD, global_extrema);

    return std::make_pair(-global_extrema[0], global_extrema[1]);
  }


  // @sect5{BoussinesqFlowProblem::compute_viscosity}

  // The function that calculates the viscosity is purely local and so needs
  // no communication at all. It is mostly the same as in step-31 but with an
  // updated formulation of the viscosity if $\alpha=2$ is chosen:
  template <int dim>
  double
  BoussinesqFlowProblem<dim>::
  compute_viscosity (const std::vector<double>          &old_temperature,
                     const std::vector<double>          &old_old_temperature,
                     const std::vector<Tensor<1,dim> >  &old_temperature_grads,
                     const std::vector<Tensor<1,dim> >  &old_old_temperature_grads,
                     const std::vector<double>          &old_temperature_laplacians,
                     const std::vector<double>          &old_old_temperature_laplacians,
                     const std::vector<Tensor<1,dim> >  &old_velocity_values,
                     const std::vector<Tensor<1,dim> >  &old_old_velocity_values,
                     const std::vector<SymmetricTensor<2,dim> >  &old_strain_rates,
                     const std::vector<SymmetricTensor<2,dim> >  &old_old_strain_rates,
                     const double                        global_u_infty,
                     const double                        global_T_variation,
                     const double                        average_temperature,
                     const double                        global_entropy_variation,
                     const double                        cell_diameter) const
  {
    if (global_u_infty == 0)
      return 5e-3 * cell_diameter;

    const unsigned int n_q_points = old_temperature.size();

    double max_residual = 0;
    double max_velocity = 0;

    for (unsigned int q=0; q < n_q_points; ++q)
      {
        const Tensor<1,dim> u = (old_velocity_values[q] +
                                 old_old_velocity_values[q]) / 2;

        const SymmetricTensor<2,dim> strain_rate = (old_strain_rates[q] +
                                                    old_old_strain_rates[q]) / 2;

        const double T = (old_temperature[q] + old_old_temperature[q]) / 2;
        const double dT_dt = (old_temperature[q] - old_old_temperature[q])
                             / old_time_step;
        const double u_grad_T = u * (old_temperature_grads[q] +
                                     old_old_temperature_grads[q]) / 2;

        const double kappa_Delta_T = EquationData::kappa
                                     * (old_temperature_laplacians[q] +
                                        old_old_temperature_laplacians[q]) / 2;
        const double gamma
          = ((EquationData::radiogenic_heating * EquationData::density(T)
              +
              2 * EquationData::eta * strain_rate * strain_rate) /
             (EquationData::density(T) * EquationData::specific_heat));

        double residual
          = std::abs(dT_dt + u_grad_T - kappa_Delta_T - gamma);
        if (parameters.stabilization_alpha == 2)
          residual *= std::abs(T - average_temperature);

        max_residual = std::max (residual,        max_residual);
        max_velocity = std::max (std::sqrt (u*u), max_velocity);
      }

    const double max_viscosity = (parameters.stabilization_beta *
                                  max_velocity * cell_diameter);
    if (timestep_number == 0)
      return max_viscosity;
    else
      {
        Assert (old_time_step > 0, ExcInternalError());

        double entropy_viscosity;
        if (parameters.stabilization_alpha == 2)
          entropy_viscosity = (parameters.stabilization_c_R *
                               cell_diameter * cell_diameter *
                               max_residual /
                               global_entropy_variation);
        else
          entropy_viscosity = (parameters.stabilization_c_R *
                               cell_diameter * global_Omega_diameter *
                               max_velocity * max_residual /
                               (global_u_infty * global_T_variation));

        return std::min (max_viscosity, entropy_viscosity);
      }
  }



  // @sect5{BoussinesqFlowProblem::project_temperature_field}

  // This function is new compared to step-31. What is does is to re-implement
  // the library function <code>VectorTools::project()</code> for an MPI-based
  // parallelization, a function we used for generating an initial vector for
  // temperature based on some initial function. The library function only
  // works with shared memory but doesn't know how to utilize multiple
  // machines coupled through MPI to compute the projected field. The details
  // of a <code>project()</code> function are not very difficult. All we do is
  // to use a mass matrix and put the evaluation of the initial value function
  // on the right hand side. The mass matrix for temperature we can simply
  // generate using the respective assembly function, so all we need to do
  // here is to create the right hand side and do a CG solve. The assembly
  // function does a loop over all cells and evaluates the function in the
  // <code>EquationData</code> namespace, and does this only on cells owned by
  // the respective processor. The implementation of this assembly differs
  // from the assembly we do for the principal assembly functions further down
  // (which include thread-based parallelization with the WorkStream
  // concept). Here we chose to keep things simple (keeping in mind that this
  // function is also only called once at the beginning of the program, not in
  // every time step), and generating the right hand side is cheap anyway so
  // we won't even notice that this part is not parallized by threads.
  //
  // Regarding the implementation of inhomogeneous Dirichlet boundary
  // conditions: Since we use the temperature ConstraintMatrix, we could apply
  // the boundary conditions directly when building the respective matrix and
  // right hand side. In this case, the boundary conditions are inhomogeneous,
  // which makes this procedure somewhat tricky since we get the matrix from
  // some other function that uses its own integration and assembly
  // loop. However, the correct imposition of boundary conditions needs the
  // matrix data we work on plus the right hand side simultaneously, since the
  // right hand side is created by Gaussian elimination on the matrix rows. In
  // order to not introduce the matrix assembly at this place, but still
  // having the matrix data available, we choose to create a dummy matrix
  // <code>matrix_for_bc</code> that we only fill with data when we need it
  // for imposing boundary conditions. These positions are exactly those where
  // we have an inhomogeneous entry in the ConstraintMatrix. There are only a
  // few such positions (on the boundary DoFs), so it is still much cheaper to
  // use this function than to create the full matrix here. To implement this,
  // we ask the constraint matrix whether the DoF under consideration is
  // inhomogeneously constrained. In that case, we generate the respective
  // matrix column that we need for creating the correct right hand side. Note
  // that this (manually generated) matrix entry needs to be exactly the entry
  // that we would fill the matrix with &mdash; otherwise, this will not work.
  template <int dim>
  void BoussinesqFlowProblem<dim>::project_temperature_field ()
  {
    assemble_temperature_matrix ();

    QGauss<dim> quadrature(parameters.temperature_degree+2);
    UpdateFlags update_flags = UpdateFlags(update_values   |
                                           update_quadrature_points |
                                           update_JxW_values);
    FEValues<dim> fe_values (mapping, temperature_fe, quadrature, update_flags);

    const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                       n_q_points    = fe_values.n_quadrature_points;

    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    Vector<double> cell_vector (dofs_per_cell);
    FullMatrix<double> matrix_for_bc (dofs_per_cell, dofs_per_cell);

    std::vector<double> rhs_values(n_q_points);

    TrilinosWrappers::MPI::Vector
    rhs (temperature_mass_matrix.row_partitioner()),
        solution (temperature_mass_matrix.row_partitioner());

    const EquationData::TemperatureInitialValues<dim> initial_temperature;

    typename DoFHandler<dim>::active_cell_iterator
    cell = temperature_dof_handler.begin_active(),
    endc = temperature_dof_handler.end();

    for (; cell!=endc; ++cell)
      if (cell->is_locally_owned())
        {
          cell->get_dof_indices (local_dof_indices);
          fe_values.reinit (cell);

          initial_temperature.value_list (fe_values.get_quadrature_points(),
                                          rhs_values);

          cell_vector = 0;
          matrix_for_bc = 0;
          for (unsigned int point=0; point<n_q_points; ++point)
            for (unsigned int i=0; i<dofs_per_cell; ++i)
              {
                cell_vector(i) += rhs_values[point] *
                                  fe_values.shape_value(i,point) *
                                  fe_values.JxW(point);
                if (temperature_constraints.is_inhomogeneously_constrained(local_dof_indices[i]))
                  {
                    for (unsigned int j=0; j<dofs_per_cell; ++j)
                      matrix_for_bc(j,i) += fe_values.shape_value(i,point) *
                                            fe_values.shape_value(j,point) *
                                            fe_values.JxW(point);
                  }
              }

          temperature_constraints.distribute_local_to_global (cell_vector,
                                                              local_dof_indices,
                                                              rhs,
                                                              matrix_for_bc);
        }

    rhs.compress (VectorOperation::add);

    // Now that we have the right linear system, we solve it using the CG
    // method with a simple Jacobi preconditioner:
    SolverControl solver_control(5*rhs.size(), 1e-12*rhs.l2_norm());
    SolverCG<TrilinosWrappers::MPI::Vector> cg(solver_control);

    TrilinosWrappers::PreconditionJacobi preconditioner_mass;
    preconditioner_mass.initialize(temperature_mass_matrix, 1.3);

    cg.solve (temperature_mass_matrix, solution, rhs, preconditioner_mass);

    temperature_constraints.distribute (solution);

    // Having so computed the current temperature field, let us set the member
    // variable that holds the temperature nodes. Strictly speaking, we really
    // only need to set <code>old_temperature_solution</code> since the first
    // thing we will do is to compute the Stokes solution that only requires
    // the previous time step's temperature field. That said, nothing good can
    // come from not initializing the other vectors as well (especially since
    // it's a relatively cheap operation and we only have to do it once at the
    // beginning of the program) if we ever want to extend our numerical
    // method or physical model, and so we initialize
    // <code>temperature_solution</code> and
    // <code>old_old_temperature_solution</code> as well. As a sidenote, while
    // the <code>solution</code> vector is strictly distributed (i.e. each
    // processor only stores a mutually exclusive subset of elements), the
    // assignment makes sure that the vectors on the left hand side (which
    // where initialized to contain ghost elements as well) also get the
    // correct ghost elements. In other words, the assignment here requires
    // communication between processors:
    temperature_solution = solution;
    old_temperature_solution = solution;
    old_old_temperature_solution = solution;
  }




  // @sect4{The BoussinesqFlowProblem setup functions}

  // The following three functions set up the Stokes matrix, the matrix used
  // for the Stokes preconditioner, and the temperature matrix. The code is
  // mostly the same as in step-31, but it has been broken out into three
  // functions of their own for simplicity.
  //
  // The main functional difference between the code here and that in step-31
  // is that the matrices we want to set up are distributed across multiple
  // processors. Since we still want to build up the sparsity pattern first
  // for efficiency reasons, we could continue to build the <i>entire</i>
  // sparsity pattern as a BlockCompressedSimpleSparsityPattern, as we did in
  // step-31. However, that would be inefficient: every processor would build
  // the same sparsity pattern, but only initialize a small part of the matrix
  // using it. It also violates the principle that every processor should only
  // work on those cells it owns (and, if necessary the layer of ghost cells
  // around it).
  //
  // Rather, we use an object of type TrilinosWrappers::BlockSparsityPattern,
  // which is (obviously) a wrapper around a sparsity pattern object provided
  // by Trilinos. The advantage is that the Trilinos sparsity pattern class
  // can communicate across multiple processors: if this processor fills in
  // all the nonzero entries that result from the cells it owns, and every
  // other processor does so as well, then at the end after some MPI
  // communication initiated by the <code>compress()</code> call, we will have
  // the globally assembled sparsity pattern available with which the global
  // matrix can be initialized.
  //
  // There is one important aspect when initializing Trilinos sparsity
  // patterns in parallel: In addition to specifying the locally owned rows
  // and columns of the matrices via the @p stokes_partitioning index set, we
  // also supply information about all the rows we are possibly going to write
  // into when assembling on a certain processor. The set of locally relevant
  // rows contains all such rows (possibly also a few unnecessary ones, but it
  // is difficult to find the exact row indices before actually getting
  // indices on all cells and resolving constraints). This additional
  // information allows to exactly determine the structure for the
  // off-processor data found during assembly. While Trilinos matrices are
  // able to collect this information on the fly as well (when initializing
  // them from some other reinit method), it is less efficient and leads to
  // problems when assembling matrices with multiple threads. In this program,
  // we pessimistically assume that only one processor at a time can write
  // into the matrix while assembly (whereas the computation is parallel),
  // which is fine for Trilinos matrices. In practice, one can do better by
  // hinting WorkStream at cells that do not share vertices, allowing for
  // parallelism among those cells (see the graph coloring algorithms and
  // WorkStream with colored iterators argument). However, that only works
  // when only one MPI processor is present because Trilinos' internal data
  // structures for accumulating off-processor data on the fly are not thread
  // safe. With the initialization presented here, there is no such problem
  // and one could safely introduce graph coloring for this algorithm.
  //
  // The only other change we need to make is to tell the
  // DoFTools::make_sparsity_pattern() function that it is only supposed to
  // work on a subset of cells, namely the ones whose
  // <code>subdomain_id</code> equals the number of the current processor, and
  // to ignore all other cells.
  //
  // This strategy is replicated across all three of the following functions.
  //
  // Note that Trilinos matrices store the information contained in the
  // sparsity patterns, so we can safely release the <code>sp</code> variable
  // once the matrix has been given the sparsity structure.
  template <int dim>
  void BoussinesqFlowProblem<dim>::
  setup_stokes_matrix (const std::vector<IndexSet> &stokes_partitioning,
                       const std::vector<IndexSet> &stokes_relevant_partitioning)
  {
    stokes_matrix.clear ();

    TrilinosWrappers::BlockSparsityPattern sp(stokes_partitioning, stokes_partitioning,
                                              stokes_relevant_partitioning,
                                              MPI_COMM_WORLD);

    Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
    for (unsigned int c=0; c<dim+1; ++c)
      for (unsigned int d=0; d<dim+1; ++d)
        if (! ((c==dim) && (d==dim)))
          coupling[c][d] = DoFTools::always;
        else
          coupling[c][d] = DoFTools::none;

    DoFTools::make_sparsity_pattern (stokes_dof_handler,
                                     coupling, sp,
                                     stokes_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(MPI_COMM_WORLD));
    sp.compress();

    stokes_matrix.reinit (sp);
  }



  template <int dim>
  void BoussinesqFlowProblem<dim>::
  setup_stokes_preconditioner (const std::vector<IndexSet> &stokes_partitioning,
                               const std::vector<IndexSet> &stokes_relevant_partitioning)
  {
    Amg_preconditioner.reset ();
    Mp_preconditioner.reset ();

    stokes_preconditioner_matrix.clear ();

    TrilinosWrappers::BlockSparsityPattern sp(stokes_partitioning, stokes_partitioning,
                                              stokes_relevant_partitioning,
                                              MPI_COMM_WORLD);

    Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
    for (unsigned int c=0; c<dim+1; ++c)
      for (unsigned int d=0; d<dim+1; ++d)
        if (c == d)
          coupling[c][d] = DoFTools::always;
        else
          coupling[c][d] = DoFTools::none;

    DoFTools::make_sparsity_pattern (stokes_dof_handler,
                                     coupling, sp,
                                     stokes_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(MPI_COMM_WORLD));
    sp.compress();

    stokes_preconditioner_matrix.reinit (sp);
  }


  template <int dim>
  void BoussinesqFlowProblem<dim>::
  setup_temperature_matrices (const IndexSet &temperature_partitioner,
                              const IndexSet &temperature_relevant_partitioner)
  {
    T_preconditioner.reset ();
    temperature_mass_matrix.clear ();
    temperature_stiffness_matrix.clear ();
    temperature_matrix.clear ();

    TrilinosWrappers::SparsityPattern sp(temperature_partitioner,
                                         temperature_partitioner,
                                         temperature_relevant_partitioner,
                                         MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern (temperature_dof_handler, sp,
                                     temperature_constraints, false,
                                     Utilities::MPI::
                                     this_mpi_process(MPI_COMM_WORLD));
    sp.compress();

    temperature_matrix.reinit (sp);
    temperature_mass_matrix.reinit (sp);
    temperature_stiffness_matrix.reinit (sp);
  }



  // The remainder of the setup function (after splitting out the three
  // functions above) mostly has to deal with the things we need to do for
  // parallelization across processors. Because setting all of this up is a
  // significant compute time expense of the program, we put everything we do
  // here into a timer group so that we can get summary information about the
  // fraction of time spent in this part of the program at its end.
  //
  // At the top as usual we enumerate degrees of freedom and sort them by
  // component/block, followed by writing their numbers to the screen from
  // processor zero. The DoFHandler::distributed_dofs() function, when applied
  // to a parallel::distributed::Triangulation object, sorts degrees of
  // freedom in such a way that all degrees of freedom associated with
  // subdomain zero come before all those associated with subdomain one,
  // etc. For the Stokes part, this entails, however, that velocities and
  // pressures become intermixed, but this is trivially solved by sorting
  // again by blocks; it is worth noting that this latter operation leaves the
  // relative ordering of all velocities and pressures alone, i.e. within the
  // velocity block we will still have all those associated with subdomain
  // zero before all velocities associated with subdomain one, etc. This is
  // important since we store each of the blocks of this matrix distributed
  // across all processors and want this to be done in such a way that each
  // processor stores that part of the matrix that is roughly equal to the
  // degrees of freedom located on those cells that it will actually work on.
  //
  // When printing the numbers of degrees of freedom, note that these numbers
  // are going to be large if we use many processors. Consequently, we let the
  // stream put a comma separator in between every three digits. The state of
  // the stream, using the locale, is saved from before to after this
  // operation. While slightly opaque, the code works because the default
  // locale (which we get using the constructor call
  // <code>std::locale("")</code>) implies printing numbers with a comma
  // separator for every third digit (i.e., thousands, millions, billions).
  template <int dim>
  void BoussinesqFlowProblem<dim>::setup_dofs ()
  {
    computing_timer.enter_section("Setup dof systems");

    std::vector<unsigned int> stokes_sub_blocks (dim+1,0);
    stokes_sub_blocks[dim] = 1;
    stokes_dof_handler.distribute_dofs (stokes_fe);
    DoFRenumbering::component_wise (stokes_dof_handler, stokes_sub_blocks);

    temperature_dof_handler.distribute_dofs (temperature_fe);

    std::vector<types::global_dof_index> stokes_dofs_per_block (2);
    DoFTools::count_dofs_per_block (stokes_dof_handler, stokes_dofs_per_block,
                                    stokes_sub_blocks);

    const unsigned int n_u = stokes_dofs_per_block[0],
                       n_p = stokes_dofs_per_block[1],
                       n_T = temperature_dof_handler.n_dofs();

    std::locale s = pcout.get_stream().getloc();
    pcout.get_stream().imbue(std::locale(""));
    pcout << "Number of active cells: "
          << triangulation.n_global_active_cells()
          << " (on "
          << triangulation.n_levels()
          << " levels)"
          << std::endl
          << "Number of degrees of freedom: "
          << n_u + n_p + n_T
          << " (" << n_u << '+' << n_p << '+'<< n_T <<')'
          << std::endl
          << std::endl;
    pcout.get_stream().imbue(s);


    // After this, we have to set up the various partitioners (of type
    // <code>IndexSet</code>, see the introduction) that describe which parts
    // of each matrix or vector will be stored where, then call the functions
    // that actually set up the matrices, and at the end also resize the
    // various vectors we keep around in this program.
    std::vector<IndexSet> stokes_partitioning, stokes_relevant_partitioning;
    IndexSet temperature_partitioning (n_T), temperature_relevant_partitioning (n_T);
    IndexSet stokes_relevant_set;
    {
      IndexSet stokes_index_set = stokes_dof_handler.locally_owned_dofs();
      stokes_partitioning.push_back(stokes_index_set.get_view(0,n_u));
      stokes_partitioning.push_back(stokes_index_set.get_view(n_u,n_u+n_p));

      DoFTools::extract_locally_relevant_dofs (stokes_dof_handler,
                                               stokes_relevant_set);
      stokes_relevant_partitioning.push_back(stokes_relevant_set.get_view(0,n_u));
      stokes_relevant_partitioning.push_back(stokes_relevant_set.get_view(n_u,n_u+n_p));

      temperature_partitioning = temperature_dof_handler.locally_owned_dofs();
      DoFTools::extract_locally_relevant_dofs (temperature_dof_handler,
                                               temperature_relevant_partitioning);
    }

    // Following this, we can compute constraints for the solution vectors,
    // including hanging node constraints and homogeneous and inhomogeneous
    // boundary values for the Stokes and temperature fields. Note that as for
    // everything else, the constraint objects can not hold <i>all</i>
    // constraints on every processor. Rather, each processor needs to store
    // only those that are actually necessary for correctness given that it
    // only assembles linear systems on cells it owns. As discussed in the
    // @ref distributed_paper "this paper", the set of constraints we need to
    // know about is exactly the set of constraints on all locally relevant
    // degrees of freedom, so this is what we use to initialize the constraint
    // objects.
    {
      stokes_constraints.clear ();
      stokes_constraints.reinit (stokes_relevant_set);

      DoFTools::make_hanging_node_constraints (stokes_dof_handler,
                                               stokes_constraints);

      FEValuesExtractors::Vector velocity_components(0);
      VectorTools::interpolate_boundary_values (stokes_dof_handler,
                                                0,
                                                ZeroFunction<dim>(dim+1),
                                                stokes_constraints,
                                                stokes_fe.component_mask(velocity_components));

      std::set<types::boundary_id> no_normal_flux_boundaries;
      no_normal_flux_boundaries.insert (1);
      VectorTools::compute_no_normal_flux_constraints (stokes_dof_handler, 0,
                                                       no_normal_flux_boundaries,
                                                       stokes_constraints,
                                                       mapping);
      stokes_constraints.close ();
    }
    {
      temperature_constraints.clear ();
      temperature_constraints.reinit (temperature_relevant_partitioning);

      DoFTools::make_hanging_node_constraints (temperature_dof_handler,
                                               temperature_constraints);
      VectorTools::interpolate_boundary_values (temperature_dof_handler,
                                                0,
                                                EquationData::TemperatureInitialValues<dim>(),
                                                temperature_constraints);
      VectorTools::interpolate_boundary_values (temperature_dof_handler,
                                                1,
                                                EquationData::TemperatureInitialValues<dim>(),
                                                temperature_constraints);
      temperature_constraints.close ();
    }

    // All this done, we can then initialize the various matrix and vector
    // objects to their proper sizes. At the end, we also record that all
    // matrices and preconditioners have to be re-computed at the beginning of
    // the next time step. Note how we initialize the vectors for the Stokes
    // and temperature right hand sides: These are writable vectors (last
    // boolean argument set to @p true) that have the correct one-to-one
    // partitioning of locally owned elements but are still given the relevant
    // partitioning for means of figuring out the vector entries that are
    // going to be set right away. As for matrices, this allows for writing
    // local contributions into the vector with multiple threads (always
    // assuming that the same vector entry is not accessed by multiple threads
    // at the same time). The other vectors only allow for read access of
    // individual elements, including ghosts, but are not suitable for
    // solvers.
    setup_stokes_matrix (stokes_partitioning, stokes_relevant_partitioning);
    setup_stokes_preconditioner (stokes_partitioning,
                                 stokes_relevant_partitioning);
    setup_temperature_matrices (temperature_partitioning,
                                temperature_relevant_partitioning);

    stokes_rhs.reinit (stokes_partitioning, stokes_relevant_partitioning,
                       MPI_COMM_WORLD, true);
    stokes_solution.reinit (stokes_relevant_partitioning, MPI_COMM_WORLD);
    old_stokes_solution.reinit (stokes_solution);

    temperature_rhs.reinit (temperature_partitioning,
                            temperature_relevant_partitioning,
                            MPI_COMM_WORLD, true);
    temperature_solution.reinit (temperature_relevant_partitioning, MPI_COMM_WORLD);
    old_temperature_solution.reinit (temperature_solution);
    old_old_temperature_solution.reinit (temperature_solution);

    rebuild_stokes_matrix              = true;
    rebuild_stokes_preconditioner      = true;
    rebuild_temperature_matrices       = true;
    rebuild_temperature_preconditioner = true;

    computing_timer.exit_section();
  }



  // @sect4{The BoussinesqFlowProblem assembly functions}
  //
  // Following the discussion in the introduction and in the @ref threads
  // module, we split the assembly functions into different parts:
  //
  // <ul> <li> The local calculations of matrices and right hand sides, given
  // a certain cell as input (these functions are named
  // <code>local_assemble_*</code> below). The resulting function is, in other
  // words, essentially the body of the loop over all cells in step-31. Note,
  // however, that these functions store the result from the local
  // calculations in variables of classes from the CopyData namespace.
  //
  // <li>These objects are then given to the second step which writes the
  // local data into the global data structures (these functions are named
  // <code>copy_local_to_global_*</code> below). These functions are pretty
  // trivial.
  //
  // <li>These two subfunctions are then used in the respective assembly
  // routine (called <code>assemble_*</code> below), where a WorkStream object
  // is set up and runs over all the cells that belong to the processor's
  // subdomain.  </ul>

  // @sect5{Stokes preconditioner assembly}
  //
  // Let us start with the functions that builds the Stokes
  // preconditioner. The first two of these are pretty trivial, given the
  // discussion above. Note in particular that the main point in using the
  // scratch data object is that we want to avoid allocating any objects on
  // the free space each time we visit a new cell. As a consequence, the
  // assembly function below only has automatic local variables, and
  // everything else is accessed through the scratch data object, which is
  // allocated only once before we start the loop over all cells:
  template <int dim>
  void
  BoussinesqFlowProblem<dim>::
  local_assemble_stokes_preconditioner (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                        Assembly::Scratch::StokesPreconditioner<dim> &scratch,
                                        Assembly::CopyData::StokesPreconditioner<dim> &data)
  {
    const unsigned int   dofs_per_cell   = stokes_fe.dofs_per_cell;
    const unsigned int   n_q_points      = scratch.stokes_fe_values.n_quadrature_points;

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    scratch.stokes_fe_values.reinit (cell);
    cell->get_dof_indices (data.local_dof_indices);

    data.local_matrix = 0;

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grad_phi_u[k] = scratch.stokes_fe_values[velocities].gradient(k,q);
            scratch.phi_p[k]      = scratch.stokes_fe_values[pressure].value (k, q);
          }

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            data.local_matrix(i,j) += (EquationData::eta *
                                       scalar_product (scratch.grad_phi_u[i],
                                                       scratch.grad_phi_u[j])
                                       +
                                       (1./EquationData::eta) *
                                       EquationData::pressure_scaling *
                                       EquationData::pressure_scaling *
                                       (scratch.phi_p[i] * scratch.phi_p[j]))
                                      * scratch.stokes_fe_values.JxW(q);
      }
  }



  template <int dim>
  void
  BoussinesqFlowProblem<dim>::
  copy_local_to_global_stokes_preconditioner (const Assembly::CopyData::StokesPreconditioner<dim> &data)
  {
    stokes_constraints.distribute_local_to_global (data.local_matrix,
                                                   data.local_dof_indices,
                                                   stokes_preconditioner_matrix);
  }


  // Now for the function that actually puts things together, using the
  // WorkStream functions.  WorkStream::run needs a start and end iterator to
  // enumerate the cells it is supposed to work on. Typically, one would use
  // DoFHandler::begin_active() and DoFHandler::end() for that but here we
  // actually only want the subset of cells that in fact are owned by the
  // current processor. This is where the FilteredIterator class comes into
  // play: you give it a range of cells and it provides an iterator that only
  // iterates over that subset of cells that satisfy a certain predicate (a
  // predicate is a function of one argument that either returns true or
  // false). The predicate we use here is IteratorFilters::LocallyOwnedCell,
  // i.e., it returns true exactly if the cell is owned by the current
  // processor. The resulting iterator range is then exactly what we need.
  //
  // With this obstacle out of the way, we call the WorkStream::run
  // function with this set of cells, scratch and copy objects, and
  // with pointers to two functions: the local assembly and
  // copy-local-to-global function. These functions need to have very
  // specific signatures: three arguments in the first and one
  // argument in the latter case (see the documentation of the
  // WorkStream::run function for the meaning of these arguments).
  // Note how we use the construct <code>std_cxx11::bind</code> to
  // create a function object that satisfies this requirement. It uses
  // placeholders <code>std_cxx11::_1, std_cxx11::_2,
  // std_cxx11::_3</code> for the local assembly function that specify
  // cell, scratch data, and copy data, as well as the placeholder
  // <code>std_cxx11::_1</code> for the copy function that expects the
  // data to be written into the global matrix (for placeholder
  // arguments, also see the discussion in step-13's
  // <code>assemble_linear_system()</code> function). On the other
  // hand, the implicit zeroth argument of member functions (namely
  // the <code>this</code> pointer of the object on which that member
  // function is to operate on) is <i>bound</i> to the
  // <code>this</code> pointer of the current function. The
  // WorkStream::run function, as a consequence, does not need to know
  // anything about the object these functions work on.
  //
  // When the WorkStream is executed, it will create several local assembly
  // routines of the first kind for several cells and let some available
  // processors work on them. The function that needs to be synchronized,
  // i.e., the write operation into the global matrix, however, is executed by
  // only one thread at a time in the prescribed order. Of course, this only
  // holds for the parallelization on a single MPI process. Different MPI
  // processes will have their own WorkStream objects and do that work
  // completely independently (and in different memory spaces). In a
  // distributed calculation, some data will accumulate at degrees of freedom
  // that are not owned by the respective processor. It would be inefficient
  // to send data around every time we encounter such a dof. What happens
  // instead is that the Trilinos sparse matrix will keep that data and send
  // it to the owner at the end of assembly, by calling the
  // <code>compress()</code> command.
  template <int dim>
  void
  BoussinesqFlowProblem<dim>::assemble_stokes_preconditioner ()
  {
    stokes_preconditioner_matrix = 0;

    const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree+1);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.end()),
         std_cxx11::bind (&BoussinesqFlowProblem<dim>::
                          local_assemble_stokes_preconditioner,
                          this,
                          std_cxx11::_1,
                          std_cxx11::_2,
                          std_cxx11::_3),
         std_cxx11::bind (&BoussinesqFlowProblem<dim>::
                          copy_local_to_global_stokes_preconditioner,
                          this,
                          std_cxx11::_1),
         Assembly::Scratch::
         StokesPreconditioner<dim> (stokes_fe, quadrature_formula,
                                    mapping,
                                    update_JxW_values |
                                    update_values |
                                    update_gradients),
         Assembly::CopyData::
         StokesPreconditioner<dim> (stokes_fe));

    stokes_preconditioner_matrix.compress(VectorOperation::add);
  }



  // The final function in this block initiates assembly of the Stokes
  // preconditioner matrix and then in fact builds the Stokes
  // preconditioner. It is mostly the same as in the serial case. The only
  // difference to step-31 is that we use a Jacobi preconditioner for the
  // pressure mass matrix instead of IC, as discussed in the introduction.
  template <int dim>
  void
  BoussinesqFlowProblem<dim>::build_stokes_preconditioner ()
  {
    if (rebuild_stokes_preconditioner == false)
      return;

    computing_timer.enter_section ("   Build Stokes preconditioner");
    pcout << "   Rebuilding Stokes preconditioner..." << std::flush;

    assemble_stokes_preconditioner ();

    std::vector<std::vector<bool> > constant_modes;
    FEValuesExtractors::Vector velocity_components(0);
    DoFTools::extract_constant_modes (stokes_dof_handler,
                                      stokes_fe.component_mask(velocity_components),
                                      constant_modes);

    Mp_preconditioner.reset  (new TrilinosWrappers::PreconditionJacobi());
    Amg_preconditioner.reset (new TrilinosWrappers::PreconditionAMG());

    TrilinosWrappers::PreconditionAMG::AdditionalData Amg_data;
    Amg_data.constant_modes = constant_modes;
    Amg_data.elliptic = true;
    Amg_data.higher_order_elements = true;
    Amg_data.smoother_sweeps = 2;
    Amg_data.aggregation_threshold = 0.02;

    Mp_preconditioner->initialize (stokes_preconditioner_matrix.block(1,1));
    Amg_preconditioner->initialize (stokes_preconditioner_matrix.block(0,0),
                                    Amg_data);

    rebuild_stokes_preconditioner = false;

    pcout << std::endl;
    computing_timer.exit_section();
  }


  // @sect5{Stokes system assembly}

  // The next three functions implement the assembly of the Stokes system,
  // again split up into a part performing local calculations, one for writing
  // the local data into the global matrix and vector, and one for actually
  // running the loop over all cells with the help of the WorkStream
  // class. Note that the assembly of the Stokes matrix needs only to be done
  // in case we have changed the mesh. Otherwise, just the
  // (temperature-dependent) right hand side needs to be calculated
  // here. Since we are working with distributed matrices and vectors, we have
  // to call the respective <code>compress()</code> functions in the end of
  // the assembly in order to send non-local data to the owner process.
  template <int dim>
  void
  BoussinesqFlowProblem<dim>::
  local_assemble_stokes_system (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                Assembly::Scratch::StokesSystem<dim> &scratch,
                                Assembly::CopyData::StokesSystem<dim> &data)
  {
    const unsigned int dofs_per_cell = scratch.stokes_fe_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.stokes_fe_values.n_quadrature_points;

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    scratch.stokes_fe_values.reinit (cell);

    typename DoFHandler<dim>::active_cell_iterator
    temperature_cell (&triangulation,
                      cell->level(),
                      cell->index(),
                      &temperature_dof_handler);
    scratch.temperature_fe_values.reinit (temperature_cell);

    if (rebuild_stokes_matrix)
      data.local_matrix = 0;
    data.local_rhs = 0;

    scratch.temperature_fe_values.get_function_values (old_temperature_solution,
                                                       scratch.old_temperature_values);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        const double old_temperature = scratch.old_temperature_values[q];

        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.phi_u[k] = scratch.stokes_fe_values[velocities].value (k,q);
            if (rebuild_stokes_matrix)
              {
                scratch.grads_phi_u[k] = scratch.stokes_fe_values[velocities].symmetric_gradient(k,q);
                scratch.div_phi_u[k]   = scratch.stokes_fe_values[velocities].divergence (k, q);
                scratch.phi_p[k]       = scratch.stokes_fe_values[pressure].value (k, q);
              }
          }

        if (rebuild_stokes_matrix == true)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              data.local_matrix(i,j) += (EquationData::eta * 2 *
                                         (scratch.grads_phi_u[i] * scratch.grads_phi_u[j])
                                         - (EquationData::pressure_scaling *
                                            scratch.div_phi_u[i] * scratch.phi_p[j])
                                         - (EquationData::pressure_scaling *
                                            scratch.phi_p[i] * scratch.div_phi_u[j]))
                                        * scratch.stokes_fe_values.JxW(q);

        const Tensor<1,dim>
        gravity = EquationData::gravity_vector (scratch.stokes_fe_values
                                                .quadrature_point(q));

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          data.local_rhs(i) += (EquationData::density(old_temperature) *
                                gravity  *
                                scratch.phi_u[i]) *
                               scratch.stokes_fe_values.JxW(q);
      }

    cell->get_dof_indices (data.local_dof_indices);
  }



  template <int dim>
  void
  BoussinesqFlowProblem<dim>::
  copy_local_to_global_stokes_system (const Assembly::CopyData::StokesSystem<dim> &data)
  {
    if (rebuild_stokes_matrix == true)
      stokes_constraints.distribute_local_to_global (data.local_matrix,
                                                     data.local_rhs,
                                                     data.local_dof_indices,
                                                     stokes_matrix,
                                                     stokes_rhs);
    else
      stokes_constraints.distribute_local_to_global (data.local_rhs,
                                                     data.local_dof_indices,
                                                     stokes_rhs);
  }



  template <int dim>
  void BoussinesqFlowProblem<dim>::assemble_stokes_system ()
  {
    computing_timer.enter_section ("   Assemble Stokes system");

    if (rebuild_stokes_matrix == true)
      stokes_matrix=0;

    stokes_rhs=0;

    const QGauss<dim> quadrature_formula(parameters.stokes_velocity_degree+1);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     stokes_dof_handler.end()),
         std_cxx11::bind (&BoussinesqFlowProblem<dim>::
                          local_assemble_stokes_system,
                          this,
                          std_cxx11::_1,
                          std_cxx11::_2,
                          std_cxx11::_3),
         std_cxx11::bind (&BoussinesqFlowProblem<dim>::
                          copy_local_to_global_stokes_system,
                          this,
                          std_cxx11::_1),
         Assembly::Scratch::
         StokesSystem<dim> (stokes_fe, mapping, quadrature_formula,
                            (update_values    |
                             update_quadrature_points  |
                             update_JxW_values |
                             (rebuild_stokes_matrix == true
                              ?
                              update_gradients
                              :
                              UpdateFlags(0))),
                            temperature_fe,
                            update_values),
         Assembly::CopyData::
         StokesSystem<dim> (stokes_fe));

    if (rebuild_stokes_matrix == true)
      stokes_matrix.compress(VectorOperation::add);
    stokes_rhs.compress(VectorOperation::add);

    rebuild_stokes_matrix = false;

    pcout << std::endl;
    computing_timer.exit_section();
  }


  // @sect5{Temperature matrix assembly}

  // The task to be performed by the next three functions is to calculate a
  // mass matrix and a Laplace matrix on the temperature system. These will be
  // combined in order to yield the semi-implicit time stepping matrix that
  // consists of the mass matrix plus a time step-dependent weight factor
  // times the Laplace matrix. This function is again essentially the body of
  // the loop over all cells from step-31.
  //
  // The two following functions perform similar services as the ones above.
  template <int dim>
  void BoussinesqFlowProblem<dim>::
  local_assemble_temperature_matrix (const typename DoFHandler<dim>::active_cell_iterator &cell,
                                     Assembly::Scratch::TemperatureMatrix<dim> &scratch,
                                     Assembly::CopyData::TemperatureMatrix<dim> &data)
  {
    const unsigned int dofs_per_cell = scratch.temperature_fe_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.temperature_fe_values.n_quadrature_points;

    scratch.temperature_fe_values.reinit (cell);
    cell->get_dof_indices (data.local_dof_indices);

    data.local_mass_matrix = 0;
    data.local_stiffness_matrix = 0;

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.grad_phi_T[k] = scratch.temperature_fe_values.shape_grad (k,q);
            scratch.phi_T[k]      = scratch.temperature_fe_values.shape_value (k, q);
          }

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          for (unsigned int j=0; j<dofs_per_cell; ++j)
            {
              data.local_mass_matrix(i,j)
              += (scratch.phi_T[i] * scratch.phi_T[j]
                  *
                  scratch.temperature_fe_values.JxW(q));
              data.local_stiffness_matrix(i,j)
              += (EquationData::kappa * scratch.grad_phi_T[i] * scratch.grad_phi_T[j]
                  *
                  scratch.temperature_fe_values.JxW(q));
            }
      }
  }



  template <int dim>
  void
  BoussinesqFlowProblem<dim>::
  copy_local_to_global_temperature_matrix (const Assembly::CopyData::TemperatureMatrix<dim> &data)
  {
    temperature_constraints.distribute_local_to_global (data.local_mass_matrix,
                                                        data.local_dof_indices,
                                                        temperature_mass_matrix);
    temperature_constraints.distribute_local_to_global (data.local_stiffness_matrix,
                                                        data.local_dof_indices,
                                                        temperature_stiffness_matrix);
  }


  template <int dim>
  void BoussinesqFlowProblem<dim>::assemble_temperature_matrix ()
  {
    if (rebuild_temperature_matrices == false)
      return;

    computing_timer.enter_section ("   Assemble temperature matrices");
    temperature_mass_matrix = 0;
    temperature_stiffness_matrix = 0;

    const QGauss<dim> quadrature_formula(parameters.temperature_degree+2);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     temperature_dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     temperature_dof_handler.end()),
         std_cxx11::bind (&BoussinesqFlowProblem<dim>::
                          local_assemble_temperature_matrix,
                          this,
                          std_cxx11::_1,
                          std_cxx11::_2,
                          std_cxx11::_3),
         std_cxx11::bind (&BoussinesqFlowProblem<dim>::
                          copy_local_to_global_temperature_matrix,
                          this,
                          std_cxx11::_1),
         Assembly::Scratch::
         TemperatureMatrix<dim> (temperature_fe, mapping, quadrature_formula),
         Assembly::CopyData::
         TemperatureMatrix<dim> (temperature_fe));

    temperature_mass_matrix.compress(VectorOperation::add);
    temperature_stiffness_matrix.compress(VectorOperation::add);

    rebuild_temperature_matrices = false;
    rebuild_temperature_preconditioner = true;

    computing_timer.exit_section();
  }


  // @sect5{Temperature right hand side assembly}

  // This is the last assembly function. It calculates the right hand side of
  // the temperature system, which includes the convection and the
  // stabilization terms. It includes a lot of evaluations of old solutions at
  // the quadrature points (which are necessary for calculating the artificial
  // viscosity of stabilization), but is otherwise similar to the other
  // assembly functions. Notice, once again, how we resolve the dilemma of
  // having inhomogeneous boundary conditions, by just making a right hand
  // side at this point (compare the comments for the <code>project()</code>
  // function above): We create some matrix columns with exactly the values
  // that would be entered for the temperature stiffness matrix, in case we
  // have inhomogeneously constrained dofs. That will account for the correct
  // balance of the right hand side vector with the matrix system of
  // temperature.
  template <int dim>
  void BoussinesqFlowProblem<dim>::
  local_assemble_temperature_rhs (const std::pair<double,double> global_T_range,
                                  const double                   global_max_velocity,
                                  const double                   global_entropy_variation,
                                  const typename DoFHandler<dim>::active_cell_iterator &cell,
                                  Assembly::Scratch::TemperatureRHS<dim> &scratch,
                                  Assembly::CopyData::TemperatureRHS<dim> &data)
  {
    const bool use_bdf2_scheme = (timestep_number != 0);

    const unsigned int dofs_per_cell = scratch.temperature_fe_values.get_fe().dofs_per_cell;
    const unsigned int n_q_points    = scratch.temperature_fe_values.n_quadrature_points;

    const FEValuesExtractors::Vector velocities (0);

    data.local_rhs = 0;
    data.matrix_for_bc = 0;
    cell->get_dof_indices (data.local_dof_indices);

    scratch.temperature_fe_values.reinit (cell);

    typename DoFHandler<dim>::active_cell_iterator
    stokes_cell (&triangulation,
                 cell->level(),
                 cell->index(),
                 &stokes_dof_handler);
    scratch.stokes_fe_values.reinit (stokes_cell);

    scratch.temperature_fe_values.get_function_values (old_temperature_solution,
                                                       scratch.old_temperature_values);
    scratch.temperature_fe_values.get_function_values (old_old_temperature_solution,
                                                       scratch.old_old_temperature_values);

    scratch.temperature_fe_values.get_function_gradients (old_temperature_solution,
                                                          scratch.old_temperature_grads);
    scratch.temperature_fe_values.get_function_gradients (old_old_temperature_solution,
                                                          scratch.old_old_temperature_grads);

    scratch.temperature_fe_values.get_function_laplacians (old_temperature_solution,
                                                           scratch.old_temperature_laplacians);
    scratch.temperature_fe_values.get_function_laplacians (old_old_temperature_solution,
                                                           scratch.old_old_temperature_laplacians);

    scratch.stokes_fe_values[velocities].get_function_values (stokes_solution,
                                                              scratch.old_velocity_values);
    scratch.stokes_fe_values[velocities].get_function_values (old_stokes_solution,
                                                              scratch.old_old_velocity_values);
    scratch.stokes_fe_values[velocities].get_function_symmetric_gradients (stokes_solution,
        scratch.old_strain_rates);
    scratch.stokes_fe_values[velocities].get_function_symmetric_gradients (old_stokes_solution,
        scratch.old_old_strain_rates);

    const double nu
      = compute_viscosity (scratch.old_temperature_values,
                           scratch.old_old_temperature_values,
                           scratch.old_temperature_grads,
                           scratch.old_old_temperature_grads,
                           scratch.old_temperature_laplacians,
                           scratch.old_old_temperature_laplacians,
                           scratch.old_velocity_values,
                           scratch.old_old_velocity_values,
                           scratch.old_strain_rates,
                           scratch.old_old_strain_rates,
                           global_max_velocity,
                           global_T_range.second - global_T_range.first,
                           0.5 * (global_T_range.second + global_T_range.first),
                           global_entropy_variation,
                           cell->diameter());

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        for (unsigned int k=0; k<dofs_per_cell; ++k)
          {
            scratch.phi_T[k]      = scratch.temperature_fe_values.shape_value (k, q);
            scratch.grad_phi_T[k] = scratch.temperature_fe_values.shape_grad (k,q);
          }


        const double T_term_for_rhs
          = (use_bdf2_scheme ?
             (scratch.old_temperature_values[q] *
              (1 + time_step/old_time_step)
              -
              scratch.old_old_temperature_values[q] *
              (time_step * time_step) /
              (old_time_step * (time_step + old_time_step)))
             :
             scratch.old_temperature_values[q]);

        const double ext_T
          = (use_bdf2_scheme ?
             (scratch.old_temperature_values[q] *
              (1 + time_step/old_time_step)
              -
              scratch.old_old_temperature_values[q] *
              time_step/old_time_step)
             :
             scratch.old_temperature_values[q]);

        const Tensor<1,dim> ext_grad_T
          = (use_bdf2_scheme ?
             (scratch.old_temperature_grads[q] *
              (1 + time_step/old_time_step)
              -
              scratch.old_old_temperature_grads[q] *
              time_step/old_time_step)
             :
             scratch.old_temperature_grads[q]);

        const Tensor<1,dim> extrapolated_u
          = (use_bdf2_scheme ?
             (scratch.old_velocity_values[q] *
              (1 + time_step/old_time_step)
              -
              scratch.old_old_velocity_values[q] *
              time_step/old_time_step)
             :
             scratch.old_velocity_values[q]);

        const SymmetricTensor<2,dim> extrapolated_strain_rate
          = (use_bdf2_scheme ?
             (scratch.old_strain_rates[q] *
              (1 + time_step/old_time_step)
              -
              scratch.old_old_strain_rates[q] *
              time_step/old_time_step)
             :
             scratch.old_strain_rates[q]);

        const double gamma
          = ((EquationData::radiogenic_heating * EquationData::density(ext_T)
              +
              2 * EquationData::eta * extrapolated_strain_rate * extrapolated_strain_rate) /
             (EquationData::density(ext_T) * EquationData::specific_heat));

        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            data.local_rhs(i) += (T_term_for_rhs * scratch.phi_T[i]
                                  -
                                  time_step *
                                  extrapolated_u * ext_grad_T * scratch.phi_T[i]
                                  -
                                  time_step *
                                  nu * ext_grad_T * scratch.grad_phi_T[i]
                                  +
                                  time_step *
                                  gamma * scratch.phi_T[i])
                                 *
                                 scratch.temperature_fe_values.JxW(q);

            if (temperature_constraints.is_inhomogeneously_constrained(data.local_dof_indices[i]))
              {
                for (unsigned int j=0; j<dofs_per_cell; ++j)
                  data.matrix_for_bc(j,i) += (scratch.phi_T[i] * scratch.phi_T[j] *
                                              (use_bdf2_scheme ?
                                               ((2*time_step + old_time_step) /
                                                (time_step + old_time_step)) : 1.)
                                              +
                                              scratch.grad_phi_T[i] *
                                              scratch.grad_phi_T[j] *
                                              EquationData::kappa *
                                              time_step)
                                             *
                                             scratch.temperature_fe_values.JxW(q);
              }
          }
      }
  }


  template <int dim>
  void
  BoussinesqFlowProblem<dim>::
  copy_local_to_global_temperature_rhs (const Assembly::CopyData::TemperatureRHS<dim> &data)
  {
    temperature_constraints.distribute_local_to_global (data.local_rhs,
                                                        data.local_dof_indices,
                                                        temperature_rhs,
                                                        data.matrix_for_bc);
  }



  // In the function that runs the WorkStream for actually calculating the
  // right hand side, we also generate the final matrix. As mentioned above,
  // it is a sum of the mass matrix and the Laplace matrix, times some time
  // step-dependent weight. This weight is specified by the BDF-2 time
  // integration scheme, see the introduction in step-31. What is new in this
  // tutorial program (in addition to the use of MPI parallelization and the
  // WorkStream class), is that we now precompute the temperature
  // preconditioner as well. The reason is that the setup of the Jacobi
  // preconditioner takes a noticeable time compared to the solver because we
  // usually only need between 10 and 20 iterations for solving the
  // temperature system (this might sound strange, as Jacobi really only
  // consists of a diagonal, but in Trilinos it is derived from more general
  // framework for point relaxation preconditioners which is a bit
  // inefficient). Hence, it is more efficient to precompute the
  // preconditioner, even though the matrix entries may slightly change
  // because the time step might change. This is not too big a problem because
  // we remesh every few time steps (and regenerate the preconditioner then).
  template <int dim>
  void BoussinesqFlowProblem<dim>::assemble_temperature_system (const double maximal_velocity)
  {
    const bool use_bdf2_scheme = (timestep_number != 0);

    if (use_bdf2_scheme == true)
      {
        temperature_matrix.copy_from (temperature_mass_matrix);
        temperature_matrix *= (2*time_step + old_time_step) /
                              (time_step + old_time_step);
        temperature_matrix.add (time_step, temperature_stiffness_matrix);
      }
    else
      {
        temperature_matrix.copy_from (temperature_mass_matrix);
        temperature_matrix.add (time_step, temperature_stiffness_matrix);
      }

    if (rebuild_temperature_preconditioner == true)
      {
        T_preconditioner.reset (new TrilinosWrappers::PreconditionJacobi());
        T_preconditioner->initialize (temperature_matrix);
        rebuild_temperature_preconditioner = false;
      }

    // The next part is computing the right hand side vectors.  To do so, we
    // first compute the average temperature $T_m$ that we use for evaluating
    // the artificial viscosity stabilization through the residual $E(T) =
    // (T-T_m)^2$. We do this by defining the midpoint between maximum and
    // minimum temperature as average temperature in the definition of the
    // entropy viscosity. An alternative would be to use the integral average,
    // but the results are not very sensitive to this choice. The rest then
    // only requires calling WorkStream::run again, binding the arguments to
    // the <code>local_assemble_temperature_rhs</code> function that are the
    // same in every call to the correct values:
    temperature_rhs = 0;

    const QGauss<dim> quadrature_formula(parameters.temperature_degree+2);
    const std::pair<double,double>
    global_T_range = get_extrapolated_temperature_range();

    const double average_temperature = 0.5 * (global_T_range.first +
                                              global_T_range.second);
    const double global_entropy_variation =
      get_entropy_variation (average_temperature);

    typedef
    FilteredIterator<typename DoFHandler<dim>::active_cell_iterator>
    CellFilter;

    WorkStream::
    run (CellFilter (IteratorFilters::LocallyOwnedCell(),
                     temperature_dof_handler.begin_active()),
         CellFilter (IteratorFilters::LocallyOwnedCell(),
                     temperature_dof_handler.end()),
         std_cxx11::bind (&BoussinesqFlowProblem<dim>::
                          local_assemble_temperature_rhs,
                          this,
                          global_T_range,
                          maximal_velocity,
                          global_entropy_variation,
                          std_cxx11::_1,
                          std_cxx11::_2,
                          std_cxx11::_3),
         std_cxx11::bind (&BoussinesqFlowProblem<dim>::
                          copy_local_to_global_temperature_rhs,
                          this,
                          std_cxx11::_1),
         Assembly::Scratch::
         TemperatureRHS<dim> (temperature_fe, stokes_fe, mapping,
                              quadrature_formula),
         Assembly::CopyData::
         TemperatureRHS<dim> (temperature_fe));

    temperature_rhs.compress(VectorOperation::add);
  }




  // @sect4{BoussinesqFlowProblem::solve}

  // This function solves the linear systems in each time step of the
  // Boussinesq problem. First, we work on the Stokes system and then on the
  // temperature system. In essence, it does the same things as the respective
  // function in step-31. However, there are a few changes here.
  //
  // The first change is related to the way we store our solution: we keep the
  // vectors with locally owned degrees of freedom plus ghost nodes on each
  // MPI node. When we enter a solver which is supposed to perform
  // matrix-vector products with a distributed matrix, this is not the
  // appropriate form, though. There, we will want to have the solution vector
  // to be distributed in the same way as the matrix, i.e. without any
  // ghosts. So what we do first is to generate a distributed vector called
  // <code>distributed_stokes_solution</code> and put only the locally owned
  // dofs into that, which is neatly done by the <code>operator=</code> of the
  // Trilinos vector.
  //
  // Next, we scale the pressure solution (or rather, the initial guess) for
  // the solver so that it matches with the length scales in the matrices, as
  // discussed in the introduction. We also immediately scale the pressure
  // solution back to the correct units after the solution is completed.  We
  // also need to set the pressure values at hanging nodes to zero. This we
  // also did in step-31 in order not to disturb the Schur complement by some
  // vector entries that actually are irrelevant during the solve stage. As a
  // difference to step-31, here we do it only for the locally owned pressure
  // dofs. After solving for the Stokes solution, each processor copies the
  // distributed solution back into the solution vector that also includes
  // ghost elements.
  //
  // The third and most obvious change is that we have two variants for the
  // Stokes solver: A fast solver that sometimes breaks down, and a robust
  // solver that is slower. This is what we already discussed in the
  // introduction. Here is how we realize it: First, we perform 30 iterations
  // with the fast solver based on the simple preconditioner based on the AMG
  // V-cycle instead of an approximate solve (this is indicated by the
  // <code>false</code> argument to the
  // <code>LinearSolvers::BlockSchurPreconditioner</code> object). If we
  // converge, everything is fine. If we do not converge, the solver control
  // object will throw an exception SolverControl::NoConvergence. Usually,
  // this would abort the program because we don't catch them in our usual
  // <code>solve()</code> functions. This is certainly not what we want to
  // happen here. Rather, we want to switch to the strong solver and continue
  // the solution process with whatever vector we got so far. Hence, we catch
  // the exception with the C++ try/catch mechanism. We then simply go through
  // the same solver sequence again in the <code>catch</code> clause, this
  // time passing the @p true flag to the preconditioner for the strong
  // solver, signaling an approximate CG solve.
  template <int dim>
  void BoussinesqFlowProblem<dim>::solve ()
  {
    computing_timer.enter_section ("   Solve Stokes system");

    {
      pcout << "   Solving Stokes system... " << std::flush;

      TrilinosWrappers::MPI::BlockVector
      distributed_stokes_solution (stokes_rhs);
      distributed_stokes_solution = stokes_solution;

      distributed_stokes_solution.block(1) /= EquationData::pressure_scaling;

      const unsigned int
      start = (distributed_stokes_solution.block(0).size() +
               distributed_stokes_solution.block(1).local_range().first),
              end   = (distributed_stokes_solution.block(0).size() +
                       distributed_stokes_solution.block(1).local_range().second);
      for (unsigned int i=start; i<end; ++i)
        if (stokes_constraints.is_constrained (i))
          distributed_stokes_solution(i) = 0;


      PrimitiveVectorMemory<TrilinosWrappers::MPI::BlockVector> mem;

      unsigned int n_iterations = 0;
      const double solver_tolerance = 1e-8 * stokes_rhs.l2_norm();
      SolverControl solver_control (30, solver_tolerance);

      try
        {
          const LinearSolvers::BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG,
                TrilinosWrappers::PreconditionJacobi>
                preconditioner (stokes_matrix, stokes_preconditioner_matrix,
                                *Mp_preconditioner, *Amg_preconditioner,
                                false);

          SolverFGMRES<TrilinosWrappers::MPI::BlockVector>
          solver(solver_control, mem,
                 SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::
                 AdditionalData(30, true));
          solver.solve(stokes_matrix, distributed_stokes_solution, stokes_rhs,
                       preconditioner);

          n_iterations = solver_control.last_step();
        }

      catch (SolverControl::NoConvergence)
        {
          const LinearSolvers::BlockSchurPreconditioner<TrilinosWrappers::PreconditionAMG,
                TrilinosWrappers::PreconditionJacobi>
                preconditioner (stokes_matrix, stokes_preconditioner_matrix,
                                *Mp_preconditioner, *Amg_preconditioner,
                                true);

          SolverControl solver_control_refined (stokes_matrix.m(), solver_tolerance);
          SolverFGMRES<TrilinosWrappers::MPI::BlockVector>
          solver(solver_control_refined, mem,
                 SolverFGMRES<TrilinosWrappers::MPI::BlockVector>::
                 AdditionalData(50, true));
          solver.solve(stokes_matrix, distributed_stokes_solution, stokes_rhs,
                       preconditioner);

          n_iterations = (solver_control.last_step() +
                          solver_control_refined.last_step());
        }


      stokes_constraints.distribute (distributed_stokes_solution);

      distributed_stokes_solution.block(1) *= EquationData::pressure_scaling;

      stokes_solution = distributed_stokes_solution;
      pcout << n_iterations  << " iterations."
            << std::endl;
    }
    computing_timer.exit_section();


    // Now let's turn to the temperature part: First, we compute the time step
    // size. We found that we need smaller time steps for 3D than for 2D for
    // the shell geometry. This is because the cells are more distorted in
    // that case (it is the smallest edge length that determines the CFL
    // number). Instead of computing the time step from maximum velocity and
    // minimal mesh size as in step-31, we compute local CFL numbers, i.e., on
    // each cell we compute the maximum velocity times the mesh size, and
    // compute the maximum of them. Hence, we need to choose the factor in
    // front of the time step slightly smaller.
    //
    // After temperature right hand side assembly, we solve the linear system
    // for temperature (with fully distributed vectors without any ghosts),
    // apply constraints and copy the vector back to one with ghosts.
    //
    // In the end, we extract the temperature range similarly to step-31 to
    // produce some output (for example in order to help us choose the
    // stabilization constants, as discussed in the introduction). The only
    // difference is that we need to exchange maxima over all processors.
    computing_timer.enter_section ("   Assemble temperature rhs");
    {
      old_time_step = time_step;

      const double scaling = (dim==3 ? 0.25 : 1.0);
      time_step = (scaling/(2.1*dim*std::sqrt(1.*dim)) /
                   (parameters.temperature_degree *
                    get_cfl_number()));

      const double maximal_velocity = get_maximal_velocity();
      pcout << "   Maximal velocity: "
            << maximal_velocity *EquationData::year_in_seconds * 100
            << " cm/year"
            << std::endl;
      pcout << "   " << "Time step: "
            << time_step/EquationData::year_in_seconds
            << " years"
            << std::endl;

      temperature_solution = old_temperature_solution;
      assemble_temperature_system (maximal_velocity);
    }
    computing_timer.exit_section ();

    computing_timer.enter_section ("   Solve temperature system");
    {
      SolverControl solver_control (temperature_matrix.m(),
                                    1e-12*temperature_rhs.l2_norm());
      SolverCG<TrilinosWrappers::MPI::Vector>   cg (solver_control);

      TrilinosWrappers::MPI::Vector
      distributed_temperature_solution (temperature_rhs);
      distributed_temperature_solution = temperature_solution;

      cg.solve (temperature_matrix, distributed_temperature_solution,
                temperature_rhs, *T_preconditioner);

      temperature_constraints.distribute (distributed_temperature_solution);
      temperature_solution = distributed_temperature_solution;

      pcout << "   "
            << solver_control.last_step()
            << " CG iterations for temperature" << std::endl;
      computing_timer.exit_section();

      double temperature[2] = { std::numeric_limits<double>::max(),
                                -std::numeric_limits<double>::max()
                              };
      double global_temperature[2];

      for (unsigned int i=distributed_temperature_solution.local_range().first;
           i < distributed_temperature_solution.local_range().second; ++i)
        {
          temperature[0] = std::min<double> (temperature[0],
                                             distributed_temperature_solution(i));
          temperature[1] = std::max<double> (temperature[1],
                                             distributed_temperature_solution(i));
        }

      temperature[0] *= -1.0;
      Utilities::MPI::max (temperature, MPI_COMM_WORLD, global_temperature);
      global_temperature[0] *= -1.0;

      pcout << "   Temperature range: "
            << global_temperature[0] << ' ' << global_temperature[1]
            << std::endl;
    }
  }


  // @sect4{BoussinesqFlowProblem::output_results}

  // Next comes the function that generates the output. The quantities to
  // output could be introduced manually like we did in step-31. An
  // alternative is to hand this task over to a class PostProcessor that
  // inherits from the class DataPostprocessor, which can be attached to
  // DataOut. This allows us to output derived quantities from the solution,
  // like the friction heating included in this example. It overloads the
  // virtual function DataPostprocessor::compute_derived_quantities_vector,
  // which is then internally called from DataOut::build_patches. We have to
  // give it values of the numerical solution, its derivatives, normals to the
  // cell, the actual evaluation points and any additional quantities. This
  // follows the same procedure as discussed in step-29 and other programs.
  template <int dim>
  class BoussinesqFlowProblem<dim>::Postprocessor : public DataPostprocessor<dim>
  {
  public:
    Postprocessor (const unsigned int partition,
                   const double       minimal_pressure);

    virtual
    void
    compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                       const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                       const std::vector<std::vector<Tensor<2,dim> > > &dduh,
                                       const std::vector<Point<dim> >                  &normals,
                                       const std::vector<Point<dim> >                  &evaluation_points,
                                       std::vector<Vector<double> >                    &computed_quantities) const;

    virtual std::vector<std::string> get_names () const;

    virtual
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    get_data_component_interpretation () const;

    virtual UpdateFlags get_needed_update_flags () const;

  private:
    const unsigned int partition;
    const double       minimal_pressure;
  };


  template <int dim>
  BoussinesqFlowProblem<dim>::Postprocessor::
  Postprocessor (const unsigned int partition,
                 const double       minimal_pressure)
    :
    partition (partition),
    minimal_pressure (minimal_pressure)
  {}


  // Here we define the names for the variables we want to output. These are
  // the actual solution values for velocity, pressure, and temperature, as
  // well as the friction heating and to each cell the number of the processor
  // that owns it. This allows us to visualize the partitioning of the domain
  // among the processors. Except for the velocity, which is vector-valued,
  // all other quantities are scalar.
  template <int dim>
  std::vector<std::string>
  BoussinesqFlowProblem<dim>::Postprocessor::get_names() const
  {
    std::vector<std::string> solution_names (dim, "velocity");
    solution_names.push_back ("p");
    solution_names.push_back ("T");
    solution_names.push_back ("friction_heating");
    solution_names.push_back ("partition");

    return solution_names;
  }


  template <int dim>
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
  BoussinesqFlowProblem<dim>::Postprocessor::
  get_data_component_interpretation () const
  {
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    interpretation (dim,
                    DataComponentInterpretation::component_is_part_of_vector);

    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    interpretation.push_back (DataComponentInterpretation::component_is_scalar);

    return interpretation;
  }


  template <int dim>
  UpdateFlags
  BoussinesqFlowProblem<dim>::Postprocessor::get_needed_update_flags() const
  {
    return update_values | update_gradients | update_q_points;
  }


  // Now we implement the function that computes the derived quantities. As we
  // also did for the output, we rescale the velocity from its SI units to
  // something more readable, namely cm/year. Next, the pressure is scaled to
  // be between 0 and the maximum pressure. This makes it more easily
  // comparable -- in essence making all pressure variables positive or
  // zero. Temperature is taken as is, and the friction heating is computed as
  // $2 \eta \varepsilon(\mathbf{u}) \cdot \varepsilon(\mathbf{u})$.
  //
  // The quantities we output here are more for illustration, rather than for
  // actual scientific value. We come back to this briefly in the results
  // section of this program and explain what one may in fact be interested
  // in.
  template <int dim>
  void
  BoussinesqFlowProblem<dim>::Postprocessor::
  compute_derived_quantities_vector (const std::vector<Vector<double> >              &uh,
                                     const std::vector<std::vector<Tensor<1,dim> > > &duh,
                                     const std::vector<std::vector<Tensor<2,dim> > > &/*dduh*/,
                                     const std::vector<Point<dim> >                  &/*normals*/,
                                     const std::vector<Point<dim> >                  &/*evaluation_points*/,
                                     std::vector<Vector<double> >                    &computed_quantities) const
  {
    const unsigned int n_quadrature_points = uh.size();
    Assert (duh.size() == n_quadrature_points,                  ExcInternalError());
    Assert (computed_quantities.size() == n_quadrature_points,  ExcInternalError());
    Assert (uh[0].size() == dim+2,                              ExcInternalError());

    for (unsigned int q=0; q<n_quadrature_points; ++q)
      {
        for (unsigned int d=0; d<dim; ++d)
          computed_quantities[q](d)
            = (uh[q](d) *  EquationData::year_in_seconds * 100);

        const double pressure = (uh[q](dim)-minimal_pressure);
        computed_quantities[q](dim) = pressure;

        const double temperature = uh[q](dim+1);
        computed_quantities[q](dim+1) = temperature;

        Tensor<2,dim> grad_u;
        for (unsigned int d=0; d<dim; ++d)
          grad_u[d] = duh[q][d];
        const SymmetricTensor<2,dim> strain_rate = symmetrize (grad_u);
        computed_quantities[q](dim+2) = 2 * EquationData::eta *
                                        strain_rate * strain_rate;

        computed_quantities[q](dim+3) = partition;
      }
  }


  // The <code>output_results()</code> function has a similar task to the one
  // in step-31. However, here we are going to demonstrate a different
  // technique on how to merge output from different DoFHandler objects. The
  // way we're going to achieve this recombination is to create a joint
  // DoFHandler that collects both components, the Stokes solution and the
  // temperature solution. This can be nicely done by combining the finite
  // elements from the two systems to form one FESystem, and let this
  // collective system define a new DoFHandler object. To be sure that
  // everything was done correctly, we perform a sanity check that ensures
  // that we got all the dofs from both Stokes and temperature even in the
  // combined system. We then combine the data vectors. Unfortunately, there
  // is no straight-forward relation that tells us how to sort Stokes and
  // temperature vector into the joint vector. The way we can get around this
  // trouble is to rely on the information collected in the FESystem. For each
  // dof on a cell, the joint finite element knows to which equation component
  // (velocity component, pressure, or temperature) it belongs – that's the
  // information we need! So we step through all cells (with iterators into
  // all three DoFHandlers moving in sync), and for each joint cell dof, we
  // read out that component using the FiniteElement::system_to_base_index
  // function (see there for a description of what the various parts of its
  // return value contain). We also need to keep track whether we're on a
  // Stokes dof or a temperature dof, which is contained in
  // joint_fe.system_to_base_index(i).first.first. Eventually, the dof_indices
  // data structures on either of the three systems tell us how the relation
  // between global vector and local dofs looks like on the present cell,
  // which concludes this tedious work. We make sure that each processor only
  // works on the subdomain it owns locally (and not on ghost or artificial
  // cells) when building the joint solution vector. The same will then have
  // to be done in DataOut::build_patches(), but that function does so
  // automatically.
  //
  // What we end up with is a set of patches that we can write using the
  // functions in DataOutBase in a variety of output formats. Here, we then
  // have to pay attention that what each processor writes is really only its
  // own part of the domain, i.e. we will want to write each processor's
  // contribution into a separate file. This we do by adding an additional
  // number to the filename when we write the solution. This is not really
  // new, we did it similarly in step-40. Note that we write in the compressed
  // format @p .vtu instead of plain vtk files, which saves quite some
  // storage.
  //
  // All the rest of the work is done in the PostProcessor class.
  template <int dim>
  void BoussinesqFlowProblem<dim>::output_results ()
  {
    computing_timer.enter_section ("Postprocessing");

    const FESystem<dim> joint_fe (stokes_fe, 1,
                                  temperature_fe, 1);

    DoFHandler<dim> joint_dof_handler (triangulation);
    joint_dof_handler.distribute_dofs (joint_fe);
    Assert (joint_dof_handler.n_dofs() ==
            stokes_dof_handler.n_dofs() + temperature_dof_handler.n_dofs(),
            ExcInternalError());

    TrilinosWrappers::MPI::Vector joint_solution;
    joint_solution.reinit (joint_dof_handler.locally_owned_dofs(), MPI_COMM_WORLD);

    {
      std::vector<types::global_dof_index> local_joint_dof_indices (joint_fe.dofs_per_cell);
      std::vector<types::global_dof_index> local_stokes_dof_indices (stokes_fe.dofs_per_cell);
      std::vector<types::global_dof_index> local_temperature_dof_indices (temperature_fe.dofs_per_cell);

      typename DoFHandler<dim>::active_cell_iterator
      joint_cell       = joint_dof_handler.begin_active(),
      joint_endc       = joint_dof_handler.end(),
      stokes_cell      = stokes_dof_handler.begin_active(),
      temperature_cell = temperature_dof_handler.begin_active();
      for (; joint_cell!=joint_endc;
           ++joint_cell, ++stokes_cell, ++temperature_cell)
        if (joint_cell->is_locally_owned())
          {
            joint_cell->get_dof_indices (local_joint_dof_indices);
            stokes_cell->get_dof_indices (local_stokes_dof_indices);
            temperature_cell->get_dof_indices (local_temperature_dof_indices);

            for (unsigned int i=0; i<joint_fe.dofs_per_cell; ++i)
              if (joint_fe.system_to_base_index(i).first.first == 0)
                {
                  Assert (joint_fe.system_to_base_index(i).second
                          <
                          local_stokes_dof_indices.size(),
                          ExcInternalError());

                  joint_solution(local_joint_dof_indices[i])
                    = stokes_solution(local_stokes_dof_indices
                                      [joint_fe.system_to_base_index(i).second]);
                }
              else
                {
                  Assert (joint_fe.system_to_base_index(i).first.first == 1,
                          ExcInternalError());
                  Assert (joint_fe.system_to_base_index(i).second
                          <
                          local_temperature_dof_indices.size(),
                          ExcInternalError());
                  joint_solution(local_joint_dof_indices[i])
                    = temperature_solution(local_temperature_dof_indices
                                           [joint_fe.system_to_base_index(i).second]);
                }
          }
    }

    joint_solution.compress(VectorOperation::insert);

    IndexSet locally_relevant_joint_dofs(joint_dof_handler.n_dofs());
    DoFTools::extract_locally_relevant_dofs (joint_dof_handler, locally_relevant_joint_dofs);
    TrilinosWrappers::MPI::Vector locally_relevant_joint_solution;
    locally_relevant_joint_solution.reinit (locally_relevant_joint_dofs, MPI_COMM_WORLD);
    locally_relevant_joint_solution = joint_solution;

    Postprocessor postprocessor (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD),
                                 stokes_solution.block(1).minimal_value());

    DataOut<dim> data_out;
    data_out.attach_dof_handler (joint_dof_handler);
    data_out.add_data_vector (locally_relevant_joint_solution, postprocessor);
    data_out.build_patches ();

    static int out_index=0;
    const std::string filename = ("solution-" +
                                  Utilities::int_to_string (out_index, 5) +
                                  "." +
                                  Utilities::int_to_string
                                  (triangulation.locally_owned_subdomain(), 4) +
                                  ".vtu");
    std::ofstream output (filename.c_str());
    data_out.write_vtu (output);


    // At this point, all processors have written their own files to disk. We
    // could visualize them individually in Visit or Paraview, but in reality
    // we of course want to visualize the whole set of files at once. To this
    // end, we create a master file in each of the formats understood by Visit
    // (<code>.visit</code>) and Paraview (<code>.pvtu</code>) on the zeroth
    // processor that describes how the individual files are defining the
    // global data set.
    if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
      {
        std::vector<std::string> filenames;
        for (unsigned int i=0; i<Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); ++i)
          filenames.push_back (std::string("solution-") +
                               Utilities::int_to_string (out_index, 5) +
                               "." +
                               Utilities::int_to_string(i, 4) +
                               ".vtu");
        const std::string
        pvtu_master_filename = ("solution-" +
                                Utilities::int_to_string (out_index, 5) +
                                ".pvtu");
        std::ofstream pvtu_master (pvtu_master_filename.c_str());
        data_out.write_pvtu_record (pvtu_master, filenames);

        const std::string
        visit_master_filename = ("solution-" +
                                 Utilities::int_to_string (out_index, 5) +
                                 ".visit");
        std::ofstream visit_master (visit_master_filename.c_str());
        data_out.write_visit_record (visit_master, filenames);
      }

    computing_timer.exit_section ();
    out_index++;
  }



  // @sect4{BoussinesqFlowProblem::refine_mesh}

  // This function isn't really new either. Since the <code>setup_dofs</code>
  // function that we call in the middle has its own timer section, we split
  // timing this function into two sections. It will also allow us to easily
  // identify which of the two is more expensive.
  //
  // One thing of note, however, is that we only want to compute error
  // indicators on the locally owned subdomain. In order to achieve this, we
  // pass one additional argument to the KellyErrorEstimator::estimate
  // function. Note that the vector for error estimates is resized to the
  // number of active cells present on the current process, which is less than
  // the total number of active cells on all processors (but more than the
  // number of locally owned active cells); each processor only has a few
  // coarse cells around the locally owned ones, as also explained in step-40.
  //
  // The local error estimates are then handed to a %parallel version of
  // GridRefinement (in namespace parallel::distributed::GridRefinement, see
  // also step-40) which looks at the errors and finds the cells that need
  // refinement by comparing the error values across processors. As in
  // step-31, we want to limit the maximum grid level. So in case some cells
  // have been marked that are already at the finest level, we simply clear
  // the refine flags.
  template <int dim>
  void BoussinesqFlowProblem<dim>::refine_mesh (const unsigned int max_grid_level)
  {
    computing_timer.enter_section ("Refine mesh structure, part 1");
    Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

    KellyErrorEstimator<dim>::estimate (temperature_dof_handler,
                                        QGauss<dim-1>(parameters.temperature_degree+1),
                                        typename FunctionMap<dim>::type(),
                                        temperature_solution,
                                        estimated_error_per_cell,
                                        ComponentMask(),
                                        0,
                                        0,
                                        triangulation.locally_owned_subdomain());

    parallel::distributed::GridRefinement::
    refine_and_coarsen_fixed_fraction (triangulation,
                                       estimated_error_per_cell,
                                       0.3, 0.1);

    if (triangulation.n_levels() > max_grid_level)
      for (typename Triangulation<dim>::active_cell_iterator
           cell = triangulation.begin_active(max_grid_level);
           cell != triangulation.end(); ++cell)
        cell->clear_refine_flag ();

    // With all flags marked as necessary, we set up the
    // parallel::distributed::SolutionTransfer object to transfer the
    // solutions for the current time level and the next older one. The syntax
    // is similar to the non-%parallel solution transfer (with the exception
    // that here a pointer to the vector entries is enough). The remainder of
    // the function is concerned with setting up the data structures again
    // after mesh refinement and restoring the solution vectors on the new
    // mesh.
    std::vector<const TrilinosWrappers::MPI::Vector *> x_temperature (2);
    x_temperature[0] = &temperature_solution;
    x_temperature[1] = &old_temperature_solution;
    std::vector<const TrilinosWrappers::MPI::BlockVector *> x_stokes (2);
    x_stokes[0] = &stokes_solution;
    x_stokes[1] = &old_stokes_solution;

    parallel::distributed::SolutionTransfer<dim,TrilinosWrappers::MPI::Vector>
    temperature_trans(temperature_dof_handler);
    parallel::distributed::SolutionTransfer<dim,TrilinosWrappers::MPI::BlockVector>
    stokes_trans(stokes_dof_handler);

    triangulation.prepare_coarsening_and_refinement();
    temperature_trans.prepare_for_coarsening_and_refinement(x_temperature);
    stokes_trans.prepare_for_coarsening_and_refinement(x_stokes);

    triangulation.execute_coarsening_and_refinement ();
    computing_timer.exit_section();

    setup_dofs ();

    computing_timer.enter_section ("Refine mesh structure, part 2");

    {
      TrilinosWrappers::MPI::Vector distributed_temp1 (temperature_rhs);
      TrilinosWrappers::MPI::Vector distributed_temp2 (temperature_rhs);

      std::vector<TrilinosWrappers::MPI::Vector *> tmp (2);
      tmp[0] = &(distributed_temp1);
      tmp[1] = &(distributed_temp2);
      temperature_trans.interpolate(tmp);

      temperature_solution     = distributed_temp1;
      old_temperature_solution = distributed_temp2;
    }

    {
      TrilinosWrappers::MPI::BlockVector distributed_stokes (stokes_rhs);
      TrilinosWrappers::MPI::BlockVector old_distributed_stokes (stokes_rhs);

      std::vector<TrilinosWrappers::MPI::BlockVector *> stokes_tmp (2);
      stokes_tmp[0] = &(distributed_stokes);
      stokes_tmp[1] = &(old_distributed_stokes);

      stokes_trans.interpolate (stokes_tmp);
      stokes_solution     = distributed_stokes;
      old_stokes_solution = old_distributed_stokes;
    }

    computing_timer.exit_section();
  }



  // @sect4{BoussinesqFlowProblem::run}

  // This is the final and controlling function in this class. It, in fact,
  // runs the entire rest of the program and is, once more, very similar to
  // step-31. We use a different mesh now (a GridGenerator::hyper_shell
  // instead of a simple cube geometry), and use the
  // <code>project_temperature_field()</code> function instead of the library
  // function <code>VectorTools::project</code>, the rest is as before.
  template <int dim>
  void BoussinesqFlowProblem<dim>::run ()
  {
    GridGenerator::hyper_shell (triangulation,
                                Point<dim>(),
                                EquationData::R0,
                                EquationData::R1,
                                (dim==3) ? 96 : 12,
                                true);
    static HyperShellBoundary<dim> boundary;
    triangulation.set_boundary (0, boundary);
    triangulation.set_boundary (1, boundary);

    global_Omega_diameter = GridTools::diameter (triangulation);

    triangulation.refine_global (parameters.initial_global_refinement);

    setup_dofs();

    unsigned int pre_refinement_step = 0;

start_time_iteration:

    project_temperature_field ();

    timestep_number           = 0;
    time_step = old_time_step = 0;

    double time = 0;

    do
      {
        pcout << "Timestep " << timestep_number
              << ":  t=" << time/EquationData::year_in_seconds
              << " years"
              << std::endl;

        assemble_stokes_system ();
        build_stokes_preconditioner ();
        assemble_temperature_matrix ();

        solve ();

        pcout << std::endl;

        if ((timestep_number == 0) &&
            (pre_refinement_step < parameters.initial_adaptive_refinement))
          {
            refine_mesh (parameters.initial_global_refinement +
                         parameters.initial_adaptive_refinement);
            ++pre_refinement_step;
            goto start_time_iteration;
          }
        else if ((timestep_number > 0)
                 &&
                 (timestep_number % parameters.adaptive_refinement_interval == 0))
          refine_mesh (parameters.initial_global_refinement +
                       parameters.initial_adaptive_refinement);

        if ((parameters.generate_graphical_output == true)
            &&
            (timestep_number % parameters.graphical_output_interval == 0))
          output_results ();

        // In order to speed up linear solvers, we extrapolate the solutions
        // from the old time levels to the new one. This gives a very good
        // initial guess, cutting the number of iterations needed in solvers
        // by more than one half. We do not need to extrapolate in the last
        // iteration, so if we reached the final time, we stop here.
        //
        // As the last thing during a time step (before actually bumping up
        // the number of the time step), we check whether the current time
        // step number is divisible by 100, and if so we let the computing
        // timer print a summary of CPU times spent so far.
        if (time > parameters.end_time * EquationData::year_in_seconds)
          break;

        TrilinosWrappers::MPI::BlockVector old_old_stokes_solution;
        old_old_stokes_solution      = old_stokes_solution;
        old_stokes_solution          = stokes_solution;
        old_old_temperature_solution = old_temperature_solution;
        old_temperature_solution     = temperature_solution;
        if (old_time_step > 0)
          {
            //Trilinos sadd does not like ghost vectors even as input. Copy
            //into distributed vectors for now:
            {
              TrilinosWrappers::MPI::BlockVector distr_solution (stokes_rhs);
              distr_solution = stokes_solution;
              TrilinosWrappers::MPI::BlockVector distr_old_solution (stokes_rhs);
              distr_old_solution = old_old_stokes_solution;
              distr_solution .sadd (1.+time_step/old_time_step, -time_step/old_time_step,
                                    distr_old_solution);
              stokes_solution = distr_solution;
            }
            {
              TrilinosWrappers::MPI::Vector distr_solution (temperature_rhs);
              distr_solution = temperature_solution;
              TrilinosWrappers::MPI::Vector distr_old_solution (temperature_rhs);
              distr_old_solution = old_old_temperature_solution;
              distr_solution .sadd (1.+time_step/old_time_step, -time_step/old_time_step,
                                    distr_old_solution);
              temperature_solution = distr_solution;
            }
          }

        if ((timestep_number > 0) && (timestep_number % 100 == 0))
          computing_timer.print_summary ();

        time += time_step;
        ++timestep_number;
      }
    while (true);

    // If we are generating graphical output, do so also for the last time
    // step unless we had just done so before we left the do-while loop
    if ((parameters.generate_graphical_output == true)
        &&
        !((timestep_number-1) % parameters.graphical_output_interval == 0))
      output_results ();
  }
}



// @sect3{The <code>main</code> function}

// The main function is short as usual and very similar to the one in
// step-31. Since we use a parameter file which is specified as an argument in
// the command line, we have to read it in here and pass it on to the
// Parameters class for parsing. If no filename is given in the command line,
// we simply use the <code>\step-32.prm</code> file which is distributed
// together with the program.
//
// Because 3d computations are simply very slow unless you throw a lot of
// processors at them, the program defaults to 2d. You can get the 3d version
// by changing the constant dimension below to 3.
int main (int argc, char *argv[])
{
  using namespace Step32;
  using namespace dealii;

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv,
                                                      numbers::invalid_unsigned_int);

  try
    {
      deallog.depth_console (0);

      std::string parameter_filename;
      if (argc>=2)
        parameter_filename = argv[1];
      else
        parameter_filename = "step-32.prm";

      const int dim = 2;
      BoussinesqFlowProblem<dim>::Parameters  parameters(parameter_filename);
      BoussinesqFlowProblem<dim> flow_problem (parameters);
      flow_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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
