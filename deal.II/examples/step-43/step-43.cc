/* Author: Chih-Che Chueh, University of Victoria, 2010 */
/*         Wolfgang Bangerth, Texas A&M University, 2010 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2010, 2011, 2012 by Chih-Che Chueh and the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


				 // @sect3{Include files}

				 // The first step, as always, is to
				 // include the functionality of a
				 // number of deal.II and C++ header
				 // files.
				 //
				 // The list includes some header
				 // files that provide vector, matrix,
				 // and preconditioner classes that
				 // implement interfaces to the
				 // respective Trilinos classes; some
				 // more information on these may be
				 // found in step-31.
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/constraint_matrix.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vectors.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>

#include <fstream>
#include <sstream>


				 // At the end of this top-matter, we
				 // open a namespace for the current
				 // project into which all the
				 // following material will go, and
				 // then import all deal.II names into
				 // this namespace:
namespace Step43
{
  using namespace dealii;


				   // @sect3{Pressure right hand side, pressure boundary values and saturation initial value classes}

				   // The following part is taken
				   // directly from step-21 so there is
				   // no need to repeat the
				   // descriptions found there.
  template <int dim>
  class PressureRightHandSide : public Function<dim>
  {
    public:
      PressureRightHandSide () : Function<dim>(1) {}

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;
  };



  template <int dim>
  double
  PressureRightHandSide<dim>::value (const Point<dim>  &/*p*/,
				     const unsigned int /*component*/) const
  {
    return 0;
  }


  template <int dim>
  class PressureBoundaryValues : public Function<dim>
  {
    public:
      PressureBoundaryValues () : Function<dim>(1) {}

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;
  };


  template <int dim>
  double
  PressureBoundaryValues<dim>::value (const Point<dim>  &p,
				      const unsigned int /*component*/) const
  {
    return 1-p[0];
  }


  template <int dim>
  class SaturationBoundaryValues : public Function<dim>
  {
    public:
      SaturationBoundaryValues () : Function<dim>(1) {}

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;
  };



  template <int dim>
  double
  SaturationBoundaryValues<dim>::value (const Point<dim> &p,
					const unsigned int /*component*/) const
  {
    if (p[0] == 0)
      return 1;
    else
      return 0;
  }


  template <int dim>
  class SaturationInitialValues : public Function<dim>
  {
    public:
      SaturationInitialValues () : Function<dim>(1) {}

      virtual double value (const Point<dim>   &p,
			    const unsigned int  component = 0) const;

      virtual void vector_value (const Point<dim> &p,
				 Vector<double>   &value) const;
  };


  template <int dim>
  double
  SaturationInitialValues<dim>::value (const Point<dim>  &/*p*/,
				       const unsigned int /*component*/) const
  {
    return 0;
  }


  template <int dim>
  void
  SaturationInitialValues<dim>::vector_value (const Point<dim> &p,
					      Vector<double>   &values) const
  {
    for (unsigned int c=0; c<this->n_components; ++c)
      values(c) = SaturationInitialValues<dim>::value (p,c);
  }


				   // @sect3{Permeability models}

				   // In this tutorial, we still use
				   // the two permeability models
				   // previously used in step-21 so we
				   // again refrain from commenting in
				   // detail about them.
  namespace SingleCurvingCrack
  {
    template <int dim>
    class KInverse : public TensorFunction<2,dim>
    {
      public:
	KInverse ()
			:
			TensorFunction<2,dim> ()
	  {}

	virtual void value_list (const std::vector<Point<dim> > &points,
				 std::vector<Tensor<2,dim> >    &values) const;
    };


    template <int dim>
    void
    KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
			       std::vector<Tensor<2,dim> >    &values) const
    {
      Assert (points.size() == values.size(),
	      ExcDimensionMismatch (points.size(), values.size()));

      for (unsigned int p=0; p<points.size(); ++p)
	{
	  values[p].clear ();

	  const double distance_to_flowline
	    = std::fabs(points[p][1]-0.5-0.1*std::sin(10*points[p][0]));

	  const double permeability = std::max(std::exp(-(distance_to_flowline*
							  distance_to_flowline)
							/ (0.1 * 0.1)),
					       0.01);

	  for (unsigned int d=0; d<dim; ++d)
	    values[p][d][d] = 1./permeability;
	}
    }
  }


  namespace RandomMedium
  {
    template <int dim>
    class KInverse : public TensorFunction<2,dim>
    {
      public:
	KInverse ()
			:
			TensorFunction<2,dim> ()
	  {}

	virtual void value_list (const std::vector<Point<dim> > &points,
				 std::vector<Tensor<2,dim> >    &values) const;

      private:
	static std::vector<Point<dim> > centers;

	static std::vector<Point<dim> > get_centers ();
    };



    template <int dim>
    std::vector<Point<dim> >
    KInverse<dim>::centers = KInverse<dim>::get_centers();


    template <int dim>
    std::vector<Point<dim> >
    KInverse<dim>::get_centers ()
    {
      const unsigned int N = (dim == 2 ?
			      40 :
			      (dim == 3 ?
			       100 :
			       throw ExcNotImplemented()));

      std::vector<Point<dim> > centers_list (N);
      for (unsigned int i=0; i<N; ++i)
	for (unsigned int d=0; d<dim; ++d)
	  centers_list[i][d] = static_cast<double>(rand())/RAND_MAX;

      return centers_list;
    }



    template <int dim>
    void
    KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
			       std::vector<Tensor<2,dim> >    &values) const
    {
      Assert (points.size() == values.size(),
	      ExcDimensionMismatch (points.size(), values.size()));

      for (unsigned int p=0; p<points.size(); ++p)
	{
	  values[p].clear ();

	  double permeability = 0;
	  for (unsigned int i=0; i<centers.size(); ++i)
	    permeability += std::exp(-(points[p]-centers[i]).square()
				     / (0.05 * 0.05));

	  const double normalized_permeability
	    = std::min (std::max(permeability, 0.01), 4.);

	  for (unsigned int d=0; d<dim; ++d)
	    values[p][d][d] = 1./normalized_permeability;
	}
    }
  }


				   // @sect3{Physical quantities}

				   // The implementations of all the
				   // physical quantities such as
				   // total mobility $\lambda_t$ and
				   // fractional flow of water $F$ are
				   // taken from step-21 so again we
				   // don't have do any comment about
				   // them. Compared to step-21 we
				   // have added checks that the
				   // saturation passed to these
				   // functions is in fact within the
				   // physically valid
				   // range. Furthermore, given that
				   // the wetting phase moves at speed
				   // $\mathbf u F'(S)$ it is clear
				   // that $F'(S)$ must be greater or
				   // equal to zero, so we assert that
				   // as well to make sure that our
				   // calculations to get at the
				   // formula for the derivative made
				   // sense.
  double mobility_inverse (const double S,
			   const double viscosity)
  {
    return 1.0 / (1.0/viscosity * S * S + (1-S) * (1-S));
  }


  double fractional_flow (const double S,
			  const double viscosity)
  {
    Assert ((S >= 0) && (S<=1),
	    ExcMessage ("Saturation is outside its physically valid range."));

    return S*S / ( S * S + viscosity * (1-S) * (1-S));
  }


  double fractional_flow_derivative (const double S,
				     const double viscosity)
  {
    Assert ((S >= 0) && (S<=1),
	    ExcMessage ("Saturation is outside its physically valid range."));

    const double temp = ( S * S + viscosity * (1-S) * (1-S) );

    const double numerator   =  2.0 * S * temp
				-
				S * S *
				( 2.0 * S - 2.0 * viscosity * (1-S) );
    const double denominator =  std::pow(temp, 2.0);

    const double F_prime = numerator / denominator;

    Assert (F_prime >= 0, ExcInternalError());

    return F_prime;
  }


				   // @sect3{Helper classes for solvers and preconditioners}

				   // In this first part we define a
				   // number of classes that we need
				   // in the construction of linear
				   // solvers and
				   // preconditioners. This part is
				   // essentially the same as that
				   // used in step-31. The only
				   // difference is that the original
				   // variable name stokes_matrix is
				   // replaced by another name
				   // darcy_matrix to match our
				   // problem.
  namespace LinearSolvers
  {
    template <class Matrix, class Preconditioner>
    class InverseMatrix : public Subscriptor
    {
      public:
	InverseMatrix (const Matrix         &m,
		       const Preconditioner &preconditioner);


	template <typename VectorType>
	void vmult (VectorType       &dst,
		    const VectorType &src) const;

      private:
	const SmartPointer<const Matrix> matrix;
	const Preconditioner &preconditioner;
    };


    template <class Matrix, class Preconditioner>
    InverseMatrix<Matrix,Preconditioner>::
    InverseMatrix (const Matrix &m,
		   const Preconditioner &preconditioner)
		    :
		    matrix (&m),
		    preconditioner (preconditioner)
    {}



    template <class Matrix, class Preconditioner>
    template <typename VectorType>
    void
    InverseMatrix<Matrix,Preconditioner>::
    vmult (VectorType       &dst,
	   const VectorType &src) const
    {
      SolverControl solver_control (src.size(), 1e-7*src.l2_norm());
      SolverCG<VectorType> cg (solver_control);

      dst = 0;

      try
	{
	  cg.solve (*matrix, dst, src, preconditioner);
	}
      catch (std::exception &e)
	{
	  Assert (false, ExcMessage(e.what()));
	}
    }

    template <class PreconditionerA, class PreconditionerMp>
    class BlockSchurPreconditioner : public Subscriptor
    {
      public:
	BlockSchurPreconditioner (
	  const TrilinosWrappers::BlockSparseMatrix     &S,
	  const InverseMatrix<TrilinosWrappers::SparseMatrix,
			      PreconditionerMp>         &Mpinv,
	  const PreconditionerA                         &Apreconditioner);

	void vmult (TrilinosWrappers::BlockVector       &dst,
		    const TrilinosWrappers::BlockVector &src) const;

      private:
	const SmartPointer<const TrilinosWrappers::BlockSparseMatrix> darcy_matrix;
	const SmartPointer<const InverseMatrix<TrilinosWrappers::SparseMatrix,
					       PreconditionerMp > > m_inverse;
	const PreconditionerA &a_preconditioner;

	mutable TrilinosWrappers::Vector tmp;
    };



    template <class PreconditionerA, class PreconditionerMp>
    BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::
    BlockSchurPreconditioner(const TrilinosWrappers::BlockSparseMatrix  &S,
			     const InverseMatrix<TrilinosWrappers::SparseMatrix,
						 PreconditionerMp>      &Mpinv,
			     const PreconditionerA                      &Apreconditioner)
		    :
		    darcy_matrix            (&S),
		    m_inverse               (&Mpinv),
		    a_preconditioner        (Apreconditioner),
		    tmp                     (darcy_matrix->block(1,1).m())
    {}


    template <class PreconditionerA, class PreconditionerMp>
    void BlockSchurPreconditioner<PreconditionerA, PreconditionerMp>::vmult (
      TrilinosWrappers::BlockVector       &dst,
      const TrilinosWrappers::BlockVector &src) const
    {
      a_preconditioner.vmult (dst.block(0), src.block(0));
      darcy_matrix->block(1,0).residual(tmp, dst.block(0), src.block(1));
      tmp *= -1;
      m_inverse->vmult (dst.block(1), tmp);
    }
  }


				   // @sect3{The TwoPhaseFlowProblem class}

				   // The definition of the class that
				   // defines the top-level logic of
				   // solving the time-dependent
				   // advection-dominated two-phase
				   // flow problem (or
				   // Buckley-Leverett problem
				   // [Buckley 1942]) is mainly based
				   // on tutorial programs step-21 and
				   // step-33, and in particular on
				   // step-31 where we have used
				   // basically the same general
				   // structure as done here. As in
				   // step-31, the key routines to
				   // look for in the implementation
				   // below are the <code>run()</code>
				   // and <code>solve()</code>
				   // functions.
				   //
				   // The main difference to step-31
				   // is that, since adaptive operator
				   // splitting is considered, we need
				   // a couple more member variables
				   // to hold the last two computed
				   // Darcy (velocity/pressure)
				   // solutions in addition to the
				   // current one (which is either
				   // computed directly, or
				   // extrapolated from the previous
				   // two), and we need to remember
				   // the last two times we computed
				   // the Darcy solution. We also need
				   // a helper function that figures
				   // out whether we do indeed need to
				   // recompute the Darcy solution.
				   //
				   // Unlike step-31, this step uses
				   // one more ConstraintMatrix object
				   // called
				   // darcy_preconditioner_constraints. This
				   // constraint object is used only
				   // for assembling the matrix for
				   // the Darcy preconditioner and
				   // includes hanging node constrants
				   // as well as Dirichlet boundary
				   // value constraints for the
				   // pressure variable. We need this
				   // because we are building a
				   // Laplace matrix for the pressure
				   // as an approximation of the Schur
				   // complement) which is only
				   // positive definite if boundary
				   // conditions are applied.
				   //
				   // The collection of member
				   // functions and variables thus
				   // declared in this class is then
				   // rather similar to those in
				   // step-31:
  template <int dim>
  class TwoPhaseFlowProblem
  {
    public:
      TwoPhaseFlowProblem (const unsigned int degree);
      void run ();

    private:
      void setup_dofs ();
      void assemble_darcy_preconditioner ();
      void build_darcy_preconditioner ();
      void assemble_darcy_system ();
      void assemble_saturation_system ();
      void assemble_saturation_matrix ();
      void assemble_saturation_rhs ();
      void assemble_saturation_rhs_cell_term (const FEValues<dim>             &saturation_fe_values,
					      const FEValues<dim>             &darcy_fe_values,
					      const double                     global_max_u_F_prime,
					      const double                     global_S_variation,
					      const std::vector<unsigned int> &local_dof_indices);
      void assemble_saturation_rhs_boundary_term (const FEFaceValues<dim>             &saturation_fe_face_values,
						  const FEFaceValues<dim>             &darcy_fe_face_values,
						  const std::vector<unsigned int>     &local_dof_indices);
      void solve ();
      void compute_refinement_indicators (const TrilinosWrappers::Vector &predicted_saturation_solution,
					  Vector<double> &refinement_indicators) const;
      void refine_mesh (const unsigned int    min_grid_level,
			const unsigned int    max_grid_level,
			const Vector<double> &indicator);
      void output_results () const;

				       // We follow with a number of
				       // helper functions that are
				       // used in a variety of places
				       // throughout the program:
      double                   get_max_u_F_prime () const;
      std::pair<double,double> get_extrapolated_saturation_range () const;
      bool                     determine_whether_to_solve_for_pressure_and_velocity () const;
      void                     project_back_saturation ();
      double                   compute_viscosity (const std::vector<double>          &old_saturation,
						  const std::vector<double>          &old_old_saturation,
						  const std::vector<Tensor<1,dim> >  &old_saturation_grads,
						  const std::vector<Tensor<1,dim> >  &old_old_saturation_grads,
						  const std::vector<Vector<double> > &present_darcy_values,
						  const double                        global_max_u_F_prime,
						  const double                        global_S_variation,
						  const double                        cell_diameter,
						  const double                        old_time_step,
						  const double                        viscosity,
						  const double                        porosity) const;


				       // This all is followed by the
				       // member variables, most of
				       // which are similar to the
				       // ones in step-31, with the
				       // exception of the ones that
				       // pertain to the macro time
				       // stepping for the
				       // velocity/pressure system:
      const unsigned int degree;

      Triangulation<dim>                   triangulation;
      double                               global_Omega_diameter;

      const unsigned int                   darcy_degree;
      FESystem<dim>                        darcy_fe;
      DoFHandler<dim>                      darcy_dof_handler;
      ConstraintMatrix                     darcy_constraints;

      ConstraintMatrix                     darcy_preconditioner_constraints;

      TrilinosWrappers::BlockSparseMatrix  darcy_matrix;
      TrilinosWrappers::BlockSparseMatrix  darcy_preconditioner_matrix;

      TrilinosWrappers::BlockVector        darcy_solution;
      TrilinosWrappers::BlockVector        darcy_rhs;

      TrilinosWrappers::BlockVector        last_computed_darcy_solution;
      TrilinosWrappers::BlockVector        second_last_computed_darcy_solution;


      const unsigned int                   saturation_degree;
      FE_Q<dim>                            saturation_fe;
      DoFHandler<dim>                      saturation_dof_handler;
      ConstraintMatrix                     saturation_constraints;

      TrilinosWrappers::SparseMatrix       saturation_matrix;


      TrilinosWrappers::Vector             saturation_solution;
      TrilinosWrappers::Vector             old_saturation_solution;
      TrilinosWrappers::Vector             old_old_saturation_solution;
      TrilinosWrappers::Vector             saturation_rhs;

      TrilinosWrappers::Vector             saturation_matching_last_computed_darcy_solution;

      const double                         saturation_refinement_threshold;

      double                               current_macro_time_step;
      double                               old_macro_time_step;

      double                               time_step;
      double                               old_time_step;
      unsigned int                         timestep_number;
      double                               viscosity;
      double                               porosity;
      double                               AOS_threshold;

      std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionIC> Amg_preconditioner;
      std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionIC> Mp_preconditioner;

      bool                                rebuild_saturation_matrix;

				       // At the very end we declare a
				       // variable that denotes the
				       // material model. Compared to
				       // step-21, we do this here as
				       // a member variable since we
				       // will want to use it in a
				       // variety of places and so
				       // having a central place where
				       // such a variable is declared
				       // will make it simpler to
				       // replace one class by another
				       // (e.g. replace
				       // RandomMedium::KInverse by
				       // SingleCurvingCrack::KInverse).
      const RandomMedium::KInverse<dim>   k_inverse;
  };


				   // @sect3{TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem}

				   // The constructor of this class is an
				   // extension of the constructors in step-21
				   // and step-31. We need to add the various
				   // variables that concern the saturation. As
				   // discussed in the introduction, we are
				   // going to use $Q_2 \times Q_1$
				   // (Taylor-Hood) elements again for the Darcy
				   // system, an element combination that fulfills
				   // the Ladyzhenskaya-Babuska-Brezzi (LBB)
				   // conditions
				   // [Brezzi and Fortin 1991, Chen 2005], and $Q_1$
				   // elements for the saturation. However, by
				   // using variables that store the polynomial
				   // degree of the Darcy and temperature finite
				   // elements, it is easy to consistently
				   // modify the degree of the elements as well
				   // as all quadrature formulas used on them
				   // downstream. Moreover, we initialize the
				   // time stepping variables related to
				   // operator splitting as well as the option
				   // for matrix assembly and preconditioning:
  template <int dim>
  TwoPhaseFlowProblem<dim>::TwoPhaseFlowProblem (const unsigned int degree)
		  :
		  degree (degree),
		  darcy_degree (degree),
		  darcy_fe (FE_Q<dim>(darcy_degree+1), dim,
			    FE_Q<dim>(darcy_degree), 1),
		  darcy_dof_handler (triangulation),

		  saturation_degree (degree),
		  saturation_fe (saturation_degree),
		  saturation_dof_handler (triangulation),

		  saturation_refinement_threshold (0.5),

		  current_macro_time_step (0),
		  old_macro_time_step (0),

		  time_step (0),
		  old_time_step (0),
		  viscosity (0.2),
                  porosity (1.0),
                  AOS_threshold (3.0),

		  rebuild_saturation_matrix (true)
  {}


				   // @sect3{TwoPhaseFlowProblem<dim>::setup_dofs}

				   // This is the function that sets up the
				   // DoFHandler objects we have here (one for
				   // the Darcy part and one for the saturation
				   // part) as well as set to the right sizes
				   // the various objects required for the
				   // linear algebra in this program. Its basic
				   // operations are similar to what
				   // step-31 did.
				   //
				   // The body of the function first enumerates
				   // all degrees of freedom for the Darcy and
				   // saturation systems. For the Darcy part,
				   // degrees of freedom are then sorted to
				   // ensure that velocities precede pressure
				   // DoFs so that we can partition the Darcy
				   // matrix into a $2 \times 2$ matrix.
				   //
				   // Then, we need to incorporate
				   // hanging node constraints and
				   // Dirichlet boundary value
				   // constraints into
				   // darcy_preconditioner_constraints.
				   // The boundary condition
				   // constraints are only set on the
				   // pressure component since the
				   // Schur complement preconditioner
				   // that corresponds to the porous
				   // media flow operator in non-mixed
				   // form, $-\nabla \cdot [\mathbf K
				   // \lambda_t(S)]\nabla$, acts only
				   // on the pressure
				   // variable. Therefore, we use a
				   // component_mask that filters out
				   // the velocity component, so that
				   // the condensation is performed on
				   // pressure degrees of freedom
				   // only.
				   //
				   // After having done so, we count
				   // the number of degrees of freedom
				   // in the various blocks. This
				   // information is then used to
				   // create the sparsity pattern for
				   // the Darcy and saturation system
				   // matrices as well as the
				   // preconditioner matrix from which
				   // we build the Darcy
				   // preconditioner. As in step-31,
				   // we choose to create the pattern
				   // not as in the first few tutorial
				   // programs, but by using the
				   // blocked version of
				   // CompressedSimpleSparsityPattern. The
				   // reason for doing this is mainly
				   // memory, that is, the
				   // SparsityPattern class would
				   // consume too much memory when
				   // used in three spatial dimensions
				   // as we intend to do for this
				   // program. So, for this, we follow
				   // the same way as step-31 did and
				   // we don't have to repeat
				   // descriptions again for the rest
				   // of the member function.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::setup_dofs ()
  {
    std::vector<unsigned int> darcy_block_component (dim+1,0);
    darcy_block_component[dim] = 1;
    {
      darcy_dof_handler.distribute_dofs (darcy_fe);
      DoFRenumbering::Cuthill_McKee (darcy_dof_handler);
      DoFRenumbering::component_wise (darcy_dof_handler, darcy_block_component);

      darcy_constraints.clear ();
      DoFTools::make_hanging_node_constraints (darcy_dof_handler, darcy_constraints);
      darcy_constraints.close ();
    }
    {
      saturation_dof_handler.distribute_dofs (saturation_fe);

      saturation_constraints.clear ();
      DoFTools::make_hanging_node_constraints (saturation_dof_handler, saturation_constraints);
      saturation_constraints.close ();
    }
    {
      darcy_preconditioner_constraints.clear ();

      std::vector<bool> component_mask (dim+1, false);
      component_mask[dim] = true;


      DoFTools::make_hanging_node_constraints (darcy_dof_handler, darcy_preconditioner_constraints);
      DoFTools::make_zero_boundary_constraints (darcy_dof_handler, darcy_preconditioner_constraints, component_mask);

      darcy_preconditioner_constraints.close ();
    }


    std::vector<unsigned int> darcy_dofs_per_block (2);
    DoFTools::count_dofs_per_block (darcy_dof_handler, darcy_dofs_per_block, darcy_block_component);
    const unsigned int n_u = darcy_dofs_per_block[0],
		       n_p = darcy_dofs_per_block[1],
		       n_s = saturation_dof_handler.n_dofs();

    std::cout << "Number of active cells: "
	      << triangulation.n_active_cells()
	      << " (on "
	      << triangulation.n_levels()
	      << " levels)"
	      << std::endl
	      << "Number of degrees of freedom: "
	      << n_u + n_p + n_s
	      << " (" << n_u << '+' << n_p << '+'<< n_s <<')'
	      << std::endl
	      << std::endl;

    {
      darcy_matrix.clear ();

      BlockCompressedSimpleSparsityPattern csp (2,2);

      csp.block(0,0).reinit (n_u, n_u);
      csp.block(0,1).reinit (n_u, n_p);
      csp.block(1,0).reinit (n_p, n_u);
      csp.block(1,1).reinit (n_p, n_p);

      csp.collect_sizes ();

      Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);

      for (unsigned int c=0; c<dim+1; ++c)
	for (unsigned int d=0; d<dim+1; ++d)
	  if (! ((c==dim) && (d==dim)))
	    coupling[c][d] = DoFTools::always;
	  else
	    coupling[c][d] = DoFTools::none;


      DoFTools::make_sparsity_pattern (darcy_dof_handler, coupling, csp,
				       darcy_constraints, false);

      darcy_matrix.reinit (csp);
    }

    {
      Amg_preconditioner.reset ();
      Mp_preconditioner.reset ();
      darcy_preconditioner_matrix.clear ();

      BlockCompressedSimpleSparsityPattern csp (2,2);

      csp.block(0,0).reinit (n_u, n_u);
      csp.block(0,1).reinit (n_u, n_p);
      csp.block(1,0).reinit (n_p, n_u);
      csp.block(1,1).reinit (n_p, n_p);

      csp.collect_sizes ();

      Table<2,DoFTools::Coupling> coupling (dim+1, dim+1);
      for (unsigned int c=0; c<dim+1; ++c)
	for (unsigned int d=0; d<dim+1; ++d)
	  if (c == d)
	    coupling[c][d] = DoFTools::always;
	  else
	    coupling[c][d] = DoFTools::none;

      DoFTools::make_sparsity_pattern (darcy_dof_handler, coupling, csp,
				       darcy_constraints, false);

      darcy_preconditioner_matrix.reinit (csp);
    }


    {
      saturation_matrix.clear ();

      CompressedSimpleSparsityPattern csp (n_s, n_s);

      DoFTools::make_sparsity_pattern (saturation_dof_handler, csp,
				       saturation_constraints, false);


      saturation_matrix.reinit (csp);
    }

    darcy_solution.reinit (2);
    darcy_solution.block(0).reinit (n_u);
    darcy_solution.block(1).reinit (n_p);
    darcy_solution.collect_sizes ();

    last_computed_darcy_solution.reinit (2);
    last_computed_darcy_solution.block(0).reinit (n_u);
    last_computed_darcy_solution.block(1).reinit (n_p);
    last_computed_darcy_solution.collect_sizes ();

    second_last_computed_darcy_solution.reinit (2);
    second_last_computed_darcy_solution.block(0).reinit (n_u);
    second_last_computed_darcy_solution.block(1).reinit (n_p);
    second_last_computed_darcy_solution.collect_sizes ();

    darcy_rhs.reinit (2);
    darcy_rhs.block(0).reinit (n_u);
    darcy_rhs.block(1).reinit (n_p);
    darcy_rhs.collect_sizes ();

    saturation_solution.reinit (n_s);
    old_saturation_solution.reinit (n_s);
    old_old_saturation_solution.reinit (n_s);

    saturation_matching_last_computed_darcy_solution.reinit (n_s);

    saturation_rhs.reinit (n_s);
  }


				   // @sect3{Assembling matrices and preconditioners}

				   // The next few functions are
				   // devoted to setting up the
				   // various system and
				   // preconditioner matrices and
				   // right hand sides that we have to
				   // deal with in this program.

				   // @sect4{TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner}

				   // This function assembles the matrix we use
				   // for preconditioning the Darcy system. What
				   // we need are a vector mass matrix weighted by
				   // $\left(\mathbf{K} \lambda_t\right)^{-1}$
				   // on the velocity components and a mass
				   // matrix weighted by $\left(\mathbf{K}
				   // \lambda_t\right)$ on the pressure
				   // component. We start by generating a
				   // quadrature object of appropriate order,
				   // the FEValues object that can give values
				   // and gradients at the quadrature points
				   // (together with quadrature weights). Next
				   // we create data structures for the cell
				   // matrix and the relation between local and
				   // global DoFs. The vectors phi_u and
				   // grad_phi_p are going to hold the values of
				   // the basis functions in order to faster
				   // build up the local matrices, as was
				   // already done in step-22. Before we start
				   // the loop over all active cells, we have to
				   // specify which components are pressure and
				   // which are velocity.
				   //
				   // The creation of the local matrix
				   // is rather simple. There are only
				   // a term weighted by
				   // $\left(\mathbf{K}
				   // \lambda_t\right)^{-1}$ (on the
				   // velocity) and a Laplace matrix
				   // weighted by $\left(\mathbf{K}
				   // \lambda_t\right)$ to be
				   // generated, so the creation of
				   // the local matrix is done in
				   // essentially two lines. Since the
				   // material model functions at the
				   // top of this file only provide
				   // the inverses of the permeability
				   // and mobility, we have to compute
				   // $\mathbf K$ and $\lambda_t$ by
				   // hand from the given values, once
				   // per quadrature point.
				   //
				   // Once the
				   // local matrix is ready (loop over
				   // rows and columns in the local
				   // matrix on each quadrature
				   // point), we get the local DoF
				   // indices and write the local
				   // information into the global
				   // matrix. We do this by directly
				   // applying the constraints
				   // (i.e. darcy_preconditioner_constraints)
				   // that takes care of hanging node
				   // and zero Dirichlet boundary
				   // condition constraints. By doing
				   // so, we don't have to do that
				   // afterwards, and we later don't
				   // have to use
				   // ConstraintMatrix::condense and
				   // MatrixTools::apply_boundary_values,
				   // both functions that would need
				   // to modify matrix and vector
				   // entries and so are difficult to
				   // write for the Trilinos classes
				   // where we don't immediately have
				   // access to individual memory
				   // locations.
  template <int dim>
  void
  TwoPhaseFlowProblem<dim>::assemble_darcy_preconditioner ()
  {
    std::cout << "   Rebuilding darcy preconditioner..." << std::endl;

    darcy_preconditioner_matrix = 0;

    const QGauss<dim> quadrature_formula(darcy_degree+2);
    FEValues<dim>     darcy_fe_values (darcy_fe, quadrature_formula,
				       update_JxW_values |
				       update_values |
				       update_gradients |
				       update_quadrature_points);
    FEValues<dim> saturation_fe_values (saturation_fe, quadrature_formula,
					update_values);

    const unsigned int   dofs_per_cell   = darcy_fe.dofs_per_cell;
    const unsigned int   n_q_points      = quadrature_formula.size();

    std::vector<Tensor<2,dim> >       k_inverse_values (n_q_points);

    std::vector<double>               old_saturation_values (n_q_points);

    FullMatrix<double>                local_matrix (dofs_per_cell, dofs_per_cell);
    std::vector<unsigned int>         local_dof_indices (dofs_per_cell);

    std::vector<Tensor<1,dim> > phi_u   (dofs_per_cell);
    std::vector<Tensor<1,dim> > grad_phi_p (dofs_per_cell);

    const FEValuesExtractors::Vector velocities (0);
    const FEValuesExtractors::Scalar pressure (dim);

    typename DoFHandler<dim>::active_cell_iterator
      cell = darcy_dof_handler.begin_active(),
      endc = darcy_dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator
      saturation_cell = saturation_dof_handler.begin_active();

    for (; cell!=endc; ++cell, ++saturation_cell)
      {
	darcy_fe_values.reinit (cell);
	saturation_fe_values.reinit (saturation_cell);

	local_matrix = 0;

	saturation_fe_values.get_function_values (old_saturation_solution, old_saturation_values);

	k_inverse.value_list (darcy_fe_values.get_quadrature_points(),
			      k_inverse_values);

	for (unsigned int q=0; q<n_q_points; ++q)
	  {
	    const double old_s = old_saturation_values[q];

	    const double        inverse_mobility = mobility_inverse(old_s,viscosity);
	    const double        mobility         = 1.0 / inverse_mobility;
	    const Tensor<2,dim> permeability     = invert(k_inverse_values[q]);

	    for (unsigned int k=0; k<dofs_per_cell; ++k)
	      {
		phi_u[k]       = darcy_fe_values[velocities].value (k,q);
		grad_phi_p[k]  = darcy_fe_values[pressure].gradient (k,q);
	      }

	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		{
		  local_matrix(i,j) += (k_inverse_values[q] * inverse_mobility *
					phi_u[i] * phi_u[j]
					+
					permeability * mobility *
					grad_phi_p[i] * grad_phi_p[j])
				       * darcy_fe_values.JxW(q);
		}
	  }

	cell->get_dof_indices (local_dof_indices);
	darcy_preconditioner_constraints.distribute_local_to_global (local_matrix,
								     local_dof_indices,
								     darcy_preconditioner_matrix);
      }
  }


				   // @sect4{TwoPhaseFlowProblem<dim>::build_darcy_preconditioner}

				   // After calling the above
				   // functions to assemble the
				   // preconditioner matrix, this
				   // function generates the inner
				   // preconditioners that are going
				   // to be used for the Schur
				   // complement block
				   // preconditioner. The
				   // preconditioners need to be
				   // regenerated at every saturation
				   // time step since they depend on
				   // the saturation $S$ that varies
				   // with time.
				   //
				   // In here, we set up the
				   // preconditioner for the
				   // velocity-velocity matrix
				   // $\mathbf{M}^{\mathbf{u}}$ and
				   // the Schur complement
				   // $\mathbf{S}$. As explained in
				   // the introduction, we are going
				   // to use an IC preconditioner
				   // based on the vector matrix
				   // $\mathbf{M}^{\mathbf{u}}$ and
				   // another based on the scalar
				   // Laplace matrix
				   // $\tilde\mathbf{S}^p$ (which is
				   // spectrally close to the Schur
				   // complement of the Darcy
				   // matrix). Usually, the
				   // TrilinosWrappers::PreconditionIC
				   // class can be seen as a good
				   // black-box preconditioner which
				   // does not need any special
				   // knowledge of the matrix
				   // structure and/or the operator
				   // that's behind it.
  template <int dim>
  void
  TwoPhaseFlowProblem<dim>::build_darcy_preconditioner ()
  {
    assemble_darcy_preconditioner ();

    Amg_preconditioner = std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionIC>
			 (new TrilinosWrappers::PreconditionIC());
    Amg_preconditioner->initialize(darcy_preconditioner_matrix.block(0,0));

    Mp_preconditioner = std_cxx1x::shared_ptr<TrilinosWrappers::PreconditionIC>
			(new TrilinosWrappers::PreconditionIC());
    Mp_preconditioner->initialize(darcy_preconditioner_matrix.block(1,1));

  }


				   // @sect4{TwoPhaseFlowProblem<dim>::assemble_darcy_system}

				   // This is the function that assembles the
				   // linear system for the Darcy system.
				   //
				   // Regarding the technical details of
				   // implementation, the procedures are similar
				   // to those in step-22 and step-31. We reset
				   // matrix and vector, create a quadrature
				   // formula on the cells, and then create the
				   // respective FEValues object.
				   //
				   // There is one thing that needs to be
				   // commented: since we have a separate
				   // finite element and DoFHandler for the
				   // saturation, we need to generate a second
				   // FEValues object for the proper evaluation
				   // of the saturation solution. This isn't too
				   // complicated to realize here: just use the
				   // saturation structures and set an update
				   // flag for the basis function values which
				   // we need for evaluation of the saturation
				   // solution. The only important part to
				   // remember here is that the same quadrature
				   // formula is used for both FEValues objects
				   // to ensure that we get matching information
				   // when we loop over the quadrature points of
				   // the two objects.
				   //
				   // The declarations proceed with some
				   // shortcuts for array sizes, the creation of
				   // the local matrix, right hand side as well
				   // as the vector for the indices of the local
				   // dofs compared to the global system.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::assemble_darcy_system ()
  {
    darcy_matrix = 0;
    darcy_rhs    = 0;

    QGauss<dim>   quadrature_formula(darcy_degree+2);
    QGauss<dim-1> face_quadrature_formula(darcy_degree+2);

    FEValues<dim> darcy_fe_values (darcy_fe, quadrature_formula,
				   update_values    | update_gradients |
				   update_quadrature_points  | update_JxW_values);

    FEValues<dim> saturation_fe_values (saturation_fe, quadrature_formula,
					update_values);

    FEFaceValues<dim> darcy_fe_face_values (darcy_fe, face_quadrature_formula,
					    update_values    | update_normal_vectors |
					    update_quadrature_points  | update_JxW_values);

    const unsigned int   dofs_per_cell   = darcy_fe.dofs_per_cell;

    const unsigned int   n_q_points      = quadrature_formula.size();
    const unsigned int   n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    const PressureRightHandSide<dim>  pressure_right_hand_side;
    const PressureBoundaryValues<dim> pressure_boundary_values;

    std::vector<double>               pressure_rhs_values (n_q_points);
    std::vector<double>               boundary_values (n_face_q_points);
    std::vector<Tensor<2,dim> >       k_inverse_values (n_q_points);

				     // Next we need a vector that
				     // will contain the values of the
				     // saturation solution at the
				     // previous time level at the
				     // quadrature points to assemble
				     // the saturation dependent
				     // coefficients in the Darcy
				     // equations.
				     //
				     // The set of vectors we create
				     // next hold the evaluations of
				     // the basis functions as well as
				     // their gradients that will be
				     // used for creating the
				     // matrices. Putting these into
				     // their own arrays rather than
				     // asking the FEValues object for
				     // this information each time it
				     // is needed is an optimization
				     // to accelerate the assembly
				     // process, see step-22 for
				     // details.
				     //
				     // The last two declarations are used to
				     // extract the individual blocks (velocity,
				     // pressure, saturation) from the total FE
				     // system.
    std::vector<double>               old_saturation_values (n_q_points);

    std::vector<Tensor<1,dim> >       phi_u (dofs_per_cell);
    std::vector<double>               div_phi_u (dofs_per_cell);
    std::vector<double>               phi_p (dofs_per_cell);

    const FEValuesExtractors::Vector  velocities (0);
    const FEValuesExtractors::Scalar  pressure (dim);

				     // Now start the loop over all
				     // cells in the problem. We are
				     // working on two different
				     // DoFHandlers for this assembly
				     // routine, so we must have two
				     // different cell iterators for
				     // the two objects in use. This
				     // might seem a bit peculiar, but
				     // since both the Darcy system
				     // and the saturation system use
				     // the same grid we can assume
				     // that the two iterators run in
				     // sync over the cells of the two
				     // DoFHandler objects.
				     //
				     // The first statements within
				     // the loop are again all very
				     // familiar, doing the update of
				     // the finite element data as
				     // specified by the update flags,
				     // zeroing out the local arrays
				     // and getting the values of the
				     // old solution at the quadrature
				     // points.  At this point we also
				     // have to get the values of the
				     // saturation function of the
				     // previous time step at the
				     // quadrature points. To this
				     // end, we can use the
				     // FEValues::get_function_values
				     // (previously already used in
				     // step-9, step-14 and step-15),
				     // a function that takes a
				     // solution vector and returns a
				     // list of function values at the
				     // quadrature points of the
				     // present cell. In fact, it
				     // returns the complete
				     // vector-valued solution at each
				     // quadrature point, i.e. not
				     // only the saturation but also
				     // the velocities and pressure.
				     //
				     // Then we are ready to loop over
				     // the quadrature points on the
				     // cell to do the
				     // integration. The formula for
				     // this follows in a
				     // straightforward way from what
				     // has been discussed in the
				     // introduction.
				     //
				     // Once this is done, we start the loop over
				     // the rows and columns of the local matrix
				     // and feed the matrix with the relevant
				     // products.
				     //
				     // The last step in the loop over all cells
				     // is to enter the local contributions into
				     // the global matrix and vector structures to
				     // the positions specified in
				     // local_dof_indices. Again, we let the
				     // ConstraintMatrix class do the insertion of
				     // the cell matrix elements to the global
				     // matrix, which already condenses the
				     // hanging node constraints.
    typename DoFHandler<dim>::active_cell_iterator
      cell = darcy_dof_handler.begin_active(),
      endc = darcy_dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator
      saturation_cell = saturation_dof_handler.begin_active();

    for (; cell!=endc; ++cell, ++saturation_cell)
      {
	darcy_fe_values.reinit (cell);
	saturation_fe_values.reinit (saturation_cell);

	local_matrix = 0;
	local_rhs = 0;

	saturation_fe_values.get_function_values (old_saturation_solution, old_saturation_values);

	pressure_right_hand_side.value_list (darcy_fe_values.get_quadrature_points(),
					     pressure_rhs_values);
	k_inverse.value_list (darcy_fe_values.get_quadrature_points(),
			      k_inverse_values);

	for (unsigned int q=0; q<n_q_points; ++q)
	  {
	    for (unsigned int k=0; k<dofs_per_cell; ++k)
	      {
		phi_u[k]     = darcy_fe_values[velocities].value (k,q);
		div_phi_u[k] = darcy_fe_values[velocities].divergence (k,q);
		phi_p[k]     = darcy_fe_values[pressure].value (k,q);
	      }
	    for (unsigned int i=0; i<dofs_per_cell; ++i)
	      {
		const double old_s = old_saturation_values[q];
		for (unsigned int j=0; j<=i; ++j)
		  {
		    local_matrix(i,j) += (phi_u[i] * k_inverse_values[q] *
					  mobility_inverse(old_s,viscosity) * phi_u[j]
					  - div_phi_u[i] * phi_p[j]
					  - phi_p[i] * div_phi_u[j])
					 * darcy_fe_values.JxW(q);
		  }

		local_rhs(i) += (-phi_p[i] * pressure_rhs_values[q])*
				darcy_fe_values.JxW(q);
	      }
	  }

	for (unsigned int face_no=0;
	     face_no<GeometryInfo<dim>::faces_per_cell;
	     ++face_no)
	  if (cell->at_boundary(face_no))
	    {
	      darcy_fe_face_values.reinit (cell, face_no);

	      pressure_boundary_values
		.value_list (darcy_fe_face_values.get_quadrature_points(),
			     boundary_values);

	      for (unsigned int q=0; q<n_face_q_points; ++q)
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  {
		    const Tensor<1,dim>
		      phi_i_u = darcy_fe_face_values[velocities].value (i, q);

		    local_rhs(i) += -(phi_i_u *
				      darcy_fe_face_values.normal_vector(q) *
				      boundary_values[q] *
				      darcy_fe_face_values.JxW(q));
		  }
	    }

	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=i+1; j<dofs_per_cell; ++j)
	    local_matrix(i,j) = local_matrix(j,i);

	cell->get_dof_indices (local_dof_indices);

	darcy_constraints.distribute_local_to_global (local_matrix,
						      local_rhs,
						      local_dof_indices,
						      darcy_matrix,
						      darcy_rhs);

      }
  }


				   // @sect4{TwoPhaseFlowProblem<dim>::assemble_saturation_system}

				   // This function is to assemble the linear
				   // system for the saturation transport
				   // equation. It calls, if necessary, two
				   // other member functions:
				   // assemble_saturation_matrix() and
				   // assemble_saturation_rhs(). The former
				   // function then assembles the saturation
				   // matrix that only needs to be changed
				   // occasionally. On the other hand, the
				   // latter function that assembles the right
				   // hand side must be called at every
				   // saturation time step.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::assemble_saturation_system ()
  {
    if (rebuild_saturation_matrix == true)
      {
	saturation_matrix = 0;
	assemble_saturation_matrix ();
      }

    saturation_rhs = 0;
    assemble_saturation_rhs ();
  }



				   // @sect4{TwoPhaseFlowProblem<dim>::assemble_saturation_matrix}

				   // This function is easily understood since
				   // it only forms a simple mass matrix for the
				   // left hand side of the saturation linear
				   // system by basis functions phi_i_s and
				   // phi_j_s only. Finally, as usual, we enter
				   // the local contribution into the global
				   // matrix by specifying the position in
				   // local_dof_indices. This is done by letting
				   // the ConstraintMatrix class do the
				   // insertion of the cell matrix elements to
				   // the global matrix, which already condenses
				   // the hanging node constraints.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::assemble_saturation_matrix ()
  {
    QGauss<dim> quadrature_formula(saturation_degree+2);

    FEValues<dim> saturation_fe_values (saturation_fe, quadrature_formula,
					update_values | update_JxW_values);

    const unsigned int dofs_per_cell = saturation_fe.dofs_per_cell;

    const unsigned int n_q_points = quadrature_formula.size();

    FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
    Vector<double>       local_rhs (dofs_per_cell);

    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    typename DoFHandler<dim>::active_cell_iterator
      cell = saturation_dof_handler.begin_active(),
      endc = saturation_dof_handler.end();
    for (; cell!=endc; ++cell)
      {
	saturation_fe_values.reinit (cell);
	local_matrix = 0;
	local_rhs    = 0;

	for (unsigned int q=0; q<n_q_points; ++q)
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      const double phi_i_s = saturation_fe_values.shape_value (i,q);
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		{
		  const double phi_j_s = saturation_fe_values.shape_value (j,q);
		  local_matrix(i,j) += porosity * phi_i_s * phi_j_s * saturation_fe_values.JxW(q);
		}
	    }
	cell->get_dof_indices (local_dof_indices);

	saturation_constraints.distribute_local_to_global (local_matrix,
							   local_dof_indices,
							   saturation_matrix);

      }
  }



				   // @sect4{TwoPhaseFlowProblem<dim>::assemble_saturation_rhs}

				   // This function is to assemble the right
				   // hand side of the saturation transport
				   // equation. Before going about it, we have to
				   // create two FEValues objects for the Darcy
				   // and saturation systems respectively and,
				   // in addition, two FEFaceValues objects for
				   // the two systems because we have a
				   // boundary integral term in the weak form of
				   // saturation equation. For the FEFaceValues
				   // object of the saturation system, we also
				   // require normal vectors, which we request
				   // using the update_normal_vectors flag.
				   //
				   // Next, before looping over all the cells,
				   // we have to compute some parameters
				   // (e.g. global_u_infty, global_S_variation,
				   // and global_Omega_diameter) that the
				   // artificial viscosity $\nu$ needs. This is
				   // largely the same as was done in
				   // step-31, so you may see there for more
				   // information.
				   //
				   // The real works starts with the loop over all the
				   // saturation and Darcy cells to put the
				   // local contributions into the global
				   // vector. In this loop, in order to simplify
				   // the implementation, we split some of the
				   // work into two helper functions:
				   // assemble_saturation_rhs_cell_term and
				   // assemble_saturation_rhs_boundary_term.
				   // We note that we insert cell or boundary
				   // contributions into the global vector in
				   // the two functions rather than in this
				   // present function.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::assemble_saturation_rhs ()
  {
    QGauss<dim>   quadrature_formula(saturation_degree+2);
    QGauss<dim-1> face_quadrature_formula(saturation_degree+2);

    FEValues<dim> saturation_fe_values                   (saturation_fe, quadrature_formula,
							  update_values    | update_gradients |
							  update_quadrature_points  | update_JxW_values);
    FEValues<dim> darcy_fe_values                        (darcy_fe, quadrature_formula,
							  update_values);
    FEFaceValues<dim> saturation_fe_face_values          (saturation_fe, face_quadrature_formula,
							  update_values    | update_normal_vectors |
							  update_quadrature_points  | update_JxW_values);
    FEFaceValues<dim> darcy_fe_face_values               (darcy_fe, face_quadrature_formula,
							  update_values);
    FEFaceValues<dim> saturation_fe_face_values_neighbor (saturation_fe, face_quadrature_formula,
							  update_values);

    const unsigned int dofs_per_cell = saturation_dof_handler.get_fe().dofs_per_cell;
    std::vector<unsigned int> local_dof_indices (dofs_per_cell);

    const double global_max_u_F_prime = get_max_u_F_prime ();
    const std::pair<double,double>
      global_S_range = get_extrapolated_saturation_range ();
    const double global_S_variation = global_S_range.second - global_S_range.first;

    typename DoFHandler<dim>::active_cell_iterator
      cell = saturation_dof_handler.begin_active(),
      endc = saturation_dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator
      darcy_cell = darcy_dof_handler.begin_active();
    for (; cell!=endc; ++cell, ++darcy_cell)
      {
	saturation_fe_values.reinit (cell);
	darcy_fe_values.reinit (darcy_cell);

	cell->get_dof_indices (local_dof_indices);

	assemble_saturation_rhs_cell_term (saturation_fe_values,
					   darcy_fe_values,
					   global_max_u_F_prime,
					   global_S_variation,
					   local_dof_indices);

	for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
	     ++face_no)
	  {

	    if (cell->at_boundary(face_no))
	      {
		darcy_fe_face_values.reinit (darcy_cell, face_no);
		saturation_fe_face_values.reinit (cell, face_no);
		assemble_saturation_rhs_boundary_term (saturation_fe_face_values,
						       darcy_fe_face_values,
						       local_dof_indices);
	      }
	  }
      }
  }



				   // @sect4{TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_cell_term}

				   // This function takes care of integrating
				   // the cell terms of the right hand side of
				   // the saturation equation, and then
				   // assembling it into the global right hand
				   // side vector. Given the discussion in the
				   // introduction, the form of these
				   // contributions is clear. The only tricky
				   // part is getting the artificial viscosity
				   // and all that is necessary to compute
				   // it. The first half of the function is
				   // devoted to this task.
				   //
				   // The last part of the function is copying
				   // the local contributions into the global
				   // vector with position specified in
				   // local_dof_indices.
  template <int dim>
  void
  TwoPhaseFlowProblem<dim>::
  assemble_saturation_rhs_cell_term (const FEValues<dim>             &saturation_fe_values,
				     const FEValues<dim>             &darcy_fe_values,
				     const double                     global_max_u_F_prime,
				     const double                     global_S_variation,
				     const std::vector<unsigned int> &local_dof_indices)
  {
    const unsigned int dofs_per_cell = saturation_fe_values.dofs_per_cell;
    const unsigned int n_q_points    = saturation_fe_values.n_quadrature_points;

    std::vector<double>          old_saturation_solution_values(n_q_points);
    std::vector<double>          old_old_saturation_solution_values(n_q_points);
    std::vector<Tensor<1,dim> >  old_grad_saturation_solution_values(n_q_points);
    std::vector<Tensor<1,dim> >  old_old_grad_saturation_solution_values(n_q_points);
    std::vector<Vector<double> > present_darcy_solution_values(n_q_points, Vector<double>(dim+1));

    saturation_fe_values.get_function_values (old_saturation_solution, old_saturation_solution_values);
    saturation_fe_values.get_function_values (old_old_saturation_solution, old_old_saturation_solution_values);
    saturation_fe_values.get_function_grads (old_saturation_solution, old_grad_saturation_solution_values);
    saturation_fe_values.get_function_grads (old_old_saturation_solution, old_old_grad_saturation_solution_values);
    darcy_fe_values.get_function_values (darcy_solution, present_darcy_solution_values);

    const double nu
      = compute_viscosity (old_saturation_solution_values,
			   old_old_saturation_solution_values,
			   old_grad_saturation_solution_values,
			   old_old_grad_saturation_solution_values,
			   present_darcy_solution_values,
			   global_max_u_F_prime,
			   global_S_variation,
			   saturation_fe_values.get_cell()->diameter(),
			   old_time_step,
			   viscosity,
			   porosity);

    Vector<double> local_rhs (dofs_per_cell);

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  const double old_s = old_saturation_solution_values[q];
	  Tensor<1,dim> present_u;
	  for (unsigned int d=0; d<dim; ++d)
	    present_u[d] = present_darcy_solution_values[q](d);

	  const double        phi_i_s      = saturation_fe_values.shape_value (i, q);
	  const Tensor<1,dim> grad_phi_i_s = saturation_fe_values.shape_grad (i, q);

	  local_rhs(i) += (time_step *
			   fractional_flow(old_s,viscosity) *
			   present_u *
			   grad_phi_i_s
			   -
			   time_step *
			   nu *
			   old_grad_saturation_solution_values[q] * grad_phi_i_s
			   +
			   porosity * old_s * phi_i_s)
			  *
			  saturation_fe_values.JxW(q);
	}

    saturation_constraints.distribute_local_to_global (local_rhs,
						       local_dof_indices,
						       saturation_rhs);
  }


				   // @sect4{TwoPhaseFlowProblem<dim>::assemble_saturation_rhs_boundary_term}

				   // The next function is responsible for the
				   // boundary integral terms in the right
				   // hand side form of the saturation
				   // equation.  For these, we have to compute
				   // the upwinding flux on the global
				   // boundary faces, i.e. we impose Dirichlet
				   // boundary conditions weakly only on
				   // inflow parts of the global boundary. As
				   // before, this has been described in
				   // step-21 so we refrain from giving more
				   // descriptions about that.
  template <int dim>
  void
  TwoPhaseFlowProblem<dim>::
  assemble_saturation_rhs_boundary_term (const FEFaceValues<dim>             &saturation_fe_face_values,
					 const FEFaceValues<dim>             &darcy_fe_face_values,
					 const std::vector<unsigned int>     &local_dof_indices)
  {
    const unsigned int dofs_per_cell      = saturation_fe_face_values.dofs_per_cell;
    const unsigned int n_face_q_points    = saturation_fe_face_values.n_quadrature_points;

    Vector<double> local_rhs (dofs_per_cell);

    std::vector<double>          old_saturation_solution_values_face(n_face_q_points);
    std::vector<Vector<double> > present_darcy_solution_values_face(n_face_q_points,
								    Vector<double>(dim+1));
    std::vector<double>          neighbor_saturation (n_face_q_points);

    saturation_fe_face_values.get_function_values (old_saturation_solution,
						   old_saturation_solution_values_face);
    darcy_fe_face_values.get_function_values (darcy_solution,
					      present_darcy_solution_values_face);

    SaturationBoundaryValues<dim> saturation_boundary_values;
    saturation_boundary_values
      .value_list (saturation_fe_face_values.get_quadrature_points(),
		   neighbor_saturation);

    for (unsigned int q=0; q<n_face_q_points; ++q)
      {
	Tensor<1,dim> present_u_face;
	for (unsigned int d=0; d<dim; ++d)
	  present_u_face[d] = present_darcy_solution_values_face[q](d);

	const double normal_flux = present_u_face *
				   saturation_fe_face_values.normal_vector(q);

	const bool is_outflow_q_point = (normal_flux >= 0);

	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  local_rhs(i) -= time_step *
			  normal_flux *
			  fractional_flow((is_outflow_q_point == true
					   ?
					   old_saturation_solution_values_face[q]
					   :
					   neighbor_saturation[q]),
					  viscosity) *
			  saturation_fe_face_values.shape_value (i,q) *
			  saturation_fe_face_values.JxW(q);
      }
    saturation_constraints.distribute_local_to_global (local_rhs,
						       local_dof_indices,
						       saturation_rhs);
  }


				   // @sect3{TwoPhaseFlowProblem<dim>::solve}

				   // This function implements the operator
				   // splitting algorithm, i.e. in each time
				   // step it either re-computes the solution
				   // of the Darcy system or extrapolates
				   // velocity/pressure from previous time
				   // steps, then determines the size of the
				   // time step, and then updates the
				   // saturation variable. The implementation
				   // largely follows similar code in
				   // step-31. It is, next to the run()
				   // function, the central one in this
				   // program.
				   //
				   // At the beginning of the function, we ask
				   // whether to solve the pressure-velocity
				   // part by evaluating the posteriori
				   // criterion (see the following
				   // function). If necessary, we will solve
				   // the pressure-velocity part using the
				   // GMRES solver with the Schur complement
				   // block preconditioner as is described in
				   // the introduction.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::solve ()
  {
    const bool
      solve_for_pressure_and_velocity = determine_whether_to_solve_for_pressure_and_velocity ();

    if (solve_for_pressure_and_velocity == true)
      {
	std::cout << "   Solving Darcy (pressure-velocity) system..." << std::endl;

	assemble_darcy_system ();
	build_darcy_preconditioner ();

	{
	  const LinearSolvers::InverseMatrix<TrilinosWrappers::SparseMatrix,
	    TrilinosWrappers::PreconditionIC>
	    mp_inverse (darcy_preconditioner_matrix.block(1,1), *Mp_preconditioner);

	  const LinearSolvers::BlockSchurPreconditioner<TrilinosWrappers::PreconditionIC,
							TrilinosWrappers::PreconditionIC>
	    preconditioner (darcy_matrix, mp_inverse, *Amg_preconditioner);

	  SolverControl solver_control (darcy_matrix.m(),
					1e-16*darcy_rhs.l2_norm());

	  SolverGMRES<TrilinosWrappers::BlockVector>
	    gmres (solver_control,
		   SolverGMRES<TrilinosWrappers::BlockVector >::AdditionalData(100));

	  for (unsigned int i=0; i<darcy_solution.size(); ++i)
	    if (darcy_constraints.is_constrained(i))
	      darcy_solution(i) = 0;

	  gmres.solve(darcy_matrix, darcy_solution, darcy_rhs, preconditioner);

	  darcy_constraints.distribute (darcy_solution);

	  std::cout << "     "
		    << solver_control.last_step()
		    << " GMRES iterations for Darcy (pressure-velocity) system."
		    << std::endl;
	}

	{
	  second_last_computed_darcy_solution              = last_computed_darcy_solution;
	  last_computed_darcy_solution                     = darcy_solution;

	  saturation_matching_last_computed_darcy_solution = saturation_solution;
	}
      }
				     // On the other hand, if we have decided
				     // that we don't want to compute the
				     // solution of the Darcy system for the
				     // current time step, then we need to
				     // simply extrapolate the previous two
				     // Darcy solutions to the same time as we
				     // would have computed the
				     // velocity/pressure at. We do a simple
				     // linear extrapolation, i.e. given the
				     // current length $dt$ of the macro time
				     // step from the time when we last
				     // computed the Darcy solution to now
				     // (given by
				     // <code>current_macro_time_step</code>),
				     // and $DT$ the length of the last macro
				     // time step (given by
				     // <code>old_macro_time_step</code>),
				     // then we get
				     // $u^\ast = u_p + dt \frac{u_p-u_{pp}}{DT}
				     // = (1+dt/DT)u_p - dt/DT u_{pp}$, where
				     // $u_p$ and $u_{pp}$ are the last two
				     // computed Darcy solutions. We can
				     // implement this formula using just
				     // two lines of code.
				     //
				     // Note that the algorithm here only
				     // works if we have at least two
				     // previously computed Darcy solutions
				     // from which we can extrapolate to the
				     // current time, and this is ensured by
				     // requiring re-computation of the Darcy
				     // solution for the first 2 time steps.
    else
      {
	darcy_solution = last_computed_darcy_solution;
	darcy_solution.sadd (1 + current_macro_time_step / old_macro_time_step,
			     -current_macro_time_step / old_macro_time_step,
			     second_last_computed_darcy_solution);
      }


				     // With the so computed velocity
				     // vector, compute the optimal
				     // time step based on the CFL
				     // criterion discussed in the
				     // introduction...
    old_time_step = time_step;
    time_step = porosity *
		GridTools::minimal_cell_diameter(triangulation) /
		get_max_u_F_prime() / 12;


				     // ...and then also update the
				     // length of the macro time steps
				     // we use while we're dealing
				     // with time step sizes. In
				     // particular, this involves: (i)
				     // If we have just recomputed the
				     // Darcy solution, then the
				     // length of the previous macro
				     // time step is now fixed and the
				     // length of the current macro
				     // time step is, up to now,
				     // simply the length of the
				     // current (micro) time
				     // step. (ii) If we have not
				     // recomputed the Darcy solution,
				     // then the length of the current
				     // macro time step has just grown
				     // by <code>time_step</code>.
    if (solve_for_pressure_and_velocity == true)
      {
	old_macro_time_step     = current_macro_time_step;
	current_macro_time_step = time_step;
      }
    else
      current_macro_time_step += time_step;

				     // The last step in this function
				     // is to recompute the saturation
				     // solution based on the velocity
				     // field we've just
				     // obtained. This naturally
				     // happens in every time step,
				     // and we don't skip any of these
				     // computations. At the end of
				     // computing the saturation, we
				     // project back into the allowed
				     // interval $[0,1]$ to make sure
				     // our solution remains physical.
    {
      std::cout << "   Solving saturation transport equation..." << std::endl;

      assemble_saturation_system ();

      SolverControl solver_control (saturation_matrix.m(),
				    1e-16*saturation_rhs.l2_norm());
      SolverCG<TrilinosWrappers::Vector> cg (solver_control);

      TrilinosWrappers::PreconditionIC preconditioner;
      preconditioner.initialize (saturation_matrix);

      cg.solve (saturation_matrix, saturation_solution,
		saturation_rhs, preconditioner);

      saturation_constraints.distribute (saturation_solution);
      project_back_saturation ();

      std::cout << "     "
		<< solver_control.last_step()
		<< " CG iterations for saturation."
		<< std::endl;
    }
  }


				   // @sect3{Tool functions}

				   // @sect4{TwoPhaseFlowProblem<dim>::determine_whether_to_solve_for_pressure_and_velocity}

				   // This function is to implement the a
				   // posteriori criterion for
				   // adaptive operator splitting. As mentioned
				   // in step-31, we use two FEValues objects
				   // initialized with two cell iterators that
				   // we walk in parallel through the two
				   // DoFHandler objects associated with the
				   // same Triangulation object; for these two
				   // FEValues objects, we use of course the
				   // same quadrature objects so that we can
				   // iterate over the same set of quadrature
				   // points, but each FEValues object will get
				   // update flags only according to what it
				   // actually needs to compute.
				   //
				   // In addition to this, if someone doesn't
				   // want to perform their simulation with
				   // operator splitting, they can lower the
				   // criterion value (default value is $5.0$)
				   // down to zero ad therefore numerical
				   // algorithm becomes the original IMPES
				   // method.
  template <int dim>
  bool
  TwoPhaseFlowProblem<dim>::determine_whether_to_solve_for_pressure_and_velocity () const
  {
    if (timestep_number <= 2)
      return true;

    const QGauss<dim> quadrature_formula(saturation_degree+2);
    const unsigned int   n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (saturation_fe, quadrature_formula,
			     update_values | update_quadrature_points);

    std::vector<double> old_saturation_after_solving_pressure (n_q_points);
    std::vector<double> present_saturation (n_q_points);

    std::vector<Tensor<2,dim> >       k_inverse_values (n_q_points);

    double max_global_aop_indicator = 0.0;

    typename DoFHandler<dim>::active_cell_iterator
      cell = saturation_dof_handler.begin_active(),
      endc = saturation_dof_handler.end();
    for (; cell!=endc; ++cell)
      {
	double max_local_mobility_reciprocal_difference = 0.0;
	double max_local_permeability_inverse_l1_norm = 0.0;

	fe_values.reinit(cell);
	fe_values.get_function_values (saturation_matching_last_computed_darcy_solution,
				       old_saturation_after_solving_pressure);
	fe_values.get_function_values (saturation_solution,
				       present_saturation);

	k_inverse.value_list (fe_values.get_quadrature_points(),
			      k_inverse_values);

	for (unsigned int q=0; q<n_q_points; ++q)
	  {
	    double mobility_reciprocal_difference = std::fabs( mobility_inverse(present_saturation[q],viscosity)
							       -
							       mobility_inverse(old_saturation_after_solving_pressure[q],viscosity) );

	    max_local_mobility_reciprocal_difference = std::max(max_local_mobility_reciprocal_difference,
								mobility_reciprocal_difference);

	    max_local_permeability_inverse_l1_norm = std::max(max_local_permeability_inverse_l1_norm,
							      k_inverse_values[q][0][0]);
	  }

	max_global_aop_indicator = std::max(max_global_aop_indicator,
					    (max_local_mobility_reciprocal_difference*max_local_permeability_inverse_l1_norm));
      }

    if ( max_global_aop_indicator > AOS_threshold )
      {
	return true;
      }
    else
      {
	std::cout << "   Activating adaptive operating splitting" << std::endl;
	return false;
      }
  }



				   // @sect3{TwoPhaseFlowProblem<dim>::compute_refinement_indicators}

				   // This function is to to compute the
				   // refinement indicator discussed in the
				   // introduction for each cell and its
				   // implementation is similar to that
				   // contained in step-33. There is no need to
				   // repeat descriptions about it.
  template <int dim>
  void
  TwoPhaseFlowProblem<dim>::
  compute_refinement_indicators (const TrilinosWrappers::Vector &predicted_saturation_solution,
				 Vector<double>                 &refinement_indicators) const
  {
    const QMidpoint<dim> quadrature_formula;
    FEValues<dim> fe_values (saturation_fe, quadrature_formula, update_gradients);
    std::vector<Tensor<1,dim> > grad_saturation (1);

    double max_refinement_indicator = 0.0;

    typename DoFHandler<dim>::active_cell_iterator
      cell = saturation_dof_handler.begin_active(),
      endc = saturation_dof_handler.end();
    for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
      {
	fe_values.reinit(cell);
	fe_values.get_function_grads (predicted_saturation_solution,
				      grad_saturation);

	refinement_indicators(cell_no) = grad_saturation[0].norm();
	max_refinement_indicator = std::max(max_refinement_indicator,
					    refinement_indicators(cell_no));
      }
  }



				   // @sect3{TwoPhaseFlowProblem<dim>::refine_mesh}

				   // This function is to decide if every cell
				   // is refined or coarsened with computed
				   // refinement indicators in the previous
				   // function and do the interpolations of the
				   // solution vectors. The main difference from
				   // the previous time-dependent tutorials is
				   // that there is no need to do the solution
				   // interpolations if we don't have any cell
				   // that is refined or coarsend, saving some
				   // additional computing time.
  template <int dim>
  void
  TwoPhaseFlowProblem<dim>::
  refine_mesh (const unsigned int    min_grid_level,
	       const unsigned int    max_grid_level,
	       const Vector<double> &refinement_indicators)
  {
				     //TODO: use a useful refinement criterion, in much the same way as we do in step-31
    {
      typename DoFHandler<dim>::active_cell_iterator
	cell = saturation_dof_handler.begin_active(),
	endc = saturation_dof_handler.end();

      for (unsigned int cell_no=0; cell!=endc; ++cell, ++cell_no)
	{
	  cell->clear_coarsen_flag();
	  cell->clear_refine_flag();

	  if ((static_cast<unsigned int>(cell->level()) < max_grid_level) &&
	      (std::fabs(refinement_indicators(cell_no)) > saturation_refinement_threshold))
	    cell->set_refine_flag();
	  else
	    if ((static_cast<unsigned int>(cell->level()) > min_grid_level) &&
		(std::fabs(refinement_indicators(cell_no)) < 0.75 * saturation_refinement_threshold))
	      cell->set_coarsen_flag();
	}
    }

    triangulation.prepare_coarsening_and_refinement ();

    unsigned int number_of_cells_refine  = 0;
    unsigned int number_of_cells_coarsen = 0;

    {
      typename DoFHandler<dim>::active_cell_iterator
	cell = saturation_dof_handler.begin_active(),
	endc = saturation_dof_handler.end();

      for (; cell!=endc; ++cell)
	if (cell->refine_flag_set())
	  ++number_of_cells_refine;
	else
	  if (cell->coarsen_flag_set())
	    ++number_of_cells_coarsen;
    }

    std::cout << "   "
	      << number_of_cells_refine
	      << " cell(s) are going to be refined."
	      << std::endl;
    std::cout << "   "
	      << number_of_cells_coarsen
	      << " cell(s) are going to be coarsened."
	      << std::endl;

    std::cout << std::endl;

    if ( number_of_cells_refine > 0 || number_of_cells_coarsen > 0 )
      {
	std::vector<TrilinosWrappers::Vector> x_saturation (3);
	x_saturation[0] = saturation_solution;
	x_saturation[1] = old_saturation_solution;
	x_saturation[2] = saturation_matching_last_computed_darcy_solution;

	std::vector<TrilinosWrappers::BlockVector> x_darcy (2);
	x_darcy[0] = last_computed_darcy_solution;
	x_darcy[1] = second_last_computed_darcy_solution;

	SolutionTransfer<dim,TrilinosWrappers::Vector> saturation_soltrans(saturation_dof_handler);

	SolutionTransfer<dim,TrilinosWrappers::BlockVector> darcy_soltrans(darcy_dof_handler);


	triangulation.prepare_coarsening_and_refinement();
	saturation_soltrans.prepare_for_coarsening_and_refinement(x_saturation);

	darcy_soltrans.prepare_for_coarsening_and_refinement(x_darcy);

	triangulation.execute_coarsening_and_refinement ();
	setup_dofs ();

	std::vector<TrilinosWrappers::Vector> tmp_saturation (3);
	tmp_saturation[0].reinit (saturation_solution);
	tmp_saturation[1].reinit (saturation_solution);
	tmp_saturation[2].reinit (saturation_solution);
	saturation_soltrans.interpolate(x_saturation, tmp_saturation);

	saturation_solution = tmp_saturation[0];
	old_saturation_solution = tmp_saturation[1];
	saturation_matching_last_computed_darcy_solution = tmp_saturation[2];

	std::vector<TrilinosWrappers::BlockVector> tmp_darcy (2);
	tmp_darcy[0].reinit (darcy_solution);
	tmp_darcy[1].reinit (darcy_solution);
	darcy_soltrans.interpolate(x_darcy, tmp_darcy);

	last_computed_darcy_solution = tmp_darcy[0];
	second_last_computed_darcy_solution = tmp_darcy[1];

	rebuild_saturation_matrix    = true;
      }
    else
      {
	rebuild_saturation_matrix    = false;

	std::vector<unsigned int> darcy_block_component (dim+1,0);
	darcy_block_component[dim] = 1;

	std::vector<unsigned int> darcy_dofs_per_block (2);
	DoFTools::count_dofs_per_block (darcy_dof_handler, darcy_dofs_per_block, darcy_block_component);
	const unsigned int n_u = darcy_dofs_per_block[0],
			   n_p = darcy_dofs_per_block[1],
			   n_s = saturation_dof_handler.n_dofs();

	std::cout << "Number of active cells: "
		  << triangulation.n_active_cells()
		  << " (on "
		  << triangulation.n_levels()
		  << " levels)"
		  << std::endl
		  << "Number of degrees of freedom: "
		  << n_u + n_p + n_s
		  << " (" << n_u << '+' << n_p << '+'<< n_s <<')'
		  << std::endl
		  << std::endl;
      }

  }



				   // @sect3{TwoPhaseFlowProblem<dim>::output_results}

				   // This function to process the output
				   // data. We only store the results when we
				   // actually solve the pressure and velocity
				   // part at the present time step. The rest of
				   // the implementation is similar to that
				   // output function in step-31, which
				   // implementations has been explained in that
				   // tutorial.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::output_results ()  const
  {
    const FESystem<dim> joint_fe (darcy_fe, 1,
				  saturation_fe, 1);
    DoFHandler<dim> joint_dof_handler (triangulation);
    joint_dof_handler.distribute_dofs (joint_fe);
    Assert (joint_dof_handler.n_dofs() ==
	    darcy_dof_handler.n_dofs() + saturation_dof_handler.n_dofs(),
	    ExcInternalError());

    Vector<double> joint_solution (joint_dof_handler.n_dofs());

    {
      std::vector<unsigned int> local_joint_dof_indices (joint_fe.dofs_per_cell);
      std::vector<unsigned int> local_darcy_dof_indices (darcy_fe.dofs_per_cell);
      std::vector<unsigned int> local_saturation_dof_indices (saturation_fe.dofs_per_cell);

      typename DoFHandler<dim>::active_cell_iterator
	joint_cell      = joint_dof_handler.begin_active(),
	joint_endc      = joint_dof_handler.end(),
	darcy_cell      = darcy_dof_handler.begin_active(),
	saturation_cell = saturation_dof_handler.begin_active();

      for (; joint_cell!=joint_endc; ++joint_cell, ++darcy_cell, ++saturation_cell)
	{
	  joint_cell->get_dof_indices (local_joint_dof_indices);
	  darcy_cell->get_dof_indices (local_darcy_dof_indices);
	  saturation_cell->get_dof_indices (local_saturation_dof_indices);

	  for (unsigned int i=0; i<joint_fe.dofs_per_cell; ++i)
	    if (joint_fe.system_to_base_index(i).first.first == 0)
	      {
		Assert (joint_fe.system_to_base_index(i).second
			<
			local_darcy_dof_indices.size(),
			ExcInternalError());
		joint_solution(local_joint_dof_indices[i])
		  = darcy_solution(local_darcy_dof_indices[joint_fe.system_to_base_index(i).second]);
	      }
	    else
	      {
		Assert (joint_fe.system_to_base_index(i).first.first == 1,
			ExcInternalError());
		Assert (joint_fe.system_to_base_index(i).second
			<
			local_darcy_dof_indices.size(),
			ExcInternalError());
		joint_solution(local_joint_dof_indices[i])
		  = saturation_solution(local_saturation_dof_indices[joint_fe.system_to_base_index(i).second]);
	      }

	}
    }
    std::vector<std::string> joint_solution_names;
    switch (dim)
      {
	case 2:
	      joint_solution_names.push_back ("u");
	      joint_solution_names.push_back ("v");
	      break;

	case 3:
	      joint_solution_names.push_back ("u");
	      joint_solution_names.push_back ("v");
	      joint_solution_names.push_back ("w");
	      break;

	default:
	      Assert (false, ExcNotImplemented());
      }
    joint_solution_names.push_back ("pressure");
    joint_solution_names.push_back ("saturation");

    std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation
      (dim, DataComponentInterpretation::component_is_part_of_vector);
    data_component_interpretation
      .push_back (DataComponentInterpretation::component_is_scalar);
    data_component_interpretation
      .push_back (DataComponentInterpretation::component_is_scalar);

    DataOut<dim> data_out;

    data_out.attach_dof_handler (joint_dof_handler);
    data_out.add_data_vector (joint_solution, joint_solution_names,
			      DataOut<dim>::type_dof_data,
			      data_component_interpretation);

    data_out.build_patches ();

    std::string filename = "solution-" +
			   Utilities::int_to_string (timestep_number, 5) + ".vtu";
    std::ofstream output (filename.c_str());
    data_out.write_vtu (output);
  }



				   // @sect3{TwoPhaseFlowProblem<dim>::THE_REMAINING_FUNCTIONS}

				   // The remaining functions that have been
				   // used in step-31 so we don't have to
				   // describe their implementations.
  template <int dim>
  void
  TwoPhaseFlowProblem<dim>::project_back_saturation ()
  {
    for (unsigned int i=0; i<saturation_solution.size(); ++i)
      if (saturation_solution(i) < 0)
	saturation_solution(i) = 0;
      else
	if (saturation_solution(i) > 1)
	  saturation_solution(i) = 1;
  }


  template <int dim>
  double
  TwoPhaseFlowProblem<dim>::get_max_u_F_prime () const
  {
    QGauss<dim>   quadrature_formula(darcy_degree+2);
    const unsigned int   n_q_points
      = quadrature_formula.size();

    FEValues<dim> darcy_fe_values (darcy_fe, quadrature_formula,
				   update_values);
    FEValues<dim> saturation_fe_values (saturation_fe, quadrature_formula,
					update_values);

    std::vector<Vector<double> > darcy_solution_values(n_q_points,
						       Vector<double>(dim+1));
    std::vector<double>          saturation_values (n_q_points);

    double max_velocity_times_dF_dS = 0;

    typename DoFHandler<dim>::active_cell_iterator
      cell = darcy_dof_handler.begin_active(),
      endc = darcy_dof_handler.end();
    typename DoFHandler<dim>::active_cell_iterator
      saturation_cell = saturation_dof_handler.begin_active();
    for (; cell!=endc; ++cell, ++saturation_cell)
      {
	darcy_fe_values.reinit (cell);
	saturation_fe_values.reinit (saturation_cell);

	darcy_fe_values.get_function_values (darcy_solution, darcy_solution_values);
	saturation_fe_values.get_function_values (old_saturation_solution, saturation_values);

	for (unsigned int q=0; q<n_q_points; ++q)
	  {
	    Tensor<1,dim> velocity;
	    for (unsigned int i=0; i<dim; ++i)
	      velocity[i] = darcy_solution_values[q](i);

	    const double dF_dS = fractional_flow_derivative(saturation_values[q],viscosity);

	    max_velocity_times_dF_dS = std::max (max_velocity_times_dF_dS,
						 velocity.norm() * dF_dS);
	  }
      }

    return max_velocity_times_dF_dS;
  }


  template <int dim>
  std::pair<double,double>
  TwoPhaseFlowProblem<dim>::get_extrapolated_saturation_range () const
  {
    const QGauss<dim>  quadrature_formula(saturation_degree+2);
    const unsigned int n_q_points = quadrature_formula.size();

    FEValues<dim> fe_values (saturation_fe, quadrature_formula,
			     update_values);
    std::vector<double> old_saturation_values(n_q_points);
    std::vector<double> old_old_saturation_values(n_q_points);

    if (timestep_number != 0)
      {
	double min_saturation = (1. + time_step/old_time_step) *
				old_saturation_solution.linfty_norm()
				+
				time_step/old_time_step *
				old_old_saturation_solution.linfty_norm(),
	       max_saturation = -min_saturation;

	typename DoFHandler<dim>::active_cell_iterator
	  cell = saturation_dof_handler.begin_active(),
	  endc = saturation_dof_handler.end();
	for (; cell!=endc; ++cell)
	  {
	    fe_values.reinit (cell);
	    fe_values.get_function_values (old_saturation_solution,
					   old_saturation_values);
	    fe_values.get_function_values (old_old_saturation_solution,
					   old_old_saturation_values);

	    for (unsigned int q=0; q<n_q_points; ++q)
	      {
		const double saturation =
		  (1. + time_step/old_time_step) * old_saturation_values[q]-
		  time_step/old_time_step * old_old_saturation_values[q];

		min_saturation = std::min (min_saturation, saturation);
		max_saturation = std::max (max_saturation, saturation);
	      }
	  }

	return std::make_pair(min_saturation, max_saturation);
      }
    else
      {
	double min_saturation = old_saturation_solution.linfty_norm(),
	       max_saturation = -min_saturation;

	typename DoFHandler<dim>::active_cell_iterator
	  cell = saturation_dof_handler.begin_active(),
	  endc = saturation_dof_handler.end();
	for (; cell!=endc; ++cell)
	  {
	    fe_values.reinit (cell);
	    fe_values.get_function_values (old_saturation_solution,
					   old_saturation_values);

	    for (unsigned int q=0; q<n_q_points; ++q)
	      {
		const double saturation = old_saturation_values[q];

		min_saturation = std::min (min_saturation, saturation);
		max_saturation = std::max (max_saturation, saturation);
	      }
	  }

	return std::make_pair(min_saturation, max_saturation);
      }
  }



  template <int dim>
  double
  TwoPhaseFlowProblem<dim>::
  compute_viscosity (const std::vector<double>          &old_saturation,
		     const std::vector<double>          &old_old_saturation,
		     const std::vector<Tensor<1,dim> >  &old_saturation_grads,
		     const std::vector<Tensor<1,dim> >  &old_old_saturation_grads,
		     const std::vector<Vector<double> > &present_darcy_values,
		     const double                        global_max_u_F_prime,
		     const double                        global_S_variation,
		     const double                        cell_diameter,
		     const double                        old_time_step,
		     const double                        viscosity,
		     const double                        porosity) const
  {
    const double beta = .35 * dim;
    const double alpha = 1;

    if (global_max_u_F_prime == 0)
      return 5e-3 * cell_diameter;

    const unsigned int n_q_points = old_saturation.size();

    double max_residual = 0;
    double max_velocity_times_dF_dS = 0;

    const bool use_dF_dS = true;

    for (unsigned int q=0; q < n_q_points; ++q)
      {
	Tensor<1,dim> u;
	for (unsigned int d=0; d<dim; ++d)
	  u[d] = present_darcy_values[q](d);

	const double dS_dt = porosity * (old_saturation[q] - old_old_saturation[q])
			     / old_time_step;

	const double dF_dS = fractional_flow_derivative ((old_saturation[q] + old_old_saturation[q]) / 2.0,viscosity);

	const double u_grad_S = u * dF_dS *
				(old_saturation_grads[q] + old_old_saturation_grads[q]) / 2.0;

	const double residual
	  = std::abs((dS_dt + u_grad_S) *
		     std::pow((old_saturation[q]+old_old_saturation[q]) / 2,
			      alpha-1.));

	max_residual = std::max (residual,        max_residual);
	max_velocity_times_dF_dS = std::max (std::sqrt (u*u) *
					     (use_dF_dS
					      ?
					      std::max(dF_dS,1.)
					      :
					      1),
					     max_velocity_times_dF_dS);
      }

    const double c_R = 0.0003;
    const double global_scaling = c_R * porosity * (global_max_u_F_prime) * global_S_variation /
				  std::pow(global_Omega_diameter, alpha - 2.);

// first order stabilization
//  return beta * max_velocity_times_dF_dS * cell_diameter / porosity;

// entropy diffusion stabilization

//  if (time > 1500.*24*3600)
//    std::cout << "xxx "
//	      << (beta *
//		  1./porosity *
//		  (max_velocity_times_dF_dS) *
//		  cell_diameter)
//	      << ' '
//	      << (beta *
//		  1./porosity *
//		  (max_velocity_times_dF_dS) *
//		  std::pow(cell_diameter,alpha) *
//		  max_residual / global_scaling)
//	      << std::endl;
    return (beta *
	    (max_velocity_times_dF_dS) *
	    std::min (cell_diameter,
		      std::pow(cell_diameter,alpha) *
		      max_residual / global_scaling));
  }


				   // @sect3{TwoPhaseFlowProblem<dim>::run}

				   // In this function, we follow the structure
				   // of the same function partly in step-21 and
				   // partly in step-31 so again there is no
				   // need to repeat it. However, since we
				   // consider the simulation with grid
				   // adaptivity, we need to compute a
				   // saturation predictor, which implementation
				   // was first used in step-33, for the
				   // function that computes the refinement
				   // indicators.
  template <int dim>
  void TwoPhaseFlowProblem<dim>::run ()
  {
    const unsigned int initial_refinement     = (dim == 2 ? 4 : 2);
    const unsigned int n_pre_refinement_steps = (dim == 2 ? 4 : 3);


    GridGenerator::hyper_cube (triangulation, 0, 1);
    triangulation.refine_global (initial_refinement);
    global_Omega_diameter = GridTools::diameter (triangulation);

    setup_dofs ();

    unsigned int pre_refinement_step = 0;

    start_time_iteration:

    VectorTools::project (saturation_dof_handler,
			  saturation_constraints,
			  QGauss<dim>(saturation_degree+2),
			  SaturationInitialValues<dim>(),
			  old_saturation_solution);

    timestep_number = 0;
    time_step = old_time_step = 0;
    current_macro_time_step = old_macro_time_step = 0;

    double time = 0;

    do
      {
	std::cout << "Timestep " << timestep_number
		  << ":  t=" << time
		  << ", dt=" << time_step
		  << std::endl;

	solve ();

	output_results ();

	{
					   // check if this already initializes the vector of if we need the next line
	  TrilinosWrappers::Vector predicted_saturation_solution (saturation_solution);
	  predicted_saturation_solution = saturation_solution;
	  predicted_saturation_solution.sadd (2.0, -1.0, old_saturation_solution);

					   // TODO: move this into refine_mesh
	  Vector<double> refinement_indicators (triangulation.n_active_cells());

	  compute_refinement_indicators(predicted_saturation_solution,
					refinement_indicators);
	  refine_mesh (initial_refinement,
		       initial_refinement + n_pre_refinement_steps,
		       refinement_indicators);
	}

	if ((timestep_number == 0) &&
	    (pre_refinement_step < n_pre_refinement_steps))
	  {
	    ++pre_refinement_step;
	    goto start_time_iteration;
	  }

	time += time_step;
	++timestep_number;

	old_old_saturation_solution = old_saturation_solution;
	old_saturation_solution = saturation_solution;
      }
    while (time <= 250);
  }
}



				 // @sect3{The <code>main</code> function}
				 //
				 // The main function looks almost the
				 // same as in all other programs. In
				 // particular, it is essentially the
				 // same as in step-31 where we also
				 // explain the need to initialize the
				 // MPI subsystem.
int main (int argc, char *argv[])
{
  try
    {
      using namespace dealii;
      using namespace Step43;

      deallog.depth_console (0);

      Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv);

      TwoPhaseFlowProblem<3> two_phase_flow_problem(1);
      two_phase_flow_problem.run ();
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
