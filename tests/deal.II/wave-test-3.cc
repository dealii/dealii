// deal_II_libraries.g=-ldeal_II_2d.g
// deal_II_libraries=-ldeal_II_2d

#include <grid/tria_boundary_lib.h>
#include <base/parameter_handler.h>
#include <base/forward_declarations.h>
#include <lac/forward_declarations.h>
#include <grid/forward_declarations.h>
#include <numerics/time_dependent.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_tools.h>

#include <fstream>
#include <string>

ofstream logfile("wave-test-3.output");

class UserMatrix;
class SweepInfo;
template <int dim> class SweepData;
template <int dim> class WaveParameters;
template <int dim> class TimeStep_Primal;
template <int dim> class TimeStep_Dual;
template <int dim> class DualFunctional;
template <int dim> class EvaluationBase;
template <int dim> class TimeStep_ErrorEstimation;
template <int dim> class TimeStep_Postprocess;



template <int dim>
class TimeStepBase_Wave :  public TimeStepBase_Tria<dim>{
  public:
    TimeStepBase_Wave ();
    TimeStepBase_Wave (const double                    time,
		       TimeStepBase_Tria<dim>::Flags   flags,
		       const WaveParameters<dim>      &parameters);
    const TimeStep_Primal<dim> & get_timestep_primal () const;
    const TimeStep_Dual<dim> &   get_timestep_dual () const;
    const TimeStep_Postprocess<dim> &   get_timestep_postprocess () const;
    string tmp_filename_base (const string &branch_signature) const;
    void attach_sweep_info (SweepInfo &sweep_info);
    void attach_sweep_data (SweepData<dim> &sweep_data);

  protected:
    const WaveParameters<dim> &parameters;
    SweepInfo *sweep_info;
    SweepData<dim> *sweep_data;
};



template <int dim>
class TimeStep_Wave :  public virtual TimeStepBase_Wave<dim>
{
  public:
    TimeStep_Wave (const string fe_name);
    ~TimeStep_Wave();
    virtual void wake_up (const unsigned int wakeup_level);
    virtual void sleep (const unsigned int sleep_level);
    virtual void end_sweep ();
    unsigned int solve (const UserMatrix       &matrix,
			Vector<double>         &solution,
			const Vector<double>   &rhs) const;
    virtual string branch_signature () const = 0;
    DeclException0 (ExcIO);
    DeclException0 (ExcCoarsestGridsDiffer);
    
    
  protected:
    struct StatisticData 
    {
	StatisticData ();
	StatisticData (const unsigned int        n_active_cells,
		       const unsigned int        n_dofs,
		       const unsigned int        n_solver_steps_helmholtz,
		       const unsigned int        n_solver_steps_projection,
		       const pair<double,double> energy);
	static void write_descriptions (ostream &out);
	void write (ostream &out) const;
	unsigned int n_active_cells;
	unsigned int n_dofs;
	unsigned int n_solver_steps_helmholtz;
	unsigned int n_solver_steps_projection;
	pair<double,double> energy;
    };	

    DoFHandler<dim>           *dof_handler;
    const FiniteElement<dim> &fe;
    const Quadrature<dim>    &quadrature;
    const Quadrature<dim-1>  &quadrature_face;
    ConstraintMatrix          constraints;
    SparsityPattern           system_sparsity;
    SparseMatrix<double>       mass_matrix, laplace_matrix;
    Vector<double>             u, v;
    StatisticData              statistic_data;
    void create_matrices ();
    void transfer_old_solutions (Vector<double> &old_u,
				 Vector<double> &old_v) const;
    void transfer_old_solutions (const typename DoFHandler<dim>::cell_iterator &old_cell,
				 const typename DoFHandler<dim>::cell_iterator &new_cell,
				 const Vector<double>  &old_grid_u,
				 const Vector<double>  &old_grid_v,
				 Vector<double>        &old_u,
				 Vector<double>        &old_v) const;
    pair<double,double> compute_energy ();
    template <int anydim> friend class DualFunctional;
    template <int anydim> friend class EvaluationBase;
    template <int anydim> friend class TimeStep_ErrorEstimation;
    template <int anydim> friend class TimeStep_Postprocess;
};



template <int dim>
class TimeStep_Primal :  public TimeStep_Wave<dim> 
{
  public:
    TimeStep_Primal (const string &primal_fe);
    void do_initial_step ();
    void do_timestep ();
    virtual void solve_primal_problem ();    
    virtual string branch_signature () const;
    virtual void wake_up (const unsigned int wakeup_level);
    
  private:
    void assemble_vectors (Vector<double> &right_hand_side1,
			   Vector<double> &right_hand_side2);
    void build_rhs (Vector<double> &right_hand_side1,
		    Vector<double> &right_hand_side2);
    void build_rhs (const typename DoFHandler<dim>::cell_iterator &old_cell,
		    const typename DoFHandler<dim>::cell_iterator &new_cell,
		    FEValues<dim>        &fe_values,
		    Vector<double>       &right_hand_side1,
		    Vector<double>       &right_hand_side2);
    unsigned int
    collect_from_children (const typename DoFHandler<dim>::cell_iterator &old_cell,
			   FEValues<dim>  &fe_values,
			   Vector<double>        &rhs1,
			   Vector<double>        &rhs2) const;
    unsigned int
    distribute_to_children (const typename DoFHandler<dim>::cell_iterator &cell,
			    FEValues<dim>  &fe_values,
			    const Vector<double>  &old_dof_values_u,
			    const Vector<double>  &old_dof_values_v,
			    Vector<double>        &right_hand_side1,
			    Vector<double>        &right_hand_side2);
};



template <int dim>
class TimeStep_Dual :  public TimeStep_Wave<dim>
{
  public:
    TimeStep_Dual (const string &dual_fe);
    void do_initial_step ();
    void do_timestep ();
    virtual void solve_dual_problem ();    
    virtual string branch_signature () const;
    virtual void wake_up (const unsigned int wakeup_level);

  private:
    void assemble_vectors (Vector<double> &right_hand_side1,
			   Vector<double> &right_hand_side2);
    void build_rhs (Vector<double> &right_hand_side1,
		    Vector<double> &right_hand_side2);
    void build_rhs (const typename DoFHandler<dim>::cell_iterator &old_cell,
		    const typename DoFHandler<dim>::cell_iterator &new_cell,
		    FEValues<dim>        &fe_values,
		    Vector<double>       &right_hand_side1,
		    Vector<double>       &right_hand_side2);
    unsigned int
    collect_from_children (const typename DoFHandler<dim>::cell_iterator &old_cell,
			   FEValues<dim>  &fe_values,
			   Vector<double>        &rhs1,
			   Vector<double>        &rhs2) const;
    unsigned int
    distribute_to_children (const typename DoFHandler<dim>::cell_iterator &cell,
			    FEValues<dim>  &fe_values,
			    const Vector<double>  &old_dof_values_u,
			    const Vector<double>  &old_dof_values_v,
			    Vector<double>        &right_hand_side1,
			    Vector<double>        &right_hand_side2);
};



#include <lac/full_matrix.h>



template <int dim>
class TimeStep_ErrorEstimation :  public virtual TimeStepBase_Wave<dim>
{
  public:
    TimeStep_ErrorEstimation ();
    virtual void estimate_error ();
    virtual void wake_up (const unsigned int wakeup_level);
    virtual void sleep (const unsigned int sleep_level);
    virtual void get_tria_refinement_criteria (Vector<float> &indicators) const;
    void get_error_indicators (Vector<float> &indicators) const;
    virtual string branch_signature () const = 0;
    
  protected:
    struct StatisticData 
    {
	StatisticData ();
	StatisticData (const double estimated_error);
	static void write_descriptions (ostream &out);
	void write (ostream &out) const;
	double estimated_error;
    };

    struct ErrorOnCell {
	double part[8];
	ErrorOnCell ();
	ErrorOnCell operator += (const ErrorOnCell &eoc);
	double sum () const;
    };


    struct CellwiseError
    {
	CellwiseError (const unsigned int n_errors);
	vector<ErrorOnCell>                    errors;
	typename vector<ErrorOnCell>::iterator next_free_slot;
    };

    Vector<float> estimated_error_per_cell;
    FullMatrix<double> embedding_matrix;
    FullMatrix<double> interpolation_matrix;
    FullMatrix<double> difference_matrix;
    StatisticData      statistic_data;
    void estimate_error_energy (const unsigned int which_variables);
    void estimate_error_dual ();
    void estimate_error_dual (const typename DoFHandler<dim>::cell_iterator &primal_cell,
			      const typename DoFHandler<dim>::cell_iterator &dual_cell,
			      const typename DoFHandler<dim>::cell_iterator &primal_cell_old,
			      const typename DoFHandler<dim>::cell_iterator &dual_cell_old,
			      CellwiseError                                 &cellwise_error,
			      FEValues<dim>                                 &fe_values) const;
    void compute_error_on_new_children (const typename DoFHandler<dim>::cell_iterator &primal_cell,
					const typename DoFHandler<dim>::cell_iterator &dual_cell,
					const Vector<double>  &local_u_old,
					const Vector<double>  &local_v_old,
					const Vector<double>  &local_u_bar_old,
					const Vector<double>  &local_v_bar_old,
					CellwiseError         &cellwise_error,
					FEValues<dim>         &fe_values) const;
    ErrorOnCell collect_error_from_children (const typename DoFHandler<dim>::cell_iterator &primal_cell_old,
					     const typename DoFHandler<dim>::cell_iterator &dual_cell_old,
					     const Vector<double>  &local_u,
					     const Vector<double>  &local_v,
					     const Vector<double>  &local_u_bar,
					     const Vector<double>  &local_v_bar,
					     const Vector<double>  &local_Ih_u_bar,
					     const Vector<double>  &local_Ih_v_bar,
					     const Vector<double>  &local_Ih_u_bar_old,
					     const Vector<double>  &local_Ih_v_bar_old,
					     FEValues<dim>  &fe_values) const;
    ErrorOnCell error_formula (const typename DoFHandler<dim>::active_cell_iterator &cell,
			       const Vector<double>  &local_u,
			       const Vector<double>  &local_v,
			       const Vector<double>  &local_u_bar,
			       const Vector<double>  &local_v_bar,
			       const Vector<double>  &local_u_old,
			       const Vector<double>  &local_v_old,
			       const Vector<double>  &local_u_bar_old,
			       const Vector<double>  &local_v_bar_old,
			       FEValues<dim>         &fe_values) const;
    ErrorOnCell error_formula (const typename DoFHandler<dim>::active_cell_iterator &cell,
			       const Vector<double>  &local_u,
			       const Vector<double>  &local_v,
			       const Vector<double>  &local_u_bar,
			       const Vector<double>  &local_v_bar,
			       const Vector<double>  &local_u_old,
			       const Vector<double>  &local_v_old,
			       const Vector<double>  &local_u_bar_old,
			       const Vector<double>  &local_v_bar_old,
			       const Vector<double>  &local_difference_u_bar,
			       const Vector<double>  &local_difference_v_bar,
			       const Vector<double>  &local_difference_u_bar_old,
			       const Vector<double>  &local_difference_v_bar_old,
			       FEValues<dim>         &fe_values) const;
    void make_interpolation_matrices ();
};


#include <vector>


template <int dim>
class TimeStep_Postprocess : public TimeStep_ErrorEstimation<dim> 
{
  public:
    virtual void postprocess_timestep ();
    virtual void wake_up (const unsigned int wakeup_level);
    virtual void sleep (const unsigned int sleep_level);
    virtual void end_sweep ();
    string branch_signature () const;

  protected:
    struct StatisticData 
    {
	static void write_descriptions (ostream &out,
					const WaveParameters<dim> &parameters);
	void write (ostream &out) const;
	vector<double> evaluation_results;
    };

    StatisticData statistic_data;

  private:
    void interpolate_dual_solution (Vector<double> &interpolated_u_bar,
				    Vector<double> &interpolated_v_bar) const;
};



template <int dim> class WaveParameters;



template <int dim>
class TimeStep :  public TimeStep_Primal<dim>, public TimeStep_Dual<dim>, public TimeStep_Postprocess<dim>
{
  public:
    TimeStep (const double               time,
	      const WaveParameters<dim> &parameters);

    virtual void wake_up (const unsigned int wakeup_level);
    virtual void sleep (const unsigned int sleep_level);
    virtual void end_sweep ();
    static void write_statistics_descriptions (ostream                   &out,
					       const WaveParameters<dim> &parameters);
    void write_statistics (ostream &out) const;
};

template <int dim> class TimeStep_Primal;
template <int dim> class TimeStep_Dual;



template <int dim>
class DualFunctional {
  public:
				     /**
				      * Constructor. Specify whether an
				      * actual functional needs the primal
				      * solution at all times or at the
				      * endtime. Default is #false# is
				      * both cases which means that the
				      * functional is linear.
				      */
    DualFunctional (const bool use_primal_problem            = false,
		    const bool use_primal_problem_at_endtime = false);

				     /**
				      * Return that part of the dual functional
				      * related to a delta function in time at
				      * the end time.
				      *
				      * The default is to return zero.
				      */
    virtual void compute_endtime_vectors (Vector<double> &final_u_bar,
					  Vector<double> &final_v_bar);

				     /**
				      * Return that part of the dual functional
				      * related to the regular time integral.
				      *
				      * The default is to return zero.
				      */
    virtual void compute_functionals (Vector<double> &j1,
				      Vector<double> &j2);

				     /**
				      * Return whether this object uses
				      * information from the primal problem
				      * (i.e. whether it is nonlinear or not).
				      * The necessary information is set in
				      * the constructor.
				      *
				      * This function refers to all times.
				      */
    bool use_primal_solutions () const;

    				     /**
				      * Return whether this object uses
				      * information from the primal problem
				      * (i.e. whether it is nonlinear or not).
				      * The necessary information is set in
				      * the constructor.
				      *
				      * This function refers to the solution
				      * at the end time. There are functionals
				      * which only evaluate at the endpoint
				      * but are nonlinear anyway. For them it
				      * is not necessary to reload the primal
				      * data at other times than the end time.
				      */
    bool use_primal_solutions_at_endtime () const;

				     /**
				      * Reset the functional to the present
				      * time level. This function needs to be
				      * called at each time level if the
				      * functional is nonlinear and at the
				      * endtime if the functional is nonlinear
				      * only at the endtime.
				      */
    virtual void reset (const TimeStep_Primal<dim> &primal_problem);

				     /**
				      * Reset the functional to the present
				      * time level. This function needs to be
				      * called at each time level. It resets
				      * pointers to the dof handler, the
				      * triangulation and several other
				      * objects which are needed to compute
				      * the dual functional.
				      */
    virtual void reset (const TimeStep_Dual<dim> &dual_problem);

				     /**
				      * Exception
				      */
    DeclException0 (ExcPrimalProblemNotRequested);
    
  protected:
    const bool use_primal_problem;
    const bool use_primal_problem_at_endtime;

    const Triangulation<dim> *tria;
    const Boundary<dim>      *boundary;
    const DoFHandler<dim>    *dof;
    const FiniteElement<dim> *fe;
    const Quadrature<dim>    *quadrature;
    const Quadrature<dim-1>  *quadrature_face;
    const Function<dim>      *density, *stiffness;

    const DoFHandler<dim>    *primal_dof;
    const FiniteElement<dim> *primal_fe;
    const Quadrature<dim>    *primal_quadrature;
    const Quadrature<dim-1>  *primal_quadrature_face;

    const Vector<double>     *u;
    const Vector<double>     *v;

    double       time;
    double       time_step;
    unsigned int step_no;
};





/**
 * Compute the dual functional which is approximately associated
 * with the end time energy in the high atmosphere above 4000km.
 * The energy in a domain $D$ is given by
 * $E_D = \int_D (v^2 + \nabla u a \nabla u)_{t=T}$ and the
 * associated functional for the error is approximately
 * $J(\Psi) = \int_D v_h(T) \psi + \nabla u_h(T) a \nabla \phi$.
 */
template <int dim>
class EndEnergy : public DualFunctional<dim> {
  public:
				     /**
				      * Constructor.
				      */
    EndEnergy (const bool use_primal_problem_at_any_time = false);

  protected:
    enum PartOfDomain { low_atmosphere, high_atmosphere };

				     /**
				      * Compute the initial values of the
				      * dual problem.
				      */
    void compute_vectors (const PartOfDomain pod,
			  Vector<double> &final_u_bar,
			  Vector<double> &final_v_bar) const;
};







/**
 * Let the point value of $u$ at the origin integrated over time
 * be the goal.
 */
template <int dim>
class IntegratedValueAtOrigin : public EndEnergy<dim> {
  public:
				     /**
				      * Evaluate the dual functionals and
				      * return the right hand side contributions
				      * thereof for the present time step.
				      */
    virtual void compute_functionals (Vector<double> &j1,
				      Vector<double> &j2);

				     /**
				      * Exception.
				      */
    DeclException0 (ExcVertexNotFound);
};





/**
 * Dual function corresponding to the #EvaluateSeismicSignal# class.
 */
template <int dim>
class SeismicSignal : public DualFunctional<dim> {
  public:
				     /**
				      * Evaluate the dual functionals and
				      * return the right hand side contributions
				      * thereof for the present time step.
				      */
    virtual void compute_functionals (Vector<double> &j1,
				      Vector<double> &j2);
};



/**
 * Compute the dual problem associated with the functional
 * $J(\Psi) = \int u ds$ with the integral being over some
 * parts of the boundary.
 */
template <int dim>
class EarthSurface : public DualFunctional<dim> {
  public:
				     /**
				      * Evaluate the dual functionals and
				      * return the right hand side contributions
				      * thereof for the present time step.
				      */
    virtual void compute_functionals (Vector<double> &j1,
				      Vector<double> &j2);
};



/**
 * Compute $J(\Psi) = \int_0^0.25 u(x=2,y,t=2.2) dy.
 */
template <int dim>
class SplitSignal : public DualFunctional<dim> {
  public:
				     /**
				      * Evaluate the dual functionals and
				      * return the right hand side contributions
				      * thereof for the present time step.
				      */
    virtual void compute_functionals (Vector<double> &j1,
				      Vector<double> &j2);
};





/**
 * 1d test case, evaluating the region (-.5,.5) at the endtime. Intended for some
 * tests on split triangulations with one fine and one coarse region.
 */
template <int dim>
class SplitLine : public DualFunctional<dim> {
  public:
				     /**
				      * Compute the initial values of the
				      * dual problem.
				      */
    void compute_endtime_vectors (Vector<double> &final_u_bar,
				  Vector<double> &final_v_bar);
};



/**
 * Compute $J(\Psi) = \int_{-0.6}^{-0.4} u(x,t=2.5) dx.
 */
template <int dim>
class OneBranch1d : public DualFunctional<dim> {
  public:
				     /**
				      * Evaluate the dual functionals and
				      * return the right hand side contributions
				      * thereof for the present time step.
				      */
    virtual void compute_functionals (Vector<double> &j1,
				      Vector<double> &j2);
};



/**
 * Compute $J(\Psi) = \int_{-0.1}^{0.1} u(x,t=2.4) dx.
 */
template <int dim>
class SecondCrossing : public DualFunctional<dim> {
  public:
				     /**
				      * Evaluate the dual functionals and
				      * return the right hand side contributions
				      * thereof for the present time step.
				      */
    virtual void compute_functionals (Vector<double> &j1,
				      Vector<double> &j2);
};



/**
 */
template <int dim>
class HuyghensWave : public DualFunctional<dim> {
  public:
				     /**
				      * Evaluate the dual functionals and
				      * return the right hand side contributions
				      * thereof for the present time step.
				      */
    virtual void compute_functionals (Vector<double> &j1,
				      Vector<double> &j2);
};


#include <grid/forward_declarations.h>
#include <base/forward_declarations.h>
#include <lac/forward_declarations.h>
#include <base/exceptions.h>
#include <string>





/**
 * This class provides a simple interface to do arbitrary evaluations of
 * the numerical solution. Concrete classes implementing evaluations
 * need access to the solution vectors #u# and #v# as well as to the
 * triangulation and the associated degrees of freedoms, which is what
 * this class provides. This way is chosen to separate the problem
 * classes which do the actual solution from the evaluation classes, since
 * they don't need to know much about the solution classes apart from
 * the solution itself. Thus, we reduce dependencies which speeds up
 * compilation and makes software engineering more simple.
 */
template <int dim>
class EvaluationBase {
  public:
				     /**
				      * Constructor. Set all pointers in this
				      * class to invalid values.
				      */
    EvaluationBase ();

				     /**
				      * Destructor. Does nothing but needs
				      * to be declared to make it virtual.
				      */
    virtual ~EvaluationBase () {};

				     /**
				      * Reset pointers to triangulation, dof
				      * handler, quadrature formulae etc.
				      * to the right values for the time level
				      * to be evaluated next. This function
				      * needs to be called each time an
				      * evaluation is to take place.
				      */
    virtual void reset_timelevel (const TimeStep_Primal<dim> &target);

				     /**
				      * Template for the evaluation functions.
				      * Return one value for the output file.
				      */
    virtual double evaluate () = 0;

    				     /**
				      * Reset the evaluator for the
				      * next sweep. This may be useful
				      * if you want to sum up the contributions
				      * of each time step and print them
				      * at the end; you then have to
				      * reset the sum at the start of
				      * the next sweep, which is done through
				      * this function.
				      *
				      * Default is: do nothing.
				      */
    virtual void reset ();

    				     /**
				      * Print the result at the end of
				      * each sweep. This function may
				      * print lines of data with four
				      * spaces at the beginning of each
				      * line.
				      *
				      * Default is: do nothing.
				      */
    virtual void print_final_result (ostream &out);

				     /**
				      * Return the final result as a number
				      * for the result file.
				      *
				      * Default is: do nothing.
				      */
    virtual double get_final_result ();

				     /**
				      * Return a brief string of description
				      * which will go into the first line
				      * of the "results" file.
				      */
    virtual string description () const = 0;

				     /**
				      * Exception.
				      */
    DeclException0 (ExcIO);
    
  protected:
				     /**
				      * Pointers to the solution vectors
				      * of the primal problem.
				      */
    const Vector<double>     *u, *v;

				     /**
				      * Underlying triangulation.
				      */
    const Triangulation<dim> *tria;

				     /**
				      * Boundary object.
				      */
    const Boundary<dim>      *boundary;

				     /**
				      * Degrees of freedom of the primal
				      * problem.
				      */
    const DoFHandler<dim>    *dof;

				     /**
				      * Primal finite element.
				      */
    const FiniteElement<dim> *fe;

				     /**
				      * Quadrature rule appropriate for
				      * the primal finite element.
				      */
    const Quadrature<dim>    *quadrature;

				     /**
				      * Same for quadrature on faces.
				      */
    const Quadrature<dim-1>  *quadrature_face;

				     /**
				      * Density and stiffness coefficients
				      * for the modell presently under
				      * investigation.
				      */
    const Function<dim>      *density, *stiffness;

				     /**
				      * Continuous time of the time step
				      * we are evaluating at present.
				      */
    double                    time;

				     /**
				      * Length of the last time step, i.e. in
				      * the backward direction in time. If
				      * this is the first timestep, the this
				      * value is set to zero.
				      */
    double                    time_step;

				     /**
				      * Number of that time step.
				      */
    unsigned int              step_no;

				     /**
				      * Base of the filenames under which
				      * we shall store our results.
				      */
    string                    base_file_name;
};





/**
 * This class is a common base class to the following two. It provides
 * for some infrastructure for evaluations computing the energy in part
 * of the domain and computing the in/outflow of energy.
 *
 * Central is the #compute_energy# function, which takes an argument
 * describing which part of the domain to take and returns the energy
 * therein.
 */
template <int dim>
class EvaluateEnergyContent : public EvaluationBase<dim> {
  public:
				     /**
				      * Constructor.
				      */
    EvaluateEnergyContent ();

				     /**
				      * Reset the accumulated energy to zero.
				      */
    virtual void reset ();

  protected:
				     /**
				      * Enum denoting for which of the two
				      * subdomains the computation is to be
				      * performed.
				      */
    enum PartOfDomain { low_atmosphere, high_atmosphere };

				     /**
				      * Compute the energy for the given
				      * subdomain.
				      */
    double compute_energy (const PartOfDomain pod) const;
      
  protected:
				     /**
				      * Energy in the domain in the previous
				      * time step. This information is needed
				      * to compute the accumulated in/outflux
				      * of energy from the domain.
				      */
    double old_energy;

				     /**
				      * Accumulated in/outflux into/from the
				      * domain integrated over time.
				      */
    double integrated_outflux;
};







/**
 * Evaluate the value of $u$ at the origin, i.e. $u(0,0)$.
 *
 * As final result, the time integrated value at the origin is computed.
 * The origin shall be a vertex in the finest grid.
 */
template <int dim>
class EvaluateIntegratedValueAtOrigin : public EvaluationBase<dim> {
  public:
    EvaluateIntegratedValueAtOrigin ():
		    integrated_value (0) {};
    
    virtual double evaluate ();
    virtual void print_final_result (ostream &out);
    virtual double get_final_result ();
    virtual string description () const;

				     /**
				      * Reset the average value to zero.
				      */
    virtual void reset ();

				     /**
				      * Exception.
				      */
    DeclException0 (ExcVertexNotFound);
    
  private:
    double integrated_value;
};






/**
 * Integrate the value of $u_h$ at the top boundary over $x$ and $t$ using a
 * highly oscillatory weight.
 */
template <int dim>
class EvaluateSeismicSignal : public EvaluationBase<dim> {
  public:
    EvaluateSeismicSignal () :
		    result (0) {};

    static inline double weight (const Point<dim> &p, const double time) {
      const double pi = 3.14159265359;
      return sin(3*pi*p(0))*sin(5*pi*time/2);
    };
    
    
    virtual double evaluate ();
    virtual void print_final_result (ostream &out);
    virtual double get_final_result ();
    virtual string description () const;

				     /**
				      * Reset the value to zero.
				      */
    virtual void reset ();

  private:
    double result;
};



/**
 * Integrate the value of $u_h$ at the top line $x=1.5, y=0..1/16$ at $t=1.6..1.8$.
 */
template <int dim>
class EvaluateSplitSignal : public EvaluationBase<dim> {
  public:
    EvaluateSplitSignal () :
		    result (0) {};
    
    
    virtual double evaluate ();
    virtual void print_final_result (ostream &out);
    virtual double get_final_result ();
    virtual string description () const;

				     /**
				      * Reset the value to zero.
				      */
    virtual void reset ();

  private:
    double result;
};

    

template <int dim>
class EvaluateOneBranch1d : public EvaluationBase<dim> {
  public:
    EvaluateOneBranch1d () :
		    result (0) {};
    
    
    virtual double evaluate ();
    virtual void print_final_result (ostream &out);
    virtual double get_final_result ();
    virtual string description () const;

				     /**
				      * Reset the value to zero.
				      */
    virtual void reset ();

  private:
    double result;
};




template <int dim>
class EvaluateSecondCrossing1d : public EvaluationBase<dim> {
  public:
    EvaluateSecondCrossing1d () :
		    result (0) {};
    
    
    virtual double evaluate ();
    virtual void print_final_result (ostream &out);
    virtual double get_final_result ();
    virtual string description () const;

				     /**
				      * Reset the value to zero.
				      */
    virtual void reset ();

  private:
    double result;
};



template <int dim>
class EvaluateHuyghensWave : public EvaluationBase<dim> {
  public:
    EvaluateHuyghensWave () :
		    integrated_value (0),
		    weighted_value (0) {};
    
    
    virtual double evaluate ();
    virtual void print_final_result (ostream &out);
    virtual double get_final_result ();
    virtual string description () const;

				     /**
				      * Reset the value to zero.
				      */
    virtual void reset ();

  private:
    double integrated_value, weighted_value;
};


    

template <int dim> class DataOutStack;



/**
 * This class has some data members which are shared between the different
 * time steps within one sweep. Unlike the #SweepInfo# class, the members
 * do not collect information for later output, but provide services to
 * the time steps.
 */
template <int dim>
class SweepData 
{
  public:
    SweepData (const bool use_data_out_stack);
    ~SweepData ();

    DataOutStack<dim> *data_out_stack;
};



#include <base/timer.h>
#include <iostream>
#include <list>



/**
 * This class provides some data members which collect information on the
 * different time steps of one sweep.
 */
class SweepInfo 
{
  public:
    struct Data 
    {
					 /**
					  * Constructor. Set all fields to
					  * their initial values.
					  */
	Data ();
	
	double accumulated_error;

	unsigned int cells;
	unsigned int primal_dofs;
	unsigned int dual_dofs;
    };


    struct Timers 
    {
	Timer grid_generation;
	Timer primal_problem;
	Timer dual_problem;
	Timer error_estimation;
	Timer postprocessing;
    };


    Data & get_data ();

    Timers & get_timers ();
    
    
    template <int dim>
    void write_summary (const list<EvaluationBase<dim>*> & eval_list,
			ostream &out) const;
    
  private:
    Data data;
    Timers timers;
};



#include <lac/sparse_matrix.h>




/**
 * Enum denoting the different possibilities to precondition a solver.
 */
enum Preconditioning {
      no_preconditioning,
      jacobi,
      sor,
      ssor
};



/**
 * Wrapper for the #SparseMatrix<double># class which handles the preconditioning.
 */
class UserMatrix :  public SparseMatrix<double> {
  public:
				     /**
				      * Constructor. The parameter specifies
				      * which way to precondition.
				      */
    UserMatrix (Preconditioning p) :
		    SparseMatrix<double> (),
		    preconditioning (p) {};

				     /**
				      * Constructor. The second parameter
				      * specifies which way to precondition.
				      * The first parameter is simply passed
				      * down to the base class.
				      */
    UserMatrix (const SparsityPattern &sparsity,
		Preconditioning       p) :
		    SparseMatrix<double>(sparsity),
		    preconditioning (p) {};

				     /**
				      * Precondition a vector #src# and write
				      * the result into #dst#. This function
				      * does not much more than delegating to
				      * the respective #precondition_*#
				      * function of the base class, according
				      * to the preconditioning method specified
				      * to the constructor of this class.
				      */
    void precondition (Vector<double> &dst, const Vector<double> &src) const;

  private:
				     /**
				      * Variable denoting the preconditioning
				      * method.
				      */
    Preconditioning preconditioning;
};


#include <base/exceptions.h>
#include <grid/forward_declarations.h>

#include <string>


string int_to_string (const unsigned int i, const unsigned int digits);


template <typename number>
inline number sqr (const number a) {
  return a*a;
};



/**
 * This is a helper class which has a collection of static elements and returns
 * the right finite element as a pointer when the name of the element is given.
 * It is also able to return the correct quadrature formula for domain and
 * boundary integrals for the specified finite element.
 */
template <int dim>
struct FEHelper {
    static const FEQ1<dim>           fe_linear;
    static const FEQ2<dim>           fe_quadratic_sub;
#if 2 < 3    
    static const FEQ3<dim>           fe_cubic_sub;
    static const FEQ4<dim>           fe_quartic_sub;
#endif

    static const QGauss2<dim>        q_gauss_2;
    static const QGauss3<dim>        q_gauss_3;
    static const QGauss4<dim>        q_gauss_4;
    static const QGauss5<dim>        q_gauss_5;
    static const QGauss6<dim>        q_gauss_6;
    static const QGauss7<dim>        q_gauss_7;
    static const QGauss8<dim>        q_gauss_8;

    static const QGauss2<dim-1>      q_gauss_2_face;
    static const QGauss3<dim-1>      q_gauss_3_face;
    static const QGauss4<dim-1>      q_gauss_4_face;
    static const QGauss5<dim-1>      q_gauss_5_face;
    static const QGauss6<dim-1>      q_gauss_6_face;
    static const QGauss7<dim-1>      q_gauss_7_face;
    static const QGauss8<dim-1>      q_gauss_8_face;

				     /**
				      * Return a reference to the finite
				      * element specified by the name
				      * #name#.
				      */
    static const FiniteElement<dim> & get_fe (const string &name);

				     /**
				      * Return the correct domain quadrature
				      * formula for the finite element denoted
				      * by the name #name#.
				      */
    static const Quadrature<dim>    & get_quadrature (const string &name);

				     /**
				      * Return the correct boundary quadrature
				      * formula for the finite element denoted
				      * by the name #name#.
				      */
    static const Quadrature<dim-1>  & get_quadrature_face (const string &name);
};



#include <base/forward_declarations.h>
#include <grid/forward_declarations.h>
#include <list>
#include <string>

template <int dim> class DualFunctional;
template <int dim> class EvaluationBase;




/**
 * This is a class holding all the input parameters to the program. It is more
 * or less a loose collection of data and the only purpose of this class is
 * to assemble all the parameters and the functions evaluating them from the
 * input file at one place without the need to scatter this functionality
 * all over the program.
 *
 *
 * \section{Description of the input parameters}
 *
 * Note that this list may not be up-tp-date at present.
 *
 * \subsection{Subsection #Grid#}
 * \begin{itemize}
 * \item #Coarse mesh#: Names a grid to be taken as a coarse grid. The following
 *    names are allowed:
 *    \begin{itemize}
 *    \item #uniform channel#: The domain is $[0,3]\times[0,1]$, triangulated
 *        by three cells. Left and right boundary are of Dirichlet type, top
 *        and bottom boundary are of homogeneous Neumann type.
 *    \item #split channel bottom#: As above, but the lower half is refined once
 *        more than the top half.
 *    \item #split channel {left | right}#: Same as #uniform channel#, but with
 *        cells on the left or right, according to the last word, more refined
 *        than on the other side.
 *    \item #square#: $[-1,1]\times[-1,1]$.
 *    \item #seismic square#: same as #square#, but with Neumann boundary
 *        at top.
 *    \item #temperature-square#: Square with size $400,000,000$ (we use the
 *        cgs system, so this amounts to 4000 km).
 *    \item #temperature-testcase#: As above, but with a sequence of
 *        continuously growing cells set atop to avoid the implementation of
 *        absorbing boundary conditions. The left boundary is of Neumann
 *        type (mirror boundary).
 *    \item #random#: Unit square, but randomly refined to test for correctness
 *        of the time stepping scheme.
 *    \item #earth#: Circle with radius 6371 (measured in km).
 *    \end{itemize}
 * \item #Initial refinement#: States how often the grid named by the above
 *    parameter shall be globally refined to form the coarse mesh.
 * \item #Maximum refinement#: maximum refinement level a cell may attain.
 *    Cells with such a refinement level are flagged as others are, but they
 *    are not refined any more; it is therefore not necessary to lower the
 *    fraction of cells to be refined in order to avoid the refinement of a
 *    similar number of cells with a lower level number.
 *
 *    The default to this value is zero, meaning no limit.
 * \item #Refinement fraction#: Upon refinement, those cells are refined which
 *    together make up for a given fraction of the total error. This parameter
 *    gives that fraction. Default is #0.95#.
 * \item #Coarsening fraction#: Similar as above, gives the fraction of the
 *    total error for which the cells shall be coarsened. Default is #0.03#.
 * \item #Top cell number deviation#: Denotes a fraction by which the number of
 *    cells on a time level may be higher than the number of cells on the
 *    previous time level. This and the next two parameters help to avoid
 *    to much differing grids on the time levels and try to smooth the numbers
 *    of cells as a function of time. The default value is #0.1#.
 * \item #Bottom cell number deviation#: Denotes the fraction by which the
 *    number of cells on a time level may be lower than on the previous time
 *    level. Default is #0.03#.
 * \item #Cell number correction steps#: Usually, the goal denoted by the two
 *    parameters above cannot be reached directly because the number of cells
 *    is modified by grid regularization etc. The goal can therefore only be
 *    reached by an iterative process. This parameter tells how many iterations
 *    of this process shall be done. Default is #2#.
 * \end{itemize}
 *
 * \subsection{Subsection #Equation data#}
 * \begin{itemize}
 * \item #Coefficient#: Names for the different coefficients for the Laplace
 *    like part of the wave operator. Allowed values are:
 *    \begin{itemize}
 *    \item #unit#: Constant one.
 *    \item #kink#: One for $y<\frac 13$, 4 otherwise.
 *    \item #gradient#: $1+8*y^2$.
 *    \item #tube#: $0.2$ for $|x|<0.2$, one otherwise.
 *    \item #temperature VAL81#: Coefficient computed from the temperature
 *        field given by Varnazza, Avrett, Loeser 1981.
 *    \item #temperature kolmogorov#: Broadened temperature spectrum.
 *    \item #temperature undisturbed#: Quiet atmosphere.
 *    \item #temperature monochromatic 20s#: Temperature as computed with
 *        shock waves with $T=20s$.
 *    \item #temperature monochromatic 40s#: Temperature as computed with
 *        shock waves with $T=40s$.
 *    \end{itemize}
 * \item #Initial u#: Names for the initial value for the amplitude. Allowed
 *    names are:
 *    \begin{itemize}
 *    \item #zero#: $u_0=0$.
 *    \item #eigenmode#: $u_0=sin(2\pi x)sin(2\pi y)$.
 *    \item #bump#: $u_0=(1-\frac{\vec x^2}{a^2})e^{-\frac{\vec x^2}{a^2}}$
 *        for $|\vec x|<a$ and $u_0=0$ otherwise. $a=0.1$
 *    \item #center-kink#: $u_0=r/a$ for $r<a$, $u_0=2-r/a$ for $a<r<2a$,
 *        $u=0$ otherwise. $a=0.1$, $r=|\vec x|$.
 *    \item #shifted bump#: Same as #bump# but the center of the bump is
 *        located at $x=0.5, y=0$.
 *    \item #tube#: $u_0=1$ for $|x|<0.2, zero otherwise.
 *    \end{itemize}
 * \item #Initial v#: Names for the initial value for the amplitude. Allowed
 *    names are the same as above.
 * \item #Boundary#: Names for the boundary functions. The boundary values
 *    for $u$ and $v$ are always set together. The boundary values apply only
 *    to those boundary parts which are of Dirichlet type. Allowed names are:
 *    \begin{itemize}
 *    \item #zero#: Homogeneous boundary values.
 *    \item #wave from left#: For $t<T=0.4$ we set $u=sin^2(\pi \frac tT)$ at
 *        the boundary where $x=0$.
 *    \item #wave from left center#: For $t<T=0.4$ and $0.4<y<0.6$ we set
 *        $u=sin^2(\pi \frac tT) (y-0.4) (0.6-y)$ at
 *        the boundary where $x=0$.
 *    \item #wave from left bottom#: For $t<T=60s$ and $r=|\vec x|<a=5000000cm=50km$
 *        let $u=(cos(\pi/2 r/a) sin(\pi t/T))^2$.
 *        This boundary condition is only suited to the temperature domains. 
 *    \end{itemize}
 * \end{itemize}
 *
 * \subsection{Subsection #Time stepping#}
 * \begin{itemize}
 * \item #Primal method#: Time stepping method for the primal problem.
 *     Allowed values are:
 *     \begin{itemize}
 *     \item #theta#: Use the $\theta$ scheme with the $\theta$-parameter
 *         as given below.
 *     \item #fractional step#: Use the fractional step $\theta$ scheme.
 *     \end{itemize}
 * \item #Dual method#: Time stepping method for the dual problem. Allowed
 *     values are the same as above. Note that the fractional step scheme
 *     is not implemented for right hand sides not equal to zero, i.e. the
 *     fractional step scheme will fail of the error functional evaluates
 *     to non-zero at times not equal to the end time.
 * \item #Theta#: $\theta$ parameter for the $\theta$ time stepping scheme.
 *     $\theta=1/2$ denotes the Crank-Nicolson scheme.
 * \item #Time step#: Selfdocumenting.
 * \item #End time#: Selfdocumenting.
 * \end{itemize}
 */
template <int dim>
class WaveParameters
{
  public:
				     /**
				      * Constructor.
				      */
    WaveParameters ();

				     /**
				      * Destructor.
				      */
    ~WaveParameters ();

				     /**
				      * Declare all the parameters to the
				      * given parameter handler.
				      */
    void declare_parameters (ParameterHandler &prm);

				     /**
				      * Extract the parameters values provided
				      * by the input file and/or the default
				      * values from the parameter handler.
				      */
    void parse_parameters (ParameterHandler &prm);

				     /**
				      * Delete the contents of this class and
				      * set up a clean state.
				      */
    void delete_parameters ();

				     /**
				      * Enum holding a list of possible coarse
				      * mesh choices.
				      */
    enum InitialMesh {
	  uniform_channel,
	  split_channel_bottom,
	  split_channel_right,
	  split_channel_left,
	  square,
	  ring,
	  seismic_square,
	  earth,
	  line,
	  split_line
    };

				     /**
				      * Enum holding a list of possible
				      * boundary condition choices.
				      */
    enum BoundaryConditions {
	  wave_from_left,
	  fast_wave_from_left,
	  wave_from_left_center,
	  wave_from_left_bottom,
	  zero
    };

				     /**
				      * Enum denoting possible strategies
				      * for output of meshes and solutions.
				      * This enum tells us, at which sweeps
				      * data is to be written.
				      */
    enum WriteStrategy {
	  never,
	  all_sweeps,
	  last_sweep_only
    };

				     /**
				      * Boundary values. Continuous function
				      * of space and time.
				      */
    Function<dim>      *boundary_values_u;

				     /**
				      * Same for the velocity variable v.
				      */
    Function<dim>      *boundary_values_v;

				     /**
				      * Initial values for u.
				      */
    Function<dim>      *initial_u;

				     /**
				      * Same for the velocity variable v.
				      */
    Function<dim>      *initial_v;

				     /**
				      * Object describing the boundary. By
				      * default the domain is polygonal made
				      * from the vertices of the coarsest
				      * triangulation. However, some of the
				      * example geometries set in
				      * #make_coarse_grid# may set this variable
				      * to another address. The object pointed
				      * will be deleted at the end of the
				      * lifetime of this object; when setting
				      * this variable to another object, you
				      * may want to delete the object pointed
				      * to previously.
				      */
    const Boundary<dim>*boundary;

				     /**
				      * Function denoting the coefficient
				      * within the generalized laplacian
				      * operator.
				      */
    Function<dim>      *density;

				     /**
				      * Same for the stiffness parameter.
				      */
    Function<dim>      *stiffness;

				     /**
				      * Store whether the density is a function
				      * that is constant in space (not
				      * necessarily in time as well, but at
				      * each fixed time).
				      */
    bool density_constant;

				     /**
				      * Same thing for the stiffness parameter.
				      */
    bool stiffness_constant;
    
				     /**
				      * Pointer to an object denoting the
				      * error functional.
				      */
    DualFunctional<dim>*dual_functional;
    
				     /**
				      * Level of initial refinement, i.e. the
				      * minimum level cells on all grids at
				      * all times need to have.
				      */
    unsigned int        initial_refinement;

				     /**
				      * Maximum refinement level a cell may
				      * have. This one defaults to zero,
				      * meaning no limit.
				      */
    unsigned int        maximum_refinement;

				     /**
				      * Define structure of initial mesh:
				      * created by regular refinement of
				      * the coarsest mesh (uniform) or
				      * refine one half once more than
				      * the other (split) or some other
				      */
    Triangulation<dim>      *coarse_grid;
    
				     /**
				      * Pair of numbers denoting the fraction
				      * of the total error for which the cells
				      * are to be refined (first) and
				      * coarsened (second).
				      */
    pair<double,double>      refinement_fraction;

				     /**
				      * Fraction by which the number of cells
				      * on a time level may differ from the
				      * number on the previous time level
				      * (first: top deviation, second: bottom
				      * deviation).
				      */
    pair<double,double>      cell_number_corridor;

				     /**
				      * Number of iterations to be performed
				      * to adjust the number of cells on a
				      * time level to those on the previous
				      * one.
				      */
    unsigned int             cell_number_correction_steps;

				     /**
				      * Shall we renumber the degrees of
				      * freedom according to the Cuthill-McKee
				      * algorithm or not.
				      */
    bool                     renumber_dofs;

				     /**
				      * Compare error indicators globally or
				      * refine each time step separately from
				      * the others.
				      */
    bool                     compare_indicators_globally;
    
    				     /**
				      * Parameters for the time discretization
				      * of the two equations using the
				      * theta scheme.
				      */
    double              theta;

				     /**
				      * Time step size.
				      */
    double              time_step;
    
				     /**
				      * Time up to which we want to compute.
				      */
    double              end_time;

				     /**
				      * Mode of preconditioning.
				      */
    Preconditioning     preconditioning;

				     /**
				      * Use extrapolated values of the old
				      * solutions as starting values for
				      * the solver on the new timestep.
				      */
    bool                extrapolate_old_solutions;
    
				     /**
				      * Directory to which we want the output
				      * written.
				      */
    string              output_directory;

				     /**
				      * Directory to which we want the temporary
				      * file to be written.
				      */
    string              tmp_directory;
    
				     /**
				      * Format in which the results on the
				      * meshes is to be written to files.
				      */
    string              output_format;

				     /**
				      * Denotes in which sweeps the solution is
				      * to be written.
				      */
    WriteStrategy       write_solution_strategy;

				     /**
				      * Denote the interval between the steps
				      * which are to be written.
				      */
    unsigned int        write_steps_interval;

				     /**
				      * Specify whether error information is
				      * to be written as cell data or node
				      * data.
				      */
    bool                write_error_as_cell_data;

				     /**
				      * Flag determining whether we shall
				      * write out the data of the different
				      * time steps stacked together for a
				      * whole sweep, and into one file for
				      * the whole sweep.
				      */
    bool                write_stacked_data;

				     /**
				      * Same as #write_steps_interval#, but
				      * for stacked output.
				      */
    unsigned int        write_stacked_interval;
    
				     /**
				      * Write statistics for the error
				      * distribution in each sweep.
				      */
    bool                produce_error_statistics;

				     /**
				      * Number of histogram intervals for
				      * the error statistics.
				      */
    unsigned int        error_statistic_intervals;

				     /**
				      * How to break the intervals: linear
				      * or logarithmic.
				      */
    string              error_statistics_scaling;
    
				     /**
				      * Names of the finite element classes to
				      * be used for the primal and dual problems.
				      */
    string              primal_fe, dual_fe;
    
				     /**
				      * Strategy for mesh refinement.
				      */
    enum { energy_estimator, dual_estimator } refinement_strategy;

				     /**
				      * Try to adjust the mesh to the error
				      * functional as well as to the dual
				      * solution. For the dual solution, an
				      * energy estimator is used.
				      */
    bool adapt_mesh_to_dual_solution;

				     /**
				      * When adapting the mesh for the dual
				      * problem as well, we have to weigh
				      * the error indicator for the dual
				      * problem with that for the primal
				      * one. This is the factor.
				      */
    double primal_to_dual_weight;

				     /**
				      * Number of sweeps at the beginning
				      * where the energy estimator is to
				      * be used rather than the dual
				      * estimator.
				      */
    unsigned int initial_energy_estimator_sweeps;

				     /**
				      * How many adaptive cycles of solving
				      * the whole problem shall be made.
				      */
    unsigned int        number_of_sweeps;

				     /**
				      * List of operations which shall be
				      * done on each time step after finishing
				      * a sweep.
				      */
    list<EvaluationBase<dim>*> eval_list;

				     /**
				      * Symbolic name of the boundary conditions
				      * (additionally to the boundary functions
				      * themselves), which may be used by some
				      * of the evaluations and other functionals
				      * in the program.
				      */
    BoundaryConditions boundary_conditions;

				     /**
				      * Exception.
				      */
    DeclException1 (ExcParameterNotInList,
		    string,
		    << "The given parameter <" << arg1 << "> is not "
		    << "recognized to be a valid one.");
    
  private:

				     /**
				      * Undefined copy constructor.
				      */
    WaveParameters (const WaveParameters &);

				     /**
				      * Undefined copy operator.
				      */
    WaveParameters & operator = (const WaveParameters &);
    
    
				     /**
				      * List of names for the initial values.
				      * Make this a member of the templated
				      * class since the supported initial
				      * values could be different from
				      * dimension to dimension.
				      */
    static const string initial_value_names;

				     /**
				      * Names of coefficient functions. The
				      * same applies as for
				      * #initial_value_names#.
				      */
    static const string coefficient_names;

    				     /**
				      * Names of boundary value functions. The
				      * same applies as for
				      * #initial_value_names#.
				      */
    static const string boundary_function_names;

    				     /**
				      * Names of error functionals. The
				      * same applies as for
				      * #initial_value_names#.
				      */
    static const string dual_functional_names;

    
				     /**
				      * Set the initial function pointers
				      * depending on the given names.
				      */
    void set_initial_functions (const string &u_name,
				const string &v_name);

				     /**
				      * Set the coefficient functions.
				      */
    void set_coefficient_functions (const string &name);

				     /**
				      * Set the boundary values.
				      */
    void set_boundary_functions (const string &name);

				     /**
				      * Make a list of evaluations to be
				      * performed after each sweep.
				      */
    void make_eval_list (const string &names);

				     /**
				      * Set the dual functional after
				      * which the dual solution will be
				      * computed.
				      */
    void set_dual_functional (const string &name);

				     /**
				      * Create the coarse grid for
				      * this run.
				      */
    void make_coarse_grid (const string &name);
};



#include <numerics/time_dependent.h>

template <int dim> class WaveParameters;
template <int dim> class DataOutStack;
class SweepInfo;



/**
 * Top-level class of the timestepping mechanism. This class manages
 * the execution and solution of primal and dual problem, of computing
 * error estimates and doing the refinement of grids.
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int dim>
class TimestepManager : public TimeDependent {
  public:
				     /**
				      * Constructor.
				      */
    TimestepManager (const WaveParameters<dim> &parameters);

				     /**
				      * Run a complete sweep, consisting
				      * of the solution of the primal problem,
				      * the solution of the dual problem if
				      * requested, computation of error
				      * quantities and refinement.
				      */
    void run_sweep (const unsigned int sweep_no);

				     /**
				      * Exception
				      */
    DeclException0 (ExcIO);
    
  private:
				     /**
				      * Reference to the global parameters
				      * object.
				      */
    const WaveParameters<dim> &parameters;

				     /**
				      * Refine the grids, or, better, find
				      * out which cells need to be refined.
				      * Refinement is done by a following
				      * sweep.
				      */
    void refine_grids ();

				     /**
				      * Write some statistics to a file.
				      */
    void write_statistics (const SweepInfo &sweep_info) const;

				     /**
				      * Write the data stacked together
				      * from all the time steps into
				      * one single file.
				      */
    void write_stacked_data (DataOutStack<dim> &data_out_stack) const;
};





/**
 * Top-level class providing the set up of a simulation. The
 * class provides an interface suitable to the #MultipleParameterLoop#
 * class to do several simulations in a row, stores global simulation
 * parameters, and so on.
 *
 * @author Wolfgang Bangerth, 1998, 1999
 */
template <int dim>
class WaveProblem :  public MultipleParameterLoop::UserClass {
  public:

				     /**
				      * Constructor.
				      */
    WaveProblem ();

				     /**
				      * Destructor.
				      */
    virtual ~WaveProblem ();

				     /**
				      * Put this object into a clean state.
				      * This function is called at the
				      * beginning of each loop by the
				      * #MultipleParameterHandler#.
				      */
    virtual void create_new (const unsigned int run_no);

				     /**
				      * Make the list of parameters known
				      * to the parameter handler. This
				      * function only delegates its work
				      * to the #parameters# sub-object.
				      */
    virtual void declare_parameters (ParameterHandler &prm);

    				     /**
				      * Parse the list of parameters given
				      * by the parameter handler. This
				      * function only delegates its work
				      * to the #parameters# sub-object.
				      */
    virtual void parse_parameters (ParameterHandler &prm);

				     /**
				      * Run a complete simulation.
				      */
    virtual void run (ParameterHandler &prm);

  private:
				     /**
				      * Object holding the parameters of
				      * the present simulation.
				      */
    WaveParameters<dim> parameters;
};





/* $Id$ */


#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_constraints.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <base/quadrature.h>
#include <base/function.h>

#include <cmath>



/*------------------------ DualFunctional --------------------------------*/

template <int dim>
DualFunctional<dim>::DualFunctional (const bool use_primal_problem,
				     const bool use_primal_problem_at_endtime) :
		use_primal_problem (use_primal_problem),
		use_primal_problem_at_endtime (use_primal_problem_at_endtime),
		tria (0),
		boundary (0),
		dof (0),
		fe(0),
		quadrature(0),
		quadrature_face(0),
		density(0),
		stiffness(0),
		primal_dof(0),
		primal_fe(0),
		primal_quadrature(0),
		primal_quadrature_face(0),
		u(0),
		v(0),
		time(0),
		time_step(0),
		step_no(0)
{};



template <int dim>
void DualFunctional<dim>::compute_functionals (Vector<double> &j1,
					       Vector<double> &j2) {
  j1.reinit (dof->n_dofs());
  j2.reinit (dof->n_dofs());
};



template <int dim>
void DualFunctional<dim>::compute_endtime_vectors (Vector<double> &final_u_bar,
						   Vector<double> &final_v_bar) {
  final_u_bar.reinit (dof->n_dofs());
  final_v_bar.reinit (dof->n_dofs());
};



template <int dim>
bool DualFunctional<dim>::use_primal_solutions () const {
  return use_primal_problem;
};



template <int dim>
bool DualFunctional<dim>::use_primal_solutions_at_endtime () const {
  return use_primal_problem_at_endtime;
};



template <int dim>
void DualFunctional<dim>::reset (const TimeStep_Primal<dim> &primal_problem) {
  Assert (use_primal_problem ||
	  (use_primal_problem_at_endtime &&
	   (primal_problem.parameters.end_time==primal_problem.time)),
	  ExcPrimalProblemNotRequested());

  primal_dof             = primal_problem.dof_handler;
  primal_fe              = &primal_problem.fe;
  primal_quadrature      = &primal_problem.quadrature;
  primal_quadrature_face = &primal_problem.quadrature_face;

  u = &primal_problem.u;
  v = &primal_problem.v;
};



template <int dim>
void DualFunctional<dim>::reset (const TimeStep_Dual<dim> &dual_problem) {
  tria            = dual_problem.tria;
  boundary        = dual_problem.parameters.boundary;
  dof             = dual_problem.dof_handler;
  fe              = &dual_problem.fe;
  quadrature      = &dual_problem.quadrature;
  quadrature_face = &dual_problem.quadrature_face;
  density         = dual_problem.parameters.density;
  stiffness       = dual_problem.parameters.stiffness;
  time            = dual_problem.time;
  time_step       = (dual_problem.next_timestep == 0 ?
		     0 :
		     dual_problem.get_forward_timestep());
  step_no         = dual_problem.timestep_no;
};









/* ----------------------- EndEnergy ------------------------------*/


template <int dim>
EndEnergy<dim>::EndEnergy (const bool use_primal_problem) :
		DualFunctional<dim> (use_primal_problem, true) {};




template <int dim>
void EndEnergy<dim>::compute_vectors (const PartOfDomain pod,
				      Vector<double> &final_u_bar,
				      Vector<double> &final_v_bar) const {
  const double y_offset = 300000000;
  const double n_q_points = quadrature->n_quadrature_points;
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  
  final_u_bar.reinit (dof->n_dofs());
  final_v_bar.reinit (dof->n_dofs());

  DoFHandler<dim>::active_cell_iterator cell, primal_cell, endc;
  cell = dof->begin_active ();
  endc = dof->end ();
  primal_cell = primal_dof->begin_active();

  FEValues<dim> fe_values (*fe, *quadrature,
			   UpdateFlags(update_values         |
				       update_gradients      |
				       update_JxW_values     |
				       update_q_points));
  FEValues<dim> fe_values_primal (*primal_fe, *quadrature,
				  UpdateFlags(update_values | update_gradients));
  
  FullMatrix<double>  cell_matrix (dofs_per_cell, dofs_per_cell);

  vector<Tensor<1,dim> > local_u_grad (n_q_points);
  vector<double>         local_v (n_q_points);
  
  vector<double> density_values(quadrature->n_quadrature_points);
  vector<double> stiffness_values(quadrature->n_quadrature_points);

  vector<int> cell_dof_indices (dofs_per_cell);

  for (; cell!=endc; ++cell, ++primal_cell)
    {
				       // only consider cells in the specified
				       // domain
      switch (pod)
	{
	  case low_atmosphere:
		if (cell->center()(1) >= y_offset)
		  continue;
		break;
	  case high_atmosphere:
		if (cell->center()(1) < y_offset)
		  continue;
		break;
	};
      
      
      fe_values.reinit (cell);
      fe_values_primal.reinit (primal_cell);
      fe_values_primal.get_function_values (*v, local_v);
      fe_values_primal.get_function_grads (*u, local_u_grad);

				       // get the coefficients at the
				       // quadrature points
      density->value_list (fe_values.get_quadrature_points(),
			   density_values);
      stiffness->value_list (fe_values.get_quadrature_points(),
			     stiffness_values);
      
				       // set up a vector of the gradients
				       // of the finite element basis
				       // functions on this face at the
				       // quadrature points
      const vector<vector<Tensor<1,dim> > > &shape_grads = fe_values.get_shape_grads ();
      const FullMatrix<double>            &shape_values = fe_values.get_shape_values ();
      const vector<double>      &JxW_values (fe_values.get_JxW_values());
      
      vector<double> local_functional1 (dofs_per_cell, 0);
      vector<double> local_functional2 (dofs_per_cell, 0);
      for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
	for (unsigned int point=0; point<n_q_points; ++point) 
	  {
	    local_functional1[shape_func] += local_u_grad[point] *
					     shape_grads[shape_func][point] *
					     stiffness_values[point] *
					     JxW_values[point];
	    local_functional2[shape_func] += local_v[point] *
					     shape_values(shape_func,point) *
					     density_values[point] *
					     JxW_values[point];
	  };

      cell->get_dof_indices (cell_dof_indices);
      for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
	{
	  final_u_bar(cell_dof_indices[shape_func]) += local_functional1[shape_func];
	  final_v_bar(cell_dof_indices[shape_func]) += local_functional2[shape_func];
	};
    };
};





/*------------------------ IntegrateValueAtOrigin --------------------------------*/


template <int dim>
void IntegratedValueAtOrigin<dim>::compute_functionals (Vector<double> &j1,
							Vector<double> &j2) {
  j1.reinit (dof->n_dofs());
  j2.reinit (dof->n_dofs());

  DoFHandler<dim>::active_cell_iterator cell = dof->begin_active(),
					endc = dof->end();

  Point<dim> origin;

  bool origin_found = false;
  for (; (cell!=endc) && !origin_found; ++cell) 
    {
      for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
	if (cell->vertex(vertex) == origin) 
	  {
	    j1(cell->vertex_dof_index(vertex,0)) = 1;
	    origin_found = true;
	  };
    };

  Assert (origin_found, ExcVertexNotFound());
};






/*------------------------ SeismicSignal --------------------------------*/


template <>
void SeismicSignal<1>::compute_functionals (Vector<double> &,
					    Vector<double> &)
{
  Assert (false, ExcNotImplemented());
};



template <int dim>
void SeismicSignal<dim>::compute_functionals (Vector<double> &j1,
					      Vector<double> &j2) {
  const double y_offset = 1.0;
  const unsigned int n_q_points = quadrature_face->n_quadrature_points;
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  
  j1.reinit (dof->n_dofs());
  j2.reinit (dof->n_dofs());

  DoFHandler<dim>::active_cell_iterator cell, endc;
  DoFHandler<dim>::face_iterator        face;
  cell = dof->begin_active();
  endc = dof->end();

  vector<int> cell_dof_indices (dofs_per_cell);

  FEFaceValues<dim> fe_face_values (*fe, *quadrature_face,
				    UpdateFlags(update_values         |
						update_JxW_values     |
						update_q_points));
  
  for (; cell!=endc; ++cell)
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
	 ++face_no)
      if (face=cell->face(face_no),
	  (face->vertex(0)(1) == y_offset) &&
	  (face->vertex(1)(1) == y_offset))
					 // this is one of the faces we
					 // are interested in, i.e. which
					 // lie on the interesting line
	{
	  fe_face_values.reinit (cell, face_no);
	  const FullMatrix<double>            &shape_values = fe_face_values.
						    get_shape_values ();
	  const vector<double>      &JxW_values (fe_face_values.
						 get_JxW_values());
	  const vector<Point<dim> > &q_points (fe_face_values.get_quadrature_points());

					   // now compute the local integral
					   // \int w(x,t) phi_i(x,y,t) ds
					   // through this line for each
					   // of the basis functions
	  vector<double> local_integral (dofs_per_cell, 0);
	  for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
	    for (unsigned int point=0; point<n_q_points; ++point)
	      local_integral[shape_func] += shape_values(shape_func,point) *
					    (EvaluateSeismicSignal<dim>
					     ::weight(q_points[point], time)) *
					    JxW_values[point];

	  cell->get_dof_indices (cell_dof_indices);
	  for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
	    j1(cell_dof_indices[shape_func]) += local_integral[shape_func];
	};
};




/*------------------------ EarthSurface --------------------------------*/


template <>
void EarthSurface<1>::compute_functionals (Vector<double> &,
					   Vector<double> &)
{
  Assert (false, ExcNotImplemented());
};



template <int dim>
void EarthSurface<dim>::compute_functionals (Vector<double> &j1,
					     Vector<double> &j2) {
  const unsigned int face_dofs = fe->dofs_per_face;
  
  j1.reinit (dof->n_dofs());
  j2.reinit (dof->n_dofs());

  DoFHandler<dim>::active_cell_iterator cell, endc;
  DoFHandler<dim>::face_iterator        face;
  cell = dof->begin_active();
  endc = dof->end();

  vector<int> face_dof_indices (face_dofs);

  for (; cell!=endc; ++cell)
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
	 ++face_no)
      if (face=cell->face(face_no),
	  face->at_boundary())
					 // this is one of the faces we
					 // may be interested in
	{
					   // find out whether it is part of
					   // the boundary portions we are
					   // looking for
	  const double x = face->center()(0),
		       y = face->center()(1);

	  if (!  (((x>0) && (fabs(y) < 500)) ||
		  ((x>0) && (y<0) && (fabs(x+y)<500))))
	    continue;

					   // doubtful for higher
					   // order elements!
	  const double h = face->measure ();
	  
	  face->get_dof_indices (face_dof_indices);
	  for (unsigned int shape_func=0; shape_func<face_dofs; ++shape_func)
					     // also doubtful!
	    j1(face_dof_indices[shape_func]) = h;
	};
};




/*------------------------ SplitSignal --------------------------------*/


template <>
void SplitSignal<1>::compute_functionals (Vector<double> &,
					  Vector<double> &)
{
  Assert (false, ExcInternalError());
};




template <int dim>
void SplitSignal<dim>::compute_functionals (Vector<double> &j1,
					    Vector<double> &j2) {
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q_points = quadrature_face->n_quadrature_points;
  
  j1.reinit (dof->n_dofs());
  j2.reinit (dof->n_dofs());

  if ((time<=1.6) || (time>1.8))
    return;
  
  DoFHandler<dim>::active_cell_iterator cell, endc;
  DoFHandler<dim>::face_iterator        face;
  cell = dof->begin_active();
  endc = dof->end();

  vector<int> dof_indices (fe->dofs_per_cell);
  FEFaceValues<dim> fe_face_values (*fe, *quadrature_face, UpdateFlags(update_values | update_JxW_values));
  
  for (; cell!=endc; ++cell)
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
	 ++face_no)
      if (cell->face(face_no)->center()(0) == 1.5)
					 // this is one of the faces we
					 // may be interested in
	{
	  face=cell->face(face_no);
					   // check whether it really is
	  bool wrong_face = face->center()(1) > 0.0625;
	  if (!wrong_face)
	    for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
	      if (face->vertex(v)(0) != 1.5)
		{
		  wrong_face=true;
		  break;
		};
	  if (wrong_face)
	    continue;

	  fe_face_values.reinit (cell, face_no);
	  const FullMatrix<double> &shape_values = fe_face_values.get_shape_values ();
	  const vector<double>     &JxW_values   = fe_face_values.get_JxW_values();
	  cell->get_dof_indices (dof_indices);

	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      double sum=0;
	      for (unsigned int j=0; j<n_q_points; ++j)
		sum += shape_values(i,j)*JxW_values[j];

					       // since we integrate over each face
					       // twice, add only half of it
	      j1(dof_indices[i]) += sum * time_step / 2;
	    };
	};
};




/* ------------------------------ Split line 1d case ----------------------------- */

template <int dim>
void SplitLine<dim>::compute_endtime_vectors (Vector<double> &,
					      Vector<double> &) {
  Assert (false, ExcNotImplemented ());
};


#if 2 == 1

template <>
void SplitLine<1>::compute_endtime_vectors (Vector<double> &final_u_bar,
					    Vector<double> &final_v_bar) {
  const unsigned int dim = 1;

  const double n_q_points = quadrature->n_quadrature_points;
  const unsigned int dofs_per_cell = fe->dofs_per_cell;

  final_u_bar.reinit (dof->n_dofs());
  final_v_bar.reinit (dof->n_dofs());

  DoFHandler<dim>::active_cell_iterator cell = dof->begin_active (),
					endc = dof->end ();

  FEValues<dim> fe_values (*fe, *quadrature, UpdateFlags(update_values | update_JxW_values));
  vector<int> cell_dof_indices (dofs_per_cell);

  for (; cell!=endc; ++cell)
    {
      if ((cell->vertex(0)(0) < -0.5) ||
	  (cell->vertex(1)(0)  > 0.5))
	continue;
      
      fe_values.reinit (cell);

      const FullMatrix<double>  &shape_values = fe_values.get_shape_values ();
      const vector<double>      &JxW_values (fe_values.get_JxW_values());

      vector<double> local_functional (dofs_per_cell, 0);
      for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
	for (unsigned int point=0; point<n_q_points; ++point) 
	  local_functional[shape_func] += shape_values(shape_func,point) *
					   JxW_values[point];

      cell->get_dof_indices (cell_dof_indices);
      for (unsigned int shape_func=0; shape_func<dofs_per_cell; ++shape_func)
	final_u_bar(cell_dof_indices[shape_func]) += local_functional[shape_func];
    };
};

#endif



/*------------------------ OneBranch1d --------------------------------*/



template <int dim>
void OneBranch1d<dim>::compute_functionals (Vector<double> &j1,
					    Vector<double> &j2) {
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q_points = quadrature->n_quadrature_points;
  
  j1.reinit (dof->n_dofs());
  j2.reinit (dof->n_dofs());
  
				   // take the time step right before 2.5
  if ((time<=2.5-time_step) || (time>2.5))
    return;

  DoFHandler<dim>::active_cell_iterator cell, endc;
  cell = dof->begin_active();
  endc = dof->end();

  vector<int> dof_indices (fe->dofs_per_cell);
  FEValues<dim> fe_values (*fe, *quadrature, UpdateFlags(update_values | update_JxW_values));
  
  for (; cell!=endc; ++cell)
    if ((cell->center()(0) > -0.6) &&
	(cell->center()(0) < -0.4))
      {
	fe_values.reinit (cell);
	const FullMatrix<double>  &shape_values = fe_values.get_shape_values ();
	const vector<double>      &JxW_values   = fe_values.get_JxW_values();
	cell->get_dof_indices (dof_indices);
	
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    double sum=0;
	    for (unsigned int j=0; j<n_q_points; ++j)
	      sum += shape_values(i,j)
		     *JxW_values[j];
	    
					     // since we integrate over each face
					     // twice, add only half of it
	    j1(dof_indices[i]) += sum;
	  };
      };
};



/*------------------------ SecondCrossing --------------------------------*/



template <int dim>
void SecondCrossing<dim>::compute_functionals (Vector<double> &j1,
					       Vector<double> &j2) {
  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q_points = quadrature->n_quadrature_points;
  
  j1.reinit (dof->n_dofs());
  j2.reinit (dof->n_dofs());
  
				   // take the time step right before 2.4
  if ((time<=2.4-time_step) || (time>2.4))
    return;

  DoFHandler<dim>::active_cell_iterator cell, endc;
  cell = dof->begin_active();
  endc = dof->end();

  vector<int> dof_indices (fe->dofs_per_cell);
  FEValues<dim> fe_values (*fe, *quadrature, UpdateFlags(update_values | update_JxW_values));
  
  for (; cell!=endc; ++cell)
    if ((cell->center()(0) > -0.03) &&
	(cell->center()(0) < 0.03))
      {
	fe_values.reinit (cell);
	const FullMatrix<double>  &shape_values = fe_values.get_shape_values ();
	const vector<double>      &JxW_values   = fe_values.get_JxW_values();
	cell->get_dof_indices (dof_indices);
	
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  {
	    double sum=0;
	    for (unsigned int j=0; j<n_q_points; ++j)
	      sum += shape_values(i,j)
		     *JxW_values[j];
	    
	    j1(dof_indices[i]) += sum / time_step;
	  };
      };
};



/*------------------------ HuyghensWave --------------------------------*/



template <int dim>
void HuyghensWave<dim>::compute_functionals (Vector<double> &j1,
					     Vector<double> &j2) {
  j1.reinit (dof->n_dofs());
  j2.reinit (dof->n_dofs());
  
  if ((time < 0.5) || (time > 0.69)) 
    return;
  
  Point<dim> p;
  p(0) = 0.75;
  const Point<dim> evaluation_point (p);

  const DoFHandler<dim>::cell_iterator endc = dof->end(3);
  bool point_found = false;
  for (DoFHandler<dim>::cell_iterator cell=dof->begin(3);
       (cell!=endc) && !point_found; ++cell)
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex) 
      if (cell->vertex(vertex) == evaluation_point)
	{
					   // step down the list of children
					   // until we find a terminal cell
	  DoFHandler<dim>::cell_iterator terminal_cell = cell;
	  while (terminal_cell->has_children())
	    terminal_cell = terminal_cell->child(vertex);
	  
					   // now terminal cell is the right one
	  j1(cell->vertex_dof_index(vertex,0)) = time*time_step;
	  point_found = true;

	  break;
	};
  
  AssertThrow (point_found, ExcInternalError());
};



// explicit specializations

template class DualFunctional<2>;
template class EndEnergy<2>;
template class IntegratedValueAtOrigin<2>;
template class SeismicSignal<2>;
template class EarthSurface<2>;
template class SplitSignal<2>;
template class SplitLine<2>;
template class OneBranch1d<2>;
template class SecondCrossing<2>;
template class HuyghensWave<2>;
/*  $Id$ */


#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_constraints.h>
#include <grid/geometry_info.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <base/quadrature.h>
#include <base/function.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>

#include <cmath>
#include <fstream>



/*--------------------------- EvaluationBase --------------------------*/

template <int dim>
EvaluationBase<dim>::EvaluationBase () :
		u (0),
		v (0),
		tria (0),
		boundary (0),
		dof (0),
		fe (0),
		quadrature (0),
		quadrature_face (0),
		density (0),
		stiffness (0),
		time (0),
		time_step (0),
		step_no (0)
{};


template <int dim>
void EvaluationBase<dim>::reset_timelevel (const TimeStep_Primal<dim> &target) {
  u               = &target.u;
  v               = &target.v;
  tria            = target.tria;
  boundary        = target.parameters.boundary;
  dof             = target.dof_handler;
  fe              = &target.fe;
  quadrature      = &target.quadrature;
  quadrature_face = &target.quadrature_face;
  density         = target.parameters.density;
  stiffness       = target.parameters.stiffness;
  time            = target.time;
  time_step       = (target.timestep_no == 0 ?
		     0 :
		     target.get_backward_timestep());
  step_no         = target.timestep_no;

  base_file_name  = target.parameters.output_directory +
		    "sweep"+int_to_string(target.sweep_no, 2) + "/evaluation/" +
		    int_to_string(step_no,4);
};



template <int dim>
void EvaluationBase<dim>::reset () {};



template <int dim>
void EvaluationBase<dim>::print_final_result (ostream &) {};


template <int dim>
double EvaluationBase<dim>::get_final_result () {
  return 0;
};






/*--------------------------- EvaluateEnergyContent ----------------------*/

template <int dim>
EvaluateEnergyContent<dim>::EvaluateEnergyContent () :
		old_energy (0),
		integrated_outflux (0) {};


template <int dim>
void EvaluateEnergyContent<dim>::reset () {
  old_energy         = 0;
  integrated_outflux = 0;
};



template <int dim>
double EvaluateEnergyContent<dim>::compute_energy (const PartOfDomain pod) const {
  const double y_offset = 300000000;

  DoFHandler<dim>::active_cell_iterator cell, endc;
  cell = dof->begin_active ();
  endc = dof->end ();

  FEValues<dim> fe_values (*fe, *quadrature,
			   UpdateFlags(update_values         |
				       update_gradients      |
				       update_JxW_values     |
				       update_q_points));
  FullMatrix<double>  cell_matrix (fe->dofs_per_cell, fe->dofs_per_cell);
  Vector<double>   local_u (fe->dofs_per_cell);
  Vector<double>   local_v (fe->dofs_per_cell);
  
  vector<double> density_values(quadrature->n_quadrature_points);
  vector<double> stiffness_values(quadrature->n_quadrature_points);

  double total_energy = 0;
  
  for (; cell!=endc; ++cell)
    {
				       // only consider cells in the specified
				       // domain
      switch (pod)
	{
	  case low_atmosphere:
		if (cell->center()(1) >= y_offset)
		  continue;
		break;
	  case high_atmosphere:
		if (cell->center()(1) < y_offset)
		  continue;
		break;
	};
      
      
      fe_values.reinit (cell);
      const FullMatrix<double>             &values    = fe_values.get_shape_values();
      const vector<vector<Tensor<1,dim> > >&gradients = fe_values.get_shape_grads ();
      const vector<double>                 &weights   = fe_values.get_JxW_values ();

      cell->get_dof_values (*u, local_u);
      cell->get_dof_values (*v, local_v);

				       // compute mass matrix
      cell_matrix.clear ();
      density->value_list (fe_values.get_quadrature_points(),
			   density_values);
      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	for (unsigned int i=0; i<fe->dofs_per_cell; ++i) 
	  for (unsigned int j=0; j<fe->dofs_per_cell; ++j)
	    cell_matrix(i,j) += (values(i,point) *
				 values(j,point)) *
				weights[point] *
				density_values[point];

      total_energy += 1./2. * cell_matrix.matrix_norm (local_v);

				       // now for the part with the laplace
				       // matrix
      cell_matrix.clear ();
      stiffness->value_list (fe_values.get_quadrature_points(),
			     stiffness_values);
      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	for (unsigned int i=0; i<fe->dofs_per_cell; ++i) 
	  for (unsigned int j=0; j<fe->dofs_per_cell; ++j)
	    cell_matrix(i,j) += (gradients[i][point] *
				 gradients[j][point]) *
				weights[point] *
				stiffness_values[point];
      total_energy += 1./2. * cell_matrix.matrix_norm (local_u);
    };

  return total_energy;
};




/* ---------------------------- EvaluateIntegratedValueAtOrigin ------------------- */


template <int dim>
void EvaluateIntegratedValueAtOrigin<dim>::print_final_result (ostream &out) {
  out << "    Integrated value of u at origin: "
      << integrated_value << endl;
};



template <int dim>
double EvaluateIntegratedValueAtOrigin<dim>::get_final_result () {
  return integrated_value;
};



template <int dim>
string EvaluateIntegratedValueAtOrigin<dim>::description () const {
  return "integrated value at origin";
};



template <int dim>
void EvaluateIntegratedValueAtOrigin<dim>::reset () {
  integrated_value = 0;
};



template <int dim>
double EvaluateIntegratedValueAtOrigin<dim>::evaluate () {
  DoFHandler<dim>::active_cell_iterator cell = dof->begin_active(),
					endc = dof->end();

  double     value_at_origin = 0;
  Point<dim> origin;

  bool origin_found = false;
  for (; (cell!=endc) && !origin_found; ++cell) 
    {
      for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
	if (cell->vertex(vertex) == origin) 
	  {
	    value_at_origin = (*u)(cell->vertex_dof_index(vertex,0));
	    origin_found = true;
	  };
    };

  Assert (origin_found, ExcVertexNotFound());

  if (time > 0)
    integrated_value += value_at_origin * time_step;
      
  return value_at_origin;
};






/*------------------------- EvaluateSeismicSignal --------------------------*/


template <int dim>
void EvaluateSeismicSignal<dim>::print_final_result (ostream &out) {
  out << "    Integrated seismic signal: " << result << endl;
};



template <int dim>
double EvaluateSeismicSignal<dim>::get_final_result () {
  return result;
};



template <int dim>
string EvaluateSeismicSignal<dim>::description () const {
  return "Integrated seismic signal at top";
};



template <int dim>
void EvaluateSeismicSignal<dim>::reset () {
  result = 0;
};



template <>
double EvaluateSeismicSignal<1>::evaluate ()
{
  Assert (false, ExcNotImplemented());
  return 0;
};



template <int dim>
double EvaluateSeismicSignal<dim>::evaluate () {
  const unsigned int n_q_points = quadrature_face->n_quadrature_points;

  ofstream out((base_file_name + ".seismic").c_str());
  AssertThrow (out, ExcIO());
  
  DoFHandler<dim>::active_cell_iterator cell = dof->begin_active(),
					endc = dof->end();
  double u_integrated=0;
  FEFaceValues<dim> face_values (*fe, *quadrature_face,
				 UpdateFlags(update_values         |
					     update_JxW_values     |
					     update_q_points));
  vector<double>    face_u (fe->dofs_per_face);
  
  for (; cell!=endc; ++cell)
    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
				       // check if face is at top boundary
      if (cell->face(face)->center()(1) == 1.0)
	{
	  face_values.reinit (cell, face);
	  face_values.get_function_values (*u, face_u);
	  const vector<double>      &JxW_values (face_values.get_JxW_values());
	  const vector<Point<dim> > &q_points   (face_values.get_quadrature_points());
	  
	  double local_integral = 0;
	  for (unsigned int point=0; point<n_q_points; ++point)
	    local_integral += face_u[point] *
			      weight (q_points[point], time) *
			      JxW_values[point];
	  u_integrated += local_integral;

					   // output the t and x coordinate
	  out << time
	      << ' '
	      << cell->face(face)->vertex(0)(0)
	      << "  "
	      << (*u)(cell->face(face)->vertex_dof_index(0,0))
	      << endl
	      << time
	      << ' '
	      << cell->face(face)->vertex(1)(0)
	      << "  "
	      << (*u)(cell->face(face)->vertex_dof_index(1,0))
	      << endl
	      << endl;
	};
  AssertThrow (out, ExcIO());
  out.close ();
  
  if (time!=0)
    result += u_integrated*time_step;
  
  return u_integrated;
};




/*------------------------- EvaluateSplitSignal --------------------------*/


template <int dim>
void EvaluateSplitSignal<dim>::print_final_result (ostream &out) {
  out << "    Integrated split signal: " << result << endl;
};



template <int dim>
double EvaluateSplitSignal<dim>::get_final_result () {
  return result;
};



template <int dim>
string EvaluateSplitSignal<dim>::description () const {
  return "Integrated split signal (exact: (2+pi)/(16-pi)=0.010229)";
};



template <int dim>
void EvaluateSplitSignal<dim>::reset () {
  result = 0;
};



template <>
double EvaluateSplitSignal<1>::evaluate ()
{
  Assert (false, ExcNotImplemented());
  return 0;
};



template <int dim>
double EvaluateSplitSignal<dim>::evaluate () {
  if ((time<=1.6) || (time>1.8))
    return 0;

  const unsigned int n_q_points = quadrature_face->n_quadrature_points;
  DoFHandler<dim>::active_cell_iterator cell = dof->begin_active(),
					endc = dof->end();
  double u_integrated=0;
  FEFaceValues<dim> face_values (*fe, *quadrature_face, UpdateFlags(update_values | update_JxW_values));
  vector<double>    face_u (fe->dofs_per_face);
  
  for (; cell!=endc; ++cell)
    for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell; ++face_no)
					 // this is one of the faces we
					 // may be interested in
      if (cell->face(face_no)->center()(0) == 1.5)
	{
	  DoFHandler<dim>::face_iterator face=cell->face(face_no);
					   // check whether it really is
	  bool wrong_face = face->center()(1) > 0.0625;
	  if (!wrong_face)
	    for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
	      if (face->vertex(v)(0) != 1.5)
		{
		  wrong_face=true;
		  break;
		};
	  if (wrong_face)
	    continue;

	  face_values.reinit (cell, face_no);
	  face_values.get_function_values (*u, face_u);
	  const vector<double>      &JxW_values (face_values.get_JxW_values());
	  
	  double local_integral = 0;
	  for (unsigned int point=0; point<n_q_points; ++point)
	    local_integral += face_u[point] *
			      JxW_values[point];
	  u_integrated += local_integral;
	};

				   // note that we integrate over each line twice, so
				   // we divide the result by two
  if (time!=0)
    result += u_integrated*time_step / 2;
  
  return u_integrated;
};




/*------------------------- EvaluateOneBranch1d --------------------------*/


template <int dim>
void EvaluateOneBranch1d<dim>::print_final_result (ostream &out) {
  out << "    One branch integrated: " << result << endl;
};



template <int dim>
double EvaluateOneBranch1d<dim>::get_final_result () {
  return result;
};



template <int dim>
string EvaluateOneBranch1d<dim>::description () const {
  return "One branch integrated (exact: 0.055735)";
};



template <int dim>
void EvaluateOneBranch1d<dim>::reset () {
  result = 0;
};



template <int dim>
double EvaluateOneBranch1d<dim>::evaluate ()
{
  Assert (false, ExcNotImplemented());
  return 0;
};


#if 2 == 1

template <>
double EvaluateOneBranch1d<1>::evaluate () {
  if ((time<=2.5-time_step) || (time>2.5))
    return 0;

  const unsigned int n_q_points = quadrature->n_quadrature_points;
  DoFHandler<1>::active_cell_iterator cell = dof->begin_active(),
				      endc = dof->end();
  double u_integrated=0;
  FEValues<1>  fe_values (*fe, *quadrature, UpdateFlags(update_values | update_JxW_values));
  vector<double> cell_u (fe->dofs_per_cell);
  
  for (; cell!=endc; ++cell)
    if ((cell->center()(0) > -0.6) &&
	(cell->center()(0) < -0.4))
      {
	fe_values.reinit (cell);
	fe_values.get_function_values (*u, cell_u);
	const vector<double>    &JxW_values (fe_values.get_JxW_values());

	double local_integral = 0;
	for (unsigned int point=0; point<n_q_points; ++point)
	  local_integral += cell_u[point] *
			    JxW_values[point];
	u_integrated += local_integral;
      };
  result += u_integrated;
  
  return u_integrated;
};

#endif




/*------------------------- EvaluateSecondCrossing1d --------------------------*/


template <int dim>
void EvaluateSecondCrossing1d<dim>::print_final_result (ostream &out) {
  out << "    Second crossing: " << result << endl;
};



template <int dim>
double EvaluateSecondCrossing1d<dim>::get_final_result () {
  return result;
};



template <int dim>
string EvaluateSecondCrossing1d<dim>::description () const {
  return "Second crossing (exact: 0.011147)";
};



template <int dim>
void EvaluateSecondCrossing1d<dim>::reset () {
  result = 0;
};



template <int dim>
double EvaluateSecondCrossing1d<dim>::evaluate ()
{
  Assert (false, ExcNotImplemented());
  return 0;
};


#if 2 == 1

template <>
double EvaluateSecondCrossing1d<1>::evaluate () {
  if ((time<=2.4-time_step) || (time>2.4))
    return 0;

  const unsigned int n_q_points = quadrature->n_quadrature_points;
  DoFHandler<1>::active_cell_iterator cell = dof->begin_active(),
				      endc = dof->end();
  double u_integrated=0;
  FEValues<1>  fe_values (*fe, *quadrature,
			  UpdateFlags(update_values | update_JxW_values | update_q_points));
  vector<double> cell_u (fe->dofs_per_cell);
  
  for (; cell!=endc; ++cell)
    if ((cell->center()(0) > -0.03) &&
	(cell->center()(0) < 0.03))
      {
	fe_values.reinit (cell);
	fe_values.get_function_values (*u, cell_u);
	const vector<double>    &JxW_values (fe_values.get_JxW_values());
	const vector<Point<1> > &quadrature_points (fe_values.get_quadrature_points());
	
	double local_integral = 0;
	for (unsigned int point=0; point<n_q_points; ++point)
	  local_integral += cell_u[point] *
			    JxW_values[point];
	u_integrated += local_integral;
      };
  result += u_integrated;
  
  return u_integrated;
};

#endif



/*------------------------- EvaluateHuyghensWave --------------------------*/


template <int dim>
void EvaluateHuyghensWave<dim>::print_final_result (ostream &out) {
  out << "    Hughens wave -- weighted time: " << weighted_value / integrated_value << endl;
  out << "                    average      : " << integrated_value << endl;
};



template <int dim>
double EvaluateHuyghensWave<dim>::get_final_result () {
  return weighted_value / integrated_value;
};



template <int dim>
string EvaluateHuyghensWave<dim>::description () const {
  return "Huyghens wave";
};



template <int dim>
void EvaluateHuyghensWave<dim>::reset () {
  integrated_value = weighted_value = 0;
};



template <int dim>
double EvaluateHuyghensWave<dim>::evaluate ()
{
  double     value_at_origin = 0;
  Point<dim> p;
  p(0) = 0.75;
  const Point<dim> evaluation_point (p);

  const DoFHandler<dim>::cell_iterator endc = dof->end(3);
  bool point_found = false;
  for (DoFHandler<dim>::cell_iterator cell=dof->begin(3);
       (cell!=endc) && !point_found; ++cell)
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
      if (cell->vertex(vertex) == evaluation_point)
	{
					   // step down the list of children
					   // until we find a terminal cell
	  DoFHandler<dim>::cell_iterator terminal_cell = cell;
	  while (terminal_cell->has_children())
	    terminal_cell = terminal_cell->child(vertex);
	  
					   // now terminal cell is the right one
	  value_at_origin = (*u)(cell->vertex_dof_index(vertex,0));
	  point_found = true;
	  
	  break;
	};

  AssertThrow (point_found, ExcInternalError());

  if ((time > 0.5) && (time < 0.69)) 
    {
      integrated_value += value_at_origin * time_step;
      weighted_value += value_at_origin * time_step * time;
    };
  
  return value_at_origin;
};


		
// explicit instantiations
template class EvaluationBase<2>;
template class EvaluateEnergyContent<2>;
template class EvaluateIntegratedValueAtOrigin<2>;
template class EvaluateSeismicSignal<2>;
template class EvaluateSplitSignal<2>;
template class EvaluateOneBranch1d<2>;
template class EvaluateSecondCrossing1d<2>;
template class EvaluateHuyghensWave<2>;

/* $Id$ */

#include <base/data_out_base.h>
#include <numerics/histogram.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_handler.h>
#include <grid/geometry_info.h>

#include <cmath>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <iomanip>

#include <numerics/data_out_stack.h>


template <int dim>
TimestepManager<dim>::TimestepManager (const WaveParameters<dim> &parameters) :
		 TimeDependent(TimeDependent::TimeSteppingData(0,1),
			       TimeDependent::TimeSteppingData(0,1),
			       TimeDependent::TimeSteppingData(0,1)),
		 parameters (parameters)
{};



template <int dim>
void TimestepManager<dim>::run_sweep (const unsigned int sweep_no)
{
  SweepInfo      sweep_info;
  SweepData<dim> sweep_data (parameters.write_stacked_data);
  if (parameters.write_stacked_data)
    {
      sweep_data.data_out_stack->declare_data_vector ("u", DataOutStack<dim>::dof_vector);
      sweep_data.data_out_stack->declare_data_vector ("v", DataOutStack<dim>::dof_vector);
      if ((parameters.refinement_strategy == WaveParameters<dim>::dual_estimator)
	  &&
	  (sweep_no >= parameters.initial_energy_estimator_sweeps))
	{
	  sweep_data.data_out_stack->declare_data_vector ("dual_u", DataOutStack<dim>::dof_vector);
	  sweep_data.data_out_stack->declare_data_vector ("dual_v", DataOutStack<dim>::dof_vector);
	};
      if ((sweep_no<parameters.number_of_sweeps-1) ||
	  (parameters.refinement_strategy == WaveParameters<dim>::dual_estimator))
	sweep_data.data_out_stack->declare_data_vector ("est_error", DataOutStack<dim>::cell_vector);
    };
  
  
  deallog << "Sweep " << setw(2) << sweep_no << ':' << endl
       << "---------" << endl;

  for (typename list<EvaluationBase<dim>*>::const_iterator i = parameters.eval_list.begin();
       i != parameters.eval_list.end(); ++i)
    (*i)->reset ();
  
  start_sweep (sweep_no);

				   // attach the present sweep_info object
				   // to all the time steps. also for
				   // the sweep_data object
  for (vector<TimeStepBase*>::iterator timestep=timesteps.begin();
       timestep!=timesteps.end(); ++timestep)
    {
      dynamic_cast<TimeStepBase_Wave<dim>*>(*timestep)->attach_sweep_info (sweep_info);
      dynamic_cast<TimeStepBase_Wave<dim>*>(*timestep)->attach_sweep_data (sweep_data);
    };
  
  solve_primal_problem ();
  deallog << endl;

  if ((parameters.refinement_strategy == WaveParameters<dim>::dual_estimator)
      &&
      (sweep_no >= parameters.initial_energy_estimator_sweeps))
    {
      solve_dual_problem ();
      deallog << endl;
    };
  
  postprocess ();

  if (parameters.write_stacked_data)
    write_stacked_data (*sweep_data.data_out_stack);

  deallog << endl;
  
  if (sweep_no != parameters.number_of_sweeps-1)
    refine_grids ();

  write_statistics (sweep_info);

  end_sweep ();
  
  deallog << endl << endl;
};




template <int dim>
void TimestepManager<dim>::refine_grids () 
{
  deallog << "  Collecting refinement data: " << endl;

  
  const unsigned int n_timesteps = timesteps.size();
      
				       // first collect all the error indicators
  vector<Vector<float> > indicators (n_timesteps);
      
  for (unsigned int i=0; i<n_timesteps; ++i)
    static_cast<TimeStepBase_Wave<dim>*>(timesteps[i])
      ->get_timestep_postprocess().get_tria_refinement_criteria (indicators[i]);


				   // count the number of cells for some
				   // statistics and other things
  unsigned int total_number_of_cells = 0;
  for (unsigned int i=0; i<timesteps.size(); ++i)
    total_number_of_cells += indicators[i].size();



  if (parameters.produce_error_statistics)
    {
      deallog << "    Generating error statistics ";

      vector<double> time_values (timesteps.size());
      for (unsigned int i=0; i<timesteps.size(); ++i)
	time_values[i] = timesteps[i]->get_time();
      
      Histogram error_statistics;
      error_statistics.evaluate (indicators,
				 time_values,
				 parameters.error_statistic_intervals,
				 Histogram::parse_interval_spacing(parameters.error_statistics_scaling));
      error_statistics.write_gnuplot (logfile);

      deallog << endl;
    };


  if (parameters.compare_indicators_globally)
    {

				       // collect all indicators in one
				       // array; delete the old data as soon
				       // as possible, i.e. right after
				       // copying
      Vector<float> all_indicators (total_number_of_cells);
      unsigned int next_index=0;
      for (unsigned int i=0; i<timesteps.size(); ++i)
	{
	  copy (indicators[0].begin(),
		indicators[0].end(),
		&all_indicators(next_index));
	  next_index += (indicators[0].end() - indicators[0].begin());
	  
	  indicators.erase (indicators.begin());
	};

      Assert (next_index==all_indicators.size(),
	      ExcInternalError());

				       /////////////////////////////////////
				       // now find the thresholds for
				       // refinement and coarsening
				       //
      				       // let #all_indicators# be the list
				       // of indicators sorted in *descending*
				       // order. #partial_sums# is the list
				       // of partial sums of #all_indicator#'s
				       // elements from the first to the present
				       // one.
      const double total_error = all_indicators.l1_norm();

      Vector<float> partial_sums(all_indicators.size());
      sort (all_indicators.begin(), all_indicators.end(), greater<double>());
      partial_sum (all_indicators.begin(), all_indicators.end(),
		   partial_sums.begin());

      const Vector<float>::const_iterator
	p = upper_bound (partial_sums.begin(), partial_sums.end(),
			 total_error*(1-parameters.refinement_fraction.second)),
	q = lower_bound (partial_sums.begin(), partial_sums.end(),
			 parameters.refinement_fraction.first*total_error);

      double bottom_threshold = all_indicators(p != partial_sums.end() ?
					       p-partial_sums.begin() :
					       all_indicators.size()-1),
	     top_threshold    = all_indicators(q-partial_sums.begin());
      
      if (bottom_threshold==top_threshold)
	  bottom_threshold = 0.999*top_threshold;

      deallog << "    " << all_indicators.size()
	   << " cells in total."
	   << endl;
      deallog << "    Thresholds are [" << bottom_threshold << "," << top_threshold << "]"
	   << " out of ["
	   << *min_element(all_indicators.begin(),all_indicators.end())
	   << ','
	   << *max_element(all_indicators.begin(),all_indicators.end())
	   << "]. "
	   << endl;
      deallog << "    Expecting "
	   << (all_indicators.size() +
	       (q-partial_sums.begin())*(GeometryInfo<dim>::children_per_cell-1) -
	       (partial_sums.end() - p)/(GeometryInfo<dim>::children_per_cell-1))
	   << " cells in next sweep."
	   << endl;
      deallog << "    Now refining...";
      do_loop (mem_fun (&TimeStepBase_Tria<dim>::init_for_refinement),
	       bind2nd (mem_fun1 (&TimeStepBase_Wave<dim>::refine_grid),
			TimeStepBase_Tria<dim>::RefinementData (top_threshold,
								bottom_threshold)),
	       TimeDependent::TimeSteppingData (0,1),
	       TimeDependent::forward);
      deallog << endl;
    }

  else
				     // refine each time step individually
    {
      deallog << "    Refining each time step separately." << endl;
      
      for (unsigned int timestep=0; timestep<timesteps.size(); ++timestep)
	static_cast<TimeStepBase_Tria<dim>*>(timesteps[timestep])->init_for_refinement();

      unsigned int total_expected_cells = 0;
      
      for (unsigned int timestep=0; timestep<timesteps.size(); ++timestep)
	{
	  TimeStepBase_Wave<dim> *this_timestep
	    = static_cast<TimeStepBase_Wave<dim>*>(timesteps[timestep]);
	    
	  this_timestep->wake_up (0);

					   // copy criteria and delete the old
					   // vector
	  Assert (indicators.size() > 0, ExcInternalError());
	  Vector<float> criteria (indicators[0]);
	  indicators.erase (indicators.begin());
	  
	  const double total_error = criteria.l1_norm();
	  
	  Vector<float> partial_sums(criteria.size());

					   // sort the largest errors to the
					   // beginning of the vector
	  sort (criteria.begin(), criteria.end(), greater<double>());
	  partial_sum (criteria.begin(), criteria.end(),
		       partial_sums.begin());

	  const Vector<float>::const_iterator
	    p = upper_bound (partial_sums.begin(), partial_sums.end(),
			     total_error*(1-parameters.refinement_fraction.second)),
	    q = lower_bound (partial_sums.begin(), partial_sums.end(),
			     parameters.refinement_fraction.first*total_error);

	  double bottom_threshold = criteria(p != partial_sums.end() ?
					     p-partial_sums.begin() :
					     criteria.size()-1),
		 top_threshold    = criteria(q != partial_sums.end() ?
					     q-partial_sums.begin() :
					     criteria.size()-1);
	  
	  if (bottom_threshold==top_threshold)
	    bottom_threshold = 0.999*top_threshold;
	  
	  total_expected_cells += (criteria.size() +
				   (q-partial_sums.begin())*(GeometryInfo<dim>::children_per_cell-1) -
				   (partial_sums.end() - p)/(GeometryInfo<dim>::children_per_cell-1));
	  
	  this_timestep->refine_grid (TimeStepBase_Tria<dim>::RefinementData (top_threshold,
									      bottom_threshold));

	  this_timestep->sleep (0);
	  if (timestep!=0)
	    static_cast<TimeStepBase_Tria<dim>*>(timesteps[timestep-1])->sleep(1);
	};
      
      if (timesteps.size() != 0)
	static_cast<TimeStepBase_Tria<dim>*>(timesteps.back())->sleep(1);


      deallog << "    Got " << total_number_of_cells << " presently, expecting "
	   << total_expected_cells << " for next sweep." << endl;
    };
};




template <int dim>
void TimestepManager<dim>::write_statistics (const SweepInfo &sweep_info) const 
{
				   // write statistics
  if (true)
    {
      deallog << "    Writing statistics for whole sweep.";
      
      deallog << "#  Description of fields" << endl
	   << "#  =====================" << endl
	   << "#  General:"              << endl
	   << "#    time"                << endl;
      
      TimeStep<dim>::write_statistics_descriptions (logfile, parameters);
      deallog << endl << endl;
      
      for (unsigned int timestep=0; timestep<timesteps.size(); ++timestep)
	{
	  deallog << timesteps[timestep]->get_time()
	       << "   ";
	  dynamic_cast<TimeStep<dim>*>
	    (static_cast<TimeStepBase_Wave<dim>*>
	     (timesteps[timestep]))->write_statistics (logfile);
	  deallog << endl;
	};

      AssertThrow (logfile, ExcIO());
      
      deallog << endl;
    };

  
				   // write summary
  if (true)
    {
      deallog << "    Writing summary.";
      
      sweep_info.write_summary (parameters.eval_list,
				logfile);
      AssertThrow (logfile, ExcIO());

      deallog << endl;
    };
};

							

template <int dim>
void TimestepManager<dim>::write_stacked_data (DataOutStack<dim> &data_out_stack) const 
{
  typename DataOutInterface<dim+1>::OutputFormat output_format
    = DataOutInterface<dim+1>::parse_output_format (parameters.output_format);
  
  deallog << "    Writing stacked time steps";
  DataOutBase::EpsFlags eps_flags;
  eps_flags.height_vector = eps_flags.color_vector = 2;
  eps_flags.draw_mesh = false;
  eps_flags.draw_cells = true;
  eps_flags.azimut_angle = 0;
  eps_flags.turn_angle = 0;
  data_out_stack.set_flags (eps_flags);
  data_out_stack.write (logfile, output_format);
  deallog << '.' << endl;
};




//explicit instantiation
template class TimestepManager<2>;
  
/* $Id$ */

#include <base/exceptions.h>
#include <base/function.h>
#include <base/parameter_handler.h>
#include <numerics/histogram.h>
#include <base/data_out_base.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <grid/grid_generator.h>

#include <map>
#include <list>
#include <cmath>


template <int dim>
const string WaveParameters<dim>::initial_value_names ("zero"
						       "|eigenmode"
						       "|bump"
						       "|small bump"
						       "|center-kink"
						       "|shifted bump"
						       "|plateau"
						       "|earthquake");
template <int dim>
const string WaveParameters<dim>::coefficient_names ("unit"
						     "|kink"
						     "|gradient"
						     "|preliminary earth model"
						     "|distorted");
template <int dim>
const string WaveParameters<dim>::boundary_function_names ("wave from left"
							   "|fast wave from left"
							   "|wave from left center"
							   "|wave from left bottom"
							   "|zero");
template <int dim>
const string WaveParameters<dim>::dual_functional_names ("none"
							 "|integrated value at origin"
							 "|seismic signature"
							 "|split signal"
							 "|earth surface"
							 "|split line"
							 "|one branch 1d"
							 "|second crossing"
							 "|Huyghens wave");



DeclException1 (ExcUnknownName,
		string,
		<< "Unknown description string " << arg1);



template <int dim>
class InitialValues {
  public:
    class EigenMode : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double pi = 3.1415926539;
	  return sin(2*pi*p(0))*sin(2*pi*p(1));
	};
    };

    class Bump : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double width = 0.1;
	  const double r2 = p.square();
	  return exp(-r2/width/width) * (r2<width*width ?
					 1-r2/width/width :
					 0);
	};
    };


    class SmallBump : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double width = 0.02;
	  const double r2 = p.square();
	  return exp(-r2/width/width) * (r2<width*width ?
					 1-r2/width/width :
					 0);
	};
    };


    class ShiftedBump : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double width = 0.1;
	  Point<dim> shift;
	  shift(0) = 0.5;
	  const double r2 = (p-shift).square();
	  return exp(-r2/width/width) * (r2<width*width ?
					 1-r2/width/width :
					 0);
	};
    };

    class CenterKink : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double width = 0.1;
	  const double r     = sqrt(p.square());
	  return (r<width ? r/width : (r<2*width ? 2-r/width : 0));
	};
    };

    class Plateau : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double width = 0.1;
	  const double r     = sqrt(p.square());
	  return (r<width ? 1 : 0);
	};
    };


    class Earthquake : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  Point<dim> earthquake_center = p;
	  earthquake_center(1) -= 5500;
	  const double r2 = earthquake_center.square();
	  
	  return (r2<300*300 ? 1-r2/300/300 : 0);
	};
    };
};



template <int dim>
class Coefficients {
  public:
    class Kink : public Function<dim> {
      public:
	inline virtual double value (const Point<dim> &p,
				     const unsigned int) const {
					   // always let the kink be
					   // in direction of the last
					   // variable
	  return 1+8*(p(dim-1)>1./5. ? 1. : 0.);
	};
	
	virtual void value_list (const vector<Point<dim> > &points,
				 vector<double>            &values,
				 const unsigned int) const {
	  Assert (values.size() == points.size(),
		  ExcVectorHasWrongSize(values.size(), points.size()));
	  for (unsigned int i=0; i<points.size(); ++i)
	    values[i]  = this->Kink::value(points[i], 0);
	};

	virtual Tensor<1,dim> gradient (const Point<dim> &p,
					const unsigned int) const {
	  Tensor<1,dim> tmp;
	  if (fabs(p(1)-1./5.) < 1./400.)
	    tmp[1] = 100;
	  return tmp;
	};
	
	virtual void gradient_list (const vector<Point<dim> > &points,
				    vector<Tensor<1,dim> >    &gradients,
				    const unsigned int) const {
	  for (unsigned int i=0; i<points.size(); ++i)
	    gradients[i] = Kink::gradient (points[i], 0);
	};
    };


    class Gradient : public Function<dim> {
      public:
	inline virtual double value (const Point<dim> &p,
				     const unsigned int) const {
	  return 1+8*p(1)*p(1);
	};
	
	virtual void value_list (const vector<Point<dim> > &points,
				 vector<double>            &values,
				 const unsigned int) const {
	  Assert (values.size() == points.size(),
		  ExcVectorHasWrongSize(values.size(), points.size()));
	  for (unsigned int i=0; i<points.size(); ++i)
	    values[i]  = this->Gradient::value(points[i], 0);
	};

	virtual Tensor<1,dim> gradient (const Point<dim> &p,
					const unsigned int) const {
	  Tensor<1,dim> tmp;
	  tmp[1] = 16*p(1);
	  return tmp;
	};
	
	virtual void gradient_list (const vector<Point<dim> > &points,
				    vector<Tensor<1,dim> >    &gradients,
				    const unsigned int) const {
	  for (unsigned int i=0; i<points.size(); ++i)
	    gradients[i] = Gradient::gradient (points[i], 0);
	};
    };



    class PreliminaryEarthModel : public Function<dim> {
      public:
	inline virtual double value (const Point<dim> &p,
				     const unsigned int) const {
	  const double r=sqrt(p.square());
					   // this data just ad hoc, not taken
					   // from the PREM
	  return 10+2.5*(2-r/6371)*(2-r/6371)+20*(r<2000 ? 1 : 0);
	};
	
	virtual void value_list (const vector<Point<dim> > &points,
				 vector<double>            &values,
				 const unsigned int) const {
	  Assert (values.size() == points.size(),
		  ExcVectorHasWrongSize(values.size(), points.size()));
	  for (unsigned int i=0; i<points.size(); ++i)
	    values[i]  = this->PreliminaryEarthModel::value(points[i], 0);
	};

	virtual Tensor<1,dim> gradient (const Point<dim> &p,
					const unsigned int) const {
					   // gradient is derivative with
					   // respect to r times a unit vector
					   // in direction of p
	  Tensor<1,dim> tmp(p);
	  const double r=sqrt(p.square());	  
	  tmp *= 1./r * 2*(10-5*r/6371);
	  return tmp;
	};
	
	virtual void gradient_list (const vector<Point<dim> > &points,
				    vector<Tensor<1,dim> >    &gradients,
				    const unsigned int) const {
	  for (unsigned int i=0; i<points.size(); ++i)
	    gradients[i] = PreliminaryEarthModel::gradient (points[i], 0);
	};
    };

    
    
    class Distorted : public Function<dim> {
      public:
	inline virtual double value (const Point<dim> &p,
				     const unsigned int) const {
	  const double x=p(0),
		       y=p(1);
	  const double pi = 3.1415926539;
	  
	  return (1+0.5*((sin(3*pi*x)>0 ? 1 : 0)+
			 (sin(3*pi*(2*x+y)/sqrt(3)))>0 ? 1 : 0));
	};
	
	virtual void value_list (const vector<Point<dim> > &points,
				 vector<double>            &values,
				 const unsigned int) const {
	  Assert (values.size() == points.size(),
		  ExcVectorHasWrongSize(values.size(), points.size()));
	  for (unsigned int i=0; i<points.size(); ++i)
	    values[i]  = this->Distorted::value(points[i], 0);
	};

	virtual Tensor<1,dim> gradient (const Point<dim> &,
					const unsigned int) const {
					   // return zero, since we don't know
					   // how to do better (regularize?)
	  return Tensor<1,dim>();
	};
	
	virtual void gradient_list (const vector<Point<dim> > &points,
				    vector<Tensor<1,dim> >    &gradients,
				    const unsigned int) const {
	  for (unsigned int i=0; i<points.size(); ++i)
	    gradients[i] = Distorted::gradient (points[i], 0);
	};
    };
};




template <int dim>
class BoundaryValues {
  public:

    class WaveFromLeft_u : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double pi = 3.1415926536;
//	  if ((get_time()<0.4) && (p(0)==0))
	  if (p(0)==0)
	    return sin(pi*get_time()/0.4)*sin(pi*get_time()/0.4);
	  else
	    return 0;
	};
    };

    class WaveFromLeft_v : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double pi = 3.1415926536;
//	  if ((get_time()<0.4) && (p(0)==0))
	  if (p(0)==0)
	    return 2*pi/0.4*sin(pi*get_time()/0.4)*cos(pi*get_time()/0.4);
	  else
	    return 0;
	};
    };


    class FastWaveFromLeft_u : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double pi = 3.1415926536;
	  if ((get_time()<0.2) && (p(0)==0))
	    return sin(pi*get_time()/0.2)*sin(pi*get_time()/0.2);
	  else
	    return 0;
	};
    };

    class FastWaveFromLeft_v : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double pi = 3.1415926536;
	  if ((get_time()<0.2) && (p(0)==0))
	    return 2*pi/0.2*sin(pi*get_time()/0.2)*cos(pi*get_time()/0.2);
	  else
	    return 0;
	};
    };

    
    class WaveFromLeftCenter_u : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double pi = 3.1415926536;
	  if ((0.4 <= p(1)) && (p(1) <= 0.6) && (p(0) <= 0.5))
	    return (p(1)-0.4)*(0.6-p(1)) * sin(pi*get_time()/0.2);
	  else
	    return 0;
	};
    };

    class WaveFromLeftCenter_v : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double pi = 3.1415926536;
	  if ((0.4 <= p(1)) && (p(1) <= 0.6) && (p(0) <= 0.5))
	    return pi/0.2*(p(1)-0.4)*(0.6-p(1)) * cos(pi*get_time()/0.2);
	  else
	    return 0;
	};
    };


    class WaveFromLeftBottom_u : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double pi = 3.1415926536;
	  const double r  = sqrt(p.square());
					   // let the radius of
					   // the excited site be
					   // 50 km
	  const double a  = 5000000;

					   // let the period be
					   // 60 seconds
	  const double period = 60;

	  if ((get_time()>=period) || (r>=a))
	    return 0;

	  const double s = cos(r/a*pi/2)*sin(pi*get_time()/period);
	  return s*s;
	};
    };

    class WaveFromLeftBottom_v : public Function<dim> {
      public:
	virtual double value (const Point<dim> &p,
			      const unsigned int) const {
	  const double pi = 3.1415926536;
	  const double r  = sqrt(p.square());
					   // let the radius of
					   // the excited site be
					   // 50 km
	  const double a  = 5000000;
					   // let the period be
					   // 60 seconds
	  const double period = 60;

	  if ((get_time()>=period) || (r>=a))
	    return 0;
	  else
	    return (2*pi/period*cos(r/a*pi/2)*cos(r/a*pi/2)*
		    sin(pi*get_time()/period)*cos(pi*get_time()/period));
	};
    };

};



template <int dim>
class Boundaries 
{
  public:
    class Ring :  public StraightBoundary<dim> 
    {
      public:
	virtual Point<dim>
	get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const {
	  Point<dim> middle = StraightBoundary<dim>::get_new_point_on_line (line);
	  middle *= sqrt(line->vertex(0).square()) / sqrt(middle.square());
	  return middle;
	};
	
	
	virtual Point<dim>
	get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const {
	  Point<dim> middle = StraightBoundary<dim>::get_new_point_on_quad (quad);
	  middle *= sqrt(quad->vertex(0).square()) / sqrt(middle.square());
	  return middle;
	};
    };
};




template <int dim>
WaveParameters<dim>::WaveParameters () :
		boundary_values_u (0),
		boundary_values_v (0),
		initial_u (0),
		initial_v (0),
		boundary (0),
		density (0),
		stiffness (0),
		dual_functional (0),
		coarse_grid (0)
{};



template <int dim>
WaveParameters<dim>::~WaveParameters ()
{
  delete_parameters ();
};



template <int dim>
void WaveParameters<dim>::delete_parameters ()
{
  if (boundary_values_u)
    delete boundary_values_u;
  boundary_values_u = 0;

  if (boundary_values_v)
    delete boundary_values_v;
  boundary_values_v = 0;

  if (initial_u)
    delete initial_u;
  initial_u = 0;

  if (initial_v)
    delete initial_v;
  initial_v = 0;

  if (boundary)
    delete boundary;
  boundary = 0;

  if (density)
    delete density;
  density = 0;

  if (stiffness)
    delete stiffness;
  stiffness = 0;

  if (dual_functional)
    delete dual_functional;
  dual_functional = 0;

  if (coarse_grid)
    delete coarse_grid;
  coarse_grid = 0;

  				   // free memory used by the evaluation
				   // objects
  for (typename list<EvaluationBase<dim>*>::iterator i=eval_list.begin();
       i!=eval_list.end(); ++i)
    delete *i;
  eval_list.erase (eval_list.begin(), eval_list.end());
};




template <int dim>
void WaveParameters<dim>::set_initial_functions (const string &u_name,
						 const string &v_name) {
  Assert (initial_u==0, ExcInternalError());
  Assert (initial_v==0, ExcInternalError());
  
  const string names[2] = {u_name, v_name};
  Function<dim> *functions[2];
  
  for (unsigned int i=0; i<2; ++i)
    {
      if (names[i]=="eigenmode")
	functions[i] = new InitialValues<dim>::EigenMode();
      else
	if (names[i]=="zero")
	  functions[i] = new ZeroFunction<dim>();
	else
	  if (names[i]=="center-kink")
	    functions[i] = new InitialValues<dim>::CenterKink();
	  else
	    if (names[i]=="bump")
	      functions[i] = new InitialValues<dim>::Bump();
	    else
	      if (names[i]=="small bump")
		functions[i] = new InitialValues<dim>::SmallBump();
	      else
		if (names[i]=="shifted bump")
		  functions[i] = new InitialValues<dim>::ShiftedBump();
		else
		  if (names[i]=="plateau")
		    functions[i] = new InitialValues<dim>::Plateau ();
		  else
		    if (names[i]=="earthquake")
		      functions[i] = new InitialValues<dim>::Earthquake ();
		    else
		      AssertThrow (false, ExcUnknownName(names[i]));
    };

  initial_u = functions[0];
  initial_v = functions[1];
};





template <int dim>
void WaveParameters<dim>::set_coefficient_functions (const string &name) {  
  Assert (density==0, ExcInternalError());
  Assert (stiffness==0, ExcInternalError());

  density = new ConstantFunction<dim>(1);
  density_constant = true;
  
  if (name=="kink")
    {
      stiffness = new Coefficients<dim>::Kink();
      stiffness_constant = false;
    }
  else
    if (name=="gradient")
      {
	stiffness = new Coefficients<dim>::Gradient();
	stiffness_constant = false;
      }
    else
      if (name=="unit")
	{
	  stiffness = new ConstantFunction<dim>(1);
	  stiffness_constant = true;
	}
      else
	if (name=="preliminary earth model")
	  {
	    stiffness = new Coefficients<dim>::PreliminaryEarthModel();
	    stiffness_constant = false;
	  }
	else
	  if (name=="distorted")
	    {
	      stiffness = new Coefficients<dim>::Distorted();
	      stiffness_constant = false;
	  }
	  else
	    AssertThrow (false, ExcUnknownName (name));
};



template <int dim>
void WaveParameters<dim>::set_boundary_functions (const string &name) {
  Assert (boundary_values_u==0, ExcInternalError());
  Assert (boundary_values_v==0, ExcInternalError());
  
  if (name=="wave from left") 
    {
      boundary_values_u = new BoundaryValues<dim>::WaveFromLeft_u ();
      boundary_values_v = new BoundaryValues<dim>::WaveFromLeft_v ();
    }
  else
    if (name=="fast wave from left") 
      {
	boundary_values_u = new BoundaryValues<dim>::FastWaveFromLeft_u ();
	boundary_values_v = new BoundaryValues<dim>::FastWaveFromLeft_v ();
      }
    else
      if (name=="wave from left center")
	{
	  boundary_values_u = new BoundaryValues<dim>::WaveFromLeftCenter_u ();
	  boundary_values_v = new BoundaryValues<dim>::WaveFromLeftCenter_v ();
	}
      else
	if (name=="wave from left bottom")
	  {
	    boundary_values_u = new BoundaryValues<dim>::WaveFromLeftBottom_u ();
	    boundary_values_v = new BoundaryValues<dim>::WaveFromLeftBottom_v ();
	  }
	else
	  if (name=="zero")
	    {
	      boundary_values_u = new ZeroFunction<dim>();
	      boundary_values_v = new ZeroFunction<dim>();
	    }
	  else
	    AssertThrow (false, ExcUnknownName (name));
};



template <int dim>
void WaveParameters<dim>::make_eval_list (const string &names) {
  Assert (eval_list.size()==0, ExcInternalError());
  string split_list = names;

  while (split_list.length())
    {
      string name;
      name = split_list;
      
      if (name.find(",") != string::npos)
	{
	  name.erase (name.find(","), string::npos);
	  split_list.erase (0, split_list.find(",")+1);
	}
      else
	split_list = "";

      while (name[0] == ' ')
	name.erase (0,1);
      while (name[name.length()-1] == ' ')
	name.erase (name.length()-1, 1);

      if (name == "integrated value at origin")
	eval_list.push_back (new EvaluateIntegratedValueAtOrigin<dim>());
      else
	if (name == "seismic signature")
	  eval_list.push_back (new EvaluateSeismicSignal<dim>());
	else
	  if (name == "split signal")
	    eval_list.push_back (new EvaluateSplitSignal<dim>());
	  else
	    if (name == "one branch 1d")
	      eval_list.push_back (new EvaluateOneBranch1d<dim>());
	    else
	      if (name == "second crossing")
		eval_list.push_back (new EvaluateSecondCrossing1d<dim>());
	      else
		if (name == "Huyghens wave")
		  eval_list.push_back (new EvaluateHuyghensWave<dim>());
		else
		  AssertThrow (false, ExcUnknownName (name));
    };
};




template <int dim>
void WaveParameters<dim>::set_dual_functional (const string &name) {
  Assert (dual_functional==0, ExcInternalError());
  if (name == "none")
    dual_functional = new DualFunctional<dim>();
  else
    if (name == "integrated value at origin")
      dual_functional = new IntegratedValueAtOrigin<dim> ();
    else
      if (name == "seismic signature")
	dual_functional = new SeismicSignal<dim> ();
      else
	if (name == "split signal")
	  dual_functional = new SplitSignal<dim> ();
	else
	  if (name == "earth surface")
	    dual_functional = new EarthSurface<dim> ();
	  else
	    if (name == "split line")
	      dual_functional = new SplitLine<dim> ();
	    else
	      if (name == "one branch 1d")
		dual_functional = new OneBranch1d<dim> ();
	      else
		if (name == "second crossing")
		  dual_functional = new SecondCrossing<dim> ();
		else
		  if (name == "Huyghens wave")
		    dual_functional = new HuyghensWave<dim> ();
		  else
		    AssertThrow (false, ExcUnknownName (name));
};




#if 2 == 1

template <>
void WaveParameters<1>::make_coarse_grid (const string &name) {
  const unsigned int dim = 1;
  
  coarse_grid = new Triangulation<dim>(MeshSmoothing(smoothing_on_refinement |
						     eliminate_refined_inner_islands));

  if (name == "line")
    GridGenerator::hyper_cube (*coarse_grid, -1, 1);
  else
    if (name == "split line") 
      {
	const Point<1> vertices[4] = { Point<1>(-1.),
				       Point<1>(-1./3.),
				       Point<1>(1./3.),
				       Point<1>(1.) };
	vector<CellData<1> > cells (3, CellData<1>());
	cells[0].vertices[0] = 0;
	cells[0].vertices[1] = 1;
	cells[0].material_id = 0;

	cells[1].vertices[0] = 1;
	cells[1].vertices[1] = 2;
	cells[1].material_id = 0;

	cells[2].vertices[0] = 2;
	cells[2].vertices[1] = 3;
	cells[2].material_id = 0;

	coarse_grid->create_triangulation (vector<Point<1> >(&vertices[0],
							     &vertices[4]),
					   cells,
					   SubCellData());

					 // refine two of the three cells
	Triangulation<dim>::active_cell_iterator cell = coarse_grid->begin_active();
	(++cell)->set_refine_flag ();
	(++cell)->set_refine_flag ();
	coarse_grid->execute_coarsening_and_refinement ();

					 // refine the level 1 cells 
					 // twice more
	for (int k=0; k<2; ++k)
	  {
	    for (cell=coarse_grid->begin_active(); cell!=coarse_grid->end(); ++cell)
	      if (cell->level() == k+1)
		cell->set_refine_flag ();
	    coarse_grid->execute_coarsening_and_refinement ();
	  };
      }
    else
      AssertThrow (false, ExcParameterNotInList(name));
  
  coarse_grid->refine_global (initial_refinement);
};

#endif



#if 2 == 2

template <>
void WaveParameters<2>::make_coarse_grid (const string &name) {
  const unsigned int dim=2;

  map<string,InitialMesh> initial_mesh_list;
  initial_mesh_list["split channel bottom"] = split_channel_bottom;
  initial_mesh_list["split channel left"]   = split_channel_left;
  initial_mesh_list["split channel right"]  = split_channel_right;
  initial_mesh_list["uniform channel"] = uniform_channel;
  initial_mesh_list["square"]          = square;
  initial_mesh_list["ring"]            = ring;
  initial_mesh_list["earth"]           = earth;
  initial_mesh_list["seismic square"]  = seismic_square;
  AssertThrow (initial_mesh_list.find(name) != initial_mesh_list.end(),
	       ExcParameterNotInList(name));

  const InitialMesh initial_mesh = initial_mesh_list[name];

  coarse_grid = new Triangulation<dim>(MeshSmoothing(smoothing_on_refinement |
						     eliminate_refined_inner_islands));
  
  switch (initial_mesh) 
    {
      case uniform_channel:
      case split_channel_bottom:
      case split_channel_left:
      case split_channel_right:
      {
	const Point<dim> vertices[8] = { Point<dim> (0,0),
					 Point<dim> (1,0),
					 Point<dim> (1,1),
					 Point<dim> (0,1),
					 Point<dim> (2,0),
					 Point<dim> (2,1),
					 Point<dim> (3,0),
					 Point<dim> (3,1)  };
	const int cell_vertices[3][4] = {{0, 1, 2, 3},
					 {1, 4, 5, 2},
					 {4, 6, 7, 5}};
	
	vector<CellData<dim> > cells (3, CellData<dim>());
	
	for (unsigned int i=0; i<3; ++i) 
	  {
	    for (unsigned int j=0; j<4; ++j)
	      cells[i].vertices[j] = cell_vertices[i][j];
	    cells[i].material_id = 0;
	  };
	
	SubCellData boundary_info;
	if ((boundary_conditions == wave_from_left) ||
	    (boundary_conditions == fast_wave_from_left))
	  {
	    for (unsigned int i=0; i<6; ++i)
	      {
		boundary_info.boundary_lines.push_back (CellData<1>());
						 // use Neumann boundary
						 // conditions at top
						 // and bottom of channel
		boundary_info.boundary_lines.back().material_id = 1;
	      };
	    
	    boundary_info.boundary_lines[0].vertices[0] = 0;
	    boundary_info.boundary_lines[0].vertices[1] = 1;
	    boundary_info.boundary_lines[1].vertices[0] = 1;
	    boundary_info.boundary_lines[1].vertices[1] = 4;
	    boundary_info.boundary_lines[2].vertices[0] = 4;
	    boundary_info.boundary_lines[2].vertices[1] = 6;
	    boundary_info.boundary_lines[3].vertices[0] = 3;
	    boundary_info.boundary_lines[3].vertices[1] = 2;
	    boundary_info.boundary_lines[4].vertices[0] = 2;
	    boundary_info.boundary_lines[4].vertices[1] = 5;
	    boundary_info.boundary_lines[5].vertices[0] = 5;
	    boundary_info.boundary_lines[5].vertices[1] = 7;      
	  };

	if (boundary_conditions == wave_from_left_bottom)
	  {
					     // use Neumann bc at left
					     // (mirror condition)
	    boundary_info.boundary_lines.push_back (CellData<1>());
	    boundary_info.boundary_lines.back().material_id = 1;
	    boundary_info.boundary_lines[0].vertices[0] = 0;
	    boundary_info.boundary_lines[0].vertices[1] = 3;
	  };
	
	coarse_grid->create_triangulation (vector<Point<dim> >(&vertices[0],
							       &vertices[8]),
					   cells, boundary_info);
	
	if (initial_refinement >= 1) 
	  {
	    coarse_grid->refine_global (1);

	    switch (initial_mesh)
	      {
		case split_channel_bottom:
		{
		  Triangulation<dim>::active_cell_iterator cell;
		  cell = coarse_grid->begin_active();
		  (cell++)->set_refine_flag ();
		  (cell++)->set_refine_flag ();
		  ++cell; ++cell;
		  (cell++)->set_refine_flag ();
		  (cell++)->set_refine_flag ();
		  ++cell; ++cell;
		  (cell++)->set_refine_flag ();
		  (cell++)->set_refine_flag ();
		  coarse_grid->execute_coarsening_and_refinement ();

		  coarse_grid->refine_global (initial_refinement-1);

		  break;
		};

		case split_channel_left:
		case split_channel_right:
		{
		  coarse_grid->refine_global (1);
		  for (unsigned int i=0; i<2; ++i)
		    {
		      Triangulation<dim>::active_cell_iterator
			cell = coarse_grid->begin_active();

		      for (; cell!=coarse_grid->end(); ++cell)
			if (((cell->center()(0) >= 1) &&
			     (initial_mesh == split_channel_right)) ||
			    ((cell->center()(0) <= 1) &&
			     (initial_mesh == split_channel_left)))
			  cell->set_refine_flag ();
		      coarse_grid->execute_coarsening_and_refinement ();
		    };

		  if (initial_refinement > 4)
		    coarse_grid->refine_global (initial_refinement-4);

		  break;
		};
		  
		
		case uniform_channel:
		{
		  coarse_grid->refine_global (initial_refinement-1);
		  break;
		};


		default:
		      Assert (false, ExcInternalError());
	      };
	  };
	break;
      };

      
      case square:
      case seismic_square:
      {
	GridGenerator::hyper_cube (*coarse_grid, -1, 1);
	if (initial_mesh==seismic_square)
	  coarse_grid->begin_active()->face(2)->set_boundary_indicator(1);

	coarse_grid->refine_global (initial_refinement);

	break;
      };

      case earth:
      {
					 // create ball
	GridGenerator::hyper_ball (*coarse_grid, Point<dim>(), 6371);

	if (boundary)
	  delete boundary;
	
					 // set all boundary to Neumann type
	Triangulation<dim>::active_face_iterator face;
	for (face=coarse_grid->begin_active_face();
	     face != coarse_grid->end_face();
	     ++face)
	  if (face->at_boundary())
	    face->set_boundary_indicator (1);

	const Point<dim> origin;
	boundary = new HyperBallBoundary<dim>(origin, 6371);
					 // set boundary. note that only
					 // id 1 is used
	coarse_grid->set_boundary (1, *boundary);

	coarse_grid->refine_global (initial_refinement);

	break;
      };

      case ring:
      {
	const double radius = 1.;
	const double a = radius/2;
	const Point<2> vertices[8] = { Point<2>(-1,-1)*(radius/sqrt(2)),
				       Point<2>(+1,-1)*(radius/sqrt(2)),
				       Point<2>(-1,-1)*(radius/sqrt(2)*a),
				       Point<2>(+1,-1)*(radius/sqrt(2)*a),
				       Point<2>(-1,+1)*(radius/sqrt(2)*a),
				       Point<2>(+1,+1)*(radius/sqrt(2)*a),
				       Point<2>(-1,+1)*(radius/sqrt(2)),
				       Point<2>(+1,+1)*(radius/sqrt(2)) };
	
	const int cell_vertices[4][4] = {{0, 1, 3, 2},
					 {0, 2, 4, 6},
					 {1, 7, 5, 3},
					 {6, 4, 5, 7}};
	
	vector<CellData<2> > cells (4, CellData<2>());
	
	for (unsigned int i=0; i<4; ++i) 
	  {
	    for (unsigned int j=0; j<4; ++j)
	      cells[i].vertices[j] = cell_vertices[i][j];
	    cells[i].material_id = 0;
	  };
  
	coarse_grid->create_triangulation (vector<Point<2> >(&vertices[0],
							     &vertices[8]),
					   cells,
					   SubCellData());
	if (boundary)
	  delete boundary;
	boundary = new Boundaries<dim>::Ring();
	coarse_grid->set_boundary (0, *boundary);

	coarse_grid->refine_global (initial_refinement);
	
	break;
      };
      
      default:
	    Assert (false, ExcInternalError());
    };
};

#endif


#if 2 == 3

template <>
void WaveParameters<3>::make_coarse_grid (const string &name) {
  const unsigned int dim=3;

  map<string,InitialMesh> initial_mesh_list;
  initial_mesh_list["square"]          = square;
  initial_mesh_list["earth"]           = earth;
  initial_mesh_list["seismic square"]  = seismic_square;
  AssertThrow (initial_mesh_list.find(name) != initial_mesh_list.end(),
	       ExcParameterNotInList(name));

  const InitialMesh initial_mesh = initial_mesh_list[name];

  coarse_grid = new Triangulation<dim>(MeshSmoothing(smoothing_on_refinement |
						     eliminate_refined_inner_islands));
  
  switch (initial_mesh) 
    {
      case square:
      case seismic_square:
      {
	GridGenerator::hyper_cube (*coarse_grid, -1, 1);
	if (initial_mesh==seismic_square)
	  coarse_grid->begin_active()->face(2)->set_boundary_indicator(1);

	coarse_grid->refine_global (initial_refinement);

	break;
      };

      case earth:
      {
					 // create ball
	GridGenerator::hyper_ball (*coarse_grid, Point<dim>(), 6371);

	if (boundary)
	  delete boundary;
	
					 // set all boundary to Neumann type
	Triangulation<dim>::active_face_iterator face;
	for (face=coarse_grid->begin_active_face();
	     face != coarse_grid->end_face();
	     ++face)
	  if (face->at_boundary())
	    face->set_boundary_indicator (1);

	const Point<dim> origin;
	boundary = new HyperBallBoundary<dim>(origin, 6371);
					 // set boundary. note that only
					 // id 1 is used
	coarse_grid->set_boundary (1, *boundary);

	coarse_grid->refine_global (initial_refinement);

	break;
      };

      default:
	    AssertThrow (false, ExcInternalError());
	    break;
    };
};

#endif




template <int dim>
void WaveParameters<dim>::declare_parameters (ParameterHandler &prm) 
{
  prm.enter_subsection ("Grid");
  if (true) {
    prm.declare_entry ("Initial refinement", "0", Patterns::Integer());
    prm.declare_entry ("Coarse mesh", "uniform channel",
		       Patterns::Selection ("uniform channel|split channel bottom|"
					   "split channel left|split channel right|"
					   "square|line|split line|ring|"
					   "seismic square|temperature-square|"
					   "temperature-testcase|random|earth"));
    prm.enter_subsection ("Refinement");
    if (true) {
      prm.declare_entry ("Refinement fraction", "0.95",
			 Patterns::Double());
      prm.declare_entry ("Coarsening fraction", "0.02",
			 Patterns::Double());
      prm.declare_entry ("Compare indicators globally", "true", Patterns::Bool());
      prm.declare_entry ("Maximum refinement", "0", Patterns::Integer());
      prm.declare_entry ("Adapt mesh to dual solution", "true",
			 Patterns::Bool());
      prm.declare_entry ("Primal to dual weight", "1.0",
			 Patterns::Double());
      prm.declare_entry ("Initial energy estimator sweeps", "0",
			 Patterns::Integer());
    };
    prm.leave_subsection ();
     
    prm.enter_subsection ("Mesh smoothing");
    if (true) {
      prm.declare_entry ("Top cell number deviation", "0.1", Patterns::Double());
      prm.declare_entry ("Bottom cell number deviation", "0.03", Patterns::Double());
      prm.declare_entry ("Cell number correction steps", "2", Patterns::Integer());
    };
    prm.leave_subsection ();
  };
  prm.declare_entry ("Renumber dofs", "false", Patterns::Bool());
  prm.leave_subsection ();

  prm.enter_subsection ("Equation data");
  if (true) {
    prm.declare_entry ("Coefficient", "unit", Patterns::Selection(coefficient_names));
    prm.declare_entry ("Initial u", "zero", Patterns::Selection (initial_value_names));
    prm.declare_entry ("Initial v", "zero", Patterns::Selection (initial_value_names));
    prm.declare_entry ("Boundary", "wave from left",
		       Patterns::Selection (boundary_function_names));
  };
  prm.leave_subsection ();

  prm.enter_subsection ("Discretization");
  prm.declare_entry ("Primal FE", "linear",
		     Patterns::Selection ("linear|quadratic|cubic|quartic"));
  prm.declare_entry ("Dual FE", "linear",
		     Patterns::Selection ("linear|quadratic|cubic|quartic"));

  prm.enter_subsection ("Time stepping");
  prm.declare_entry ("Primal method", "fractional step",
		     Patterns::Selection ("theta|fractional step"));
  prm.declare_entry ("Dual method", "fractional step",
		     Patterns::Selection ("theta|fractional step"));
  prm.declare_entry ("Theta", "0.5", Patterns::Double());
  prm.declare_entry ("Time step", "0.1", Patterns::Double());
  prm.declare_entry ("End time", "1",  Patterns::Double());
  prm.leave_subsection ();
  prm.leave_subsection ();

  prm.enter_subsection ("Solver");
  prm.declare_entry ("Preconditioning", "none",
		     Patterns::Selection ("none|jacobi|sor|ssor"));
  prm.declare_entry ("Extrapolate old solutions", "true",
		     Patterns::Bool());
  prm.leave_subsection ();

  prm.enter_subsection ("Output");
  prm.declare_entry ("Format", "gnuplot",
		     Patterns::Selection(DataOutInterface<dim>::get_output_format_names()));
  prm.declare_entry ("Directory", "data");
  prm.declare_entry ("Directory for temporaries", "data/tmp");
  prm.declare_entry ("Write solutions", "all sweeps",
		     Patterns::Selection ("never|all sweeps|last sweep only"));
  prm.declare_entry ("Write stacked time steps", "false", Patterns::Bool());
  prm.declare_entry ("Write stacked interval", "1", Patterns::Integer());
  prm.declare_entry ("Write steps interval", "1", Patterns::Integer());
  prm.declare_entry ("Write error as cell data", "true", Patterns::Bool());
  prm.enter_subsection ("Error statistics");
  prm.declare_entry ("Produce error statistics", "false", Patterns::Bool());
  prm.declare_entry ("Number of intervals", "10", Patterns::Integer());
  prm.declare_entry ("Interval spacing", "linear",
		     Patterns::Selection(Histogram::get_interval_spacing_names()));
  prm.leave_subsection ();
  prm.leave_subsection ();


  prm.enter_subsection ("Goal");
  prm.declare_entry ("Goal", "none",
		     Patterns::Selection (dual_functional_names));
  prm.declare_entry ("Evaluate", "");
  prm.leave_subsection ();  


  prm.declare_entry ("Refinement criterion", "energy estimator",
 		     Patterns::Selection ("energy estimator|dual estimator"));
  prm.declare_entry ("Sweeps", "3", Patterns::Integer());
};



template <int dim>
void WaveParameters<dim>::parse_parameters (ParameterHandler &prm) {
				   // declare some maps for convenience,
				   // to avoid those annoying if then else
				   // clauses...
  map<string,BoundaryConditions> boundary_conditions_list;
  boundary_conditions_list["wave from left"]        = wave_from_left;
  boundary_conditions_list["fast wave from left"]   = fast_wave_from_left;
  boundary_conditions_list["wave from left center"] = wave_from_left_center;
  boundary_conditions_list["wave from left bottom"] = wave_from_left_bottom;
  boundary_conditions_list["zero"] = zero;
  
  map<string,Preconditioning> preconditioning_list;
  preconditioning_list["jacobi"] = jacobi;
  preconditioning_list["sor"]    = sor;
  preconditioning_list["ssor"]   = ssor;
  preconditioning_list["none"]   = no_preconditioning;
  
  map<string,WriteStrategy> write_strategy_list;
  write_strategy_list["never"] = never;
  write_strategy_list["all sweeps"] = all_sweeps;
  write_strategy_list["last sweep only"] = last_sweep_only;
  
  
  prm.enter_subsection ("Grid");
  initial_refinement = prm.get_integer ("Initial refinement");
				   // don't make the grid here already, since
				   // it may depend on the chosen boundary
				   // conditions (which need some boundary
				   // flags to be set), etc.

  prm.enter_subsection ("Refinement");
  {
    refinement_fraction.first   = prm.get_double ("Refinement fraction");
    refinement_fraction.second  = prm.get_double ("Coarsening fraction");
    compare_indicators_globally = prm.get_bool ("Compare indicators globally");
    maximum_refinement          = prm.get_integer ("Maximum refinement");
    adapt_mesh_to_dual_solution = prm.get_bool ("Adapt mesh to dual solution");
    primal_to_dual_weight       = prm.get_double ("Primal to dual weight");
    initial_energy_estimator_sweeps = prm.get_integer("Initial energy estimator sweeps");
  };
  prm.leave_subsection ();

  prm.enter_subsection ("Mesh smoothing");
  {
    cell_number_corridor.first  = prm.get_double ("Top cell number deviation");
    cell_number_corridor.second = prm.get_double ("Bottom cell number deviation");
    cell_number_correction_steps= prm.get_integer ("Cell number correction steps");
  };
  prm.leave_subsection ();

  renumber_dofs               = prm.get_bool ("Renumber dofs");
  prm.leave_subsection ();

  prm.enter_subsection ("Equation data");
  set_coefficient_functions (prm.get("Coefficient"));
  set_initial_functions (prm.get("Initial u"), prm.get("Initial v"));
  boundary_conditions = boundary_conditions_list[prm.get("Boundary")];
  set_boundary_functions (prm.get("Boundary"));  
  Assert (boundary_conditions_list.find(prm.get("Boundary")) !=
	  boundary_conditions_list.end(),
	  ExcParameterNotInList(prm.get("Boundary")));
  prm.leave_subsection ();
  
  prm.enter_subsection ("Discretization");
  primal_fe = prm.get("Primal FE");
  dual_fe = prm.get("Dual FE");
  prm.enter_subsection ("Time stepping");
  theta    = prm.get_double ("Theta");
  time_step= prm.get_double ("Time step");
  end_time = prm.get_double ("End time");
  prm.leave_subsection ();
  prm.leave_subsection ();

  prm.enter_subsection ("Solver");
  preconditioning = preconditioning_list[prm.get("Preconditioning")];
  Assert (preconditioning_list.find(prm.get("Preconditioning")) !=
	  preconditioning_list.end(),
	  ExcParameterNotInList(prm.get("Preconditioning")));
  extrapolate_old_solutions = prm.get_bool ("Extrapolate old solutions");
  prm.leave_subsection ();
  
  prm.enter_subsection ("Output");
  output_format = prm.get("Format");
  output_directory = prm.get("Directory");
  if (output_directory[output_directory.size()-1] != '/')
    output_directory += '/';
  tmp_directory = prm.get ("Directory for temporaries");
  if (tmp_directory[tmp_directory.size()-1] != '/')
    tmp_directory += '/';
  write_solution_strategy = write_strategy_list[prm.get("Write solutions")];
  Assert (write_strategy_list.find(prm.get("Write solutions")) !=
	  write_strategy_list.end(),
	  ExcParameterNotInList(prm.get("Write solutions")));
  write_stacked_data       = prm.get_bool ("Write stacked time steps");
  write_stacked_interval   = prm.get_integer ("Write stacked interval");
  write_steps_interval     = prm.get_integer ("Write steps interval");
  write_error_as_cell_data = prm.get_bool ("Write error as cell data");
  prm.enter_subsection ("Error statistics");
  produce_error_statistics = prm.get_bool ("Produce error statistics");
  error_statistic_intervals= prm.get_integer ("Number of intervals");
  error_statistics_scaling = prm.get ("Interval spacing");
  prm.leave_subsection ();
  prm.leave_subsection ();


  prm.enter_subsection ("Goal");
  set_dual_functional (prm.get("Goal"));
  make_eval_list (prm.get("Evaluate"));
  prm.leave_subsection ();

  
  
  if (prm.get("Refinement criterion")=="energy estimator")
    refinement_strategy = energy_estimator;
  else
    refinement_strategy = dual_estimator;

  number_of_sweeps = prm.get_integer ("Sweeps");

				   // now that we know everything, we can make
				   // the grid
  prm.enter_subsection ("Grid");
  make_coarse_grid (prm.get("Coarse mesh"));
  prm.leave_subsection ();
};




// explicit instantiations
template class WaveParameters<2>;
/* $Id$ */

#include <numerics/data_out_stack.h>
#include <dofs/dof_handler.h>  //??
#include <lac/vector.h>


template <int dim>
SweepData<dim>::SweepData (const bool use_data_out_stack) 
{
  if (use_data_out_stack)
    data_out_stack = new DataOutStack<dim>();
  else
    data_out_stack = 0;
};



template <int dim>
SweepData<dim>::~SweepData () 
{
  if (data_out_stack != 0)
    delete data_out_stack;
  data_out_stack = 0;
};




// explicit instantiations
template class SweepData<2>;
/* $Id$ */


#include <iomanip>
#include <ctime>


SweepInfo::Data &
SweepInfo::get_data () 
{
  return data;
};



SweepInfo::Timers &
SweepInfo::get_timers () 
{
  return timers;
};



template <int dim>
void
SweepInfo::write_summary (const list<EvaluationBase<dim>*> &eval_list,
			  ostream &out) const
{
  out << "Summary of this sweep:" << endl
      << "======================" << endl
      << endl;

  out << "  Accumulated number of cells: " << data.cells       << endl
      << "  Acc. number of primal dofs : " << data.primal_dofs << endl
      << "  Acc. number of dual dofs   : " << data.dual_dofs   << endl
      << "  Accumulated error          : " << data.accumulated_error       << endl;
  
  if (eval_list.size() != 0)
    {
      out << endl;
      out << "  Evaluations:" << endl
	  << "  ------------" << endl;
      
      for (typename list<EvaluationBase<dim>*>::const_iterator i = eval_list.begin();
	   i != eval_list.end(); ++i)

      (*i)->print_final_result (out);
    };
  
  time_t  time1= time (0);
  tm     *time = localtime(&time1); 
  out << "  Time tag: "
      << time->tm_year+1900 << "/"
      << time->tm_mon+1 << "/"
      << time->tm_mday << ' '
      << int_to_string (time->tm_hour, 2) << ":"
      << int_to_string (time->tm_min, 2) << ":"
      << int_to_string (time->tm_sec, 2) << endl;
};


  




SweepInfo::Data::Data () :
		accumulated_error (0),
		cells (0),
		primal_dofs (0),
		dual_dofs (0)
{};







// explicit instantiations
template 
void SweepInfo::write_summary (const list<EvaluationBase<2>*> &eval_list,
			       ostream &out) const;

/* $Id$ */


#include <base/quadrature.h>
#include <base/function.h>
#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/vector_memory.h>
#include <lac/precondition.h>
#include <grid/geometry_info.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <fe/fe_update_flags.h>
#include <numerics/matrices.h>
#include <numerics/dof_renumbering.h>


#include <fstream>
#include <iomanip>




static const pair<unsigned int, double> relaxations[3]
= { make_pair(100,5), make_pair(300,3), make_pair(500,2) };


static const TimeStepBase_Tria<2>::RefinementFlags::CorrectionRelaxations
wave_correction_relaxations (1,
			     vector<pair<unsigned int,double> > (&relaxations[0],
								 &relaxations[3]));



template <int dim>
TimeStepBase_Wave<dim>::TimeStepBase_Wave ():
		TimeStepBase_Tria<dim> (),
		parameters (*static_cast<WaveParameters<dim>*>(0))
{};



template <int dim>
TimeStepBase_Wave<dim>::TimeStepBase_Wave (const double                    time,
					   TimeStepBase_Tria<dim>::Flags   flags,
					   const WaveParameters<dim>      &parameters)
		:
		TimeStepBase_Tria<dim> (time,
					*parameters.coarse_grid,
					flags,
					typename TimeStepBase_Wave<dim>::RefinementFlags
					(parameters.maximum_refinement,
					 1,
					 0,
					 parameters.cell_number_corridor.first,
					 parameters.cell_number_corridor.first,
					 wave_correction_relaxations,
					 parameters.cell_number_correction_steps,
					 (parameters.refinement_strategy ==
					  WaveParameters<dim>::dual_estimator),
					 true)),
		parameters (parameters)
{};



template <int dim>
const TimeStep_Primal<dim> &
TimeStepBase_Wave<dim>::get_timestep_primal () const
{
  return dynamic_cast<const TimeStep_Primal<dim> &> (*this);
};



template <int dim>
const TimeStep_Dual<dim> &
TimeStepBase_Wave<dim>::get_timestep_dual () const
{
  return dynamic_cast<const TimeStep_Dual<dim> &> (*this);
};



template <int dim>
const TimeStep_Postprocess<dim> &
TimeStepBase_Wave<dim>::get_timestep_postprocess () const
{
  return dynamic_cast<const TimeStep_Postprocess<dim> &> (*this);
};



template <int dim>
string TimeStepBase_Wave<dim>::tmp_filename_base (const string &branch_signature) const
{
  return (parameters.tmp_directory +
	  branch_signature + 's' +
	  int_to_string (sweep_no, 2) + 't' +
	  int_to_string (timestep_no, 4));
};



template <int dim>
void TimeStepBase_Wave<dim>::attach_sweep_info (SweepInfo &si)
{
  sweep_info = &si;
};



template <int dim>
void TimeStepBase_Wave<dim>::attach_sweep_data (SweepData<dim> &sd)
{
  sweep_data = &sd;
};






/* --------------------------------------------------------------*/


template <int dim>
TimeStep_Wave<dim>::TimeStep_Wave (const string fe_name) :
		dof_handler (0),
		fe (FEHelper<dim>::get_fe(fe_name)),
		quadrature (FEHelper<dim>::get_quadrature(fe_name)),
		quadrature_face (FEHelper<dim>::get_quadrature_face(fe_name)),
		statistic_data()
{};



template <int dim>
TimeStep_Wave<dim>::~TimeStep_Wave ()
{
  Assert (dof_handler == 0, ExcInternalError());
  Assert (constraints.n_constraints() == 0, ExcInternalError());
  Assert (system_sparsity.empty(), ExcInternalError());
  Assert (mass_matrix.empty(), ExcInternalError());
  Assert (laplace_matrix.empty(), ExcInternalError());
  Assert (u.size() ==0, ExcInternalError());
  Assert (v.size() ==0, ExcInternalError());
};

				   

template <int dim>
void TimeStep_Wave<dim>::wake_up (const unsigned int wakeup_level) 
{
				   // only do something if we are
				   // right at the beginning of a
				   // time level
  if (wakeup_level==0)
    {
				       // first make the dof handler
      Assert (dof_handler==0, ExcInternalError());

      sweep_info->get_timers().grid_generation.start();

      dof_handler = new DoFHandler<dim>(*tria);
      dof_handler->distribute_dofs (fe);

      if (parameters.renumber_dofs)
	DoFRenumbering::Cuthill_McKee (*dof_handler);
      

      constraints.clear ();
      DoFTools::make_hanging_node_constraints (*dof_handler, constraints);
      constraints.close ();

      sweep_info->get_timers().grid_generation.stop();
      
      Assert (u.size()==0, ExcInternalError ());
      Assert (v.size()==0, ExcInternalError ());

      switch (next_action)
	{
	  case primal_problem:
	  case dual_problem:
	  {
					     // assert that this function only
					     // wakes up data members in the right
					     // branch of the multiple inheritance
					     // lattice, i.e. the dual problem
					     // branch may only be woken up if the
					     // dual problem is solved and vica
					     // versa
	    Assert (((next_action == primal_problem) &&
		     (static_cast<const TimeStep_Wave<dim>*>(&get_timestep_primal())
		      == this))
		    ||
		    ((next_action == dual_problem) &&
		     (static_cast<const TimeStep_Wave<dim>*>(&get_timestep_dual())
		      == this)),
		    ExcInternalError());
	    
					     // if we are to extrapolate the old
					     // solutions, we overwrite the previous
					     // content of the vectors anyway, so
					     // we can use the fast initialization
	    u.reinit (dof_handler->n_dofs(),
		      parameters.extrapolate_old_solutions && (timestep_no!=0));
	    v.reinit (dof_handler->n_dofs(),
		      parameters.extrapolate_old_solutions && (timestep_no!=0));
	    break;
	  };
	   
	  case postprocess:
	  {
	    sweep_info->get_timers().postprocessing.start();
					     // reload data vectors from disk
	    ifstream tmp_in(tmp_filename_base(branch_signature()).c_str());
	    u.block_read (tmp_in);
	    v.block_read (tmp_in);
	    tmp_in.close ();

	    sweep_info->get_timers().postprocessing.stop();
		    
	    break;
	  };
	   
	  default:
		Assert (false, ExcInternalError());
	};
    };
};



template <int dim>
void TimeStep_Wave<dim>::sleep (const unsigned int sleep_level) 
{
  switch (sleep_level)
    {
      case 1:
      {
	Assert (dof_handler!=0, ExcInternalError());
      
	delete dof_handler;
	dof_handler = 0;

	Assert (u.size() != 0, ExcInternalError());
	Assert (v.size() != 0, ExcInternalError());

	ofstream tmp_out(tmp_filename_base(branch_signature()).c_str());
	u.block_write (tmp_out);
	v.block_write (tmp_out);
	tmp_out.close ();
	
	u.reinit (0);
	v.reinit (0);
	
	Assert (constraints.n_constraints() == 0, ExcInternalError());
	Assert (system_sparsity.empty(), ExcInternalError());
	Assert (mass_matrix.empty(), ExcInternalError());
	Assert (laplace_matrix.empty(), ExcInternalError());

	break;
      };

      case 0:
      {
					 // these are the data we don't need
					 // any more right after the time step
					 // do this action for the derived classes
	constraints.clear ();
	system_sparsity.reinit (0,0,0);
	mass_matrix.reinit (system_sparsity);
	laplace_matrix.reinit (system_sparsity);

	break;
      };

      default:
	    Assert (false, ExcInternalError());
    };
};



template <int dim>
void TimeStep_Wave<dim>::end_sweep ()
{
  string tmp_filename = tmp_filename_base(branch_signature());
  remove (tmp_filename.c_str());
};


    
template <int dim>
unsigned int TimeStep_Wave<dim>::solve (const UserMatrix       &matrix,
					Vector<double>         &solution,
					const Vector<double>   &rhs) const {
  SolverControl            control(2000, 1.e-12);
  PrimitiveVectorMemory<>  memory;
  SolverCG<UserMatrix>     pcg(control,memory);

				   // solve
  pcg.solve (matrix, solution, rhs,
	     PreconditionUseMatrix<UserMatrix>
	     (matrix,
	      &UserMatrix::precondition));
				   // distribute solution
  constraints.distribute (solution);

  return control.last_step();
};



template <int dim>
void TimeStep_Wave<dim>::create_matrices () 
{        
				   // reinitialize sparsity and vector size
  system_sparsity.reinit (dof_handler->n_dofs(), dof_handler->n_dofs(),
			  dof_handler->max_couplings_between_dofs());
				   // build sparsity pattern and condense
				   // with hanging nodes
  DoFTools::make_sparsity_pattern (*dof_handler, system_sparsity);
  constraints.condense (system_sparsity);
  system_sparsity.compress ();
      
				       // reinit matrices
  laplace_matrix.reinit (system_sparsity);
  mass_matrix.reinit (system_sparsity);

				   // now actually assemble the matrices
  const unsigned int dofs_per_cell       = fe.dofs_per_cell,
		     n_q_points       = quadrature.n_quadrature_points;

  const bool   density_constant = parameters.density_constant,
	     stiffness_constant = parameters.stiffness_constant;

  vector<double> density_values   (n_q_points, 1.);
  vector<double> stiffness_values (n_q_points, 1.);

				   // if a coefficient is constant, get
				   // its value
  if (density_constant)
    fill_n (density_values.begin(), n_q_points,
	    parameters.density->value(Point<dim>()));
  if (stiffness_constant)
    fill_n (stiffness_values.begin(), n_q_points,
	    parameters.stiffness->value(Point<dim>()));
  
  
  FEValues<dim>  fe_values (fe, quadrature,
			    UpdateFlags(update_values |
					update_gradients  |
					update_JxW_values |
					(!density_constant || !stiffness_constant ?
					 update_q_points :
					 0)));

				   // indices of all the dofs on this
				   // cell
  vector<int>    dof_indices_on_cell (dofs_per_cell);
  FullMatrix<double> cell_mass_matrix (dofs_per_cell, dofs_per_cell);
  FullMatrix<double> cell_laplace_matrix (dofs_per_cell, dofs_per_cell);

  
  for (typename DoFHandler<dim>::active_cell_iterator cell=dof_handler->begin_active();
       cell != dof_handler->end(); ++cell)
    {
      fe_values.reinit (cell);
      cell_mass_matrix.clear ();
      cell_laplace_matrix.clear ();
      cell->get_dof_indices (dof_indices_on_cell);

      const FullMatrix<double>              &shape_values = fe_values.get_shape_values ();
      const vector<vector<Tensor<1,dim> > > &shape_grads  = fe_values.get_shape_grads ();
      const vector<double>                  &JxW_values   = fe_values.get_JxW_values ();

				       // if necessary: get the values of any
				       // of the coefficients at the quadrature
				       // points
      if (!density_constant || !stiffness_constant)
	{
	  const vector<Point<dim> > &quadrature_points = fe_values.get_quadrature_points ();
	  if (!density_constant)
	    parameters.density->value_list (quadrature_points,
					    density_values);
	  if (!stiffness_constant)
	    parameters.stiffness->value_list (quadrature_points,
					      stiffness_values);
	};
      
				       // now do the loop
      for (unsigned int q_point=0; q_point<fe_values.n_quadrature_points; ++q_point)
	for (unsigned int i=0; i<dofs_per_cell; ++i)
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    {
	      cell_mass_matrix(i,j) += (shape_values(i, q_point) *
					shape_values(j, q_point) *
					JxW_values[q_point]      *
					density_values[q_point]);
	      cell_laplace_matrix(i,j) += (shape_grads[i][q_point] *
					   shape_grads[j][q_point] *
					   JxW_values[q_point]      *
					   stiffness_values[q_point]);
	    };

				       // now transfer to global matrices
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	for (unsigned int j=0; j<dofs_per_cell; ++j)
	  {
	    mass_matrix.add(dof_indices_on_cell[i],
			    dof_indices_on_cell[j],
			    cell_mass_matrix(i,j));
	    laplace_matrix.add(dof_indices_on_cell[i],
			       dof_indices_on_cell[j],
			       cell_laplace_matrix(i,j));
	  };
    };
};



template <int dim>
void TimeStep_Wave<dim>::transfer_old_solutions (Vector<double> &old_u,
						 Vector<double> &old_v) const 
{
  const DoFHandler<dim> *present_dof_handler = dof_handler,
			*    old_dof_handler = 0;
  const Vector<double>  *old_grid_u = 0,
			*old_grid_v = 0;
  
  switch (next_action)
    {
      case primal_problem:
	    Assert (previous_timestep != 0, ExcInternalError());
	    
	    old_dof_handler = (static_cast<const TimeStepBase_Wave<dim>*>
			       (previous_timestep)->get_timestep_primal()).dof_handler;
	    old_grid_u      = &(static_cast<const TimeStepBase_Wave<dim>*>
				(previous_timestep)->get_timestep_primal()).u;
	    old_grid_v      = &(static_cast<const TimeStepBase_Wave<dim>*>
				(previous_timestep)->get_timestep_primal()).v;
	    
	    break;

      case dual_problem:
	    Assert (next_timestep != 0, ExcInternalError());
	    
	    old_dof_handler = (static_cast<const TimeStepBase_Wave<dim>*>
			       (next_timestep)->get_timestep_dual()).dof_handler;
	    old_grid_u      = &(static_cast<const TimeStepBase_Wave<dim>*>
				(next_timestep)->get_timestep_dual()).u;
	    old_grid_v      = &(static_cast<const TimeStepBase_Wave<dim>*>
				(next_timestep)->get_timestep_dual()).v;

	    break;
    };
  
  Assert (old_dof_handler != 0, ExcInternalError());

  DoFHandler<dim>::cell_iterator old_cell = old_dof_handler->begin(),
				 new_cell = present_dof_handler->begin();
  for (; old_cell != (old_dof_handler->get_tria().n_levels() == 1  ?
		      static_cast<DoFHandler<dim>::cell_iterator>(old_dof_handler->end()) :
		      old_dof_handler->begin(1));
       ++old_cell, new_cell)
    transfer_old_solutions (old_cell, new_cell,
			    *old_grid_u, *old_grid_v,
			    old_u, old_v);
};



template <int dim>
void
TimeStep_Wave<dim>::transfer_old_solutions (const typename DoFHandler<dim>::cell_iterator &old_cell,
					    const typename DoFHandler<dim>::cell_iterator &new_cell,
					    const Vector<double>  &old_grid_u,
					    const Vector<double>  &old_grid_v,
					    Vector<double>        &old_u,
					    Vector<double>        &old_v) const 
{
  if (!old_cell->has_children() && !new_cell->has_children()) 
    {
				       // none of the children are active, so
				       // recurse into the triangulation
      for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	transfer_old_solutions (old_cell->child(c),
				new_cell->child(c),
				old_grid_u, old_grid_v,
				old_u, old_v);
    }
  else
				     // one of the cells is active
    {
				       // get values from
				       // old cell and set on the new one
      Vector<double> cell_data (fe.dofs_per_cell);

      old_cell->get_interpolated_dof_values (old_grid_u, cell_data);
      new_cell->set_dof_values_by_interpolation (cell_data, old_u);
      
      old_cell->get_interpolated_dof_values (old_grid_v, cell_data);
      new_cell->set_dof_values_by_interpolation (cell_data, old_v);
    };
};



template <int dim>
pair<double,double>
TimeStep_Wave<dim>::compute_energy () {
  pair<double,double> energy;
  
  switch (next_action)
    {
      case primal_problem:
	    energy.first = 0.5*laplace_matrix.matrix_norm (u);
	    energy.second = 0.5*mass_matrix.matrix_norm(v);
	    break;

      case dual_problem:
	    energy.first = 0.5*laplace_matrix.matrix_norm (v);
	    energy.second = 0.5*mass_matrix.matrix_norm(u);
	    break;

      default:
	    Assert (false, ExcInternalError());
    };

  return energy;
};



template <int dim>
TimeStep_Wave<dim>::StatisticData::
StatisticData () :
		n_active_cells (0),
		n_dofs (0),
		n_solver_steps_helmholtz (0),
		n_solver_steps_projection (0),
		energy (make_pair(0.0, 0.0))
{};



template <int dim>
TimeStep_Wave<dim>::StatisticData::
StatisticData (const unsigned int        n_active_cells,
	       const unsigned int        n_dofs,
	       const unsigned int        n_solver_steps_helmholtz,
	       const unsigned int        n_solver_steps_projection,
	       const pair<double,double> energy) :
		n_active_cells (n_active_cells),
		n_dofs (n_dofs),
		n_solver_steps_helmholtz (n_solver_steps_helmholtz),
		n_solver_steps_projection (n_solver_steps_projection),
		energy (energy)
{};



template <int dim>
void
TimeStep_Wave<dim>::StatisticData::write_descriptions (ostream &out) 
{
  out << "#    number of active cells"                 << endl
      << "#    number of degrees of freedom"           << endl
      << "#    iterations for the helmholtz equation"  << endl
      << "#    iterations for the projection equation" << endl
      << "#    elastic energy"                         << endl
      << "#    kinetic energy"                         << endl
      << "#    total energy"                           << endl;
};



template <int dim>
void TimeStep_Wave<dim>::StatisticData::write (ostream &out) const
{
  out << n_active_cells             << ' '
      << n_dofs                     << ' '
      << n_solver_steps_helmholtz   << ' '
      << n_solver_steps_projection  << ' '
      << energy.first               << ' '
      << energy.second              << ' '
      << energy.first+energy.second;
};







// explicit instantiations
template class TimeStepBase_Wave<2>;
template class TimeStep_Wave<2>;
/* $Id$ */

#include <base/function.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe_values.h>
#include <fe/fe.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>


#include <iomanip>



template <int dim>
TimeStep_Dual<dim>::TimeStep_Dual (const string &dual_fe)
		:
		TimeStep_Wave<dim> (dual_fe)
{};



template <int dim>
void TimeStep_Dual<dim>::do_initial_step () {
  deallog << "  Dual problem: time="
       << time
       << ", step=" << timestep_no
       << ", sweep=" << sweep_no
       << ". "
       << tria->n_active_cells() << " cells, "
       << dof_handler->n_dofs() << " dofs";
  
				   // add up sweep-accumulated data. count
				   // u and v as separate dofs
				   //
				   // do not add up cells, since this is already
				   // done in the primal problem
  sweep_info->get_data().dual_dofs += dof_handler->n_dofs() * 2;

  Vector<double> tmp_u_bar, tmp_v_bar;

				   // get evaluation of dual functional
				   // at end time
  parameters.dual_functional->reset (*this);
  parameters.dual_functional->
    compute_endtime_vectors (tmp_u_bar, tmp_v_bar);
				   // compute final values for the dual
				   // problem by projection, i.e. by
				   // inversion of the mass matrix; don't
				   // do so if the solution will be zero
				   // (inversion would not take long, but
				   // assembling the matrices is expensive)
  u.reinit (tmp_u_bar.size());
  v.reinit (tmp_v_bar.size());
  if ((tmp_u_bar.linfty_norm() > 0) || (tmp_v_bar.linfty_norm() > 0))
    {
      UserMatrix system_matrix (system_sparsity,
				parameters.preconditioning);
      system_matrix.copy_from (mass_matrix);
      constraints.condense (system_matrix);
      const unsigned int
	solver_steps1 = solve (system_matrix, u, tmp_u_bar),
	solver_steps2 = solve (system_matrix, v, tmp_v_bar);

      statistic_data = StatisticData (tria->n_active_cells(),
				      dof_handler->n_dofs(),
				      solver_steps1, solver_steps2,
				      compute_energy ());
    }
  else
    statistic_data = StatisticData (tria->n_active_cells(),
				    dof_handler->n_dofs(),
				    0, 0,
				    make_pair (0.0, 0.0));
  deallog << "." << endl;
};



template <int dim>
void TimeStep_Dual<dim>::do_timestep ()
{
  deallog << "  Dual problem: time="
       << time
       << ", step=" << timestep_no
       << ", sweep=" << sweep_no
       << ". "
       << tria->n_active_cells() << " cells, "
       << dof_handler->n_dofs() << " dofs";

				   // add up sweep-accumulated data. count
				   // u and v as separate dofs
				   //
				   // do not add up cells, since this is already
				   // done in the primal problem
  sweep_info->get_data().dual_dofs += dof_handler->n_dofs() * 2;

  const double time_step = get_forward_timestep ();

 				   // Vectors holding the right hand sides of
 				   // the two equations.
  Vector<double> right_hand_side1 (dof_handler->n_dofs());
  Vector<double> right_hand_side2 (dof_handler->n_dofs());
  
				   // Vector holding a the values for
				   // u and v of the previous time step.
				   // these are used in case we want to
				   // use extrapolation from the previous
				   // time step to the present one
  Vector<double> old_u, old_v;
  if (parameters.extrapolate_old_solutions)
    {
      old_u.reinit (dof_handler->n_dofs());
      old_v.reinit (dof_handler->n_dofs());

      transfer_old_solutions (old_u, old_v);
    };
    
  assemble_vectors (right_hand_side1, right_hand_side2);

  UserMatrix system_matrix (system_sparsity, parameters.preconditioning);
  system_matrix.copy_from (mass_matrix);
  system_matrix.add_scaled (time_step * time_step *
			    parameters.theta *
			    parameters.theta,
			    laplace_matrix);
  constraints.condense (system_matrix);
	
  if (parameters.extrapolate_old_solutions)
				     // solve with a hopefully good guess
				     // as start vector
    {
      v  = old_v;
      v.add (time_step, old_u);
    };
				   // in the other case, the wake_up
				   // function of the base class has set
				   // the solution vector's values to
				   // zero already.


				   // in 1d, do not set boundary conditions
				   // at all
				   //
				   // note: in boundary_value_map, all entries
				   // for dirichlet boundary nodes are set to
				   // zero. we re-use them later, and because
				   // zero is such a universal constant, we
				   // don't even need to recompute the values!
  map<int,double> boundary_value_list;
  if (dim != 1)
    {
      static const ZeroFunction<dim> boundary_values;
      
      VectorTools::interpolate_boundary_values (*dof_handler, 0, boundary_values,
						     boundary_value_list);
      MatrixTools<dim>::apply_boundary_values (boundary_value_list,
					       system_matrix, v,
					       right_hand_side1);
    };
  
  const unsigned int solver_steps1 = solve (system_matrix, v, right_hand_side1);
	
  system_matrix.copy_from (mass_matrix);
  constraints.condense (system_matrix);
  if (true)
    {
      Vector<double> tmp (right_hand_side2.size());
      laplace_matrix.vmult (tmp, v);
      right_hand_side2.add (-parameters.theta*time_step, tmp);
    };
  constraints.condense (right_hand_side2);
				   ///////////////////////////
				   // This is not ok here, for two reasons:
				   // 1. it assumes that for v the same
				   //    bc hold as for u; build the list
				   //    of bc for v separately, this way
				   //    it only holds for u=v=0
				   // 2. v has no boundary conditions at
				   //    all!
				   ///////////////////////////
  if (dim != 1)
				     // note: the values in boundary_value_map
				     // are already set for the first component
				     // and have not been touched since.
    MatrixTools<dim>::apply_boundary_values (boundary_value_list,
					     system_matrix, u,
					     right_hand_side2);
  
  if (parameters.extrapolate_old_solutions)
				     // solve with a hopefully good guess
				     // as start vector
    {
      u  = v;
      u -= old_v;
      u.scale (2./time_step);
      u -= old_u;
    };
  
  const unsigned int solver_steps2 = solve (system_matrix, u, right_hand_side2);

  statistic_data = StatisticData (tria->n_active_cells(),
				  dof_handler->n_dofs(),
				  solver_steps1,
				  solver_steps2,
				  compute_energy ());
  
  deallog << "." << endl;
};




template <int dim>
void TimeStep_Dual<dim>::solve_dual_problem ()
{
  sweep_info->get_timers().dual_problem.start();
  if (next_timestep == 0)
    do_initial_step ();
  else
    do_timestep ();
  sweep_info->get_timers().dual_problem.stop();
};



template <int dim>
string TimeStep_Dual<dim>::branch_signature () const 
{
  return "d";
};



template <int dim>
void TimeStep_Dual<dim>::wake_up (const unsigned int wakeup_level)
{
  TimeStep_Wave<dim>::wake_up (wakeup_level);
  
  sweep_info->get_timers().dual_problem.start();
  if ((wakeup_level==0) && (next_action==dual_problem))
    {
      Assert (system_sparsity.empty(), ExcInternalError());
      
      create_matrices ();
    };
  sweep_info->get_timers().dual_problem.stop();
};



template <int dim>
void TimeStep_Dual<dim>::assemble_vectors (Vector<double> &right_hand_side1,
					   Vector<double> &right_hand_side2) {
				   // don't do some things for the initial
				   // step since we don't need them there
  Assert (next_timestep != 0, ExcInternalError());
  
				   // construct right hand side
  build_rhs (right_hand_side1, right_hand_side2);

				   // compute contributions of error
				   // functional to right hand sides
  Vector<double> dual1, dual2;
  parameters.dual_functional->reset (*this);
  parameters.dual_functional->compute_functionals (dual1, dual2);

  const double timestep = get_forward_timestep();
  right_hand_side1.add (timestep, dual2);
  right_hand_side1.add (parameters.theta * timestep * timestep, dual1);

  right_hand_side2.add (timestep, dual1);

				   // condense right hand side in-place
  constraints.condense (right_hand_side1);
};



template <int dim>
void TimeStep_Dual<dim>::build_rhs (Vector<double> &right_hand_side1,
				    Vector<double> &right_hand_side2) {
				   // select the TimeStep_Wave part in the
				   // TimeStep_Primal branch
  const TimeStep_Dual<dim> &previous_time_level
    = static_cast<const TimeStepBase_Wave<dim>*>(next_timestep)->get_timestep_dual();

  Assert (previous_time_level.tria->n_cells(0) == tria->n_cells(0),
	  ExcCoarsestGridsDiffer());

				   // convenience typedef
  typedef DoFHandler<dim>::cell_iterator cell_iterator;

				   // create this here and pass it to
				   // the cellwise function since it
				   // is expensive to create it for
				   // every cell
  FEValues<dim> fe_values (fe, quadrature,
			   UpdateFlags(update_values |
				       update_gradients |
				       update_JxW_values |
				       update_q_points));

  
  cell_iterator old_cell = previous_time_level.dof_handler->begin(),
		new_cell = dof_handler->begin(),
		end_cell = (tria->n_levels() == 1                  ?
			    static_cast<cell_iterator>(dof_handler->end()) :
			    dof_handler->begin(1));  
  for (; new_cell!=end_cell; ++new_cell, ++old_cell)
    build_rhs (old_cell, new_cell,
	       fe_values,
	       right_hand_side1, right_hand_side2);
};



template <int dim>
void
TimeStep_Dual<dim>::build_rhs (const DoFHandler<dim>::cell_iterator &old_cell,
			       const DoFHandler<dim>::cell_iterator &new_cell,
			       FEValues<dim>        &fe_values,
			       Vector<double>       &right_hand_side1,
			       Vector<double>       &right_hand_side2) {
				   // declare this type for convenience
  typedef DoFHandler<dim>::cell_iterator cell_iterator;

				   // both cells have children, so
				   // recurse into the tree
  if (old_cell->has_children() && new_cell->has_children()) 
    {
      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
	build_rhs (old_cell->child(child),
		   new_cell->child(child),
		   fe_values,
		   right_hand_side1,
		   right_hand_side2);
      return;
    };


				   // select the TimeStep_Wave part in the
				   // TimeStep_Dual branch
  const TimeStep_Dual<dim> &previous_time_level
    = static_cast<const TimeStepBase_Wave<dim>*>(next_timestep)->get_timestep_dual();

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const double time_step = get_forward_timestep();

				   // both cells are on the same refinement
				   // level
  if (!old_cell->has_children() && !new_cell->has_children()) 
    {
      fe_values.reinit (old_cell);
      const FullMatrix<double>             &values    = fe_values.get_shape_values ();
      const vector<vector<Tensor<1,dim> > >&gradients = fe_values.get_shape_grads ();
      const vector<double>                 &weights   = fe_values.get_JxW_values ();

      FullMatrix<double>    cell_matrix (dofs_per_cell, dofs_per_cell);

      vector<double> density_values(fe_values.n_quadrature_points);
      parameters.density->value_list (fe_values.get_quadrature_points(),
				      density_values);
      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	for (unsigned int i=0; i<dofs_per_cell; ++i) 
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    cell_matrix(i,j) += (values(i,point) *
				 values(j,point)) *
				weights[point] *
				density_values[point];

      Vector<double> tmp (dofs_per_cell);
				       // this is the right hand side of the
				       // first equation
				       // for the theta scheme:
				       //    rhs1 := Mv^1 + kMu^1
				       //           -(1-theta)theta k^2 Av^1
      Vector<double> rhs1 (dofs_per_cell);

				       // this is the part of the right hand side
				       // of the second equation which depends
				       // on the solutions of the previous time
				       // step.
				       // for the theta scheme:
				       //    rhs2 := Mu^1-(1-theta)kAv^1
      Vector<double> rhs2 (dofs_per_cell);
	    
      				       // vector of values of the function on the
				       // old grid restricted to one cell
      Vector<double> old_dof_values_v (dofs_per_cell);
				       // vector of old u and v times the local
				       // mass matrix
      Vector<double> local_M_u (dofs_per_cell);
      Vector<double> local_M_v (dofs_per_cell);
      Vector<double> local_A_v (dofs_per_cell);      
				       // transfer v+k*u. Note that we need
				       // old_dof_values_u again below
      old_cell->get_dof_values (previous_time_level.v, old_dof_values_v);
      cell_matrix.vmult (local_M_v, old_dof_values_v);
      
      old_cell->get_dof_values (previous_time_level.u, tmp);
      cell_matrix.vmult (local_M_u, tmp);

				       // now for the part with the laplace
				       // matrix
      cell_matrix.clear ();
      vector<double> stiffness_values(fe_values.n_quadrature_points);
      parameters.stiffness->value_list (fe_values.get_quadrature_points(),
					stiffness_values);
      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	for (unsigned int i=0; i<dofs_per_cell; ++i) 
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    cell_matrix(i,j) += (gradients[i][point] *
				 gradients[j][point]) *
				weights[point] *
				stiffness_values[point];
      cell_matrix.vmult (local_A_v, old_dof_values_v);

      rhs1 = local_M_v;
      rhs1.add (time_step, local_M_u);
      rhs1.add ((-time_step*time_step*
		 parameters.theta*
		 (1-parameters.theta)),
		local_A_v);
      rhs2 = local_M_u;
      rhs2.add (-(1-parameters.theta)*
		time_step,
		local_A_v);

				       // transfer into the global
				       // right hand side
      vector<int> new_dof_indices (dofs_per_cell, -1);
      new_cell->get_dof_indices (new_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  right_hand_side1(new_dof_indices[i]) += rhs1(i);
	  right_hand_side2(new_dof_indices[i]) += rhs2(i);
	};
      
      return;
    };

				   // only old cell is refined
  if (old_cell->has_children() && !new_cell->has_children())
    {
				       // this is the right hand side of the
				       // first equation
				       // for the theta scheme:
				       //    rhs1 := Mv^0 + kMu^1
				       //           -(1-theta)theta k^2 Av^1
      Vector<double> rhs1 (dofs_per_cell);
				       // this is the part of the right hand side
				       // of the second equation which depends
				       // on the solutions of the previous time
				       // step.
				       // for the theta scheme:
				       //    rhs2 := Mu^1-(1-theta)kAv^1
      Vector<double> rhs2 (dofs_per_cell);

				       // collect the contributions of the
				       // child cells (and possibly their
				       // children as well)
      collect_from_children (old_cell, fe_values, rhs1, rhs2);
      
				       // transfer into the global
				       // right hand side
      vector<int> new_dof_indices (dofs_per_cell);
      new_cell->get_dof_indices (new_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i) 
	{
	  right_hand_side1(new_dof_indices[i]) += rhs1(i);
	  right_hand_side2(new_dof_indices[i]) += rhs2(i);
	};

      return;
    };

  				   // only new cell is refined
  if (!old_cell->has_children() && new_cell->has_children())
    {
				       // vector of values of the function
				       // on the old grid restricted to
				       // the large (old) cell
      Vector<double>  old_dof_values_u (dofs_per_cell);
      Vector<double>  old_dof_values_v (dofs_per_cell);
      old_cell->get_dof_values (previous_time_level.u, old_dof_values_u);
      old_cell->get_dof_values (previous_time_level.v, old_dof_values_v);

				       // distribute the contribution of the
				       // large old cell to the children on
				       // the new cell
      distribute_to_children (new_cell, fe_values,
			      old_dof_values_u, old_dof_values_v,
 			      right_hand_side1, right_hand_side2);

      return;
    };

  Assert (false, ExcInternalError());
};



template <int dim>
unsigned int
TimeStep_Dual<dim>::collect_from_children (const DoFHandler<dim>::cell_iterator &old_cell,
					   FEValues<dim>  &fe_values,
					   Vector<double> &rhs1,
					   Vector<double> &rhs2) const {
				   // maximal difference of levels between the
				   // cell to which we write and the cells from
				   // which we read. Default is one, but this is
				   // increased with each level of recursion
  unsigned int level_difference = 1;  
  
				   // select the TimeStep_Wave part in the
				   // TimeStep_Primal branch
  const TimeStep_Dual<dim> &previous_time_level
    = static_cast<const TimeStepBase_Wave<dim>*>(next_timestep)->get_timestep_dual();

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const double time_step = get_forward_timestep();

  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);

				       // these will hold the values of the
				       // solution on the old grid, i.e. on
				       // the small cells
  Vector<double>  local_old_dof_values_u (dofs_per_cell);
  Vector<double>  local_old_dof_values_v (dofs_per_cell);

				       // same for the contributions to the
				       // right hand sides of the projection
  Vector<double>  local_M_u (dofs_per_cell);
  Vector<double>  local_M_v (dofs_per_cell);
  Vector<double>  local_A_v (dofs_per_cell);
      
				   // this is the right hand side of the
				   // first equation
				   // for the theta scheme:
				   //    rhs1 := Mv^0 + kMu^1
				   //           -(1-theta)theta k^2 Av^1
  Vector<double> child_rhs1 (dofs_per_cell);
				   // this is the part of the right hand side
				   // of the second equation which depends
				   // on the solutions of the previous time
				   // step.
				   // for the theta scheme:
				   //    rhs2 := Mu^1-(1-theta)kAv^1
  Vector<double> child_rhs2 (dofs_per_cell);
      
  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c) 
    {
      const DoFHandler<dim>::cell_iterator old_child = old_cell->child(c);

      child_rhs1.clear ();
      child_rhs2.clear ();
      
				       // if this child is further subdivided:
				       // collect the contributions of the
				       // children
      if (old_child->has_children())
	{
	  const unsigned int l = collect_from_children (old_child, fe_values,
							child_rhs1, child_rhs2);
	  level_difference = max (l+1, level_difference);
	}
      else
	{
	  fe_values.reinit (old_child);
	  const FullMatrix<double>             &values    = fe_values.get_shape_values();
	  const vector<vector<Tensor<1,dim> > >&gradients = fe_values.get_shape_grads ();
	  const vector<double>                 &weights   = fe_values.get_JxW_values ();

					   // get solutions restricted to small
					   // cell
	  old_child->get_dof_values (previous_time_level.u, local_old_dof_values_u);
	  old_child->get_dof_values (previous_time_level.v, local_old_dof_values_v);

					   // compute M*(v+ku) on the small cell
	  cell_matrix.clear ();
	  vector<double> density_values(fe_values.n_quadrature_points);
	  parameters.density->value_list (fe_values.get_quadrature_points(),
					  density_values);
	  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		cell_matrix(i,j) += (values(i,point) *
				     values(j,point)) *
				    weights[point] *
				    density_values[point];

	  cell_matrix.vmult (local_M_u, local_old_dof_values_u);
	  cell_matrix.vmult (local_M_v, local_old_dof_values_v);

					   // now for the part with the laplace
					   // matrix
	  cell_matrix.clear ();
	  vector<double> stiffness_values(fe_values.n_quadrature_points);
	  parameters.stiffness->value_list (fe_values.get_quadrature_points(),
						stiffness_values);
	  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		cell_matrix(i,j) += (gradients[i][point] *
				     gradients[j][point]) *
				    weights[point] *
				    stiffness_values[point];
	  cell_matrix.vmult (local_A_v, local_old_dof_values_v);
	  
	  child_rhs1 = local_M_v;
	  child_rhs1.add (time_step, local_M_u);
	  child_rhs1.add ((-time_step*time_step*
			   parameters.theta*
			   (1-parameters.theta)),
			  local_A_v);
	  child_rhs2 = local_M_u;
	  child_rhs2.add (-(1-parameters.theta)*
			  time_step,
			  local_A_v);
	};
      
				       // transfer the contribution of this
				       // child cell to its parent cell
				       // (#true# means: add up)
      fe.prolongate(c).Tvmult (rhs1, child_rhs1, true);
      fe.prolongate(c).Tvmult (rhs2, child_rhs2, true);
    };

  return level_difference;
};



template <int dim>
unsigned int
TimeStep_Dual<dim>::distribute_to_children (const DoFHandler<dim>::cell_iterator &new_cell,
					    FEValues<dim>         &fe_values,
					    const Vector<double>  &old_dof_values_u,
					    const Vector<double>  &old_dof_values_v,
					    Vector<double>        &right_hand_side1,
					    Vector<double>        &right_hand_side2) {
				   // maximal difference of levels between the
				   // cell to which we write and the cells from
				   // which we read. Default is one, but this is
				   // increased with each level of recursion
  unsigned int level_difference = 1;  
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const double time_step = get_forward_timestep();

  FullMatrix<double>    cell_matrix(dofs_per_cell, dofs_per_cell);
				   // set up a vector which will hold the
				   // restriction of the old
				   // functions (u,v) to a childcell
  Vector<double> local_old_dof_values_u (dofs_per_cell);
  Vector<double> local_old_dof_values_v (dofs_per_cell);

				   // vector of old u and v times the local
				   // mass matrix (on the small cells
				   // respectively)
  Vector<double> local_M_u (dofs_per_cell);
  Vector<double> local_M_v (dofs_per_cell);
  Vector<double> local_A_v (dofs_per_cell);

				   // this is the right hand side of the
				   // first equation
				   // for the theta scheme:
				   //    rhs1 := Mv^1 + kMu^1
				   //           -(1-theta)theta k^2 Av^1
  Vector<double> rhs1 (dofs_per_cell);

				   // this is the part of the right hand side
				   // of the second equation which depends
				   // on the solutions of the previous time
				   // step.
				   // for the theta scheme:
				   //    rhs2 := Mu^1-(1-theta)kAv^1
  Vector<double> rhs2 (dofs_per_cell);
	    
				   // indices of the dofs of a cell on
				   // the new grid
  vector<int> new_dof_indices (dofs_per_cell, -1);

      
				       // loop over the child cells
  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c) 
    {
      const DoFHandler<dim>::cell_iterator new_child = new_cell->child(c);

				       // get u and v on the childcells
      fe.prolongate(c).vmult (local_old_dof_values_u,
				   old_dof_values_u);
      fe.prolongate(c).vmult (local_old_dof_values_v,
				   old_dof_values_v);

      if (new_child->has_children())
					 // cell on new grid is further refined
					 // distribute data on this local cell
					 // to its children
	{
	  const unsigned int l = distribute_to_children (new_child, fe_values,
							 local_old_dof_values_u,
							 local_old_dof_values_v,
 							 right_hand_side1,
 							 right_hand_side2);
	  level_difference = max (l+1, level_difference);
	}
      else
					 // child is not further refined
					 // -> directly distribute data
	{
	  fe_values.reinit (new_child);
	  const FullMatrix<double>             &values    = fe_values.get_shape_values();
	  const vector<vector<Tensor<1,dim> > >&gradients = fe_values.get_shape_grads ();
	  const vector<double>                 &weights   = fe_values.get_JxW_values ();

					   // transfer v+ku
	  cell_matrix.clear ();
	  vector<double> density_values(fe_values.n_quadrature_points);
	  parameters.density->value_list (fe_values.get_quadrature_points(),
					  density_values);
	  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		cell_matrix(i,j) += (values(i,point) *
				     values(j,point)) *
				    weights[point] *
				    density_values[point];

	  cell_matrix.vmult (local_M_u, local_old_dof_values_u);
	  cell_matrix.vmult (local_M_v, local_old_dof_values_v);

					   // now for the part with the laplace
					   // matrix
	  cell_matrix.clear ();
	  vector<double> stiffness_values(fe_values.n_quadrature_points);
	  parameters.stiffness->value_list (fe_values.get_quadrature_points(),
					    stiffness_values);
	  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		cell_matrix(i,j) += (gradients[i][point] *
				     gradients[j][point]) *
				    weights[point] *
				    stiffness_values[point];
	  cell_matrix.vmult (local_A_v, local_old_dof_values_v);

	  rhs1 = local_M_v;
	  rhs1.add (time_step, local_M_u);
	  rhs1.add ((-time_step*time_step*
		     parameters.theta*
		     (1-parameters.theta)),
		    local_A_v);
	  rhs2 = local_M_u;
	  rhs2.add (-(1-parameters.theta)*
		    time_step,
		    local_A_v);
	  
					   // transfer into the global
					   // right hand side
	  new_child->get_dof_indices (new_dof_indices);
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      right_hand_side1(new_dof_indices[i]) += rhs1(i);
	      right_hand_side2(new_dof_indices[i]) += rhs2(i);
	    };
	};
    };

  return level_difference;
};




// explicit instantiations
template class TimeStep_Dual<2>;
/* $Id$ */


#include <base/tensor.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_constraints.h>
#include <fe/fe.h>
#include <fe/fe_values.h>
#include <numerics/error_estimator.h>


#include <fstream>
#include <iomanip>
#include <cmath>
#include <numeric>



template <int dim>
TimeStep_ErrorEstimation<dim>::TimeStep_ErrorEstimation () 
{};



template <int dim>
void TimeStep_ErrorEstimation<dim>::estimate_error ()
{
  sweep_info->get_timers().error_estimation.start();

  deallog << "[ee]";
  
  if ((parameters.refinement_strategy == WaveParameters<dim>::energy_estimator)
      ||
      (sweep_no < parameters.initial_energy_estimator_sweeps))
    estimate_error_energy (0);

  else
    {
				       // can't estimate error
				       // this way for the initial
				       // time level
      if (timestep_no != 0)
	estimate_error_dual ();
    };

  const double accumulated_error = accumulate (estimated_error_per_cell.begin(),
					       estimated_error_per_cell.end(),
					       0.0);
  statistic_data = StatisticData (accumulated_error);
  sweep_info->get_data().accumulated_error += accumulated_error;

  sweep_info->get_timers().error_estimation.stop();
};



template <int dim>
void TimeStep_ErrorEstimation<dim>::wake_up (const unsigned int wakeup_level)
{
  Assert (next_action==postprocess, ExcInternalError());

  if (wakeup_level==0)
    {
      Assert (estimated_error_per_cell.size()==0,
	      ExcInternalError());
      
      estimated_error_per_cell.reinit (tria->n_active_cells());
    };
};



template <int dim>
void TimeStep_ErrorEstimation<dim>::sleep (const unsigned int sleep_level)
{
  Assert (next_action==postprocess, ExcInternalError());

  if (sleep_level==0)
    {
      Assert (estimated_error_per_cell.size()!=0,
	      ExcInternalError());

      ofstream tmp_out(tmp_filename_base(branch_signature()).c_str());
      estimated_error_per_cell.block_write (tmp_out);
      tmp_out.close ();

      estimated_error_per_cell.reinit (0);
    };
};



template <int dim>
void
TimeStep_ErrorEstimation<dim>::get_tria_refinement_criteria (Vector<float> &indicators) const 
{
  get_error_indicators (indicators);
  for (Vector<float>::iterator i=indicators.begin(); i!=indicators.end(); ++i)
    *i = fabs(*i);
};


template <int dim>
void
TimeStep_ErrorEstimation<dim>::get_error_indicators (Vector<float> &indicators) const 
{
  ifstream in (tmp_filename_base(branch_signature()).c_str());
  indicators.block_read (in);
};



template <int dim>
void TimeStep_ErrorEstimation<dim>::estimate_error_energy (const unsigned int which_variables) {
  Assert (which_variables<=1, ExcInternalError());
  
  KellyErrorEstimator<dim>::FunctionMap neumann_boundary;
  static ZeroFunction<dim> homogeneous_neumann_bc;
  neumann_boundary[1] = &homogeneous_neumann_bc;

  const TimeStep_Wave<dim> &target = (which_variables==0 ?
				      static_cast<const TimeStep_Wave<dim>&>(get_timestep_primal()) :
				      static_cast<const TimeStep_Wave<dim>&>(get_timestep_dual ()));

  KellyErrorEstimator<dim>::estimate (*target.dof_handler,
				      target.quadrature_face,
				      neumann_boundary,
				      (which_variables==0 ?
				       target.u :
				       target.v),
				      estimated_error_per_cell,
				      vector<bool>(),
				      parameters.stiffness);

				   // if we are at the first time step, we
				   // try to adapt the mesh to the variable
				   // v also, since in some cases only v.neq.0
				   // and then the error indicator results in
				   // zero on all cells
  if (((previous_timestep == 0) && (which_variables==0)) ||
      ((next_timestep     == 0) && (which_variables==1)  ))
    {
      Vector<float> v_estimator(estimated_error_per_cell.size());
      KellyErrorEstimator<dim>::estimate (*target.dof_handler,
					  target.quadrature_face,
					  neumann_boundary,
					  (which_variables==0 ?
					   target.v :
					   target.u),
					  v_estimator,
					  vector<bool>(),
					  parameters.density);
      estimated_error_per_cell += v_estimator;
    };
};



template <int dim>
void TimeStep_ErrorEstimation<dim>::estimate_error_dual () {
  CellwiseError cellwise_error (tria->n_active_cells());

  const TimeStep_Primal<dim> &primal_problem     = get_timestep_primal(),
			     &primal_problem_old = static_cast<const TimeStepBase_Wave<dim>*>
						   (previous_timestep)->get_timestep_primal();
  const TimeStep_Dual<dim>   &dual_problem     = get_timestep_dual(),
			     &dual_problem_old = static_cast<const TimeStepBase_Wave<dim>*>
						 (previous_timestep)->get_timestep_dual();


				   // first clear the user pointers of
				   // the cells we need
  if (true)
    {
      DoFHandler<dim>::active_cell_iterator
	cell = primal_problem.dof_handler->begin_active();
      const DoFHandler<dim>::active_cell_iterator
	endc = primal_problem.dof_handler->end();
      for (; cell!=endc; ++cell)
	cell->clear_user_pointer();
    };
  
				   // set up some matrices used by the
				   // functions called in the sequel
  make_interpolation_matrices ();

				   // then go recursively through the two
				   // grids and collect the data
  if (true)
    {
      FEValues<dim> fe_values (dual_problem.fe,
			       dual_problem.quadrature,
			       UpdateFlags(update_values |
					   update_gradients |
					   update_second_derivatives |
					   update_JxW_values |
					   update_q_points));

				       // get dof iterators for the primal
				       // and dual dof handlers for the
				       // present and the last time level.
				       // since the coarse grids are the
				       // same and since we only loop
				       // over coarse grid cells here,
				       // the cells over which we loop
				       // match each other
      DoFHandler<dim>::cell_iterator
	primal_cell     = primal_problem.dof_handler->begin(),
	dual_cell       = dual_problem.dof_handler->begin(),
	primal_cell_old = primal_problem_old.dof_handler->begin(),
	dual_cell_old   = dual_problem_old.dof_handler->begin();
				       // get last cell to loop over. note that
				       // we only loop over the coarsest mesh
				       // in this function
      const DoFHandler<dim>::cell_iterator
	endc            = primal_problem.dof_handler->end(0);

				       // loop over all corse grid cells, since
				       // they are the same on the two time
				       // levels
      for (; primal_cell!=endc; (++primal_cell, ++dual_cell,
				 ++primal_cell_old, ++dual_cell_old))
	estimate_error_dual (primal_cell, dual_cell,
			     primal_cell_old, dual_cell_old,
			     cellwise_error,
			     fe_values);

      Assert (cellwise_error.next_free_slot == cellwise_error.errors.end(),
	      ::ExcInternalError());
    };

				   // compute the sum of the errors
				   // on the cells
  ErrorOnCell total_estimated_error;
  

				   // now fill the data we collected to the
				   // error_per_cell array
  Vector<float>::iterator i = estimated_error_per_cell.begin();
  DoFHandler<dim>::active_cell_iterator
    cell = primal_problem.dof_handler->begin_active();
  const DoFHandler<dim>::active_cell_iterator
    endc = primal_problem.dof_handler->end();
  for (; cell!=endc; ++cell, ++i)
    {
      const typename vector<ErrorOnCell>::iterator
	error_on_this_cell = static_cast<typename vector<ErrorOnCell>::iterator>(cell->user_pointer());
      Assert (error_on_this_cell != 0, ::ExcInternalError());

      cell->clear_user_pointer ();
      
      *i = error_on_this_cell->sum();
      total_estimated_error += *error_on_this_cell;
    };
};



template <int dim>
void
TimeStep_ErrorEstimation<dim>::estimate_error_dual (const DoFHandler<dim>::cell_iterator &primal_cell,
						    const DoFHandler<dim>::cell_iterator &dual_cell,
						    const DoFHandler<dim>::cell_iterator &primal_cell_old,
						    const DoFHandler<dim>::cell_iterator &dual_cell_old,
						    CellwiseError  &cellwise_error,
						    FEValues<dim>  &fe_values) const {
  
				   // if both of the two cells have children:
				   // recurse into the grid
  if (primal_cell->has_children() && primal_cell_old->has_children())
    {
      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
	estimate_error_dual (primal_cell->child(child),
			     dual_cell->child(child),
			     primal_cell_old->child(child),
			     dual_cell_old->child(child),
			     cellwise_error,
			     fe_values);
      return;
    };

  

  const TimeStep_Primal<dim> &primal_problem     = get_timestep_primal(),
			     &primal_problem_old = static_cast<const TimeStepBase_Wave<dim>*>
						   (previous_timestep)->get_timestep_primal();
  const TimeStep_Dual<dim>   &dual_problem     = get_timestep_dual(),
			     &dual_problem_old = static_cast<const TimeStepBase_Wave<dim>*>
						 (previous_timestep)->get_timestep_dual();

  const FiniteElement<dim> &primal_fe = get_timestep_primal().fe,
			   &dual_fe   = get_timestep_dual().fe;

  const unsigned int        dofs_per_cell_primal = primal_fe.dofs_per_cell,
			    dofs_per_cell_dual   = dual_fe.dofs_per_cell;  


				   // none of the two cells has children
  if (!primal_cell->has_children() && !primal_cell_old->has_children())
    {
				       // vector holding the solutions on
				       // this time level. u and v will
				       // hold the solution interpolated
				       // up to the ansatz degree of the
				       // dual problem.
      Vector<double> local_u(dofs_per_cell_dual), local_v(dofs_per_cell_dual);
      Vector<double> local_u_bar(dofs_per_cell_dual), local_v_bar(dofs_per_cell_dual);

				       // same thing for old solutions
      Vector<double> local_u_old(dofs_per_cell_dual), local_v_old(dofs_per_cell_dual);
      Vector<double> local_u_bar_old(dofs_per_cell_dual), local_v_bar_old(dofs_per_cell_dual);

				       // vectors to hold dof values on
				       // the primal/dual cell (temporary)
      Vector<double> primal_tmp(dofs_per_cell_primal);
      
				       // fill local solution vectors
      primal_cell->get_dof_values (primal_problem.u, primal_tmp);
      embedding_matrix.vmult (local_u, primal_tmp);

      primal_cell->get_dof_values (primal_problem.v, primal_tmp);
      embedding_matrix.vmult (local_v, primal_tmp);

      dual_cell->get_dof_values (dual_problem.u, local_u_bar);
      dual_cell->get_dof_values (dual_problem.v, local_v_bar);


      				       // fill local old solution vectors.
				       // no problems here, since the two
				       // cells are both unrefined
      primal_cell_old->get_dof_values (primal_problem_old.u, primal_tmp);
      embedding_matrix.vmult (local_u_old, primal_tmp);

      primal_cell_old->get_dof_values (primal_problem_old.v, primal_tmp);
      embedding_matrix.vmult (local_v_old, primal_tmp);

      dual_cell_old->get_dof_values (dual_problem_old.u, local_u_bar_old);
      dual_cell_old->get_dof_values (dual_problem_old.v, local_v_bar_old);

				       // store the error on this cell
      primal_cell->set_user_pointer (cellwise_error.next_free_slot);
      *cellwise_error.next_free_slot = error_formula (dual_cell,
						      local_u,     local_v,
						      local_u_bar, local_v_bar,
						      local_u_old,     local_v_old,
						      local_u_bar_old, local_v_bar_old,
						      fe_values);
      ++cellwise_error.next_free_slot;
      
      return;
    };



				   // only new cell has children. handle this
				   // case by prolonging the solutions on the
				   // old cell to its children and recursing
				   // thereon
  if (!primal_cell_old->has_children() && primal_cell->has_children())
    {
      Vector<double> local_u_old(dofs_per_cell_dual), local_v_old(dofs_per_cell_dual);
      Vector<double> local_u_bar_old(dofs_per_cell_dual), local_v_bar_old(dofs_per_cell_dual);

				       // vectors to hold dof values on
				       // the primal/dual cell (temporary)
      Vector<double> primal_tmp(dofs_per_cell_primal);
      
				       // fill local old solution vectors.
				       // no problems here, since the two
				       // cells are both unrefined
      primal_cell_old->get_dof_values (primal_problem_old.u, primal_tmp);
      embedding_matrix.vmult (local_u_old, primal_tmp);
  
      primal_cell_old->get_dof_values (primal_problem_old.v, primal_tmp);
      embedding_matrix.vmult (local_v_old, primal_tmp);
  
      dual_cell_old->get_dof_values (dual_problem_old.u, local_u_bar_old);
      dual_cell_old->get_dof_values (dual_problem_old.v, local_v_bar_old);


      compute_error_on_new_children (primal_cell, dual_cell,
				     local_u_old,
				     local_v_old,
				     local_u_bar_old,
				     local_v_bar_old,
				     cellwise_error,
				     fe_values);

      return;
    };


				   // last possibility: new cell is not
				   // refined, but old one is. in this case:
				   // collect error on this cell from the
				   // smaller ones on the old grid
				   //
				   // note that we have to perform the
				   // interpolation of the dual solution
				   // on the large cell of the new grid
				   // and have to pass the interpolant
				   // down to the children (which are
				   // taken from the old grid)
  if (primal_cell_old->has_children() && !primal_cell->has_children())
    {
				       // vector holding the solutions on
				       // this time level. u and v will
				       // hold the solution interpolated
				       // up to the ansatz degree of the
				       // dual problem.
      Vector<double> local_u(dofs_per_cell_dual), local_v(dofs_per_cell_dual);
      Vector<double> local_u_bar(dofs_per_cell_dual), local_v_bar(dofs_per_cell_dual);
      Vector<double> local_u_bar_old(dofs_per_cell_dual), local_v_bar_old(dofs_per_cell_dual);
      Vector<double> local_Ih_u_bar(dofs_per_cell_dual), local_Ih_v_bar(dofs_per_cell_dual);
      Vector<double> local_Ih_u_bar_old(dofs_per_cell_dual), local_Ih_v_bar_old(dofs_per_cell_dual);
      
				       // vectors to hold dof values on
				       // the primal/dual cell (temporary)
      Vector<double> primal_tmp(embedding_matrix.n());
      
				       // fill local solution vectors
      primal_cell->get_dof_values (primal_problem.u, primal_tmp);
      embedding_matrix.vmult (local_u, primal_tmp);
      
      primal_cell->get_dof_values (primal_problem.v, primal_tmp);
      embedding_matrix.vmult (local_v, primal_tmp);

				       // get the dual solution on the new
				       // time level to allow its interpolation
      dual_cell->get_dof_values (dual_problem.u, local_u_bar);
      dual_cell->get_dof_values (dual_problem.v, local_v_bar);

				       // now we have to get the interpolant
				       // of the dual solution on the old
				       // time level. Originally I wanted
				       // to do the following
				       //      dual_cell_old->get_dof_values
				       //	(previous_time_level->u_bar,
				       //	 local_u_bar_old
				       //	);
				       //      dual_cell_old->get_dof_values
				       //	(previous_time_level->v_bar,
				       //	 local_v_bar_old
				       //	);
				       //
				       // However, this must fail since
				       // dual_cell_old has children and
				       // we can't access data values on
				       // nonterminal cells...
				       //
				       // therefore, we use a new function
				       // which does exactly this interpolation
      dual_cell_old->get_interpolated_dof_values (dual_problem_old.u,
						  local_u_bar_old);
      dual_cell_old->get_interpolated_dof_values (dual_problem_old.v,
						  local_v_bar_old);

				       // compute the interpolant of w_bar and
				       // w_bar_old on the large cell
      interpolation_matrix.vmult (local_Ih_u_bar, local_u_bar);
      interpolation_matrix.vmult (local_Ih_v_bar, local_v_bar);
      interpolation_matrix.vmult (local_Ih_u_bar_old, local_u_bar_old);
      interpolation_matrix.vmult (local_Ih_v_bar_old, local_v_bar_old);

      primal_cell->set_user_pointer (cellwise_error.next_free_slot);
      *cellwise_error.next_free_slot
	= collect_error_from_children (primal_cell_old,
				       dual_cell_old,
				       local_u,            local_v,
				       local_u_bar,        local_v_bar,
				       local_Ih_u_bar,     local_Ih_v_bar,
				       local_Ih_u_bar_old, local_Ih_v_bar_old,
				       fe_values);
      ++cellwise_error.next_free_slot;
      
      return;
    };


  Assert (false, ExcInternalError());
};




template <int dim>
void TimeStep_ErrorEstimation<dim>::
compute_error_on_new_children (const DoFHandler<dim>::cell_iterator &primal_cell,
			       const DoFHandler<dim>::cell_iterator &dual_cell,
			       const Vector<double>  &local_u_old,
			       const Vector<double>  &local_v_old,
			       const Vector<double>  &local_u_bar_old,
			       const Vector<double>  &local_v_bar_old,
			       CellwiseError         &cellwise_error,
			       FEValues<dim>  &fe_values) const {
  const TimeStep_Primal<dim> &primal_problem = get_timestep_primal();
  const TimeStep_Dual<dim>   &dual_problem   = get_timestep_dual();

  const FiniteElement<dim> &dual_fe  = get_timestep_dual().fe;
  const unsigned int dofs_per_cell_dual = dual_fe.dofs_per_cell;
  
  
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    {
      				       // we have the solutions on the
				       // old (large) cell, we restrict it to
				       // each of the small cells
      Vector<double> child_u_old(dofs_per_cell_dual), child_v_old(dofs_per_cell_dual);
      Vector<double> child_u_bar_old(dofs_per_cell_dual), child_v_bar_old(dofs_per_cell_dual);

      dual_fe.prolongate(child).vmult (child_u_old, local_u_old);
      dual_fe.prolongate(child).vmult (child_v_old, local_v_old);
      dual_fe.prolongate(child).vmult (child_u_bar_old, local_u_bar_old);
      dual_fe.prolongate(child).vmult (child_v_bar_old, local_v_bar_old);

      const DoFHandler<dim>::cell_iterator
	new_primal_child = primal_cell->child(child),
	new_dual_child   = dual_cell->child(child);

      if (new_primal_child->has_children())
					 // cell on new grid is further refined
					 // distribute data on this local cell
					 // to its children
	compute_error_on_new_children (new_primal_child, new_dual_child,
				       child_u_old,
				       child_v_old,
				       child_u_bar_old,
				       child_v_bar_old,
				       cellwise_error,
				       fe_values);
      else
					 // we have reached the final level
					 // -> gather the information from
					 // the new cell and compute the
					 // error
	{
					   // vector holding the solutions on
					   // this time level. u and v will
					   // hold the solution interpolated
					   // up to the ansatz degree of the
					   // dual problem.
	  Vector<double> local_u(dofs_per_cell_dual), local_v(dofs_per_cell_dual);
	  Vector<double> local_u_bar(dofs_per_cell_dual), local_v_bar(dofs_per_cell_dual);
	  
					   // vectors to hold dof values on
					   // the primal/dual cell (temporary)
	  Vector<double> primal_tmp(embedding_matrix.n());
	  
					   // fill local solution vectors
	  new_primal_child->get_dof_values (primal_problem.u, primal_tmp);
	  embedding_matrix.vmult (local_u, primal_tmp);
	  
	  new_primal_child->get_dof_values (primal_problem.v, primal_tmp);
	  embedding_matrix.vmult (local_v, primal_tmp);
	  
	  new_dual_child->get_dof_values (dual_problem.u, local_u_bar);
	  new_dual_child->get_dof_values (dual_problem.v, local_v_bar);

	  new_primal_child->set_user_pointer (cellwise_error.next_free_slot);
	  *cellwise_error.next_free_slot
	    = error_formula (new_dual_child,
			     local_u,     local_v,
			     local_u_bar, local_v_bar,
			     child_u_old,     child_v_old,
			     child_u_bar_old, child_v_bar_old,
			     fe_values);
	  ++cellwise_error.next_free_slot;
	};
    };
};



template <int dim>
typename TimeStep_ErrorEstimation<dim>::ErrorOnCell
TimeStep_ErrorEstimation<dim>::collect_error_from_children (const DoFHandler<dim>::cell_iterator &primal_cell_old,
							    const DoFHandler<dim>::cell_iterator &dual_cell_old,
							    const Vector<double>  &local_u,
							    const Vector<double>  &local_v,
							    const Vector<double>  &local_u_bar,
							    const Vector<double>  &local_v_bar,
							    const Vector<double>  &local_Ih_u_bar,
							    const Vector<double>  &local_Ih_v_bar,
							    const Vector<double>  &local_Ih_u_bar_old,
							    const Vector<double>  &local_Ih_v_bar_old,
							    FEValues<dim>  &fe_values) const {
  const TimeStep_Primal<dim> &primal_problem_old = static_cast<const TimeStepBase_Wave<dim>*>
						   (previous_timestep)->get_timestep_primal();
  const TimeStep_Dual<dim>   &dual_problem_old = static_cast<const TimeStepBase_Wave<dim>*>
						 (previous_timestep)->get_timestep_dual();
  const FiniteElement<dim>   &dual_fe          = dual_problem_old.fe;
  
  ErrorOnCell error_sum;

  const unsigned int dofs_per_cell_dual = local_u_bar.size();
  
  for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
    {
      				       // we have the solutions on the
				       // new (large) cell, we restrict it to
				       // each of the small cells
      Vector<double> child_u(dofs_per_cell_dual), child_v(dofs_per_cell_dual);
      Vector<double> child_u_bar(dofs_per_cell_dual), child_v_bar(dofs_per_cell_dual);
      Vector<double> child_Ih_u_bar(dofs_per_cell_dual), child_Ih_v_bar(dofs_per_cell_dual);
      Vector<double> child_Ih_u_bar_old(dofs_per_cell_dual), child_Ih_v_bar_old(dofs_per_cell_dual);      

      dual_fe.prolongate(child).vmult (child_u, local_u);
      dual_fe.prolongate(child).vmult (child_v, local_v);
      dual_fe.prolongate(child).vmult (child_u_bar, local_u_bar);
      dual_fe.prolongate(child).vmult (child_v_bar, local_v_bar);
      dual_fe.prolongate(child).vmult (child_Ih_u_bar, local_Ih_u_bar);
      dual_fe.prolongate(child).vmult (child_Ih_v_bar, local_Ih_v_bar);
      dual_fe.prolongate(child).vmult (child_Ih_u_bar_old, local_Ih_u_bar_old);
      dual_fe.prolongate(child).vmult (child_Ih_v_bar_old, local_Ih_v_bar_old);

      const DoFHandler<dim>::cell_iterator
	old_primal_child = primal_cell_old->child(child),
	old_dual_child   = dual_cell_old->child(child);

      if (old_primal_child->has_children())
					 // the old cell was further
					 // refined -> recurse into the tree
	error_sum += collect_error_from_children (old_primal_child,
						  old_dual_child,
						  child_u,            child_v,
						  child_u_bar,        child_v_bar,
						  child_Ih_u_bar,     child_Ih_v_bar,
						  child_Ih_u_bar_old, child_Ih_v_bar_old,
						  fe_values);
      else
					 // the old cell was not further
					 // refined -> go on here directly
	{
	  Vector<double> local_u_old(dofs_per_cell_dual), local_v_old(dofs_per_cell_dual);
	  Vector<double> local_u_bar_old(dofs_per_cell_dual), local_v_bar_old(dofs_per_cell_dual);

					   // vectors to hold dof values on
					   // the primal/dual cell (temporary)
	  Vector<double> primal_tmp(embedding_matrix.n());
      
				       // fill local old solution vectors.
				       // no problems here, since the two
				       // cells are both unrefined
	  old_primal_child->get_dof_values (primal_problem_old.u, primal_tmp);
	  embedding_matrix.vmult (local_u_old, primal_tmp);
  
	  old_primal_child->get_dof_values (primal_problem_old.v, primal_tmp);
	  embedding_matrix.vmult (local_v_old, primal_tmp);

	  Vector<double> child_difference_u_bar (dofs_per_cell_dual);
	  Vector<double> child_difference_v_bar (dofs_per_cell_dual);
	  Vector<double> local_difference_u_bar_old (dofs_per_cell_dual);
	  Vector<double> local_difference_v_bar_old (dofs_per_cell_dual);

	  child_difference_u_bar =  child_u_bar;
	  child_difference_u_bar -= child_Ih_u_bar;
	  child_difference_v_bar =  child_v_bar;
	  child_difference_v_bar -= child_Ih_v_bar;

	  local_difference_u_bar_old =  local_u_bar_old;
	  local_difference_u_bar_old -= local_Ih_u_bar_old;
	  local_difference_v_bar_old =  local_v_bar_old;
	  local_difference_v_bar_old -= local_Ih_v_bar_old;
	  
	  
	  error_sum += error_formula (old_dual_child,
				      child_u,            child_v,
				      child_u_bar,        child_v_bar,
				      local_u_old,        local_v_old,
				      local_u_bar_old,    local_v_bar_old,
				      fe_values);
	};
    };

  return error_sum;
};



template <int dim>
typename TimeStep_ErrorEstimation<dim>::ErrorOnCell
TimeStep_ErrorEstimation<dim>::error_formula (const DoFHandler<dim>::active_cell_iterator &cell,
					      const Vector<double>  &local_u,
					      const Vector<double>  &local_v,
					      const Vector<double>  &local_u_bar,
					      const Vector<double>  &local_v_bar,
					      const Vector<double>  &local_u_old,
					      const Vector<double>  &local_v_old,
					      const Vector<double>  &local_u_bar_old,
					      const Vector<double>  &local_v_bar_old,
					      FEValues<dim>  &fe_values) const {
  Vector<double> local_difference_u_bar(local_u_bar.size());
  Vector<double> local_difference_v_bar(local_u_bar.size());
  Vector<double> local_difference_u_bar_old(local_u_bar.size());
  Vector<double> local_difference_v_bar_old(local_u_bar.size());
  
  difference_matrix.vmult (local_difference_u_bar, local_u_bar);
  difference_matrix.vmult (local_difference_v_bar, local_v_bar);
  difference_matrix.vmult (local_difference_u_bar_old, local_u_bar_old);
  difference_matrix.vmult (local_difference_v_bar_old, local_v_bar_old);

  return error_formula (cell,
			local_u,            local_v,
			local_u_bar,        local_v_bar,
			local_u_old,        local_v_old,
			local_u_bar_old,    local_v_bar_old,
			local_difference_u_bar,
			local_difference_v_bar,
			local_difference_u_bar_old,
			local_difference_v_bar_old,
			fe_values);			
};



template <int dim>
typename TimeStep_ErrorEstimation<dim>::ErrorOnCell
TimeStep_ErrorEstimation<dim>::error_formula (const DoFHandler<dim>::active_cell_iterator &cell,
					      const Vector<double>  &local_u,
					      const Vector<double>  &local_v,
					      const Vector<double>  &local_u_bar,
					      const Vector<double>  &local_v_bar,
					      const Vector<double>  &local_u_old,
					      const Vector<double>  &local_v_old,
					      const Vector<double>  &local_u_bar_old,
					      const Vector<double>  &local_v_bar_old,
					      const Vector<double>  &local_difference_u_bar,
					      const Vector<double>  &local_difference_v_bar,
					      const Vector<double>  &local_difference_u_bar_old,
					      const Vector<double>  &local_difference_v_bar_old,
					      FEValues<dim>         &fe_values) const {

				   // this will be used to sum up the
				   // different parts of the error
				   // identity on this cell
  ErrorOnCell error_on_cell;

  const unsigned int dofs_per_cell = get_timestep_dual().fe.dofs_per_cell;
  
				   // two temporaries needed for the
				   // calculation of the scalar products
  Vector<double> tmp1(dofs_per_cell);
  Vector<double> tmp2(dofs_per_cell);


  vector<double> stiffness(fe_values.n_quadrature_points);
  parameters.stiffness->value_list (fe_values.get_quadrature_points(),
				    stiffness);
  vector<Tensor<1,dim> > grad_stiffness(fe_values.n_quadrature_points);
  parameters.stiffness->gradient_list (fe_values.get_quadrature_points(),
				       grad_stiffness);

				   // matrix for (phi_i, phi_j)
  FullMatrix<double> mass_matrix (tmp1.size(), tmp1.size());
				   // matrix for (a\Delta phi_i, phi_j)
//  FullMatrix<double> delta_matrix (tmp1.size(), tmp1.size());
				   // matrix for (grad a . grad phi_i, phi_j)
//  FullMatrix<double> grad_grad_matrix (tmp1.size(), tmp1.size());
  FullMatrix<double> laplace_matrix (tmp1.size(), tmp1.size());
  
				   // first task: create matrices
  fe_values.reinit (cell);
  const FullMatrix<double>             &values    = fe_values.get_shape_values();
  const vector<vector<Tensor<1,dim> > >&gradients = fe_values.get_shape_grads ();
  const vector<vector<Tensor<2,dim> > >&second_derivatives
    = fe_values.get_shape_2nd_derivatives ();
  const vector<double>                 &weights   = fe_values.get_JxW_values ();

  vector<double> density_values(fe_values.n_quadrature_points);
  parameters.density->value_list (fe_values.get_quadrature_points(),
				  density_values);
  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
    for (unsigned int i=0; i<dofs_per_cell; ++i) 
      for (unsigned int j=0; j<dofs_per_cell; ++j)
	{
	  mass_matrix(i,j) += (values(i,point) *
			       values(j,point)) *
			      weights[point] *
			      density_values[point];

					   // compute laplacian of phi_i
					   // by summing the trace of the
					   // tensor of second derivatives
	  double laplace_phi_i = 0;
	  for (unsigned int t=0; t<dim; ++t)
	    laplace_phi_i += second_derivatives[i][point][t][t];
	  
// 	  delta_matrix(i,j) += stiffness[point] *
// 			       laplace_phi_i *
// 			       values(j,point) *
// 			       weights[point];

// 	  grad_grad_matrix(i,j) += (grad_stiffness[point] *
// 				    gradients[i][point]) *
// 				   values(j,point) *
// 				   weights[point];

	  laplace_matrix(i,j) += (gradients[i][point] *
				  gradients[j][point]) *
				 weights[point] *
				 stiffness[point];
	};
  
			       
				   
				   // ////////////////////////////////////
				   // Compute the different contributions
				   // separately. Note that the parts
				   // 1 and 2a+2b together should give
				   // a small quantity, since they form
				   // the first domain residual, which is
				   // small for elements of odd order.
  
      
				   // ////////////////////////////////////
				   // PART 1
				   //
				   // let #tmp_dual2# hold the contribution
				   // 1/2 (1-I)(u_bar^n + u_bar^(n-1))
				   // with I the interpolation operator
  tmp2  = local_difference_u_bar;
  tmp2 += local_difference_u_bar_old;
  tmp2.scale (1./2.);
  
				   // let #tmp_dual1# hold
				   // u^n - u^(n-1)
  tmp1  = local_u;
  tmp1 -= local_u_old;
  
  error_on_cell.part[0] = mass_matrix.matrix_scalar_product (tmp1, tmp2);
  
  
				   // same thing for the second part
				   // with v instead of u
  tmp2  = local_difference_v_bar;
  tmp2 += local_difference_v_bar_old;
  tmp2.scale (1./2.);
  
  tmp1  = local_v;
  tmp1 -= local_v_old;
  
  error_on_cell.part[1] = mass_matrix.matrix_scalar_product (tmp1, tmp2);
  

  
				   // ////////////////////////////////
				   // PART 2a
				   //
				   // let tmp2=(1-I)(u_bar^n+u_bar^(n-1))
				   // let tmp1 = v^n+v^(n-1)
  tmp2  = local_difference_u_bar;
  tmp2 += local_difference_u_bar_old;

  tmp1  = local_v;
  tmp1 += local_v_old;

  error_on_cell.part[2] = -(get_backward_timestep() / 4 *
			    mass_matrix.matrix_scalar_product (tmp1, tmp2));
  

				   // ////////////////////////////////
				   // PART 2b
				   //
				   // let tmp1 = v^n-v^(n-1)
				   // let tmp2=u_bar^n - u_bar^(n-1)
  tmp1  = local_v;
  tmp1 -= local_v_old;

  tmp2  = local_u_bar;
  tmp2 -= local_u_bar_old;

  error_on_cell.part[3] = -(get_backward_timestep() / 12 *
			    mass_matrix.matrix_scalar_product (tmp1, tmp2));


				   // ////////////////////////////////
				   // PART 3a
				   //
				   // let tmp2=(1-I)(v_bar^n+v_bar^(n-1))
				   // let tmp1 = u^n+u^(n-1)
  tmp2  = local_difference_v_bar;
  tmp2 += local_difference_v_bar_old;

  tmp1  = local_u;
  tmp1 += local_u_old;

  error_on_cell.part[4] = (get_backward_timestep() / 4 *
			   laplace_matrix.matrix_scalar_product (tmp1, tmp2));


				   // ////////////////////////////////
				   // PART 3b
				   //
				   // let tmp1 = u^n-u^(n-1)
				   // let tmp2 = (v_bar^n - v_bar^(n-1))
  tmp1  = local_u;
  tmp1 -= local_u_old;

  tmp2  = local_v_bar;
  tmp2 -= local_v_bar_old;

  error_on_cell.part[5] = (get_backward_timestep() / 12 *
			   laplace_matrix.matrix_scalar_product (tmp1, tmp2));



// 				   // ///////////////////////////
// 				   // PART 0:
// 				   // tmp1 = u^n-u^(n-1)
// 				   // tmp2 = 1/2 (1-I) (u_bar^n + u_bar^(n-1)
//   tmp1  = local_u;
//   tmp1 -= local_u_old;

//   tmp2  = local_difference_u_bar;
//   tmp2 += local_difference_u_bar_old;
//   tmp2.scale (1./2.);

//   error_on_cell.part[0] = -1. * mass_matrix.matrix_scalar_product (tmp1, tmp2);

// 				   // ////////////////////////////
// 				   // PART 1:
// 				   // same as above, but with u and
// 				   // v interchanged
//   tmp1  = local_v;
//   tmp1 -= local_v_old;

//   tmp2  = local_difference_v_bar;
//   tmp2 += local_difference_v_bar_old;
//   tmp2.scale (1./2.);

//   error_on_cell.part[1] = -1. * mass_matrix.matrix_scalar_product (tmp1, tmp2);


// 				   // /////////////////////////////
// 				   // PART 2:
// 				   // tmp1 = v^n+v^(n-1)
// 				   // tmp2 = (1-I) (u_bar^n + u_bar^(n-1))
//   tmp1  = local_v;
//   tmp1 += local_v_old;

//   tmp2  = local_difference_u_bar;
//   tmp2 += local_difference_u_bar_old;

//   error_on_cell.part[2]  = mass_matrix.matrix_scalar_product (tmp1, tmp2);
//   error_on_cell.part[2] *= get_backward_timestep() / 4;

// 				   // //////////////////////////////
// 				   // PART 3:
// 				   // tmp1 = v^n-v^(n-1)
// 				   // tmp2 = u_bar^n - u_bar^(n-1)
//   tmp1  = local_v;
//   tmp1 -= local_v_old;

//   tmp2  = local_u_bar;
//   tmp2 -= local_u_bar_old;


//   error_on_cell.part[3]  = mass_matrix.matrix_scalar_product (tmp1, tmp2);
//   error_on_cell.part[3] *= get_backward_timestep() / 12;
  

// 				   // /////////////////////////////
// 				   // PART 4 and 6:
// 				   // tmp1 = u^n+u^(n-1)
// 				   // tmp2 = (1-I) (v_bar^n + v_bar^(n-1))
//   tmp1  = local_u;
//   tmp1 += local_u_old;

//   tmp2  = local_difference_v_bar;
//   tmp2 += local_difference_v_bar_old;

//   error_on_cell.part[4]  = delta_matrix.matrix_scalar_product (tmp1, tmp2);
//   error_on_cell.part[4] *= get_backward_timestep() / 4;

//   error_on_cell.part[6]  = grad_grad_matrix.matrix_scalar_product (tmp1, tmp2);
//   error_on_cell.part[6] *= get_backward_timestep() / 12;

// 				   // //////////////////////////////
// 				   // PART 5 and 7:
// 				   // tmp1 = u^n-u^(n-1)
// 				   // tmp2 = v_bar^n - v_bar^(n-1)
//   tmp1  = local_u;
//   tmp1 -= local_u_old;

//   tmp2  = local_v_bar;
//   tmp2 -= local_v_bar_old;


//   error_on_cell.part[5]  = delta_matrix.matrix_scalar_product (tmp1, tmp2);
//   error_on_cell.part[5] *= get_backward_timestep() / 12;

//   error_on_cell.part[7]  = grad_grad_matrix.matrix_scalar_product (tmp1, tmp2);
//   error_on_cell.part[7] *= get_backward_timestep() / 12;

  
  return error_on_cell;
};






template <int dim>
void TimeStep_ErrorEstimation<dim>::make_interpolation_matrices () {
  const FiniteElement<dim> &primal_fe = get_timestep_primal().fe,
			   &dual_fe   = get_timestep_dual().fe;
  
  embedding_matrix.reinit (dual_fe.dofs_per_cell,
			   primal_fe.dofs_per_cell);

  vector<Point<dim> > unit_support_points (dual_fe.dofs_per_cell);
  dual_fe.get_unit_support_points (unit_support_points);
  
  for (unsigned int i=0; i<dual_fe.dofs_per_cell; ++i)
    for (unsigned int j=0; j<primal_fe.dofs_per_cell; ++j)
      embedding_matrix(i,j) = primal_fe.shape_value (j, unit_support_points[i]);


				   // construct the difference between the
				   // identity and the interpolation operator
				   // to the primal ansatz space. The
				   // interpolation operator is to act from
				   // and to the dual space (not as above
				   // where it acted from one space into
				   // the other), so we construct it by
				   // first interpolating down to the
				   // primal space and then back to the
				   // dual space (by injection, using
				   // the matrix constructed above)
  FullMatrix<double> inverse_interpolation (primal_fe.dofs_per_cell,
					    dual_fe.dofs_per_cell);
  unit_support_points.resize (primal_fe.dofs_per_cell);
  primal_fe.get_unit_support_points (unit_support_points);
  
  for (unsigned int i=0; i<primal_fe.dofs_per_cell; ++i)
    for (unsigned int j=0; j<dual_fe.dofs_per_cell; ++j)
      inverse_interpolation(i,j) = dual_fe.shape_value (j, unit_support_points[i]);

  interpolation_matrix.reinit (dual_fe.dofs_per_cell, dual_fe.dofs_per_cell);
  embedding_matrix.mmult (interpolation_matrix, inverse_interpolation);
  
  difference_matrix.reinit (dual_fe.dofs_per_cell, dual_fe.dofs_per_cell);
				   // initialize with the unit matrix
  for (unsigned int i=0; i<dual_fe.dofs_per_cell; ++i)
    difference_matrix(i,i) = 1.;
				   // compute difference
  difference_matrix.add (-1, interpolation_matrix);
};





template <int dim>
TimeStep_ErrorEstimation<dim>::StatisticData::StatisticData () :
		estimated_error (0)
{};



template <int dim>
TimeStep_ErrorEstimation<dim>::StatisticData::StatisticData (const double estimated_error) :
		estimated_error (estimated_error)
{};



template <int dim>
void TimeStep_ErrorEstimation<dim>::StatisticData::write_descriptions (ostream &out)
{
  out << "#    total estimated error in this timestep" << endl;
};



template <int dim>
void TimeStep_ErrorEstimation<dim>::StatisticData::write (ostream &out) const
{
  out << estimated_error;
};



template <int dim>
TimeStep_ErrorEstimation<dim>::ErrorOnCell::ErrorOnCell () {
  for (unsigned int i=0; i<sizeof(part)/sizeof(part[0]); ++i)
    part[i] = 0;
};



template <int dim>
typename TimeStep_ErrorEstimation<dim>::ErrorOnCell
TimeStep_ErrorEstimation<dim>::ErrorOnCell::operator += (const ErrorOnCell &eoc) {
  for (unsigned int i=0; i<sizeof(part)/sizeof(part[0]); ++i)
    part[i] += eoc.part[i];
  return *this;
};



template <int dim>
double TimeStep_ErrorEstimation<dim>::ErrorOnCell::sum () const {
  double x=0;
  for (unsigned int i=0; i<sizeof(part)/sizeof(part[0]); ++i)
    x += part[i];
  return x;
};



template <int dim>
TimeStep_ErrorEstimation<dim>::CellwiseError::CellwiseError (const unsigned int n_errors) :
		errors (n_errors),
		next_free_slot (errors.begin())
{};



// explicit instantiations
template class TimeStep_ErrorEstimation<2>;
/* $Id$ */


#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_handler.h>




template <int dim>
TimeStep<dim>::TimeStep (const double               time,
			 const WaveParameters<dim> &parameters):
		TimeStepBase_Wave<dim> (time,
					TimeStepBase_Tria<dim>::Flags(true, 0, 1),
					parameters),
		TimeStep_Primal<dim>(parameters.primal_fe),
		TimeStep_Dual<dim>  (parameters.dual_fe)
{};




template <int dim>
void TimeStep<dim>::wake_up (const unsigned int wakeup_level)
{
  sweep_info->get_timers().grid_generation.start();
  TimeStepBase_Wave<dim>::wake_up (wakeup_level);
  sweep_info->get_timers().grid_generation.stop();

  switch (next_action)
    {
      case primal_problem:
	    TimeStep_Primal<dim>::wake_up (wakeup_level);
	    break;
	    
      case dual_problem:
	    TimeStep_Dual<dim>::wake_up (wakeup_level);
	    break;
	    
      case postprocess:	    
	    TimeStep_Primal<dim>::wake_up (wakeup_level);

	    if ((parameters.refinement_strategy == WaveParameters<dim>::dual_estimator)
		&&
		(sweep_no >= parameters.initial_energy_estimator_sweeps))
	      TimeStep_Dual<dim>::wake_up (wakeup_level);
	    
	    TimeStep_Postprocess<dim>::wake_up (wakeup_level);
	    
	    break;

      case grid_refinement:
					     // do nothing except for waking
					     // up the grid
	    break;

      default:
	    Assert (false, ExcInternalError());
    };
};



template <int dim>
void TimeStep<dim>::sleep (const unsigned int sleep_level)
{  
  switch (next_action)
    {
      case primal_problem:
	    TimeStep_Primal<dim>::sleep (sleep_level);
	    break;
	    
      case dual_problem:
	    TimeStep_Dual<dim>::sleep (sleep_level);
	    break;

      case postprocess:
	    TimeStep_Primal<dim>::sleep (sleep_level);
	    
	    if ((parameters.refinement_strategy == WaveParameters<dim>::dual_estimator)
		&&
		(sweep_no >= parameters.initial_energy_estimator_sweeps))
	      TimeStep_Dual<dim>::sleep (sleep_level);
	    
	    TimeStep_Postprocess<dim>::sleep (sleep_level);
	    break;

      case grid_refinement:
					     // save the flags since the grid
					     // will be deleted next along with
					     // the flags
	    if (sleep_level == 1)
	      save_refine_flags ();
	    break;

      default:
	    Assert (false, ExcInternalError());
    };

  sweep_info->get_timers().grid_generation.start();
  TimeStepBase_Wave<dim>::sleep (sleep_level);
  sweep_info->get_timers().grid_generation.stop();
};



template <int dim>
void TimeStep<dim>::end_sweep ()
{
  TimeStep_Primal<dim>::end_sweep ();
  TimeStep_Dual<dim>::end_sweep ();
  TimeStep_Postprocess<dim>::end_sweep ();
};



template <int dim>
void TimeStep<dim>::write_statistics_descriptions (ostream                   &out,
						   const WaveParameters<dim> &parameters)
{
  out << "#  Primal problem:" << endl;
  TimeStep_Primal<dim>::StatisticData::write_descriptions (out);

  out << "#  Dual problem:" << endl;
  TimeStep_Dual<dim>::StatisticData::write_descriptions (out);

  out << "#  Error estimation:" << endl;
  TimeStep_ErrorEstimation<dim>::StatisticData::write_descriptions (out);

  if (parameters.eval_list.size() != 0)
    {
      out << "#  Postprocessing:" << endl;
      TimeStep_Postprocess<dim>::StatisticData::write_descriptions (out, parameters);
    };
};



template <int dim>
void TimeStep<dim>::write_statistics (ostream &out) const 
{
  TimeStep_Primal<dim>::statistic_data.write (out);
  out << "    ";
  TimeStep_Dual<dim>::statistic_data.write (out);
  out << "    ";
  TimeStep_ErrorEstimation<dim>::statistic_data.write (out);
  out << "    ";
  TimeStep_Postprocess<dim>::statistic_data.write (out);
};



// explicit instantiations
template class TimeStep<2>;
/* $Id$ */


#include <lac/vector.h>
#include <numerics/data_out.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/geometry_info.h>
#include <numerics/data_out_stack.h>

#include <fstream>
#include <iomanip>



template <int dim>
void TimeStep_Postprocess<dim>::postprocess_timestep () 
{
  deallog << "  Postprocessing: time="
       << time
       << ", step=" << timestep_no
       << ", sweep=" << sweep_no
       << ". ";

  if ((sweep_no < parameters.number_of_sweeps-1) ||
      (parameters.refinement_strategy == WaveParameters<dim>::dual_estimator))
    estimate_error ();

				   // the error estimator has its own timer,
				   // so start the postprocessing timer
				   // only here
  sweep_info->get_timers().postprocessing.start();

    				   // do the user evaluations
  statistic_data.evaluation_results.clear();
  for (typename list<EvaluationBase<dim>*>::const_iterator i = parameters.eval_list.begin();
       i != parameters.eval_list.end(); ++i)
    {
      (*i)->reset_timelevel (get_timestep_primal());
      statistic_data.evaluation_results.push_back ((*i)->evaluate());
    };

				   // write data if requested
  if (((parameters.write_solution_strategy == WaveParameters<dim>::all_sweeps) ||
       ((parameters.write_solution_strategy == WaveParameters<dim>::last_sweep_only) &&
	(sweep_no == parameters.number_of_sweeps-1)))
      &&
      (((timestep_no % parameters.write_steps_interval) == 0) ||
       (next_timestep == 0)))
    {
      deallog << "[o]";

      DataOut<dim>::OutputFormat output_format
	= DataOut<dim>::parse_output_format (parameters.output_format);
      
      string data_filename	= (parameters.output_directory +
				   "sweep" + int_to_string(sweep_no,2) +
				   "/" + int_to_string(timestep_no,4) +
				   DataOut<dim>::default_suffix (output_format));
      DataOut<dim> out;
      out.attach_dof_handler (*get_timestep_primal().dof_handler);
      out.add_data_vector (get_timestep_primal().u, "u");
      out.add_data_vector (get_timestep_primal().v, "v");

				       // vectors holding the dual variables,
				       // if needed
      Vector<double> u_bar, v_bar;
      
				       // if dual problem was computed
      if ((parameters.refinement_strategy == WaveParameters<dim>::dual_estimator)
	  &&
	  (sweep_no >= parameters.initial_energy_estimator_sweeps))
	{
	  u_bar.reinit (get_timestep_primal().u.size());
	  v_bar.reinit (get_timestep_primal().u.size());
	  
	  if (parameters.primal_fe == parameters.dual_fe)
					     // if primal and dual solution
					     // were computed using the same
					     // ansatz, we may add the dual
					     // solutions "as is"
	    {
	      u_bar = get_timestep_dual().u;
	      v_bar = get_timestep_dual().v;
	    }
	  else
					     // otherwise: first interpolate
					     // the dual solutions to the
					     // same degree
	    interpolate_dual_solution (u_bar, v_bar);
	  
	  out.add_data_vector (u_bar, "dual_u");
	  out.add_data_vector (v_bar, "dual_v");
	};

				       // add error vector if error
				       // was computed
      Vector<double> estimated_error;
      if ((sweep_no<parameters.number_of_sweeps-1) ||
	  (parameters.refinement_strategy == WaveParameters<dim>::dual_estimator))
	{
	  if (parameters.write_error_as_cell_data) 
	    {
	      estimated_error.reinit (estimated_error_per_cell.size());
	      copy_n (estimated_error_per_cell.begin(),
		      estimated_error_per_cell.size(),
		      estimated_error.begin());
	    }
	  else
	    {
	      estimated_error.reinit (get_timestep_primal().dof_handler->n_dofs());
	      DoFTools::distribute_cell_to_dof_vector (*get_timestep_primal().dof_handler,
						       estimated_error_per_cell,
						       estimated_error);
	    };
      
	  out.add_data_vector (estimated_error, "est_error");
	};

      out.build_patches ();
      
      out.write (logfile, output_format);

      deallog << ".";
    };

  if (parameters.write_stacked_data &&
      (timestep_no % parameters.write_stacked_interval == 0))
    {
      deallog << "[st]";

      sweep_data->data_out_stack->new_parameter_value (time,
						       (timestep_no == 0 ?
							0 :
							get_backward_timestep() *
							parameters.write_stacked_interval));
      sweep_data->data_out_stack->attach_dof_handler (*get_timestep_primal().dof_handler);
      sweep_data->data_out_stack->add_data_vector (get_timestep_primal().u, "u");
      sweep_data->data_out_stack->add_data_vector (get_timestep_primal().v, "v");

				       // if dual problem was computed
      if ((parameters.refinement_strategy == WaveParameters<dim>::dual_estimator)
	  &&
	  (sweep_no >= parameters.initial_energy_estimator_sweeps))
	{
	  if (parameters.primal_fe == parameters.dual_fe)
					     // if primal and dual solution
					     // were computed using the same
					     // ansatz, we may add the dual
					     // solutions "as is"
	    {
	      sweep_data->data_out_stack->add_data_vector (get_timestep_dual().u, "dual_u");
	      sweep_data->data_out_stack->add_data_vector (get_timestep_dual().v, "dual_v");
	    }
	  else
					     // otherwise: first interpolate
					     // the dual solutions to the
					     // same degree
	    {
	      Vector<double> u_bar(get_timestep_primal().dof_handler->n_dofs());
	      Vector<double> v_bar(get_timestep_primal().dof_handler->n_dofs());
	      
	      interpolate_dual_solution (u_bar, v_bar);
	      
	      sweep_data->data_out_stack->add_data_vector (u_bar, "dual_u");
	      sweep_data->data_out_stack->add_data_vector (v_bar, "dual_v");
	    };
	};

				       // add error estimator if that was
				       // computed
      if ((sweep_no < parameters.number_of_sweeps-1) ||
	  (parameters.refinement_strategy == WaveParameters<dim>::dual_estimator))
	sweep_data->data_out_stack->add_data_vector (estimated_error_per_cell, "est_error");

      sweep_data->data_out_stack->build_patches ();
      sweep_data->data_out_stack->finish_parameter_value ();
    };
  
	
  deallog << endl;
  sweep_info->get_timers().postprocessing.stop();
};



template <int dim>
void TimeStep_Postprocess<dim>::wake_up (const unsigned int wakeup_level) 
{
  TimeStep_ErrorEstimation<dim>::wake_up (wakeup_level);
};



template <int dim>
void TimeStep_Postprocess<dim>::sleep (const unsigned int sleep_level) 
{
  TimeStep_ErrorEstimation<dim>::sleep (sleep_level);
};



template <int dim>
string TimeStep_Postprocess<dim>::branch_signature () const 
{
  return "o";
};



template <int dim>
void TimeStep_Postprocess<dim>::end_sweep ()
{
  string tmp_filename = tmp_filename_base(branch_signature());
  remove (tmp_filename.c_str());
};


    



template <int dim>
void TimeStep_Postprocess<dim>::interpolate_dual_solution (Vector<double> &interpolated_u_bar,
							   Vector<double> &interpolated_v_bar) const {
  const unsigned int n_primal_dofs = get_timestep_primal().dof_handler->n_dofs();
  
  interpolated_u_bar.reinit (n_primal_dofs);
  interpolated_v_bar.reinit (n_primal_dofs);

  const TimeStep_Wave<dim> &target = get_timestep_dual ();
  
  typename DoFHandler<dim>::active_cell_iterator primal_cell, dual_cell, endc;
  primal_cell = get_timestep_primal().dof_handler->begin_active();
  endc        = get_timestep_primal().dof_handler->end();
  dual_cell   = target.dof_handler->begin_active();
  
				   // loop over all cells and set the vertex
				   // values of the interpolated vector to
				   // the vertex values of the dual solutions.
				   // don't care that we set these values
				   // more than once...
  for (; primal_cell != endc; ++primal_cell, ++dual_cell)
    for (unsigned int vertex=0; vertex<GeometryInfo<dim>::vertices_per_cell; ++vertex)
      {
	const unsigned int  primal_vertex_index = primal_cell->vertex_dof_index(vertex,0),
			    dual_vertex_index   = dual_cell->vertex_dof_index(vertex,0);
	interpolated_u_bar(primal_vertex_index) = target.u(dual_vertex_index);
	interpolated_v_bar(primal_vertex_index) = target.v(dual_vertex_index);
      };      
};




template <int dim>
void TimeStep_Postprocess<dim>::StatisticData::
write_descriptions (ostream &out,
		    const WaveParameters<dim> &parameters) 
{
  for (typename list<EvaluationBase<dim>*>::const_iterator i = parameters.eval_list.begin();
       i != parameters.eval_list.end(); ++i)
    out << "#    " << (*i)->description() << endl;
};



template <int dim>
void TimeStep_Postprocess<dim>::StatisticData::write (ostream &out) const
{
  for (unsigned int i=0; i<evaluation_results.size(); ++i)
    out << evaluation_results[i] << ' ';
};




// explicit instantiations
template class TimeStep_Postprocess<2>;
/* $Id$ */

#include <base/function.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <lac/sparse_matrix.h>
#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <fe/fe_values.h>
#include <fe/fe.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>


#include <iomanip>



template <int dim>
TimeStep_Primal<dim>::TimeStep_Primal (const string &primal_fe)
		:
		TimeStep_Wave<dim> (primal_fe)
{};

		


template <int dim>
void TimeStep_Primal<dim>::do_initial_step ()
{
  deallog << "  Primal problem: time="
       << time
       << ", step=" << timestep_no
       << ", sweep=" << sweep_no
       << ". "
       << tria->n_active_cells() << " cells, "
       << dof_handler->n_dofs() << " dofs";


				   // add up sweep-accumulated data. count
				   // u and v as separate dofs
  sweep_info->get_data().cells       += tria->n_active_cells();
  sweep_info->get_data().primal_dofs += dof_handler->n_dofs() * 2;

				   // use L2-projection for u0 and v0
#if 2 == 1
  VectorTools::interpolate (*dof_handler, *parameters.initial_u, u);
  VectorTools::interpolate (*dof_handler, *parameters.initial_v, v);
#else
  VectorTools::project (*dof_handler, constraints,
			     quadrature, *parameters.initial_u, u,
			     false, quadrature_face, (dim==2 ? true : false));
  VectorTools::project (*dof_handler, constraints,
			     quadrature, *parameters.initial_v, v,
			     false, quadrature_face, (dim==2 ? true : false));
#endif
				   // set energy to zero since we
				   // don't want to assemble the matrices
				   // needed for this
  statistic_data = StatisticData (tria->n_active_cells(),
				  dof_handler->n_dofs(),
				  0,
				  0,
				  make_pair (0.0, 0.0));

  deallog << "." << endl;
};




template <int dim>
void TimeStep_Primal<dim>::do_timestep ()
{
  deallog << "  Primal problem: time="
       << time
       << ", step=" << timestep_no
       << ", sweep=" << sweep_no
       << ". "
       << tria->n_active_cells() << " cells, "
       << dof_handler->n_dofs() << " dofs";
  
				   // add up sweep-accumulated data. count
				   // u and v as separate dofs
  sweep_info->get_data().cells       += tria->n_active_cells();
  sweep_info->get_data().primal_dofs += dof_handler->n_dofs() * 2;
  
  
  const double time_step = get_backward_timestep ();
  
				   // Vectors holding the right hand sides of
				   // the two equations.
  Vector<double> right_hand_side1 (dof_handler->n_dofs());
  Vector<double> right_hand_side2 (dof_handler->n_dofs());

				   // Vector holding a the values for
				   // u and v of the previous time step.
				   // these are used in case we want to
				   // use extrapolation from the previous
				   // time step to the present one
  Vector<double> old_u, old_v;
  if (parameters.extrapolate_old_solutions)
    {
      old_u.reinit (dof_handler->n_dofs());
      old_v.reinit (dof_handler->n_dofs());

      transfer_old_solutions (old_u, old_v);
    };

  
  assemble_vectors (right_hand_side1, right_hand_side2);

  UserMatrix system_matrix (system_sparsity, parameters.preconditioning);
  system_matrix.copy_from (mass_matrix);
  system_matrix.add_scaled (time_step * time_step *
			    parameters.theta *
			    parameters.theta,
			    laplace_matrix);
  constraints.condense (system_matrix);
	
  if (parameters.extrapolate_old_solutions)
				     // solve with a hopefully good guess
				     // as start vector
    {
      u  = old_u;
      u.add (time_step, old_v);
    };

				   // in 1d, do not set boundary conditions
				   // at all
  if (dim!=1) 
    {
				       // in the other case, the wake_up
				       // function of the base class has set
				       // the solution vector's values to
				       // zero already.
      parameters.boundary_values_u->set_time (time);
      parameters.boundary_values_v->set_time (time);
      
      map<int,double> boundary_value_list;
      VectorTools::interpolate_boundary_values (*dof_handler, 0,
						     *(parameters.boundary_values_u),
						     boundary_value_list);
      MatrixTools<dim>::apply_boundary_values (boundary_value_list,
					       system_matrix, u,
					       right_hand_side1);
    };

  const unsigned int solver_steps1 = solve (system_matrix, u, right_hand_side1);
		
  system_matrix.copy_from (mass_matrix);
  constraints.condense (system_matrix);
  if (true)
    { 
      Vector<double> tmp (right_hand_side2.size());
      laplace_matrix.vmult (tmp, u);
      right_hand_side2.add (-parameters.theta*time_step, tmp);
    };
  constraints.condense (right_hand_side2);


				   // in 1d, do not set boundary conditions
				   // at all
  if (dim != 1)
    {
      map<int,double> boundary_value_list;
      VectorTools::interpolate_boundary_values (*dof_handler, 0,
						     *(parameters.boundary_values_v),
						     boundary_value_list);
      MatrixTools<dim>::apply_boundary_values (boundary_value_list,
					       system_matrix, v,
					       right_hand_side2);
    };
  

  if (parameters.extrapolate_old_solutions)
				     // solve with a hopefully good guess
				     // as start vector
    {
      v  = u;
      v -= old_u;
      v.scale (2./time_step);
      v -= old_v;
    };

  const unsigned int solver_steps2 = solve (system_matrix, v, right_hand_side2);

  statistic_data = StatisticData (tria->n_active_cells(),
				  dof_handler->n_dofs(),
				  solver_steps1,
				  solver_steps2,
				  compute_energy ());
  
  deallog << "." << endl;
};




template <int dim>
void TimeStep_Primal<dim>::solve_primal_problem ()
{
  sweep_info->get_timers().primal_problem.start();
  if (timestep_no == 0)
    do_initial_step ();
  else
    do_timestep ();
  sweep_info->get_timers().primal_problem.stop();
};



template <int dim>
string TimeStep_Primal<dim>::branch_signature () const 
{
  return "p";
};



template <int dim>
void TimeStep_Primal<dim>::wake_up (const unsigned int wakeup_level)
{
  TimeStep_Wave<dim>::wake_up (wakeup_level);
  
  sweep_info->get_timers().primal_problem.start();
  if ((wakeup_level==0) && (next_action==primal_problem))
    {
      Assert (system_sparsity.empty(), ExcInternalError());
      
      create_matrices ();
    };
  sweep_info->get_timers().primal_problem.stop();
};



template <int dim>
void TimeStep_Primal<dim>::assemble_vectors (Vector<double> &right_hand_side1,
					     Vector<double> &right_hand_side2) {
				   // don't do some things for the initial
				   // step since we don't need them there
  Assert (timestep_no>=1, ExcInternalError());
  
				   // construct right hand side
  build_rhs (right_hand_side1, right_hand_side2);
				   // condense right hand side in-place
  constraints.condense (right_hand_side1);
};








template <int dim>
void TimeStep_Primal<dim>::build_rhs (Vector<double> &right_hand_side1,
				      Vector<double> &right_hand_side2) {
				   // select the TimeStep_Wave part in the
				   // TimeStep_Primal branch
  const TimeStep_Primal<dim> &previous_time_level
    = static_cast<const TimeStepBase_Wave<dim>*>(previous_timestep)->get_timestep_primal();

  Assert (previous_time_level.tria->n_cells(0) == tria->n_cells(0),
	  ExcCoarsestGridsDiffer());

				   // convenience typedef
  typedef DoFHandler<dim>::cell_iterator cell_iterator;

				   // create this here and pass it to
				   // the cellwise function since it
				   // is expensive to create it for
				   // every cell
  FEValues<dim> fe_values (fe, quadrature,
			   UpdateFlags(update_values |
				       update_gradients |
				       update_JxW_values |
				       update_q_points));

  
  cell_iterator old_cell = previous_time_level.dof_handler->begin(),
		new_cell = dof_handler->begin(),
		end_cell = (tria->n_levels() == 1                  ?
			    static_cast<cell_iterator>(dof_handler->end()) :
			    dof_handler->begin(1));  
  for (; new_cell!=end_cell; ++new_cell, ++old_cell)
    build_rhs (old_cell, new_cell,
	       fe_values,
	       right_hand_side1, right_hand_side2);
};



template <int dim>
void
TimeStep_Primal<dim>::build_rhs (const DoFHandler<dim>::cell_iterator &old_cell,
				 const DoFHandler<dim>::cell_iterator &new_cell,
				 FEValues<dim>        &fe_values,
				 Vector<double>       &right_hand_side1,
				 Vector<double>       &right_hand_side2) {
				   // declare this type for convenience
  typedef DoFHandler<dim>::cell_iterator cell_iterator;
  
				   // both cells have children, so
				   // recurse into the tree
  if (old_cell->has_children() && new_cell->has_children()) 
    {
      for (unsigned int child=0; child<GeometryInfo<dim>::children_per_cell; ++child)
	build_rhs (old_cell->child(child),
		   new_cell->child(child),
		   fe_values,
		   right_hand_side1,
		   right_hand_side2);
      return;
    };


				   // select the TimeStep_Wave part in the
				   // TimeStep_Primal branch
  const TimeStep_Primal<dim> &previous_time_level
    = static_cast<const TimeStepBase_Wave<dim>*>(previous_timestep)->get_timestep_primal();

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const double time_step = get_backward_timestep();

				   // both cells are on the same refinement
				   // level
  if (!old_cell->has_children() && !new_cell->has_children()) 
    {
      fe_values.reinit (old_cell);

      FullMatrix<double>                    cell_matrix (dofs_per_cell, dofs_per_cell);
      const FullMatrix<double>             &values    = fe_values.get_shape_values ();
      const vector<vector<Tensor<1,dim> > >&gradients = fe_values.get_shape_grads ();
      const vector<double>                 &weights   = fe_values.get_JxW_values ();

				       // assemble mass matrix
      vector<double> density_values(fe_values.n_quadrature_points);
      parameters.density->value_list (fe_values.get_quadrature_points(),
				      density_values);
      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	for (unsigned int i=0; i<dofs_per_cell; ++i) 
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    cell_matrix(i,j) += (values(i,point) *
				 values(j,point)) *
				weights[point] *
				density_values[point];

      Vector<double> tmp (dofs_per_cell);
				       // this is the right hand side of the
				       // first equation
				       // for the theta scheme:
				       //    rhs1 := Mu^0 + kMv^0
				       //           -(1-theta)theta k^2 Au^0
      Vector<double> rhs1 (dofs_per_cell);

				       // this is the part of the right hand side
				       // of the second equation which depends
				       // on the solutions of the previous time
				       // step.
				       // for the theta scheme:
				       //    rhs2 := Mv^0-(1-theta)kA^0
      Vector<double> rhs2 (dofs_per_cell);
	    
      				       // vector of values of the function on the
				       // old grid restricted to one cell
      Vector<double>     old_dof_values_u (dofs_per_cell);
				       // vector of old u and v times the local
				       // mass matrix
      Vector<double> local_M_u (dofs_per_cell);
      Vector<double> local_M_v (dofs_per_cell);
      Vector<double> local_A_u (dofs_per_cell);      
				       // transfer u+k*v. Note that we need
				       // old_dof_values_u again below
      old_cell->get_dof_values (previous_time_level.u, old_dof_values_u);
      cell_matrix.vmult (local_M_u, old_dof_values_u);
      
      old_cell->get_dof_values (previous_time_level.v, tmp);
      cell_matrix.vmult (local_M_v, tmp);

				       // now for the part with the laplace
				       // matrix
      cell_matrix.clear ();
      vector<double> stiffness_values(fe_values.n_quadrature_points);
      parameters.stiffness->value_list (fe_values.get_quadrature_points(),
					    stiffness_values);
      for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	for (unsigned int i=0; i<dofs_per_cell; ++i) 
	  for (unsigned int j=0; j<dofs_per_cell; ++j)
	    cell_matrix(i,j) += (gradients[i][point] *
				 gradients[j][point]) *
				weights[point] *
				stiffness_values[point];
      cell_matrix.vmult (local_A_u, old_dof_values_u);


      rhs1 = local_M_u;
      rhs1.add (time_step, local_M_v);
      rhs1.add ((-time_step*time_step*
		 parameters.theta*
		 (1-parameters.theta)),
		local_A_u);
      rhs2 = local_M_v;
      rhs2.add (-(1-parameters.theta)*
		time_step,
		local_A_u);

				       // transfer into the global
				       // right hand side
      vector<int> new_dof_indices (dofs_per_cell, -1);
      new_cell->get_dof_indices (new_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  right_hand_side1(new_dof_indices[i]) += rhs1(i);
	  right_hand_side2(new_dof_indices[i]) += rhs2(i);
	};
      
      return;
    };
  
				   // only old cell is refined
  if (old_cell->has_children() && !new_cell->has_children())
    {
				       // this is the right hand side of the
				       // first equation
				       // for the theta scheme:
				       //    rhs1 := Mu^0 + kMv^0
				       //           -(1-theta)theta k^2 Au^0
      Vector<double> rhs1 (dofs_per_cell);

				       // this is the part of the right hand side
				       // of the second equation which depends
				       // on the solutions of the previous time
				       // step.
				       // for the theta scheme:
				       //    rhs2 := Mv^0-(1-theta)kA^0
      Vector<double> rhs2 (dofs_per_cell);

				       // collect the contributions of the
				       // child cells (and possibly their
				       // children as well)
      collect_from_children (old_cell, fe_values, rhs1, rhs2);

				       // transfer into the global
				       // right hand side
      vector<int> new_dof_indices (dofs_per_cell);
      new_cell->get_dof_indices (new_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i) 
	{
	  right_hand_side1(new_dof_indices[i]) += rhs1(i);
	  right_hand_side2(new_dof_indices[i]) += rhs2(i);
	};

      return;
    };

  				   // only new cell is refined
  if (!old_cell->has_children() && new_cell->has_children())
    {
				       // vector of values of the function
				       // on the old grid restricted to
				       // the large (old) cell
      Vector<double>  old_dof_values_u (dofs_per_cell);
      Vector<double>  old_dof_values_v (dofs_per_cell);
      old_cell->get_dof_values (previous_time_level.u, old_dof_values_u);
      old_cell->get_dof_values (previous_time_level.v, old_dof_values_v);

				       // distribute the contribution of the
				       // large old cell to the children on
				       // the new cell
      distribute_to_children (new_cell, fe_values,
			      old_dof_values_u, old_dof_values_v,
			      right_hand_side1, right_hand_side2);

      return;
    };

  Assert (false, ExcInternalError());
};


     
  
template <int dim>
unsigned int
TimeStep_Primal<dim>::collect_from_children (const DoFHandler<dim>::cell_iterator &old_cell,
					     FEValues<dim>  &fe_values,
					     Vector<double> &rhs1,
					     Vector<double> &rhs2) const {
				   // maximal difference of levels between the
				   // cell to which we write and the cells from
				   // which we read. Default is one, but this is
				   // increased with each level of recursion
  unsigned int level_difference = 1;  
  
				   // select the TimeStep_Wave part in the
				   // TimeStep_Primal branch
  const TimeStep_Primal<dim> &previous_time_level
    = static_cast<const TimeStepBase_Wave<dim>*>(previous_timestep)->get_timestep_primal();

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const double time_step = get_backward_timestep();

  
  FullMatrix<double>   cell_matrix (dofs_per_cell, dofs_per_cell);
  
				   // these will hold the values of the
				   // solution on the old grid, i.e. on
				   // the small cells
  Vector<double>  local_old_dof_values_u (dofs_per_cell);
  Vector<double>  local_old_dof_values_v (dofs_per_cell);
  
				   // same for the contributions to the
				   // right hand sides of the projection
  Vector<double>  local_M_u (dofs_per_cell);
  Vector<double>  local_M_v (dofs_per_cell);
  Vector<double>  local_A_u (dofs_per_cell);
				   // this is the right hand side of the
				   // first equation
				   // for the theta scheme:
				   //    rhs1 := Mu^0 + kMv^0
				   //           -(1-theta)theta k^2 Au^0
  Vector<double> child_rhs1 (dofs_per_cell);

				   // this is the part of the right hand side
				   // of the second equation which depends
				   // on the solutions of the previous time
				   // step.
				   // for the theta scheme:
				   //    rhs2 := Mv^0-(1-theta)kA^0
  Vector<double> child_rhs2 (dofs_per_cell);

  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c) 
    {
      const DoFHandler<dim>::cell_iterator old_child = old_cell->child(c);

      child_rhs1.clear ();
      child_rhs2.clear ();
      
				       // if this child is further subdivided:
				       // collect the contributions of the
				       // children
      if (old_child->has_children())
	{
	  const unsigned int l = collect_from_children (old_child, fe_values,
							child_rhs1, child_rhs2);
	  level_difference = max (l+1, level_difference);
	}
      else
	{
	  fe_values.reinit (old_child);
	  const FullMatrix<double>             &values    = fe_values.get_shape_values ();
	  const vector<vector<Tensor<1,dim> > >&gradients = fe_values.get_shape_grads ();
	  const vector<double>                 &weights   = fe_values.get_JxW_values ();

					   // get solutions restricted to small
					   // cell
	  old_child->get_dof_values (previous_time_level.u, local_old_dof_values_u);
	  old_child->get_dof_values (previous_time_level.v, local_old_dof_values_v);
      
					   // compute M*(u+kv) on the small cell
	  cell_matrix.clear ();
	  vector<double> density_values(fe_values.n_quadrature_points);
	  parameters.density->value_list (fe_values.get_quadrature_points(),
					  density_values);
	  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		cell_matrix(i,j) += (values(i,point) *
				     values(j,point)) *
				    weights[point] *
				    density_values[point];

	  cell_matrix.vmult (local_M_u, local_old_dof_values_u);
	  cell_matrix.vmult (local_M_v, local_old_dof_values_v);
      
					   // now for the part with the laplace
					   // matrix
	  cell_matrix.clear ();

	  vector<double> stiffness_values(fe_values.n_quadrature_points);
	  parameters.stiffness->value_list (fe_values.get_quadrature_points(),
					    stiffness_values);
	  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		cell_matrix(i,j) += (gradients[i][point] *
				     gradients[j][point]) *
				    weights[point] *
				    stiffness_values[point];
	  cell_matrix.vmult (local_A_u, local_old_dof_values_u);
      
	  child_rhs1 = local_M_u;
	  child_rhs1.add (time_step, local_M_v);
	  child_rhs1.add ((-time_step*time_step*
			   parameters.theta*
			   (1-parameters.theta)),
			  local_A_u);
	  child_rhs2 = local_M_v;
	  child_rhs2.add (-(1-parameters.theta)*
			  time_step,
			  local_A_u);
	};
      
				       // transfer the contribution of this
				       // child cell to its parent cell
				       // (#true# means: add up)
      fe.prolongate(c).Tvmult (rhs1, child_rhs1, true);
      fe.prolongate(c).Tvmult (rhs2, child_rhs2, true);
    };

  return level_difference;
};




template <int dim>
unsigned int
TimeStep_Primal<dim>::distribute_to_children (const DoFHandler<dim>::cell_iterator &new_cell,
					      FEValues<dim>         &fe_values,
					      const Vector<double>  &old_dof_values_u,
					      const Vector<double>  &old_dof_values_v,
					      Vector<double>        &right_hand_side1,
					      Vector<double>        &right_hand_side2) {
				   // maximal difference of levels between the
				   // cell to which we write and the cells from
				   // which we read. Default is one, but this is
				   // increased with each level of recursion
  unsigned int level_difference = 1;  
  
  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  const double time_step = get_backward_timestep();

  FullMatrix<double>    cell_matrix(dofs_per_cell, dofs_per_cell);
				   // set up a vector which will hold the
				   // restriction of the old
				   // functions (u,v) to a childcell
  Vector<double> local_old_dof_values_u (dofs_per_cell);
  Vector<double> local_old_dof_values_v (dofs_per_cell);

				   // vector of old u and v times the local
				   // mass matrix (on the small cells
				   // respectively)
  Vector<double> local_M_u (dofs_per_cell);
  Vector<double> local_M_v (dofs_per_cell);
  Vector<double> local_A_u (dofs_per_cell);

				   // this is the right hand side of the
				   // first equation
				   // for the theta scheme:
				   //    rhs1 := Mu^0 + kMv^0
				   //           -(1-theta)theta k^2 Au^0
  Vector<double> rhs1 (dofs_per_cell);

				   // this is the part of the right hand side
				   // of the second equation which depends
				   // on the solutions of the previous time
				   // step.
				   // for the theta scheme:
				   //    rhs2 := Mv^0-(1-theta)kA^0
  Vector<double> rhs2 (dofs_per_cell);
	    
				   // indices of the dofs of a cell on
				   // the new grid
  vector<int> new_dof_indices (dofs_per_cell, -1);

      
				       // loop over the child cells
  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c) 
    {
      const DoFHandler<dim>::cell_iterator new_child = new_cell->child(c);

				       // get u and v on the childcells
      fe.prolongate(c).vmult (local_old_dof_values_u,
			      old_dof_values_u);
      fe.prolongate(c).vmult (local_old_dof_values_v,
			      old_dof_values_v);

      if (new_child->has_children())
					 // cell on new grid is further refined
					 // distribute data on this local cell
					 // to its children
	{
	  const unsigned int l = distribute_to_children (new_child, fe_values,
							 local_old_dof_values_u,
							 local_old_dof_values_v,
							 right_hand_side1,
							 right_hand_side2);
	  level_difference = max (l+1, level_difference);
	}
      else
					 // child is not further refined
					 // -> directly distribute data
	{
	  fe_values.reinit (new_child);
	  const FullMatrix<double>             &values = fe_values.get_shape_values ();
	  const vector<vector<Tensor<1,dim> > >&gradients = fe_values.get_shape_grads ();
	  const vector<double>                 &weights   = fe_values.get_JxW_values ();

					   // transfer u+kv
	  cell_matrix.clear ();
	  vector<double> density_values(fe_values.n_quadrature_points);
	  parameters.density->value_list (fe_values.get_quadrature_points(),
					      density_values);
	  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		cell_matrix(i,j) += (values(i,point) *
				     values(j,point)) *
				    weights[point] *
				    density_values[point];

	  cell_matrix.vmult (local_M_u, local_old_dof_values_u);
	  cell_matrix.vmult (local_M_v, local_old_dof_values_v);

					   // now for the part with the laplace
					   // matrix
	  cell_matrix.clear ();
	  vector<double> stiffness_values(fe_values.n_quadrature_points);
	  parameters.stiffness->value_list (fe_values.get_quadrature_points(),
						stiffness_values);
	  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
	    for (unsigned int i=0; i<dofs_per_cell; ++i) 
	      for (unsigned int j=0; j<dofs_per_cell; ++j)
		cell_matrix(i,j) += (gradients[i][point] *
				     gradients[j][point]) *
				    weights[point] *
				    stiffness_values[point];
	  cell_matrix.vmult (local_A_u, local_old_dof_values_u);

	  rhs1 = local_M_u;
	  rhs1.add (time_step, local_M_v);
	  rhs1.add ((-time_step*time_step*
		     parameters.theta*
		     (1-parameters.theta)),
		    local_A_u);
	  rhs2 = local_M_v;
	  rhs2.add (-(1-parameters.theta)*
		    time_step,
		    local_A_u);
	  
					   // transfer into the global
					   // right hand side
	  new_child->get_dof_indices (new_dof_indices);
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {
	      right_hand_side1(new_dof_indices[i]) += rhs1(i);
	      right_hand_side2(new_dof_indices[i]) += rhs2(i);
	    };
	};
    };

  return level_difference;
};







// explicit instantiations
template class TimeStep_Primal<2>;
/* $Id$ */

#include <lac/vector.h>



void UserMatrix::precondition (Vector<double> &dst,
			       const Vector<double> &src) const {
  switch (preconditioning)
    {
      case jacobi:
	    precondition_Jacobi (dst, src);
	    return;
      case sor:
	    precondition_SOR (dst, src);
	    return;
      case ssor:
	    precondition_SSOR (dst, src);
	    return;
      default:
	    dst = src;
	    return;
    };
};


/* $Id$ */



#include <fe/fe_lib.lagrange.h>
#include <base/quadrature_lib.h>







// static objects

const FEQ1<2>   FEHelper<2>::fe_linear;
const FEQ2<2>   FEHelper<2>::fe_quadratic_sub;
#if 2 < 3
const FEQ3<2>   FEHelper<2>::fe_cubic_sub;
const FEQ4<2>   FEHelper<2>::fe_quartic_sub;
#endif
    
const QGauss2<2> FEHelper<2>::q_gauss_2;
const QGauss3<2> FEHelper<2>::q_gauss_3;
const QGauss4<2> FEHelper<2>::q_gauss_4;
const QGauss5<2> FEHelper<2>::q_gauss_5;
const QGauss6<2> FEHelper<2>::q_gauss_6;
const QGauss7<2> FEHelper<2>::q_gauss_7;
const QGauss8<2> FEHelper<2>::q_gauss_8;

#if 2 > 1
const QGauss2<2-1> FEHelper<2>::q_gauss_2_face;
const QGauss3<2-1> FEHelper<2>::q_gauss_3_face;
const QGauss4<2-1> FEHelper<2>::q_gauss_4_face;
const QGauss5<2-1> FEHelper<2>::q_gauss_5_face;
const QGauss6<2-1> FEHelper<2>::q_gauss_6_face;
const QGauss7<2-1> FEHelper<2>::q_gauss_7_face;
const QGauss8<2-1> FEHelper<2>::q_gauss_8_face;
#endif


template <int dim>
const FiniteElement<dim> & FEHelper<dim>::get_fe (const string &name) {
  if (name=="linear")
    return fe_linear;
  else
    if (name=="quadratic")
      return fe_quadratic_sub;
#if 2 < 3
    else
      if (name=="cubic")
	return fe_cubic_sub;
      else
	if (name=="quartic")
	  return fe_quartic_sub;
#endif
  
  Assert (false, ExcInternalError());

  return fe_linear;
};



template <int dim>
const Quadrature<dim> &FEHelper<dim>::get_quadrature (const string &name) {
  if (name=="linear")
    return q_gauss_2;
  else
    if (name=="quadratic")
      return q_gauss_3;
#if 2 < 3
    else
      if (name=="cubic")
	return q_gauss_4;
      else
	if (name=="quartic")
	  return q_gauss_5;
#endif
  
  Assert (false, ExcInternalError());

  return q_gauss_2;
};



template <>
const Quadrature<0> &FEHelper<1>::get_quadrature_face (const string &) {
  static const Quadrature<0> dummy_quadrature(1);
  return dummy_quadrature;
};



template <int dim>
const Quadrature<dim-1> &FEHelper<dim>::get_quadrature_face (const string &name) {
  if (name=="linear")
    return q_gauss_2_face;
  else
    if (name=="quadratic")
      return q_gauss_3_face;
#if 2 < 3
    else
      if (name=="cubic")
	return q_gauss_4_face;
      else
	if (name=="quartic")
	  return q_gauss_5_face;
#endif
  
  Assert (false, ExcInternalError());

  return q_gauss_2_face;
};



string int_to_string (const unsigned int i, const unsigned int digits) {
  string s;
  switch (digits) 
    {
      case 4:
	    s += '0' + i/1000;
      case 3:
	    s += '0' + (i%1000)/100;
      case 2:
	    s += '0' + (i%100)/10;
      case 1:
	    s += '0' + i%10;
	    break;
      default:
	    s += "invalid digits information";
    };
  return s;
};





// explicit instantiations
template class FEHelper<2>;


/* $Id$ */



#include <base/logstream.h>
#include <lac/sparse_matrix.h>
#include <lac/vector.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>




template <int dim>
WaveProblem<dim>::WaveProblem ()
{};



template <int dim>
WaveProblem<dim>::~WaveProblem ()
{};



template <int dim>
void WaveProblem<dim>::declare_parameters (ParameterHandler &prm)
{
  parameters.declare_parameters (prm);
};



template <int dim>
void WaveProblem<dim>::parse_parameters (ParameterHandler &prm)
{
  parameters.parse_parameters (prm);
};



template <int dim>
void WaveProblem<dim>::create_new (const unsigned int)
{
  parameters.delete_parameters ();
};



template <int dim>
void WaveProblem<dim>::run (ParameterHandler &prm) 
{
  parse_parameters (prm);
//  prm.print_parameters (logfile, Text);


				   ////////////////////////////////
				   // Set up the time step objects
  TimestepManager<dim> timestep_manager (parameters);
  if (true) {
				     // push back initial level
    timestep_manager.add_timestep (new TimeStep<dim>(0, parameters));
    double       time = 0;
    unsigned int step_no = 0;
    double       local_time_step;
    
    while (time<parameters.end_time)
      {
	++step_no;
	
					 // if on last time step
					 // allow last time step to
					 // be at most 10% longer than
					 // initially wanted
	if (time+parameters.time_step*1.1 >= parameters.end_time)
	  local_time_step = parameters.end_time-time;
	else
					   // equilibrate time step size
					   // of the two last time steps
	  if (time+2*parameters.time_step >= parameters.end_time)
	    local_time_step = (parameters.end_time-time)/2;
	  else
					     // regular time step
	    local_time_step = parameters.time_step;
	
	time += local_time_step;
	
	timestep_manager.add_timestep (new TimeStep<dim>(time, parameters));
      };
  };


				   ////////////////////////////////////
				   // actually do the work (or rather:
				   // let the work be done)
  for (unsigned int sweep=0; sweep<parameters.number_of_sweeps; ++sweep)
    timestep_manager.run_sweep (sweep);
};




int main ()
{
				   // no additional output to console
  deallog.attach(logfile);
  deallog.depth_console(0);

  logfile.precision(4);
  
  WaveProblem<2> waves;
  MultipleParameterLoop input_data;

  waves.declare_parameters(input_data);

  try 
    {
      input_data.read_input ("wave-test-3.prm");
    }
  catch (exception &e)
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Exception on input: " << e.what() << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
				       // abort
      return 1;
    };

  try
    {
      input_data.loop (waves);
    }
  catch (exception &e)
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Exception on processing: " << e.what() << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
				       // abort
      return 2;
    }
  catch (...) 
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Unknown exception!" << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
				       // abort
      return 3;
    };
  
  
  return 0;
};



