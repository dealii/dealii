/*----------------------------   poisson.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __poisson_H
#define __poisson_H
/*----------------------------   poisson.h     ---------------------------*/


#include <grid/tria.h>
#include <dofs/dof_handler.h>
#include <grid/tria_accessor.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <dofs/dof_constraints.h>
#include <numerics/data_out.h>
#include <base/function.h>
#include <base/parameter_handler.h>
#include <fe/fe_lib.lagrange.h>
#include <base/quadrature_lib.h>
#include "../problem_base.h"
#include <numerics/assembler.h>
#include <lac/sparse_matrix.h>


#include <map>
#include <fstream>
#include <cmath>
#include <cstdlib>







template <int dim>
class PoissonEquation :  public Equation<dim> {
  public:
    PoissonEquation (const Function<dim> &rhs) :
		    Equation<dim>(1),
		    right_hand_side (rhs)  {};

    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;
    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;
    virtual void assemble (Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;
  protected:
    const Function<dim> &right_hand_side;
};






template <int dim>
class PoissonProblem : public ProblemBase<dim>,
		       public MultipleParameterLoop::UserClass {
  public:
    PoissonProblem ();
    virtual ~PoissonProblem();
    
    void clear ();
    
    virtual void create_new (const unsigned int run_no);
    virtual void declare_parameters (ParameterHandler &prm);
    virtual void run (ParameterHandler &prm);


    bool make_grid (ParameterHandler &prm);
    void make_zoom_in_grid ();
    void make_random_grid ();

    bool set_right_hand_side (ParameterHandler &prm);
    bool set_boundary_values (ParameterHandler &prm);
    
  protected:
    Triangulation<dim> *tria;
    DoFHandler<dim>    *dof;
    
    Function<dim>      *rhs;
    Function<dim>      *boundary_values;

    Boundary<dim>      *boundary;
};





/*----------------------------   poisson.h     ---------------------------*/
/* end of #ifndef __poisson_H */
#endif
/*----------------------------   poisson.h     ---------------------------*/
