/*----------------------------   poisson.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __poisson_H
#define __poisson_H
/*----------------------------   poisson.h     ---------------------------*/


#include <grid/tria.h>
#include <grid/dof.h>
#include <grid/tria_accessor.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <grid/dof_constraints.h>
#include <basic/data_io.h>
#include <basic/function.h>
#include <basic/parameter_handler.h>
#include <fe/fe_lib.lagrange.h>
#include <fe/quadrature_lib.h>
#include <numerics/base.h>
#include <numerics/assembler.h>
#include <lac/dsmatrix.h>


#include <map.h>
#include <fstream.h>
#include <cmath>
extern "C" {
#  include <stdlib.h>
}






template <int dim>
class PoissonEquation :  public Equation<dim> {
  public:
    PoissonEquation (const Function<dim> &rhs) :
		    Equation<dim>(1),
		    right_hand_side (rhs)  {};

    virtual void assemble (dFMatrix            &cell_matrix,
			   dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const Triangulation<dim>::cell_iterator &cell) const;
    virtual void assemble (dFMatrix            &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const Triangulation<dim>::cell_iterator &cell) const;
    virtual void assemble (dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const Triangulation<dim>::cell_iterator &cell) const;
  protected:
    const Function<dim> &right_hand_side;
};






template <int dim>
class PoissonProblem : public ProblemBase<dim>,
		       public MultipleParameterLoop::UserClass {
  public:
    PoissonProblem ();

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
