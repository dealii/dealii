/* $Id$ */


#include <grid/tria.h>
#include <grid/dof.h>
#include <grid/tria_accessor.h>
#include <grid/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary.h>
#include <grid/dof_constraints.h>
#include <basic/data_io.h>
#include <fe/fe_lib.h>
#include <fe/quadrature_lib.h>
#include <numerics/base.h>
#include <numerics/assembler.h>
#include <lac/dsmatrix.h>

#include <fstream.h>
#include <cmath>
extern "C" {
#  include <stdlib.h>
}

extern TriaActiveIterator<1,CellAccessor<1> > x;
extern TriaActiveIterator<2,CellAccessor<2> > y;
extern TriaRawIterator<1,DoFLineAccessor<1,LineAccessor<1> > > z;



template <int dim>
class PoissonEquation :  public Equation<dim> {
  public:
    PoissonEquation () :
		    Equation<dim>(1) {};

    virtual void assemble (dFMatrix            &cell_matrix,
			   vector<dVector>     &rhs,
			   const FEValues<dim> &fe_values,
			   const Triangulation<dim>::cell_iterator &cell) const;
    double right_hand_side (const Point<dim> &) const;
};



template <int dim>
inline
double PoissonEquation<dim>::right_hand_side (const Point<dim> &p) const {
  switch (dim) 
    {
      case 1:
//	    return ((1-4*3.1415926536*3.1415926536) *
//		    cos(2*3.1415926536*p(0)));
	    return p(0)*p(0)*p(0)-3./2.*p(0)*p(0)-6*p(0)+3;
      case 2:
//	    return ((1-3.1415926536*3.1415926536) *
//		    cos(3.1415926536*p(0)) *
//		    cos(3.1415926536*p(1)));
	    return (p(0)*p(0)*p(0)+p(1)*p(1)*p(1)
		    - 3./2.*(p(0)*p(0)+p(1)*p(1))
		    - 6*(p(0)+p(1))
		    + 6);
      default:
	    return 0;
    };
};



void PoissonEquation<1>::assemble (dFMatrix            &cell_matrix,
				   vector<dVector>     &rhs,
				   const FEValues<1>   &fe_values,
				   const Triangulation<1>::cell_iterator &) const {
  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
      {
	for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	  cell_matrix(i,j) += (fe_values.shape_grad(i,point) *
			       fe_values.shape_grad(j,point) +
			       fe_values.shape_value(i,point) *
			       fe_values.shape_value(j,point)) *
			      fe_values.JxW(point);
	rhs[0](i) += fe_values.shape_value(i,point) *
		     right_hand_side(fe_values.quadrature_point(point)) *
		     fe_values.JxW(point);
      };
};



void PoissonEquation<2>::assemble (dFMatrix            &cell_matrix,
				   vector<dVector>     &rhs,
				   const FEValues<2>   &fe_values,
				   const Triangulation<2>::cell_iterator &) const {
  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
      {
	for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	  cell_matrix(i,j) += (fe_values.shape_grad(i,point) *
			       fe_values.shape_grad(j,point) +
			       fe_values.shape_value(i,point) *
			       fe_values.shape_value(j,point)) *
			      fe_values.JxW(point);
	rhs[0](i) += fe_values.shape_value(i,point) *
		     right_hand_side(fe_values.quadrature_point(point)) *
		     fe_values.JxW(point);
      };
};

  


int main () {
  Triangulation<2>   tria;
  DoFHandler<2>      dof(&tria);
  FELinear<2>        fe;
  ProblemBase<2>     problem(&tria, &dof);
  PoissonEquation<2> equation;
  QGauss4<2>         quadrature;

				   
//  HyperBallBoundary<2> boundary(Point<2>(2,3), 4);

  tria.create_hypercube ();
//  tria.create_hyper_ball(Point<2>(2,3),4);
//  tria.set_boundary (&boundary);
  
/*  tria.refine_global (1);
  tria.begin_active()->set_refine_flag();
  tria.execute_refinement ();
  tria.refine_global (3);
*/

  cout << "Making grid..." << endl;
  
  const unsigned int dim=2;
  tria.refine_global (1);
	
  Triangulation<dim>::active_cell_iterator cell, endc;
  for (int i=0; i<8; ++i) 
    {
      int n_levels = tria.n_levels();
      cell = tria.begin_active();
      endc = tria.end();
      
      for (; cell!=endc; ++cell) 
	{
	  double r      = rand()*1.0/RAND_MAX,
		 weight = 1.*
			  (cell->level()*cell->level()) /
			  (n_levels*n_levels);
	  
	  if (r <= 0.5*weight)
	    cell->set_refine_flag ();
	};
      
      tria.execute_refinement ();
    };
  tria.refine_global (1);
  
  cout << "Distributing dofs..." << endl; 
  dof.distribute_dofs (fe);

  cout << "Assembling matrices..." << endl;
  problem.assemble (equation, quadrature, fe);

  cout << "Solving..." << endl;
  problem.solve ();

  cout << "Printing..." << endl;
  DataOut<2> out;
  ofstream gnuplot("gnuplot.out.5");
  problem.fill_data (out); 
  out.write_gnuplot (gnuplot);
  gnuplot.close ();
  
  return 0;
};
