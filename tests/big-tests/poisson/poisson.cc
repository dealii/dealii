/* $Id$ */


#include <grid/dof.h>
#include <grid/tria.h>
#include <fe/fe_lib.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <lac/dsmatrix.h>
#include <grid/dof_constraints.h>
#include <basic/data_io.h>
#include <numerics/base.h>
#include <numerics/assembler.h>

#include <fstream.h>
#include <cmath>
extern "C" {
#  include <stdlib.h>
}

extern TriaActiveIterator<1,CellAccessor<1> > x;
extern TriaActiveIterator<2,CellAccessor<2> > y;



template <int dim>
class PoissonEquation :  public Equation<dim> {
  public:
    PoissonEquation () :
		    Equation<dim>(1) {};

    virtual void assemble (dFMatrix            &cell_matrix,
			   vector<dVector>     &rhs,
			   const FEValues<dim> &fe_values,
			   const Triangulation<dim>::cell_iterator &cell) const;
};



void PoissonEquation<1>::assemble (dFMatrix            &cell_matrix,
				   vector<dVector>     &rhs,
				   const FEValues<1>   &fe_values,
				   const Triangulation<1>::cell_iterator &) const {
  for (unsigned int point=0; point<fe_values.n_quadrature_points; ++point)
    for (unsigned int i=0; i<fe_values.total_dofs; ++i) 
      {
	for (unsigned int j=0; j<fe_values.total_dofs; ++j)
	  cell_matrix(i,j) = fe_values.shape_grad(i,point) *
			     fe_values.shape_grad(j,point) *
			     fe_values.JxW(point);
	rhs[0] = fe_values.shape_value(i,point) *
		 right_hand_side(fe_values.quadrature_point(point)) *
		 fe_values.JxW(point);
      };
};

  


void main () {
  Triangulation<1> tria;
  DoFHandler<1>    dof(tria);
  FELinear<1>      fe;
  ProblemBase<1>   problem(&tria, &dof);
  Quadrature<1>    q;
  
  tria.create_hypercube();
  dof.distribute_dofs (fe);
  problem.assemble (q);
};
