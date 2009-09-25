//----------------------------  bem_integration.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  bem_integration.cc  ---------------------------


// 


#include "../tests.h"
#include <fstream>
#include <base/logstream.h>

// all include files you need here

#include <base/convergence_table.h>
#include <base/quadrature_lib.h>
#include <base/table.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_dgp.h>
//#include <fe/fe_q.h>
#include <fe/fe_tools.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>
#include <grid/grid_generator.h>
#include <grid/grid_in.h>
#include <grid/grid_out.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria_iterator.h>
#include <lac/full_matrix.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/vector.h>
#include <numerics/data_out.h>
#include <base/smartpointer.h>

#include <cmath>
#include <iostream>
#include <math.h>
#include <string>

using namespace std;

ofstream logfile("bem_integration/output");


template <int dim>
class LaplaceKernelIntegration
{
public:

    LaplaceKernelIntegration();
    ~LaplaceKernelIntegration();
  
    void run();
    
    void compute_SD_integral_on_cell(vector<double> &dst,
				     typename DoFHandler<dim,dim+1>::active_cell_iterator &cell,
				     const Point<dim+1> &point);

private:
    double term_S(const Point<3> &r,
		  const Point<3> &a1,
		  const Point<3> &a2,
		  const Point<3> &n,
		  const double &rn_c);

    double term_D(const Point<3> &r,
		  const Point<3> &a1,
		  const Point<3> &a2);

    SmartPointer<FEValues<dim,dim+1> > fe_values;
};

template <>
LaplaceKernelIntegration<2>::LaplaceKernelIntegration()
{
    static FE_DGP<2,3> fe(0);
    vector<Point<2> > qps(5);
    qps[0] = Point<2>(0,0);
    qps[1] = Point<2>(0,1);
    qps[2] = Point<2>(1,0);
    qps[3] = Point<2>(1,1);
    qps[4] = Point<2>(.5,.5);
    vector<double> ws(5,1.);
    static Quadrature<2> quadrature(qps, ws);
    fe_values = new FEValues<2,3>(fe,quadrature,
				  update_values | 
				  update_jacobians );
}

template <int dim>
LaplaceKernelIntegration<dim>::~LaplaceKernelIntegration() {
    FEValues<dim,dim+1> *fp = fe_values;
    fe_values = 0;
    delete fp;
}


template <>
void
LaplaceKernelIntegration<2>::compute_SD_integral_on_cell(vector<double> &dst,
							 DoFHandler<2,3>::active_cell_iterator &cell,
							 const Point<3> &point)
{
    Assert(dst.size() == 2,
	   ExcDimensionMismatch(dst.size(), 2));
    fe_values->reinit(cell);
    vector<Tensor<2,3> > jacobians = fe_values->get_jacobians();
    Point<3> r,a1,a2,n,r_c,n_c;
    r_c = point-cell->center();
    n_c = jacobians[4][2];
    double rn_c = r_c*n_c;
    vector<double> i_S(4);
    vector<double> i_D(4);
    for (unsigned int q_point=0; q_point < 4; ++q_point)
    {
	r = point-cell->vertex(q_point);
	a1 = jacobians[q_point][0];
	a2 = jacobians[q_point][1];
	n =  jacobians[q_point][2];
	i_S[q_point]=term_S(r,a1,a2,n,rn_c);
	i_D[q_point]=term_D(r,a1,a2);
    }
    dst[0] = (i_S[3]-i_S[1]-i_S[2]+i_S[0]);
    dst[1] = (i_D[3]-i_D[1]-i_D[2]+i_D[0]);

}

template <int dim>
double
LaplaceKernelIntegration<dim>::term_S (const Point<3> &r,
				       const Point<3> &a1,
				       const Point<3> &a2,
				       const Point<3> &n,
				       const double &rn_c)
{
    Point<3> ra1, ra2, a12;

    cross_product(ra1,r,a1);
    cross_product(ra2,r,a2);
    cross_product(a12,a1,a2);    
    
    double integral =
	-1./2./numbers::PI
	*(
	    - ra1*n/a1.norm() * asinh( r*a1/ra1.norm() )
	    + ra2*n/a2.norm() * asinh( r*a2/ra2.norm() )
	    + rn_c * atan( ra1*ra2 / (r.norm()* (r*(a12))))
	    );

    return integral;

}

template <int dim>
double
LaplaceKernelIntegration<dim>::term_D (const Point<3> &r,
				       const Point<3> &a1,
				       const Point<3> &a2)    
{
    Point<3> ra1, ra2, a12;

    cross_product(ra1,r,a1);
    cross_product(ra2,r,a2);
    cross_product(a12,a1,a2);    
    
    double integral = 1./2./numbers::PI
	*atan( ra1*ra2 / (r.norm()* (r*(a12))));

    return integral;

}


double
integration(Point<3> point)
{
    Triangulation<2,3> square;
    GridGenerator::hyper_cube<2,3>(square,0,2);
    DoFHandler<2,3> dof_handler(square);
    static const FE_DGP<2,3> fe(0);
    dof_handler.distribute_dofs(fe);

    DoFHandler<2,3>::active_cell_iterator
	cell = dof_handler.begin_active();

    LaplaceKernelIntegration<2> laplace;
    vector<double> integrals(2);
    laplace.compute_SD_integral_on_cell(integrals,
					cell,
					point);
    return integrals[0];
     
}



int main()
{
    
    deallog.attach(logfile);
    deallog.depth_console(0);
    deallog<<std::fixed;
    deallog<<std::setprecision(5);

    Point<3> point(.5, .5, 0);
    double true_result= -3.163145629/numbers::PI;
    deallog<< "Error on  " << point 
	<< " : " << integration(point)-true_result <<endl;

    point= Point<3>(3, 3, 0);
    true_result= -.2306783616;
    deallog<< "Error on  " << point 
	<< " : " << integration(point)-true_result <<endl;

    point= Point<3>(1.5, .5, 0);
    true_result= -1.006860525;
    deallog<< "Error on  " << point 
	<< " : " << integration(point)-true_result <<endl;

    return 0;
}
