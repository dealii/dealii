//----------------------------  bem_integration.cc  ---------------------------
//    $Id: testsuite.html 13373 2006-07-13 13:12:08Z manigrasso $
//    Version: $Name$ 
//
//    Copyright (C) 2009 by the deal.II authors 
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//    Authors: Luca Heltai, Cataldo Manigrasso
//
//----------------------------  bem_integration.cc  ---------------------------


// 

#include <fstream>
#include <base/logstream.h>

#include <base/convergence_table.h>
#include <base/quadrature_lib.h>
#include <base/table.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_renumbering.h>
#include <fe/fe_dgp.h>
#include <fe/fe_system.h>
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
#include <lac/sparse_direct.h>
#include <lac/lapack_full_matrix.h>
#include <numerics/data_out.h>
#include <base/smartpointer.h>

#include <cmath>
#include <iostream>
#include <math.h>
#include <string>

using namespace std;
using namespace dealii;

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

template <int dim> 
class BEMProblem 
{
public:
    BEMProblem();
    ~BEMProblem();
    
    /** Starts the Boundary Element Method Computation. */
    void run();
    
    /** Initialize mesh and vector space. */
    void read_domain();

    /** Refine and resize all vectors for the active step. */
    void refine_and_resize();
    
    /** Assemble the two system matrices as well as the system right
     * hands side. */
    void assemble_system();

    /** Solve the system. */
    void solve_system();

    /** Output results for the given cycle. */
    void output_results(unsigned int cycle);
    
private:
    /** The boundary element method triangulation. */
    Triangulation<dim-1, dim> tria;

    /** The finite element space for the potential. */
    FE_DGP<2,3> fe;
    
    /** The finite element space for the velocity. */
    FESystem<2,3> fev;
    
    /** The potential degrees of freedom. */
    DoFHandler<dim-1,dim> dh;
    
    /** The velocity degrees of freedom. */
    DoFHandler<dim-1,dim> dhv;

    /** The system matrix. This is I-C. Since the LAPACKFullMatrix
     * does not have a reinit method, we need to work around this a
     * little. */
    SmartPointer<LAPACKFullMatrix<double> > system_matrix;
    
    /** Single layer potential matrix. */
    FullMatrix<double> single_layer_matrix;
    
    /** Normal component of the wind field. */
    Vector<double> Vn;
    
    /** System rhs. */
    Vector<double> rhs;
    
    /** Potential. */
    Vector<double> phi;
    
    /** Something else. */
    Vector<double> vs;
};

template <int dim>
BEMProblem<dim>::BEMProblem() :
    fe(0),
    fev(FE_DGP<dim-1,dim>(0), dim),
    dh(tria),
    dhv(tria)
{}


template <int dim>
BEMProblem<dim>::~BEMProblem() {
    LAPACKFullMatrix<double> * p = system_matrix;
    system_matrix = 0;
    delete p;
}


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
				  update_jacobians |
				  update_quadrature_points );
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
    vector<Point<3> > quad_points = fe_values->get_quadrature_points();
    Point<3> r,a1,a2,n,r_c,n_c;
    r_c = point-cell->center();
    n_c = jacobians[4][2];
    double rn_c = r_c*n_c;
    vector<double> i_S(4);
    vector<double> i_D(4);
    for (unsigned int q_point=0; q_point < 4; ++q_point)
    {
	r = point-quad_points[q_point];
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
    
    double integral = -1./2./numbers::PI
	*atan( ra1*ra2 / (r.norm()* (r*(a12))));

    return integral;

}

template <int dim>
void BEMProblem<dim>::read_domain() {
    // Center of the ball. It is the origin by default.
    Point<dim> p;
    static HyperBallBoundary<dim-1, dim> boundary(p,1.);    

    // Read the sphere from
    GridIn<dim-1, dim> gi;
    gi.attach_triangulation (tria);
    if(dim == 3) {
	std::ifstream in ("coarse_sphere.inp");
	gi.read_ucd (in);
    } else if(dim == 2) {
	std::ifstream in ("coarse_circle.inp");
	gi.read_ucd (in);
    }
    tria.set_boundary(1, boundary);
}



template <int dim>
void BEMProblem<dim>::refine_and_resize() {
    tria.refine_global(1);
    
    dh.distribute_dofs(fe);
    dhv.distribute_dofs(fev);
    
    const unsigned int ndofs =  dh.n_dofs();
    const unsigned int nvdofs =  dhv.n_dofs();
    
    deallog << "Levels: " << tria.n_levels()
	    << ", potential dofs: " << ndofs 
	    << ", velocity dofs: " << nvdofs << endl;
    
    if(system_matrix) {
	LAPACKFullMatrix<double> * p = system_matrix;
	system_matrix = 0;
	delete p;
    }
    
    system_matrix = new LAPACKFullMatrix<double>(ndofs, ndofs);

    Vn.reinit(ndofs);
    rhs.reinit(ndofs);
    phi.reinit(ndofs);
    vs.reinit(nvdofs);
}    

template <int dim>
void BEMProblem<dim>::assemble_system() {
    
    Point<dim> wind;
    wind[0] = 1.;
    
    typename DoFHandler<dim-1,dim>::active_cell_iterator
	celli = dh.begin_active(),
	cellj = dh.begin_active(),
	cellv = dhv.begin_active(),
	endc = dh.end();
    
    QMidpoint<dim-1> midpoint;
    FEValues<dim-1,dim> fe_mid(fe, midpoint,
			       update_values |
			       update_cell_normal_vectors |
			       update_quadrature_points);
    
    vector<unsigned int> dofsi(fe_mid.n_quadrature_points);
    vector<unsigned int> dofsj(fe_mid.n_quadrature_points);
    vector<unsigned int> dofsv(dim);
    
    // Temporary matrix 
    FullMatrix<double> _B(dh.n_dofs(), dh.n_dofs());
    
    // The kernel.
    LaplaceKernelIntegration<dim-1> kernel;
    
    // i runs on points, j runs on cells.
    for(; celli != endc; ++celli, ++cellv) {
	    fe_mid.reinit(celli);
	    Point<dim> ci = celli->center();
	    Point<dim> ni = fe_mid.cell_normal_vector(0);

	    celli->get_dof_indices(dofsi);
	    cellv->get_dof_indices(dofsv);

	    // Vn vector:
	    Vn(dofsi[0]) = wind*ni;
	    vs(dofsv[0]) = ni[0];
	    vs(dofsv[1]) = ni[1];
	    vs(dofsv[2]) = ni[2];
        
	    // Now the two matrices.
	    for(cellj = dh.begin_active(); cellj != endc; ++cellj) {
		vector<double> SD(2,0.); // Single and Double layers.
		cellj->get_dof_indices(dofsj);
		kernel.compute_SD_integral_on_cell(SD,
						   cellj,
						   ci);
		_B (dofsi[0], dofsj[0]) += -SD[0];
		if(dofsi[0] != dofsj[0])
		    (*system_matrix)(dofsi[0], dofsj[0]) += -SD[1];
		if(dofsi[0] == dofsj[0])
		    (*system_matrix)(dofsi[0], dofsj[0]) += 1.;
	    }
	}
    _B.vmult(rhs, Vn);
}

template <int dim>
void BEMProblem<dim>::solve_system() {
    phi.swap(rhs);
    system_matrix->compute_lu_factorization();
    system_matrix->apply_lu_factorization(phi, false);
}


template <int dim>
void BEMProblem<dim>::output_results(unsigned int cycle) {
    
    DataOut<dim-1, DoFHandler<dim-1, dim> > dataout;
    
    dataout.attach_dof_handler(dh);
    dataout.add_data_vector(phi, "phi");
    dataout.build_patches();
    
    char fname[100];
    sprintf(fname, "test_%02d.vtk", cycle);
    std::ofstream file(fname);
    
    dataout.write_vtk(file);
}

template <int dim>
void BEMProblem<dim>::run() {
    read_domain();
    
    const unsigned int number_of_cycles = 4;
    
    for(unsigned int cycle=0; cycle<number_of_cycles; ++cycle) {
	refine_and_resize();
	assemble_system();
	solve_system();
	output_results(cycle);
    }
}


int main () 
{
  try
  {
      deallog.depth_console (3);
      
      BEMProblem<3> laplace_problem;
      laplace_problem.run();
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
