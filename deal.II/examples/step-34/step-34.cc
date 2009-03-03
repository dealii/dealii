//----------------------------  step-34.cc  ---------------------------
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
//----------------------------  step-34.cc  ---------------------------


// 

#include <fstream>
#include <base/logstream.h>

#include <base/convergence_table.h>
#include <base/quadrature_lib.h>
#include <base/quadrature_selector.h>
#include <base/table.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_renumbering.h>
#include <base/parsed_function.h>
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

    LaplaceKernelIntegration(FiniteElement<dim-1,dim> &fe);
    ~LaplaceKernelIntegration();
  
    void run();
    
    // This functions computes the integral of the single and double
    // layer potentials on the cell given as a parameter, at the
    // quadrature points @p q. In practice this function produces the objects
    // 
    // \f[
    // \text{dst}_{ik0} := \int_{\text{cell}} G(y - \text[q]_k) \phi_i dy
    // \f]
    // 
    // and 
    //
    // \f[
    // \text{dst}_{ik1} := \int_{\text{cell}} \frac{\partial
    // G}{\partial \textbf{n}} (y - \text[q]_k) \phi_i dy
    // \f]
    void compute_SD_integral_on_cell(vector<vector<vector<double> > > &dst,
				     typename DoFHandler<dim-1,dim>::active_cell_iterator &cell,
				     const vector<Point<dim> > &q,
				     const Function<dim> &rhs);

    // The following two functions are the actual calculations of the
    // single and double layer potential kernels, with a minus sign in
    // front of them. They are well defined only if the vector $r =
    // x-y$ is different from zero.
    double nS(const Point<dim> &R);
    Point<dim> nD(const Point<dim> &R);
    
private:
    // The following two helper functions should only be called when
    // dim=3. If this is not the case, the default implementation is
    // to throw an exception. When the dimension is equal to two, it
    // is possible to compute the singular integrals using the
    // GaussLog quadrature formulas.

    double term_S(const Point<3> &r,
		  const Point<3> &a1,
		  const Point<3> &a2,
		  const Point<3> &n,
		  const double &rn_c) {
	AssertThrow(false, ExcImpossibleInDim());
	return 0;
    };

    double term_D(const Point<3> &r,
		  const Point<3> &a1,
		  const Point<3> &a2) {
	AssertThrow(false, ExcImpossibleInDim());
	return 0;
    };

    SmartPointer<FEValues<dim-1,dim> > fe_values;
};

template <int dim> 
class BEMProblem 
{
public:
    BEMProblem();
    ~BEMProblem();
    
    // Read parameters.
    void read_parameters(std::string filename);
    
    // Starts the Boundary Element Method Computation.
    void run();
    
    // Initialize mesh and vector space. 
    void read_domain();

    // Refine and resize all vectors for the active step. 
    void refine_and_resize();
    
    // Assemble the two system matrices as well as the system right
    // hands side.
    void assemble_system();

    // Solve the system. 
    void solve_system();

    // Output results for the given cycle. 
    void output_results(unsigned int cycle);
    
private:
    // The boundary element method triangulation. 
    Triangulation<dim-1, dim> tria;

    // The finite element spaces for the potential and the velocity. 
    FE_DGP<dim-1,dim> fe;
    FESystem<dim-1,dim> fev;
    
    // Finite element space used to smoothen the potential solution
    // (from piecewise constant to continuous piecewise quadratic)
    FE_Q<dim-1, dim>  fe_q;

    // And the relevant degrees of freedom.
    DoFHandler<dim-1,dim> dh;
    DoFHandler<dim-1,dim> dhv;
    DoFHandler<dim-1,dim> dhq;

    // The system matrix. This is I-C. Since the LAPACKFullMatrix does
    // not have a reinit method, we need to work around this a little.
    SmartPointer<LAPACKFullMatrix<double> > system_matrix;
    
    // The right hand side, the potential and its smoothed version
    Vector<double> system_rhs;
    Vector<double> phi;
    Vector<double> smooth_phi;
    
    // These are the parameters that we read in from a parameter file.
    // In particular we define the wind function and the outer
    // quadrature.  We use a parsed function, for its ease of
    // definition, and the quadrature formula
    Functions::ParsedFunction<dim> wind;
    SmartPointer<Quadrature<dim-1> > outer_quadrature_pointer;
    SmartPointer<Quadrature<dim-1> > inner_quadrature_pointer;
    unsigned int n_cycles;
};

template <int dim>
BEMProblem<dim>::BEMProblem() :
    fe(0),
    fev(FE_DGP<dim-1,dim>(0), dim),
    fe_q(FE_Q<dim-1,dim>(2)),
    dh(tria),
    dhv(tria),
    dhq(tria),
    wind(dim)
{}

template <int dim>
BEMProblem<dim>::~BEMProblem() {
    LAPACKFullMatrix<double> * p = system_matrix;
    system_matrix = 0;
    delete p;
}


template <int dim> 
void BEMProblem<dim>::read_parameters(std::string filename) {
    ParameterHandler prm;
    
    prm.declare_entry("Number of cycles", "4", Patterns::Integer());
    
    prm.enter_subsection("Outer quadrature rule");
    prm.declare_entry("Quadrature type", "midpoint", 
		      Patterns::Selection(QuadratureSelector<(dim-1)>::get_quadrature_names()));
    prm.declare_entry("Quadrature order", "0", Patterns::Integer());
    prm.leave_subsection();


    prm.enter_subsection("Inner quadrature rule");
    prm.declare_entry("Quadrature type", "midpoint", 
		      Patterns::Selection(QuadratureSelector<(dim-1)>::get_quadrature_names()));
    prm.declare_entry("Quadrature order", "0", Patterns::Integer());
    prm.leave_subsection();
    
    prm.enter_subsection("Wind function");
    Functions::ParsedFunction<dim>::declare_parameters(prm, dim);
    prm.leave_subsection();
    
    prm.read_input(filename);
    
    n_cycles = prm.get_integer("Number of cycles");		      
		      
    prm.enter_subsection("Outer quadrature rule");
    static QuadratureSelector<dim-1> outer_quadrature
		      (prm.get("Quadrature type"),
		       prm.get_integer("Quadrature order"));
    prm.leave_subsection();


    prm.enter_subsection("Inner quadrature rule");
    static QuadratureSelector<dim-1> inner_quadrature
		      (prm.get("Quadrature type"),
		       prm.get_integer("Quadrature order"));
    prm.leave_subsection();
    
    prm.enter_subsection("Wind function");
    wind.parse_parameters(prm);
    prm.leave_subsection();

    outer_quadrature_pointer = &outer_quadrature;
    inner_quadrature_pointer = &inner_quadrature;
}



template <int dim>
double LaplaceKernelIntegration<dim>::nS(const Point<dim> &R) {
    if(dim == 2)
	return (-std::log(R.norm()) / numbers::PI);
    else if(dim == 3)
	return (1./(R.norm()*numbers::PI) );
    else {
	Assert(false, ExcInternalError());
    }
    return 0.;
}
	


template <int dim>
Point<dim> LaplaceKernelIntegration<dim>::nD(const Point<dim> &R) {
    Point<dim> D(R);
    if(dim == 2)
	D /= -numbers::PI * R.square();
    else if(dim == 3)
	D /= -2*numbers::PI * R.square() * R.norm();
    else {
	Assert(false, ExcInternalError());
    }
    return D;
}
	


template <>
LaplaceKernelIntegration<3>::LaplaceKernelIntegration(FiniteElement<2,3> &fe)
{
    // In order to perform the two dimensional singular integration on
    // the given cell, we use standard formulas derived by Morino and
    // Chu, as explained in the introduction. In order to do so, we
    // generate a custom quadrature point with the four vertices and
    // the middle point. We won't use the weights, and we set them to
    // 1.

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
				  update_cell_normal_vectors |
				  update_quadrature_points );
}

template <int dim>
LaplaceKernelIntegration<dim>::~LaplaceKernelIntegration() {
    // We delete the pointer. Since this was created via the new
    // operator, we need to destroy it using delete. But delete does
    // not take smart pointers, which implies we need to first remove
    // detach the smart pointer from the fe_values object, and then
    // delete it.
    FEValues<dim-1,dim> *fp = fe_values;
    fe_values = 0;
    delete fp;
}


template <>
double
LaplaceKernelIntegration<3>::term_S (const Point<3> &r,
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

template <>
double
LaplaceKernelIntegration<3>::term_D (const Point<3> &r,
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

template <>
void
LaplaceKernelIntegration<3>::compute_SD_integral_on_cell(vector<vector<vector<double> > > &dstvv,
							 DoFHandler<2,3>::active_cell_iterator &cell,
							 const vector<Point<3> > &q_points, 
							 const Function<3> &rhs)
{
    fe_values->reinit(cell);
    const vector<Tensor<2,3> > &jacobians = fe_values->get_jacobians();
    const vector<Point<3> > &quad_points = fe_values->get_quadrature_points();
    const vector<Point<3> > &normals = fe_values->get_cell_normal_vectors();
    
    static vector<Vector<double> > cell_wind
	( (*fe_values).n_quadrature_points, Vector<double>(3) );
    static vector<double> normal_wind(quad_points.size());
    
    rhs.vector_value_list(quad_points, cell_wind);
    
    for(unsigned int q=0; q<quad_points.size(); ++q) {
	normal_wind[q] = 0;
	for(unsigned int d=0; d<3; ++d)
	    normal_wind[q] += normals[q][d] * cell_wind[q](d);
    }
    Point<3> r,a1,a2,n,r_c,n_c;
    
    Assert(dstvv.size() == fe_values->dofs_per_cell,
	   ExcDimensionMismatch(dstvv.size(), fe_values->dofs_per_cell));
	   
    for(unsigned int i=0; i<fe_values->dofs_per_cell; ++i) {
	vector<vector<double> > & dstv = dstvv[i];
	Assert(dstv.size() == q_points.size(),
	       ExcDimensionMismatch(dstv.size(), q_points.size()));
    
	/* Check only the first size. */
	Assert(dstv[0].size() == 2,
	       ExcDimensionMismatch(dstv[0].size(), 2));
	
    
	n_c = jacobians[4][2];
	
	for(unsigned int outer_q=0; outer_q<q_points.size(); ++outer_q) {
	    const Point<3> &point = q_points[outer_q];
	    vector<double> &dst = dstv[outer_q];
	    r_c = point-cell->center();
	    double rn_c = r_c*n_c;
	    vector<double> i_S(4);
	    vector<double> i_D(4);
	    for (unsigned int inner_q_point=0; inner_q_point < 4; ++inner_q_point)
	    {
		r = point-quad_points[inner_q_point];
		a1 = jacobians[inner_q_point][0];
		a2 = jacobians[inner_q_point][1];
		n =  jacobians[inner_q_point][2];
		i_S[inner_q_point]= term_S(r,a1,a2,n,rn_c) * normal_wind[inner_q_point];
		i_D[inner_q_point]= term_D(r,a1,a2) * fe_values->shape_value(i,inner_q_point);
	    }
	    dst[0] = (i_S[3]-i_S[1]-i_S[2]+i_S[0]);
	    dst[1] = (i_D[3]-i_D[1]-i_D[2]+i_D[0]);
	}
    }
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

    system_rhs.reinit(ndofs);
    phi.reinit(ndofs);
}    

template <int dim>
void BEMProblem<dim>::assemble_system() {
    
    typename DoFHandler<dim-1,dim>::active_cell_iterator
	celli = dh.begin_active(),
	cellj = dh.begin_active(),
	endc = dh.end();
    
    // Outer quadrature rule. If we choose midpoint quadrature rule,
    // then this is a collocation method. If we choose any other
    // Quadrature rule, then this is Galerkin method.
    Quadrature<dim-1> &outer_quadrature = *outer_quadrature_pointer;
    Quadrature<dim-1> &inner_quadrature = *inner_quadrature_pointer;
    
    FEValues<dim-1,dim> fe_outer(fe, outer_quadrature,
				 update_values |
				 update_cell_normal_vectors |
				 update_quadrature_points |
				 update_JxW_values);
    
    FEValues<dim-1,dim> fe_inner(fe, inner_quadrature,
				 update_values |
				 update_cell_normal_vectors |
				 update_quadrature_points |
				 update_JxW_values);
    
    const unsigned int n_q_points_outer = fe_outer.n_quadrature_points;
    const unsigned int n_q_points_inner = fe_inner.n_quadrature_points;
    
    vector<unsigned int> dofs_i(fe.dofs_per_cell);
    vector<unsigned int> dofs_j(fe.dofs_per_cell);

    vector<Vector<double> > inner_cell_wind(n_q_points_inner, Vector<double>(dim) );
    double inner_normal_wind;
    
    Vector<double>	local_rhs(fe.dofs_per_cell);
    FullMatrix<double>  local_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
    
    // The kernel.
    LaplaceKernelIntegration<dim> kernel(fe);
    
    vector<vector<vector<double> > > single_double_layer_potentials
	(fe.dofs_per_cell, vector<vector<double> >
	 (n_q_points_outer, vector<double> (2, 0.) ) ); 
    
    Point<dim> R;

    
    // The index i runs on outer integration, while j runs on inner integration.
    for(; celli != endc; ++celli) {
	fe_outer.reinit(celli);
	
	const vector<Point<dim> > &q_points_outer = fe_outer.get_quadrature_points();
	const vector<Point<dim> > &normals_outer = fe_outer.get_cell_normal_vectors();

	celli->get_dof_indices(dofs_i);
        
	for(cellj = dh.begin_active(); cellj != endc; ++cellj) {

	    // If we are on the same cell, then the integrals we are
	    // performing are singular, and they require a special
	    // treatment, as explained in the introduction.
	    //
	    // In all other cases, standard Gauss quadrature rules can
	    // be used.
	    bool is_singular = (cellj->index() == celli->index());
	    
	    local_rhs = 0;
	    local_matrix = 0;
	    
	    fe_inner.reinit(cellj);
	    cellj->get_dof_indices(dofs_j);
	    
	    const vector<Point<dim> > &q_points_inner = fe_inner.get_quadrature_points();
	    const vector<Point<dim> > &normals_inner = fe_inner.get_cell_normal_vectors();
	    wind.vector_value_list(q_points_inner, inner_cell_wind);
	    
	    if(is_singular == false) {
		for(unsigned int q_inner=0; q_inner<n_q_points_inner; ++q_inner) {
		    inner_normal_wind = 0;
		    for(unsigned int d=0; d<dim; ++d) 
			inner_normal_wind += normals_inner[q_inner][d]*inner_cell_wind[q_inner](d);
		    
		    for(unsigned int q_outer=0; q_outer<n_q_points_outer; ++q_outer) {
		    
			R = q_points_outer[q_outer]-q_points_inner[q_inner];
			
			for(unsigned int i=0; i<fe.dofs_per_cell; ++i) {
			    local_rhs(i) += ( fe_outer.shape_value(i,q_outer)   *
					      fe_outer.JxW(q_outer)		*
					      //
					      kernel.nS(R)			* 
					      inner_normal_wind			*
					      fe_inner.JxW(q_inner) );
				
			    for(unsigned int j=0; j<fe.dofs_per_cell; ++j) {
				
				local_matrix(i,j) += ( fe_outer.shape_value(i,q_outer)  *
						       fe_outer.JxW(q_outer)		*
						       //
						       ( kernel.nD(R)			* 
							 normals_inner[q_inner] )	*
						       fe_inner.shape_value(j,q_inner)  *
						       fe_inner.JxW(q_inner)	);
			    }
			}
		    }
		}
	    } else {
		// Now we treat the more delicate case. If we are
		// here, it means that the cell that runs on the j
		// index and the one that runs on the i index are the
		// same. In this case both the single and the double
		// layer potential are singular, and they require a
		// special treatment, as explained in the
		// introduction. 

		kernel.compute_SD_integral_on_cell(single_double_layer_potentials, 
						   cellj, q_points_outer, wind);
		
		for(unsigned int i=0; i<fe.dofs_per_cell; ++i) {
		    for(unsigned int q_outer=0; q_outer<n_q_points_outer; ++q_outer) {
			local_rhs(i) += ( - single_double_layer_potentials[0][q_outer][0]  * 
					  fe_outer.shape_value(i,q_outer)	           *
					  fe_outer.JxW(q_outer) );
			
			for(unsigned int j=0; j<fe.dofs_per_cell; ++j) {
			    
			    // When the indices are the same, we
			    // assemble also the mass matrix.
			    local_matrix(i,j) += ( fe_outer.shape_value(i,q_outer)	     * 
						   fe_outer.shape_value(j,q_outer)	     *
						   fe_outer.JxW(q_outer) );
			    
			    local_matrix(i,j) += ( -single_double_layer_potentials[j][q_outer][1] * 
						   fe_outer.shape_value(i,q_outer)		     *
						   fe_outer.JxW(q_outer) );
			}
		    }
		}
	    }
	    // Move the local matrix and rhs to the global one.
	    for(unsigned int i=0; i<fe.dofs_per_cell; ++i) {
		system_rhs(dofs_i[i]) += local_rhs(i);
		for(unsigned int j=0; j<fe.dofs_per_cell; ++j) 
		    (*system_matrix)(dofs_i[i],dofs_j[j]) += local_matrix(i,j);
	    }
	}
    }
}

template <int dim>
void BEMProblem<dim>::solve_system() {
    phi.swap(system_rhs);
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
    
    read_parameters("parameters.prm");
    read_domain();
    
    for(unsigned int cycle=0; cycle<n_cycles; ++cycle) {
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
