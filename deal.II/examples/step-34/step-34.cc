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

#include <base/logstream.h>
#include <base/smartpointer.h>
#include <base/convergence_table.h>
#include <base/quadrature_lib.h>
#include <base/quadrature_selector.h>
#include <base/table.h>
#include <base/parsed_function.h>
#include <base/utilities.h>

#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>
#include <grid/grid_generator.h>
#include <grid/grid_in.h>
#include <grid/grid_out.h>
#include <grid/tria_boundary_lib.h>

#include <fe/fe_dgp.h>
#include <fe/fe_system.h>
#include <fe/fe_tools.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>

#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <lac/full_matrix.h>
#include <lac/precondition.h>
#include <lac/solver_cg.h>
#include <lac/vector.h>
#include <lac/sparse_direct.h>
#include <lac/lapack_full_matrix.h>

#include <numerics/data_out.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

using namespace std;
using namespace dealii;


template <int dim>
class LaplaceKernelIntegration;


template <int dim> 
class BEMProblem 
{
public:
    BEMProblem(const unsigned int degree = 0);
    ~BEMProblem();
    
    // The structure of a boundary element method code is very similar
    // to the structure of a finite element code. By now you should be
    // familiar with reading paramaters from an external file, and
    // with the splitting of the different tasks into different
    // modules. The same applyes to boundary element methods, and we
    // won't comment too much on them, except on the differences.

    void read_parameters(std::string filename);
    
    void run();
    
    void read_domain();

    void refine_and_resize();
    
    // The only really different function that we find here is the
    // assembly routine. We wrote this function in the most possible
    // general way, in order to allow for easy generalization to
    // higher order methods and to different fundamental solutions
    // (e.g., Stokes or Maxwell).
    //
    // The most noticeable difference is the fact that the final
    // matrix is full, and that we have two nested loops on cells
    // instead of the usual one we have in finite element method.
    //
    // The reason for this is that while the basis functions have a
    // compact support, their convolution with the fundamental
    // solution of the laplace equation is global, and needs to be
    // integrated against all other basis functions.
    //
    // The practical consequence is that we have two sets of
    // quadrature formulas, finite element values and temporary
    // elements, one for the inner integration and one for the outer
    // integration. We allow for different quadrature rules to be used
    // in the two integrations to preserve generality and to allow,
    // for example, the use of collocation method (by specifying midpoint 
    // quadrature formula on the outer integration).
    void assemble_system();

    // The only difference in the solution of the system is that the
    // matrix is a LAPACKFullMatrix, which requires a different
    // treatment with respect to what we saw in most of the other
    // examples. Besides from this detail, things proceeds pretty much
    // in the same way as usual.
    void solve_system();
    
    // Once we obtained a solution on the codimension one domain, we
    // want to interpolate it to the rest of the
    // space. This is done by performing again the convolution of the
    // solution with the kernel in the interpolate() function.
    //
    // We would like to plot the velocity variable which is the
    // gradient of the potential solution. The potential solution is
    // only known on the boundary, but we use the convolution with the
    // fundamental solution to interpolate it on a standard dim
    // dimensional continuous finite element space. The plot of the
    // gradient of the extrapolated solution will give us the velocity
    // we want.
    void interpolate();
    
    void output_results(unsigned int cycle);
    
private:
    // The usual deal.II classes can be used for boundary element
    // methods by specifying the "codimension" of the problem. This is
    // done by setting the optional template arguments to
    // Triangulation, FiniteElement and DoFHandler to the dimension of
    // the embedding space. In our case we generate either 1 or 2
    // dimensional meshes embedded in 2 or 3 dimensional spaces.
    //
    // The optional argument by default is equal to the first
    // argument, and produces the usual finite element classes that we
    // saw in all previous examples.

    Triangulation<dim-1, dim>	tria;
    FE_DGP<dim-1,dim>		fe;
    DoFHandler<dim-1,dim>	dh;

    // In BEM methods, the matrix that is generated is
    // dense. Depending on the size of the problem, the final system
    // might be solved by direct LU decomposition, or by iterative
    // methods. Just for the purpose of illustrating the use of the
    // LAPACK classes, we opt for LU decomposition of the final
    // system. Note that this will be very inefficient when the number
    // of dofs grows, since it is of order $n^3$.

    SmartPointer<LAPACKFullMatrix<double> >	system_matrix;    
    Vector<double>				system_rhs;
    Vector<double>				phi;
    
    // The reconstruction of the solution in the entire space is done
    // on a continuous finite element grid of dimension dim. These are
    // the usual ones, and we don't comment any further on them.
    
    Triangulation<dim>	external_tria;
    FE_Q<dim>		external_fe;
    DoFHandler<dim>	external_dh;
    Vector<double>	external_phi;
    
    // The following variables are the one that we fill through a
    // parameter file.
    // The new objects that we use in this example are the
    // ParsedFunction object and the QuadratureSelector object.
    //
    // The ParsedFunction class allows us to easily and quickly define
    // new function objects via parameter files, with custom
    // definitions which can be very
    // complex (see the documentation of that class for all the
    // available options).
    //
    // The QuadratureSelector class allows us to generate quadrature
    // formulas based on an identifying string and on the possible
    // degree of the formula itself. We used this to allow custom
    // selection of quadrature formulas for the inner as well as the
    // outer integration in the calculation of the boundary element
    // matrix.
    //
    // Notice that selecting the midpoint rule as the outer
    // integration formula on uniformly refined meshes is equivalent
    // (up to a scaling factor) to solving the boundary element method
    // via collocation instead of Galerkin technique.
    Functions::ParsedFunction<dim> wind;
    SmartPointer<Quadrature<dim-1> > outer_quadrature_pointer;
    SmartPointer<Quadrature<dim-1> > inner_quadrature_pointer;
    unsigned int n_cycles;
    unsigned int external_refinement;
};



template <int dim>
class LaplaceKernelIntegration
{
public:

    LaplaceKernelIntegration(const FiniteElement<dim-1,dim> &fe);
    ~LaplaceKernelIntegration();

    // This functions computes the integral of the single and double
    // layer potentials on the cell given as a parameter, at the
    // quadrature points @p q. In practice this function produces the objects
    // 
    // \f[
    // \text{dst}_{ik0} := \int_{\text{cell}} G(y - \text[q]_k) rhs(y) dy
    // \f]
    // 
    // and 
    //
    // \f[
    // \text{dst}_{ik1} := \int_{\text{cell}} \frac{\partial
    // G}{\partial \textbf{n}} (y - \text[q]_k) \phi_i(y) dy
    // \f]
    void compute_SD_integral_on_cell(vector<vector<vector<double> > > &dst,
				     typename DoFHandler<dim-1,dim>::active_cell_iterator &cell,
				     const vector<Point<dim> > &q,
				     const Function<dim> &rhs);

    // The following two functions are the actual calculations of the
    // single and double layer potential kernels, with a minus sign in
    // front of them. They are well defined only if the vector $R =
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
    
    SmartPointer<const FiniteElement<dim-1, dim> > fe;
    SmartPointer<FEValues<dim-1,dim> > fe_values;
};



template <int dim>
BEMProblem<dim>::BEMProblem(const unsigned int degree) :
    fe(degree),
    dh(tria),
    external_fe(1),
    external_dh(external_tria),
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
    prm.declare_entry("External refinement", "5", Patterns::Integer());
    
    prm.enter_subsection("Outer quadrature rule");
    prm.declare_entry("Quadrature type", "midpoint", 
		      Patterns::Selection(QuadratureSelector<(dim-1)>::get_quadrature_names()));
    prm.declare_entry("Quadrature order", "0", Patterns::Integer());
    prm.leave_subsection();


    prm.enter_subsection("Inner quadrature rule");
    prm.declare_entry("Quadrature type", "gauss", 
		      Patterns::Selection(QuadratureSelector<(dim-1)>::get_quadrature_names()));
    prm.declare_entry("Quadrature order", "2", Patterns::Integer());
    prm.leave_subsection();
    
    prm.enter_subsection("Wind function 2d");
    Functions::ParsedFunction<2>::declare_parameters(prm, 2);
    prm.leave_subsection();

    prm.enter_subsection("Wind function 3d");
    Functions::ParsedFunction<3>::declare_parameters(prm, 3);
    prm.leave_subsection();
    
    prm.read_input(filename);
    
    n_cycles = prm.get_integer("Number of cycles");		      
    external_refinement = prm.get_integer("External refinement");		      
		      
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
    
    prm.enter_subsection(std::string("Wind function ")+
			 Utilities::int_to_string(dim)+std::string("d"));
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
LaplaceKernelIntegration<3>::LaplaceKernelIntegration(const FiniteElement<2,3> &fe) :
    fe(&fe)
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


// The one dimensional singular integration can be calculated
// exploiting QGaussLogR quadrature formula. The quadrature formula
// is constructed in each step, so the constructor is empty.
template <>
LaplaceKernelIntegration<2>::LaplaceKernelIntegration(const FiniteElement<1,2> &fe) :
    fe(&fe)
{}

template <int dim>
LaplaceKernelIntegration<dim>::~LaplaceKernelIntegration() {
    // We delete the pointer. Since this was created via the new
    // operator, we need to destroy it using delete. But delete does
    // not take smart pointers, which implies we need to first remove
    // detach the smart pointer from the fe_values object, and then
    // delete it.
    if(fe_values) {
	FEValues<dim-1,dim> *fp = fe_values;
	fe_values = 0;
	delete fp;
    }
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
    
    // A boundary element method triangulation is basically the same
    // as a (dim-1) triangulation, with the difference that the
    // vertices belong to a (dim) dimensional space.
    //
    // Some of the mesh formats supported in deal.II use by default
    // three dimensional points to describe meshes. These are the
    // formats which are compatible with the boundary element method
    // capabilities of deal.II. In particular we can use either UCD or
    // GMSH formats. In both cases, we have to be particularly careful
    // with
    // the orientation of the mesh, because, unlike in the standard
    // finite element case, no reordering or compatibility check is
    // performed here.
    //
    // All meshes are considered as oriented, because they are
    // embedded in a higher dimensional space. See the documentation
    // of the GridIn and of the Triangulation for further details on
    // the orientation.
    //
    // The other detail that is required for appropriate refinement of
    // the boundary element mesh, is an accurate description of the
    // manifold that the mesh is approximating. We already saw this
    // several times for the boundary of standard finite element
    // meshes, and here the principle and usage is the same, except
    // that the Boundary description class takes an additional
    // template parameter that specifies the embedding space
    // dimension. 
    
    Point<dim> p;
    static HyperBallBoundary<dim-1, dim> boundary(p,1.);    

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
    
    const unsigned int ndofs =  dh.n_dofs();
    
    deallog << "Levels: " << tria.n_levels()
	    << ", potential dofs: " << ndofs <<  endl;
    
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
		if(dim == 3) {
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
		} else {
		    // In the two dimensional case we only need a
		    // QGaussLogR quadrature formula to correctly
		    // integrate the single layer potential.
		    for(unsigned int q_outer=0; q_outer<n_q_points_outer; ++q_outer) {
			QGaussLogR<1> singular_quad(inner_quadrature.size(),
						    outer_quadrature.point(q_outer),
						    1./cellj->measure());
			FEValues<1,2> fe_v_singular(fe, singular_quad, 
						    update_jacobians |
						    update_cell_normal_vectors |
						    update_quadrature_points );
			fe_v_singular.reinit(cellj);
			
			static vector<Vector<double> > singular_cell_wind(singular_quad.size(), 
									  Vector<double>(dim) );
			
			const vector<Point<dim> > &singular_normals = fe_v_singular.get_cell_normal_vectors();
			const vector<Point<dim> > &singular_q_points = fe_v_singular.get_quadrature_points();
			
			wind.vector_value_list(singular_q_points, singular_cell_wind);
			
			for(unsigned int i=0; i<fe.dofs_per_cell; ++i) {
			    for(unsigned int q_inner=0; q_inner<singular_quad.size(); ++q_inner) {
				double normal_wind = 0;
				for(unsigned int d=0; d<dim; ++d)
				    normal_wind += (singular_cell_wind[q_inner](d)*
						    singular_normals[q_inner][d]);
				
				local_rhs(i) -= ( normal_wind *
						  fe_v_singular.JxW(q_inner)		/
						  numbers::PI				*
						  fe_outer.shape_value(i,q_outer)	*
						  fe_outer.JxW(q_outer) );
				
				for(unsigned int j=0; j<fe.dofs_per_cell; ++j) {
				    
				    // When the indices are the same, we
				    // assemble also the mass matrix.
				    local_matrix(i,j) += ( fe_outer.shape_value(i,q_outer)	     * 
							   fe_outer.shape_value(j,q_outer)	     *
							   fe_outer.JxW(q_outer) );
				}
			    }
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


// We assume here that the boundary element domain is contained in the
// box $[-2,2]^{\text{dim}}$, and we extrapolate the actual solution
// inside this box using the convolution with the fundamental solution.
template <int dim>
void BEMProblem<dim>::interpolate() {
    // Generate the mesh, refine it and distribute dofs on it.
    GridGenerator::hyper_cube(external_tria, -2, 2);
    external_tria.refine_global(external_refinement);
    external_dh.distribute_dofs(external_fe);
    external_phi.reinit(external_dh.n_dofs());
    
    typename DoFHandler<dim-1,dim>::active_cell_iterator
	cell = dh.begin_active(),
	endc = dh.end();

    
    Quadrature<dim-1> &quadrature = *inner_quadrature_pointer;
    
    FEValues<dim-1,dim> fe_v(fe, quadrature,
			     update_values |
			     update_cell_normal_vectors |
			     update_quadrature_points |
			     update_JxW_values);
    
    const unsigned int n_q_points = fe_v.n_quadrature_points;
    
    vector<unsigned int> dofs(fe.dofs_per_cell);
    
    vector<double> local_phi(n_q_points);
    vector<Vector<double> > local_wind(n_q_points, Vector<double>(dim) );
    double normal_wind;
    
    LaplaceKernelIntegration<dim> kernel(fe);
    Point<dim> R;


    typename DoFHandler<dim>::active_cell_iterator
	external_cell = external_dh.begin_active(),
	external_endc = external_dh.end();
    
    vector<unsigned int> external_dofs(external_fe.dofs_per_cell);
    vector<bool> dof_is_treated(external_dh.n_dofs(), false);
    

    for(; external_cell != external_endc; ++external_cell) {
	external_cell->get_dof_indices(external_dofs);
	
	for(unsigned int i=0; i<external_fe.dofs_per_cell; ++i)
	    if(dof_is_treated[external_dofs[i]] == false) {
		
		dof_is_treated[external_dofs[i]] = true;
		
		external_phi(external_dofs[i]) = 0;
		
		for(cell = dh.begin_active(); cell != endc; ++cell) {
		    fe_v.reinit(cell);
		    
		    const vector<Point<dim> > &q_points = fe_v.get_quadrature_points();
		    const vector<Point<dim> > &normals = fe_v.get_cell_normal_vectors();
		    
		    cell->get_dof_indices(dofs);
		    fe_v.get_function_values(phi, local_phi);
		    
		    wind.vector_value_list(q_points, local_wind);
		    
		    for(unsigned int q=0; q<n_q_points; ++q) {
			normal_wind = 0;
			for(unsigned int d=0; d<dim; ++d) 
			    normal_wind += normals[q][d]*local_wind[q](d);
		    
			R =  external_cell->vertex(i) - q_points[q];
			
			external_phi(external_dofs[i]) += ( ( - kernel.nS(R)	* 
							      normal_wind 	-
							      //
							      ( kernel.nD(R)	* 
								normals[q] )	*
							      local_phi[q] )	*
							    fe_v.JxW(q) );
		    }
		}
	    }
    }
    DataOut<dim, DoFHandler<dim> > dataout;
    
    dataout.attach_dof_handler(external_dh);
    dataout.add_data_vector(external_phi, "external_phi");
    dataout.build_patches();
    
    std::string filename = Utilities::int_to_string(dim) + "d_external.vtk";
    std::ofstream file(filename.c_str());
    dataout.write_vtk(file);
}


template <int dim>
void BEMProblem<dim>::output_results(unsigned int cycle) {
    
    DataOut<dim-1, DoFHandler<dim-1, dim> > dataout;
    
    dataout.attach_dof_handler(dh);
    dataout.add_data_vector(phi, "phi");
    dataout.build_patches();
    
    std::string filename = ( Utilities::int_to_string(dim) + 
			     "d_boundary_solution_" +
			     Utilities::int_to_string(cycle) +
			     ".vtk" );
    std::ofstream file(filename.c_str());
    
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
    
    interpolate();
}


int main () 
{
  try
  {
      deallog.depth_console (3);
      BEMProblem<2> laplace_problem_2d;
      // BEMProblem<3> laplace_problem_3d;      

      laplace_problem_2d.run();
      // laplace_problem_3d.run();
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
