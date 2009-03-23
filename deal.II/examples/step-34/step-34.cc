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
#include <lac/matrix_lib.h>

#include <numerics/data_out.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

using namespace std;
using namespace dealii;


template <int dim>
class LaplaceKernel;


template <int dim> 
class BEMProblem 
{
public:
    BEMProblem();
    
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
    FE_Q<dim-1,dim>		fe;
    DoFHandler<dim-1,dim>	dh;

    // In BEM methods, the matrix that is generated is
    // dense. Depending on the size of the problem, the final system
    // might be solved by direct LU decomposition, or by iterative
    // methods. Just for the purpose of illustrating the use of the
    // LAPACK classes, we opt for LU decomposition of the final
    // system. Note that this will be very inefficient when the number
    // of dofs grows, since it is of order $n^3$.

    SparsityPattern		sparsity;
    SparseMatrix<double> 	system_matrix;    
    Vector<double>		system_rhs;
    Vector<double>		phi;
    
    // The reconstruction of the solution in the entire space is done
    // on a continuous finite element grid of dimension dim. These are
    // the usual ones, and we don't comment any further on them.
    
    Triangulation<dim>	external_tria;
    FE_Q<dim>		external_fe;
    DoFHandler<dim>	external_dh;
    Vector<double>	external_phi;
    
    // The following variables are the one that we fill through a
    // parameter file.  The new objects that we use in this example
    // are the ParsedFunction object and the QuadratureSelector
    // object.
    //
    // The ParsedFunction class allows us to easily and quickly define
    // new function objects via parameter files, with custom
    // definitions which can be very complex (see the documentation of
    // that class for all the available options).
    //
    // The QuadratureSelector class allows us to generate quadrature
    // formulas based on an identifying string and on the possible
    // degree of the formula itself. We used this to allow custom
    // selection of the quadrature formulas for the inner integration.
    //
    // Notice that the pointer given below for the quadrature rule is
    // only used for non singular integrals. Whenever the integral is
    // singular, then only the degree of the quadrature pointer is
    // used, and the integration is a special one (see the
    // assemble_matrix below for further details).
    //
    // We also define a couple of parameters which are used in case we
    // wanted to extend the solution to the entire domain. 
    Functions::ParsedFunction<dim> wind;
    SmartPointer<Quadrature<dim-1> > quadrature_pointer;
    unsigned int n_cycles;
    unsigned int external_refinement;
    bool extend_solution;
};



template <int dim>
class LaplaceKernel
{
public:
    // The following two functions are the actual calculations of the
    // single and double layer potential kernels, that is G and Grad
    // G. They are well defined only if the vector $R = x-y$ is
    // different from zero.
    // 
    // Whenever the integration is performed with the singularity
    // inside the given cell, then a special quadrature formula is
    // used that allows one to integrate arbitrary functions against a
    // singular weight on the reference cell.
    //
    // In order to do so, it is necessary to provide a
    // "desingularized" single and double layer potentials which can
    // then be integrated on the given cell. When the @p
    // factor_out_singularity parameter is set to true, then the
    // computed kernels do not conatain the singular factor, which is
    // included in the quadrature formulas as a weighting function.
    //
    // Notice that the QGaussLog quadrature formula is made to
    // integrate f(x)ln|x-x0|, but the kernel for two dimensional
    // problems has the opposite sign. This is taken care of by
    // switching the sign of the two dimensional desingularized
    // kernel.
    double single_layer(const Point<dim> &R, 
			bool factor_out_singularity = false);
    Point<dim> double_layer(const Point<dim> &R, 
			    bool factor_out_singularity = false);
};



template <int dim>
BEMProblem<dim>::BEMProblem() :
    fe(1),
    dh(tria),
    external_fe(1),
    external_dh(external_tria),
    wind(dim)
{}

template <int dim> 
void BEMProblem<dim>::read_parameters(std::string filename) {
    ParameterHandler prm;
    
    prm.declare_entry("Number of cycles", "4", Patterns::Integer());
    prm.declare_entry("External refinement", "5", Patterns::Integer());
    prm.declare_entry("Extend solution on the -2,2 box", "false", Patterns::Bool());
    
    prm.enter_subsection("Quadrature rule");
    prm.declare_entry("Quadrature type", "gauss", 
		      Patterns::Selection(QuadratureSelector<(dim-1)>::get_quadrature_names()));
    prm.declare_entry("Quadrature order", "5", Patterns::Integer());
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
    extend_solution = prm.get_bool("Extend solution on the -2,2 box");

    prm.enter_subsection("Quadrature rule");
    static QuadratureSelector<dim-1> quadrature
		      (prm.get("Quadrature type"),
		       prm.get_integer("Quadrature order"));
    prm.leave_subsection();
    
    prm.enter_subsection(std::string("Wind function ")+
			 Utilities::int_to_string(dim)+std::string("d"));
    wind.parse_parameters(prm);
    prm.leave_subsection();

    quadrature_pointer = &quadrature;
}


template <int dim>
double LaplaceKernel<dim>::single_layer(const Point<dim> &R, 
					bool factor_out_singularity) {
    if(factor_out_singularity == true) 
	return (dim == 2 ? -1. : 1.)/(2*(dim-1)*numbers::PI);
    else
	if(dim == 2)
	    return (-std::log(R.norm()) / (2*numbers::PI) );
	else if(dim == 3)
	    return (1./( R.norm()*4*numbers::PI ) );
	else {
	    Assert(false, ExcInternalError());
	    return 0.;
	}
    return 0.;
}
	


template <int dim>
Point<dim> LaplaceKernel<dim>::double_layer(const Point<dim> &R,
					    bool factor_out_singularity) {
    Point<dim> D(R);
    switch(dim) {
    case 2:
	factor_out_singularity ?  D *= 0  : D /=  -2*numbers::PI * R.square();
	break;
    case 3:
	D /= ( -4*numbers::PI * R.square() * 
	       ( factor_out_singularity ? 1. : R.norm() ) );
	break;
    default:
	Assert(false, ExcInternalError());
	break;
    }
    return D;
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
    
    // The matrix is a full matrix. Notwithstanding this fact, the
    // SparseMatrix class coupled with the SparseDirectUMFPACK solver
    // are still faster than Lapack solvers. The drawback is that we
    // need to assemble a full SparsityPattern.
    system_matrix.clear();
    sparsity.reinit(ndofs, ndofs, ndofs);
    for(unsigned int i=0; i<ndofs;++i)
	for(unsigned int j=0; j<ndofs; ++j)
	    sparsity.add(i,j);
    sparsity.compress();
    system_matrix.reinit(sparsity);
    
    system_rhs.reinit(ndofs);
    phi.reinit(ndofs);
}    

template <int dim>
void BEMProblem<dim>::assemble_system() {
    
    typename DoFHandler<dim-1,dim>::active_cell_iterator
	celli = dh.begin_active(),
	cellj = dh.begin_active(),
	endc = dh.end();
    
    // Quadrature formula for the integration of the kernel in non
    // singular cells. This quadrature is selected with the parameter
    // file, and should be quite precise, since the functions we are
    // integrating are not polynomial functions.
    Quadrature<dim-1> &quadrature = *quadrature_pointer;
    
    // We create initially the singular quadratures for the
    // threedimensional problem, since in this case it is only
    // dependent on the reference element. This quadrature is a
    // standard Gauss quadrature formula reparametrized in such a way
    // that allows one to integrate singularities of the kind 1/R
    // centered at one of the vertices. Here we define a vector of
    // four such quadratures that will be used later on.
    vector<QGaussOneOverR<2> > sing_quadratures_3d; 
    for(unsigned int i=0; i<4; ++i)
	sing_quadratures_3d.push_back(QGaussOneOverR<2>(quadrature.size(), i));
	

    FEValues<dim-1,dim> fe_v(fe, quadrature,
			     update_values |
			     update_cell_normal_vectors |
			     update_quadrature_points |
			     update_JxW_values);
    
    const unsigned int n_q_points = fe_v.n_quadrature_points;
    
    vector<unsigned int> dofs_i(fe.dofs_per_cell);
    vector<unsigned int> dofs_j(fe.dofs_per_cell);

    vector<Vector<double> > cell_wind(n_q_points, Vector<double>(dim) );
    double normal_wind;
    
    Vector<double>	local_rhs(fe.dofs_per_cell);
    FullMatrix<double>  local_matrix(fe.dofs_per_cell, fe.dofs_per_cell);
    
    // The kernel.
    LaplaceKernel<dim> kernel;
    
    Point<dim> R;

    // The index i runs on the collocation points, which are the
    // support of the ith basis function, while j runs on inner
    // integration. We perform this check here to ensure that we are
    // not trying to use this code for high order elements. It will
    // only work with Q1 elements, that is, for fe_dofs_per_cell =
    // GeometryInfo<dim>::vertices_per_cell.
    AssertThrow(fe.dofs_per_cell == GeometryInfo<dim-1>::vertices_per_cell,
		ExcDimensionMismatch(fe.dofs_per_cell, GeometryInfo<dim-1>::vertices_per_cell));
    
    for(; celli != endc; ++celli) {
	
	// On the outer cell, we only need to know how to go from
	// local numbering to global numbering. Each degree of freedom
	// is associated with its support point, which is the ith
	// vertex of the cell.
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
	    
	    fe_v.reinit(cellj);
	    cellj->get_dof_indices(dofs_j);
	    
	    const vector<Point<dim> > &q_points = fe_v.get_quadrature_points();
	    const vector<Point<dim> > &normals = fe_v.get_cell_normal_vectors();
	    wind.vector_value_list(q_points, cell_wind);
	    
	    if(is_singular == false) {
		for(unsigned int q=0; q<n_q_points; ++q) {
		    normal_wind = 0;
		    for(unsigned int d=0; d<dim; ++d) 
			normal_wind += normals[q][d]*cell_wind[q](d);
			
		    for(unsigned int i=0; i<fe.dofs_per_cell; ++i) {
			
			// Distance between the external support point
			// and the quadrature point on the internal
			// cell.
			R = celli->vertex(i)-q_points[q];
			    
			local_rhs(i) += ( kernel.single_layer(R)	* 
					  normal_wind			*
					  fe_v.JxW(q) );
				
			for(unsigned int j=0; j<fe.dofs_per_cell; ++j) {
			    
			    local_matrix(i,j) += ( ( kernel.double_layer(R)	* 
						     normals[q] )		*
						   fe_v.shape_value(j,q)	*
						   fe_v.JxW(q)	);
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
		// 
		// In the two dimensional case we perform the
		// integration using a QGaussLogR quadrature formula,
		// which is specifically designed to integrate
		// logarithmic singularities on the unit interval,
		// while in three dimensions we use the
		// QGaussOneOverR, which allows us to integrate 1/R
		// singularities on the vertices of the reference
		// element. Since we don't want to rebuild the two
		// dimensional quadrature formula at each singular
		// integration, we built them outside the loop on the
		// cells, and we only use a pointer to that quadrature
		// here.
		//
		// Notice that in one dimensional integration this is
		// not possible, since we need to know the scaling
		// parameter for the quadrature, which is not known a
		// priori.
		//
		// Dimension independent programming here is a little
		// tricky, but can be achieved via dynamic casting. We
		// check that everything went ok with an assertion at
		// the end of this block. Notice that the dynamic cast
		// will only work when the dimension is the correct
		// one, in which case it is possible to cast a
		// QGaussLogR and QGaussOneOverR to a Quadrature<1>
		// and Quadrature<2> object.
		//
		// In the other cases this won't be called, and even
		// if it was, the dynamic_cast function would just
		// return a null pointer. We check that this is not
		// the case with the Assert at the end.
		//
		// Notice that in two dimensions the singular
		// quadrature rule depends also on the size of the
		// current cell. For this reason, it is necessary to
		// create a new quadrature for each singular
		// integration. Since we create it using the new
		// operator of C++, we also need to destroy it using
		// the dual of new: delete. This is done at the end,
		// and only if dim == 2.
		Quadrature<dim-1> * singular_quadrature;
		for(unsigned int i=0; i<fe.dofs_per_cell; ++i) {
		    if(dim == 2) {
			singular_quadrature = dynamic_cast<Quadrature<dim-1> *>(
			    new QGaussLogR<1>(quadrature.size(),
					      Point<1>((double)i),
					      1./cellj->measure()));
		    } else {
			singular_quadrature = dynamic_cast<Quadrature<dim-1> *>(
			    & sing_quadratures_3d[i]);
		    }
		    
		    Assert(singular_quadrature, ExcInternalError());

		    
		    FEValues<dim-1,dim> fe_v_singular(fe, *singular_quadrature, 
						      update_jacobians |
						      update_values |
						      update_cell_normal_vectors |
						      update_quadrature_points );
		    fe_v_singular.reinit(cellj);
		    
		    static vector<Vector<double> > singular_cell_wind( (*singular_quadrature).size(), 
								       Vector<double>(dim) );
			
		    const vector<Point<dim> > &singular_normals = fe_v_singular.get_cell_normal_vectors();
		    const vector<Point<dim> > &singular_q_points = fe_v_singular.get_quadrature_points();
			
		    wind.vector_value_list(singular_q_points, singular_cell_wind);
		    
		    for(unsigned int q=0; q<singular_quadrature->size(); ++q) {
			R = celli->vertex(i)-singular_q_points[q];
			double normal_wind = 0;
			for(unsigned int d=0; d<dim; ++d)
			    normal_wind += (singular_cell_wind[q](d)*
					    singular_normals[q][d]);
			
			local_rhs(i) += ( kernel.single_layer(R, is_singular) *
					  normal_wind			      *
					  fe_v_singular.JxW(q) );

			for(unsigned int j=0; j<fe.dofs_per_cell; ++j) {
			    local_matrix(i,j) += (( kernel.double_layer(R, is_singular)	 *
						    singular_normals[q])		 *
						  fe_v_singular.shape_value(j,q)	 *
						  fe_v_singular.JxW(q)	);
			}
		    }
		    if(dim==2) delete singular_quadrature;
		}
	    }
	    // Move the local matrix and rhs to the global one.
	    for(unsigned int i=0; i<fe.dofs_per_cell; ++i) {
		system_rhs(dofs_i[i]) += local_rhs(i);
		for(unsigned int j=0; j<fe.dofs_per_cell; ++j) 
		    system_matrix.add(dofs_i[i],dofs_j[j], local_matrix(i,j));
	    }
	}
    }
    // One quick way to compute the matrix of the solid angles, is to
    // use the Neumann matrix itself. It is enough to multiply the
    // matrix with the vector of ones, to get the diagonal matrix of
    // the alpha solid angles.
    Vector<double> ones(dh.n_dofs()), alpha(dh.n_dofs());
    for(unsigned int i=0; i<dh.n_dofs(); ++i) 
	ones(i) = 1.;
    system_matrix.vmult(alpha, ones);
    for(unsigned int i=0; i<dh.n_dofs(); ++i) 
	system_matrix.add(i,i,-alpha(i));
}

template <int dim>
void BEMProblem<dim>::solve_system() {
    SparseDirectUMFPACK LU;
    LU.initialize(system_matrix);
    LU.vmult(phi, system_rhs);
    
    // Since we are solving a purely Neumann problem, the solution is
    // only known up to a constant potential. We filter out the mean
    // value using the MeanValueFilter class. 
    MeanValueFilter mean_filter;
    mean_filter.filter(phi);
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

    
    Quadrature<dim-1> &quadrature = *quadrature_pointer;
    
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
    
    LaplaceKernel<dim> kernel;
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
			
			external_phi(external_dofs[i]) += ( ( kernel.single_layer(R)	* 
							      normal_wind 	-
							      //
							      ( kernel.double_layer(R)	* 
								normals[q] )		*
							      local_phi[q] )		*
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
    if(extend_solution == true)
	interpolate();
}


int main () 
{
  try
  {
      deallog.depth_console (3);
      BEMProblem<2> laplace_problem_2d;
      laplace_problem_2d.run();

      BEMProblem<3> laplace_problem_3d;      
      laplace_problem_3d.run();
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
