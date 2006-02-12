/* $Id$ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2005, 2006 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2005, 2006 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */

                                 // @sect3{Include files}

				 // Since this program is only an
				 // adaptation of step-4, there is not
				 // much new stuff in terms of header
				 // files. In deal.II, we usually list
				 // include files in the order
				 // base-lac-grid-dofs-fe-numerics,
				 // followed by C++ standard include
				 // files:
#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

#include <fstream>
#include <iostream>

				 // This is the only significant new
				 // header, namely the one in which
				 // the Raviart-Thomas finite element
				 // is declared:
#include <fe/fe_raviart_thomas.h>

				 // Finally, as a bonus in this
				 // program, we will use a tensorial
				 // coefficient. Since it may have a
				 // spatial dependence, we consider it
				 // a tensor-valued function. The
				 // following include file provides
				 // the ``TensorFunction'' class that
				 // offers such functionality:
#include <base/tensor_function.h>


                                 // @sect3{The ``MixedLaplaceProblem'' class template}

				 // Again, since this is an adaptation
				 // of step-6, the main class is
				 // almost the same as the one in that
				 // tutorial program. In terms of
				 // member functions, the main
				 // differences are that the
				 // constructor takes the degree of
				 // the Raviart-Thomas element as an
				 // argument (and that there is a
				 // corresponding member variable to
				 // store this value) and the addition
				 // of the ``compute_error'' function
				 // in which, no surprise, we will
				 // compute the difference between the
				 // exact and the numerical solution
				 // to determine convergence of our
				 // computations:
template <int dim>
class MixedLaplaceProblem 
{
  public:
    MixedLaplaceProblem (const unsigned int degree);
    void run ();
    
  private:
    void make_grid_and_dofs ();
    void assemble_system ();
    void solve ();
    void compute_errors () const;
    void output_results () const;

    const unsigned int   degree;
    
    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    DoFHandler<dim>      dof_handler;

				     // The second difference is that
				     // the sparsity pattern, the
				     // system matrix, and solution
				     // and right hand side vectors
				     // are now blocked. What this
				     // means and what one can do with
				     // such objects is explained in
				     // the introduction to this
				     // program as well as further
				     // down below when we explain the
				     // linear solvers and
				     // preconditioners for this
				     // problem:
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    BlockVector<double>       solution;
    BlockVector<double>       system_rhs;
};


				 // @sect3{Right hand side, boundary values, and exact solution}

				 // Our next task is to define the
				 // right hand side of our problem
				 // (i.e., the scalar right hand side
				 // for the pressure in the original
				 // Laplace equation), boundary values
				 // for the pressure, as well as a
				 // function that describes both the
				 // pressure and the velocity of the
				 // exact solution for later
				 // computations of the error. Note
				 // that these functions have one,
				 // one, and ``dim+1'' components,
				 // respectively, and that we pass the
				 // number of components down to the
				 // ``Function<dim>'' base class. For
				 // the exact solution, we only
				 // declare the function that actually
				 // returns the entire solution vector
				 // (i.e. all components of it) at
				 // once. Here are the respective
				 // declarations:
template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    RightHandSide () : Function<dim>(1) {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};



template <int dim>
class PressureBoundaryValues : public Function<dim> 
{
  public:
    PressureBoundaryValues () : Function<dim>(1) {};
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
};


template <int dim>
class ExactSolution : public Function<dim> 
{
  public:
    ExactSolution () : Function<dim>(dim+1) {};
    
    virtual void vector_value (const Point<dim> &p, 
			       Vector<double>   &value) const;
};


				 // And then we also have to define
				 // these respective functions, of
				 // course. Given our discussion in
				 // the introduction of how the
				 // solution should look like, the
				 // following computations should be
				 // straightforward:
template <int dim>
double RightHandSide<dim>::value (const Point<dim>  &/*p*/,
				  const unsigned int /*component*/) const 
{
  return 0;
}



template <int dim>
double PressureBoundaryValues<dim>::value (const Point<dim>  &p,
					   const unsigned int /*component*/) const 
{
  const double alpha = 0.3;
  const double beta = 1;
  return -(alpha*p[0]*p[1]*p[1]/2 + beta*p[0] - alpha*p[0]*p[0]*p[0]/6);
}



template <int dim>
void
ExactSolution<dim>::vector_value (const Point<dim> &p,
				  Vector<double>   &values) const 
{
  Assert (values.size() == dim+1,
	  ExcDimensionMismatch (values.size(), dim+1));

  const double alpha = 0.3;
  const double beta = 1;

  values(0) = alpha*p[1]*p[1]/2 + beta - alpha*p[0]*p[0]/2;
  values(1) = alpha*p[0]*p[1];
  values(2) = -(alpha*p[0]*p[1]*p[1]/2 + beta*p[0] - alpha*p[0]*p[0]*p[0]/6);
}



				 // @sect3{The inverse permability tensor}

                                 // In addition to the other equation
                                 // data, we also want to use a
                                 // permeability tensor, or better --
                                 // because this is all that appears
                                 // in the weak form -- the inverse of
                                 // the permeability tensor,
                                 // ``KInverse''. For the purpose of
                                 // verifying the exactness of the
                                 // solution and determining
                                 // convergence orders, this tensor is
                                 // more in the way than helpful. We
                                 // will therefore simply set it to
                                 // the identity matrix.
                                 //
                                 // However, a spatially varying
                                 // permeability tensor is
                                 // indispensable in real-life porous
                                 // media flow simulations, and we
                                 // would like to use the opportunity
                                 // to demonstrate the technique to
                                 // use tensor valued functions.
                                 //
                                 // Possibly unsurprising, deal.II
                                 // also has a base class not only for
                                 // scalar and generally vector-valued
                                 // functions (the ``Function'' base
                                 // class) but also for functions that
                                 // return tensors of fixed dimension
                                 // and rank, the ``TensorFunction''
                                 // template. Here, the function under
                                 // consideration returns a dim-by-dim
                                 // matrix, i.e. a tensor of rank 2
                                 // and dimension ``dim''. We then
                                 // choose the template arguments of
                                 // the base class appropriately.
                                 //
                                 // The interface that the
                                 // ``TensorFunction'' class provides
                                 // is essentially equivalent to the
                                 // ``Function'' class. In particular,
                                 // there exists a ``value_list''
                                 // function that takes a list of
                                 // points at which to evaluate the
                                 // function, and returns the values
                                 // of the function in the second
                                 // argument, a list of tensors:
template <int dim>
class KInverse : public TensorFunction<2,dim>
{
  public:
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<Tensor<2,dim> >    &values) const;
};


                                 // The implementation is less
                                 // interesting. As in previous
                                 // examples, we add a check to the
                                 // beginning of the class to make
                                 // sure that the sizes of input and
                                 // output parameters are the same
                                 // (see step-5 for a discussion of
                                 // this technique). Then we loop over
                                 // all evaluation points, and for
                                 // each one first clear the output
                                 // tensor and then set all its
                                 // diagonal elements to one
                                 // (i.e. fill the tensor with the
                                 // identity matrix):
template <int dim>
void
KInverse<dim>::value_list (const std::vector<Point<dim> > &points,
			      std::vector<Tensor<2,dim> >    &values) const
{
  Assert (points.size() == values.size(),
	  ExcDimensionMismatch (points.size(), values.size()));

  for (unsigned int p=0; p<points.size(); ++p)
    {
      values[p].clear ();

      for (unsigned int d=0; d<dim; ++d)
	values[p][d][d] = 1.;
    }
}


                                 // @sect3{extract_u and friends}

                                 // The next three functions are
                                 // needed for matrix and right hand
                                 // side assembly. They are described
                                 // in detail in the introduction to
                                 // this program, so that we do not
                                 // need to discuss them here again:
template <int dim>
Tensor<1,dim>
extract_u (const FEValuesBase<dim> &fe_values,
	   const unsigned int i,
	   const unsigned int q)
{
  Tensor<1,dim> tmp;

  for (unsigned int d=0; d<dim; ++d)
    tmp[d] = fe_values.shape_value_component (i,q,d);

  return tmp;
}



template <int dim>
double
extract_div_u (const FEValuesBase<dim> &fe_values,
	       const unsigned int i,
	       const unsigned int q)
{
  double divergence = 0;
  for (unsigned int d=0; d<dim; ++d)
    divergence += fe_values.shape_grad_component (i,q,d)[d];

  return divergence;
}


  
template <int dim>
double extract_p (const FEValuesBase<dim> &fe_values,
                  const unsigned int i,
                  const unsigned int q)
{
  return fe_values.shape_value_component (i,q,dim);
}



                                 // @sect3{MixedLaplaceProblem class implementation}

                                 // @sect4{MixedLaplaceProblem::MixedLaplaceProblem}

                                 // In the constructor of this class,
                                 // we first store the value that was
                                 // passed in concerning the degree of
                                 // the finite elements we shall use
                                 // (a degree of zero, for example,
                                 // means to use RT(0) and DG(0)), and
                                 // then construct the vector valued
                                 // element belonging to the space X_h
                                 // described in the introduction. The
                                 // rest of the constructor is as in
                                 // the early tutorial programs.
                                 //
                                 // The only thing worth describing
                                 // here is the constructor call of
                                 // the ``fe'' variable. The
                                 // ``FESystem'' class to which this
                                 // variable belongs has a number of
                                 // different constructors that all
                                 // refer to binding simpler elements
                                 // together into one larger
                                 // element. In the present case, we
                                 // want to couple a single RT(degree)
                                 // element with a single DQ(degree)
                                 // element. The constructor to
                                 // ``FESystem'' that does this
                                 // requires us to specity first the
                                 // first base element (the
                                 // ``FE_RaviartThomas'' object of
                                 // given degree) and then the number
                                 // of copies for this base element,
                                 // and then similarly the kind and
                                 // number of ``FE_DGQ''
                                 // elements. Note that the Raviart
                                 // Thomas element already has ``dim''
                                 // vector components, so that the
                                 // coupled element will have
                                 // ``dim+1'' vector components, the
                                 // first ``dim'' of which correspond
                                 // to the velocity variable whereas the
                                 // last one corresponds to the
                                 // pressure.
                                 //
                                 // It is also worth comparing the way
                                 // we constructed this element from
                                 // its base elements, with the way we
                                 // have done so in step-8: there, we
                                 // have built it as ``fe
                                 // (FE_Q<dim>(1), dim)'', i.e. we
                                 // have simply used ``dim'' copies of
                                 // the ``FE_Q(1)'' element, one copy
                                 // for the displacement in each
                                 // coordinate direction.
template <int dim>
MixedLaplaceProblem<dim>::MixedLaplaceProblem (const unsigned int degree)
		:
		degree (degree),
                fe (FE_RaviartThomas<dim>(degree), 1,
                    FE_DGQ<dim>(degree), 1),
		dof_handler (triangulation)
{}



                                 // @sect4{MixedLaplaceProblem::make_grid_and_dofs}

                                 // This next function starts out with
                                 // well-known functions calls that
                                 // create and refine a mesh, and then
                                 // associate degrees of freedom with
                                 // it:
template <int dim>
void MixedLaplaceProblem<dim>::make_grid_and_dofs ()
{
  GridGenerator::hyper_cube (triangulation, -1, 1);
  triangulation.refine_global (3);
  
  dof_handler.distribute_dofs (fe);

                                   // However, then things become
                                   // different. As mentioned in the
                                   // introduction, we want to
                                   // subdivide the matrix into blocks
                                   // corresponding to the two
                                   // different kinds of variables,
                                   // velocity and pressure. To this end,
                                   // we first have to make sure that
                                   // the indices corresponding to
                                   // velocities and pressures are not
                                   // intermingled: First all velocity
                                   // degrees of freedom, then all
                                   // pressure DoFs. This way, the
                                   // global matrix separates nicely
                                   // into a 2x2 system. To achieve
                                   // this, we have to renumber
                                   // degrees of freedom base on their
                                   // vector component, an operation
                                   // that conveniently is already
                                   // implemented:
  DoFRenumbering::component_wise (dof_handler);

                                   // The next thing is that we want
                                   // to figure out the sizes of these
                                   // blocks, so that we can allocate
                                   // an appropriate amount of
                                   // space. To this end, we call the
                                   // ``DoFTools::count_dofs_per_component''
                                   // function that counts how many
                                   // shape functions are non-zero for
                                   // a particular vector
                                   // component. We have ``dim+1''
                                   // vector components, and we have
                                   // to use the knowledge that for
                                   // Raviart-Thomas elements all
                                   // shape functions are nonzero in
                                   // all components. In other words,
                                   // the number of velocity shape
                                   // functions equals the number of
                                   // overall shape functions that are
                                   // nonzero in the zeroth vector
                                   // component. On the other hand,
                                   // the number of pressure variables
                                   // equals the number of shape
                                   // functions that are nonzero in
                                   // the dim-th component. Let us
                                   // compute these numbers and then
                                   // create some nice output with
                                   // that:
  std::vector<unsigned int> dofs_per_component (dim+1);
  DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);  
  const unsigned int n_u = dofs_per_component[0],
                     n_p = dofs_per_component[dim];

  std::cout << "   Number of active cells: "
	    << triangulation.n_active_cells()
	    << std::endl
	    << "   Total number of cells: "
	    << triangulation.n_cells()
	    << std::endl
            << "   Number of degrees of freedom: "
	    << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << ')'
	    << std::endl;

                                   // The next task is to allocate a
                                   // sparsity pattern for the matrix
                                   // that we will create. The way
                                   // this works is that we first
                                   // obtain a guess for the maximal
                                   // number of nonzero entries per
                                   // row (this could be done more
                                   // efficiently in this case, but we
                                   // only want to solve relatively
                                   // small problems for which this is
                                   // not so important). In the second
                                   // step, we allocate a 2x2 block
                                   // pattern and then reinitialize
                                   // each of the blocks to its
                                   // correct size using the ``n_u''
                                   // and ``n_p'' variables defined
                                   // above that hold the number of
                                   // velocity and pressure
                                   // variables. In this second step,
                                   // we only operate on the
                                   // individual blocks of the
                                   // system. In the third step, we
                                   // therefore have to instruct the
                                   // overlying block system to update
                                   // its knowledge about the sizes of
                                   // the blocks it manages; this
                                   // happens with the
                                   // ``sparsity_pattern.collect_sizes()''
                                   // call:
  const unsigned int
    n_couplings = dof_handler.max_couplings_between_dofs();
  
  sparsity_pattern.reinit (2,2);
  sparsity_pattern.block(0,0).reinit (n_u, n_u, n_couplings);
  sparsity_pattern.block(1,0).reinit (n_p, n_u, n_couplings);
  sparsity_pattern.block(0,1).reinit (n_u, n_p, n_couplings);
  sparsity_pattern.block(1,1).reinit (n_p, n_p, n_couplings);
  sparsity_pattern.collect_sizes();

                                   // Now that the sparsity pattern
                                   // and its blocks have the correct
                                   // sizes, we actually need to
                                   // construct the content of this
                                   // pattern, and as usual compress
                                   // it, before we also initialize a
                                   // block matrix with this block
                                   // sparsity pattern:
  DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
  sparsity_pattern.compress();

  system_matrix.reinit (sparsity_pattern);

                                   // Then we have to resize the
                                   // solution and right hand side
                                   // vectors in exactly the same way:
  solution.reinit (2);
  solution.block(0).reinit (n_u);
  solution.block(1).reinit (n_p);
  solution.collect_sizes ();
  
  system_rhs.reinit (2);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_p);
  system_rhs.collect_sizes ();
}


                                 // @sect4{MixedLaplaceProblem::assemble_system}
                                 // Similarly, the function that
                                 // assembles the linear system has
                                 // mostly been discussed already in
                                 // the introduction to this
                                 // example. At its top, what happens
                                 // are all the usual steps, with the
                                 // addition that we do not only
                                 // allocate quadrature and
                                 // ``FEValues'' objects for the cell
                                 // terms, but also for face
                                 // terms. After that, we define the
                                 // usual abbreviations for variables,
                                 // and the allocate space for the
                                 // local matrix and right hand side
                                 // contributions, and the array that
                                 // holds the global numbers of the
                                 // degrees of freedom local to the
                                 // present cell.
template <int dim>
void MixedLaplaceProblem<dim>::assemble_system () 
{  
  QGauss<dim>   quadrature_formula(degree+2);
  QGauss<dim-1> face_quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula, 
			   update_values    | update_gradients |
                           update_q_points  | update_JxW_values);
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
				    update_values    | update_normal_vectors |
				    update_q_points  | update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.n_quadrature_points;
  const unsigned int   n_face_q_points = face_quadrature_formula.n_quadrature_points;

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
                                   // The next step is to declare
                                   // objects that represent the
                                   // source term, pressure boundary
                                   // value, and coefficient in the
                                   // equation. In addition to these
                                   // objects that represent
                                   // continuous functions, we also
                                   // need arrays to hold their values
                                   // at the quadrature points of
                                   // individual cells (or faces, for
                                   // the boundary values). Note that
                                   // in the case of the coefficient,
                                   // the array has to be one of
                                   // matrices.
  const RightHandSide<dim>          right_hand_side;
  const PressureBoundaryValues<dim> pressure_boundary_values;
  const KInverse<dim>               k_inverse;
  
  std::vector<double> rhs_values (n_q_points);
  std::vector<double> boundary_values (n_face_q_points);
  std::vector<Tensor<2,dim> > k_inverse_values (n_q_points);


                                   // With all this in place, we can
                                   // go on with the loop over all
                                   // cells. The body of this loop has
                                   // been discussed in the
                                   // introduction, and will not be
                                   // commented any further here:
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      right_hand_side.value_list (fe_values.get_quadrature_points(),
                                  rhs_values);
      k_inverse.value_list (fe_values.get_quadrature_points(),
                            k_inverse_values);
      
      for (unsigned int q=0; q<n_q_points; ++q) 
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const Tensor<1,dim> phi_i_u = extract_u (fe_values, i, q);
            const double div_phi_i_u = extract_div_u (fe_values, i, q);
            const double phi_i_p = extract_p (fe_values, i, q);
            
            for (unsigned int j=0; j<dofs_per_cell; ++j)
              {
                const Tensor<1,dim> phi_j_u = extract_u (fe_values, j, q);
                const double div_phi_j_u = extract_div_u (fe_values, j, q);
                const double phi_j_p = extract_p (fe_values, j, q);
                
                local_matrix(i,j) += (phi_i_u * k_inverse_values[q] * phi_j_u
                                      - div_phi_i_u * phi_j_p
                                      - phi_i_p * div_phi_j_u)
                                     * fe_values.JxW(q);
              }

            local_rhs(i) += -phi_i_p *
                            rhs_values[q] *
                            fe_values.JxW(q);
          }

      for (unsigned int face_no=0;
	   face_no<GeometryInfo<dim>::faces_per_cell;
	   ++face_no)
	if (cell->at_boundary(face_no))
	  {
	    fe_face_values.reinit (cell, face_no);
	    
	    pressure_boundary_values
	      .value_list (fe_face_values.get_quadrature_points(),
			   boundary_values);

	    for (unsigned int q=0; q<n_face_q_points; ++q) 
	      for (unsigned int i=0; i<dofs_per_cell; ++i)
		{
		  const Tensor<1,dim>
		    phi_i_u = extract_u (fe_face_values, i, q);

		  local_rhs(i) += -(phi_i_u *
				    fe_face_values.normal_vector(q) *
				    boundary_values[q] *
				    fe_face_values.JxW(q));
		}
	  }

                                       // The final step in the loop
                                       // over all cells is to
                                       // transfer local contributions
                                       // into the global matrix and
                                       // right hand side vector. Note
                                       // that we use exactly the same
                                       // interface as in previous
                                       // examples, although we now
                                       // use block matrices and
                                       // vectors instead of the
                                       // regular ones. In other
                                       // words, to the outside world,
                                       // block objects have the same
                                       // interface as matrices and
                                       // vectors, but they
                                       // additionally allow to access
                                       // individual blocks.
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        for (unsigned int j=0; j<dofs_per_cell; ++j)
          system_matrix.add (local_dof_indices[i],
                             local_dof_indices[j],
                             local_matrix(i,j));
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += local_rhs(i);
    }
}


                                 // @sect3{Linear solvers and preconditioners}

                                 // The linear solvers and
                                 // preconditioners we use in this
                                 // example have been discussed in
                                 // significant detail already in the
                                 // introduction. We will therefore
                                 // not discuss the rationale for
                                 // these classes here any more, but
                                 // rather only comment on
                                 // implementational aspects.

                                 // @sect4{The ``InverseMatrix'' class template}

                                 // The first component of our linear
                                 // solver scheme was the creation of
                                 // a class that acts like the inverse
                                 // of a matrix, i.e. which has a
                                 // ``vmult'' function that multiplies
                                 // a vector with an inverse matrix by
                                 // solving a linear system.
                                 //
                                 // While most of the code below
                                 // should be obvious given the
                                 // purpose of this class, two
                                 // comments are in order. First, the
                                 // class is derived from the
                                 // ``Subscriptor'' class so that we
                                 // can use the ``SmartPointer'' class
                                 // with inverse matrix objects. The
                                 // use of the ``Subscriptor'' class
                                 // has been explained before in
                                 // step-7 and step-20. The present
                                 // class also sits on the receiving
                                 // end of this
                                 // ``Subscriptor''/``SmartPointer''
                                 // pair: it holds its pointer to the
                                 // matrix it is supposed to be the
                                 // inverse of through a
                                 // ``SmartPointer'' to make sure that
                                 // this matrix is not destroyed while
                                 // we still have a pointer to it.
                                 //
                                 // Secondly, we realize that we will
                                 // probably perform many
                                 // matrix-vector products with
                                 // inverse matrix objects. Now, every
                                 // time we do so, we have to call the
                                 // CG solver to solve a linear
                                 // system. To work, the CG solver
                                 // needs to allocate four temporary
                                 // vectors that it will release again
                                 // at the end of its operation. What
                                 // this means is that through
                                 // repeated calls to the ``vmult''
                                 // function of this class we have to
                                 // allocate and release vectors over
                                 // and over again.
                                 //
                                 // The natural question is then:
                                 // Wouldn't it be nice if we could
                                 // avoid this, and allocate vectors
                                 // only once? In fact, deal.II offers
                                 // a way to do exactly this. What all
                                 // the linear solvers do is not to
                                 // allocate memory using ``new'' and
                                 // ``delete'', but rather to allocate
                                 // them from an object derived from
                                 // the ``VectorMemory'' class (see
                                 // the module on Vector memory
                                 // management in the API reference
                                 // manual). By default, the linear
                                 // solvers use a derived class
                                 // ``PrimitiveVectorMemory'' that,
                                 // ever time a vector is requested,
                                 // allocates one using ``new'', and
                                 // calls ``delete'' on it again once
                                 // the solver returns it to the
                                 // ``PrimitiveVectorMemory''
                                 // object. This is the appropriate
                                 // thing to do if we do not
                                 // anticipate that the vectors may be
                                 // reused any time soon.
                                 //
                                 // On the other hand, for the present
                                 // case, we would like to have a
                                 // vector memory object that
                                 // allocates vectors when asked by a
                                 // linear solver, but when the linear
                                 // solver returns the vectors, the
                                 // vector memory object holds on to
                                 // them for later requests by linear
                                 // solvers. The
                                 // ``GrowingVectorMemory'' class does
                                 // exactly this: when asked by a
                                 // linear solver for a vector, it
                                 // first looks whether it has unused
                                 // ones in its pool and if so offers
                                 // this vector. If it doesn't, it
                                 // simply grows its pool. Vectors are
                                 // only returned to the C++ runtime
                                 // memory system once the
                                 // ``GrowingVectorMemory'' object is
                                 // destroyed itself.
                                 //
                                 // What we therefore need to do is
                                 // have the present matrix have an
                                 // object of type
                                 // ``GrowingVectorMemory'' as a
                                 // member variable and use it
                                 // whenever we create a linear solver
                                 // object. There is a slight
                                 // complication here: Since the
                                 // ``vmult'' function is marked as
                                 // ``const'' (it doesn't change the
                                 // state of the object, after all,
                                 // and simply operates on its
                                 // arguments), it can only pass an
                                 // unchanging vector memory object to
                                 // the solvers. The solvers, however,
                                 // do change the state of the vector
                                 // memory object, even though this
                                 // has no impact on the actual state
                                 // of the inverse matrix object. The
                                 // compiler would therefore flag any
                                 // such attempt as an error, if we
                                 // didn't make use of a rarely used
                                 // feature of C++: we mark the
                                 // variable as ``mutable''. What this
                                 // does is to allow us to change a
                                 // member variable even from a
                                 // ``const'' member function.
template <class Matrix>
class InverseMatrix : public Subscriptor
{
  public:
    InverseMatrix (const Matrix &m);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const Matrix> matrix;

    mutable GrowingVectorMemory<> vector_memory;    
};


template <class Matrix>
InverseMatrix<Matrix>::InverseMatrix (const Matrix &m)
                :
                matrix (&m)
{}


                                 // Here now is the function that
                                 // implements multiplication with the
                                 // inverse matrix by calling a CG
                                 // solver. Note how we pass the
                                 // vector memory object discussed
                                 // above to the linear solver. Note
                                 // also that we set the solution
                                 // vector to zero before starting the
                                 // solve, since we do not want to use
                                 // the possible previous and unknown
                                 // content of that variable as
                                 // starting vector for the linear
                                 // solve:
template <class Matrix>
void InverseMatrix<Matrix>::vmult (Vector<double>       &dst,
                                   const Vector<double> &src) const
{
  SolverControl solver_control (src.size(), 1e-8*src.l2_norm());
  SolverCG<> cg (solver_control, vector_memory);

  dst = 0;
  
  cg.solve (*matrix, dst, src, PreconditionIdentity());        
}



class SchurComplement : public Subscriptor
{
  public:
    SchurComplement (const BlockSparseMatrix<double> &A,
                     const InverseMatrix<SparseMatrix<double> > &Minv);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double> > > m_inverse;
    
    mutable Vector<double> tmp1, tmp2;
};


SchurComplement::SchurComplement (const BlockSparseMatrix<double> &A,
                                  const InverseMatrix<SparseMatrix<double> > &Minv)
                :
                system_matrix (&A),
                m_inverse (&Minv),
                tmp1 (A.block(0,0).m()),
                tmp2 (A.block(0,0).m())
{}


void SchurComplement::vmult (Vector<double>       &dst,
                             const Vector<double> &src) const
{
  system_matrix->block(0,1).vmult (tmp1, src);
  m_inverse->vmult (tmp2, tmp1);
  system_matrix->block(1,0).vmult (dst, tmp2);
}



class ApproximateSchurComplement : public Subscriptor
{
  public:
    ApproximateSchurComplement (const BlockSparseMatrix<double> &A);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    
    mutable Vector<double> tmp1, tmp2;
};


ApproximateSchurComplement::ApproximateSchurComplement (const BlockSparseMatrix<double> &A)
                :
                system_matrix (&A),
                tmp1 (A.block(0,0).m()),
                tmp2 (A.block(0,0).m())
{}


void ApproximateSchurComplement::vmult (Vector<double>       &dst,
                                   const Vector<double> &src) const
{
  system_matrix->block(0,1).vmult (tmp1, src);
  system_matrix->block(0,0).precondition_Jacobi (tmp2, tmp1);
  system_matrix->block(1,0).vmult (dst, tmp2);
}





template <int dim>
void MixedLaplaceProblem<dim>::solve () 
{
  const InverseMatrix<SparseMatrix<double> >
    m_inverse (system_matrix.block(0,0));
  Vector<double> tmp (solution.block(0).size());
  
  {
    Vector<double> schur_rhs (solution.block(1).size());

    m_inverse.vmult (tmp, system_rhs.block(0));
    system_matrix.block(1,0).vmult (schur_rhs, tmp);
    schur_rhs -= system_rhs.block(1);

    SchurComplement
      schur_complement (system_matrix, m_inverse);
    
    ApproximateSchurComplement
      approximate_schur_complement (system_matrix);
      
    InverseMatrix<ApproximateSchurComplement>
      preconditioner (approximate_schur_complement);
    
    SolverControl solver_control (system_matrix.block(0,0).m(),
				  1e-6*schur_rhs.l2_norm());
    SolverCG<>    cg (solver_control);

    cg.solve (schur_complement, solution.block(1), schur_rhs,
              preconditioner);
  
    std::cout << "   " << solver_control.last_step()
              << " CG Schur complement iterations needed to obtain convergence."
              << std::endl;
  }
  {
    system_matrix.block(0,1).vmult (tmp, solution.block(1));
    tmp *= -1;
    tmp += system_rhs.block(0);
    
    m_inverse.vmult (solution.block(0), tmp);
  }
}



template <int dim>
void MixedLaplaceProblem<dim>::compute_errors () const
{
  Vector<double> tmp (triangulation.n_active_cells());
  ExactSolution<dim> exact_solution;

                                   // do NOT use QGauss here!
  QTrapez<1> q_trapez;
  QIterated<dim> quadrature (q_trapez, degree+2);
  {
    const ComponentSelectFunction<dim> mask (dim, 1., dim+1);
    VectorTools::integrate_difference (dof_handler, solution, exact_solution,
				       tmp, quadrature,
				       VectorTools::L2_norm,
				       &mask);
  }
  const double p_l2_error = tmp.l2_norm();
  
  double u_l2_error = 0;
  for (unsigned int d=0; d<dim; ++d)
    {
      const ComponentSelectFunction<dim> mask(d, 1., dim+1);
      VectorTools::integrate_difference (dof_handler, solution, exact_solution,
					 tmp, quadrature,
					 VectorTools::L2_norm,
					 &mask);
      u_l2_error = std::sqrt (u_l2_error * u_l2_error +
			      tmp.l2_norm() * tmp.l2_norm());
    }
  
  std::cout << "Errors: ||e_p||_L2 = " << p_l2_error
	    << ",   ||e_u||_L2 = " << u_l2_error
	    << std::endl;
}


template <int dim>
void MixedLaplaceProblem<dim>::output_results () const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector (solution, "solution");

  data_out.build_patches (degree+1);

  std::ofstream output (dim == 2 ?
			"solution-2d.gmv" :
			"solution-3d.gmv");
  data_out.write_gmv (output);
}



template <int dim>
void MixedLaplaceProblem<dim>::run () 
{
  std::cout << "Solving problem in " << dim
            << " space dimensions." << std::endl;
  
  make_grid_and_dofs();
  assemble_system ();
  solve ();
  compute_errors ();
  output_results ();
}

    
                                 // @sect3{The ``main'' function}

				 // The main function we stole from
				 // step-6 instead of step-4. It is
				 // almost equal to the one in step-6
				 // (apart from the changed class
				 // names, of course), the only
				 // exception is that we pass the
				 // degree of the finite element space
				 // to the constructor of the mixed
				 // laplace problem.
int main () 
{
  try
    {
      deallog.depth_console (0);

      MixedLaplaceProblem<2> mixed_laplace_problem(0);
      mixed_laplace_problem.run ();
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
