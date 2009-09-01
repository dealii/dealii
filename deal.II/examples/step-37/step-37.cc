/* $Id$ */
/* Author: Martin Kronbichler, Uppsala University, 2009 */

/*    $Id$       */
/*                                                                */
/*    Copyright (C) 2009 by the deal.II authors                   */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */


				 // The include files are mostly similar to
				 // the ones in step-16.
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>

#include <fe/fe_q.h>
#include <fe/fe_values.h>

#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_generator.h>

#include <multigrid/multigrid.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_accessor.h>
#include <multigrid/mg_transfer.h>
#include <multigrid/mg_tools.h>
#include <multigrid/mg_coarse.h>
#include <multigrid/mg_smoother.h>
#include <multigrid/mg_matrix.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>

#include <fstream>
#include <sstream>

using namespace dealii;



				 // @sect3{Equation data.}

				 // We define a variable coefficient
				 // function for the Poisson problem. It is
				 // similar to the function in step-5. As a
				 // difference, we use the formulation
				 // $\frac{1}{0.1 + \|\bf x\|^2}$ instead of
				 // a discontinuous one. It is merely to
				 // demonstrate the possibilities of this
				 // implemenation, rather than being
				 // physically reasonable.
template <int dim>
class Coefficient : public Function<dim> 
{
  public:
    Coefficient ()  : Function<dim>() {}
    
    virtual double value (const Point<dim>   &p,
			  const unsigned int  component = 0) const;
    
    virtual void value_list (const std::vector<Point<dim> > &points,
			     std::vector<double>            &values,
			     const unsigned int              component = 0) const;
};



template <int dim>
double Coefficient<dim>::value (const Point<dim> &p,
				const unsigned int /*component*/) const 
{
  return 1./(0.1+p.square());
}



template <int dim>
void Coefficient<dim>::value_list (const std::vector<Point<dim> > &points,
				   std::vector<double>            &values,
				   const unsigned int              component) const 
{
  Assert (values.size() == points.size(), 
	  ExcDimensionMismatch (values.size(), points.size()));
  Assert (component == 0, 
	  ExcIndexRange (component, 0, 1));

  const unsigned int n_points = points.size();

  for (unsigned int i=0; i<n_points; ++i)
    values[i] = 1./(0.1+points[i].square());
}



				 // @sect3{Matrix-free implementation.}

				 // Next comes the implemenation of the
				 // Matrix-Free class. It provides standard
				 // information we expect for matrices (like
				 // the size of the matrix), and it is able
				 // to perform matrix-vector
				 // multiplications.
				 // 
				 // TODO: Use WorkStream for parallelization
				 // instead of apply_to_subranges, once we
				 // have realized the best way for doing
				 // that.
template <typename number, class Transformation>
class MatrixFree : public Subscriptor
{
public:
  MatrixFree ();

  void reinit (const unsigned int        n_dofs,
	       const unsigned int        n_cells,
	       const FullMatrix<double> &cell_matrix,
	       const unsigned int        n_points_per_cell);

  void clear();

  unsigned int m () const;
  unsigned int n () const;

  ConstraintMatrix & get_constraints ();

  void set_local_dof_indices (const unsigned int               cell_no,
			      const std::vector<unsigned int> &local_dof_indices);

  void set_derivative_data (const unsigned int    cell_no,
			    const unsigned int    quad_point,
			    const Transformation &trans_in);

  template <typename number2>
  void vmult_on_subrange (const unsigned int first_cell,
			  const unsigned int last_cell,
			  Vector<number2> &dst,
			  const Vector<number2> &src) const;

  template <typename number2>
  void vmult (Vector<number2> &dst,
	      const Vector<number2> &src) const;

  template <typename number2>
  void Tvmult (Vector<number2> &dst,
	       const Vector<number2> &src) const;

  template <typename number2>
  void vmult_add (Vector<number2> &dst,
		  const Vector<number2> &src) const;

  template <typename number2>
  void Tvmult_add (Vector<number2> &dst,
		   const Vector<number2> &src) const;

  number el (const unsigned int row, const unsigned int col) const;

  std::size_t memory_consumption () const;

private:
  FullMatrix<number>      small_matrix;
  ConstraintMatrix        constraints;

  Table<2,unsigned int>   indices_local_to_global;
  Table<2,Transformation> derivatives;

  mutable Vector<number>  diagonal_values;
  mutable bool            diagonal_is_calculated;

  struct SmallMatrixData
  {
    unsigned int m;
    unsigned int n;
    unsigned int n_points;
    unsigned int n_comp;
  };
  unsigned int            n_dofs, n_cols, n_cells;
  SmallMatrixData         matrix_data;
};



template <typename number, class Transformation>
MatrixFree<number,Transformation>::MatrixFree ()
:
  Subscriptor()
{}



template <typename number, class Transformation>
void MatrixFree<number,Transformation>::
reinit (const unsigned int        n_dofs_in,
	const unsigned int        n_cells_in,
	const FullMatrix<double> &small_matrix_in,
	const unsigned int        n_points_per_cell)
{
  n_dofs = n_dofs_in;
  n_cells = n_cells_in;
  small_matrix = small_matrix_in;
  matrix_data.m = small_matrix.m();
  matrix_data.n = small_matrix.n();
  matrix_data.n_points = n_points_per_cell;
  matrix_data.n_comp   = small_matrix.n()/matrix_data.n_points;

  Assert(matrix_data.n_comp * n_points_per_cell == small_matrix.n(),
	 ExcInternalError());

  derivatives.reinit (n_cells, n_points_per_cell);
  indices_local_to_global.reinit (n_cells, small_matrix.m());
  diagonal_is_calculated = false;
}



template <typename number, class Transformation>
void
MatrixFree<number,Transformation>::clear ()
{
  n_dofs = 0;
  n_cells = 0;
  small_matrix.reinit(0,0);
  derivatives.reinit (0,0);
  indices_local_to_global.reinit(0,0);
  diagonal_values.reinit (0);
  constraints.clear();
  diagonal_is_calculated = false;
}



template <typename number, class Transformation>
unsigned int
MatrixFree<number,Transformation>::m () const
{
  return n_dofs;
}



template <typename number, class Transformation>
unsigned int
MatrixFree<number,Transformation>::n () const
{
  return n_dofs;
}



template <typename number, class Transformation>
ConstraintMatrix &
MatrixFree<number,Transformation>::get_constraints ()
{
  return constraints;
}



template <typename number, class Transformation>
void MatrixFree<number,Transformation>::
set_local_dof_indices (const unsigned int               cell_no,
		       const std::vector<unsigned int> &local_dof_indices)
{
  Assert (local_dof_indices.size() == matrix_data.m,
	  ExcDimensionMismatch(local_dof_indices.size(),
			       matrix_data.m));
  for (unsigned int i=0; i<matrix_data.m; ++i)
    {
      Assert (local_dof_indices[i] < n_dofs, ExcInternalError());
      indices_local_to_global(cell_no,i) = local_dof_indices[i];
    }
  diagonal_is_calculated = false;
}



template <typename number, class Transformation>
void MatrixFree<number,Transformation>::
set_derivative_data (const unsigned int cell_no,
		     const unsigned int quad_point,
		     const Transformation &trans_in)
{
  derivatives(cell_no,quad_point) = trans_in;
  diagonal_is_calculated = false;
}



template <typename number, class Transformation>
template <typename number2>
void 
MatrixFree<number,Transformation>::
vmult_on_subrange (const unsigned int    first_cell,
		   const unsigned int    last_cell,
		   Vector<number2>       &dst,
		   const Vector<number2> &src) const
{
  FullMatrix<number> solution_cells, solution_points;

  const unsigned int n_chunks = (last_cell-first_cell)/100 + 1;
  const unsigned int chunk_size = 
    (last_cell-first_cell)/n_chunks + ((last_cell-first_cell)%n_chunks>0);

  for (unsigned int k=first_cell; k<last_cell; k+=chunk_size)
    {
      const unsigned int current_chunk_size = 
	k+chunk_size>last_cell ? last_cell-k : chunk_size;

      solution_cells.reinit (current_chunk_size,matrix_data.m, true);
      solution_points.reinit (current_chunk_size,matrix_data.n, true);

      for (unsigned int i=0; i<current_chunk_size; ++i)
	for (unsigned int j=0; j<matrix_data.m; ++j)
	  solution_cells(i,j) = (number)src(indices_local_to_global(i+k,j));

      solution_cells.mmult (solution_points, small_matrix);

      for (unsigned int i=0; i<current_chunk_size; ++i)
	for (unsigned int j=0; j<matrix_data.n_points; ++j)
	  derivatives(i+k,j).transform(&solution_points(i, j*matrix_data.n_comp));

      solution_points.mTmult (solution_cells, small_matrix);

      static Threads::Mutex mutex;
      Threads::Mutex::ScopedLock lock (mutex);
      for (unsigned int i=0; i<current_chunk_size; ++i)
	for (unsigned int j=0; j<matrix_data.m; ++j)
	  dst(indices_local_to_global(i+k,j)) += (number2)solution_cells(i,j);
    }
}



template <typename number, class Transformation>
template <typename number2>
void 
MatrixFree<number,Transformation>::vmult (Vector<number2>       &dst,
					  const Vector<number2> &src) const
{
  dst = 0;
  vmult_add (dst, src);
}



template <typename number, class Transformation>
template <typename number2>
void 
MatrixFree<number,Transformation>::Tvmult (Vector<number2>       &dst,
					   const Vector<number2> &src) const
{
  dst = 0;
  Tvmult_add (dst,src);
}



template <typename number, class Transformation>
template <typename number2>
void 
MatrixFree<number,Transformation>::vmult_add (Vector<number2>       &dst,
					      const Vector<number2> &src) const
{
  Vector<number2> src_copy (src);
  constraints.distribute(src_copy);
  
  vmult_on_subrange (0, n_cells, dst, src_copy);
  constraints.condense (dst);

				 // Need to do this in order to be
				 // consistent even at constrained
				 // dofs. Need to find a better solution in
				 // the future (e.g. by switching to smaller
				 // vectors that do not contain any
				 // constrained entries).
  for (unsigned int i=0; i<n_dofs; ++i)
    if (constraints.is_constrained(i) == true)
      dst(i) = el(i,i) * src(i);
}



template <typename number, class Transformation>
template <typename number2>
void 
MatrixFree<number,Transformation>::Tvmult_add (Vector<number2>       &dst,
					       const Vector<number2> &src) const
{
  vmult_add (dst,src);
}



template <typename number, class Transformation>
number
MatrixFree<number,Transformation>::el (const unsigned int row,
				       const unsigned int col) const
{
  Assert (row == col, ExcNotImplemented());

  if (diagonal_is_calculated == false)
    {
      diagonal_values.reinit (n_dofs);
      std::vector<number> calculation (matrix_data.n_comp);
      for (unsigned int cell=0; cell<n_cells; ++cell)
	for (unsigned int dof=0; dof<matrix_data.m; ++dof)
	  {
	    double diag_value = 0;
	    for (unsigned int j=0; j<matrix_data.n_points; ++j)
	      {
		for (unsigned int d=0; d<matrix_data.n_comp; ++d)
		  calculation[d] = small_matrix(dof,j*matrix_data.n_comp+d);
		derivatives(cell,j).transform(&calculation[0]);
		for (unsigned int d=0; d<matrix_data.n_comp; ++d)
		  diag_value += calculation[d]*small_matrix(dof,j*matrix_data.n_comp+d);
	      }
	    diagonal_values(indices_local_to_global(cell,dof)) += diag_value;
	  }
      diagonal_is_calculated = true;
    }

  return diagonal_values(row);
}



template <typename number, class Transformation>
std::size_t MatrixFree<number,Transformation>::memory_consumption () const
{
  std::size_t glob_size = derivatives.memory_consumption() + 
    indices_local_to_global.memory_consumption() + 
    constraints.memory_consumption() +
    small_matrix.memory_consumption() + sizeof(*this);
  return glob_size;
}



				 // @sect3{Laplace operator.}

				 // This implements the local action of a
				 // Laplace preconditioner.
template <int dim,typename number>
class LaplaceOperator
{
public:
  LaplaceOperator ();

  LaplaceOperator (const Tensor<2,dim> &tensor);

  void transform (number * result) const;

  LaplaceOperator<dim,number>&
  operator = (const Tensor<2,dim> &tensor);

  number transformation[dim*(dim+1)/2];
};

template<int dim,typename number>
LaplaceOperator<dim,number>::LaplaceOperator()
{}

template<int dim,typename number>
LaplaceOperator<dim,number>::LaplaceOperator(const Tensor<2,dim> &tensor)
{
  *this = tensor;
}

template <int dim, typename number>
void LaplaceOperator<dim,number>::transform (number* result) const
{
  number temp_result[dim];
  for (unsigned int d=0; d<dim; ++d)
    temp_result[d] = result[d];
  for (unsigned int d=0; d<dim; ++d)
    {
      number output = transformation[d]*temp_result[d];
      if (dim == 2)
	output += transformation[2]*temp_result[1-d];
      else if (dim == 3)
	{
	  if (d==0)
	    output += transformation[3]*temp_result[1] + transformation[4]*temp_result[2];
	  else if (d==1)
	    output += transformation[3]*temp_result[0] + transformation[5]*temp_result[2];
	  else
	    output += transformation[4]*temp_result[0] + transformation[5]*temp_result[1];
	}
      result[d] = output;
    }
}

template <int dim, typename number>
LaplaceOperator<dim,number>&
LaplaceOperator<dim,number>::operator=(const Tensor<2,dim> &tensor) 
{
  if (dim == 2)
    {
      transformation[0] = tensor[0][0];
      transformation[1] = tensor[1][1];
      transformation[2] = tensor[0][1];
      Assert (std::fabs(tensor[1][0]-tensor[0][1])<1e-15,
	      ExcInternalError());
    }
  else if (dim == 3)
    {
      transformation[0] = tensor[0][0];
      transformation[1] = tensor[1][1];
      transformation[2] = tensor[2][2];
      transformation[3] = tensor[0][1];
      transformation[4] = tensor[0][2];
      transformation[5] = tensor[1][2];
      Assert (std::fabs(tensor[1][0]-tensor[0][1])<1e-15,
	      ExcInternalError());
      Assert (std::fabs(tensor[2][0]-tensor[0][2])<1e-15,
	      ExcInternalError());
      Assert (std::fabs(tensor[2][1]-tensor[1][2])<1e-15,
	      ExcInternalError());
    }
  return *this;
}



				 // @sect3{LaplaceProblem class.}

				 // This class is based on the same class in
				 // step-16. We replaced the
				 // SparseMatrix<double> class by our
				 // matrix-free implementation, which means
				 // that we can skip the sparsity patterns.
template <int dim>
class LaplaceProblem 
{
  public:
    LaplaceProblem (const unsigned int degree);
    void run ();
    
  private:
    void setup_system ();
    void assemble_system ();
    void assemble_multigrid ();
    void solve ();
    void output_results (const unsigned int cycle) const;

    Triangulation<dim>   triangulation;
    FE_Q<dim>            fe;
    MGDoFHandler<dim>    mg_dof_handler;

    MatrixFree<double,LaplaceOperator<dim,double> > system_matrix;
    typedef MatrixFree<float,LaplaceOperator<dim,float> > MatrixFreeType;
    MGLevelObject<MatrixFreeType> mg_matrices;
    FullMatrix<float>             coarse_matrix;

    Vector<double>       solution;
    Vector<double>       system_rhs;
};



template <int dim>
LaplaceProblem<dim>::LaplaceProblem (const unsigned int degree) :
                fe (degree),
		mg_dof_handler (triangulation)
{}



				 // This is the function of step-16 with
				 // relevant changes due to the MatrixFree
				 // class.
template <int dim>
void LaplaceProblem<dim>::setup_system ()
{
  system_matrix.clear();
  mg_matrices.clear();

  mg_dof_handler.distribute_dofs (fe);

  std::cout << "Number of degrees of freedom: "
	    << mg_dof_handler.n_dofs()
	    << std::endl;

  const unsigned int nlevels = triangulation.n_levels();
  mg_matrices.resize(0, nlevels-1);

  QGauss<dim>  quadrature_formula(fe.degree+1);
  FEValues<dim> fe_values2 (fe, quadrature_formula, 
			    update_gradients);
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria, 0, 1);
  fe_values2.reinit (tria.begin());
  FullMatrix<double> data_matrix (fe.dofs_per_cell, 
				  quadrature_formula.size()*dim);
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    {
      for (unsigned int j=0; j<quadrature_formula.size(); ++j)
	{
	  for (unsigned int d=0; d<dim; ++d)
	    data_matrix(i,j*dim+d) = fe_values2.shape_grad(i,j)[d];
	}
    }
  system_matrix.reinit (mg_dof_handler.n_dofs(), triangulation.n_active_cells(),
			data_matrix, quadrature_formula.size());
  VectorTools::interpolate_boundary_values (mg_dof_handler,
					    0,
					    ZeroFunction<dim>(),
					    system_matrix.get_constraints());
  system_matrix.get_constraints().close();
  std::cout.precision(4);
  std::cout << "System matrix memory consumption: " 
	    << (double)system_matrix.memory_consumption()*std::pow(2.,-20.) << " MBytes." 
	    << std::endl;

  solution.reinit (mg_dof_handler.n_dofs());
  system_rhs.reinit (mg_dof_handler.n_dofs());

				 // Initialize the matrices for the
				 // multigrid method on all the levels.
  typename FunctionMap<dim>::type dirichlet_boundary;
  ZeroFunction<dim>               homogeneous_dirichlet_bc (1);
  dirichlet_boundary[0] = &homogeneous_dirichlet_bc;
  std::vector<std::set<unsigned int> > boundary_indices(triangulation.n_levels());
  MGTools::make_boundary_list (mg_dof_handler,
			       dirichlet_boundary,
			       boundary_indices);
  for (unsigned int level=0;level<nlevels;++level)
    {
      mg_matrices[level].reinit(mg_dof_handler.n_dofs(level),
				triangulation.n_cells(level),
				data_matrix,
				quadrature_formula.size());
      std::set<unsigned int>::iterator bc_it = boundary_indices[level].begin();
      for ( ; bc_it != boundary_indices[level].end(); ++bc_it)
	mg_matrices[level].get_constraints().add_line(*bc_it);
      mg_matrices[level].get_constraints().close();
    }
  coarse_matrix.reinit (mg_dof_handler.n_dofs(0),
			mg_dof_handler.n_dofs(0));
}



template <int dim>
void LaplaceProblem<dim>::assemble_system () 
{
  QGauss<dim>  quadrature_formula(fe.degree+1);
  MappingQ<dim> mapping (fe.degree);
  FEValues<dim> fe_values (mapping, fe, quadrature_formula, 
			   update_values   | update_inverse_jacobians |
			   update_gradients |
                           update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  const Coefficient<dim> coefficient;
  std::vector<double>    coefficient_values (n_q_points);

  unsigned int cell_no = 0;

  typename DoFHandler<dim>::active_cell_iterator cell = mg_dof_handler.begin_active(),
						 endc = mg_dof_handler.end();
  for (; cell!=endc; ++cell, ++cell_no)
    {
      cell->get_dof_indices (local_dof_indices);
      fe_values.reinit (cell);
      coefficient.value_list (fe_values.get_quadrature_points(),
			      coefficient_values);

      system_matrix.set_local_dof_indices (cell_no, local_dof_indices);
      for (unsigned int q=0; q<n_q_points; ++q)
	system_matrix.set_derivative_data 
	  (cell_no, q,
	   (transpose(fe_values.inverse_jacobian(q)) * 
	    fe_values.inverse_jacobian(q)) * 
	   fe_values.JxW(q) * coefficient_values[q]);

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  double rhs_val = 0;
	  for (unsigned int q=0; q<n_q_points; ++q)
	    rhs_val += (fe_values.shape_value(i,q) * 1.0 *
			fe_values.JxW(q));
	  system_rhs(local_dof_indices[i]) += rhs_val;
	}
    }
  system_matrix.get_constraints().condense(system_rhs);
}


				 // Here is another assemble
				 // function. The integration core is
				 // the same as above. Only the loop
				 // goes over all existing cells now
				 // and the results must be entered
				 // into the correct matrix.

				 // Since we only do multi-level
				 // preconditioning, no right-hand
				 // side is assembled here.
template <int dim>
void LaplaceProblem<dim>::assemble_multigrid () 
{
  coarse_matrix = 0;
  QGauss<dim>  quadrature_formula(fe.degree+1);
  MappingQ<dim> mapping (fe.degree);
  FEValues<dim> fe_values (mapping, fe, quadrature_formula, 
			   update_gradients  | update_inverse_jacobians |
                           update_quadrature_points | update_JxW_values);

  const unsigned int   dofs_per_cell = fe.dofs_per_cell;
  const unsigned int   n_q_points    = quadrature_formula.size();

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  const Coefficient<dim> coefficient;
  std::vector<double>    coefficient_values (n_q_points);

  std::vector<unsigned int> cell_no(triangulation.n_levels());
  typename MGDoFHandler<dim>::cell_iterator cell = mg_dof_handler.begin(),
					    endc = mg_dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      const unsigned int level = cell->level();
      cell->get_mg_dof_indices (local_dof_indices);
      fe_values.reinit (cell);
      coefficient.value_list (fe_values.get_quadrature_points(),
			      coefficient_values);

      mg_matrices[level].set_local_dof_indices (cell_no[level], 
						local_dof_indices);
      for (unsigned int q=0; q<n_q_points; ++q)
	mg_matrices[level].set_derivative_data 
	  (cell_no[level], q,
	   (transpose(fe_values.inverse_jacobian(q)) * 
	    fe_values.inverse_jacobian(q)) * 
	   fe_values.JxW(q) * coefficient_values[q]);

      ++cell_no[level];
      if (level == 0)
	{
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      {
		double add_value = 0;
		for (unsigned int q=0; q<n_q_points; ++q)
		  add_value += (fe_values.shape_grad(i,q) *
				fe_values.shape_grad(j,q) *
				coefficient_values[q] *
				fe_values.JxW(q));
		coarse_matrix(local_dof_indices[i],
			      local_dof_indices[j]) += add_value;
	      }
	}
    }
  for (unsigned int i=0; i<coarse_matrix.m(); ++i)
    if (mg_matrices[0].get_constraints().is_constrained(i))
      for (unsigned int j=0; j<coarse_matrix.n(); ++j)
	if (i!=j)
	  {
	    coarse_matrix(i,j) = 0;
	    coarse_matrix(j,i) = 0;
	  }
}



				 // The solution process again looks like
				 // step-16. We now use a Chebyshev smoother
				 // instead of SSOR (which is difficult to
				 // implement if we do not have the matrix
				 // elements explicitly available).
template <int dim>
void LaplaceProblem<dim>::solve () 
{
  GrowingVectorMemory<>   vector_memory;

  MGTransferPrebuilt<Vector<double> > mg_transfer;
  mg_transfer.build_matrices(mg_dof_handler);

  MGCoarseGridHouseholder<float, Vector<double> > mg_coarse;
  mg_coarse.initialize(coarse_matrix);

  typedef PreconditionChebyshev<MatrixFreeType,Vector<double> > SMOOTHER;
  MGSmootherPrecondition<MatrixFreeType, SMOOTHER, Vector<double> >
    mg_smoother(vector_memory);

				   // Initialize the smoother with our level
				   // matrices and the required, additional
				   // data for the Chebyshev smoother. Use a
				   // higher polynomial degree for higher
				   // order elements, since smoothing gets
				   // more difficult then. Smooth out a
				   // range of
				   // $[\lambda_{\max}/8,\lambda_{\max}]$.
  typename SMOOTHER::AdditionalData smoother_data;
  smoother_data.smoothing_range = 8.;
  smoother_data.degree = fe.degree+1;
  mg_smoother.initialize(mg_matrices, smoother_data);

  MGMatrix<MatrixFreeType, Vector<double> >
    mg_matrix(&mg_matrices);

  Multigrid<Vector<double> > mg(mg_dof_handler,
				mg_matrix,
				mg_coarse,
				mg_transfer,
				mg_smoother,
				mg_smoother);
  PreconditionMG<dim, Vector<double>,
    MGTransferPrebuilt<Vector<double> > >
    preconditioner(mg_dof_handler, mg, mg_transfer);

  double multigrid_memory = 
    (double)mg_matrices.memory_consumption() +
    (double)mg_transfer.memory_consumption() +
    (double)coarse_matrix.memory_consumption();

  std::cout << "Multigrid objects memory consumption: " 
	    << multigrid_memory*std::pow(2.,-20.) 
	    << " MBytes." 
	    << std::endl;

				   // Finally, create the solver
				   // object and solve the system
  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  cg.solve (system_matrix, solution, system_rhs,
  	    preconditioner);
  
  std::cout << "Convergence in " << solver_control.last_step() 
	    << " CG iterations." << std::endl;
}



				 // Here is the data output, which is
				 // a simplified version of step-5. We
				 // do a standard vtk output for
				 // each grid produced in the
				 // refinement process.
template <int dim>
void LaplaceProblem<dim>::output_results (const unsigned int cycle) const
{
  DataOut<dim> data_out;

  data_out.attach_dof_handler (mg_dof_handler);
  data_out.add_data_vector (solution, "solution");
  data_out.build_patches ();

  std::ostringstream filename;
  filename << "solution-"
	   << cycle
	   << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
}



template <int dim>
void LaplaceProblem<dim>::run () 
{
  for (unsigned int cycle=0; cycle<8-dim; ++cycle)
    {
      std::cout << "Cycle " << cycle << std::endl;

      if (cycle == 0)
	{
					   // Generate a simple hyperball grid.
	  GridGenerator::hyper_ball(triangulation);
	  static const HyperBallBoundary<dim> boundary;
	  triangulation.set_boundary (0, boundary);
	  triangulation.refine_global (3-dim);
	}
      triangulation.refine_global (1);
      setup_system ();
      assemble_system ();
      assemble_multigrid ();
      solve ();
      output_results (cycle);
      std::cout << std::endl;
    };
}

    

int main () 
{
  deallog.depth_console (0);
  LaplaceProblem<2> laplace_problem (2);
  laplace_problem.run ();
  
  return 0;
}
