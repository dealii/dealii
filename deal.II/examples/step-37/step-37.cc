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


				 // To start with the include files are more
				 // or less the same as in step-16:
#include <base/quadrature_lib.h>
#include <base/function.h>
#include <base/logstream.h>
#include <base/work_stream.h>

#include <lac/vector.h>
#include <lac/full_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>

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

#include <numerics/data_out.h>
#include <numerics/vectors.h>

#include <fstream>
#include <sstream>

using namespace dealii;



				 // @sect3{Equation data}

				 // We define a variable coefficient function
				 // for the Poisson problem. It is similar to
				 // the function in step-5 but we use the form
				 // $a(\mathbf x)=\frac{1}{0.1 + \|\bf x\|^2}$
				 // instead of a discontinuous one. It is
				 // merely to demonstrate the possibilities of
				 // this implementation, rather than making
				 // much sense physically.
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



				 // @sect3{Matrix-free implementation}

				 // In this program, we want to make
				 // use of the ability of deal.II to
				 // runs things in parallel if compute
				 // resources are available. We will
				 // follow the general framework laid
				 // out in the @ref threads module and
				 // use the WorkStream class to do
				 // operations on the range of all
				 // cells.
				 //
				 // To this end, we first have to have
				 // a few declarations that we use for
				 // defining the %parallel layout of
				 // the vector multiplication function
				 // with the WorkStream concept in the
				 // Matrix-free class. These comprise
				 // so-called scratch data that we use
				 // for calculating cell-related
				 // information, and copy data that is
				 // eventually used in a separate
				 // function for writing local data
				 // into the global vector. The reason
				 // for this split-up definition is
				 // that many threads at a time can
				 // execute the local multiplications
				 // (and filling up the copy data),
				 // but than that copy data needs to
				 // be worked on by one process at a
				 // time.
namespace WorkStreamData
{
  template <typename number>
  struct ScratchData
  {
      ScratchData ();
      ScratchData (const ScratchData &scratch);
      FullMatrix<number> solutions;
  };

  template<typename number>
  ScratchData<number>::ScratchData ()
		  :
		  solutions ()
  {}

  template<typename number>
  ScratchData<number>::ScratchData (const ScratchData &scratch)
		  :
		  solutions ()
  {}

  template <typename number>
  struct CopyData : public ScratchData<number>
  {
      CopyData ();
      CopyData (const CopyData &scratch);
      unsigned int first_cell;
      unsigned int n_dofs;
  };

  template <typename number>
  CopyData<number>::CopyData ()
		  :
		  ScratchData<number> ()
  {}

  template <typename number>
  CopyData<number>::CopyData (const CopyData &scratch)
		  :
		  ScratchData<number> ()
  {}

}



				 // Next comes the implementation of the
				 // matrix-free class. It provides some
				 // standard information we expect for
				 // matrices (like returning the dimensions
				 // of the matrix), it implements
				 // matrix-vector multiplications in several
				 // forms, and it provides functions for
				 // filling the matrix with data.
				 //
				 // We choose to make this class generic,
				 // i.e., we do not implement the actual
				 // differential operator (here: Laplace
				 // operator) directly in this class.  We
				 // instead let the actual transformation
				 // (which happens on the level of quadrature
				 // points, see the discussion in the
				 // introduction) be a template parameter that
				 // is implemented by another class. We then
				 // only have to store a list of these objects
				 // for each quadrature point on each cell in
				 // a big list &ndash; we choose a
				 // <code>Table<2,Transformation></code> data
				 // format) &ndash; and call a transform
				 // command of the <code>Transformation</code>
				 // class. This template magic makes it easy
				 // to reuse this MatrixFree class for other
				 // problems that are based on a symmetric
				 // operation without the need for further
				 // changes.
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

    number el (const unsigned int row,
	       const unsigned int col) const;
    void calculate_diagonal () const;

    std::size_t memory_consumption () const;

				     // The private member variables of the
				     // <code>MatrixFree</code> class are a
				     // small matrix that does the
				     // transformation from solution values to
				     // quadrature points, a list with the
				     // mapping between local degrees of freedom
				     // and global degrees of freedom for each
				     // cell (stored as a two-dimensional array,
				     // where each row corresponds to one
				     // cell, and the columns within individual
				     // cells are the local degrees of freedom),
				     // the transformation variable for
				     // implementing derivatives, a constraint
				     // matrix for handling boundary conditions
				     // as well as a few other variables that
				     // store matrix properties.
  private:
    typedef std::vector<std::pair<unsigned int,unsigned int> >::const_iterator
    CellChunkIterator;
    template <typename number2>
    void local_vmult (CellChunkIterator                    cell_range,
		      WorkStreamData::ScratchData<number> &scratch,
		      WorkStreamData::CopyData<number>    &copy,
		      const Vector<number2>               &src) const;

    template <typename number2>
    void
    copy_local_to_global (const WorkStreamData::CopyData<number> &copy,
			  Vector<number2>                        &dst) const;

    FullMatrix<number>      B_ref_cell;
    Table<2,unsigned int>   indices_local_to_global;
    Table<2,Transformation> derivatives;

    ConstraintMatrix        constraints;

    mutable Vector<number>  diagonal_values;
    mutable bool            diagonal_is_calculated;

    struct MatrixSizes
    {
	unsigned int n_dofs, n_cells;
	unsigned int m, n;
	unsigned int n_points, n_comp;
	std::vector<std::pair<unsigned int,unsigned int> > chunks;
    }  matrix_sizes;
};



				 // This is the constructor of the
				 // <code>MatrixFree</code> class. All it does
				 // is to subscribe to the general deal.II
				 // <code>Subscriptor</code> scheme that makes
				 // sure that we do not delete an object of
				 // this class as long as it used somewhere
				 // else, e.g. in a preconditioner.
template <typename number, class Transformation>
MatrixFree<number,Transformation>::MatrixFree ()
		:
		Subscriptor()
{}



				 // The next functions return the
				 // number of rows and columns of the
				 // global matrix (i.e. the dimensions
				 // of the operator this class
				 // represents, the point of this
				 // tutorial program was, after all,
				 // that we don't actually store the
				 // elements of the rows and columns
				 // of this operator). Since the
				 // matrix is square, the returned
				 // numbers are the same.
template <typename number, class Transformation>
unsigned int
MatrixFree<number,Transformation>::m () const
{
  return matrix_sizes.n_dofs;
}



template <typename number, class Transformation>
unsigned int
MatrixFree<number,Transformation>::n () const
{
  return matrix_sizes.n_dofs;
}



				 // One more function that just returns an
				 // %internal variable. Note that the user
				 // will need to change this variable, so it
				 // returns a non-constant reference to the
				 // ConstraintMatrix.
template <typename number, class Transformation>
ConstraintMatrix &
MatrixFree<number,Transformation>::get_constraints ()
{
  return constraints;
}



				 // The following function takes a
				 // vector of local dof indices on
				 // cell level and writes the data
				 // into the
				 // <code>indices_local_to_global</code>
				 // field in order to have fast access
				 // to it. It performs a few sanity
				 // checks like whether the sizes in
				 // the matrix are set correctly. One
				 // tiny thing: Whenever we enter this
				 // function, we probably make some
				 // modification to the matrix. This
				 // means that the diagonal of the
				 // matrix, which we might have
				 // computed to have fast access to
				 // those elements, is invalidated. We
				 // set the respective flag to
				 // <code>false</code>.
template <typename number, class Transformation>
void MatrixFree<number,Transformation>::
set_local_dof_indices (const unsigned int               cell_no,
		       const std::vector<unsigned int> &local_dof_indices)
{
  Assert (local_dof_indices.size() == matrix_sizes.m,
	  ExcDimensionMismatch(local_dof_indices.size(),
			       matrix_sizes.m));
  for (unsigned int i=0; i<matrix_sizes.m; ++i)
    {
      Assert (local_dof_indices[i] < matrix_sizes.n_dofs, ExcInternalError());
      indices_local_to_global(cell_no,i) = local_dof_indices[i];
    }
  diagonal_is_calculated = false;
}



				 // Next a function that writes the
				 // derivative data on a certain cell
				 // and a certain quadrature point to
				 // the array that keeps the data
				 // around. Even though the array
				 // <code>derivatives</code> stands
				 // for the majority of the matrix
				 // memory consumption, it still pays
				 // off to have that data around since
				 // it would be quite expensive to
				 // manually compute it every time we
				 // make a matrix-vector product.
template <typename number, class Transformation>
void MatrixFree<number,Transformation>::
set_derivative_data (const unsigned int cell_no,
		     const unsigned int quad_point,
		     const Transformation &trans_in)
{
  Assert (quad_point < matrix_sizes.n_points, ExcInternalError());
  derivatives(cell_no,quad_point) = trans_in;
  diagonal_is_calculated = false;
}



				 // Now finally to the central function of the
				 // matrix-free class, implementing the
				 // multiplication of the matrix with a
				 // vector. This function does not actually
				 // work on all cells of a mesh, but only the subset
				 // of cells specified by the first argument
				 // <code>cell_range</code>. Since this
				 // function operates similarly irrespective
				 // on which cell chunk we are sitting, we can
				 // call it simultaneously on many processors,
				 // but with different cell range data.
				 //
				 // The goal of this function is to provide
				 // the multiplication of a vector with the
				 // local contributions of a set of cells. As
				 // mentioned in the introduction, if we were
				 // to deal with a single cell, this would
				 // amount to performing the product 
				 // @f{eqnarray*}
				 // P^T_\mathrm{cell,local-global} A_\mathrm{cell} 
				 // P_\mathrm{cell,local-global} x
				 // @f}
				 // where 
				 // @f{eqnarray*}
				 // A_\mathrm{cell} = 
				 // B_\mathrm{ref\_cell}^T J_\mathrm{cell}^T 
				 // D_\mathrm{cell} 
				 // J_\mathrm{cell} B_\mathrm{ref\_cell}
				 // @f}
				 // and $P_\mathrm{cell,local-global}$ is the
				 // transformation from local to global
				 // indices.
				 //
				 // To do this, we would have to do the
				 // following steps:
				 // <ol>
				 //   <li> Form $x_\mathrm{cell} =
				 //   P_\mathrm{cell,local-global} x$. This is
				 //   done using the command
				 //   ConstraintMatrix::get_dof_values.
				 //   <li> Form $x_1 = B_\mathrm{ref\_cell}
				 //   x_\mathrm{cell}$. The vector $x_1$
				 //   contains the reference cell gradient to
				 //   the local cell vector.
				 //   <li> Form $x_2 = J_\mathrm{cell}^T
				 //   D_\mathrm{cell} J_\mathrm{cell}
				 //   x_1$. This is a block-diagonal
				 //   operation, with the block size equal to
				 //   <code>dim</code>. The blocks just
				 //   correspond to the individual quadrature
				 //   points. The operation on each quadrature
				 //   point is implemented by the
				 //   Transformation class object that this
				 //   class is equipped with. Compared to the
				 //   introduction, the matrix
				 //   $D_\mathrm{cell}$ now contains the
				 //   <code>JxW</code> values and the
				 //   inhomogeneous coefficient.
				 //   <li> Form $y_\mathrm{cell} =
				 //   B_\mathrm{ref\_cell}^T x_2$. This gives
				 //   the local result of the matrix-vector
				 //   product.
				 //   <li> Form $y +=
				 //   P_\mathrm{cell,local-global}^T
				 //   y_\mathrm{cell}$. This adds the local
				 //   result to the global vector, which is
				 //   realized using the method
				 //   ConstraintMatrix::distribute_local_to_global.
				 //   Note that we do this in an extra
				 //   function called
				 //   <code>copy_local_to_global</code>
				 //   because that operation must not be done
				 //   in %parallel, in order to avoid two or
				 //   more processes trying to add to the same
				 //   positions in the result vector $y$.
				 //   </ol>
				 // The steps 1 to 4 can be done in %parallel
				 // by multiple processes.

				 // Now, it turns out that the most expensive
				 // part of the above is the multiplication
				 // $B_\mathrm{ref\_cell} x_\mathrm{cell}$ in
				 // the second step and the transpose
				 // operation in step 4. Note that the matrix
				 // $J^T D J$ is block-diagonal, and hence,
				 // its application is cheaper. Since the
				 // matrix $B_\mathrm{ref\_cell}$ is the same
				 // for all cells, all that changes is the
				 // vector $x_\mathrm{cell}$. Hence, nothing
				 // prevents us from collecting several cell
				 // vectors to a (rectangular) matrix, and
				 // then perform a matrix-matrix
				 // product. These matrices are both full, but
				 // not very large, having of the order
				 // <code>dofs_per_cell</code> rows and
				 // columns. This is an operation that can be
				 // much better optimized than matrix-vector
				 // products. The functions
				 // <code>FullMatrix<number>::mmult</code> and
				 // <code>FullMatrix<number>::mTmult</code>
				 // use the BLAS dgemm function (as long as
				 // BLAS has been detected in deal.II
				 // configuration), which provides optimized
				 // kernels for doing this product. In our
				 // case, a matrix-matrix product is between
				 // three and five times faster than doing the
				 // matrix-vector product on one cell after
				 // the other. The variables that hold the
				 // solution on the respective cell's support
				 // points and the quadrature points are thus
				 // full matrices. The number of rows is given
				 // by the number of cells they work on, and
				 // the number of columns is the number of
				 // degrees of freedom per cell for the first
				 // and the number of quadrature points times
				 // the number of components per point for the
				 // latter.
template <typename number, class Transformation>
template <typename number2>
void
MatrixFree<number,Transformation>::
local_vmult (CellChunkIterator                    cell_range,
	     WorkStreamData::ScratchData<number> &scratch,
	     WorkStreamData::CopyData<number>    &copy,
	     const Vector<number2>               &src) const
{
  const unsigned int chunk_size = cell_range->second - cell_range->first;

  copy.solutions.reinit    (chunk_size,matrix_sizes.m, true);
  copy.first_cell         = cell_range->first;
  copy.n_dofs             = chunk_size*matrix_sizes.m;
  scratch.solutions.reinit (chunk_size,matrix_sizes.n, true);

  constraints.get_dof_values(src, &indices_local_to_global(copy.first_cell,0),
			     &copy.solutions(0,0),
			     &copy.solutions(0,0)+copy.n_dofs);

  copy.solutions.mmult (scratch.solutions, B_ref_cell);

  for (unsigned int i=0, k = copy.first_cell; i<chunk_size; ++i, ++k)
    for (unsigned int j=0; j<matrix_sizes.n_points; ++j)
      derivatives(k,j).transform(&scratch.solutions(i, j*matrix_sizes.n_comp));

  scratch.solutions.mTmult (copy.solutions, B_ref_cell);
}



template <typename number, class Transformation>
template <typename number2>
void
MatrixFree<number,Transformation>::
copy_local_to_global (const WorkStreamData::CopyData<number> &copy,
		      Vector<number2>                        &dst) const
{
  constraints.distribute_local_to_global (&copy.solutions(0,0),
					  &copy.solutions(0,0)+copy.n_dofs,
					  &indices_local_to_global(copy.first_cell,0),
					  dst);
}



				 // Now to the <code>vmult</code> function
				 // that is called externally: It is very
				 // similar to the <code>vmult_add</code>
				 // function, so just set the destination to
				 // zero first, and then go to the other
				 // function.
template <typename number, class Transformation>
template <typename number2>
void
MatrixFree<number,Transformation>::vmult (Vector<number2>       &dst,
					  const Vector<number2> &src) const
{
  dst = 0;
  vmult_add (dst, src);
}



				 // Transposed matrix-vector products: do
				 // the same. Since we implement a symmetric
				 // operation, we can refer to the vmult_add
				 // operation.
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
MatrixFree<number,Transformation>::Tvmult_add (Vector<number2>       &dst,
					       const Vector<number2> &src) const
{
  vmult_add (dst,src);
}



				 // This is the <code>vmult_add</code>
				 // function that multiplies the
				 // matrix with vector
				 // <code>src</code> and adds the
				 // result to vector <code>dst</code>.
				 // We include a few sanity checks to
				 // make sure that the size of the
				 // vectors is the same as the
				 // dimension of the matrix. We call a
				 // %parallel function that applies
				 // the multiplication on a chunk of
				 // cells at once using the WorkStream
				 // module (cf. also the @ref threads
				 // module). The subdivision into
				 // chunks will be performed in the
				 // reinit function and is stored in
				 // the field
				 // <code>matrix_sizes.chunks</code>. What
				 // the rather cryptic command to
				 // <code>std_cxx1x::bind</code> does
				 // is to transform a function that
				 // has several arguments (source
				 // vector, chunk information) into a
				 // function which has three arguments
				 // (in the first case) or one
				 // argument (in the second), which is
				 // what the WorkStream::run function
				 // expects. The placeholders
				 // <code>_1, _2, _3</code> in the
				 // local vmult specify variable input
				 // values, given by the chunk
				 // information, scratch data and copy
				 // data that the WorkStream::run
				 // function will provide, whereas the
				 // other arguments to the
				 // <code>local_vmult</code> function
				 // are bound: to <code>this</code>
				 // and a constant reference to the
				 // <code>src</code> in the first
				 // case, and <code>this</code> and a
				 // reference to the output vector in
				 // the second. Similarly, the
				 // placeholder <code>_1</code> in the
				 // <code>copy_local_to_global</code>
				 // function sets the first explicit
				 // argument of that function, which
				 // is of class
				 // <code>CopyData</code>. We need to
				 // abstractly specify these arguments
				 // because the tasks defined by
				 // different cell chunks will be
				 // scheduled by the WorkStream class,
				 // and we will reuse available
				 // scratch and copy data.
template <typename number, class Transformation>
template <typename number2>
void
MatrixFree<number,Transformation>::vmult_add (Vector<number2>       &dst,
					      const Vector<number2> &src) const
{
  Assert (src.size() == n(), ExcDimensionMismatch(src.size(), n()));
  Assert (dst.size() == m(), ExcDimensionMismatch(dst.size(), m()));

  WorkStream::run (matrix_sizes.chunks.begin(), matrix_sizes.chunks.end(),
		   std_cxx1x::bind(&MatrixFree<number,Transformation>::
				   template local_vmult<number2>,
				   this, _1, _2, _3, boost::cref(src)),
		   std_cxx1x::bind(&MatrixFree<number,Transformation>::
				   template copy_local_to_global<number2>,
				   this, _1, boost::ref(dst)),
		   WorkStreamData::ScratchData<number>(),
		   WorkStreamData::CopyData<number>(),
		   2*multithread_info.n_default_threads,1);

				   // One thing to be cautious about:
				   // The deal.II classes expect that
				   // the matrix still contains a
				   // diagonal entry for constrained
				   // dofs (otherwise, the matrix
				   // would be singular, which is not
				   // what we want). Since the
				   // <code>distribute_local_to_global</code>
				   // command of the constraint matrix
				   // which we used for adding the
				   // local elements into the global
				   // vector does not do anything with
				   // constrained elements, we have to
				   // circumvent that problem by
				   // artificially setting the
				   // diagonal to some non-zero value
				   // and adding the source values. We
				   // simply set it to one, which
				   // corresponds to copying the
				   // respective elements of the
				   // source vector into the matching
				   // entry of the destination vector.
  for (unsigned int i=0; i<matrix_sizes.n_dofs; ++i)
    if (constraints.is_constrained(i) == true)
      dst(i) += 1.0 * src(i);
}



				 // The next function initializes the structures
				 // of the matrix. It writes the number of
				 // total degrees of freedom in the problem
				 // as well as the number of cells to the
				 // MatrixSizes struct and copies the small
				 // matrix that transforms the solution from
				 // support points to quadrature points. It
				 // uses the small matrix for determining
				 // the number of degrees of freedom per
				 // cell (number of rows in
				 // <code>B_ref_cell</code>). The number
				 // of quadrature points needs to be passed
				 // through the last variable
				 // <code>n_points_per_cell</code>, since
				 // the number of columns in the small
				 // matrix is
				 // <code>dim*n_points_per_cell</code> for
				 // the Laplace problem (the Laplacian is a
				 // tensor and has <code>dim</code>
				 // components). In this function, we also
				 // give the fields containing the
				 // derivative information and the local dof
				 // indices the correct sizes. They will be
				 // filled by calling the respective set
				 // function defined above.
template <typename number, class Transformation>
void MatrixFree<number,Transformation>::
reinit (const unsigned int        n_dofs_in,
	const unsigned int        n_cells_in,
	const FullMatrix<double> &B_ref_cell_in,
	const unsigned int        n_points_per_cell)
{
  B_ref_cell = B_ref_cell_in;

  derivatives.reinit (n_cells_in, n_points_per_cell);
  indices_local_to_global.reinit (n_cells_in, B_ref_cell.m());

  diagonal_is_calculated = false;

  matrix_sizes.n_dofs = n_dofs_in;
  matrix_sizes.n_cells = n_cells_in;
  matrix_sizes.m = B_ref_cell.m();
  matrix_sizes.n = B_ref_cell.n();
  matrix_sizes.n_points = n_points_per_cell;
  matrix_sizes.n_comp   = B_ref_cell.n()/matrix_sizes.n_points;
  Assert(matrix_sizes.n_comp * n_points_per_cell == B_ref_cell.n(),
	 ExcInternalError());

				   // One thing to make the matrix-vector
				   // product with this class efficient is to
				   // decide how many cells should be combined
				   // to one chunk, which will determine the
				   // size of the full matrix that we work
				   // on. If we choose too few cells, then the
				   // gains from using the matrix-matrix product
				   // will not be fully utilized (dgemm tends to
				   // provide more efficiency the larger the
				   // matrix dimensions get). If we choose too
				   // many, we will firstly degrade
				   // parallelization (we need to have
				   // sufficiently independent tasks), and
				   // secondly introduce an inefficiency that
				   // comes from the computer architecture: In
				   // the actual working function above, right
				   // after the first matrix-matrix
				   // multiplication, we transform the solution
				   // on quadrature points by using
				   // derivatives. Obviously, we want to have
				   // fast access to that data, so it should
				   // still be present in processor cache and
				   // not needed to be fetched from main
				   // memory. The total memory usage of the data
				   // on quadrature points should not be more
				   // than about a third of the cache size of
				   // the processor in order to be on the safe
				   // side. Since most of today's processors
				   // provide 512 kB or more cache memory per
				   // core, we choose about 150 kB as a size to
				   // leave some room for other things to be
				   // stored in the CPU. Clearly, this is an
				   // architecture-dependent value and the
				   // interested user can squeeze out some extra
				   // performance by hand-tuning this
				   // parameter. Once we have chosen the number
				   // of cells we collect in one chunk, we
				   // determine how many chunks we have on the
				   // given cell range and recalculate the
				   // actual chunk size in order to evenly
				   // distribute the chunks.
  const unsigned int divisor = 150000/(matrix_sizes.n*sizeof(double));
  const unsigned int n_chunks = std::max (matrix_sizes.n_cells/divisor + 1,
					  2*multithread_info.n_default_threads);

  const unsigned int chunk_size = (matrix_sizes.n_cells/n_chunks +
				   (matrix_sizes.n_cells%n_chunks>0));

  std::pair<unsigned int, unsigned int> chunk;
  for (unsigned int i=0; i<n_chunks; ++i)
    {
      chunk.first = i*chunk_size;
      if ((i+1)*chunk_size > matrix_sizes.n_cells)
	chunk.second = matrix_sizes.n_cells;
      else
	chunk.second = (i+1)*chunk_size;

      if (chunk.second > chunk.first)
	matrix_sizes.chunks.push_back(chunk);
      else
	break;
    }
}



				 // Then we need a function if we want to
				 // delete the content of the matrix,
				 // e.g. when we are finished with one grid
				 // level and continue to the next one. Just
				 // set all the field sizes to 0.
template <typename number, class Transformation>
void
MatrixFree<number,Transformation>::clear ()
{
  B_ref_cell.reinit(0,0);
  derivatives.reinit (0,0);
  indices_local_to_global.reinit(0,0);

  constraints.clear();

  diagonal_values.reinit (0);
  diagonal_is_calculated = false;

  matrix_sizes.n_dofs = 0;
  matrix_sizes.n_cells = 0;
  matrix_sizes.chunks.clear();
}



				 // The next function returns the entries of the
				 // matrix. Since this class is intended not
				 // to store the matrix entries, it would make
				 // no sense to provide all those
				 // elements. However, diagonal entries are
				 // explicitly needed for the implementation
				 // of the Chebyshev smoother that we intend
				 // to use in the multigrid
				 // preconditioner. This matrix is equipped
				 // with a vector that stores the diagonal,
				 // and we compute it when this function is
				 // called for the first time.
template <typename number, class Transformation>
number
MatrixFree<number,Transformation>::el (const unsigned int row,
				       const unsigned int col) const
{
  Assert (row == col, ExcNotImplemented());
  if (diagonal_is_calculated == false)
    calculate_diagonal();

  return diagonal_values(row);
}



				 // Regarding the calculation of the diagonal,
				 // remember that this is as simple (or
				 // complicated) as assembling a right hand
				 // side in deal.II. Well, it is a bit easier
				 // to do this within this class since we have
				 // all the derivative information
				 // available. What we do is to go through all
				 // the cells (now in serial, since this
				 // function should not be called very often
				 // anyway), then all the degrees of
				 // freedom. At this place, we first copy the
				 // first basis functions in all the
				 // quadrature points to a temporary array,
				 // apply the derivatives from the Jacobian
				 // matrix, and finally multiply with the
				 // second basis function. This is exactly the
				 // value that would be written into the
				 // diagonal of a sparse matrix. Note that we
				 // need to condense hanging node constraints
				 // and set the constrained diagonals to one.
template <typename number, class Transformation>
void
MatrixFree<number,Transformation>::calculate_diagonal() const
{
  diagonal_values.reinit (matrix_sizes.n_dofs);
  std::vector<number> calculation (matrix_sizes.n);
  for (unsigned int cell=0; cell<matrix_sizes.n_cells; ++cell)
    for (unsigned int dof=0; dof<matrix_sizes.m; ++dof)
      {
	memcpy (&calculation[0],&B_ref_cell(dof,0),
		matrix_sizes.n*sizeof(number));
	for (unsigned int q=0; q<matrix_sizes.n_points; ++q)
	  derivatives(cell,q).transform(&calculation[q*matrix_sizes.n_comp]);
	double diag_value = 0;
	for (unsigned int q=0; q<matrix_sizes.n; ++q)
	  diag_value += calculation[q] * B_ref_cell(dof,q);
	diagonal_values(indices_local_to_global(cell,dof)) += diag_value;
      }
  constraints.condense (diagonal_values);
  for (unsigned int i=0; i<matrix_sizes.n_dofs; ++i)
    if (constraints.is_constrained(i) == true)
      diagonal_values(i) = 1.0;

  diagonal_is_calculated = true;
}



				 // Eventually, we provide a function that
				 // calculates how much memory this class
				 // uses. We just need to sum up the memory
				 // consumption of the arrays, the
				 // constraints, the small matrix and of the
				 // local variables. Just as a remark: In 2D
				 // and with data type <code>double</code>,
				 // about 80 per cent of the memory
				 // consumption is due to the
				 // <code>derivatives</code> array, while in 3D
				 // this number is even 85 per cent.
template <typename number, class Transformation>
std::size_t MatrixFree<number,Transformation>::memory_consumption () const
{
  std::size_t glob_size = derivatives.memory_consumption() +
			  indices_local_to_global.memory_consumption() +
			  constraints.memory_consumption() +
			  B_ref_cell.memory_consumption() +
			  diagonal_values.memory_consumption() +
			  matrix_sizes.chunks.size()*2*sizeof(unsigned int) +
			  sizeof(*this);
  return glob_size;
}



				 // @sect3{Laplace operator implementation}

				 // This class implements the local action
				 // of a Laplace operator on a quadrature
				 // point. This is a very basic class
				 // implementation, providing functions for
				 // initialization with a Tensor of rank 2
				 // and implementing the
				 // <code>transform</code> operation needed
				 // by the <code>MatrixFree</code>
				 // class. There is one point worth noting:
				 // The quadrature-point related action of
				 // the Laplace operator is a tensor of rank
				 // two. It is symmetric since it is the
				 // product of the inverse Jacobian
				 // transformation between unit and real
				 // cell with its transpose (times
				 // quadrature weights and a coefficient,
				 // which are scalar), so we can just save
				 // the diagonal and upper diagonal part. We
				 // could use the SymmetricTensor<2,dim>
				 // class for doing this, however, that
				 // class is only based on
				 // <code>double</code> %numbers. Since we
				 // also want to use <code>float</code>
				 // %numbers for the multigrid
				 // preconditioner (in order to save memory
				 // and computing time), we manually
				 // implement this operator. Note that
				 // <code>dim</code> is a template argument
				 // and hence known at compile-time, so the
				 // compiler knows that this symmetric
				 // rank-2 tensor has 3 entries if used in
				 // 2D and 6 entries if used in 3D.
template <int dim,typename number>
class LaplaceOperator
{
  public:
    LaplaceOperator ();

    LaplaceOperator (const Tensor<2,dim> &tensor);

    void transform (number * result) const;

    LaplaceOperator<dim,number>&
    operator = (const Tensor<2,dim> &tensor);

  private:
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

				 // Now implement the transformation, which is
				 // just a so-called contraction
				 // operation between a tensor of rank two and a
				 // tensor of rank one. Unfortunately, we
				 // need to implement this by hand, since we
				 // chose not to use the
				 // SymmetricTensor<2,dim> class (note that
				 // the resulting values are entries in a full
				 // matrix that consists of doubles or
				 // floats). It feels a bit unsafe to operate
				 // on a pointer to the data, but that is the
				 // only possibility if we do not want to copy
				 // data back and forth, which is expensive
				 // since this is the innermost position of
				 // the loop in the <code>vmult</code>
				 // operation of the MatrixFree class. We need
				 // to pay attention to the fact that we only
				 // saved half of the (symmetric) rank-two
				 // tensor.
				 //
				 // At first sight, it seems inefficient that
				 // we have an <code>if</code> clause at this
				 // position in the code at the innermost
				 // loop, but note once again that
				 // <code>dim</code> is known when this piece
				 // of code is compiled, so the compiler can
				 // optimize away the <code>if</code> statement
				 // (and actually even inline these few lines
				 // of code into the <code>MatrixFree</code>
				 // class).
template <int dim, typename number>
void LaplaceOperator<dim,number>::transform (number* result) const
{
  if (dim == 2)
    {
      const number temp = result[0];
      result[0] = transformation[0] * temp + transformation[1] * result[1];
      result[1] = transformation[1] * temp + transformation[2] * result[1];
    }
  else if (dim == 3)
    {
      const number temp1 = result[0];
      const number temp2 = result[1];
      result[0] = transformation[0] * temp1 + transformation[1] * temp2 +
		  transformation[2] * result[2];
      result[1] = transformation[1] * temp1 + transformation[3] * temp2 +
		  transformation[4] * result[2];
      result[2] = transformation[2] * temp1 + transformation[4] * temp2 +
		  transformation[5] * result[2];
    }
  else
    ExcNotImplemented();
}

				 // The final function in this group
				 // takes the content of a rank-2
				 // tensor and writes it to the field
				 // <code>transformation</code> of
				 // this class. We save the upper part
				 // of the symmetric tensor row-wise:
				 // we first take the (0,0)-entry,
				 // then the (0,1)-entry, and so
				 // on. We only implement this for
				 // dimensions two and three, which
				 // for the moment should do just
				 // fine:
template <int dim, typename number>
LaplaceOperator<dim,number>&
LaplaceOperator<dim,number>::operator=(const Tensor<2,dim> &tensor)
{
  if (dim == 2)
    {
      transformation[0] = tensor[0][0];
      transformation[1] = tensor[0][1];
      transformation[2] = tensor[1][1];
      Assert (std::fabs(tensor[1][0]-tensor[0][1])<1e-15,
	      ExcInternalError());
    }
  else if (dim == 3)
    {
      transformation[0] = tensor[0][0];
      transformation[1] = tensor[0][1];
      transformation[2] = tensor[0][2];
      transformation[3] = tensor[1][1];
      transformation[4] = tensor[1][2];
      transformation[5] = tensor[2][2];
      Assert (std::fabs(tensor[1][0]-tensor[0][1])<1e-15,
	      ExcInternalError());
      Assert (std::fabs(tensor[2][0]-tensor[0][2])<1e-15,
	      ExcInternalError());
      Assert (std::fabs(tensor[2][1]-tensor[1][2])<1e-15,
	      ExcInternalError());
    }
  else
    ExcNotImplemented();
  return *this;
}



				 // @sect3{LaplaceProblem class}

				 // This class is based on the same
				 // class in step-16. However, we
				 // replaced the SparseMatrix<double>
				 // class by our matrix-free
				 // implementation, which means that
				 // we can also skip the sparsity
				 // patterns.
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
LaplaceProblem<dim>::LaplaceProblem (const unsigned int degree)
		:
                fe (degree),
		mg_dof_handler (triangulation)
{}



				 // @sect4{LaplaceProblem::setup_system}

				 // This is the function of step-16 with
				 // relevant changes due to the MatrixFree
				 // class. What we need to do is to somehow
				 // create a local gradient matrix that does
				 // not contain any cell-related data
				 // (gradient on the reference cell). The
				 // way to get to this matrix is to create
				 // an FEValues object with gradient
				 // information on a cell that corresponds
				 // to the reference cell, which is a cube
				 // with side length 1. So we create a
				 // pseudo triangulation, initialize the
				 // FEValues to the only cell of that
				 // triangulation, and read off the
				 // gradients (which we put in a
				 // FullMatrix). That full matrix is then
				 // passed to the reinit function of the
				 // MatrixFree class used as a system matrix
				 // and, further down, as multigrid matrices
				 // on the individual levels. We need to
				 // implement Dirichlet boundary conditions
				 // here, which is done with the
				 // ConstraintMatrix function as shown,
				 // e.g., in step-22.
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
  FEValues<dim> fe_values_reference (fe, quadrature_formula,
				     update_gradients);
  Triangulation<dim> reference_cell;
  GridGenerator::hyper_cube (reference_cell, 0, 1);
  fe_values_reference.reinit (reference_cell.begin());
  FullMatrix<double> data_matrix (fe.dofs_per_cell,
				  quadrature_formula.size()*dim);
  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    {
      for (unsigned int j=0; j<quadrature_formula.size(); ++j)
	{
	  for (unsigned int d=0; d<dim; ++d)
	    data_matrix(i,j*dim+d) = fe_values_reference.shape_grad(i,j)[d];
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
	    << system_matrix.memory_consumption()/std::pow(2.,20.)
	    << " MBytes."
	    << std::endl;

  solution.reinit (mg_dof_handler.n_dofs());
  system_rhs.reinit (mg_dof_handler.n_dofs());

				   // Next, initialize the matrices for the
				   // multigrid method on all the
				   // levels. Unfortunately, the function
				   // MGTools::make_boundary_list cannot write
				   // Dirichlet boundary conditions into a
				   // ConstraintMatrix object directly, so we
				   // first have to make the boundary list and
				   // then manually fill the boundary
				   // conditions using the command
				   // ConstraintMatrix::add_line. Once this is
				   // done, we close the ConstraintMatrix so
				   // it can be used for matrix-vector
				   // products.
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



				 // @sect4{LaplaceProblem::assemble_system}

				 // The assemble function is significantly
				 // reduced compared to step-16. All we need
				 // to do is to assemble the right hand side
				 // and to calculate the cell-dependent part
				 // of the Laplace operator. The first task
				 // is standard. The second is also not too
				 // hard given the discussion in the
				 // introduction: We need to take the
				 // inverse of the Jacobian of the
				 // transformation from unit to real cell,
				 // multiply it with its transpose and
				 // multiply the resulting rank-2 tensor
				 // with the quadrature weights and the
				 // coefficient values at the quadrature
				 // points. To make this work, we add the
				 // update flag
				 // <code>update_inverse_jacobians</code> to
				 // the FEValues constructor, and query the
				 // inverse of the Jacobian in a loop over
				 // the quadrature points (note that the
				 // Jacobian is not related to any kind of
				 // degrees of freedom directly). In the
				 // end, we condense the constraints from
				 // Dirichlet boundary conditions away from
				 // the right hand side.
template <int dim>
void LaplaceProblem<dim>::assemble_system ()
{
  QGauss<dim>  quadrature_formula(fe.degree+1);
  MappingQ<dim> mapping (fe.degree);
  FEValues<dim> fe_values (mapping, fe, quadrature_formula,
			   update_values   | update_inverse_jacobians |
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

      for (unsigned int i=0; i<dofs_per_cell; ++i)
	{
	  double rhs_val = 0;
	  for (unsigned int q=0; q<n_q_points; ++q)
	    rhs_val += (fe_values.shape_value(i,q) * 1.0 *
			fe_values.JxW(q));
	  system_rhs(local_dof_indices[i]) += rhs_val;
	}

      system_matrix.set_local_dof_indices (cell_no, local_dof_indices);
      for (unsigned int q=0; q<n_q_points; ++q)
	system_matrix.set_derivative_data (cell_no, q,
					   (transpose
					    (fe_values.inverse_jacobian(q)) *
					    fe_values.inverse_jacobian(q)) *
					   fe_values.JxW(q) *
					   coefficient_values[q]);
    }
  system_matrix.get_constraints().condense(system_rhs);
}


				 // @sect4{LaplaceProblem::assemble_multigrid}

				 // Here is another assemble
				 // function. The integration core is
				 // the same as above. Only the loop
				 // goes over all existing cells now
				 // and the results must be entered
				 // into the correct matrix.

				 // Since we only do multilevel
				 // preconditioning, no right-hand side is
				 // assembled here. Compared to step-16, there
				 // is one new thing here: we manually
				 // calculate the matrix on the coarsest
				 // level. In step-16, we could simply copy
				 // the entries from the respective sparse
				 // matrix, but this is obviously not possible
				 // here. We could have integrated this into the
				 // MatrixFree class as well, but it is simple
				 // enough, so calculate it here instead.
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
	mg_matrices[level].set_derivative_data (cell_no[level], q,
						(transpose
						 (fe_values.inverse_jacobian(q)) *
						 fe_values.inverse_jacobian(q)) *
						fe_values.JxW(q) *
						coefficient_values[q]);

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

				   // In a final step, we need to
				   // condense the boundary conditions
				   // on the coarse matrix. There is
				   // no built-in function for doing
				   // this on a full matrix, so
				   // manually delete the rows and
				   // columns of the matrix that are
				   // constrained.
  for (unsigned int i=0; i<coarse_matrix.m(); ++i)
    if (mg_matrices[0].get_constraints().is_constrained(i))
      for (unsigned int j=0; j<coarse_matrix.n(); ++j)
	if (i!=j)
	  {
	    coarse_matrix(i,j) = 0;
	    coarse_matrix(j,i) = 0;
	  }
}



				 // @sect4{LaplaceProblem::solve}

				 // The solution process again looks like
				 // step-16. We now use a Chebyshev smoother
				 // instead of SOR (SOR would be very
				 // difficult to implement because we do not
				 // have the matrix elements available
				 // explicitly, and it is difficult to make it
				 // work efficiently in %parallel). The
				 // multigrid classes provide a simple
				 // interface for using the Chebyshev smoother
				 // which is defined in a preconditioner
				 // class: MGSmootherPrecondition.
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

				   // Then, we initialize the smoother
				   // with our level matrices and the
				   // required, additional data for
				   // the Chebyshev smoother. In
				   // particular, we use a higher
				   // polynomial degree for higher
				   // order elements, since smoothing
				   // gets more difficult for
				   // these. Smooth out a range of
				   // $[\lambda_{\max}/10,\lambda_{\max}]$. In
				   // order to compute the maximum
				   // eigenvalue of the corresponding
				   // matrix, the Chebyshev
				   // initializations performs a few
				   // steps of a CG algorithm. Since
				   // all we need is a rough estimate,
				   // we choose some eight iterations
				   // (more if the finite element
				   // polynomial degree is larger,
				   // less if it is smaller than
				   // quadratic).
  typename SMOOTHER::AdditionalData smoother_data;
  smoother_data.smoothing_range = 10.;
  smoother_data.degree = fe.degree;
  smoother_data.eig_cg_n_iterations = 4+2*fe.degree;
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

				 // Finally, write out the memory
				 // consumption of the Multigrid object
				 // (or rather, of its most significant
				 // components, since there is no built-in
				 // function for the total multigrid
				 // object), then create the solver object
				 // and solve the system. This is very
				 // easy, and we didn't even see any
				 // difference in the solve process
				 // compared to step-16. The magic is all
				 // hidden behind the implementation of
				 // the MatrixFree::vmult operation.
  const unsigned int multigrid_memory
    = (mg_matrices.memory_consumption() +
       mg_transfer.memory_consumption() +
       coarse_matrix.memory_consumption());
  std::cout << "Multigrid objects memory consumption: "
	    << multigrid_memory/std::pow(2.,20.)
	    << " MBytes."
	    << std::endl;

  SolverControl           solver_control (1000, 1e-12);
  SolverCG<>              cg (solver_control);

  cg.solve (system_matrix, solution, system_rhs,
	    preconditioner);

  std::cout << "Convergence in " << solver_control.last_step()
	    << " CG iterations." << std::endl;
}



				 // @sect4{LaplaceProblem::output_results}

				 // Here is the data output, which is a
				 // simplified version of step-5. We use the
				 // standard VTK output for each grid
				 // produced in the refinement process.
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



				 // @sect4{LaplaceProblem::run}

				 // The function that runs the program is
				 // very similar to the one in step-16. We
				 // make less refinement steps in 3D
				 // compared to 2D, but that's it.
template <int dim>
void LaplaceProblem<dim>::run ()
{
  for (unsigned int cycle=0; cycle<8-dim; ++cycle)
    {
      std::cout << "Cycle " << cycle << std::endl;

      if (cycle == 0)
	{
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



				 // @sect3{The <code>main</code> function}

				 // This is as in all other programs:
int main ()
{
  try
    {
      deallog.depth_console (0);
      LaplaceProblem<2> laplace_problem (2);
      laplace_problem.run ();
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
