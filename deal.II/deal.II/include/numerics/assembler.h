//----------------------------  assembler.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  assembler.h  ---------------------------
#ifndef __deal2__assembler_h
#define __deal2__assembler_h


/*----------------------------   problem_assembler.h     ---------------------------*/

#include <base/exceptions.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_values.h>
#include <vector>


/**
 * The use of this class is now deprecated!
 *
 * This is the base class for equation objects. Equations objects describe the
 * finite element discretisation of one or more equations.
 *
 * Equation objects need only provide functions which set up the cell
 * matrices and the cell right hand side. These are then automatically inserted
 * into the global matrices and vectors.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class Equation
{
  public:
				     /**
				      * Constructor. You have to pass the number
				      * of equations you want to discretize, which
				      * equals the number of solution functions.
				      */
    Equation (const unsigned int n_equations);

				     /**
				      * Virtual function which assembles the
				      * cell matrix and the right hand side
				      * on a given cell.
				      *
				      * This function assumes the cell matrix
				      * and right hand side to have the right
				      * size and to be empty. Functions of
				      * derived classes should check for
				      * this.
				      * For that purpose, the two exceptions
				      * @p{ExcWrongSize} and @p{ExcObjectNotEmpty}
				      * are declared.
				      */
    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;

				     /**
				      * Virtual function which only assembles
				      * the cell matrix on a given cell.
				      *
				      * This function assumes the cell matrix
				      * and right hand side to have the right
				      * size and to be empty. Functions of
				      * derived classes should check for
				      * this.
				      * For that purpose, the two exceptions
				      * @p{ExcWrongSize} and @p{ExcObjectNotEmpty}
				      * are declared.
				      */
    virtual void assemble (FullMatrix<double>  &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;

				     /**
				      * Virtual function which only assembles
				      * the right hand side on a given cell.
				      *
				      * This function assumes the cell matrix
				      * and right hand side to have the right
				      * size and to be empty. Functions of
				      * derived classes should check for
				      * this.
				      * For that purpose, the two exceptions
				      * @p{ExcWrongSize} and @p{ExcObjectNotEmpty}
				      * are declared.
				      */
    virtual void assemble (Vector<double>      &rhs,
			   const FEValues<dim> &fe_values,
			   const DoFHandler<dim>::cell_iterator &cell) const;

				     /**
				      * Return number of equations for this
				      * equation object. This equals the number
				      * of solution functions.
				      */
    unsigned int n_equations () const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcPureVirtualFunctionCalled);
				     /**
				      * Exception
				      */
    DeclException2 (ExcWrongSize,
		    int, int,
		    << "Object has wrong size " << arg1
		    << ", but should have " << arg2 << ".");
				     /**
				      * Exception
				      */
    DeclException0 (ExcObjectNotEmpty);
    
  protected:
				     /**
				      * Store the number of solution functions,
				      * which is the same as the number of
				      * equations.
				      */
    const unsigned int n_eq;
};



/**
 * The use of this class is now deprecated!
 *
 * An @p{Assembler} is a specialized version of a @p{DoFCellAccessor} which adds
 * functionality to assemble global matrices and vectors from cell base ones.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class Assembler : public DoFCellAccessor<dim>
{
  public:

				     /**
				      * Structure to be passed upon
				      * construction of an assembler object
				      * through the iterator object. See
				      * \Ref{TriaRawIterator} for a discussion
				      * of this mechanism.
				      */
    struct AssemblerData {
					 /**
					  * Constructor.
					  */
	AssemblerData (const DoFHandler<dim>    &dof,
		       const bool                assemble_matrix,
		       const bool                assemble_rhs,
		       SparseMatrix<double>     &matrix,
		       Vector<double>           &rhs_vector,
		       const Quadrature<dim>    &quadrature,
		       const UpdateFlags        &update_flags);
	
					 /**
					  * Pointer to the dof handler object
					  * to be used to iterate on.
					  */
	const DoFHandler<dim>  &dof;
	
					 /**
					  * Flags to assemble the matrix.
					  */
	const bool              assemble_matrix;

					 /**
					  * Flags whether to assemble the right hand sides.
					  */
	const bool              assemble_rhs;
	
					 /**
					  * Pointer to the matrix to be assembled
					  * by this object. Elements are summed
					  * up by the assembler, so you may want
					  * to clear this object (set all entries
					  * to zero) before use.
					  */
	SparseMatrix<double>   &matrix;
	
					 /**
					  * Pointer to the vector to be assembled
					  * by this object. Elements are summed
					  * up by the assembler, so you may want
					  * to clear this object (set all entries
					  * to zero) before use.
					  */
	Vector<double>         &rhs_vector;
	
					 /**
					  * Pointer to a quadrature object to be
					  * used for this assemblage process.
					  */
	const Quadrature<dim>  &quadrature;
	
					 /**
					  * Store which of the fields of the
					  * FEValues object need to be reinitialized
					  * on each cell.
					  */
	const UpdateFlags       update_flags;
    };


				     /**
				      * Declare the data type that this accessor
				      * class expects to get passed from the
				      * iterator classes.
				      */
    typedef AssemblerData AccessorData;
    
				     /**
				      * Default constructor, unused thus not
				      * implemented.
				      */
    Assembler ();
    
				     /**
				      * Constructor. The @p{local_data}
				      * argument is assumed to be a pointer
				      * to an @p{AssemblerData} object. The data
				      * is copied, so the object need not live
				      * longer than the constructor call.
				      */
    Assembler (Triangulation<dim> *tria,
	       const int           level,
	       const int           index,
	       const AccessorData *local_data);
    
				     /**
				      * Assemble on the present cell using
				      * the given equation objectand the data
				      * passed to the constructor. The elements
				      * of the local matrix and right hand side
				      * are added to the global matrix and
				      * vector so you may want to clear the
				      * matrix before use.
				      */
    void assemble (const Equation<dim> &);

				     /**
				      * Exception.
				      */
    DeclException0 (ExcNoAssemblingRequired);
				     /**
				      * Exception.
				      */
    DeclException0 (ExcInvalidData);
				     /**
				      * Exception.
				      */
				     /**
				      * Exception.
				      */
  private:
				     /**
				      * Store a local cell matrix.
				      */
    FullMatrix<double>  cell_matrix;

				     /**
				      * Right hand side local to cell.
				      */
    Vector<double>    cell_vector;

				     /**
				      * Store whether to assemble the
				      * global matrix.
				      */
    bool              assemble_matrix;

				     /**
				      * Store whether to assemble the
				      * right hand side.
				      */
    bool              assemble_rhs;

				     /**
				      * Pointer to the matrix to be assembled
				      * by this object.
				      */
    SparseMatrix<double>    &matrix;

				     /**
				      * Pointer to the vector to be assembled
				      * by this object.
				      */
    Vector<double>          &rhs_vector;

				     /**
				      * The finite element evaluated at the
				      * quadrature points.
				      */
    FEValues<dim>     fe_values;
};


/*----------------------------   problem_assembler.h     ---------------------------*/

#endif
/*----------------------------   problem_assembler.h     ---------------------------*/
