/*----------------------------   problem_assembler.h     ---------------------------*/
/*      $Id$                 */
#ifndef __problem_assembler_H
#define __problem_assembler_H
/*----------------------------   problem_assembler.h     ---------------------------*/

#include <base/exceptions.h>
#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <fe/fe_values.h>
#include <vector>


// forward declarations
template <int dim> class FiniteElement;
template <int dim> class Quadrature;

class dFMatrix;
class dVector;




/**
  This is the base class for equation objects. Equations objects describe the
  finite element discretisation of one or more equations.

  Equation objects need only provide functions which set up the cell
  matrices and the cell right hand side. These are then automatically inserted
  into the global matrices and vectors.
  */
template <int dim>
class Equation {
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
				      */
    virtual void assemble (dFMatrix            &cell_matrix,
			   dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const Triangulation<dim>::cell_iterator &cell) const;

				     /**
				      * Virtual function which only assembles
				      * the cell matrix on a given cell.
				      */
    virtual void assemble (dFMatrix            &cell_matrix,
			   const FEValues<dim> &fe_values,
			   const Triangulation<dim>::cell_iterator &cell) const;

				     /**
				      * Virtual function which only assembles
				      * the right hand side on a given cell.
				      */
    virtual void assemble (dVector             &rhs,
			   const FEValues<dim> &fe_values,
			   const Triangulation<dim>::cell_iterator &cell) const;

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
    
  protected:
				     /**
				      * Store the number of solution functions,
				      * which is the same as the number of
				      * equations.
				      */
    const unsigned int n_eq;
};





/**
 * Structure to be passed upon
 * construction of an assembler object
 * through the iterator object. See
 * \Ref{TriaRawIterator} for a discussion
 * of this mechanism.
 */
template <int dim>
struct AssemblerData {
				     /**
				      * Constructor.
				      */
    AssemblerData (const DoFHandler<dim>    &dof,
		   const bool                assemble_matrix,
		   const bool                assemble_rhs,
		   dSMatrix                 &matrix,
		   dVector                  &rhs_vector,
		   const Quadrature<dim>    &quadrature,
		   const FiniteElement<dim> &fe,
		   const UpdateFields       &update_flags,
		   const Boundary<dim>      &boundary);
    
				     /**
				      * Pointer to the dof handler object
				      * to be used to iterate on.
				      */
    const DoFHandler<dim>  &dof;
    
				     /**
				      * Flags whether to assemble the matrix
				      * and right hand sides.
				      */
    const bool              assemble_matrix, assemble_rhs;

				     /**
				      * Pointer to the matrix to be assembled
				      * by this object.
				      */
    dSMatrix               &matrix;

				     /**
				      * Pointer to the vector to be assembled
				      * by this object.
				      */
    dVector                &rhs_vector;
    
				     /**
				      * Pointer to a quadrature object to be
				      * used for this assemblage process.
				      */
    const Quadrature<dim>  &quadrature;
    
				     /**
				      * Use this FE type for the assemblage
				      * process. The FE object must match that
				      * used to construct the #DoFHandler#
				      * object.
				      */
    const FiniteElement<dim> &fe;
    
				     /**
				      * Store which of the fields of the
				      * FEValues object need to be reinitialized
				      * on each cell.
				      */
    const UpdateFields       update_flags;

				     /**
				      * Store a pointer to the object describing
				      * the boundary of the domain. This is
				      * necessary, since we may want to use
				      * curved faces of cells at the boundary
				      * when using higher order elements.
				      */
    const Boundary<dim>     &boundary;
};




/**
  An #Assembler# is a specialized version of a #DoFCellAccessor# which adds
  functionality to assemble global matrices and vectors from cell base ones.
  */
template <int dim>
class Assembler : public DoFCellAccessor<dim> {
  public:
				     /**
				      * Default constructor, unused thus not
				      * implemented.
				      */
    Assembler ();
    
				     /**
				      * Constructor. The #local_data#
				      * argument is assumed to be a pointer
				      * to an #AssemblerData# object. The data
				      * is copied, so the object need not live
				      * longer than the constructor call.
				      */
    Assembler (Triangulation<dim> *tria,
	       const int           level,
	       const int           index,
	       const void         *local_data);
    
				     /**
				      * Assemble on the present cell using
				      * the given equation object.
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
    dFMatrix          cell_matrix;

				     /**
				      * Right hand side local to cell.
				      */
    dVector           cell_vector;

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
    dSMatrix         &matrix;

				     /**
				      * Pointer to the vector to be assembled
				      * by this object.
				      */
    dVector          &rhs_vector;

				     /**
				      * Pointer to the finite element used for
				      * the assemblage process. We store a
				      * pointer to the finite element, not a
				      * copy, but this is not a problem, since
				      * #Assembler# objects are usually created
				      * and deleted by the
				      * #ProblemBase::assemble# routine, which
				      * itself gets passed the finite element;
				      * the lifetime of the #Assembler# object
				      * is thus less than the lifetime of the
				      * finite element object and storing
				      * a pointer is not risky.
				      */
    const FiniteElement<dim> &fe;

				     /**
				      * The finite element evaluated at the
				      * quadrature points.
				      */
    FEValues<dim>     fe_values;

				     /**
				      * Store a pointer to the object describing
				      * the boundary of the domain. This is
				      * necessary, since we may want to use
				      * curved faces of cells at the boundary
				      * when using higher order elements.
				      */
    const Boundary<dim> &boundary;
};

    
    


/*----------------------------   problem_assembler.h     ---------------------------*/
/* end of #ifndef __problem_assembler_H */
#endif
/*----------------------------   problem_assembler.h     ---------------------------*/
