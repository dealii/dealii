/*----------------------------   problem_assembler.h     ---------------------------*/
/*      $Id$                 */
#ifndef __problem_assembler_H
#define __problem_assembler_H
/*----------------------------   problem_assembler.h     ---------------------------*/

#include <grid/dof.h>
#include <grid/dof_accessor.h>
#include <fe/fe.h>
#include <base/exceptions.h>
#include <vector.h>


// forward declarations
class dFMatrix;
class dVector;
template <int dim> class ProblemBase;




/**
  This is the base class for equation objects. Equations objects describe the
  finite element discretisation of one or more equations.

  Equation objects need only provide functions which set up the cell
  matrices and the right hand side(s). These are then automatically inserted
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
				      * cell matrix and a specific (user
				      * selectable) number of right hand sides
				      * on a given cell.
				      */
    virtual void assemble (dFMatrix            &cell_matrix,
			   vector<dVector>     &rhs,
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
				      * the right hand side(s) on a given cell.
				      */
    virtual void assemble (vector<dVector>     &rhs,
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
  An #Assembler# is a specialized version of a #DoFCellAccessor# which adds
  functionality to assemble global matrices and vectors from cell base ones.
  */
template <int dim>
class Assembler : public DoFCellAccessor<dim> {
  public:
				     /**
				      * Structure to be passed upon
				      * construction of an assembler object.
				      */
    struct AssemblerData {
					 /**
					  * Pointer to the dof handler object
					  * to be used to iterate on.
					  */
	DoFHandler<dim>  *dof;
					 /**
					  * Flags whether to assemble the matrix
					  * and right hand sides.
					  */
	bool              assemble_matrix, assemble_rhs;
					 /**
					  * How many right hand sides are there.
					  */
	unsigned int      n_rhs;
					 /**
					  * Pointer to the #ProblemBase# object
					  * of which the global matrices and
					  * vectors are to be assembled.
					  */
	ProblemBase<dim> *problem;
					 /**
					  * Pointer to a quadrature object to be
					  * used for this assemblage process.
					  */
	Quadrature<dim>  *quadrature;
    };

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
				      * Have room for a certain number of
				      * right hand sides.
				      */
    vector<dVector>   cell_vectors;

				     /**
				      * Store whether to assemble the
				      * global matrix.
				      */
    bool              assemble_matrix;

				     /**
				      * Store whether to assemble the
				      * right hand sides.
				      */
    bool              assemble_rhs;

				     /**
				      * Pointer to the problem class object
				      * of which we are to use the matrices
				      * and vectors.
				      */
    ProblemBase<dim> *problem;

				     /**
				      * The finite element evaluated at the
				      * quadrature points.
				      */
    FEValues<dim>     fe_values;
};

    
    


/*----------------------------   problem_assembler.h     ---------------------------*/
/* end of #ifndef __problem_assembler_H */
#endif
/*----------------------------   problem_assembler.h     ---------------------------*/
