/*----------------------------   fe_lib.dg.h     ---------------------------*/
/*      $Id$ */
/*   Ralf Hartmann, University of Heidelberg, Dez 98                        */
#ifndef __fe_lib_dg_H 
#define __fe_lib_dg_H 
/*----------------------------   fe_lib.dg.h     ---------------------------*/


#include <fe/fe_lib.lagrange.h>


/**
 * Define a constant discontinuous finite element in #dim#
 * space dimensions, along with (bi-, tri-)linear (therefore isoparametric)
 * transforms from the unit cell to the real cell.
 * @author Ralf Hartmann, 1998
 */
template <int dim>
class FEDGConstant : public FELinearMapping<dim> {
  public:
				     /**
				      * Constructor
				      */
    FEDGConstant ();
    
				     /**
				      * Return the value of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual double shape_value(const unsigned int i,
			       const Point<dim>& p) const;
    
				     /**
				      * Return the gradient of the #i#th shape
				      * function at point #p# on the unit cell.
				      */
    virtual Tensor<1,dim> shape_grad(const unsigned int i,
				     const Point<dim>& p) const;

				     /**
				      * Return the tensor of second derivatives
				      * of the #i#th shape function at
				      * point #p# on the unit cell.
				      *
				      * For linear elements, all second
				      * derivatives on the unit cell are zero.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim>   &p) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_unit_support_points (vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_support_points (const DoFHandler<dim>::cell_iterator &cell,
				     vector<Point<dim> > &support_points) const;

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;

    				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_local_mass_matrix (const DoFHandler<dim>::cell_iterator &cell,
					FullMatrix<double> &local_mass_matrix) const;
  
				     /**
				      * Return a readonly reference to the
				      * matrix which describes the transfer of a
				      * child with the given number to the
				      * mother cell. See the #restriction# array
				      * for more information.
				      *
				      * This function returns an error since the
				      * correct use of the restriction
				      * matrices is not yet finally decided
				      * about.
				      */
    const FullMatrix<double> & restrict (const unsigned int) const;
  
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotImplemented);
};



/**
 * Define a (bi-, tri-, etc)linear finite element in #dim# space dimensions,
 * along with (bi-, tri-)linear (therefore isoparametric) transforms from the
 * unit cell to the real cell allowing discontinuous Galerkin methods.
 *
 * 
 * This class is derived from and provides substantially the same 
 * as the #FELinear# class. The only difference is the new constructor that
 * calls #FELinear::FELinear(const int)#, the protected constructor of the
 * #FELinear# class using a #FiniteElement# with no dofs in the vertices and 
 * $2^d$ dofs per cell. As now the cells do not share any vertex-dof with
 * a neighboring cell the $2^d$ dofs per cell can be choosen independently not
 * needing any constraints and allowing the use of discontinuous Galerkin
 * methods. Although the basis functions now are not longer associated 
 * with the vertices but with the cell they retain their shape. As already
 * explained no constraint matrices needed to be implemented.
 * To use this element you need to think about the jump terms in your 
 * weak formulation of your discontinuous Galerkin scheme.
 * @author Ralf Hartmann, 1998
 */
template <int dim>
class FEDGLinear : public FELinear<dim>{
  public:
				     /**
				      * Constructor
				      */
    FEDGLinear();

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;


				     /**
				      * This function returns an error since the
				      * correct use of the restriction
				      * matrices is not yet finally decided
				      * about.
				      */
    const FullMatrix<double> & restrict (const unsigned int) const;
};



/**
 * Define a (bi-, tri-, etc)quadratic finite element in #dim# space dimensions,
 * along with (bi-, tri-)quadratic (therefore isoparametric) transforms from the
 * unit cell to the real cell allowing discontinuous Galerkin methods.
 * 
 * This class is derived from and provides substantially the same 
 * as the #FEQuadraticSub# class. The only difference is the new constructor that
 * calls #FEQuadraticSub::FEQuadraticSub(const int)#, the protected constructor of the
 * #FEQuadraticSub# class using a #FiniteElement# with no dofs in the vertices, no dofs on the lines and 
 * $3^d$ dofs per cell. As now the cells do not share any vertex-dof with
 * a neighboring cell the $3^d$ dofs per cell can be choosen independently not
 * needing any constraints and allowing the use of discontinuous Galerkin
 * methods. Although the basis functions now are not longer associated 
 * with the vertices but with the cell they retain their shape. As already
 * explained no constraint matrices needed to be implemented.
 * To use this element you need to think about the jump terms in your 
 * weak formulation of your discontinuous Galerkin scheme.
 * @author Ralf Hartmann, 1998
 */
template <int dim>
class FEDGQuadraticSub : public FEQuadraticSub<dim>{
  public:
				     /**
				      * Constructor
				      */
    FEDGQuadraticSub();
				 
				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;

				     /**
				      * This function returns an error since the
				      * correct use of the restriction
				      * matrices is not yet finally decided
				      * about.
				      */
    const FullMatrix<double> & restrict (const unsigned int) const;
};




/**
 * Define a (bi-, tri-, etc)cubic finite element in #dim# space dimensions,
 * along with (bi-, tri-)cubic (therefore isoparametric) transforms from the
 * unit cell to the real cell allowing discontinuous Galerkin methods.
 * 
 * This class is derived from and provides substantially the same 
 * as the #FECubicSub# class. The only difference is the new constructor that
 * calls #FECubicSub::FECubicSub(const int)#, the protected constructor of the
 * #FECubicSub# class using a #FiniteElement# with no dofs in the vertices, no dofs on the lines and 
 * $4^d$ dofs per cell. As now the cells do not share any vertex-dof with
 * a neighboring cell the $4^d$ dofs per cell can be choosen independently not
 * needing any constraints and allowing the use of discontinuous Galerkin
 * methods. Although the basis functions now are not longer associated 
 * with the vertices but with the cell they retain their shape. As already
 * explained no constraint matrices needed to be implemented.
 * To use this element you need to think about the jump terms in your 
 * weak formulation of your discontinuous Galerkin scheme.
 * @author Ralf Hartmann, 1998
 */
template <int dim>
class FEDGCubicSub : public FECubicSub<dim>{
  public:
				     /**
				      * Constructor
				      */
    FEDGCubicSub();

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;

				     /**
				      * This function returns an error since the
				      * correct use of the restriction
				      * matrices is not yet finally decided
				      * about.
				      */
    const FullMatrix<double> & restrict (const unsigned int) const;
};



/**
 * Define a (bi-, tri-, etc)quartic finite element in #dim# space dimensions,
 * along with (bi-, tri-)quartic (therefore isoparametric) transforms from the
 * unit cell to the real cell allowing discontinuous Galerkin methods.
 * 
 * This class is derived from and provides substantially the same 
 * as the #FEQuarticSub# class. The only difference is the new constructor that
 * calls #FEQuarticSub::FEQuarticSub(const int)#, the protected constructor of the
 * #FEQuarticSub# class using a #FiniteElement# with no dofs in the vertices, no dofs on the lines and 
 * $5^d$ dofs per cell. As now the cells do not share any vertex-dof with
 * a neighboring cell the $5^d$ dofs per cell can be choosen independently not
 * needing any constraints and allowing the use of discontinuous Galerkin
 * methods. Although the basis functions now are not longer associated 
 * with the vertices but with the cell they retain their shape. As already
 * explained no constraint matrices needed to be implemented.
 * To use this element you need to think about the jump terms in your 
 * weak formulation of your discontinuous Galerkin scheme.
 * @author Ralf Hartmann, 1998
 */
template <int dim>
class FEDGQuarticSub : public FEQuarticSub<dim>{
  public:
				     /**
				      * Constructor
				      */
    FEDGQuarticSub();

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;
				     
				     /**
				      * This function returns an error since the
				      * correct use of the restriction
				      * matrices is not yet finally decided
				      * about.
				      */
    const FullMatrix<double> & restrict (const unsigned int) const;
};





/*----------------------------   fe_lib.dg.h     ---------------------------*/
/* end of #ifndef __fe_lib_dg_H */
#endif
/*----------------------------   fe_lib.dg.h     ---------------------------*/
