//----------------------------  fe_lib.dg.h  ---------------------------
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
//----------------------------  fe_lib.dg.h  ---------------------------
#ifndef __deal2__fe_lib_dg_h
#define __deal2__fe_lib_dg_h


#include <fe/fe_lib.lagrange.h>

#define FEDGConstant FEDG_Q0
#define FEDGLinear FEDG_Q1
#define FEDGQuadraticSub FEDG_Q2
#define FEDGCubicSub FEDG_Q3
#define FEDGQuarticSub FEDG_Q4

/**
 * Define a constant discontinuous finite element in #dim#
 * space dimensions, along with (bi-, tri-)linear (therefore isoparametric)
 * transforms from the unit cell to the real cell.
 * @author Ralf Hartmann, 1998
 */
template <int dim>
class FEDG_Q0 : public FEQ1Mapping<dim> {
  public:
				     /**
				      * Constructor
				      */
    FEDG_Q0 ();
    
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
};


/**
 * Define a (bi-, tri-, etc)linear finite element in #dim# space dimensions,
 * along with (bi-, tri-)linear (therefore isoparametric) transforms from the
 * unit cell to the real cell allowing discontinuous Galerkin methods.
 *
 * 
 * This class is derived from and provides substantially the same 
 * as the #FEQ1# class. The only difference is the new constructor that
 * calls #FEQ1::FEQ1(const int)#, the protected constructor of the
 * #FEQ1# class using a #FiniteElement# with no dofs in the vertices and 
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
class FEDG_Q1 : public FEQ1<dim>{
  public:
				     /**
				      * Constructor
				      */
    FEDG_Q1();

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;
};


/**
 * Define a (bi-, tri-, etc)quadratic finite element in #dim# space dimensions,
 * along with (bi-, tri-)quadratic (therefore isoparametric) transforms from the
 * unit cell to the real cell allowing discontinuous Galerkin methods.
 * 
 * This class is derived from and provides substantially the same 
 * as the #FEQ2# class. The only difference is the new constructor that
 * calls #FEQ2::FEQ2(const int)#, the protected constructor of the
 * #FEQ2# class using a #FiniteElement# with no dofs in the vertices, no dofs on the lines and 
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
class FEDG_Q2 : public FEQ2<dim>{
  public:
				     /**
				      * Constructor
				      */
    FEDG_Q2();
				 
				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;
};


/**
 * Define a (bi-, tri-, etc)cubic finite element in #dim# space dimensions,
 * along with (bi-, tri-)cubic (therefore isoparametric) transforms from the
 * unit cell to the real cell allowing discontinuous Galerkin methods.
 * 
 * This class is derived from and provides substantially the same 
 * as the #FEQ3# class. The only difference is the new constructor that
 * calls #FEQ3::FEQ3(const int)#, the protected constructor of the
 * #FEQ3# class using a #FiniteElement# with no dofs in the vertices, no dofs on the lines and 
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
class FEDG_Q3 : public FEQ3<dim>{
  public:
				     /**
				      * Constructor
				      */
    FEDG_Q3();

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;
};


/**
 * Define a (bi-, tri-, etc)quartic finite element in #dim# space dimensions,
 * along with (bi-, tri-)quartic (therefore isoparametric) transforms from the
 * unit cell to the real cell allowing discontinuous Galerkin methods.
 * 
 * This class is derived from and provides substantially the same 
 * as the #FEQ4# class. The only difference is the new constructor that
 * calls #FEQ4::FEQ4(const int)#, the protected constructor of the
 * #FEQ4# class using a #FiniteElement# with no dofs in the vertices, no dofs on the lines and 
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
class FEDG_Q4 : public FEQ4<dim>{
  public:
				     /**
				      * Constructor
				      */
    FEDG_Q4();

				     /**
				      * Refer to the base class for detailed
				      * information on this function.
				      */
    virtual void get_face_support_points (const DoFHandler<dim>::face_iterator &face,
					  vector<Point<dim> > &support_points) const;
};


#endif
