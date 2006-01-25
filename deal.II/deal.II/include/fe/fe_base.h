//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_base_h
#define __deal2__fe_base_h

#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/point.h>
#include <base/tensor.h>
#include <base/table.h>
#include <base/vector_slice.h>
#include <base/geometry_info.h>
#include <lac/full_matrix.h>
#include <fe/fe_update_flags.h>
#include <fe/mapping.h>

#include <string>
#include <vector>

template<int dim> class FESystem;

/**
 * Dimension independent data for finite elements. See the derived
 * class FiniteElement class for information on its use. All
 * its data are available to the implementation in a concrete finite
 * element class.
 *
 * Remark on a change in implementation: it is now wrong to cast a
 * pointer to FiniteElement to a pointer to FiniteElementData and
 * delete it. The virtual destructor has been moved up. In a later
 * version, FiniteElementData and FiniteElement should be private
 * base classes of FiniteElement.
 *
 * @ingroup febase
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999, 2000, 2001, 2003, 2005
 */
template <int dim>
class FiniteElementData
{
  public:
				     /**
				      * Enumerator for the different
				      * types of continuity a finite
				      * element may have. Continuity
				      * is measured by the Sobolev
				      * space containing the
				      * constructed finite element
				      * space and is also called this
				      * way.
				      *
				      * Note that certain continuities
				      * may imply others. For
				      * instance, a function in
				      * <i>H<sup>1</sup></i> is in
				      * <i>H<sup>curl</sup></i> and
				      * <i>H<sup>div</sup></i> as
				      * well.
				      *
				      * If you are interested in
				      * continuity in the classical
				      * sense, then the following
				      * relations hold:
				      *
				      * <ol>
				      *
				      * <li> <i>H<sup>1</sup></i>
				      * implies that the function is
				      * continuous over cell
				      * boundaries.
				      *
				      * <li> <i>H<sup>2</sup></i>
				      * implies that the function is
				      * continuously differentiable
				      * over cell boundaries.
				      *
				      * <li> <i>L<sup>2</sup></i>
				      * indicates that the element is
				      * discontinuous. Since
				      * discontinuous elements have no
				      * topological couplings between
				      * grid cells and code may
				      * actually depend on this
				      * property, <i>L<sup>2</sup></i>
				      * conformity is handled in a
				      * special way in the sense that
				      * it is <b>not</b> implied by
				      * any higher conformity.
				      * </ol>
				      *
				      * In order to test if a finite
				      * element conforms to a certain
				      * space, use
				      * FiniteElementData<dim>::conforms().
				      */
    enum Conformity
    {
					   /**
					    * Indicates incompatible
					    * continuities of a
					    * system.
					    */
	  unknown = 0x00,
	  
					   /**
					    * Discontinuous
					    * elements. See above!
					    */
	  L2 = 0x01,
	  
					   /**
					    * Conformity with the
					    * space
					    * <i>H<sup>curl</sup></i>
					    * (continuous tangential
					    * component of a vector
					    * field)
					    */
	  Hcurl = 0x02,
	  
					   /**
					    * Conformity with the
					    * space
					    * <i>H<sup>div</sup></i>
					    * (continuous normal
					    * component of a vector
					    * field)
					    */
	  Hdiv = 0x04,
	  
					   /**
					    * Conformity with the
					    * space
					    * <i>H<sup>1</sup></i>
					    * (continuous)
					    */
	  H1 = Hcurl | Hdiv,
	  
					   /**
					    * Conformity with the
					    * space
					    * <i>H<sup>2</sup></i>
					    * (continuously
					    * differentiable)
					    */
	  H2 = 0x0e
    };
    
				     /**
				      * Number of degrees of freedom on
				      * a vertex.
				      */
    const unsigned int dofs_per_vertex;

				     /** Number of degrees of freedom
				      *  in a line; not including the
				      *  degrees of freedom on the
				      *  vertices of the line.
				      */
    const unsigned int dofs_per_line;

				     /** Number of degrees of freedom
				      *  in a quadrilateral; not
				      *  including the degrees of
				      *  freedom on the lines and
				      *  vertices of the
				      *  quadrilateral.
				      */
    const unsigned int dofs_per_quad;

				     /** Number of degrees of freedom
				      *  in a hexahedron; not
				      *  including the degrees of
				      *  freedom on the
				      *  quadrilaterals, lines and
				      *  vertices of the hecahedron.
				      */
    const unsigned int dofs_per_hex;

				     /**
				      * First index of dof on a line.
				      */
    const unsigned int first_line_index;
    
				     /**
				      * First index of dof on a quad.
				      */
    const unsigned int first_quad_index;
    
				     /**
				      * First index of dof on a hexahedron.
				      */
    const unsigned int first_hex_index;
    
				     /**
				      * First index of dof on a line for face data.
				      */
    const unsigned int first_face_line_index;
    
				     /**
				      * First index of dof on a quad for face data.
				      */
    const unsigned int first_face_quad_index;

				     /**
				      * Number of degrees of freedom
				      * on a face. This is the
				      * accumulated number of degrees
				      * of freedom on all the objects
				      * of dimension up to
				      * <tt>dim-1</tt> constituting a
				      * face.
				      */
    const unsigned int dofs_per_face;
    
				     /**
				      * Total number of degrees of freedom
				      * on a cell. This is the
				      * accumulated number of degrees
				      * of freedom on all the objects
				      * of dimension up to
				      * <tt>dim</tt> constituting a
				      * cell.
				      */
    const unsigned int dofs_per_cell;

				     /**
				      * Number of vector components of
				      * this finite element, and
				      * dimension of the image
				      * space. For vector-valued
				      * finite elements (i.e. when
				      * this number is greater than
				      * one), the number of vector
				      * components is in many cases
				      * equal to the number of base
				      * elements glued together with
				      * the help of the FESystem
				      * class. However, for elements
				      * like the Nedelec element, the
				      * number is greater than one
				      * even though we only have one
				      * base element.
				      */
    const unsigned int components;
				     /**
				      * The number of vector blocks a
				      * BlockVector for this element
				      * should contain. For primitive
				      * elements, this is equal to
				      * #components, but for vector
				      * valued base elements it may
				      * differ. Actually, in systems
				      * this is the sum of the base
				      * element multiplicities.
				      */
    const unsigned int blocks;
    
				     /**
				      * Maximal polynomial degree of a
				      * shape function in a single
				      * coordinate direction.
				      */
    const unsigned int degree;

				     /**
				      * Indicate the space this element conforms to.
				      */
    const Conformity conforming_space;
    
    				     /**
				      * Default
				      * constructor. Constructs an
				      * element with no dofs. Checking
				      * n_dofs_per_cell() is therefore
				      * a good way to check if
				      * something went wrong.
				      */
    FiniteElementData ();

				     /**
				      * Constructor, computing all
				      * necessary values from the
				      * distribution of dofs to
				      * geometrcal objects.
				      *
				      * @param dofs_per_object Number
				      * of dofs on geometrical objects
				      * for each dimension. In this
				      * vector, entry 0 refers to dofs
				      * on vertices, entry 1 on lines
				      * and so on. Its length must be
				      * <i>dim+1</i>.
				      * @param n_components Number of
				      * vector components of the
				      * element.
				      * @param degree
				      * Maximal polynomial degree in a
				      * single direction.
				      * @param conformity The finite
				      * element space has continuity
				      * of this Sobolev space.
				      * @param n_blocks Number of
				      * vector blocks.
				      */
    FiniteElementData (const std::vector<unsigned int> &dofs_per_object,
		       const unsigned int n_components,
		       const unsigned int degree,
		       const Conformity conformity = unknown,
		       const unsigned int n_blocks = deal_II_numbers::invalid_unsigned_int);

				     /**
				      * Number of dofs per vertex.
				      */
    unsigned int n_dofs_per_vertex () const;

				     /**
				      * Number of dofs per line. Not
				      * including dofs on lower
				      * dimensional objects.
				      */
    unsigned int n_dofs_per_line () const;

    				     /**
				      * Number of dofs per quad. Not
				      * including dofs on lower
				      * dimensional objects.
				      */
    unsigned int n_dofs_per_quad () const;

    				     /**
				      * Number of dofs per hex. Not
				      * including dofs on lower
				      * dimensional objects.
				      */
    unsigned int n_dofs_per_hex () const;

    				     /**
				      * Number of dofs per face,
				      * accumulating degrees of
				      * freedom of all lower
				      * dimensional objects.
				      */
    unsigned int n_dofs_per_face () const;

    				     /**
				      * Number of dofs per cell,
				      * accumulating degrees of
				      * freedom of all lower
				      * dimensional objects.
				      */
    unsigned int n_dofs_per_cell () const;

    				     /**
				      * Number of components.
				      */
    unsigned int n_components () const;

    				     /**
				      * Number of blocks.
				      */
    unsigned int n_blocks () const;

				     /**
				      * Maximal polynomial degree of a
				      * shape function in a single
				      * coordinate direction.
				      *
				      * This function can be used to
				      * determine the optimal
				      * quadrature rule.
				      */
    unsigned int tensor_degree () const;

				     /**
				      * Test whether a finite element
				      * space conforms to a certain
				      * Sobolev space.
				      *
				      * @note This function will
				      * return a true value even if
				      * the finite element space has
				      * higher regularity than asked
				      * for.
				      */
    bool conforms (const Conformity) const;
    
				     /**
				      * Comparison operator.
				      */
    bool operator == (const FiniteElementData &) const;
};


// --------- inline and template functions ---------------

template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_vertex () const
{
  return dofs_per_vertex;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_line () const
{
  return dofs_per_line;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_quad () const
{
  return dofs_per_quad;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_hex () const
{
  return dofs_per_hex;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_face () const
{
  return dofs_per_face;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_dofs_per_cell () const
{
  return dofs_per_cell;
}


template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_components () const
{
  return components;
}



template <int dim>
inline
unsigned int 
FiniteElementData<dim>::n_blocks () const
{
  return blocks;
}



template <int dim>
inline
unsigned int 
FiniteElementData<dim>::tensor_degree () const
{
  return degree;
}


template <int dim>
inline
bool
FiniteElementData<dim>::conforms (Conformity space) const
{
  return ((space & conforming_space) != 0);
}




#endif
