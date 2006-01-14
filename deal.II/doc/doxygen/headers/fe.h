//-------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-------------------------------------------------------------------------

/**
 * @defgroup feall Finite elements
 *
 * All classes related to shape functions and to access to shape
 * functions.  This concerns the actual values of finite elements. For
 * the numbering of degrees of freedom refer to the module on @ref dofs.
 *
 * The classes and functions of this module fall into several sub-groups that
 * are discussed in their respective sub-modules listed above. In addition,
 * the FETools class provides functions that provide information on finite
 * elements, transformations between elements, etc.
 */


/**
 * @defgroup febase Base classes
 *
 * The members of this sub-module describe the implementational mechanics of
 * finite element classes, without actually implementing a concrete
 * element. For example, the FiniteElement base class declares the virtual
 * functions a derived class has to implement if it wants to describe a finite
 * element space. Likewise, the FiniteElementData holds variables that
 * describe certain values characterizing a finite element, such as the number
 * of degrees of freedom per vertex, line, or face.
 *
 * On the other hand, classes like FE_Poly and FE_PolyTensor are higher
 * abstractions. They describe finite elements that are built atop polynomial
 * descriptions of the shape functions on the unit cell. Classes derived from
 * them then only have to provide a description of the particular polynomial
 * from which a finite element is built. For example, the FE_Q class that
 * implements the usual Lagrange elements uses the FE_Poly base class to
 * generate a finite element by providing it with a set of Lagrange
 * interpolation polynomials corresponding to an equidistant subdivision of
 * interpolation points.
 *
 * Finally, the FESystem class is used for vector-valued problems. There, one
 * may want to couple a number of scalar (or also vector-valued) base elements
 * together to form the joint finite element of a vector-valued operator. As
 * an example, for 3d Navier-Stokes flow, one may want to use three Q1
 * elements for the three components of the velocity, and a piecewise constant
 * Q0 element for the pressure. The FESystem class can be used to couple these
 * four base elements together into a single, vector-valued element with 4
 * vector components. The step-8, step-17, and step-18 tutorial programs give
 * an introduction into the use of this class in the context of the
 * vector-valued elasticity (Lam&eacute;) equations. step-20 discusses a mixed
 * Laplace discretization that also uses vector-valued elements.
 * 
 * @ingroup feall
 */


/**
 * @defgroup feaccess Finite element access/FEValues classes
 *
 * @ingroup feall
 */


/**
 * @defgroup fe Finite element space descriptions
 *
 * The classes here describe finite element spaces, such as the simplest Q1
 * (bi-/trilinear) spaces, and higher order Lagrangian spaces Qp, but also
 * more specialized spaces such as Nedelec or Raviart-Thomas ones. Concrete
 * implementations are derived from the abstract FiniteElement base class.
 *
 * In essence, the functions these classes have to implement provide the
 * ability to query the value or derivatives of a shape function at a given
 * point on the unit cell. To be useful in integrating matrix and right hand
 * side entries, one has to have the ability to map these shape functions and
 * gradients to the real cell. This is done using classes derived from the
 * Mapping base class (see the @ref mapping module) in conjunction with the
 * FEValues class (see the @ref feaccess module).
 *
 * The FESystem class is different since it doesn't describe shape functions
 * itself, but assembles a vector-valued finite element from other finite
 * element objects. This functionality is described step-8, step-17 and other
 * tutorial programs after that.
 * 
 * @ingroup feall
 */


/**
 * @defgroup mapping Mappings between reference and real cell
 *
 * The classes in this module are used to map from unit coordinates to the
 * coordinates of a cell in real cell. Most commonly, one uses the MappingQ1
 * class that provides a Q1 (bi-/trilinear) mapping (i.e. a mapping that is
 * isoparametric for the usual Q1 elements). However, there are other classes
 * that implement higher-order mappings as well to provide for curvilinear
 * elements. These are discussed in the step-11 and step-12 tutorial programs.
 *
 * The MappingQ1Eulerian class is an extension to the MappingQ1 class in that
 * it accepts a vector that describes a displacement field for each position
 * of the domain. This is used in Eulerian computations without the need to
 * actually move vertices after each time step.
 * 
 * In addition, the MappingC1 class provides for a boundary of the
 * computational domain that is not only curved, but also has a continuous
 * derivative at the interface between two cells on the boundary.
 * 
 * Finally, the MappingCartesian class is an optimization for elements that
 * are brick-shaped and with edges parallel to the coordinate axes.
 * 
 * @ingroup feall
 */
