// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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
 * The classes in this module are used when one wants to assemble matrices or
 * vectors. They link finite elements, quadrature objects, and mappins: the
 * finite element classes describe a finite element space on a unit cell
 * (i.e. the unit line segment, square, or cube <tt>[0,1]^d</tt>), the
 * quadrature classes describe where quadrature points are located and what
 * weight they have, and the mapping classes describe how to map a point from
 * the unit cell to a real cell and back. Since integration happens at
 * quadrature points on the real cell, and needs to know their location as
 * well as the values and gradients of finite element shape functions at these
 * points. The FEValues class coordinates getting this information. For
 * integrations on faces (for example for integration on the boundary, or
 * interfaces between cells), the FEFaceValues class offers similar
 * functionality as the FEValues class does for cells. Finally, the
 * FESubfaceValues class offers the possibility to ingrate on parts of faces
 * if the neighboring cell is refined and the present cell shares only a part
 * of its face with the neighboring cell. If vector-valued elements are used,
 * the FEValues and related classes allow access to all vector components; if
 * one wants to pick individual components, there are extractor classes that
 * make this task simpler, as described in the @ref vector_valued module.
 *
 * The last member of this group, the UpdateFlags enumeration, is used as an
 * optimization: instead of letting the FEValues class compute every possible
 * piece of data relating to a given finite element on a cell, you have to
 * specify up front which information you are actually interested in. The
 * UpdateFlags enumeration is used to offer symbolic names denoting what you
 * want the FEValues class to compute.
 * 
 * All these classes are used in all @ref Tutorial "tutorial programs" from
 * step-3 onward, and are described there in significant detail.
 *
 * The actual workings of the FEValues class and friends is
 * complicated because it has to be general yet efficient. The page on
 * @ref UpdateFlagsEssay attempts to give an overview of how this
 * works.
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
 * <h3>Vector-valued finite elements</h3>
 *
 * deal.II provides two different kinds of vector valued
 * elements. First, there is a group of genuine vector elements,
 * usually distinguished by the fact, that each vector component
 * consists of a different set of anisotropic polynomials. These
 * elements are typically associated with differential
 * forms. Currently, they are
 *
 * <ul>
 * <li> FE_ABF
 * <li> FE_BDM, FE_DGBDM
 * <li> FE_Nedelec, FE_DGNedelec
 * <li> FE_RaviartThomas, FE_DGRaviartThomas
 * </ul>
 *
 * Additionally, deal.II offers a mechanism to create a vector element
 * from existing scalar or vector elements. The FESystem class is
 * responsible for this: it doesn't describe shape functions itself, but
 * assembles a vector-valued finite element from other finite element
 * objects. This functionality is described step-8, step-17 and other
 * tutorial programs after that.
 *
 * @note Support  for the implementation of  vector-valued elements is
 * provided  by  the  class  FE_PolyTensor. Typically,  a  new  vector
 * element should be derived from this class.
 *
 * <h3>Discontinuous Galerkin</h3>
 *
 * For each finite element conforming to any space of weakly
 * differentiable functions like <i>H<sup>1</sup></i> or
 * <i>H<sup>curl</sup></i>, we can define an analogue DG space by
 * simply assigning all degrees of freedom on vertices, edges or faces
 * to the interior of the cell. This is to be understood in the
 * topological sense. The interpolation operator for such a degree of
 * freedom would still be on the boundary.  While not done so
 * consistently, we provide quite a few of these elements, plus those,
 * which have no conforming counterparts, like FE_DGP. Here is a list of the current DG elements:
 * <ul>
 * <li> scalar: FE_DGP, FE_DGQ
 * <li> scalar, different shape functions: FE_DGPMonomial, FE_DGPNonparametric, FE_DGQArbitraryNodes
 * <li> vector-valued:  FE_DGBDM, FE_DGNedelec, FE_DGRaviartThomas
 * </ul> 
 *
 * @note The implementation of vector valued DG elements is supported
 * by the class FE_DGVector, in the way, that only the vector
 * polynomial space has to be provided. The actual class derived from
 * this only has to implement a constructor and
 * FiniteElement::get_name().
 *
 *  @ingroup feall
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
