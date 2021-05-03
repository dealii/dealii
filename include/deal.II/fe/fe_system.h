// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_fe_system_h
#  define dealii_fe_system_h


/*----------------------------   fe_system.h     ---------------------------*/


#  include <deal.II/base/config.h>

#  include <deal.II/base/thread_management.h>

#  include <deal.II/fe/fe.h>
#  include <deal.II/fe/fe_tools.h>

#  include <memory>
#  include <type_traits>
#  include <utility>
#  include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#  ifndef DOXYGEN
template <int dim, int spacedim>
class FE_Enriched;
#  endif

/**
 * This class provides an interface to group several elements together into
 * one, vector-valued element. As example, consider the Taylor-Hood element
 * that is used for the solution of the Stokes and Navier-Stokes equations:
 * There, the velocity (of which there are as many components as the dimension
 * $d$ of the domain) is discretized with $Q_2$ elements and the pressure with
 * $Q_1$ elements. Mathematically, the finite element space for the coupled
 * problem is then often written as $V_h = Q_2^d \times Q_1$ where the
 * exponentiation is understood to be the tensor product of spaces -- i.e.,
 * in 2d, we have $V_h=Q_2\times Q_2\times Q_1$ -- and tensor products
 * lead to vectors where each component of the vector-valued function
 * space corresponds to a scalar function in one of the $Q_2$ or $Q_1$
 * spaces. Using the FESystem class, this space is created using
 * @code
 *   FESystem<dim> taylor_hood_fe (FE_Q<dim>(2)^dim,   // velocity components
 *                                 FE_Q<dim>(1));      // pressure component
 * @endcode
 * The creation of this element here corresponds to taking tensor-product
 * powers of the $Q_2$ element in the first line of the list of arguments
 * to the FESystem constructor, and then concatenation via another tensor
 * product with the element in the second line. This kind of construction
 * is used, for example, in the step-22 tutorial program.
 *
 * Similarly, step-8 solves an elasticity equation where we need to solve
 * for the displacement of a solid object. The displacement again has
 * $d$ components if the domain is $d$-dimensional, and so the combined
 * finite element is created using
 * @code
 *   FESystem<dim> displacement_fe (FE_Q<dim>(1)^dim);
 * @endcode
 * where now each (vector) component of the combined element corresponds to
 * a $Q_1$ space.
 *
 * To the outside world, FESystem objects look just like a usual
 * finite element object, they just happen to be composed of several other
 * finite elements that are possibly of different type. These "base elements"
 * can themselves have multiple components and, in particular, could
 * also be vector-valued -- for example, if one of the base elements
 * is an FESystem itself (see also below). An example is given in the
 * documentation of namespace FETools::Compositing, when using the
 * "tensor product" strategy.
 *
 * %Vector valued elements are discussed in a number of
 * tutorial programs, for example step-8, step-20, step-21, step-22, and in
 * particular in the
 * @ref vector_valued
 * module.
 *
 * @dealiiVideoLecture{19,20}
 *
 *
 * <h3>FESystem, components and blocks</h3>
 *
 * An FESystem, except in the most trivial case, produces a vector-valued
 * finite element with several components. The number of components
 * n_components() corresponds to the dimension of the solution function in the
 * PDE system, and correspondingly also to the number of equations your PDE
 * system has. For example, the mixed Laplace system covered in step-20 has
 * $d+1$ components in $d$ space dimensions: the scalar pressure and the $d$
 * components of the velocity vector. Similarly, the elasticity equation
 * covered in step-8 has $d$ components in $d$ space dimensions. In general,
 * the number of components of a FESystem element is the accumulated number of
 * components of all base elements times their multiplicities. A bit more on
 * components is also given in the
 * @ref GlossComponent "glossary entry on components".
 *
 * While the concept of components is important from the viewpoint of a
 * partial differential equation, the finite element side looks a bit
 * different Since not only FESystem, but also vector-valued elements like
 * FE_RaviartThomas, have several components. The concept needed here is a
 * @ref GlossBlock "block".
 * Each block encompasses the set of degrees of freedom associated with a
 * single base element of an FESystem, where base elements with multiplicities
 * count multiple times. These blocks are usually addressed using the
 * information in DoFHandler::block_info(). The number of blocks of a FESystem
 * object is simply the sum of all multiplicities of base elements and is
 * given by n_blocks().
 *
 * For example, the FESystem for the Taylor-Hood element for the
 * three-dimensional Stokes problem can be built using the code
 * @code
 * const FE_Q<3> u(2);
 * const FE_Q<3> p(1);
 * FESystem<3> sys1(u,3, p,1);
 * @endcode
 * or more concisely via
 * @code
 * FESystem<3> sys1(FE_Q<3>(2),3,
 *                  FE_Q<3>(1),1);
 * @endcode
 * or even shorter (mimicking the mathematical notation that we are dealing
 * with a $Q_2^3 \times Q_1$ element):
 * @code
 * FESystem<3> sys1(FE_Q<3>(2)^3,
 *                  FE_Q<3>(1));
 * @endcode
 *
 * This example creates an FESystem @p sys1 with four components, three for
 * the velocity components and one for the pressure, and also four blocks with
 * the degrees of freedom of each of the velocity components and the pressure
 * in a separate block each. The number of blocks is four since the first base
 * element is repeated three times.
 *
 * On the other hand, a Taylor-Hood element can also be constructed using
 *
 * @code
 * FESystem<3> U(u,3);
 * FESystem<3> sys2(U, p);
 * @endcode
 *
 * The FESystem @p sys2 created here has the same four components, but the
 * degrees of freedom are distributed into only two blocks. The first block
 * has all velocity degrees of freedom from @p U, while the second block
 * contains the pressure degrees of freedom. Note that while @p U itself has 3
 * blocks, the FESystem @p sys2 does not attempt to split @p U into its base
 * elements but considers it a block of its own. By blocking all velocities
 * into one system first as in @p sys2, we achieve the same block structure
 * that would be generated if instead of using a $Q_2^3$ element for the
 * velocities we had used vector-valued base elements, for instance like using
 * a mixed discretization of Darcy's law using
 *
 * @code
 * FE_RaviartThomas<3> u(1);
 * FE_DGQ<3> p(1);
 * FESystem<3> sys3(u, p);
 * @endcode
 *
 * This example also produces a system with four components, but only two
 * blocks.
 *
 * In most cases, the composed element behaves as if it were a usual element.
 * It just has more degrees of freedom than most of the "common" elements.
 * However the underlying structure is visible in the restriction,
 * prolongation and interface constraint matrices, which do not couple the
 * degrees of freedom of the base elements. E.g. the continuity requirement is
 * imposed for the shape functions of the subobjects separately; no
 * requirement exist between shape functions of different subobjects, i.e. in
 * the above example: on a hanging node, the respective value of the @p u
 * velocity is only coupled to @p u at the vertices and the line on the larger
 * cell next to this vertex, but there is no interaction with @p v and @p w of
 * this or the other cell.
 *
 *
 * <h3>Internal information on numbering of degrees of freedom</h3>
 *
 * The overall numbering of degrees of freedom is as follows: for each
 * subobject (vertex, line, quad, or hex), the degrees of freedom are numbered
 * such that we run over all subelements first, before turning for the next
 * dof on this subobject or for the next subobject. For example, for an
 * element of three components in one space dimension, the first two
 * components being cubic lagrange elements and the third being a quadratic
 * lagrange element, the ordering for the system <tt>s=(u,v,p)</tt> is:
 *
 * <ul>
 * <li> First vertex: <tt>u0, v0, p0 = s0, s1, s2</tt>
 * <li> Second vertex: <tt>u1, v1, p1 = s3, s4, s5</tt>
 * <li> First component on the line: <tt>u2, u3 = s4, s5</tt>
 * <li> Second component on the line: <tt>v2, v3 = s6, s7</tt>.
 * <li> Third component on the line: <tt>p2 = s8</tt>.
 * </ul>
 * That said, you should not rely on this numbering in your application as
 * these %internals might change in future. Rather use the functions
 * system_to_component_index() and component_to_system_index().
 *
 * For more information on the template parameter <tt>spacedim</tt> see the
 * documentation of Triangulation.
 *
 * @ingroup febase fe vector_valued
 *
 */
template <int dim, int spacedim = dim>
class FESystem : public FiniteElement<dim, spacedim>
{
public:
  /**
   * Delete default constructor so that `FESystem(FEPairs &&... fe_pairs)` is
   * not accidentally picked if no FiniteElement is provided.
   */
  FESystem() = delete;

  /**
   * Constructor. Take a finite element and the number of elements you want to
   * group together using this class.
   *
   * The object @p fe is not actually used for anything other than creating a
   * copy that will then be owned by the current object. In other words, it is
   * completely fine to call this constructor with a temporary object for the
   * finite element, as in this code snippet:
   * @code
   *   FESystem<dim> fe (FE_Q<dim>(2), 2);
   * @endcode
   * Here, <code>FE_Q@<dim@>(2)</code> constructs an unnamed, temporary object
   * that is passed to the FESystem constructor to create a finite element
   * that consists of two components, both of which are quadratic FE_Q
   * elements. The temporary is destroyed again at the end of the code that
   * corresponds to this line, but this does not matter because FESystem
   * creates its own copy of the FE_Q object.
   *
   * This constructor (or its variants below) is used in essentially all
   * tutorial programs that deal with vector valued problems. See step-8,
   * step-20, step-22 and others for use cases. Also see the module on
   * @ref vector_valued "Handling vector valued problems".
   *
   * @dealiiVideoLecture{19,20}
   *
   * @param[in] fe The finite element that will be used to represent the
   * components of this composed element.
   * @param[in] n_elements An integer denoting how many copies of @p fe this
   * element should consist of.
   */
  FESystem(const FiniteElement<dim, spacedim> &fe,
           const unsigned int                  n_elements);

  /**
   * Constructor for mixed discretizations with two base elements.
   *
   * See the other constructor above for an explanation of the general idea of
   * composing elements.
   */
  FESystem(const FiniteElement<dim, spacedim> &fe1,
           const unsigned int                  n1,
           const FiniteElement<dim, spacedim> &fe2,
           const unsigned int                  n2);

  /**
   * Constructor for mixed discretizations with three base elements.
   *
   * See the other constructor above for an explanation of the general idea of
   * composing elements.
   */
  FESystem(const FiniteElement<dim, spacedim> &fe1,
           const unsigned int                  n1,
           const FiniteElement<dim, spacedim> &fe2,
           const unsigned int                  n2,
           const FiniteElement<dim, spacedim> &fe3,
           const unsigned int                  n3);

  /**
   * Constructor for mixed discretizations with four base elements.
   *
   * See the first of the other constructors above for an explanation of the
   * general idea of composing elements.
   */
  FESystem(const FiniteElement<dim, spacedim> &fe1,
           const unsigned int                  n1,
           const FiniteElement<dim, spacedim> &fe2,
           const unsigned int                  n2,
           const FiniteElement<dim, spacedim> &fe3,
           const unsigned int                  n3,
           const FiniteElement<dim, spacedim> &fe4,
           const unsigned int                  n4);

  /**
   * Constructor for mixed discretizations with five base elements.
   *
   * See the first of the other constructors above for an explanation of the
   * general idea of composing elements.
   */
  FESystem(const FiniteElement<dim, spacedim> &fe1,
           const unsigned int                  n1,
           const FiniteElement<dim, spacedim> &fe2,
           const unsigned int                  n2,
           const FiniteElement<dim, spacedim> &fe3,
           const unsigned int                  n3,
           const FiniteElement<dim, spacedim> &fe4,
           const unsigned int                  n4,
           const FiniteElement<dim, spacedim> &fe5,
           const unsigned int                  n5);

  /**
   * Same as above but for any number of base elements. Pointers to the base
   * elements and their multiplicities are passed as vectors to this
   * constructor. The lengths of these vectors are assumed to be equal.
   *
   * As above, the finite element objects pointed to by the first argument are
   * not actually used other than to create copies internally. Consequently,
   * you can delete these pointers immediately again after calling this
   * constructor.
   *
   * <h4>How to use this constructor</h4>
   *
   * Using this constructor is a bit awkward at times because you need to pass
   * two vectors in a place where it may not be straightforward to construct
   * such a vector -- for example, in the member initializer list of a class
   * with an FESystem member variable. For example, if your main class looks
   * like this:
   * @code
   *   template <int dim>
   *   class MySimulator {
   *   public:
   *     MySimulator (const unsigned int polynomial_degree);
   *   private:
   *     FESystem<dim> fe;
   *   };
   *
   *   template <int dim>
   *   MySimulator<dim>::MySimulator (const unsigned int polynomial_degree)
   *     :
   *     fe (...)  // what to pass here???
   *   {}
   * @endcode
   *
   * Using the C++11 language standard (or later) you could do something like
   * this to create an element with four base elements and multiplicities 1,
   * 2, 3 and 4:
   * @code
   *   template <int dim>
   *   MySimulator<dim>::MySimulator (const unsigned int polynomial_degree)
   *     :
   *     fe (std::vector<const FiniteElement<dim>*> { new FE_Q<dim>(1),
   *                                                  new FE_Q<dim>(2),
   *                                                  new FE_Q<dim>(3),
   *                                                  new FE_Q<dim>(4) },
   *         std::vector<unsigned int> { 1, 2, 3, 4 })
   *   {}
   * @endcode
   * This creates two vectors in place and initializes them using the
   * initializer list enclosed in braces <code>{ ... }</code>.
   *
   * This code has a problem: it creates four memory leaks because the first
   * vector above is created with pointers to elements that are allocated with
   * <code>new</code> but never destroyed.
   *
   * The solution to the second of these problems is to create two static
   * member functions that can create vectors. Here is an example:
   * @code
   *   template <int dim>
   *   class MySimulator {
   *   public:
   *     MySimulator (const unsigned int polynomial_degree);
   *
   *   private:
   *     FESystem<dim> fe;
   *
   *     static std::vector<const FiniteElement<dim>*>
   *     create_fe_list (const unsigned int polynomial_degree);
   *
   *     static std::vector<unsigned int>
   *     create_fe_multiplicities ();
   *   };
   *
   *   template <int dim>
   *   std::vector<const FiniteElement<dim>*>
   *   MySimulator<dim>::create_fe_list (const unsigned int polynomial_degree)
   *   {
   *     std::vector<const FiniteElement<dim>*> fe_list;
   *     fe_list.push_back (new FE_Q<dim>(1));
   *     fe_list.push_back (new FE_Q<dim>(2));
   *     fe_list.push_back (new FE_Q<dim>(3));
   *     fe_list.push_back (new FE_Q<dim>(4));
   *     return fe_list;
   *   }
   *
   *   template <int dim>
   *   std::vector<unsigned int>
   *   MySimulator<dim>::create_fe_multiplicities ()
   *   {
   *     std::vector<unsigned int> multiplicities;
   *     multiplicities.push_back (1);
   *     multiplicities.push_back (2);
   *     multiplicities.push_back (3);
   *     multiplicities.push_back (4);
   *     return multiplicities;
   *   }
   *
   *   template <int dim>
   *   MySimulator<dim>::MySimulator (const unsigned int polynomial_degree)
   *     :
   *     fe (create_fe_list (polynomial_degree),
   *         create_fe_multiplicities ())
   *   {}
   * @endcode
   *
   * The way this works is that we have two static member functions that
   * create the necessary vectors to pass to the constructor of the member
   * variable <code>fe</code>. They need to be static because they are called
   * during the constructor of <code>MySimulator</code> at a time when the
   * <code>*this</code> object isn't fully constructed and, consequently,
   * regular member functions cannot be called yet.
   *
   * The code above does not solve the problem with the memory leak yet,
   * though: the <code>create_fe_list()</code> function creates a vector of
   * pointers, but nothing destroys these. This is the solution:
   * @code
   * template <int dim>
   * class MySimulator
   * {
   * public:
   *   MySimulator (const unsigned int polynomial_degree);
   *
   * private:
   *   FESystem<dim> fe;
   *
   *   struct VectorElementDestroyer
   *   {
   *     const std::vector<const FiniteElement<dim>*> data;
   *
   *     VectorElementDestroyer(
   *       const std::vector<const FiniteElement<dim>*> &pointers);
   *
   *      // destructor to delete the pointers
   *     ~VectorElementDestroyer ();
   *
   *     const std::vector<const FiniteElement<dim>*> & get_data () const;
   *   };
   *
   *   static std::vector<const FiniteElement<dim>*>
   *   create_fe_list (const unsigned int polynomial_degree);
   *
   *   static std::vector<unsigned int>
   *   create_fe_multiplicities ();
   * };
   *
   * template <int dim>
   * MySimulator<dim>::VectorElementDestroyer::
   * VectorElementDestroyer(
   *   const std::vector<const FiniteElement<dim>*> &pointers)
   *   :
   *   data(pointers)
   * {}
   *
   * template <int dim>
   * MySimulator<dim>::VectorElementDestroyer::
   * ~VectorElementDestroyer ()
   * {
   *   for (unsigned int i=0; i<data.size(); ++i)
   *     delete data[i];
   * }
   *
   * template <int dim>
   * const std::vector<const FiniteElement<dim>*> &
   * MySimulator<dim>::VectorElementDestroyer::
   * get_data () const
   * {
   *   return data;
   * }
   *
   * template <int dim>
   * MySimulator<dim>::MySimulator (const unsigned int polynomial_degree)
   * :
   * fe (VectorElementDestroyer(create_fe_list (polynomial_degree)).get_data(),
   *     create_fe_multiplicities ())
   * {}
   * @endcode
   *
   * In other words, the vector we receive from the
   * <code>create_fe_list()</code> is packed into a temporary object of type
   * <code>VectorElementDestroyer</code>; we then get the vector from this
   * temporary object immediately to pass it to the constructor of
   * <code>fe</code>; and finally, the <code>VectorElementDestroyer</code>
   * destructor is called at the end of the entire expression (after the
   * constructor of <code>fe</code> has finished) and destroys the elements of
   * the temporary vector. Voila: not short nor elegant, but it works!
   */
  FESystem(const std::vector<const FiniteElement<dim, spacedim> *> &fes,
           const std::vector<unsigned int> &multiplicities);

#  if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1900
  /**
   * Constructor taking an arbitrary number of parameters of type
   * <code>std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned
   * int></code>. In combination with FiniteElement::operator^, this allows to
   * construct FESystem objects as follows:
   * @code
   *   FiniteElementType1<dim,spacedim> fe_1;
   *   FiniteElementType1<dim,spacedim> fe_2;
   *   FESystem<dim,spacedim> fe_system ( fe_1^dim, fe_2 );
   * @endcode
   *
   * The `fe_1` and `fe_2` objects are not actually used for anything other than
   * creating a copy that will then be owned by the current object. In other
   * words, it is completely fine to call this constructor with a temporary
   * object for the finite element, as in this code snippet:
   * @code
   *   FESystem<dim> fe (FE_Q<dim>(2)^2);
   * @endcode
   * Here, <code>FE_Q@<dim@>(2)</code> constructs an unnamed, temporary object
   * that is passed to the FESystem constructor to create a finite element
   * that consists of two components, both of which are quadratic FE_Q
   * elements. The temporary is destroyed again at the end of the code that
   * corresponds to this line, but this does not matter because FESystem
   * creates its own copy of the FE_Q object.
   *
   * As a shortcut, this constructor also allows calling
   * @code
   *   FESystem<dim> fe (FE_Q<dim>(2)^dim, FE_Q<dim>(1));
   * @endcode
   * instead of the more explicit
   * @code
   *   FESystem<dim> fe (FE_Q<dim>(2)^dim, FE_Q<dim>(1)^1);
   * @endcode
   * In other words, if no multiplicity for an element is explicitly specified
   * via the exponentiation operation, then it is assumed to be one (as one
   * would have expected).
   *
   * @warning This feature is not available for Intel compilers
   * prior to version 19.0. Defining this
   * constructor leads to internal compiler errors for Intel compilers prior
   * to 18.0.
   */
  template <
    class... FEPairs,
    typename = typename enable_if_all<
      (std::is_same<typename std::decay<FEPairs>::type,
                    std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                              unsigned int>>::value ||
       std::is_base_of<FiniteElement<dim, spacedim>,
                       typename std::decay<FEPairs>::type>::value)...>::type>
  FESystem(FEPairs &&... fe_pairs);

  /**
   * Same as above allowing the following syntax:
   * @code
   *   FiniteElementType1<dim,spacedim> fe_1;
   *   FiniteElementType1<dim,spacedim> fe_2;
   *   FESystem<dim,spacedim> fe_system = { fe_1^dim, fe_2^1 };
   * @endcode
   *
   * @warning This feature is not available for Intel compilers
   * prior to version 19.0. The constructor is just not selected for overload
   * resolution.
   */
  FESystem(
    const std::initializer_list<
      std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
      &fe_systems);
#  endif

  /**
   * Copy constructor. This constructor is deleted, i.e., copying
   * FESystem objects is not allowed.
   */
  FESystem(const FESystem<dim, spacedim> &) = delete;

  /**
   * Move constructor.
   */
  FESystem(FESystem<dim, spacedim> &&other_fe_system) noexcept
    : FiniteElement<dim, spacedim>(std::move(other_fe_system))
  {
    base_elements = std::move(other_fe_system.base_elements);
    generalized_support_points_index_table =
      std::move(other_fe_system.generalized_support_points_index_table);
  }

  /**
   * Destructor.
   */
  virtual ~FESystem() override = default;

  /**
   * Return a string that uniquely identifies a finite element. This element
   * returns a string that is composed of the strings @p name1...@p nameN
   * returned by the basis elements. From these, we create a sequence
   * <tt>FESystem<dim>[name1^m1-name2^m2-...-nameN^mN]</tt>, where @p mi are
   * the multiplicities of the basis elements. If a multiplicity is equal to
   * one, then the superscript is omitted.
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  // make variant with ComponentMask also available:
  using FiniteElement<dim, spacedim>::get_sub_fe;

  /**
   * @copydoc FiniteElement<dim,spacedim>::get_sub_fe()
   */
  virtual const FiniteElement<dim, spacedim> &
  get_sub_fe(const unsigned int first_component,
             const unsigned int n_selected_components) const override;

  /**
   * Return the value of the @p ith shape function at the point @p p.  @p p is
   * a point on the reference element. Since this finite element is always
   * vector-valued, we return the value of the only non-zero component of the
   * vector value of this shape function. If the shape function has more than
   * one non-zero component (which we refer to with the term non-primitive),
   * then throw an exception of type @p ExcShapeFunctionNotPrimitive.
   *
   * An @p ExcUnitShapeValuesDoNotExist is thrown if the shape values of the
   * @p FiniteElement (corresponding to the @p ith shape function) depend on
   * the shape of the cell in real space.
   */
  virtual double
  shape_value(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Return the value of the @p componentth vector component of the @p ith
   * shape function at the point @p p. See the FiniteElement base class for
   * more information about the semantics of this function.
   *
   * Since this element is vector valued in general, it relays the computation
   * of these values to the base elements.
   */
  virtual double
  shape_value_component(const unsigned int i,
                        const Point<dim> & p,
                        const unsigned int component) const override;

  /**
   * Return the gradient of the @p ith shape function at the point @p p. @p p
   * is a point on the reference element, and likewise the gradient is the
   * gradient on the unit cell with respect to unit cell coordinates. Since
   * this finite element is always vector-valued, we return the value of the
   * only non-zero component of the vector value of this shape function. If
   * the shape function has more than one non-zero component (which we refer
   * to with the term non-primitive), then throw an exception of type @p
   * ExcShapeFunctionNotPrimitive.
   *
   * An @p ExcUnitShapeValuesDoNotExist is thrown if the shape values of the
   * @p FiniteElement (corresponding to the @p ith shape function) depend on
   * the shape of the cell in real space.
   */
  virtual Tensor<1, dim>
  shape_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Return the gradient of the @p componentth vector component of the @p ith
   * shape function at the point @p p. See the FiniteElement base class for
   * more information about the semantics of this function.
   *
   * Since this element is vector valued in general, it relays the computation
   * of these values to the base elements.
   */
  virtual Tensor<1, dim>
  shape_grad_component(const unsigned int i,
                       const Point<dim> & p,
                       const unsigned int component) const override;

  /**
   * Return the tensor of second derivatives of the @p ith shape function at
   * point @p p on the unit cell. The derivatives are derivatives on the unit
   * cell with respect to unit cell coordinates. Since this finite element is
   * always vector-valued, we return the value of the only non-zero component
   * of the vector value of this shape function. If the shape function has
   * more than one non-zero component (which we refer to with the term non-
   * primitive), then throw an exception of type @p
   * ExcShapeFunctionNotPrimitive.
   *
   * An @p ExcUnitShapeValuesDoNotExist is thrown if the shape values of the
   * @p FiniteElement (corresponding to the @p ith shape function) depend on
   * the shape of the cell in real space.
   */
  virtual Tensor<2, dim>
  shape_grad_grad(const unsigned int i, const Point<dim> &p) const override;

  /**
   * Return the second derivatives of the @p componentth vector component of
   * the @p ith shape function at the point @p p. See the FiniteElement base
   * class for more information about the semantics of this function.
   *
   * Since this element is vector valued in general, it relays the computation
   * of these values to the base elements.
   */
  virtual Tensor<2, dim>
  shape_grad_grad_component(const unsigned int i,
                            const Point<dim> & p,
                            const unsigned int component) const override;

  /**
   * Return the tensor of third derivatives of the @p ith shape function at
   * point @p p on the unit cell. The derivatives are derivatives on the unit
   * cell with respect to unit cell coordinates. Since this finite element is
   * always vector-valued, we return the value of the only non-zero component
   * of the vector value of this shape function. If the shape function has
   * more than one non-zero component (which we refer to with the term non-
   * primitive), then throw an exception of type @p
   * ExcShapeFunctionNotPrimitive.
   *
   * An @p ExcUnitShapeValuesDoNotExist is thrown if the shape values of the
   * @p FiniteElement (corresponding to the @p ith shape function) depend on
   * the shape of the cell in real space.
   */
  virtual Tensor<3, dim>
  shape_3rd_derivative(const unsigned int i,
                       const Point<dim> & p) const override;

  /**
   * Return the third derivatives of the @p componentth vector component of
   * the @p ith shape function at the point @p p. See the FiniteElement base
   * class for more information about the semantics of this function.
   *
   * Since this element is vector valued in general, it relays the computation
   * of these values to the base elements.
   */
  virtual Tensor<3, dim>
  shape_3rd_derivative_component(const unsigned int i,
                                 const Point<dim> & p,
                                 const unsigned int component) const override;

  /**
   * Return the tensor of fourth derivatives of the @p ith shape function at
   * point @p p on the unit cell. The derivatives are derivatives on the unit
   * cell with respect to unit cell coordinates. Since this finite element is
   * always vector-valued, we return the value of the only non-zero component
   * of the vector value of this shape function. If the shape function has
   * more than one non-zero component (which we refer to with the term non-
   * primitive), then throw an exception of type @p
   * ExcShapeFunctionNotPrimitive.
   *
   * An @p ExcUnitShapeValuesDoNotExist is thrown if the shape values of the
   * @p FiniteElement (corresponding to the @p ith shape function) depend on
   * the shape of the cell in real space.
   */
  virtual Tensor<4, dim>
  shape_4th_derivative(const unsigned int i,
                       const Point<dim> & p) const override;

  /**
   * Return the fourth derivatives of the @p componentth vector component of
   * the @p ith shape function at the point @p p. See the FiniteElement base
   * class for more information about the semantics of this function.
   *
   * Since this element is vector valued in general, it relays the computation
   * of these values to the base elements.
   */
  virtual Tensor<4, dim>
  shape_4th_derivative_component(const unsigned int i,
                                 const Point<dim> & p,
                                 const unsigned int component) const override;

  /**
   * Return the matrix interpolating from the given finite element to the
   * present one. The size of the matrix is then @p dofs_per_cell times
   * <tt>source.n_dofs_per_cell()</tt>.
   *
   * These matrices are available if source and destination element are both
   * @p FESystem elements, have the same number of base elements with same
   * element multiplicity, and if these base elements also implement their @p
   * get_interpolation_matrix functions. Otherwise, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &matrix) const override;

  /**
   * Access to a composing element. The index needs to be smaller than the
   * number of base elements. Note that the number of base elements may in
   * turn be smaller than the number of components of the system element, if
   * the multiplicities are greater than one.
   */
  virtual const FiniteElement<dim, spacedim> &
  base_element(const unsigned int index) const override;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * Projection from a fine grid space onto a coarse grid space. Overrides the
   * respective method in FiniteElement, implementing lazy evaluation
   * (initialize when requested).
   *
   * If this projection operator is associated with a matrix @p P, then the
   * restriction of this matrix @p P_i to a single child cell is returned
   * here.
   *
   * The matrix @p P is the concatenation or the sum of the cell matrices @p
   * P_i, depending on the #restriction_is_additive_flags. This distinguishes
   * interpolation (concatenation) and projection with respect to scalar
   * products (summation).
   *
   * Row and column indices are related to coarse grid and fine grid spaces,
   * respectively, consistent with the definition of the associated operator.
   *
   * If projection matrices are not implemented in the derived finite element
   * class, this function aborts with an exception of type
   * FiniteElement::ExcProjectionVoid. You can check whether this would happen
   * by first calling the restriction_is_implemented() or the
   * isotropic_restriction_is_implemented() function.
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * Embedding matrix between grids. Overrides the respective method in
   * FiniteElement, implementing lazy evaluation (initialize when queried).
   *
   * The identity operator from a coarse grid space into a fine grid space is
   * associated with a matrix @p P. The restriction of this matrix @p P_i to a
   * single child cell is returned here.
   *
   * The matrix @p P is the concatenation, not the sum of the cell matrices @p
   * P_i. That is, if the same non-zero entry <tt>j,k</tt> exists in two
   * different child matrices @p P_i, the value should be the same in both
   * matrices and it is copied into the matrix @p P only once.
   *
   * Row and column indices are related to fine grid and coarse grid spaces,
   * respectively, consistent with the definition of the associated operator.
   *
   * These matrices are used by routines assembling the prolongation matrix
   * for multi-level methods.  Upon assembling the transfer matrix between
   * cells using this matrix array, zero elements in the prolongation matrix
   * are discarded and will not fill up the transfer matrix.
   *
   * If prolongation matrices are not implemented in one of the base finite
   * element classes, this function aborts with an exception of type
   * FiniteElement::ExcEmbeddingVoid. You can check whether this would happen
   * by first calling the prolongation_is_implemented() or the
   * isotropic_prolongation_is_implemented() function.
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * Given an index in the natural ordering of indices on a face, return the
   * index of the same degree of freedom on the cell.
   *
   * To explain the concept, consider the case where we would like to know
   * whether a degree of freedom on a face, for example as part of an FESystem
   * element, is primitive. Unfortunately, the is_primitive() function in the
   * FiniteElement class takes a cell index, so we would need to find the cell
   * index of the shape function that corresponds to the present face index.
   * This function does that.
   *
   * Code implementing this would then look like this:
   * @code
   * for (i=0; i<dofs_per_face; ++i)
   *  if (fe.is_primitive(fe.face_to_cell_index(i, some_face_no)))
   *   ... do whatever
   * @endcode
   * The function takes additional arguments that account for the fact that
   * actual faces can be in their standard ordering with respect to the cell
   * under consideration, or can be flipped, oriented, etc.
   *
   * @param face_dof_index The index of the degree of freedom on a face. This
   * index must be between zero and dofs_per_face.
   * @param face The number of the face this degree of freedom lives on. This
   * number must be between zero and GeometryInfo::faces_per_cell.
   * @param face_orientation One part of the description of the orientation of
   * the face. See
   * @ref GlossFaceOrientation.
   * @param face_flip One part of the description of the orientation of the
   * face. See
   * @ref GlossFaceOrientation.
   * @param face_rotation One part of the description of the orientation of
   * the face. See
   * @ref GlossFaceOrientation.
   * @return The index of this degree of freedom within the set of degrees of
   * freedom on the entire cell. The returned value will be between zero and
   * dofs_per_cell.
   */
  virtual unsigned int
  face_to_cell_index(const unsigned int face_dof_index,
                     const unsigned int face,
                     const bool         face_orientation = true,
                     const bool         face_flip        = false,
                     const bool         face_rotation = false) const override;

  /**
   * Implementation of the respective function in the base class.
   */
  virtual Point<dim>
  unit_support_point(const unsigned int index) const override;

  /**
   * Implementation of the respective function in the base class.
   */
  virtual Point<dim - 1>
  unit_face_support_point(const unsigned int index,
                          const unsigned int face_no = 0) const override;

  /**
   * Return a list of constant modes of the element. The returns table has as
   * many rows as there are components in the element and dofs_per_cell
   * columns. To each component of the finite element, the row in the returned
   * table contains a basis representation of the constant function 1 on the
   * element. Concatenates the constant modes of each base element.
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  /**
   * @name Functions to support hp
   * @{
   */

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp-compatible".
   *
   * This function returns @p true if and only if all its base elements return
   * @p true for this function.
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Base elements of this element will have to implement this function. They
   * may only provide interpolation matrices for certain source finite
   * elements, for example those from the same family. If they don't implement
   * interpolation from a given element, then they must throw an exception of
   * type FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented, which
   * will get propagated out from this element.
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double> &                matrix,
                                const unsigned int face_no = 0) const override;


  /**
   * Return the matrix interpolating from a face of one element to the
   * subface of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Base elements of this element will have to implement this function. They
   * may only provide interpolation matrices for certain source finite
   * elements, for example those from the same family. If they don't implement
   * interpolation from a given element, then they must throw an exception of
   * type FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented, which
   * will get propagated out from this element.
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double> &                matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * If, on a vertex, several finite elements are active, the hp-code first
   * assigns the degrees of freedom of each of these FEs different global
   * indices. It then calls this function to find out which of them should get
   * identical values, and consequently can receive the same global DoF index.
   * This function therefore returns a list of identities between DoFs of the
   * present finite element object with the DoFs of @p fe_other, which is a
   * reference to a finite element object representing one of the other finite
   * elements active on this particular vertex. The function computes which of
   * the degrees of freedom of the two finite element objects are equivalent,
   * both numbered between zero and the corresponding value of
   * n_dofs_per_vertex() of the two finite elements. The first index of each
   * pair denotes one of the vertex dofs of the present element, whereas the
   * second is the corresponding index of the other finite element.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * @copydoc FiniteElement::compare_for_domination()
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  //@}

  /**
   * Implementation of the
   * FiniteElement::convert_generalized_support_point_values_to_dof_values()
   * function.
   *
   * This function simply calls
   * FiniteElement::convert_generalized_support_point_values_to_dof_values
   * of the base elements and re-assembles everything into the output
   * argument. If a base element is non-interpolatory the corresponding dof
   * values are filled with "signaling" NaNs instead.
   *
   * The function fails if none of the base elements of the FESystem are
   * interpolatory.
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double> &              dof_values) const override;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   *
   * This function is made virtual, since finite element objects are usually
   * accessed through pointers to their base class, rather than the class
   * itself.
   */
  virtual std::size_t
  memory_consumption() const override;

protected:
  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim> &       quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<dim, spacedim>::get_face_data;

  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_face_data(
    const UpdateFlags               update_flags,
    const Mapping<dim, spacedim> &  mapping,
    const hp::QCollection<dim - 1> &quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual std::unique_ptr<
    typename FiniteElement<dim, spacedim>::InternalDataBase>
  get_subface_data(
    const UpdateFlags             update_flags,
    const Mapping<dim, spacedim> &mapping,
    const Quadrature<dim - 1> &   quadrature,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim> &                                     quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  using FiniteElement<dim, spacedim>::fill_fe_face_values;

  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1> &                            quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const Quadrature<dim - 1> &                                 quadrature,
    const Mapping<dim, spacedim> &                              mapping,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                       spacedim>
      &                                                            mapping_data,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
    dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                       spacedim>
      &output_data) const override;

  /**
   * Do the work for the three <tt>fill_fe*_values</tt> functions.
   *
   * Calls (among other things) <tt>fill_fe_([sub]face)_values</tt> of the
   * base elements. Calls @p fill_fe_values if
   * <tt>face_no==invalid_face_no</tt> and <tt>sub_no==invalid_face_no</tt>;
   * calls @p fill_fe_face_values if <tt>face_no==invalid_face_no</tt> and
   * <tt>sub_no!=invalid_face_no</tt>; and calls @p fill_fe_subface_values if
   * <tt>face_no!=invalid_face_no</tt> and <tt>sub_no!=invalid_face_no</tt>.
   */
  template <int dim_1>
  void
  compute_fill(
    const Mapping<dim, spacedim> &                              mapping,
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          sub_no,
    const hp::QCollection<dim_1> &                              quadrature,
    const CellSimilarity::Similarity                            cell_similarity,
    const typename Mapping<dim, spacedim>::InternalDataBase &mapping_internal,
    const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_data,
    const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &mapping_data,
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
      &output_data) const;

private:
  /**
   * Value to indicate that a given face or subface number is invalid.
   */
  static const unsigned int invalid_face_number = numbers::invalid_unsigned_int;

  /**
   * Pointers to underlying finite element objects.
   *
   * This object contains a pointer to each contributing element of a mixed
   * discretization and its multiplicity. It is created by the constructor and
   * constant afterwards.
   */
  std::vector<std::pair<std::unique_ptr<const FiniteElement<dim, spacedim>>,
                        unsigned int>>
    base_elements;

  /**
   * An index table that maps generalized support points of a base element
   * to the vector of generalized support points of the FE System.
   * It holds true that
   * @code
   *   auto n = generalized_support_points_index_table[i][j];
   *   generalized_support_points[n] ==
   *           base_elements[i].generalized_support_points[j];
   * @endcode
   * for each base element (indexed by i) and each g. s. point of the base
   * element (index by j).
   */
  std::vector<std::vector<std::size_t>> generalized_support_points_index_table;

  /**
   * This function is simply singled out of the constructors since there are
   * several of them. It sets up the index table for the system as well as @p
   * restriction and @p prolongation matrices.
   */
  void
  initialize(const std::vector<const FiniteElement<dim, spacedim> *> &fes,
             const std::vector<unsigned int> &multiplicities);

  /**
   * Used by @p initialize.
   */
  void
  build_interface_constraints();

  /**
   * A function that computes the hp_vertex_dof_identities(),
   * hp_line_dof_identities(), or hp_quad_dof_identities(), depending on the
   * value of the template parameter.
   */
  template <int structdim>
  std::vector<std::pair<unsigned int, unsigned int>>
  hp_object_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                           const unsigned int face_no = 0) const;

  /**
   * Usually: Fields of cell-independent data.
   *
   * However, here, this class does not itself store the data but only
   * pointers to @p InternalData objects for each of the base elements.
   */
  class InternalData : public FiniteElement<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * Constructor. Is called by the @p get_data function. Sets the size of
     * the @p base_fe_datas vector to @p n_base_elements.
     */
    InternalData(const unsigned int n_base_elements);

    /**
     * Destructor. Deletes all @p InternalDatas whose pointers are stored by
     * the @p base_fe_datas vector.
     */
    ~InternalData() override;

    /**
     * Give write-access to the pointer to a @p InternalData of the @p
     * base_noth base element.
     */
    void
    set_fe_data(
      const unsigned int base_no,
      std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>);

    /**
     * Give read-access to the pointer to a @p InternalData of the @p
     * base_noth base element.
     */
    typename FiniteElement<dim, spacedim>::InternalDataBase &
    get_fe_data(const unsigned int base_no) const;

    /**
     * Give read-access to the pointer to an object to which into which the
     * <code>base_no</code>th base element will write its output when calling
     * FiniteElement::fill_fe_values() and similar functions.
     */
    internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &
    get_fe_output_object(const unsigned int base_no) const;

  private:
    /**
     * Pointers to @p InternalData objects for each of the base elements. They
     * are accessed to by the @p set_ and @p get_fe_data functions.
     *
     * The size of this vector is set to @p n_base_elements by the
     * InternalData constructor.  It is filled by the @p get_data function.
     * Note that since the data for each instance of a base class is
     * necessarily the same, we only need as many of these objects as there
     * are base elements, irrespective of their multiplicity.
     */
    typename std::vector<
      std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>>
      base_fe_datas;

    /**
     * A collection of objects to which the base elements will write their
     * output when we call FiniteElement::fill_fe_values() and related
     * functions on them.
     *
     * The size of this vector is set to @p n_base_elements by the
     * InternalData constructor.
     */
    mutable std::vector<
      internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>>
      base_fe_output_objects;
  };

  /**
   * Mutex for protecting initialization of restriction and embedding matrix.
   */
  mutable std::mutex mutex;

  friend class FE_Enriched<dim, spacedim>;
};

//------------------------variadic template constructor------------------------

#  ifndef DOXYGEN
namespace internal
{
  namespace FESystemImplementation
  {
    template <int dim, int spacedim>
    unsigned int
    count_nonzeros(
      const std::initializer_list<
        std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
        &fe_systems)
    {
      return std::count_if(
        fe_systems.begin(),
        fe_systems.end(),
        [](const std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                           unsigned int> &fe_system) {
          return fe_system.second > 0;
        });
    }



    template <int dim, int spacedim>
    std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>
    promote_to_fe_pair(const FiniteElement<dim, spacedim> &fe)
    {
      return std::make_pair(std::move(fe.clone()), 1u);
    }



    template <int dim, int spacedim>
    auto
    promote_to_fe_pair(std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                                 unsigned int> &&p)
      -> decltype(
        std::forward<std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                               unsigned int>>(p))
    {
      return std::forward<
        std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>(
        p);
    }
  } // namespace FESystemImplementation
} // namespace internal



#    if !defined(__INTEL_COMPILER) || __INTEL_COMPILER >= 1900
// We are just forwarding/delegating to the constructor taking a
// std::initializer_list. If we decide to remove the deprecated constructors, we
// might just use the variadic constructor with a suitable static_assert instead
// of the std::enable_if.
template <int dim, int spacedim>
template <class... FEPairs, typename>
FESystem<dim, spacedim>::FESystem(FEPairs &&... fe_pairs)
  : FESystem<dim, spacedim>(
      {internal::FESystemImplementation::promote_to_fe_pair<dim, spacedim>(
        std::forward<FEPairs>(fe_pairs))...})
{}



template <int dim, int spacedim>
FESystem<dim, spacedim>::FESystem(
  const std::initializer_list<
    std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>, unsigned int>>
    &fe_systems)
  : FiniteElement<dim, spacedim>(
      FETools::Compositing::multiply_dof_numbers<dim, spacedim>(fe_systems),
      FETools::Compositing::compute_restriction_is_additive_flags<dim,
                                                                  spacedim>(
        fe_systems),
      FETools::Compositing::compute_nonzero_components<dim, spacedim>(
        fe_systems))
  , base_elements(internal::FESystemImplementation::count_nonzeros(fe_systems))
{
  std::vector<const FiniteElement<dim, spacedim> *> fes;
  std::vector<unsigned int>                         multiplicities;

  const auto extract =
    [&fes, &multiplicities](
      const std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
                      unsigned int> &fe_system) {
      fes.push_back(fe_system.first.get());
      multiplicities.push_back(fe_system.second);
    };

  for (const auto &p : fe_systems)
    extract(p);

  initialize(fes, multiplicities);
}
#    endif

#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
/*----------------------------  fe_system.h  ---------------------------*/
