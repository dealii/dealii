// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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

#ifndef dealii__quadrature_h
#define dealii__quadrature_h


#include <deal.II/base/config.h>
#include <deal.II/base/point.h>
#include <deal.II/base/subscriptor.h>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup Quadrature */
/*@{*/

/**
 * Base class for quadrature formulae in arbitrary dimensions. This class
 * stores quadrature points and weights on the unit line [0,1], unit square
 * [0,1]x[0,1], etc.
 *
 * There are a number of derived classes, denoting concrete integration
 * formulae. Their names names prefixed by <tt>Q</tt>. Refer to the list of
 * derived classes for more details.
 *
 * The schemes for higher dimensions are typically tensor products of the one-
 * dimensional formulae, but refer to the section on implementation detail
 * below.
 *
 * In order to allow for dimension independent programming, a quadrature
 * formula of dimension zero exists. Since an integral over zero dimensions is
 * the evaluation at a single point, any constructor of such a formula
 * initializes to a single quadrature point with weight one. Access to the
 * weight is possible, while access to the quadrature point is not permitted,
 * since a Point of dimension zero contains no information. The main purpose
 * of these formulae is their use in QProjector, which will create a useful
 * formula of dimension one out of them.
 *
 * <h3>Mathematical background</h3>
 *
 * For each quadrature formula we denote by <tt>m</tt>, the maximal degree of
 * polynomials integrated exactly. This number is given in the documentation
 * of each formula. The order of the integration error is <tt>m+1</tt>, that
 * is, the error is the size of the cell to the <tt>m+1</tt> by the Bramble-
 * Hilbert Lemma. The number <tt>m</tt> is to be found in the documentation of
 * each concrete formula. For the optimal formulae QGauss we have $m = 2N-1$,
 * where N is the constructor parameter to QGauss. The tensor product formulae
 * are exact on tensor product polynomials of degree <tt>m</tt> in each space
 * direction, but they are still only of <tt>m+1</tt>st order.
 *
 * <h3>Implementation details</h3>
 *
 * Most integration formulae in more than one space dimension are tensor
 * products of quadrature formulae in one space dimension, or more generally
 * the tensor product of a formula in <tt>(dim-1)</tt> dimensions and one in
 * one dimension. There is a special constructor to generate a quadrature
 * formula from two others.  For example, the QGauss@<dim@> formulae include
 * <i>N<sup>dim</sup></i> quadrature points in <tt>dim</tt> dimensions, where
 * N is the constructor parameter of QGauss.
 *
 * @note Instantiations for this template are provided for dimensions 0, 1, 2,
 * and 3 (see the section on
 * @ref Instantiations).
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999, 2000, 2005, 2009
 */
template <int dim>
class Quadrature : public Subscriptor
{
public:
  /**
   * Define a typedef for a quadrature that acts on an object of one dimension
   * less. For cells, this would then be a face quadrature.
   */
  typedef Quadrature<dim-1> SubQuadrature;

  /**
   * Constructor.
   *
   * This constructor is marked as explicit to avoid involuntary accidents
   * like in <code>hp::QCollection@<dim@> q_collection(3)</code> where
   * <code>hp::QCollection@<dim@> q_collection(QGauss@<dim@>(3))</code> was
   * meant.
   */
  explicit Quadrature (const unsigned int n_quadrature_points = 0);

  /**
   * Build this quadrature formula as the tensor product of a formula in a
   * dimension one less than the present and a formula in one dimension.
   *
   * <tt>SubQuadrature<dim>::type</tt> expands to <tt>Quadrature<dim-1></tt>.
   */
  Quadrature (const SubQuadrature &,
              const Quadrature<1> &);

  /**
   * Build this quadrature formula as the <tt>dim</tt>-fold tensor product of
   * a formula in one dimension.
   *
   * Assuming that the points in the one-dimensional rule are in ascending
   * order, the points of the resulting rule are ordered lexicographically
   * with <i>x</i> running fastest.
   *
   * In order to avoid a conflict with the copy constructor in 1d, we let the
   * argument be a 0d quadrature formula for dim==1, and a 1d quadrature
   * formula for all other space dimensions.
   */
  explicit Quadrature (const Quadrature<dim != 1 ? 1 : 0> &quadrature_1d);

  /**
   * Copy constructor.
   */
  Quadrature (const Quadrature<dim> &q);

#ifdef DEAL_II_WITH_CXX11
  /**
   * Move constructor. Construct a new quadrature object by transferring the
   * internal data of another quadrature object.
   *
   * @note this constructor is only available if deal.II is configured with
   * C++11 support.
   */
  Quadrature (Quadrature<dim> &&) = default;
#endif

  /**
   * Construct a quadrature formula from given vectors of quadrature points
   * (which should really be in the unit cell) and the corresponding weights.
   * You will want to have the weights sum up to one, but this is not checked.
   */
  Quadrature (const std::vector<Point<dim> > &points,
              const std::vector<double>      &weights);

  /**
   * Construct a dummy quadrature formula from a list of points, with weights
   * set to infinity. The resulting object is therefore not meant to actually
   * perform integrations, but rather to be used with FEValues objects in
   * order to find the position of some points (the quadrature points in this
   * object) on the transformed cell in real space.
   */
  Quadrature (const std::vector<Point<dim> > &points);

  /**
   * Constructor for a one-point quadrature. Sets the weight of this point to
   * one.
   */
  Quadrature (const Point<dim> &point);

  /**
   * Virtual destructor.
   */
  virtual ~Quadrature ();

  /**
   * Assignment operator. Copies contents of #weights and #quadrature_points
   * as well as size.
   */
  Quadrature &operator = (const Quadrature<dim> &);

  /**
   * Test for equality of two quadratures.
   */
  bool operator == (const Quadrature<dim> &p) const;

  /**
   * Set the quadrature points and weights to the values provided in the
   * arguments.
   */
  void initialize(const std::vector<Point<dim> > &points,
                  const std::vector<double>      &weights);

  /**
   * Number of quadrature points.
   */
  unsigned int size () const;

  /**
   * Return the <tt>i</tt>th quadrature point.
   */
  const Point<dim> &point (const unsigned int i) const;

  /**
   * Return a reference to the whole array of quadrature points.
   */
  const std::vector<Point<dim> > &get_points () const;

  /**
   * Return the weight of the <tt>i</tt>th quadrature point.
   */
  double weight (const unsigned int i) const;

  /**
   * Return a reference to the whole array of weights.
   */
  const std::vector<double> &get_weights () const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t memory_consumption () const;

  /**
   * Write or read the data of this object to or from a stream for the purpose
   * of serialization.
   */
  template <class Archive>
  void serialize (Archive &ar, const unsigned int version);

protected:
  /**
   * List of quadrature points. To be filled by the constructors of derived
   * classes.
   */
  std::vector<Point<dim> > quadrature_points;

  /**
   * List of weights of the quadrature points.  To be filled by the
   * constructors of derived classes.
   */
  std::vector<double>      weights;
};


/**
 * Quadrature formula implementing anisotropic distributions of quadrature
 * points on the reference cell. To this end, the tensor product of
 * <tt>dim</tt> one-dimensional quadrature formulas is generated.
 *
 * @note Each constructor can only be used in the dimension matching the
 * number of arguments.
 *
 * @author Guido Kanschat, 2005
 */
template <int dim>
class QAnisotropic : public Quadrature<dim>
{
public:
  /**
   * Constructor for a one-dimensional formula. This one just copies the given
   * quadrature rule.
   */
  QAnisotropic(const Quadrature<1> &qx);

  /**
   * Constructor for a two-dimensional formula.
   */
  QAnisotropic(const Quadrature<1> &qx,
               const Quadrature<1> &qy);

  /**
   * Constructor for a three-dimensional formula.
   */
  QAnisotropic(const Quadrature<1> &qx,
               const Quadrature<1> &qy,
               const Quadrature<1> &qz);
};


/**
 * Quadrature formula constructed by iteration of another quadrature formula
 * in each direction. In more than one space dimension, the resulting
 * quadrature formula is constructed in the usual way by building the tensor
 * product of the respective iterated quadrature formula in one space
 * dimension.
 *
 * In one space dimension, the given base formula is copied and scaled onto a
 * given number of subintervals of length <tt>1/n_copies</tt>. If the
 * quadrature formula uses both end points of the unit interval, then in the
 * interior of the iterated quadrature formula there would be quadrature
 * points which are used twice; we merge them into one with a weight which is
 * the sum of the weights of the left- and the rightmost quadrature point.
 *
 * Since all dimensions higher than one are built up by tensor products of one
 * dimensional and <tt>dim-1</tt> dimensional quadrature formulae, the
 * argument given to the constructor needs to be a quadrature formula in one
 * space dimension, rather than in <tt>dim</tt> dimensions.
 *
 * The aim of this class is to provide a low order formula, where the error
 * constant can be tuned by increasing the number of quadrature points. This
 * is useful in integrating non-differentiable functions on cells.
 *
 * @author Wolfgang Bangerth 1999
 */
template <int dim>
class QIterated : public Quadrature<dim>
{
public:
  /**
   * Constructor. Iterate the given quadrature formula <tt>n_copies</tt> times
   * in each direction.
   */
  QIterated (const Quadrature<1> &base_quadrature,
             const unsigned int   n_copies);

  /**
   * Exception
   */
  DeclExceptionMsg (ExcInvalidQuadratureFormula,
                    "The quadrature formula you provided cannot be used "
                    "as the basis for iteration.");
private:
  /**
   * Check whether the given quadrature formula has quadrature points at the
   * left and right end points of the interval.
   */
  static bool
  uses_both_endpoints (const Quadrature<1> &base_quadrature);
};



/*@}*/

#ifndef DOXYGEN

// -------------------  inline and template functions ----------------


template<int dim>
inline
unsigned int
Quadrature<dim>::size () const
{
  return weights.size();
}


template <int dim>
inline
const Point<dim> &
Quadrature<dim>::point (const unsigned int i) const
{
  AssertIndexRange(i, size());
  return quadrature_points[i];
}



template <int dim>
double
Quadrature<dim>::weight (const unsigned int i) const
{
  AssertIndexRange(i, size());
  return weights[i];
}



template <int dim>
inline
const std::vector<Point<dim> > &
Quadrature<dim>::get_points () const
{
  return quadrature_points;
}



template <int dim>
inline
const std::vector<double> &
Quadrature<dim>::get_weights () const
{
  return weights;
}



template <int dim>
template <class Archive>
inline
void
Quadrature<dim>::serialize (Archive &ar, const unsigned int)
{
  // forward to serialization
  // function in the base class.
  ar   &static_cast<Subscriptor &>(*this);

  ar &quadrature_points &weights;
}



/* -------------- declaration of explicit specializations ------------- */

template <>
Quadrature<0>::Quadrature (const unsigned int);
template <>
Quadrature<0>::Quadrature (const Quadrature<-1> &,
                           const Quadrature<1> &);
template <>
Quadrature<0>::Quadrature (const Quadrature<1> &);
template <>
Quadrature<0>::~Quadrature ();

template <>
Quadrature<1>::Quadrature (const Quadrature<0> &,
                           const Quadrature<1> &);

template <>
Quadrature<1>::Quadrature (const Quadrature<0> &);

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
