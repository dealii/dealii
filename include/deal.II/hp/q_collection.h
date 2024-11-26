// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_q_collection_h
#define dealii_q_collection_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/hp/collection.h>

#include <memory>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  /**
   * This class implements a collection of quadrature objects in the same way
   * as the hp::FECollection implements a collection of finite element
   * classes.
   *
   * It implements the concepts stated in the
   * @ref hpcollection
   * topic described in the doxygen documentation.
   *
   * @ingroup hp hpcollection
   */
  template <int dim>
  class QCollection : public Collection<Quadrature<dim>>
  {
  public:
    /**
     * Default constructor. Leads to an empty collection that can later be
     * filled using push_back().
     */
    QCollection() = default;

    /**
     * Copy constructor.
     */
    template <int dim_in>
    QCollection(const QCollection<dim_in> &other);

    /**
     * Conversion constructor. This constructor creates a QCollection from a
     * single quadrature rule. More quadrature formulas can be added with
     * push_back(), if desired, though it would probably be clearer to add all
     * mappings the same way.
     */
    template <int dim_in>
    explicit QCollection(const Quadrature<dim_in> &quadrature);

    /**
     * Constructor. This constructor creates a QCollection from one or
     * more quadrature objects passed to the constructor. For this
     * call to be valid, all arguments need to be of types derived
     * from class Quadrature<dim>.
     */
    template <class... QTypes>
    explicit QCollection(const QTypes &...quadrature_objects);

    /**
     * Add a new quadrature rule to the QCollection. In most cases, you will
     * want to add quadrature rules in the same order as the elements were
     * added to the hp::FECollection for which this quadrature rule collection
     * is meant. If done this way, the hp::FEValues objects with which you will
     * use both hp::FECollection and hp::QCollection objects will automatically
     * choose corresponding elements and quadrature formulas. On the other hand,
     * it is possible to use arbitrary combinations of elements and quadrature
     * formulas in hp::FECollection and hp::QCollection objects when
     * specifically specifying appropriate indices in calls to
     * hp::FEValues::reinit() or hp::FEFaceValues::reinit(). In those cases,
     * there need not be a correspondence between elements of the
     * hp::FECollection and hp::QCollection objects; they need not even be of
     * the same size in this case.
     *
     * The same arguments about the order of elements of collections can, by
     * the way, also be made about the elements of hp::MappingCollection
     * objects.
     *
     * This class creates a copy of the given quadrature object, i.e., you can
     * do things like <tt>push_back(QGauss<dim>(3));</tt>. The internal copy
     * is later destroyed by this object upon destruction of the entire
     * collection.
     */
    template <int dim_in>
    void
    push_back(const Quadrature<dim_in> &new_quadrature);

    /**
     * Equality comparison operator. All stored Quadrature objects are compared
     * in order.
     */
    bool
    operator==(const QCollection<dim> &q_collection) const;

    /**
     * Return the maximum number of quadrature points over all the elements of
     * the collection. This is mostly useful to initialize arrays to allocate
     * the maximum amount of memory that may be used when re-sizing later on
     * to a articular quadrature formula from within this collection.
     */
    unsigned int
    max_n_quadrature_points() const;

    /**
     * Exception
     */
    DeclException0(ExcNoQuadrature);
  };



  /* --------------- inline functions ------------------- */

  template <int dim>
  template <int dim_in>
  QCollection<dim>::QCollection(const QCollection<dim_in> &other)
  {
    for (unsigned int i = 0; i < other.size(); ++i)
      push_back(other[i]);
  }



  template <int dim>
  template <class... QTypes>
  QCollection<dim>::QCollection(const QTypes &...quadrature_objects)
  {
    // loop over all of the given arguments and add the quadrature objects to
    // this collection. Inlining the definition of q_pointers causes internal
    // compiler errors on GCC 7.1.1 so we define it separately:
    if (is_base_of_all<Quadrature<dim>, QTypes...>::value)
      {
        const auto q_pointers = {
          (reinterpret_cast<const Quadrature<dim> *>(&quadrature_objects))...};
        for (const auto p : q_pointers)
          push_back(*p);
      }
    else if (is_base_of_all<Quadrature<1>, QTypes...>::value)
      {
        const auto q_pointers = {
          (reinterpret_cast<const Quadrature<1> *>(&quadrature_objects))...};
        for (const auto p : q_pointers)
          push_back(*p);
      }
    else
      {
        DEAL_II_NOT_IMPLEMENTED();
      }
  }



  template <int dim>
  inline unsigned int
  QCollection<dim>::max_n_quadrature_points() const
  {
    Assert(this->size() > 0,
           ExcMessage("You can't call this function for an empty collection"));

    unsigned int max = 0;
    for (unsigned int i = 0; i < this->size(); ++i)
      max = std::max(max, this->operator[](i).size());

    return max;
  }



  template <int dim>
  inline bool
  QCollection<dim>::operator==(const QCollection<dim> &q_collection) const
  {
    const unsigned int n_quadratures = this->size();
    if (n_quadratures != q_collection.size())
      return false;

    for (unsigned int i = 0; i < n_quadratures; ++i)
      if ((this->operator[](i) == q_collection[i]) == false)
        return false;

    return true;
  }



  template <int dim>
  template <int dim_in>
  inline QCollection<dim>::QCollection(const Quadrature<dim_in> &quadrature)
  {
    this->push_back(quadrature);
  }


  template <int dim>
  template <int dim_in>
  inline void
  QCollection<dim>::push_back(const Quadrature<dim_in> &new_quadrature)
  {
    Collection<Quadrature<dim>>::push_back(
      std::make_shared<const Quadrature<dim>>(new_quadrature));
  }

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif
