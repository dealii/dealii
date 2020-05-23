// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2020 by the deal.II authors
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

#ifndef dealii_q_collection_h
#define dealii_q_collection_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/fe/fe.h>

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
   * module described in the doxygen documentation.
   *
   * @ingroup hp hpcollection
   *
   * @author Oliver Kayser-Herold, 2005
   */
  template <int dim>
  class QCollection : public Subscriptor
  {
  public:
    /**
     * Default constructor. Leads to an empty collection that can later be
     * filled using push_back().
     */
    QCollection() = default;

    /**
     * Conversion constructor. This constructor creates a QCollection from a
     * single quadrature rule. More quadrature formulas can be added with
     * push_back(), if desired, though it would probably be clearer to add all
     * mappings the same way.
     */
    explicit QCollection(const Quadrature<dim> &quadrature);

    /**
     * Constructor. This constructor creates a QCollection from one or
     * more quadrature objects passed to the constructor. For this
     * call to be valid, all arguments need to be of types derived
     * from class Quadrature<dim>.
     */
    template <class... QTypes>
    explicit QCollection(const QTypes &... quadrature_objects);

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
    void
    push_back(const Quadrature<dim> &new_quadrature);

    /**
     * Return a reference to the quadrature rule specified by the argument.
     *
     * @pre @p index must be between zero and the number of elements of the
     * collection.
     */
    const Quadrature<dim> &operator[](const unsigned int index) const;

    /**
     * Equality comparison operator. All stored Quadrature objects are compared
     * in order.
     */
    bool
    operator==(const QCollection<dim> &q_collection) const;

    /**
     * Return the number of quadrature pointers stored in this object.
     */
    unsigned int
    size() const;

    /**
     * Return the maximum number of quadrature points over all the elements of
     * the collection. This is mostly useful to initialize arrays to allocate
     * the maximum amount of memory that may be used when re-sizing later on
     * to a articular quadrature formula from within this collection.
     */
    unsigned int
    max_n_quadrature_points() const;

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const;

    /**
     * Exception
     */
    DeclException0(ExcNoQuadrature);

  private:
    /**
     * The real container, which stores pointers to the different quadrature
     * objects.
     */
    std::vector<std::shared_ptr<const Quadrature<dim>>> quadratures;
  };



  /* --------------- inline functions ------------------- */

  template <int dim>
  template <class... QTypes>
  QCollection<dim>::QCollection(const QTypes &... quadrature_objects)
  {
    static_assert(is_base_of_all<Quadrature<dim>, QTypes...>::value,
                  "Not all of the input arguments of this function "
                  "are derived from Quadrature<dim>!");

    // loop over all of the given arguments and add the quadrature objects to
    // this collection. Inlining the definition of q_pointers causes internal
    // compiler errors on GCC 7.1.1 so we define it separately:
    const auto q_pointers = {
      (static_cast<const Quadrature<dim> *>(&quadrature_objects))...};
    for (const auto p : q_pointers)
      push_back(*p);
  }



  template <int dim>
  inline unsigned int
  QCollection<dim>::size() const
  {
    return quadratures.size();
  }



  template <int dim>
  inline unsigned int
  QCollection<dim>::max_n_quadrature_points() const
  {
    Assert(quadratures.size() > 0,
           ExcMessage("You can't call this function for an empty collection"));

    unsigned int m = 0;
    for (unsigned int i = 0; i < quadratures.size(); ++i)
      if (quadratures[i]->size() > m)
        m = quadratures[i]->size();

    return m;
  }



  template <int dim>
  inline const Quadrature<dim> &QCollection<dim>::
                                operator[](const unsigned int index) const
  {
    AssertIndexRange(index, quadratures.size());
    return *quadratures[index];
  }



  template <int dim>
  inline bool
  QCollection<dim>::operator==(const QCollection<dim> &q_collection) const
  {
    const unsigned int n_quadratures = size();
    if (n_quadratures != q_collection.size())
      return false;

    for (unsigned int i = 0; i < n_quadratures; ++i)
      if (!(*quadratures[i] == q_collection[i]))
        return false;

    return true;
  }



  template <int dim>
  inline QCollection<dim>::QCollection(const Quadrature<dim> &quadrature)
  {
    quadratures.push_back(std::make_shared<const Quadrature<dim>>(quadrature));
  }



  template <int dim>
  inline std::size_t
  QCollection<dim>::memory_consumption() const
  {
    return (sizeof(*this) + MemoryConsumption::memory_consumption(quadratures));
  }


  template <int dim>
  inline void
  QCollection<dim>::push_back(const Quadrature<dim> &new_quadrature)
  {
    quadratures.push_back(
      std::make_shared<const Quadrature<dim>>(new_quadrature));
  }

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif
