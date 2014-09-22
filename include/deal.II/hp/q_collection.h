// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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

#ifndef __deal2__q_collection_h
#define __deal2__q_collection_h

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/fe/fe.h>

#include <vector>
#include <deal.II/base/std_cxx11/shared_ptr.h>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  /**
   * This class implements a collection of quadrature objects in the same way as
   * the hp::FECollection implements a collection of finite element classes.
   *
   * It implements the concepts stated in the @ref hpcollection module described
   * in the doxygen documentation.
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
     * Default constructor. Leads
     * to an empty collection that
     * can later be filled using
     * push_back().
     */
    QCollection ();

    /**
     * Conversion constructor. This
     * constructor creates a
     * QCollection from a single
     * quadrature rule. More
     * quadrature formulas can be
     * added with push_back(), if
     * desired, though it would
     * probably be clearer to add
     * all mappings the same way.
     */
    explicit QCollection (const Quadrature<dim> &quadrature);

    /**
     * Copy constructor.
     */
    QCollection (const QCollection<dim> &q_collection);

    /**
     * Adds a new quadrature rule
     * to the QCollection.  The
     * quadrature rules have to be
     * added in the same order as
     * for the FECollection for
     * which this quadrature rule
     * collection is meant. Thus
     * the reference to the
     * quadrature rule for
     * active_fe_index 0 has to be
     * added first, followed by the
     * quadrature rule for
     * active_fe_index 1.
     *
     * This class creates a copy of
     * the given quadrature object,
     * i.e. you can do things like
     * <tt>push_back(QGauss<dim>(3));</tt>. The
     * internal copy is later
     * destroyed by this object
     * upon destruction of the
     * entire collection.
     */
    void push_back (const Quadrature<dim> &new_quadrature);

    /**
     * Returns a reference to the
     * quadrature rule specified by the
     * argument.
     *
     * @pre @p index must be between zero
     * and the number of elements of the
     * collection.
     */
    const Quadrature<dim> &
    operator[] (const unsigned int index) const;

    /**
     * Returns the number of
     * quadrature pointers stored
     * in this object.
     */
    unsigned int size () const;

    /**
     * Return the maximum number of
     * quadrature points over all
     * the elements of the
     * collection. This is mostly
     * useful to initialize arrays
     * to allocate the maxmimum
     * amount of memory that may be
     * used when re-sizing later on
     * to a articular quadrature
     * formula from within this
     * collection.
     */
    unsigned int max_n_quadrature_points () const;

    /**
     * Determine an estimate for the
     * memory consumption (in bytes)
     * of this object.
     */
    std::size_t memory_consumption () const;

    /**
     * Exception
     */
    DeclException0 (ExcNoQuadrature);

  private:
    /**
     * The real container, which stores
     * pointers to the different quadrature
     * objects.
     */
    std::vector<std_cxx11::shared_ptr<const Quadrature<dim> > > quadratures;
  };



  /* --------------- inline functions ------------------- */

  template <int dim>
  inline
  unsigned int
  QCollection<dim>::size () const
  {
    return quadratures.size();
  }



  template <int dim>
  inline
  unsigned int
  QCollection<dim>::max_n_quadrature_points () const
  {
    Assert (quadratures.size() > 0,
            ExcMessage ("You can't call this function for an empty collection"));

    unsigned int m = 0;
    for (unsigned int i=0; i<quadratures.size(); ++i)
      if (quadratures[i]->size() > m)
        m = quadratures[i]->size();

    return m;
  }



  template <int dim>
  inline
  const Quadrature<dim> &
  QCollection<dim>::operator[] (const unsigned int index) const
  {
    Assert (index < quadratures.size (),
            ExcIndexRange (index, 0, quadratures.size ()));
    return *quadratures[index];
  }



  template <int dim>
  inline
  QCollection<dim>::QCollection ()
  {}



  template <int dim>
  inline
  QCollection<dim>::QCollection (const Quadrature<dim> &quadrature)
  {
    quadratures
    .push_back (std_cxx11::shared_ptr<const Quadrature<dim> >(new Quadrature<dim>(quadrature)));
  }



  template <int dim>
  inline
  QCollection<dim>::
  QCollection (const QCollection<dim> &q_collection)
    :
    Subscriptor (),
    // copy the array
    // of shared
    // pointers. nothing
    // bad should
    // happen -- they
    // simply all point
    // to the same
    // objects, and the
    // last one to die
    // will delete the
    // mappings
    quadratures (q_collection.quadratures)
  {}



  template <int dim>
  inline
  std::size_t
  QCollection<dim>::memory_consumption () const
  {
    return (sizeof(*this) +
            MemoryConsumption::memory_consumption (quadratures));
  }


  template <int dim>
  inline
  void
  QCollection<dim>::push_back (const Quadrature<dim> &new_quadrature)
  {
    quadratures
    .push_back (std_cxx11::shared_ptr<const Quadrature<dim> >(new Quadrature<dim>(new_quadrature)));
  }

} // namespace hp


DEAL_II_NAMESPACE_CLOSE

#endif
