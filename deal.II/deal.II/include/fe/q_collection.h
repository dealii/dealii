//----------------------------  q_collection.h  ----------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  q_collection.h  ----------------------------
#ifndef __deal2__q_collection_h
#define __deal2__q_collection_h

#include <base/config.h>
#include <base/subscriptor.h>
#include <base/quadrature.h>
#include <base/smartpointer.h>
#include <fe/fe.h>

#include <vector>
#include <boost/shared_ptr.hpp>


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
                                        * Default constructor. Initialises
                                        * this QCollection.
                                        */
      QCollection ();

                                       /**
                                        * Conversion constructor. This
                                        * constructor creates a QCollection
                                        * from a single quadrature rule. In
                                        * the newly created QCollection, this
                                        * quadrature is used for all active_fe
                                        * indices.
                                        */
      QCollection (const Quadrature<dim> &quadrature);

                                       /**
                                        * Returns the number of quadrature
                                        * pointers stored in this object.
                                        */
      unsigned int size () const;

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
                                        * Adds a new quadrature rule to the
                                        * QCollection.  The quadrature rules
                                        * have to be added in the same order
                                        * as for the FECollection for which
                                        * this quadrature rule collection is
                                        * meant. Thus the reference to the
                                        * quadrature rule for active_fe_index
                                        * 0 has to be added first, followed by
                                        * the quadrature rule for
                                        * active_fe_index 1.
                                        *
                                        * This class creates a copy of the
                                        * given quadrature object, i.e. you
                                        * can do things like
                                        * <tt>add_quadrature(QGauss<dim>(3));</tt>. The
                                        * internal copy is later destroyed by
                                        * this object upon destruction of the
                                        * entire collection.
                                        */
      unsigned int add_quadrature (const Quadrature<dim> &new_quadrature);
    
                                       /**
                                        * Determine an estimate for the
                                        * memory consumption (in bytes)
                                        * of this object.
                                        */
      unsigned int memory_consumption () const;

                                       /**
                                        * Exception
                                        */
      DeclException0 (ExcNoQuadrature);
    
    private:
                                       /**
                                        * Upon construction of a
                                        * <tt>QCollection</tt> the later
                                        * functionality of the class is
                                        * specified.  Either it is a real
                                        * collection, which provides different
                                        * quadrature rules for each
                                        * active_fe_index or its a "unreal"
                                        * collection, which returns a the same
                                        * quadrature rule for all
                                        * active_fe_indices.  This boolean
                                        * remembers which type this object is.
                                        */
      const bool single_quadrature;

                                       /**
                                        * The real container, which stores
                                        * pointers to the different quadrature
                                        * objects.
                                        */
      std::vector<boost::shared_ptr<const Quadrature<dim> > > quadratures;
  };



/* --------------- inline functions ------------------- */

  template <int dim>
  inline
  unsigned int
  QCollection<dim>::size () const 
  {
    return quadratures.size();
  }
  
} // namespace hp


#endif
