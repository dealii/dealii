//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__operator_h
#define __deal2__operator_h

#include <base/config.h>
#include <base/named_data.h>
#include <numerics/event.h>

DEAL_II_NAMESPACE_OPEN

namespace Algorithms
{
/**
 * The abstract base class of all algorithms in this library. An
 * operator is an object with an #operator(), which transforms a set
 * of named vectors into another set of named vectors.
 *
 * Furthermore, an operator can be notified of parameter changes by
 * the calling routine. The outer iteration can notify() the Operator
 * of an Event, which could be for instance a change of mesh, a
 * different time step size or too slow convergence of Newton's
 * method, which would then trigger reassembling of a matrix or
 * similar things.
 *
 * @author Guido Kanschat, 2010
 */
  template <class VECTOR>
  class Operator : public Subscriptor
  {
    public:
				       /**
					* The virtual destructor.
					*/
      ~Operator();

				       /**
					* The actual operation, which
					* is implemented in a derived class.
					*/
      virtual void operator() (NamedData<VECTOR*>& out, const NamedData<VECTOR*>& in) = 0;

				       /**
					* Register an event triggered
					* by an outer iteration.
					*/
      virtual void notify(const Event&);
				       /**
					* Clear all #notifications.
					*/
      void clear_events();
    protected:
				       /**
					* Accumulate reasons for
					* reassembling here. If any of
					* those is set, the function
					* solve() of a terminal
					* application must take care
					* of reassembling the matrix.
					*/
      Event notifications;      
  };

}


DEAL_II_NAMESPACE_CLOSE

#endif
