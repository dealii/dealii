//----------------------------  fe_collection.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_collection.h  ---------------------------
#ifndef __deal2__fe_collection_h
#define __deal2__fe_collection_h

#include <base/config.h>
#include <base/subscriptor.h>
#include <fe/fe.h>


/**
 * This class acts as a collection of finite element objects used in the
 * hp::DoFHandler(). It is thus to a hp::DoFHandler() what a
 * FiniteElement is to a DoFHandler. This collection stores copies
 * of the original elements added to it, and is therefore handling its memory
 * manegement itself.
 *
 * In addition to offering access to the elements of the collection, this
 * class provides access to the maximal number of degrees of freedom per
 * vertex, line, etc, to allow allocation of as much memory as is necessary in
 * the worst case when using the finite elements associated with the cells of
 * a triangulation.
 * 
 * @author Wolfgang Bangerth, 2003
 */
template <int dim>
class FECollection : public Subscriptor
{
  public:
                                     /**
                                      * Destructor. Destroy the elements of
                                      * the collection.
                                      */
    ~FECollection ();
    
                                     /**
                                      * Add a finite element. This function
                                      * generates a copy of the given element,
                                      * i.e. you can do things like
                                      * <tt>add_fe(FE_Q<dim>(1));</tt>. The
                                      * internal copy is later destroyed by
                                      * this object upon destruction of the
                                      * entire collection.
                                      *
                                      * The returned value is the new number
                                      * of finite elements in the collection.
                                      *
                                      * When a new element is added, it needs
                                      * to have the same number of vector
                                      * components as all other elements
                                      * already in the collection.
                                      */
    unsigned int
    add_fe (const FiniteElement<dim> &new_fe);

                                     /**
                                      * Get a reference to the given element
                                      * in this collection.
                                      */
    const FiniteElement<dim> &
    get_fe (const unsigned int index) const;

                                     /**
                                      * Return the number of finite element
                                      * objects stored in this collection.
                                      */
    unsigned int n_finite_elements () const;

                                     /**
                                      * Return the number of vector components
                                      * of the finite elements in this
                                      * collection.  This number must be the
                                      * same for all elements in the
                                      * collection.
                                      */
    unsigned int n_components () const;

                                     /**
                                      * Return the maximal number of degrees
                                      * of freedom per vertex over all
                                      * elements of this collection.
                                      */
    unsigned int max_dofs_per_vertex () const;

                                     /**
                                      * Return the maximal number of degrees
                                      * of freedom per line over all elements
                                      * of this collection.
                                      */
    unsigned int max_dofs_per_line () const;

                                     /**
                                      * Return the maximal number of degrees
                                      * of freedom per quad over all elements
                                      * of this collection.
                                      */
    unsigned int max_dofs_per_quad () const;

                                     /**
                                      * Return the maximal number of degrees
                                      * of freedom per hex over all elements
                                      * of this collection.
                                      */
    unsigned int max_dofs_per_hex () const;

                                     /**
                                      * Return the maximal number of degrees
                                      * of freedom per face over all elements
                                      * of this collection.
                                      */
    unsigned int max_dofs_per_face () const;

                                     /**
                                      * Return the maximal number of degrees
                                      * of freedom per cell over all elements
                                      * of this collection.
                                      */
    unsigned int max_dofs_per_cell () const;

                                     /**
                                      * Return an estimate for the memory
                                      * allocated for this object.
                                      */
    unsigned int memory_consumption () const;

                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcNoFiniteElements);    
    
  private:
                                     /**
                                      * Array of pointers to the finite
                                      * elements stored by this collection.
                                      */
    std::vector<SmartPointer<const FiniteElement<dim> > > finite_elements;
};



/* --------------- inline functions ------------------- */

template <int dim>
inline
unsigned int
FECollection<dim>::n_finite_elements () const 
{
  return finite_elements.size();
}


template <int dim>
inline
unsigned int
FECollection<dim>::n_components () const
{
  Assert (finite_elements.size () > 0, ExcNoFiniteElements());
  return finite_elements[0]->n_components ();
}


template <int dim>
inline
const FiniteElement<dim> &
FECollection<dim>::get_fe (const unsigned int index) const 
{
  Assert (index < finite_elements.size(),
          ExcIndexRange (index, 0, finite_elements.size()));
  return *finite_elements[index];
}



template <int dim>
unsigned int
FECollection<dim>::max_dofs_per_vertex () const 
{
  Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
  unsigned int max = 0;
  for (unsigned int i=0; i<finite_elements.size(); ++i)
    if (finite_elements[i]->dofs_per_vertex > max)
      max = finite_elements[i]->dofs_per_vertex;

  return max;
}



template <int dim>
unsigned int
FECollection<dim>::max_dofs_per_line () const 
{
  Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
  unsigned int max = 0;
  for (unsigned int i=0; i<finite_elements.size(); ++i)
    if (finite_elements[i]->dofs_per_line > max)
      max = finite_elements[i]->dofs_per_line;

  return max;
}



template <int dim>
unsigned int
FECollection<dim>::max_dofs_per_quad () const 
{
  Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
  unsigned int max = 0;
  for (unsigned int i=0; i<finite_elements.size(); ++i)
    if (finite_elements[i]->dofs_per_quad > max)
      max = finite_elements[i]->dofs_per_quad;

  return max;
}



template <int dim>
unsigned int
FECollection<dim>::max_dofs_per_hex () const 
{
  Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
  unsigned int max = 0;
  for (unsigned int i=0; i<finite_elements.size(); ++i)
    if (finite_elements[i]->dofs_per_hex > max)
      max = finite_elements[i]->dofs_per_hex;

  return max;
}



template <int dim>
unsigned int
FECollection<dim>::max_dofs_per_face () const 
{
  Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
  unsigned int max = 0;
  for (unsigned int i=0; i<finite_elements.size(); ++i)
    if (finite_elements[i]->dofs_per_face > max)
      max = finite_elements[i]->dofs_per_face;

  return max;
}



template <int dim>
unsigned int
FECollection<dim>::max_dofs_per_cell () const 
{
  Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
  unsigned int max = 0;
  for (unsigned int i=0; i<finite_elements.size(); ++i)
    if (finite_elements[i]->dofs_per_cell > max)
      max = finite_elements[i]->dofs_per_cell;

  return max;
}


#endif
