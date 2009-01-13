//----------------------------  fe_collection.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2004, 2006 by the deal.II authors
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
#include <base/std_cxx0x/shared_ptr.h>
#include <fe/fe.h>

DEAL_II_NAMESPACE_OPEN

namespace hp
{
  
/**
 * This class acts as a collection of finite element objects used in the
 * hp::DoFHandler. It is thus to a hp::DoFHandler what a
 * FiniteElement is to a ::DoFHandler.
 *
 * It implements the concepts stated in the @ref hpcollection module described
 * in the doxygen documentation.
 *
 * In addition to offering access to the elements of the collection, this
 * class provides access to the maximal number of degrees of freedom per
 * vertex, line, etc, to allow allocation of as much memory as is necessary in
 * the worst case when using the finite elements associated with the cells of
 * a triangulation.
 *
 * This class has not yet been implemented for the use in the codimension
 * one case (<tt>spacedim != dim </tt>).
 * 
 * @ingroup hp hpcollection
 * 
 * @author Wolfgang Bangerth, 2003
 */
  template <int dim, int spacedim=dim>
  class FECollection : public Subscriptor
  {
    public:
                                       /**
                                        * Default constructor. Leads
                                        * to an empty collection that
                                        * can later be filled using
                                        * push_back().
                                        */
      FECollection ();

                                       /**
                                        * Conversion constructor. This
                                        * constructor creates a
                                        * FECollection from a single
                                        * finite element. More finite
                                        * element objects can be added
                                        * with push_back(), if
                                        * desired, though it would
                                        * probably be clearer to add
                                        * all mappings the same way.
                                        */
      explicit FECollection (const FiniteElement<dim,spacedim> &fe);

                                       /**
                                        * Copy constructor.
                                        */
      FECollection (const FECollection<dim,spacedim> &fe_collection);

                                       /**
                                        * Add a finite element. This
                                        * function generates a copy of
                                        * the given element, i.e. you
                                        * can do things like
                                        * <tt>push_back(FE_Q<dim>(1));</tt>. The
                                        * internal copy is later
                                        * destroyed by this object
                                        * upon destruction of the
                                        * entire collection.
                                        *
                                        * When a new element is added,
                                        * it needs to have the same
                                        * number of vector components
                                        * as all other elements
                                        * already in the collection.
                                        */
      void push_back (const FiniteElement<dim,spacedim> &new_fe);

                                       /**
                                        * Get a reference to the given
                                        * element in this collection.
                                        *
                                        * @pre @p index must be
                                        * between zero and the number
                                        * of elements of the
                                        * collection.
                                        */
      const FiniteElement<dim,spacedim> &
      operator[] (const unsigned int index) const;

                                       /**
                                        * Return the number of finite
                                        * element objects stored in
                                        * this collection.
                                        */
      unsigned int size () const;

                                       /**
                                        * Return the number of vector
                                        * components of the finite elements in
                                        * this collection.  This number must
                                        * be the same for all elements in the
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
					* Return whether all elements
					* in this collection
					* implement the hanging node
					* constraints in the new way,
					* which has to be used to make
					* elements "hp compatible".
					* If this is not the case,
					* the function returns false,
					* which implies, that at least
					* one element in the FECollection
					* does not support the new face
					* interface constraints.
					* On the other hand, if this
					* method does return
					* true, this does not imply
					* that the hp method will work!
					*
					* This behaviour is related to
					* the fact, that FiniteElement
					* classes, which provide the
					* new style hanging node constraints
					* might still not provide
					* them for all possible cases.
					* If FE_Q and FE_RaviartThomas
					* elements are included in the
					* FECollection and both properly implement
					* the get_face_interpolation_matrix
					* method, this method will return
					* true. But the get_face_interpolation_matrix
					* might still fail to find an interpolation
					* matrix between these two elements.
                                        */
      bool hp_constraints_are_implemented () const;


                                       /**
                                        * Exception
                                        */
      DeclException0 (ExcNoFiniteElements);    
    
    private:
                                       /**
                                        * Array of pointers to the finite
                                        * elements stored by this collection.
                                        */
      std::vector<std_cxx0x::shared_ptr<const FiniteElement<dim,spacedim> > > finite_elements;
  };



/* --------------- inline functions ------------------- */

  template <int dim, int spacedim>
  inline
  unsigned int
  FECollection<dim,spacedim>::size () const 
  {
    return finite_elements.size();
  }


  template <int dim, int spacedim>
  inline
  unsigned int
  FECollection<dim,spacedim>::n_components () const
  {
    Assert (finite_elements.size () > 0, ExcNoFiniteElements());
    return finite_elements[0]->n_components ();
  }


  template <int dim, int spacedim>
  inline
  const FiniteElement<dim,spacedim> &
  FECollection<dim,spacedim>::operator[] (const unsigned int index) const 
  {
    Assert (index < finite_elements.size(),
            ExcIndexRange (index, 0, finite_elements.size()));
    return *finite_elements[index];
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_vertex () const 
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_vertex > max)
        max = finite_elements[i]->dofs_per_vertex;

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_line () const 
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_line > max)
        max = finite_elements[i]->dofs_per_line;

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_quad () const 
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_quad > max)
        max = finite_elements[i]->dofs_per_quad;

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_hex () const 
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_hex > max)
        max = finite_elements[i]->dofs_per_hex;

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_face () const 
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_face > max)
        max = finite_elements[i]->dofs_per_face;

    return max;
  }



  template <int dim, int spacedim>
  unsigned int
  FECollection<dim,spacedim>::max_dofs_per_cell () const 
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
    unsigned int max = 0;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      if (finite_elements[i]->dofs_per_cell > max)
        max = finite_elements[i]->dofs_per_cell;

    return max;
  }
  

  template <int dim, int spacedim>
  bool
  FECollection<dim,spacedim>::hp_constraints_are_implemented () const 
  {
    Assert (finite_elements.size() > 0, ExcNoFiniteElements());
  
    bool hp_constraints = true;
    for (unsigned int i=0; i<finite_elements.size(); ++i)
      hp_constraints = hp_constraints &&
		       finite_elements[i]->hp_constraints_are_implemented();
    
    return hp_constraints;
  }
  

} // namespace hp

DEAL_II_NAMESPACE_CLOSE
  
#endif
