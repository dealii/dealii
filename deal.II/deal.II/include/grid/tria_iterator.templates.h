/*----------------------------   tria_iterator.templates.h     ---------------------------*/
/*      Copyright W. Bangerth, University of Heidelberg, 1998                             */
/*      $Id$                 */
#ifndef __tria_iterator_templates_H
#define __tria_iterator_templates_H
/*----------------------------   tria_iterator.templates.h     ---------------------------*/

#include <grid/tria.h>
#include <grid/tria_iterator.h>


/* Note: This file only contains template definitions and will thus
   not produce an object file. It is rather thought to be included
   into the *_accessor.cc files.
*/



/*------------------------ Functions: TriaRawIterator ------------------*/


template <int dim, typename Accessor>
inline
TriaRawIterator<dim,Accessor>::TriaRawIterator () :
		accessor (0, -2, -2, 0) {};



template <int dim, typename Accessor>
inline
TriaRawIterator<dim,Accessor>::TriaRawIterator (const TriaRawIterator<dim,Accessor> &i) :
		accessor (i.accessor) {};



template <int dim, typename Accessor>
inline
TriaRawIterator<dim,Accessor>::TriaRawIterator (Triangulation<dim> *parent,
						const int           level,
						const int           index,
						const typename Accessor::AccessorData *local_data) :
		accessor (parent, level, index, local_data) {};



template <int dim, typename Accessor>
inline
TriaRawIterator<dim,Accessor> &
TriaRawIterator<dim,Accessor>::operator = (const TriaRawIterator<dim,Accessor> &i) {
  accessor.copy_from (i.accessor);
  
  return *this;
};



template <int dim, typename Accessor>
inline
bool
TriaRawIterator<dim,Accessor>::operator == (const TriaRawIterator<dim,Accessor> &i) const {
  return accessor == i.accessor;
};



template <int dim, typename Accessor>
inline
bool
TriaRawIterator<dim,Accessor>::operator != (const TriaRawIterator<dim,Accessor> &i) const {
  return accessor != i.accessor;
};



template <int dim, typename Accessor>
inline
TriaRawIterator<dim,Accessor>
TriaRawIterator<dim,Accessor>::operator ++ (int) {
  TriaRawIterator<dim,Accessor> tmp(*this);
  operator++ ();
  
  return tmp;
};



template <int dim, typename Accessor>
inline
TriaRawIterator<dim,Accessor>
TriaRawIterator<dim,Accessor>::operator -- (int) {
  TriaRawIterator<dim,Accessor> tmp(*this);
  operator-- ();
  
  return tmp;
};




/*-----------------------  functions: TriaIterator ---------------*/


template <int dim, typename Accessor>
inline
TriaIterator<dim,Accessor>::TriaIterator () :
		TriaRawIterator<dim,Accessor> () {};


template <int dim, typename Accessor>
inline
TriaIterator<dim,Accessor>::
TriaIterator (const TriaIterator<dim,Accessor> &i) :
		TriaRawIterator<dim,Accessor> (static_cast<TriaRawIterator<dim,Accessor> >(i)) {};



template <int dim, typename Accessor>
inline
TriaIterator<dim,Accessor>::
TriaIterator (const TriaRawIterator<dim,Accessor> &i) :
		TriaRawIterator<dim,Accessor> (i)
{
#ifdef DEBUG
				   // do this like this, because:
				   // if we write
				   // "Assert (past_the_end || used)"
				   // used() is called anyway, even if
				   // state==past_the_end, and will then
				   // throw the exception!
  if (state() != past_the_end)
    Assert (accessor.used(),
	    ExcAssignmentOfUnusedObject());
#endif  
};



template <int dim, typename Accessor>
inline
TriaIterator<dim,Accessor>::TriaIterator (Triangulation<dim> *parent,
					  const int           level,
					  const int           index,
					  const typename Accessor::AccessorData *local_data) :
		TriaRawIterator<dim,Accessor> (parent, level, index, local_data)
{
#ifdef DEBUG
				   // do this like this, because:
				   // if we write
				   // "Assert (past_the_end || used)"
				   // used() is called anyway, even if
				   // state==past_the_end, and will then
				   // throw the exception!
  if (state() != past_the_end)
    Assert (accessor.used(),
	    ExcAssignmentOfUnusedObject());
#endif  
};



template <int dim, typename Accessor>
inline
TriaIterator<dim,Accessor> &
TriaIterator<dim,Accessor>::operator = (const TriaIterator<dim,Accessor> &i) {
  accessor.copy_from (i.accessor);
  return *this;
};


template <int dim, typename Accessor>
inline
TriaIterator<dim,Accessor> &
TriaIterator<dim,Accessor>::operator = (const TriaRawIterator<dim,Accessor> &i) {
  accessor.copy_from (i.accessor);
#ifdef DEBUG
				   // do this like this, because:
				   // if we write
				   // "Assert (past_the_end || used)"
				   // used() is called anyway, even if
				   // state==past_the_end, and will then
				   // throw the exception!
  if (state() != past_the_end) 
    Assert (accessor.used(),
	    ExcAssignmentOfUnusedObject());
#endif  
  return *this;
};



template <int dim, typename Accessor>
inline
TriaIterator<dim,Accessor> & TriaIterator<dim,Accessor>::operator ++ () {
  while (TriaRawIterator<dim,Accessor>::operator++(),
	 (state() == valid))
    if (accessor.used() == true)
      return *this;
  return *this;
};



template <int dim, typename Accessor>
inline
TriaIterator<dim,Accessor>  TriaIterator<dim,Accessor>::operator ++ (int) {
  TriaIterator<dim,Accessor> tmp(*this);
  operator++ ();
  
  return tmp;
};



template <int dim, typename Accessor>
inline
TriaIterator<dim,Accessor> &
TriaIterator<dim,Accessor>::operator -- () {
  while (TriaRawIterator<dim,Accessor>::operator--(),
	 (state() == valid))
    if (accessor.used() == true)
      return *this;
  return *this;
};


template <int dim, typename Accessor>
inline
TriaIterator<dim,Accessor>
TriaIterator<dim,Accessor>::operator -- (int) {
  TriaIterator<dim,Accessor> tmp(*this);
  operator-- ();
  
  return tmp;
};



/*-----------------------  functions: TriaActiveIterator ---------------*/


template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor>::TriaActiveIterator () :
		TriaIterator<dim,Accessor> () {};



template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor>::
TriaActiveIterator (const TriaActiveIterator<dim,Accessor> &i) :
		TriaIterator<dim,Accessor> (static_cast<TriaIterator<dim,Accessor> >(i)) {};



template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor>::
TriaActiveIterator (const TriaRawIterator<dim,Accessor> &i) :
		TriaIterator<dim,Accessor> (i)
{
#ifdef DEBUG
				   // do this like this, because:
				   // if we write
				   // "Assert (past_the_end || used)"
				   // has_children() is called anyway, even if
				   // state==past_the_end, and will then
				   // throw the exception!
  if (state() != past_the_end) 
    Assert (accessor.has_children()==false,
	    ExcAssignmentOfInactiveObject());
#endif  
};


template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor>::
TriaActiveIterator (const TriaIterator<dim,Accessor> &i) :
		TriaIterator<dim,Accessor> (i)
{
#ifdef DEBUG
				   // do this like this, because:
				   // if we write
				   // "Assert (past_the_end || used)"
				   // has_children() is called anyway, even if
				   // state==past_the_end, and will then
				   // throw the exception!
  if (state() != past_the_end) 
    Assert (accessor.has_children()==false,
	    ExcAssignmentOfInactiveObject());
#endif  
};



template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor>::TriaActiveIterator (Triangulation<dim> *parent,
						      const int           level,
						      const int           index,
						      const typename Accessor::AccessorData *local_data) :
		TriaIterator<dim,Accessor> (parent, level, index, local_data)
{
#ifdef DEBUG
				   // do this like this, because:
				   // if we write
				   // "Assert (past_the_end || used)"
				   // has_children() is called anyway, even if
				   // state==past_the_end, and will then
				   // throw the exception!
  if (state() != past_the_end) 
    Assert (accessor.has_children()==false,
	    ExcAssignmentOfInactiveObject());
#endif  
};



template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor> &
TriaActiveIterator<dim,Accessor>::operator = (const TriaActiveIterator<dim,Accessor> &i) {
  accessor.copy_from (i.accessor);
  return *this;
};



template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor> &
TriaActiveIterator<dim,Accessor>::operator = (const TriaRawIterator<dim,Accessor> &i) {
  accessor.copy_from (i.accessor);
#ifdef DEBUG
				   // do this like this, because:
				   // if we write
				   // "Assert (past_the_end || used)"
				   // has_chidlren() is called anyway, even if
				   // state==past_the_end, and will then
				   // throw the exception!
  if (state() != past_the_end) 
    Assert (accessor.used() && accessor.has_children()==false,
	    ExcAssignmentOfInactiveObject());
#endif  
  return *this;
};



template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor> &
TriaActiveIterator<dim,Accessor>::operator = (const TriaIterator<dim,Accessor> &i) {
  accessor.copy_from (i.accessor);
#ifdef DEBUG
				   // do this like this, because:
				   // if we write
				   // "Assert (past_the_end || used)"
				   // has_children() is called anyway, even if
				   // state==past_the_end, and will then
				   // throw the exception!
  if (state() != past_the_end) 
    Assert (accessor.has_children()==false,
	    ExcAssignmentOfInactiveObject());
#endif  
  return *this;
};



template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor> &
TriaActiveIterator<dim,Accessor>::operator ++ () {
  while (TriaIterator<dim,Accessor>::operator++(),
	 (state() == valid))
    if (accessor.has_children() == false)
      return *this;
  return *this;
};



template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor>
TriaActiveIterator<dim,Accessor>::operator ++ (int) {
  TriaActiveIterator<dim,Accessor> tmp(*this);
  operator++ ();
  
  return tmp;
};



template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor> &
TriaActiveIterator<dim,Accessor>::operator -- () {
  while (TriaIterator<dim,Accessor>::operator--(),
	 (state() == valid))
    if (accessor.has_children() == false)
      return *this;
  return *this;
};



template <int dim, typename Accessor>
inline
TriaActiveIterator<dim,Accessor> TriaActiveIterator<dim,Accessor>::operator -- (int) {
  TriaActiveIterator<dim,Accessor> tmp(*this);
  operator-- ();
  
  return tmp;
};




/*----------------------------   tria_iterator.templates.h     ---------------------------*/
/* end of #ifndef __tria_iterator_templates_H */
#endif
/*----------------------------   tria_iterator.templates.h     ---------------------------*/
