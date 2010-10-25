//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006, 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/memory_consumption.h>
#include <grid/tria.h>
#include <grid/tria_objects.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>

    

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace Triangulation
  {    
    template <>
    template <int dim, int spacedim>
    typename dealii::Triangulation<dim,spacedim>::raw_line_iterator
    TriaObjects<TriaObject<1> >::next_free_single_line (const dealii::Triangulation<dim,spacedim> &tria)
    {
				       // TODO: Think of a way to ensure that we are using the correct triangulation, i.e. the one containing *this.
      
      int pos=next_free_single,
	 last=used.size()-1;
      if (!reverse_order_next_free_single)
	{
					   // first sweep forward, only use
					   // really single slots, do not use
					   // pair slots
	  for (; pos<last; ++pos)
	    if (!used[pos])
	      if (used[++pos])
		{
						   // this was a single slot
		  pos-=1;
		  break;
		}
	  if (pos>=last)
	    {
	      reverse_order_next_free_single=true;
	      next_free_single=used.size()-1;
	      pos=used.size()-1;
	    }
	  else
	    next_free_single=pos+1;
	}
      
      if(reverse_order_next_free_single)
	{
					   // second sweep, use all slots, even
					   // in pairs
	  for(;pos>=0;--pos)
	    if (!used[pos])
	      break;
	  if (pos>0)
	    next_free_single=pos-1;
	  else
					     // no valid single line anymore
	    return tria.end_line();
	}

      return typename dealii::Triangulation<dim,spacedim>::raw_line_iterator(&tria,0,pos);
    }
    
	  

    template <>
    template <int dim, int spacedim>
    typename dealii::Triangulation<dim,spacedim>::raw_line_iterator
    TriaObjects<TriaObject<1> >::next_free_pair_line (const dealii::Triangulation<dim,spacedim> &tria)
    {
				       // TODO: Think of a way to ensure that we are using the correct triangulation, i.e. the one containing *this.
      
      int pos=next_free_pair,
	 last=used.size()-1;
      for (; pos<last; ++pos)
	if (!used[pos])
	  if (!used[++pos])
	    {
					       // this was a pair slot
	      pos-=1;
	      break;
	    }
      if (pos>=last)
					 // no free slot
	return tria.end_line();
      else
	next_free_pair=pos+2;
      
      return typename dealii::Triangulation<dim,spacedim>::raw_line_iterator(&tria,0,pos);
    }


    
    template <>
    template <int dim, int spacedim>
    typename dealii::Triangulation<dim,spacedim>::raw_quad_iterator
    TriaObjects<TriaObject<1> >::next_free_single_quad (const dealii::Triangulation<dim,spacedim> &tria)
    {
      Assert(false, ExcWrongIterator("quad","line"));
      return tria.end_quad();
    }



    template <>
    template <int dim, int spacedim>
    typename dealii::Triangulation<dim,spacedim>::raw_quad_iterator
    TriaObjects<TriaObject<1> >::next_free_pair_quad (const dealii::Triangulation<dim,spacedim> &tria)
    {
      Assert(false, ExcWrongIterator("quad","line"));
      return tria.end_quad();
    }


    
    template <>
    template <int dim, int spacedim>
    typename dealii::Triangulation<dim,spacedim>::raw_quad_iterator
    TriaObjects<TriaObject<2> >::next_free_single_quad (const dealii::Triangulation<dim,spacedim> &tria)
    {
				       // TODO: Think of a way to ensure that we are using the correct triangulation, i.e. the one containing *this.
      
      int pos=next_free_single,
	 last=used.size()-1;
      if (!reverse_order_next_free_single)
	{
					   // first sweep forward, only use
					   // really single slots, do not use
					   // pair slots
	  for (; pos<last; ++pos)
	    if (!used[pos])
	      if (used[++pos])
		{
						   // this was a single slot
		  pos-=1;
		  break;
		}
	  if (pos>=last)
	    {
	      reverse_order_next_free_single=true;
	      next_free_single=used.size()-1;
	      pos=used.size()-1;
	    }
	  else
	    next_free_single=pos+1;
	}
      
      if(reverse_order_next_free_single)
	{
					   // second sweep, use all slots, even
					   // in pairs
	  for(;pos>=0;--pos)
	    if (!used[pos])
	      break;
	  if (pos>0)
	    next_free_single=pos-1;
	  else
					     // no valid single quad anymore
	    return tria.end_quad();
	}

      return typename dealii::Triangulation<dim,spacedim>::raw_quad_iterator(&tria,0,pos);
    }
    
	  

    template <>
    template <int dim, int spacedim>
    typename dealii::Triangulation<dim,spacedim>::raw_quad_iterator
    TriaObjects<TriaObject<2> >::next_free_pair_quad (const dealii::Triangulation<dim,spacedim> &tria)
    {
				       // TODO: Think of a way to ensure that we are using the correct triangulation, i.e. the one containing *this.
      
      int pos=next_free_pair,
	 last=used.size()-1;
      for (; pos<last; ++pos)
	if (!used[pos])
	  if (!used[++pos])
	    {
					       // this was a pair slot
	      pos-=1;
	      break;
	    }
      if (pos>=last)
					 // no free slot
	return tria.end_quad();
      else
	next_free_pair=pos+2;
      
      return typename dealii::Triangulation<dim,spacedim>::raw_quad_iterator(&tria,0,pos);
    }


    template <>
    template <int dim, int spacedim>
    typename dealii::Triangulation<dim,spacedim>::raw_line_iterator
    TriaObjects<TriaObject<2> >::next_free_single_line (const dealii::Triangulation<dim,spacedim> &tria)
    {
      Assert(false, ExcWrongIterator("line","quad"));
      return tria.end_line();
    }



    template <>
    template <int dim, int spacedim>
    typename dealii::Triangulation<dim,spacedim>::raw_line_iterator
    TriaObjects<TriaObject<2> >::next_free_pair_line (const dealii::Triangulation<dim,spacedim> &tria)
    {
      Assert(false, ExcWrongIterator("line","quad"));
      return tria.end_line();
    }



    template <>
    template <int dim, int spacedim>
    typename dealii::Triangulation<dim,spacedim>::raw_hex_iterator
    TriaObjects<TriaObject<3> >::next_free_hex (const dealii::Triangulation<dim,spacedim> &tria,
					    const unsigned int               level)
    {
				       // TODO: Think of a way to ensure that we are using the correct triangulation, i.e. the one containing *this.
      
      int pos=next_free_pair,
	 last=used.size()-1;
      for (; pos<last; ++pos)
	if (!used[pos])
	  {
					     // this should be a pair slot
	      Assert(!used[pos+1], ExcInternalError());
	      break;
	    }
      if (pos>=last)
					 // no free slot
	return tria.end_hex();
      else
	next_free_pair=pos+2;
      
      return typename dealii::Triangulation<dim,spacedim>::raw_hex_iterator(&tria,level,pos);
    }


    


				     // explicit instantiations

#if deal_II_dimension > 1
    template dealii::Triangulation<deal_II_dimension>::raw_line_iterator
    TriaObjects<TriaObject<1> >::next_free_single_line(const dealii::Triangulation<deal_II_dimension> &);
    template dealii::Triangulation<deal_II_dimension>::raw_line_iterator
    TriaObjects<TriaObject<1> >::next_free_pair_line(const dealii::Triangulation<deal_II_dimension> &);
    template dealii::Triangulation<deal_II_dimension>::raw_quad_iterator
    TriaObjects<TriaObject<1> >::next_free_single_quad(const dealii::Triangulation<deal_II_dimension> &);
    template dealii::Triangulation<deal_II_dimension>::raw_quad_iterator
    TriaObjects<TriaObject<1> >::next_free_pair_quad(const dealii::Triangulation<deal_II_dimension> &);

    template dealii::Triangulation<deal_II_dimension>::raw_line_iterator
    TriaObjects<TriaObject<2> >::next_free_single_line(const dealii::Triangulation<deal_II_dimension> &);
    template dealii::Triangulation<deal_II_dimension>::raw_line_iterator
    TriaObjects<TriaObject<2> >::next_free_pair_line(const dealii::Triangulation<deal_II_dimension> &);
    template dealii::Triangulation<deal_II_dimension>::raw_quad_iterator
    TriaObjects<TriaObject<2> >::next_free_single_quad(const dealii::Triangulation<deal_II_dimension> &);
    template dealii::Triangulation<deal_II_dimension>::raw_quad_iterator
    TriaObjects<TriaObject<2> >::next_free_pair_quad(const dealii::Triangulation<deal_II_dimension> &);
#endif
#if deal_II_dimension == 3
    template dealii::Triangulation<deal_II_dimension>::raw_hex_iterator
    TriaObjects<TriaObject<3> >::next_free_hex(const dealii::Triangulation<deal_II_dimension> &, const unsigned int);
#endif
  }
}

DEAL_II_NAMESPACE_CLOSE

