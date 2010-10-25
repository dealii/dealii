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
#include <grid/tria_objects.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/tria_accessor.h>

#include <algorithm>
#include <functional>

    

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace Triangulation
  {    
    template<>
    void
    TriaObjects<TriaObject<1> >::reserve_space (const unsigned int new_lines_in_pairs,
				      const unsigned int new_lines_single)
    {
      Assert(new_lines_in_pairs%2==0, ExcInternalError());

      next_free_single=0;
      next_free_pair=0;
      reverse_order_next_free_single=false;
      
				       // count the number of lines, of
				       // unused single lines and of
				       // unused pairs of lines
      unsigned int n_lines=0;
      unsigned int n_unused_pairs=0;
      unsigned int n_unused_singles=0;
      for (unsigned int i=0; i<used.size(); ++i)
	{
	  if (used[i])
	    ++n_lines;
	  else if (i+1<used.size())
	    {
	      if (used[i+1])
		{
		  ++n_unused_singles;
		  if (next_free_single==0)
		    next_free_single=i;
		}
	      else
		{
		  ++n_unused_pairs;
		  if (next_free_pair==0)
		    next_free_pair=i;
		  ++i;
		}
	    }
	  else
	    ++n_unused_singles;
	}
      Assert(n_lines+2*n_unused_pairs+n_unused_singles==used.size(),
	     ExcInternalError());

				       // how many single lines are needed in
				       // addition to n_unused_singles?
      const int additional_single_lines=
	new_lines_single-n_unused_singles;

      unsigned int new_size=
	used.size() + new_lines_in_pairs - 2*n_unused_pairs;
      if (additional_single_lines>0)
	new_size+=additional_single_lines;

                                       // only allocate space if necessary
      if (new_size>cells.size()) 
        {
          cells.reserve (new_size);
          cells.insert (cells.end(),
			new_size-cells.size(),
			TriaObject<1> ());
  
          used.reserve (new_size);
          used.insert (used.end(),
		       new_size-used.size(),
		       false);
  
          user_flags.reserve (new_size);
          user_flags.insert (user_flags.end(),
			     new_size-user_flags.size(),
			     false);
      
          children.reserve (new_size);
          children.insert (children.end(),
			   new_size-children.size(),
			   -1);

          material_id.reserve (new_size);
          material_id.insert (material_id.end(),
			      new_size-material_id.size(),
			      255);

          user_data.reserve (new_size);
          user_data.insert (user_data.end(),
			    new_size-user_data.size(),
			    UserData());
	}

      if (n_unused_singles==0)
	{
	  next_free_single=new_size-1;
	  reverse_order_next_free_single=true;
	}
    }


    template<>
    void
    TriaObjects<TriaObject<2> >::reserve_space (const unsigned int new_quads_in_pairs,
				      const unsigned int new_quads_single)
    {
      Assert(new_quads_in_pairs%2==0, ExcInternalError());
      
      next_free_single=0;
      next_free_pair=0;
      reverse_order_next_free_single=false;

				       // count the number of lines, of
				       // unused single lines and of
				       // unused pairs of lines
      unsigned int n_quads=0;
      unsigned int n_unused_pairs=0;
      unsigned int n_unused_singles=0;
      for (unsigned int i=0; i<used.size(); ++i)
	{
	  if (used[i])
	    ++n_quads;
	  else if (i+1<used.size())
	    {
	      if (used[i+1])
		{
		  ++n_unused_singles;
		  if (next_free_single==0)
		    next_free_single=i;
		}
	      else
		{
		  ++n_unused_pairs;
		  if (next_free_pair==0)
		    next_free_pair=i;
		  ++i;
		}
	    }
	  else
	    ++n_unused_singles;
	}
      Assert(n_quads+2*n_unused_pairs+n_unused_singles==used.size(),
	     ExcInternalError());

				       // how many single quads are needed in
				       // addition to n_unused_quads?
      const int additional_single_quads=
	new_quads_single-n_unused_singles;

      unsigned int new_size=
	used.size() + new_quads_in_pairs - 2*n_unused_pairs;
      if (additional_single_quads>0)
	new_size+=additional_single_quads;

                                       // only allocate space if necessary
      if (new_size>cells.size())
        {
          cells.reserve (new_size);
          cells.insert (cells.end(),
			new_size-cells.size(),
			TriaObject<2> ());
      
          used.reserve (new_size);
          used.insert (used.end(),
		       new_size-used.size(),
		       false);
  
          user_flags.reserve (new_size);
          user_flags.insert (user_flags.end(),
			     new_size-user_flags.size(),
			     false);
  
          children.reserve (2*new_size);
          children.insert (children.end(),
			   2*new_size-children.size(),
			   -1);

	  refinement_cases.reserve (new_size);
	  refinement_cases.insert (refinement_cases.end(),
			       new_size - refinement_cases.size(),
			       RefinementCase<2>::no_refinement);
	  

          material_id.reserve (new_size);
          material_id.insert (material_id.end(),
			      new_size-material_id.size(),
			      255);

          user_data.reserve (new_size);
          user_data.insert (user_data.end(),
			    new_size-user_data.size(),
			    UserData());
        }

      if (n_unused_singles==0)
	{
	  next_free_single=new_size-1;
	  reverse_order_next_free_single=true;
	}
    }


    void
    TriaObjectsHex::reserve_space (const unsigned int new_hexes)
    {
      const unsigned int new_size = new_hexes +
                                    std::count_if (used.begin(),
                                                   used.end(),
                                                   std::bind2nd (std::equal_to<bool>(), true));

                                       // see above...
      if (new_size>cells.size())
        {
          cells.reserve (new_size);
          cells.insert (cells.end(),
			new_size-cells.size(),
			TriaObject<3> ());
      
          used.reserve (new_size);
          used.insert (used.end(),
		       new_size-used.size(),
		       false);
  
          user_flags.reserve (new_size);
          user_flags.insert (user_flags.end(),
			     new_size-user_flags.size(),
			     false);
  
          children.reserve (4*new_size);
          children.insert (children.end(),
			   4*new_size-children.size(),
			   -1);

          material_id.reserve (new_size);
          material_id.insert (material_id.end(),
			      new_size-material_id.size(),
			      255);

          user_data.reserve (new_size);
          user_data.insert (user_data.end(),
			    new_size-user_data.size(),
			    UserData());

          face_orientations.reserve (new_size * GeometryInfo<3>::faces_per_cell);
          face_orientations.insert (face_orientations.end(),
				    new_size * GeometryInfo<3>::faces_per_cell
				    - face_orientations.size(),
				    true);

	  refinement_cases.reserve (new_size);
	  refinement_cases.insert (refinement_cases.end(),
			       new_size-refinement_cases.size(),
			       RefinementCase<3>::no_refinement);

          face_flips.reserve (new_size * GeometryInfo<3>::faces_per_cell);
          face_flips.insert (face_flips.end(),
				    new_size * GeometryInfo<3>::faces_per_cell
				    - face_flips.size(),
				    false);
          face_rotations.reserve (new_size * GeometryInfo<3>::faces_per_cell);
          face_rotations.insert (face_rotations.end(),
				    new_size * GeometryInfo<3>::faces_per_cell
				    - face_rotations.size(),
				    false);
        }
      next_free_single=next_free_pair=0;
    }


    void
    TriaObjectsQuad3D::reserve_space (const unsigned int new_quads_in_pairs,
				      const unsigned int new_quads_single)
    {
      Assert(new_quads_in_pairs%2==0, ExcInternalError());

      next_free_single=0;
      next_free_pair=0;
      reverse_order_next_free_single=false;

				       // count the number of lines, of unused
				       // single lines and of unused pairs of
				       // lines
      unsigned int n_quads=0;
      unsigned int n_unused_pairs=0;
      unsigned int n_unused_singles=0;
      for (unsigned int i=0; i<used.size(); ++i)
	{
	  if (used[i])
	    ++n_quads;
	  else if (i+1<used.size())
	    {
	      if (used[i+1])
		{
		  ++n_unused_singles;
		  if (next_free_single==0)
		    next_free_single=i;
		}
	      else
		{
		  ++n_unused_pairs;
		  if (next_free_pair==0)
		    next_free_pair=i;
		  ++i;
		}
	    }
	  else
	    ++n_unused_singles;
	}
      Assert(n_quads+2*n_unused_pairs+n_unused_singles==used.size(),
	     ExcInternalError());

				       // how many single quads are needed in
				       // addition to n_unused_quads?
      const int additional_single_quads=
	new_quads_single-n_unused_singles;

      unsigned int new_size=
	used.size() + new_quads_in_pairs - 2*n_unused_pairs;
      if (additional_single_quads>0)
	new_size+=additional_single_quads;
      
                                       // see above...
      if (new_size>cells.size())
        {
					   // reseve space for the base class
	  TriaObjects<TriaObject<2> >::reserve_space(new_quads_in_pairs,new_quads_single);
					   // reserve the field of the derived
					   // class
          line_orientations.reserve (new_size * GeometryInfo<2>::lines_per_cell);
          line_orientations.insert (line_orientations.end(),
				    new_size * GeometryInfo<2>::lines_per_cell
				    - line_orientations.size(),
				    true);
        }

      if (n_unused_singles==0)
	{
	  next_free_single=new_size-1;
	  reverse_order_next_free_single=true;
	}
    }
    

    template<>
    void
    TriaObjects<TriaObject<1> >::monitor_memory (const unsigned int) const
    {
                                       // check that we have not allocated
                                       // too much memory. note that bool
                                       // vectors allocate their memory in
                                       // chunks of whole integers, so
                                       // they may over-allocate by up to
                                       // as many elements as an integer
                                       // has bits
      Assert (cells.size() <=
	      cells.capacity() + DEAL_II_MIN_VECTOR_CAPACITY,
              ExcMemoryWasted ("lines",
                               cells.size(), cells.capacity()));
      Assert (children.size() <=
	      children.capacity() + DEAL_II_MIN_VECTOR_CAPACITY,
              ExcMemoryWasted ("children",
                               children.size(), children.capacity()));
      Assert (used.size() <= used.capacity() + sizeof(int)*8 ||
              used.size()<DEAL_II_MIN_BOOL_VECTOR_CAPACITY,
              ExcMemoryWasted ("used",
                               used.size(), used.capacity()));
      Assert (user_flags.size() <= user_flags.capacity() + sizeof(int)*8 ||
              user_flags.size()<DEAL_II_MIN_BOOL_VECTOR_CAPACITY,
              ExcMemoryWasted ("user_flags",
                               user_flags.size(), user_flags.capacity()));
      Assert (cells.size() == used.size(),
              ExcMemoryInexact (cells.size(), used.size()));
      Assert (cells.size() == user_flags.size(),
              ExcMemoryInexact (cells.size(), user_flags.size()));
      Assert (cells.size() == children.size(),
              ExcMemoryInexact (cells.size(), children.size()));
      Assert (cells.size() == material_id.size(),
              ExcMemoryInexact (cells.size(), material_id.size()));
      Assert (cells.size() == user_data.size(),
              ExcMemoryInexact (cells.size(), user_data.size()));
    }


    template<>
    void
    TriaObjects<TriaObject<2> >::monitor_memory (const unsigned int) const
    {
                                       // check that we have not allocated
                                       // too much memory. note that bool
                                       // vectors allocate their memory in
                                       // chunks of whole integers, so
                                       // they may over-allocate by up to
                                       // as many elements as an integer
                                       // has bits
      Assert (cells.size() <=
	      cells.capacity() + DEAL_II_MIN_VECTOR_CAPACITY,
              ExcMemoryWasted ("quads",
                               cells.size(), cells.capacity()));
      Assert (children.size() <=
	      children.capacity() + DEAL_II_MIN_VECTOR_CAPACITY,
              ExcMemoryWasted ("children",
                               children.size(), children.capacity()));
      Assert (used.size() <= used.capacity() + sizeof(int)*8 ||
              used.size()<DEAL_II_MIN_BOOL_VECTOR_CAPACITY,
              ExcMemoryWasted ("used",
                               used.size(), used.capacity()));
      Assert (user_flags.size() <= user_flags.capacity() + sizeof(int)*8 ||
              user_flags.size()<DEAL_II_MIN_BOOL_VECTOR_CAPACITY,
              ExcMemoryWasted ("user_flags",
                               user_flags.size(), user_flags.capacity()));
      Assert (cells.size() == used.size(),
              ExcMemoryInexact (cells.size(), used.size()));
      Assert (cells.size() == user_flags.size(),
              ExcMemoryInexact (cells.size(), user_flags.size()));
      Assert (2*cells.size() == children.size(),
              ExcMemoryInexact (cells.size(), children.size()));
      Assert (cells.size() == refinement_cases.size(),
              ExcMemoryInexact (cells.size(), refinement_cases.size()));
      Assert (cells.size() == material_id.size(),
              ExcMemoryInexact (cells.size(), material_id.size()));
      Assert (cells.size() == user_data.size(),
              ExcMemoryInexact (cells.size(), user_data.size()));
    }


    void
    TriaObjectsHex::monitor_memory (const unsigned int) const
    {
                                       // check that we have not allocated
                                       // too much memory. note that bool
                                       // vectors allocate their memory in
                                       // chunks of whole integers, so
                                       // they may over-allocate by up to
                                       // as many elements as an integer
                                       // has bits
      Assert (cells.size() <=
	      cells.capacity() + DEAL_II_MIN_VECTOR_CAPACITY,
              ExcMemoryWasted ("hexes",
                               cells.size(), cells.capacity()));
      Assert (children.size() <=
	      children.capacity() + DEAL_II_MIN_VECTOR_CAPACITY,
              ExcMemoryWasted ("children",
                               children.size(), children.capacity()));
      Assert (used.size() <= used.capacity() + sizeof(int)*8 ||
              used.size()<DEAL_II_MIN_BOOL_VECTOR_CAPACITY,
              ExcMemoryWasted ("used",
                               used.size(), used.capacity()));
      Assert (user_flags.size() <= user_flags.capacity() + sizeof(int)*8 ||
              user_flags.size()<DEAL_II_MIN_BOOL_VECTOR_CAPACITY,
              ExcMemoryWasted ("user_flags",
                               user_flags.size(), user_flags.capacity()));
      Assert (cells.size() == used.size(),
              ExcMemoryInexact (cells.size(), used.size()));
      Assert (cells.size() == user_flags.size(),
              ExcMemoryInexact (cells.size(), user_flags.size()));
      Assert (4*cells.size() == children.size(),
              ExcMemoryInexact (cells.size(), children.size()));
      Assert (cells.size() == material_id.size(),
              ExcMemoryInexact (cells.size(), material_id.size()));
      Assert (cells.size() == user_data.size(),
              ExcMemoryInexact (cells.size(), user_data.size()));
      Assert (cells.size() * GeometryInfo<3>::faces_per_cell
              == face_orientations.size(),
              ExcMemoryInexact (cells.size() * GeometryInfo<3>::faces_per_cell,
                                face_orientations.size()));
      Assert (cells.size() * GeometryInfo<3>::faces_per_cell
              == face_flips.size(),
              ExcMemoryInexact (cells.size() * GeometryInfo<3>::faces_per_cell,
                                face_flips.size()));
      Assert (cells.size() * GeometryInfo<3>::faces_per_cell
              == face_rotations.size(),
              ExcMemoryInexact (cells.size() * GeometryInfo<3>::faces_per_cell,
                                face_rotations.size()));
    }


    void
    TriaObjectsQuad3D::monitor_memory (const unsigned int) const
    {
                                       // check that we have not allocated
                                       // too much memory. note that bool
                                       // vectors allocate their memory in
                                       // chunks of whole integers, so
                                       // they may over-allocate by up to
                                       // as many elements as an integer
                                       // has bits
      Assert (cells.size() * GeometryInfo<2>::lines_per_cell
              == line_orientations.size(),
              ExcMemoryInexact (cells.size() * GeometryInfo<2>::lines_per_cell,
                                line_orientations.size()));
      TriaObjects<TriaObject<2> >::monitor_memory (3);
      
    }


    template <typename G>
    void
    TriaObjects<G>::clear()
    {
      cells.clear();
      children.clear();
      refinement_cases.clear();
      used.clear();
      user_flags.clear();
      material_id.clear();
      user_data.clear();
      user_data_type = data_unknown;
    }
    

    void
    TriaObjectsHex::clear()
    {
      TriaObjects<TriaObject<3> >::clear();
      face_orientations.clear();
      face_flips.clear();
      face_rotations.clear();
    }


    void
    TriaObjectsQuad3D::clear()
    {
      TriaObjects<TriaObject<2> >::clear();
      line_orientations.clear();
    }
    
    
    template<typename G>
    unsigned int
    TriaObjects<G>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (cells) +
              MemoryConsumption::memory_consumption (children) +
              MemoryConsumption::memory_consumption (used) +
              MemoryConsumption::memory_consumption (user_flags) +
              MemoryConsumption::memory_consumption (material_id) +
              MemoryConsumption::memory_consumption (refinement_cases) +
	      user_data.capacity() * sizeof(UserData) + sizeof(user_data));
    }
  

    unsigned int
    TriaObjectsHex::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (face_orientations) +
	      MemoryConsumption::memory_consumption (face_flips) +
	      MemoryConsumption::memory_consumption (face_rotations) +
	      TriaObjects<TriaObject<3> >::memory_consumption() );
    }


    unsigned int
    TriaObjectsQuad3D::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (line_orientations) +
	      this->TriaObjects<TriaObject<2> >::memory_consumption() );
    }

    

				     // explicit instantiations

    template class TriaObjects<TriaObject<1> >;
    template class TriaObjects<TriaObject<2> >;
  }
}

DEAL_II_NAMESPACE_CLOSE

