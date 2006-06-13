//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/memory_consumption.h>
#include <grid/tria_objects.h>

#include <algorithm>
#include <functional>


namespace internal
{
  namespace Triangulation
  {
    
    template<>
    void
    TriaObjects<Line>::reserve_space (const unsigned int new_lines)
    {
      const unsigned int new_size = new_lines +
                                    std::count_if (used.begin(),
                                                   used.end(),
                                                   std::bind2nd (std::equal_to<bool>(), true));

                                       // only allocate space if necessary
      if (new_size>cells.size()) 
        {
          cells.reserve (new_size);
          cells.insert (cells.end(),
			new_size-cells.size(),
			Line());
  
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

          user_pointers.reserve (new_size);
          user_pointers.insert (user_pointers.end(),
				new_size-user_pointers.size(),
				0);
        };
    }


    template<>
    void
    TriaObjects<Quad>::reserve_space (const unsigned int new_quads)
    {
      const unsigned int new_size = new_quads +
                                    std::count_if (used.begin(),
                                                   used.end(),
                                                   std::bind2nd (std::equal_to<bool>(), true));

                                       // see above...
      if (new_size>cells.size())
        {
          cells.reserve (new_size);
          cells.insert (cells.end(),
			new_size-cells.size(),
			Quad());
      
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

          user_pointers.reserve (new_size);
          user_pointers.insert (user_pointers.end(),
				new_size-user_pointers.size(),
				0);
        };
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
			Hexahedron());
      
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

          user_pointers.reserve (new_size);
          user_pointers.insert (user_pointers.end(),
				new_size-user_pointers.size(),
				0);

          face_orientations.reserve (new_size * GeometryInfo<3>::faces_per_cell);
          face_orientations.insert (face_orientations.end(),
				    new_size * GeometryInfo<3>::faces_per_cell
				    - face_orientations.size(),
				    true);
        };
    }


    template<>
    void
    TriaObjects<Line>::monitor_memory (const unsigned int) const
    {
                                       // check that we have not allocated
                                       // too much memory. note that bool
                                       // vectors allocate their memory in
                                       // chunks of whole integers, so
                                       // they may over-allocate by up to
                                       // as many elements as an integer
                                       // has bits
      Assert (cells.size() == cells.capacity() ||
              cells.size()<DEAL_II_MIN_VECTOR_CAPACITY,
              ExcMemoryWasted ("lines",
                               cells.size(), cells.capacity()));
      Assert (children.size() == children.capacity() ||
              children.size()<DEAL_II_MIN_VECTOR_CAPACITY,
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
      Assert (cells.size() == user_pointers.size(),
              ExcMemoryInexact (cells.size(), user_pointers.size()));
    }


    template<>
    void
    TriaObjects<Quad>::monitor_memory (const unsigned int) const
    {
                                       // check that we have not allocated
                                       // too much memory. note that bool
                                       // vectors allocate their memory in
                                       // chunks of whole integers, so
                                       // they may over-allocate by up to
                                       // as many elements as an integer
                                       // has bits
      Assert (cells.size() == cells.capacity() ||
              cells.size()<DEAL_II_MIN_VECTOR_CAPACITY,
              ExcMemoryWasted ("quads",
                               cells.size(), cells.capacity()));
      Assert (children.size() == children.capacity() ||
              children.size()<DEAL_II_MIN_VECTOR_CAPACITY,
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
      Assert (cells.size() == user_pointers.size(),
              ExcMemoryInexact (cells.size(), user_pointers.size()));
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
      Assert (cells.size() == cells.capacity() ||
              cells.size()<DEAL_II_MIN_VECTOR_CAPACITY,
              ExcMemoryWasted ("hexes",
                               cells.size(), cells.capacity()));
      Assert (children.size() == children.capacity() ||
              children.size()<DEAL_II_MIN_VECTOR_CAPACITY,
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
      Assert (cells.size() == user_pointers.size(),
              ExcMemoryInexact (cells.size(), user_pointers.size()));
      Assert (cells.size() * GeometryInfo<3>::faces_per_cell
              == face_orientations.size(),
              ExcMemoryInexact (cells.size() * GeometryInfo<3>::faces_per_cell,
                                face_orientations.size()));
    }


    template <typename G>
    void
    TriaObjects<G>::clear()
    {
      cells.clear();
      children.clear();
      used.clear();
      user_flags.clear();
      material_id.clear();
      user_pointers.clear();
    }
    

    void
    TriaObjectsHex::clear()
    {
      TriaObjects<Hexahedron>::clear();
      face_orientations.clear();
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
              MemoryConsumption::memory_consumption (user_pointers));
    }
  

    unsigned int
    TriaObjectsHex::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (face_orientations) +
	      this->TriaObjects<Hexahedron>::memory_consumption() );
    }

				     // explicit instantiations

    template class TriaObjects<Line>;

#if deal_II_dimension >=2
    template class TriaObjects<Quad>;
#endif
    
  }
}
