/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1999 */


#include <grid/tria.h>
#include <grid/tria_levels.h>
#include <grid/tria_boundary.h>

#include <algorithm>
#include <numeric>


/**
 * Single out some functions which are needed by all dimensions, but
 * which are not template. They thus have the same name and when we
 * try to link with the libraries for different dimensions at the same
 * time, we get linker errors for functions defined more than once. By
 * putting these functions in a single file, the linker is allowed to
 * use it only once and throw away all other versions of this file in
 * the other libraries.
 */




bool SubCellData::check_consistency (const unsigned int dim) const
{
  switch (dim) 
    {
      case 1:
	    return ((boundary_lines.size() == 0) &&
		    (boundary_quads.size() == 0));
      case 2:
	    return (boundary_quads.size() == 0);
    };
  return true;
};

		    


void TriangulationLevel<0>::reserve_space (const unsigned int total_cells,
					   const unsigned int dimension)
{
				   // we need space for total_cells
				   // cells. Maybe we have more already
				   // with those cells which are unused,
				   // so only allocate new space if needed.
				   //
				   // note that all arrays should have equal
				   // sizes (checked by #monitor_memory#
  if (total_cells > refine_flags.size()) 
    {
      refine_flags.reserve (total_cells);
      refine_flags.insert (refine_flags.end(),
			   total_cells - refine_flags.size(),
			   false);
      
      coarsen_flags.reserve (total_cells);
      coarsen_flags.insert (coarsen_flags.end(),
			    total_cells - coarsen_flags.size(),
			    false);
      
      neighbors.reserve (total_cells*(2*dimension));
      neighbors.insert (neighbors.end(),
			total_cells*(2*dimension) - neighbors.size(),
			make_pair(-1,-1));
    };
};




void TriangulationLevel<0>::monitor_memory (const unsigned int true_dimension) const
{
//  Assert (refine_flags.size() == refine_flags.capacity() ||
//	  refine_flags.size()<256,
//	  ExcMemoryWasted ("refine_flags",
//			   refine_flags.size(), refine_flags.capacity()));
//  Assert (coarsen_flags.size() == coarsen_flags.capacity() ||
//	  coarsen_flags.size()<256,
//	  ExcMemoryWasted ("coarsen_flags",
//			   coarsen_flags.size(), coarsen_flags.capacity()));
//  Assert (neighbors.size() ==  neighbors.capacity() ||
//	  neighbors.size()<256,
//	  ExcMemoryWasted ("neighbors",
//			   neighbors.size(), neighbors.capacity()));
  Assert (2*true_dimension*refine_flags.size() == neighbors.size(),
	  ExcMemoryInexact (refine_flags.size(), neighbors.size()));
  Assert (2*true_dimension*coarsen_flags.size() == neighbors.size(),
	  ExcMemoryInexact (coarsen_flags.size(), neighbors.size()));
};



void TriangulationLevel<1>::reserve_space (const unsigned int new_lines)
{
  const unsigned int new_size = new_lines +
				count_if (lines.used.begin(),
					  lines.used.end(),
					  bind2nd (equal_to<bool>(), true));

				   // same as in #reserve_space<0>#: only
				   // allocate space if necessary
  if (new_size>lines.lines.size()) 
    {
//  cout << "  lines: pre: siz=" << lines.lines.size() << ", cap=" << lines.lines.capacity();
      lines.lines.reserve (new_size);
//  cout << " inter: siz=" << lines.lines.size() << ", cap=" << lines.lines.capacity()
//       << " (newsize=" << new_size << ")";
      lines.lines.insert (lines.lines.end(), new_size-lines.lines.size(), Line());
//  cout << " post: siz=" << lines.lines.size() << ", cap=" << lines.lines.capacity() << endl;
  
//  cout << "  used : pre: siz=" << lines.used.size() << ", cap=" << lines.used.capacity();
      lines.used.reserve (new_size);
//  cout << " inter: siz=" << lines.used.size() << ", cap=" << lines.used.capacity()
//       << " (newsize=" << new_size << ")";
      lines.used.insert (lines.used.end(), new_size-lines.used.size(), false);
//  cout << " post: siz=" << lines.used.size() << ", cap=" << lines.used.capacity() << endl;
  
      lines.user_flags.reserve (new_size);
      lines.user_flags.insert (lines.user_flags.end(),
			       new_size-lines.user_flags.size(), false);
      
      lines.children.reserve (new_size);
      lines.children.insert (lines.children.end(), new_size-lines.children.size(),
			     -1);

      lines.material_id.reserve (new_size);
      lines.material_id.insert (lines.material_id.end(),
				new_size-lines.material_id.size(),
				255);

      lines.user_pointers.reserve (new_size);
      lines.user_pointers.insert (lines.user_pointers.end(),
				  new_size-lines.user_pointers.size(), 0);
    };
};




void TriangulationLevel<1>::monitor_memory (const unsigned int true_dimension) const
{
//  Assert (lines.lines.size() == lines.lines.capacity() ||
//	  lines.lines.size()<256,
//	  ExcMemoryWasted ("lines",
//			   lines.lines.size(), lines.lines.capacity()));
//  Assert (lines.children.size() == lines.children.capacity() ||
//	  lines.children.size()<256,
//	  ExcMemoryWasted ("children",
//			   lines.children.size(), lines.children.capacity()));
//  Assert (lines.used.size() == lines.used.capacity() ||
//	  lines.used.size()<256,
//	  ExcMemoryWasted ("used",
//			   lines.used.size(), lines.used.capacity()));
//  Assert (lines.user_flags.size() == lines.user_flags.capacity() ||
//	  lines.user_flags.size()<256,
//	  ExcMemoryWasted ("user_flags",
//			   lines.user_flags.size(), lines.user_flags.capacity()));
  Assert (lines.lines.size() == lines.used.size(),
	  ExcMemoryInexact (lines.lines.size(), lines.used.size()));
  Assert (lines.lines.size() == lines.user_flags.size(),
	  ExcMemoryInexact (lines.lines.size(), lines.user_flags.size()));
  Assert (lines.lines.size() == lines.children.size(),
	  ExcMemoryInexact (lines.lines.size(), lines.children.size()));
  Assert (lines.lines.size() == lines.material_id.size(),
	  ExcMemoryInexact (lines.lines.size(), lines.material_id.size()));
  Assert (lines.lines.size() == lines.user_pointers.size(),
	  ExcMemoryInexact (lines.lines.size(), lines.user_pointers.size()));

  TriangulationLevel<0>::monitor_memory (true_dimension);
};



void TriangulationLevel<2>::reserve_space (const unsigned int new_quads)
{
  const unsigned int new_size = new_quads +
				count_if (quads.used.begin(),
					  quads.used.end(),
					  bind2nd (equal_to<bool>(), true));

				   // see above...
  if (new_size>quads.quads.size())
    {
      quads.quads.reserve (new_size);
      quads.quads.insert (quads.quads.end(), new_size-quads.quads.size(), Quad());
      
      quads.used.reserve (new_size);
      quads.used.insert (quads.used.end(), new_size-quads.used.size(), false);
  
      quads.user_flags.reserve (new_size);
      quads.user_flags.insert (quads.user_flags.end(),
			       new_size-quads.user_flags.size(), false);
  
      quads.children.reserve (new_size);
      quads.children.insert (quads.children.end(), new_size-quads.children.size(),
			     -1);

      quads.material_id.reserve (new_size);
      quads.material_id.insert (quads.material_id.end(),
				new_size-quads.material_id.size(),
				255);

      quads.user_pointers.reserve (new_size);
      quads.user_pointers.insert (quads.user_pointers.end(),
				  new_size-quads.user_pointers.size(), 0);
    };
};



void TriangulationLevel<2>::monitor_memory (const unsigned int true_dimension) const
{
//  Assert (quads.quads.size() == quads.quads.capacity() ||
//	  quads.quads.size()<256,
//	  ExcMemoryWasted ("quads",
//			   quads.quads.size(), quads.quads.capacity()));
//  Assert (quads.children.size() == quads.children.capacity() ||
//	  quads.children.size()<256,
//	  ExcMemoryWasted ("children",
//			   quads.children.size(), quads.children.capacity()));
//  Assert (quads.used.size() == quads.used.capacity() ||
//	  quads.used.size()<256,
//	  ExcMemoryWasted ("used",
//			   quads.used.size(), quads.used.capacity()));
//  Assert (quads.user_flags.size() == quads.user_flags.capacity() ||
//	  quads.user_flags.size()<256,
//	  ExcMemoryWasted ("user_flags",
//			   quads.user_flags.size(), quads.user_flags.capacity()));
  Assert (quads.quads.size() == quads.used.size(),
	  ExcMemoryInexact (quads.quads.size(), quads.used.size()));
  Assert (quads.quads.size() == quads.user_flags.size(),
	  ExcMemoryInexact (quads.quads.size(), quads.user_flags.size()));
  Assert (quads.quads.size() == quads.children.size(),
	  ExcMemoryInexact (quads.quads.size(), quads.children.size()));
  Assert (quads.quads.size() == quads.material_id.size(),
	  ExcMemoryInexact (quads.quads.size(), quads.material_id.size()));
  Assert (quads.quads.size() == quads.user_pointers.size(),
	  ExcMemoryInexact (quads.quads.size(), quads.user_pointers.size()));

  TriangulationLevel<1>::monitor_memory (true_dimension);
};



void TriangulationLevel<3>::reserve_space (const unsigned int new_hexes)
{
  const unsigned int new_size = new_hexes +
				count_if (hexes.used.begin(),
					  hexes.used.end(),
					  bind2nd (equal_to<bool>(), true));

				   // see above...
  if (new_size>hexes.hexes.size())
    {
      hexes.hexes.reserve (new_size);
      hexes.hexes.insert (hexes.hexes.end(), new_size-hexes.hexes.size(), Hexahedron());
      
      hexes.used.reserve (new_size);
      hexes.used.insert (hexes.used.end(), new_size-hexes.used.size(), false);
  
      hexes.user_flags.reserve (new_size);
      hexes.user_flags.insert (hexes.user_flags.end(),
			       new_size-hexes.user_flags.size(), false);
  
      hexes.children.reserve (new_size);
      hexes.children.insert (hexes.children.end(), new_size-hexes.children.size(),
			     -1);

      hexes.material_id.reserve (new_size);
      hexes.material_id.insert (hexes.material_id.end(),
				new_size-hexes.material_id.size(),
				255);

      hexes.user_pointers.reserve (new_size);
      hexes.user_pointers.insert (hexes.user_pointers.end(),
				  new_size-hexes.user_pointers.size(), 0);
    };
};



void TriangulationLevel<3>::monitor_memory (const unsigned int true_dimension) const
{
//  Assert (hexes.hexes.size() == hexes.hexes.capacity() ||
//	  hexes.hexes.size()<256,
//	  ExcMemoryWasted ("hexes",
//			   hexes.hexes.size(), hexes.hexes.capacity()));
//  Assert (hexes.children.size() == hexes.children.capacity() ||
//	  hexes.children.size()<256,
//	  ExcMemoryWasted ("children",
//			   hexes.children.size(), hexes.children.capacity()));
//  Assert (hexes.used.size() == hexes.used.capacity() ||
//	  hexes.used.size()<256,
//	  ExcMemoryWasted ("used",
//			   hexes.used.size(), hexes.used.capacity()));
//  Assert (hexes.user_flags.size() == hexes.user_flags.capacity() ||
//	  hexes.user_flags.size()<256,
//	  ExcMemoryWasted ("user_flags",
//			   hexes.user_flags.size(), hexes.user_flags.capacity()));
  Assert (hexes.hexes.size() == hexes.used.size(),
	  ExcMemoryInexact (hexes.hexes.size(), hexes.used.size()));
  Assert (hexes.hexes.size() == hexes.user_flags.size(),
	  ExcMemoryInexact (hexes.hexes.size(), hexes.user_flags.size()));
  Assert (hexes.hexes.size() == hexes.children.size(),
	  ExcMemoryInexact (hexes.hexes.size(), hexes.children.size()));
  Assert (hexes.hexes.size() == hexes.material_id.size(),
	  ExcMemoryInexact (hexes.hexes.size(), hexes.material_id.size()));
  Assert (hexes.hexes.size() == hexes.user_pointers.size(),
	  ExcMemoryInexact (hexes.hexes.size(), hexes.user_pointers.size()));

  TriangulationLevel<2>::monitor_memory (true_dimension);
};



TriaNumberCache<1>::TriaNumberCache () :
		n_lines (0),
		n_active_lines (0) 
						 // all other fields are
						 // default constructed
{};




TriaNumberCache<2>::TriaNumberCache () :
		n_quads (0),
		n_active_quads (0) 
						 // all other fields are
						 // default constructed
{};




TriaNumberCache<3>::TriaNumberCache () :
		n_hexes (0),
		n_active_hexes (0) 
						 // all other fields are
						 // default constructed
{};
