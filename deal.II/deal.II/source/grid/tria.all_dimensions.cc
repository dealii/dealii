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




bool SubCellData::check_consistency (const unsigned int dim) const {
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
					   const unsigned int dimension) {
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




void TriangulationLevel<0>::monitor_memory (const unsigned int true_dimension) const {
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



template <typename T>
void TriangulationLevel<0>::write_raw_vector (const unsigned int  magic_number1,
					      const vector<T>    &v,
					      const unsigned int  magic_number2,
					      ostream            &out) 
{
  AssertThrow (out, ExcIO());
  
  out << magic_number1 << ' ' << v.size() << '[';
  out.write (reinterpret_cast<const char*>(v.begin()),
	     reinterpret_cast<const char*>(v.end())
	     - reinterpret_cast<const char*>(v.begin()));
  out << ']' << magic_number2 << endl;

  AssertThrow (out, ExcIO());
};

  
  
void TriangulationLevel<0>::write_raw_vector (const unsigned int  magic_number1,
					      const vector<bool> &v,
					      const unsigned int  magic_number2,
					      ostream            &out) {
  const unsigned int N = v.size();
  unsigned char *flags = new unsigned char[N/8+1];
  for (unsigned int i=0; i<N/8+1; ++i) flags[i]=0;
  
  for (unsigned int position=0; position<N; ++position)
    flags[position/8] |= (v[position] ? (1<<(position%8)) : 0);

  AssertThrow (out, ExcIO());
  
				   // format:
				   // 0. magic number
				   // 1. number of flags
				   // 2. the flags
				   // 3. magic number
  out << magic_number1 << ' ' << N << '[';
  for (unsigned int i=0; i<N/8+1; ++i) 
    out << static_cast<unsigned int>(flags[i]) << " ";
  
  out << ']' << magic_number2 << endl;
  
  delete[] flags;

  AssertThrow (out, ExcIO());
};



template <typename T>
void TriangulationLevel<0>::write_rle_vector (const unsigned int  magic_number1,
					      const vector<T>    &v,
					      const unsigned int  magic_number2,
					      ostream            &out) 
{
  AssertThrow (out, ExcIO());

				   // store the position of the last
				   // break in the data stream from
				   // which on nothing was yet written
				   // to #out#
  typename vector<T>::const_iterator last_major_break = v.begin();
				   // pointer to an element against
				   // which we want to compare following
				   // elements
  typename vector<T>::const_iterator comparator;
				   // name says all
  typename vector<T>::const_iterator present_element;
				   // same here
  typename vector<T>::const_iterator end_of_vector = v.end();
  
  out << magic_number1 << ' ' << v.size() << '[';

  while (last_major_break < end_of_vector)
    {
				       // from the present position onward:
				       // find how many elements are equal
      comparator = last_major_break;

      while (true)
	{
	  present_element = comparator+1;
	  while ((present_element != end_of_vector) &&
		 (*present_element == *comparator)    &&
		 (present_element-comparator < 255))
	    ++present_element;

					   // now present_element points to the
					   // first element which is not equal
					   // to #comparator# or alternatively
					   // to #end_of_vector# or to
					   // #comparator+255#
					   //
					   // if #present_element# is
					   // #comparator+1#, i.e. the
					   // sequence of equal
					   // elements consisted of
					   // only one elements then
					   // discard this sequence
					   // and find the next one
					   //
					   // otherwise leave the loop
	  if ((present_element != end_of_vector) &&
	      (present_element < comparator+1))
	    comparator = present_element;
	  else
	    break;
	};

				       // if the loop broke because we have
				       // reached the end of the vector:
				       // maybe set comparator to
				       // present_element
      if ((present_element == end_of_vector) &&
	  (present_element < comparator+1))
	comparator = present_element;
      
				       // now comparator points to the start
				       // of the sequence of equal elements
				       // and present_element points past its
				       // end; they may be the same if the
				       // length of the sequence of unequal
				       // elements was 255
				       //
				       // now do the following: write out
				       // the elements of
				       // last_major_break...comparator,
				       // then write
				       // comparator...present_element
				       // if necessary; note that we need
				       // only write *one* element (that is why
				       // we do all this here)
      out << '<' << comparator-last_major_break << '>';
      if (comparator != last_major_break)
	out.write (reinterpret_cast<const char*>(last_major_break),
		   reinterpret_cast<const char*>(comparator)
		   - reinterpret_cast<const char*>(last_major_break));
      out << '<' << present_element-comparator << '>';
      if (present_element != comparator)
	out.write (reinterpret_cast<const char*>(comparator),
		   reinterpret_cast<const char*>(comparator+1)
		   - reinterpret_cast<const char*>(comparator));

      last_major_break = present_element;
    };
      
  out << ']' << magic_number2 << endl;

  AssertThrow (out, ExcIO());
};



void TriangulationLevel<0>::block_write (ostream &out) const
{
  write_raw_vector (0, refine_flags, 0, out);
  write_raw_vector (0, coarsen_flags, 0, out);
  write_raw_vector (0, neighbors, 0, out);
};





void TriangulationLevel<1>::reserve_space (const unsigned int new_lines) {
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




void TriangulationLevel<1>::monitor_memory (const unsigned int true_dimension) const {
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



void TriangulationLevel<1>::block_write (ostream &out) const
{
  TriangulationLevel<0>::block_write (out);

  write_raw_vector (0, lines.lines, 0, out);
  write_rle_vector (0, lines.children, 0, out);
  write_raw_vector (0, lines.used, 0, out);
  write_raw_vector (0, lines.user_flags, 0, out);
  write_rle_vector (0, lines.material_id, 0, out);
				   // note: user_pointers are not written
};




void TriangulationLevel<2>::reserve_space (const unsigned int new_quads) {
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



void TriangulationLevel<2>::monitor_memory (const unsigned int true_dimension) const {
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



void TriangulationLevel<2>::block_write (ostream &out) const
{
  TriangulationLevel<1>::block_write (out);

  write_raw_vector (0, quads.quads, 0, out);
  write_rle_vector (0, quads.children, 0, out);
  write_raw_vector (0, quads.used, 0, out);
  write_raw_vector (0, quads.user_flags, 0, out);
  write_rle_vector (0, quads.material_id, 0, out);
				   // note: user_pointers are not written
};






void TriangulationLevel<3>::reserve_space (const unsigned int new_hexes) {
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



void TriangulationLevel<3>::monitor_memory (const unsigned int true_dimension) const {
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



void TriangulationLevel<3>::block_write (ostream &out) const
{
  TriangulationLevel<2>::block_write (out);

  write_raw_vector (0, hexes.hexes, 0, out);
  write_rle_vector (0, hexes.children, 0, out);
  write_raw_vector (0, hexes.used, 0, out);
  write_raw_vector (0, hexes.user_flags, 0, out);
  write_rle_vector (0, hexes.material_id, 0, out);
				   // note: user_pointers are not written
};


