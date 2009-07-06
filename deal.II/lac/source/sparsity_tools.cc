//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/exceptions.h>
#include <lac/exceptions.h>
#include <lac/sparsity_pattern.h>
#include <lac/sparsity_tools.h>

#include <algorithm>

#ifdef DEAL_II_USE_METIS
// This is sorta stupid. what we really would like to do here is this:
//   extern "C" {
//   #  include <metis.h>
//   }
// The problem with this is that (i) metis.h declares a couple of
// functions with different exception specifications than are
// declared in standard header files (in particular the drand48
// function every one who's worked with Metis seems to know about);
// and (ii) metis.h may include <mpi.h> which if we run a C++
// compiler may include mpicxx.h -- but that leads to trouble
// because we have wrapped an
//   extern "C" {...}
// around the whole block that now also includes a C++ header :-(
//
// The correct solution, of course, would be if Metis would produce
// a new dot release after 10 years that would remove the unnecessary
// function prototypes, and that would wrap the contents into
//   #ifdef __C_PLUSPLUS__
//   extern "C" {
//   #endif
// But that appears to not be happening, and so we have to do the
// following hack where we just forward declare everything ourselves :-(

extern "C" 
{
// the following is from Metis-4.0/Lib/struct.h:
  typedef int idxtype;
// and this is from proto.h
  void
  METIS_PartGraphKway(int *, idxtype *, idxtype *,
		      idxtype *, idxtype *, int *,
		      int *, int *, int *, int *, idxtype *);
  void
  METIS_PartGraphRecursive(int *, idxtype *, idxtype *,
			   idxtype *, idxtype *, int *,
			   int *, int *, int *, int *, idxtype *); 
}
#endif


DEAL_II_NAMESPACE_OPEN

namespace SparsityTools
{
  
  void partition (const SparsityPattern     &sparsity_pattern,
		  const unsigned int         n_partitions,
		  std::vector<unsigned int> &partition_indices)
  {
    Assert (sparsity_pattern.n_rows()==sparsity_pattern.n_cols(),
	    ExcNotQuadratic());
    Assert (sparsity_pattern.is_compressed(),
	    SparsityPattern::ExcNotCompressed());

    Assert (n_partitions > 0, ExcInvalidNumberOfPartitions(n_partitions));
    Assert (partition_indices.size() == sparsity_pattern.n_rows(),
	    ExcInvalidArraySize (partition_indices.size(),
				 sparsity_pattern.n_rows()));

				     // check for an easy return
    if (n_partitions == 1)
      {
	std::fill_n (partition_indices.begin(), partition_indices.size(), 0U);
	return;
      }

				     // Make sure that METIS is actually
				     // installed and detected
#ifndef DEAL_II_USE_METIS
    AssertThrow (false, ExcMETISNotInstalled());
#else

				     // generate the data structures for
				     // METIS. Note that this is particularly
				     // simple, since METIS wants exactly our
				     // compressed row storage format. we only
				     // have to set up a few auxiliary arrays
    int
      n = static_cast<signed int>(sparsity_pattern.n_rows()),
      wgtflag = 0,                          // no weights on nodes or edges
      numflag = 0,                          // C-style 0-based numbering
      nparts  = static_cast<int>(n_partitions), // number of subdomains to create
      dummy;                                // the numbers of edges cut by the
				     //   resulting partition

				     // use default options for METIS
    int options[5] = { 0,0,0,0,0 };

				     // one more nuisance: we have to copy our
				     // own data to arrays that store signed
				     // integers :-(
    std::vector<idxtype> int_rowstart (sparsity_pattern.get_rowstart_indices(),
				       sparsity_pattern.get_rowstart_indices() +
				       sparsity_pattern.n_rows()+1);
    std::vector<idxtype> int_colnums (sparsity_pattern.get_column_numbers(),
				      sparsity_pattern.get_column_numbers()+
				      int_rowstart[sparsity_pattern.n_rows()]);

    std::vector<idxtype> int_partition_indices (sparsity_pattern.n_rows());

				     // Select which type of partitioning to
				     // create

				     // Use recursive if the number of
				     // partitions is less than or equal to 8
    if (n_partitions <= 8)
      METIS_PartGraphRecursive(&n, &int_rowstart[0], &int_colnums[0],
			       NULL, NULL,
			       &wgtflag, &numflag, &nparts, &options[0],
			       &dummy, &int_partition_indices[0]);
  
				     // Otherwise  use kway
    else
      METIS_PartGraphKway(&n, &int_rowstart[0], &int_colnums[0],
			  NULL, NULL,
			  &wgtflag, &numflag, &nparts, &options[0],
			  &dummy, &int_partition_indices[0]);

				     // now copy back generated indices into the
				     // output array
    std::copy (int_partition_indices.begin(),
	       int_partition_indices.end(),
	       partition_indices.begin());
#endif
  }



  void
  reorder_Cuthill_McKee (const SparsityPattern     &sparsity,
			 std::vector<unsigned int> &new_indices,
			 const std::vector<unsigned int> &starting_indices)
  {
    Assert (sparsity.n_rows() == sparsity.n_cols(),
	    ExcDimensionMismatch (sparsity.n_rows(), sparsity.n_cols()));
    Assert (sparsity.n_rows() == new_indices.size(),
	    ExcDimensionMismatch (sparsity.n_rows(), new_indices.size()));
    Assert (starting_indices.size() <= sparsity.n_rows(),
	    ExcMessage ("You can't specify more starting indices than there are rows"));
    for (unsigned int i=0; i<starting_indices.size(); ++i)
      Assert (starting_indices[i] < sparsity.n_rows(),
	      ExcMessage ("Invalid starting index"));
    
				     // store the indices of the dofs renumbered
				     // in the last round. Default to starting
				     // points
    std::vector<unsigned int> last_round_dofs (starting_indices);

				     // initialize the new_indices array with
				     // invalid values
    std::fill (new_indices.begin(), new_indices.end(),
	       numbers::invalid_unsigned_int);
    
				     // delete disallowed elements
    for (unsigned int i=0; i<last_round_dofs.size(); ++i)
      if ((last_round_dofs[i]==numbers::invalid_unsigned_int) ||
	  (last_round_dofs[i]>=sparsity.n_rows()))
	last_round_dofs[i] = numbers::invalid_unsigned_int;
  
    std::remove_if (last_round_dofs.begin(), last_round_dofs.end(),
		    std::bind2nd(std::equal_to<unsigned int>(),
				 numbers::invalid_unsigned_int));
  
				     // now if no valid points remain:
				     // find dof with lowest coordination
				     // number
    if (last_round_dofs.size() == 0)
      {
	unsigned int starting_point   = numbers::invalid_unsigned_int;
	unsigned int min_coordination = sparsity.n_rows();
	for (unsigned int row=0; row<sparsity.n_rows(); ++row) 
	  {
	    unsigned int j;

					     // loop until we hit the end
					     // of this row's entries
	    for (j=sparsity.get_rowstart_indices()[row];
		 j<sparsity.get_rowstart_indices()[row+1]; ++j)
	      if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
		break;
					     // post-condition after loop:
					     // coordination, i.e. the number
					     // of entries in this row is now
					     // j-rowstart[row]
	    if (j-sparsity.get_rowstart_indices()[row] <  min_coordination)
	      {
		min_coordination = j-sparsity.get_rowstart_indices()[row];
		starting_point   = row;
	      }
	  }
      
					 // now we still have to care for the
					 // case that no dof has a coordination
					 // number less than sparsity.n_rows(). this rather
					 // exotic case only happens if we only
					 // have one cell, as far as I can see,
					 // but there may be others as well.
					 //
					 // if that should be the case, we can
					 // chose an arbitrary dof as starting
					 // point, e.g. the one with number zero
	if (starting_point == numbers::invalid_unsigned_int)
	  starting_point = 0;
      
					 // initialize the first dof
	last_round_dofs.push_back (starting_point);
      }


				     // store next free dof index
    unsigned int next_free_number = 0;

				     // enumerate the first round dofs
    for (unsigned int i=0; i!=last_round_dofs.size(); ++i)
      new_indices[last_round_dofs[i]] = next_free_number++;

    bool all_dofs_renumbered = false;

				     // now do as many steps as needed to
				     // renumber all dofs
    while (!all_dofs_renumbered) 
      {
					 // store the indices of the dofs to be
					 // renumbered in the next round
	std::vector<unsigned int> next_round_dofs;

					 // find all neighbors of the
					 // dofs numbered in the last
					 // round
	for (unsigned int i=0; i<last_round_dofs.size(); ++i)
	  for (unsigned int j=sparsity.get_rowstart_indices()[last_round_dofs[i]];
	       j<sparsity.get_rowstart_indices()[last_round_dofs[i]+1]; ++j)
	    if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
	      break;
	    else
	      next_round_dofs.push_back (sparsity.get_column_numbers()[j]);
      
					 // sort dof numbers
	std::sort (next_round_dofs.begin(), next_round_dofs.end());

					 // delete multiple entries
	std::vector<unsigned int>::iterator end_sorted;
	end_sorted = std::unique (next_round_dofs.begin(), next_round_dofs.end());
	next_round_dofs.erase (end_sorted, next_round_dofs.end());

					 // eliminate dofs which are
					 // already numbered
	for (int s=next_round_dofs.size()-1; s>=0; --s)
	  if (new_indices[next_round_dofs[s]] != numbers::invalid_unsigned_int)
	    next_round_dofs.erase (next_round_dofs.begin() + s);

					 // check whether there are any new
					 // dofs
	all_dofs_renumbered = (next_round_dofs.size() == 0);
	if (all_dofs_renumbered)
					   // end loop if possible
	  continue;


					 // store for each coordination
					 // number the dofs with these
					 // coordination number
	std::multimap<unsigned int, int> dofs_by_coordination;
      
					 // find coordination number for
					 // each of these dofs
	for (std::vector<unsigned int>::iterator s=next_round_dofs.begin();
	     s!=next_round_dofs.end(); ++s) 
	  {
	    unsigned int coordination = 0;
	    for (unsigned int j=sparsity.get_rowstart_indices()[*s];
		 j<sparsity.get_rowstart_indices()[*s+1]; ++j)
	      if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
		break;
	      else
		++coordination;

					     // insert this dof at its
					     // coordination number
	    const std::pair<const unsigned int, int> new_entry (coordination, *s);
	    dofs_by_coordination.insert (new_entry);
	  }
      
					 // assign new DoF numbers to
					 // the elements of the present
					 // front:
	std::multimap<unsigned int, int>::iterator i;
	for (i = dofs_by_coordination.begin(); i!=dofs_by_coordination.end(); ++i) 
	  new_indices[i->second] = next_free_number++;

					 // after that: copy this round's
					 // dofs for the next round
	last_round_dofs = next_round_dofs;
      }

				     // test for all indices
				     // numbered. this mostly tests
				     // whether the
				     // front-marching-algorithm (which
				     // Cuthill-McKee actually is) has
				     // reached all points. it should
				     // usually do so, but might not for
				     // two reasons:
				     //
				     // - The algorithm above has a bug, or
				     // - The domain is not connected and
				     // consists of separate parts.
				     //
				     // In any case, if not all DoFs
				     // have been reached, renumbering
				     // will not be possible
    Assert ((std::find (new_indices.begin(), new_indices.end(), numbers::invalid_unsigned_int)
	     ==
	     new_indices.end())
	    &&
	    (next_free_number == sparsity.n_rows()),
	    ExcMessage("Some rows have not been renumbered. Maybe "
		       "the connectivity graph has two or more disconnected "
		       "subgraphs?"));
  }
}

DEAL_II_NAMESPACE_CLOSE
