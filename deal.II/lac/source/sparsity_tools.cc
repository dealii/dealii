//---------------------------------------------------------------------------
//    $Id: sparsity_pattern.cc 16733 2008-09-03 14:22:26Z heister $
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
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

#ifdef DEAL_II_USE_METIS
extern "C" {
#  include <metis.h>
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
				       sparsity_pattern.n_cols()+1);
    std::vector<idxtype> int_colnums (sparsity_pattern.get_colnums(),
				      sparsity_pattern.get_colnums()+max_vec_len+1);

    std::vector<idxtype> int_partition_indices (sparsity_pattern.n_rows());

				     // Select which type of partitioning to create

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

}

DEAL_II_NAMESPACE_CLOSE
