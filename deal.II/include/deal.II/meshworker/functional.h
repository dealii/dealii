//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mesh_worker_functional_h
#define __deal2__mesh_worker_functional_h

#include <deal.II/base/named_data.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>


DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  namespace Assembler
  {
/**
 * The class assembling local contributions to a functional into
 * global functionals.
 *
 * 
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
    template <typename number = double>
    class Functional
    {
      public:
					 /**
					  * Initialize local data to
					  * store functionals. The
					  * number <tt>n</tt> is the
					  * number of functionals to
					  * be computed.
					  */
	void initialize(unsigned int n);
					 /**
					  * Initialize the local data
					  * in the DoFInfo object used
					  * later for assembling.
					  *
					  * The info object refers to
					  * a cell if
					  * <code>!face</code>, or
					  * else to an interior or
					  * boundary face.
					  */
	template <class DOFINFO>
	void initialize_info(DOFINFO& info, bool face);
	
					 /**
					  * Assemble the local values
					  * into the global vectors.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info);
	
					 /**
					  * Assemble both local values
					  * into the global vectors.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info1,
		      const DOFINFO& info2);

					 /**
					  * The value of the ith entry
					  * in #results.
					  */
	number operator() (unsigned int i) const;
      private:
					 /**
					  * The values into which the
					  * results are added.
					  */
	std::vector<double> results;
    };

/**
 * Compute cell and face contributions of one or several functionals,
 * typically for error estimates.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
    template <typename number = double>
    class CellsAndFaces
    {
      public:
					 /**
					  * The data type for
					  * communicating the cell and
					  * face vectors.
					  */
	typedef NamedData<BlockVector<number>*> DataVectors;
	
					 /**
					  * The initialization
					  * function, specifying the
					  * #results vectors and
					  * whether face data should
					  * be collected separately.
					  *
					  * #results should contain
					  * two block vectors named
					  * "cells" and "faces" (the
					  * latter only if
					  * #separate_faces is
					  * true). In each of the two,
					  * each block should have
					  * equal size and be large
					  * enough to accomodate all
					  * user indices set in the
					  * cells and faces covered by
					  * the loop it is used
					  * in. Typically, for
					  * estimators, this is
					  * Triangulation::n_active_cells()
					  * and
					  * Triangulation::n_faces(),
					  * respectively.
					  *
					  * The use of BlockVector may
					  * seem cumbersome, but it
					  * allows us to assemble
					  * several functionals at the
					  * same time, one in each
					  * block. The typical
					  * situation for error
					  * estimate is just having a
					  * single block in each vector.
					  */
	void initialize(DataVectors& results,
			bool separate_faces = true);
					 /**
					  * Initialize the local data
					  * in the
					  * DoFInfo
					  * object used later for
					  * assembling.
					  *
					  * The info object refers to
					  * a cell if
					  * <code>!face</code>, or
					  * else to an interior or
					  * boundary face.
					  */
	template <class DOFINFO>
	void initialize_info(DOFINFO& info, bool face) const;
	
					 /**
					  * Assemble the local values
					  * into the global vectors.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info);
	
					 /**
					  * Assemble both local values
					  * into the global vectors.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info1,
		      const DOFINFO& info2);

					 /**
					  * The value of the ith entry
					  * in #results.
					  */
	number operator() (unsigned int i) const;
      private:
	DataVectors results;
	bool separate_faces;
    };
//----------------------------------------------------------------------//
    
    template <typename number>
    inline void
    Functional<number>::initialize(unsigned int n)
    {
      results.resize(n);
      std::fill(results.begin(), results.end(), 0.);
    }
    
    
    template <typename number>
    template <class DOFINFO>
    inline void
    Functional<number>::initialize_info(DOFINFO& info, bool)
    {
      info.initialize_numbers(results.size());
    }
    
    
    template <typename number>
    template <class DOFINFO>
    inline void
    Functional<number>::assemble(const DOFINFO& info)
    {
      for (unsigned int i=0;i<results.size();++i)
	results[i] += info.value(i);
    }
    

    template <typename number>
    template <class DOFINFO>
    inline void
    Functional<number>::assemble(const DOFINFO& info1,
				 const DOFINFO& info2)
    {
      for (unsigned int i=0;i<results.size();++i)
	{
	  results[i] += info1.value(i);
	  results[i] += info2.value(i);
	}
    }

    
    template <typename number>
    inline number
    Functional<number>::operator() (unsigned int i) const
    {
      AssertIndexRange(i, results.size());
      return results[i];
    }
    
//----------------------------------------------------------------------//

    template <typename number>
    inline void
    CellsAndFaces<number>::initialize(DataVectors& r, bool sep)
    {
      Assert(r.name(0) == "cells", typename DataVectors::ExcNameMismatch(0, "cells"));
      if (sep)
	{
	  Assert(r.name(1) == "faces", typename DataVectors::ExcNameMismatch(1, "faces"));
	  AssertDimension(r(0)->n_blocks(), r(1)->n_blocks());
	}
      
      results = r;
      separate_faces = sep;
    }

    
    template <typename number>
    template <class DOFINFO>
    inline void
    CellsAndFaces<number>::initialize_info(DOFINFO& info, bool) const
    {
      info.initialize_numbers(results(0)->n_blocks());
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    CellsAndFaces<number>::assemble(const DOFINFO& info)
    {
      for (unsigned int i=0;i<info.n_values();++i)
	{
	  if (separate_faces &&
	      info.face_number != deal_II_numbers::invalid_unsigned_int)
	    results(1)->block(i)(info.face->user_index()) += info.value(i);
	  else
	    results(0)->block(i)(info.cell->user_index()) += info.value(i);
	}
    }
    

    template <typename number>
    template <class DOFINFO>
    inline void
    CellsAndFaces<number>::assemble(const DOFINFO& info1,
				    const DOFINFO& info2)
    {
      for (unsigned int i=0;i<info1.n_values();++i)
	{
	  if (separate_faces)
	    {
	      const double J = info1.value(i) + info2.value(i);
	      results(1)->block(i)(info1.face->user_index()) += J;
	      if (info2.face != info1.face)
		results(1)->block(i)(info2.face->user_index()) += J;
	    }
	  else
	    {
	      results(0)->block(i)(info1.cell->user_index()) += .5*info1.value(i);
	      results(0)->block(i)(info2.cell->user_index()) += .5*info2.value(i);
	    }
	}
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
