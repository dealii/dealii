//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <numerics/mesh_worker_assembler.h>

DEAL_II_NAMESPACE_OPEN


namespace MeshWorker
{
  namespace Assembler
  {
//----------------------------------------------------------------------//
    
    template <typename number>
    inline void
    Functional<number>::initialize(unsigned int n)
    {
      results.resize(n);
    }
    
    
    template <typename number>
    template<int dim>
    inline void
    Functional<number>::assemble(const DoFInfo<dim>& info)
    {
      for (unsigned int i=0;i<results.size();++i)
	results[i] += info.J[i];
    }
    

    template <typename number>
    template<int dim>
    inline void
    Functional<number>::assemble(const DoFInfo<dim>& info1,
				 const DoFInfo<dim>& info2)
    {
      for (unsigned int i=0;i<results.size();++i)
	{
	  results[i] += info1.J[i];
	  results[i] += info2.J[i];
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
      Assert(r.name(0) == "cells", DataVectors::ExcNameMismatch(0, "cells"));
      if (sep)
	Assert(r.name(1) == "faces", DataVectors::ExcNameMismatch(1, "faces"));
      AssertDimension(r(0).n_blocks(), r(1).n_blocks());
      
      results = r;
      separate_faces = sep;
    }
    
    
    template <typename number>
    template<int dim>
    inline void
    CellsAndFaces<number>::assemble(const DoFInfo<dim>& info)
    {
      for (unsigned int i=0;i<info.J.size();++i)
	{
	  if (separate_faces &&
	      info.face_number != deal_II_numbers::invalid_unsigned_int)
	    results.vector(1).block(i)(info.face->user_index()) += info.J[i];
	  else
	    results.vector(0).block(i)(info.cell->user_index()) += info.J[i];
	}
    }
    

    template <typename number>
    template<int dim>
    inline void
    CellsAndFaces<number>::assemble(const DoFInfo<dim>& info1,
				    const DoFInfo<dim>& info2)
    {
      for (unsigned int i=0;i<info1.J.size();++i)
	{
	  if (separate_faces)
	    {
	      const double J = info1.J[i] + info2.J[i];
	      results.vector(1).block(i)(info1.face->user_index()) += J;
	      if (info2.face != info1.face)
		results.vector(1).block(i)(info2.face->user_index()) += J;
	    }
	  else
	    {
	      results.vector(0).block(i)(info1.cell->user_index()) += info1.J[i];
	      results.vector(0).block(i)(info2.cell->user_index()) += info2.J[i];
	    }
	}
    }
    

//----------------------------------------------------------------------//



  }
}


DEAL_II_NAMESPACE_CLOSE
