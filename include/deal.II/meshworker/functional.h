// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#ifndef __deal2__mesh_worker_functional_h
#define __deal2__mesh_worker_functional_h

#include <deal.II/base/named_data.h>
#include <deal.II/algorithms/any_data.h>
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
      void initialize(const unsigned int n);
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
      void initialize_info(DOFINFO &info, bool face);

      /**
       * Assemble the local values
       * into the global vectors.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble both local values
       * into the global vectors.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);

      /**
       * The value of the ith entry
       * in #results.
       */
      number operator() (const unsigned int i) const;
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
       * The initialization function, specifying the @p results
       * vectors and whether face data should be collected separately.
       *
       * @p results should contain two block vectors named "cells" and
       * "faces" (the latter only if @p separate_faces is true). In
       * each of the two, each block should have equal size and be
       * large enough to accommodate all user indices set in the cells
       * and faces covered by the loop it is used in. Typically, for
       * estimators, this is Triangulation::n_active_cells() and
       * Triangulation::n_faces(), respectively.
       *
       * The use of BlockVector may seem cumbersome, but it allows us
       * to assemble several functionals at the same time, one in each
       * block. The typical situation for error estimate is just
       * having a single block in each vector.
       */
      void initialize(AnyData &results, bool separate_faces = true);

      /**
       * @deprecated
       */
      void initialize(NamedData<BlockVector<number>*> &results,
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
      void initialize_info(DOFINFO &info, bool face) const;

      /**
       * Assemble the local values
       * into the global vectors.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble both local values
       * into the global vectors.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);

      /**
       * The value of the ith entry
       * in @p results.
       */
      number operator() (const unsigned int i) const;
    private:
      AnyData results;
      bool separate_faces;
    };
//----------------------------------------------------------------------//

    template <typename number>
    inline void
    Functional<number>::initialize(const unsigned int n)
    {
      results.resize(n);
      std::fill(results.begin(), results.end(), 0.);
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    Functional<number>::initialize_info(DOFINFO &info, bool)
    {
      info.initialize_numbers(results.size());
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    Functional<number>::assemble(const DOFINFO &info)
    {
      for (unsigned int i=0; i<results.size(); ++i)
        results[i] += info.value(i);
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    Functional<number>::assemble(const DOFINFO &info1,
                                 const DOFINFO &info2)
    {
      for (unsigned int i=0; i<results.size(); ++i)
        {
          results[i] += info1.value(i);
          results[i] += info2.value(i);
        }
    }


    template <typename number>
    inline number
    Functional<number>::operator() (const unsigned int i) const
    {
      AssertIndexRange(i, results.size());
      return results[i];
    }

//----------------------------------------------------------------------//

    template <typename number>
    inline void
    CellsAndFaces<number>::initialize(AnyData &r, bool sep)
    {
      Assert(r.name(0) == "cells", AnyData::ExcNameMismatch(0, "cells"));
      if (sep)
        {
          Assert(r.name(1) == "faces", AnyData::ExcNameMismatch(1, "faces"));
          AssertDimension(r.entry<BlockVector<double>*>(0)->n_blocks(),
                          r.entry<BlockVector<double>*>(1)->n_blocks());
        }

      results = r;
      separate_faces = sep;
    }

    template <typename number>
    inline void
    CellsAndFaces<number>::initialize(NamedData<BlockVector<number>*> &r, bool sep)
    {
      Assert(r.name(0) == "cells", AnyData::ExcNameMismatch(0, "cells"));
      if (sep)
        {
          Assert(r.name(1) == "faces", AnyData::ExcNameMismatch(1, "faces"));
          AssertDimension(r(0)->n_blocks(),r(1)->n_blocks());
        }

      results = r;
      separate_faces = sep;
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    CellsAndFaces<number>::initialize_info(DOFINFO &info, bool) const
    {
      info.initialize_numbers(results.entry<BlockVector<double>*>(0)->n_blocks());
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    CellsAndFaces<number>::assemble(const DOFINFO &info)
    {
      BlockVector<double> *v;
      if (separate_faces &&
          info.face_number != deal_II_numbers::invalid_unsigned_int)
        v = results.entry<BlockVector<double>*>(1);
      else
        v = results.entry<BlockVector<double>*>(0);

      for (unsigned int i=0; i<info.n_values(); ++i)
        v->block(i)(info.cell->user_index()) += info.value(i);
    }


    template <typename number>
    template <class DOFINFO>
    inline void
    CellsAndFaces<number>::assemble(const DOFINFO &info1,
                                    const DOFINFO &info2)
    {
      for (unsigned int i=0; i<info1.n_values(); ++i)
        {
          if (separate_faces)
            {
              BlockVector<double> *v1 = results.entry<BlockVector<double>*>(1);
              const double J = info1.value(i) + info2.value(i);
              v1->block(i)(info1.face->user_index()) += J;
              if (info2.face != info1.face)
                v1->block(i)(info2.face->user_index()) += J;
            }
          else
            {
              BlockVector<double> *v0 = results.entry<BlockVector<double>*>(0);
              v0->block(i)(info1.cell->user_index()) += .5*info1.value(i);
              v0->block(i)(info2.cell->user_index()) += .5*info2.value(i);
            }
        }
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
