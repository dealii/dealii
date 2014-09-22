// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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

#ifndef __deal2__mg_sparse_matrix_collection_h
#define __deal2__mg_sparse_matrix_collection_h

#include <deal.II/lac/vector.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_tools.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/std_cxx11/shared_ptr.h>

DEAL_II_NAMESPACE_OPEN

namespace mg
{
  /**
   * Handler and storage for all five SparseMatrix object involved in
   * using multigrid with local refinement.
   *
   * @author Baerbel Janssen, Guido Kanschat
   * @date 2013
   */
  template <typename number>
  class SparseMatrixCollection : public Subscriptor
  {
  public:
    void resize(const unsigned int minlevel, const unsigned  int maxlevel);

    template <class DH>
    void reinit(const DH &dof_handler);

    void set_zero();

    MGLevelObject<SparsityPattern> sparsity;
    MGLevelObject<SparsityPattern> sparsity_edge;

    MGLevelObject<SparseMatrix<number> > matrix;
    MGLevelObject<SparseMatrix<number> > matrix_down;
    MGLevelObject<SparseMatrix<number> > matrix_up;
    MGLevelObject<SparseMatrix<number> > matrix_in;
    MGLevelObject<SparseMatrix<number> > matrix_out;
  };


  template <typename number>
  void
  SparseMatrixCollection<number>::resize(const unsigned int minlevel, const unsigned  int maxlevel)
  {
    matrix.resize(minlevel, maxlevel);
    matrix.clear();
    matrix_up.resize(minlevel+1, maxlevel);
    matrix_up.clear();
    matrix_down.resize(minlevel+1, maxlevel);
    matrix_down.clear();
    matrix_in.resize(minlevel, maxlevel);
    matrix_in.clear();
    matrix_out.resize(minlevel, maxlevel);
    matrix_out.clear();
    sparsity.resize(minlevel, maxlevel);
    sparsity_edge.resize(minlevel, maxlevel);
  }


  template <typename number>
  template <class DH>
  void
  SparseMatrixCollection<number>::reinit(const DH &dof_handler)
  {
    AssertIndexRange(sparsity.max_level(), dof_handler.get_tria().n_levels());

    for (unsigned int level=sparsity.min_level();
         level<=sparsity.max_level(); ++level)
      {
        CompressedSparsityPattern c_sparsity(dof_handler.n_dofs(level));
        MGTools::make_flux_sparsity_pattern(dof_handler, c_sparsity, level);
        sparsity[level].copy_from(c_sparsity);
        matrix[level].reinit(sparsity[level]);
        matrix_in[level].reinit(sparsity[level]);
        matrix_out[level].reinit(sparsity[level]);
        if (level>0)
          {
            CompressedSparsityPattern ci_sparsity;
            ci_sparsity.reinit(dof_handler.n_dofs(level-1), dof_handler.n_dofs(level));
            MGTools::make_flux_sparsity_pattern_edge(dof_handler, ci_sparsity, level);
            sparsity_edge[level].copy_from(ci_sparsity);
            matrix_up[level].reinit(sparsity_edge[level]);
            matrix_down[level].reinit(sparsity_edge[level]);
          }
      }
  }

  template <typename number>
  void
  SparseMatrixCollection<number>::set_zero()
  {
    matrix = 0.;
    matrix_in = 0.;
    matrix_out = 0.;
    matrix_up = 0.;
    matrix_down = 0.;
  }

}

DEAL_II_NAMESPACE_CLOSE

#endif
