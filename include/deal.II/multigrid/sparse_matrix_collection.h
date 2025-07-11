// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mg_sparse_matrix_collection_h
#define dealii_mg_sparse_matrix_collection_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_tools.h>

#include <memory>

DEAL_II_NAMESPACE_OPEN

namespace mg
{
  /**
   * Handler and storage for all five SparseMatrix object involved in using
   * multigrid with local refinement.
   */
  template <typename number>
  class SparseMatrixCollection : public EnableObserverPointer
  {
  public:
    void
    resize(const unsigned int minlevel, const unsigned int maxlevel);

    template <int dim, int spacedim>
    void
    reinit(const DoFHandler<dim, spacedim> &dof_handler);

    void
    set_zero();

    MGLevelObject<SparsityPattern> sparsity;
    MGLevelObject<SparsityPattern> sparsity_edge;

    MGLevelObject<SparseMatrix<number>> matrix;
    MGLevelObject<SparseMatrix<number>> matrix_down;
    MGLevelObject<SparseMatrix<number>> matrix_up;
    MGLevelObject<SparseMatrix<number>> matrix_in;
    MGLevelObject<SparseMatrix<number>> matrix_out;
  };


  template <typename number>
  void
  SparseMatrixCollection<number>::resize(const unsigned int minlevel,
                                         const unsigned int maxlevel)
  {
    matrix.resize(minlevel, maxlevel);
    matrix.clear_elements();
    matrix_up.resize(minlevel + 1, maxlevel);
    matrix_up.clear_elements();
    matrix_down.resize(minlevel + 1, maxlevel);
    matrix_down.clear_elements();
    matrix_in.resize(minlevel, maxlevel);
    matrix_in.clear_elements();
    matrix_out.resize(minlevel, maxlevel);
    matrix_out.clear_elements();
    sparsity.resize(minlevel, maxlevel);
    sparsity_edge.resize(minlevel, maxlevel);
  }


  template <typename number>
  template <int dim, int spacedim>
  void
  SparseMatrixCollection<number>::reinit(
    const DoFHandler<dim, spacedim> &dof_handler)
  {
    AssertIndexRange(sparsity.max_level(),
                     dof_handler.get_triangulation().n_levels());

    for (unsigned int level = sparsity.min_level();
         level <= sparsity.max_level();
         ++level)
      {
        DynamicSparsityPattern dsp(dof_handler.n_dofs(level));
        MGTools::make_flux_sparsity_pattern(dof_handler, dsp, level);
        sparsity[level].copy_from(dsp);
        matrix[level].reinit(sparsity[level]);
        matrix_in[level].reinit(sparsity[level]);
        matrix_out[level].reinit(sparsity[level]);
        if (level > 0)
          {
            DynamicSparsityPattern ci_sparsity;
            ci_sparsity.reinit(dof_handler.n_dofs(level - 1),
                               dof_handler.n_dofs(level));
            MGTools::make_flux_sparsity_pattern_edge(dof_handler,
                                                     ci_sparsity,
                                                     level);
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
    matrix      = 0.;
    matrix_in   = 0.;
    matrix_out  = 0.;
    matrix_up   = 0.;
    matrix_down = 0.;
  }

} // namespace mg

DEAL_II_NAMESPACE_CLOSE

#endif
