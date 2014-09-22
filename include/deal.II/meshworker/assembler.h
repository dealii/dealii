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


#ifndef __deal2__mesh_worker_assembler_h
#define __deal2__mesh_worker_assembler_h

#include <deal.II/base/named_data.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/functional.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>


DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  /**
   * The namespace containing objects that can be used to assemble data
   * computed on cells and faces into global objects. This can reach
   * from collecting the total error estimate from cell and face
   * contributions to assembling matrices and multilevel matrices.
   *
   * <h3>Data models</h3>
   *
   * The class chosen from this namespace determines which data model is
   * used. For the local as well as the global objects, we have the
   * choice between two models:
   *
   * <h4>The comprehensive data model</h4>
   *
   * This is the structure set up by the FESystem class. Globally, this
   * means, data is assembled into one residual vector and into one
   * matrix. These objects may be block vectors and block matrices, but
   * the process of assembling ignores this fact.
   *
   * Similarly, there is only a single cell vector and cell matrix,
   * respectively, which is indexed by all degrees of freedom of the
   * FESystem. When building the cell matrix, it is necessary to
   * distinguish between the different components of the system and
   * select the right operator for each pair.
   *
   * <h4>The blocked data model</h4>
   *
   * Here, all the blocks are treated separately (in spite of using
   * FESystem for its convenience in other places). For instance, no
   * block matrix is assembled, but a list of blocks, which can be
   * combined later by BlockMatrixArray. Locally, this means, that each
   * matrix block of a system is generated separately and assembled into
   * the corresponding global block.
   *
   * This approach is advantageous, if the number of matrices for each
   * block position in the global system is different. For instance,
   * block preconditioners for the Oseen problem require 3 pressure
   * matrices, but only one divergence and one advection-diffusion
   * operator for velocities.
   *
   * Additionally, this approach enables the construction of a system of
   * equations from building blocks for each equation and coupling
   * operator.
   *
   * Nevertheless, since a separate FEValues object must be created for
   * each base element, it is not quite clear a priori, which data model
   * is more efficient.
   *
   * @ingroup MeshWorker
   * @author Guido Kanschat, 2009
   */
  namespace Assembler
  {
    /**
     * Assemble local residuals into global residuals.
     *
     * The global residuals are expected as an FEVectors object.
     * The local residuals are block vectors.
     *
     * Depending on whether the BlockInfo object was initialize with
     * BlockInfo::initialize_local(), the comprehensive or block data
     * model is used locally.
     *
     * In the block model, each of the blocks of the local vectors
     * corresponds to the restriction of a single block of the system to
     * this cell (@ref GlossBlock). Thus, the size of this local block is
     * the number of degrees of freedom of the corresponding base element
     * of the FESystem.
     *
     * @todo Comprehensive model currently not implemented.
     *
     * @ingroup MeshWorker
     * @author Guido Kanschat, 2009
     */
    template <class VECTOR>
    class ResidualLocalBlocksToGlobalBlocks
    {
    public:
      /**
       * Copy the BlockInfo and the
       * matrix pointers into local
       * variables.
       */
      void initialize(const BlockInfo *block_info,
                      NamedData<VECTOR *> &residuals);
      /**
       * Initialize the constraints.
       */
      void initialize(const ConstraintMatrix &constraints);
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
       * Assemble the local residuals
       * into the global residuals.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble both local residuals
       * into the global residuals.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);
    private:
      /**
       * Assemble a single local
       * residual into the global.
       */
      void assemble(VECTOR &global,
                    const BlockVector<double> &local,
                    const std::vector<types::global_dof_index> &dof);

      /**
       * The global matrices,
       * stored as a vector of
       * pointers.
       */
      NamedData<SmartPointer<VECTOR,
                ResidualLocalBlocksToGlobalBlocks<VECTOR> > > residuals;

      /**
       * A pointer to the object containing the block structure.
       */
      SmartPointer<const BlockInfo,
                   ResidualLocalBlocksToGlobalBlocks<VECTOR> > block_info;
      /**
       * A pointer to the object containing constraints.
       */
      SmartPointer<const ConstraintMatrix,
                   ResidualLocalBlocksToGlobalBlocks<VECTOR> > constraints;
    };


    /**
     * A helper class assembling local matrices into global matrices.
     *
     * The global matrices are expected as a vector of MatrixBlock
     * objects, each containing a matrix object with a function
     * corresponding to SparseMatrix::add() and information on the block
     * row and column this matrix represents in a block system.
     *
     * The local matrices are expected as a similar vector of MatrixBlock
     * objects, but containing a FullMatrix.
     *
     * Like with ResidualLocalBlocksToGlobalBlocks, the initialization of
     * the BlockInfo object decides whether the comprehensive data model
     * or the block model is used.
     *
     * In the comprehensive model, each of the LocalMatrixBlocks has
     * coordinates (0,0) and dimensions equal to the number of degrees of
     * freedom of the FESystem.
     *
     * In the comprehensive model, each block has its own block
     * coordinates and the size depends on the associated
     * FESystem::base_element(). These blocks can be generated separately
     * and will be assembled into the correct matrix block by this object.
     *
     * @ingroup MeshWorker
     * @author Guido Kanschat, 2009
     */
    template <class MATRIX, typename number = double>
    class MatrixLocalBlocksToGlobalBlocks
    {
    public:
      /**
       * Constructor, initializing
       * the #threshold, which
       * limits how small numbers
       * may be to be entered into
       * the matrix.
       */
      MatrixLocalBlocksToGlobalBlocks(double threshold = 1.e-12);

      /**
       * Copy the BlockInfo and the
       * matrix pointers into local
       * variables and initialize
       * cell matrix vectors.
       */
      void initialize(const BlockInfo *block_info,
                      MatrixBlockVector<MATRIX> &matrices);

      /**
       * Initialize the constraints.
       */
      void initialize(const ConstraintMatrix &constraints);
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
       * Assemble the local matrices
       * into the global matrices.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble all local matrices
       * into the global matrices.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);

    private:
      /**
       * Assemble a single local
       * matrix into a global one.
       */
      void assemble(
        MatrixBlock<MATRIX> &global,
        const FullMatrix<number> &local,
        const unsigned int block_row,
        const unsigned int block_col,
        const std::vector<types::global_dof_index> &dof1,
        const std::vector<types::global_dof_index> &dof2);

      /**
       * The global matrices,
       * stored as a vector of
       * pointers.
       */
      SmartPointer<MatrixBlockVector<MATRIX>,
                   MatrixLocalBlocksToGlobalBlocks<MATRIX, number> > matrices;

      /**
       * A pointer to the object containing the block structure.
       */
      SmartPointer<const BlockInfo,
                   MatrixLocalBlocksToGlobalBlocks<MATRIX, number> > block_info;
      /**
       * A pointer to the object containing constraints.
       */
      SmartPointer<const ConstraintMatrix,
                   MatrixLocalBlocksToGlobalBlocks<MATRIX,number> > constraints;

      /**
       * The smallest positive
       * number that will be
       * entered into the global
       * matrix. All smaller
       * absolute values will be
       * treated as zero and will
       * not be assembled.
       */
      const double threshold;

    };

    /**
     * A helper class assembling local matrices into global multilevel
     * matrices. This class is the multilevel equivalent of
     * MatrixLocalBlocksToGlobalBlocks and documentation of that class
     * applies here to a large extend.
     *
     * The global matrices are expected as a vector of pointers to MatrixBlock
     * objects, each containing a MGLevelObject with matrices with a function
     * corresponding to SparseMatrix::add() and information on the block
     * row and column this matrix represents in a block system.
     *
     * The local matrices are a similar vector of MatrixBlock objects, but
     * containing a FullMatrix.
     *
     * If local refinement occurs, the Multigrid method needs more
     * matrices, two for continuous elements and another two if numerical
     * fluxes are computed on interfaces. The second set can be added
     * using initialize_edge_flux(). Once added, the contributions in all
     * participating matrices will be assembled from the cell and face
     * matrices automatically.
     *
     * @ingroup MeshWorker
     * @author Guido Kanschat, 2009
     */
    template <class MATRIX, typename number = double>
    class MGMatrixLocalBlocksToGlobalBlocks
    {
    public:
      typedef MGMatrixBlockVector<MATRIX> MatrixPtrVector;
      typedef SmartPointer<MatrixPtrVector, MGMatrixLocalBlocksToGlobalBlocks<MATRIX,number> >
      MatrixPtrVectorPtr;

      /**
       * Constructor, initializing
       * the #threshold, which
       * limits how small numbers
       * may be to be entered into
       * the matrix.
       */
      MGMatrixLocalBlocksToGlobalBlocks(double threshold = 1.e-12);

      /**
       * Copy the BlockInfo and the
       * matrix pointers into local
       * variables and initialize
       * cell matrix vectors.
       */
      void initialize(const BlockInfo *block_info,
                      MatrixPtrVector &matrices);

      /**
       * Initialize the multilevel
       * constraints.
       */
      void initialize(const MGConstrainedDoFs &mg_constrained_dofs);

      /**
       * Multigrid methods on
       * locally refined meshes
       * need additional
       * matrices. For
       * discontinuous Galerkin
       * methods, these are two
       * flux matrices across the
       * refinement edge, which are
       * set by this method.
       */
      void initialize_edge_flux(MatrixPtrVector &up, MatrixPtrVector &down);

      /**
       * Multigrid methods on
       * locally refined meshes
       * need additional
       * matrices. For
       * discontinuous Galerkin
       * methods, these are two
       * flux matrices across the
       * refinement edge, which are
       * set by this method.
       */
      void initialize_interfaces (MatrixPtrVector &interface_in, MatrixPtrVector &interface_out);
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
       * Assemble the local matrices
       * into the global matrices.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble all local matrices
       * into the global matrices.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);

    private:
      /**
       * Assemble a single local
       * matrix into a global one.
       */
      void assemble(
        MATRIX &global,
        const FullMatrix<number> &local,
        const unsigned int block_row,
        const unsigned int block_col,
        const std::vector<types::global_dof_index> &dof1,
        const std::vector<types::global_dof_index> &dof2,
        const unsigned int level1,
        const unsigned int level2,
        bool transpose = false);

      /**
       * Assemble a single local
       * matrix into a global one.
       */
      void assemble_fluxes(
        MATRIX &global,
        const FullMatrix<number> &local,
        const unsigned int block_row,
        const unsigned int block_col,
        const std::vector<types::global_dof_index> &dof1,
        const std::vector<types::global_dof_index> &dof2,
        const unsigned int level1,
        const unsigned int level2);

      /**
       * Assemble a single local
       * matrix into a global one.
       */
      void assemble_up(
        MATRIX &global,
        const FullMatrix<number> &local,
        const unsigned int block_row,
        const unsigned int block_col,
        const std::vector<types::global_dof_index> &dof1,
        const std::vector<types::global_dof_index> &dof2,
        const unsigned int level1,
        const unsigned int level2);

      /**
       * Assemble a single local
       * matrix into a global one.
       */
      void assemble_down(
        MATRIX &global,
        const FullMatrix<number> &local,
        const unsigned int block_row,
        const unsigned int block_col,
        const std::vector<types::global_dof_index> &dof1,
        const std::vector<types::global_dof_index> &dof2,
        const unsigned int level1,
        const unsigned int level2);

      /**
       * Assemble a single local
       * matrix into a global one.
       */
      void assemble_in(
        MATRIX &global,
        const FullMatrix<number> &local,
        const unsigned int block_row,
        const unsigned int block_col,
        const std::vector<types::global_dof_index> &dof1,
        const std::vector<types::global_dof_index> &dof2,
        const unsigned int level1,
        const unsigned int level2);

      /**
       * Assemble a single local
       * matrix into a global one.
       */
      void assemble_out(
        MATRIX &global,
        const FullMatrix<number> &local,
        const unsigned int block_row,
        const unsigned int block_col,
        const std::vector<types::global_dof_index> &dof1,
        const std::vector<types::global_dof_index> &dof2,
        const unsigned int level1,
        const unsigned int level2);

      /**
       * The level matrices,
       * stored as a vector of
       * pointers.
       */
      MatrixPtrVectorPtr matrices;

      /**
       * The flux matrix between
       * the fine and the coarse
       * level at refinement edges.
       */
      MatrixPtrVectorPtr flux_down;

      /**
       * The flux matrix between
       * the coarse and the fine
       * level at refinement edges.
       */
      MatrixPtrVectorPtr flux_up;

      /**
       * The interface matrix between
       * the fine and the coarse
       * level at refinement edges.
       */
      MatrixPtrVectorPtr interface_out;

      /**
       * The interface matrix between
       * the coarse and the fine
       * level at refinement edges.
       */
      MatrixPtrVectorPtr interface_in;

      /**
       * A pointer to the object containing the block structure.
       */
      SmartPointer<const BlockInfo, MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number> > block_info;

      /**
       * A pointer to the object containing constraints.
       */
      SmartPointer<const MGConstrainedDoFs,MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number> > mg_constrained_dofs;


      /**
       * The smallest positive
       * number that will be
       * entered into the global
       * matrix. All smaller
       * absolute values will be
       * treated as zero and will
       * not be assembled.
       */
      const double threshold;

    };

//----------------------------------------------------------------------//

    template <class VECTOR>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VECTOR>::initialize(const BlockInfo *b,
                                                          NamedData<VECTOR *> &m)
    {
      block_info = b;
      residuals = m;
    }

    template <class VECTOR>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VECTOR>::initialize(
      const ConstraintMatrix &c)
    {
      constraints = &c;
    }


    template <class VECTOR>
    template <class DOFINFO>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VECTOR>::initialize_info(
      DOFINFO &info, bool) const
    {
      info.initialize_vectors(residuals.size());
    }

    template <class VECTOR>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VECTOR>::assemble(
      VECTOR &global,
      const BlockVector<double> &local,
      const std::vector<types::global_dof_index> &dof)
    {
      if (constraints == 0)
        {
          for (unsigned int b=0; b<local.n_blocks(); ++b)
            for (unsigned int j=0; j<local.block(b).size(); ++j)
              {
                // The coordinates of
                // the current entry in
                // DoFHandler
                // numbering, which
                // differs from the
                // block-wise local
                // numbering we use in
                // our local vectors
                const unsigned int jcell = this->block_info->local().local_to_global(b, j);
                global(dof[jcell]) += local.block(b)(j);
              }
        }
      else
        constraints->distribute_local_to_global(local, dof, global);
    }


    template <class VECTOR>
    template <class DOFINFO>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VECTOR>::assemble(
      const DOFINFO &info)
    {
      for (unsigned int i=0; i<residuals.size(); ++i)
        assemble(*residuals(i), info.vector(i), info.indices);
    }


    template <class VECTOR>
    template <class DOFINFO>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VECTOR>::assemble(
      const DOFINFO &info1,
      const DOFINFO &info2)
    {
      for (unsigned int i=0; i<residuals.size(); ++i)
        {
          assemble(*residuals(i), info1.vector(i), info1.indices);
          assemble(*residuals(i), info2.vector(i), info2.indices);
        }
    }


//----------------------------------------------------------------------//

    template <class MATRIX, typename number>
    inline
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::MatrixLocalBlocksToGlobalBlocks(
      double threshold)
      :
      threshold(threshold)
    {}


    template <class MATRIX, typename number>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize(
      const BlockInfo *b,
      MatrixBlockVector<MATRIX> &m)
    {
      block_info = b;
      matrices = &m;
    }



    template <class MATRIX, typename number>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize(
      const ConstraintMatrix &c)
    {
      constraints = &c;
    }



    template <class MATRIX ,typename number>
    template <class DOFINFO>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize_info(
      DOFINFO &info,
      bool face) const
    {
      info.initialize_matrices(*matrices, face);
    }



    template <class MATRIX, typename number>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(
      MatrixBlock<MATRIX> &global,
      const FullMatrix<number> &local,
      const unsigned int block_row,
      const unsigned int block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2)
    {
      if (constraints == 0)
        {
          for (unsigned int j=0; j<local.n_rows(); ++j)
            for (unsigned int k=0; k<local.n_cols(); ++k)
              if (std::fabs(local(j,k)) >= threshold)
                {
                  // The coordinates of
                  // the current entry in
                  // DoFHandler
                  // numbering, which
                  // differs from the
                  // block-wise local
                  // numbering we use in
                  // our local matrices
                  const unsigned int jcell = this->block_info->local().local_to_global(block_row, j);
                  const unsigned int kcell = this->block_info->local().local_to_global(block_col, k);

                  global.add(dof1[jcell], dof2[kcell], local(j,k));
                }
        }
      else
        {
          const BlockIndices &bi = this->block_info->local();
          std::vector<types::global_dof_index> sliced_row_indices (bi.block_size(block_row));
          for (unsigned int i=0; i<sliced_row_indices.size(); ++i)
            sliced_row_indices[i] = dof1[bi.block_start(block_row)+i];

          std::vector<types::global_dof_index> sliced_col_indices (bi.block_size(block_col));
          for (unsigned int i=0; i<sliced_col_indices.size(); ++i)
            sliced_col_indices[i] = dof2[bi.block_start(block_col)+i];

          constraints->distribute_local_to_global(local,
                                                  sliced_row_indices, sliced_col_indices, global);
        }
    }


    template <class MATRIX, typename number>
    template <class DOFINFO>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(
      const DOFINFO &info)
    {
      for (unsigned int i=0; i<matrices->size(); ++i)
        {
          // Row and column index of
          // the block we are dealing with
          const types::global_dof_index row = matrices->block(i).row;
          const types::global_dof_index col = matrices->block(i).column;

          assemble(matrices->block(i), info.matrix(i,false).matrix, row, col, info.indices, info.indices);
        }
    }


    template <class MATRIX, typename number>
    template <class DOFINFO>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(
      const DOFINFO &info1,
      const DOFINFO &info2)
    {
      for (unsigned int i=0; i<matrices->size(); ++i)
        {
          // Row and column index of
          // the block we are dealing with
          const types::global_dof_index row = matrices->block(i).row;
          const types::global_dof_index col = matrices->block(i).column;

          assemble(matrices->block(i), info1.matrix(i,false).matrix, row, col, info1.indices, info1.indices);
          assemble(matrices->block(i), info1.matrix(i,true).matrix, row, col, info1.indices, info2.indices);
          assemble(matrices->block(i), info2.matrix(i,false).matrix, row, col, info2.indices, info2.indices);
          assemble(matrices->block(i), info2.matrix(i,true).matrix, row, col, info2.indices, info1.indices);
        }
    }


// ----------------------------------------------------------------------//

    template <class MATRIX, typename number>
    inline
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::MGMatrixLocalBlocksToGlobalBlocks(
      double threshold)
      :
      threshold(threshold)
    {}


    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize(
      const BlockInfo *b,
      MatrixPtrVector &m)
    {
      block_info = b;
      AssertDimension(block_info->local().size(), block_info->global().size());
      matrices = &m;
    }


    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize(
      const MGConstrainedDoFs &mg_c)
    {
      mg_constrained_dofs = &mg_c;
    }


    template <class MATRIX ,typename number>
    template <class DOFINFO>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize_info(
      DOFINFO &info,
      bool face) const
    {
      info.initialize_matrices(*matrices, face);
    }



    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize_edge_flux(
      MatrixPtrVector &up,
      MatrixPtrVector &down)
    {
      flux_up = up;
      flux_down = down;
    }


    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize_interfaces(
      MatrixPtrVector &in,
      MatrixPtrVector &out)
    {
      interface_in = in;
      interface_out = out;
    }


    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(
      MATRIX &global,
      const FullMatrix<number> &local,
      const unsigned int block_row,
      const unsigned int block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int level1,
      const unsigned int level2,
      bool transpose)
    {
      for (unsigned int j=0; j<local.n_rows(); ++j)
        for (unsigned int k=0; k<local.n_cols(); ++k)
          if (std::fabs(local(j,k)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell = this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell = this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1).global_to_local(dof1[jcell]).second;
              const unsigned int kglobal = this->block_info->level(level2).global_to_local(dof2[kcell]).second;

              if (mg_constrained_dofs == 0)
                {
                  if (transpose)
                    global.add(kglobal, jglobal, local(j,k));
                  else
                    global.add(jglobal, kglobal, local(j,k));
                }
              else
                {
                  if (!mg_constrained_dofs->at_refinement_edge(level1, jglobal) &&
                      !mg_constrained_dofs->at_refinement_edge(level2, kglobal))
                    {
                      if (mg_constrained_dofs->set_boundary_values())
                        {
                          if ((!mg_constrained_dofs->is_boundary_index(level1, jglobal) &&
                               !mg_constrained_dofs->is_boundary_index(level2, kglobal))
                              ||
                              (mg_constrained_dofs->is_boundary_index(level1, jglobal) &&
                               mg_constrained_dofs->is_boundary_index(level2, kglobal) &&
                               jglobal == kglobal))
                            {
                              if (transpose)
                                global.add(kglobal, jglobal, local(j,k));
                              else
                                global.add(jglobal, kglobal, local(j,k));
                            }
                        }
                      else
                        {
                          if (transpose)
                            global.add(kglobal, jglobal, local(j,k));
                          else
                            global.add(jglobal, kglobal, local(j,k));
                        }
                    }
                }
            }
    }


    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble_fluxes(
      MATRIX &global,
      const FullMatrix<number> &local,
      const unsigned int block_row,
      const unsigned int block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int level1,
      const unsigned int level2)
    {
      for (unsigned int j=0; j<local.n_rows(); ++j)
        for (unsigned int k=0; k<local.n_cols(); ++k)
          if (std::fabs(local(j,k)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell = this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell = this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1).global_to_local(dof1[jcell]).second;
              const unsigned int kglobal = this->block_info->level(level2).global_to_local(dof2[kcell]).second;

              if (mg_constrained_dofs == 0)
                global.add(jglobal, kglobal, local(j,k));
              else
                {
                  if (!mg_constrained_dofs->non_refinement_edge_index(level1, jglobal) &&
                      !mg_constrained_dofs->non_refinement_edge_index(level2, kglobal))
                    {
                      if (!mg_constrained_dofs->at_refinement_edge(level1, jglobal) &&
                          !mg_constrained_dofs->at_refinement_edge(level2, kglobal))
                        global.add(jglobal, kglobal, local(j,k));
                    }
                }
            }
    }

    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble_up(
      MATRIX &global,
      const FullMatrix<number> &local,
      const unsigned int block_row,
      const unsigned int block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int level1,
      const unsigned int level2)
    {
      for (unsigned int j=0; j<local.n_rows(); ++j)
        for (unsigned int k=0; k<local.n_cols(); ++k)
          if (std::fabs(local(j,k)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell = this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell = this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1).global_to_local(dof1[jcell]).second;
              const unsigned int kglobal = this->block_info->level(level2).global_to_local(dof2[kcell]).second;

              if (mg_constrained_dofs == 0)
                global.add(jglobal, kglobal, local(j,k));
              else
                {
                  if (!mg_constrained_dofs->non_refinement_edge_index(level1, jglobal) &&
                      !mg_constrained_dofs->non_refinement_edge_index(level2, kglobal))
                    {
                      if (!mg_constrained_dofs->at_refinement_edge(level1, jglobal) &&
                          !mg_constrained_dofs->at_refinement_edge(level2, kglobal))
                        global.add(jglobal, kglobal, local(j,k));
                    }
                }
            }
    }

    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble_down(
      MATRIX &global,
      const FullMatrix<number> &local,
      const unsigned int block_row,
      const unsigned int block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int level1,
      const unsigned int level2)
    {
      for (unsigned int j=0; j<local.n_rows(); ++j)
        for (unsigned int k=0; k<local.n_cols(); ++k)
          if (std::fabs(local(k,j)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell = this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell = this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1).global_to_local(dof1[jcell]).second;
              const unsigned int kglobal = this->block_info->level(level2).global_to_local(dof2[kcell]).second;

              if (mg_constrained_dofs == 0)
                global.add(jglobal, kglobal, local(k,j));
              else
                {
                  if (!mg_constrained_dofs->non_refinement_edge_index(level1, jglobal) &&
                      !mg_constrained_dofs->non_refinement_edge_index(level2, kglobal))
                    {
                      if (!mg_constrained_dofs->at_refinement_edge(level1, jglobal) &&
                          !mg_constrained_dofs->at_refinement_edge(level2, kglobal))
                        global.add(jglobal, kglobal, local(k,j));
                    }
                }
            }
    }

    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble_in(
      MATRIX &global,
      const FullMatrix<number> &local,
      const unsigned int block_row,
      const unsigned int block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int level1,
      const unsigned int level2)
    {
//      AssertDimension(local.n(), dof1.size());
//      AssertDimension(local.m(), dof2.size());

      for (unsigned int j=0; j<local.n_rows(); ++j)
        for (unsigned int k=0; k<local.n_cols(); ++k)
          if (std::fabs(local(j,k)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell = this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell = this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1).global_to_local(dof1[jcell]).second;
              const unsigned int kglobal = this->block_info->level(level2).global_to_local(dof2[kcell]).second;

              if (mg_constrained_dofs == 0)
                global.add(jglobal, kglobal, local(j,k));
              else
                {
                  if (mg_constrained_dofs->at_refinement_edge(level1, jglobal) &&
                      !mg_constrained_dofs->at_refinement_edge(level2, kglobal))
                    {
                      if (mg_constrained_dofs->set_boundary_values())
                        {
                          if ((!mg_constrained_dofs->is_boundary_index(level1, jglobal) &&
                               !mg_constrained_dofs->is_boundary_index(level2, kglobal))
                              ||
                              (mg_constrained_dofs->is_boundary_index(level1, jglobal) &&
                               mg_constrained_dofs->is_boundary_index(level2, kglobal) &&
                               jglobal == kglobal))
                            global.add(jglobal, kglobal, local(j,k));
                        }
                      else
                        global.add(jglobal, kglobal, local(j,k));
                    }
                }
            }
    }

    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble_out(
      MATRIX &global,
      const FullMatrix<number> &local,
      const unsigned int block_row,
      const unsigned int block_col,
      const std::vector<types::global_dof_index> &dof1,
      const std::vector<types::global_dof_index> &dof2,
      const unsigned int level1,
      const unsigned int level2)
    {
//      AssertDimension(local.n(), dof1.size());
//      AssertDimension(local.m(), dof2.size());

      for (unsigned int j=0; j<local.n_rows(); ++j)
        for (unsigned int k=0; k<local.n_cols(); ++k)
          if (std::fabs(local(k,j)) >= threshold)
            {
              // The coordinates of
              // the current entry in
              // DoFHandler
              // numbering, which
              // differs from the
              // block-wise local
              // numbering we use in
              // our local matrices
              const unsigned int jcell = this->block_info->local().local_to_global(block_row, j);
              const unsigned int kcell = this->block_info->local().local_to_global(block_col, k);

              // The global dof
              // indices to assemble
              // in. Since we may
              // have face matrices
              // coupling two
              // different cells, we
              // provide two sets of
              // dof indices.
              const unsigned int jglobal = this->block_info->level(level1).global_to_local(dof1[jcell]).second;
              const unsigned int kglobal = this->block_info->level(level2).global_to_local(dof2[kcell]).second;

              if (mg_constrained_dofs == 0)
                global.add(jglobal, kglobal, local(k,j));
              else
                {
                  if (mg_constrained_dofs->at_refinement_edge(level1, jglobal) &&
                      !mg_constrained_dofs->at_refinement_edge(level2, kglobal))
                    {
                      if (mg_constrained_dofs->set_boundary_values())
                        {
                          if ((!mg_constrained_dofs->is_boundary_index(level1, jglobal) &&
                               !mg_constrained_dofs->is_boundary_index(level2, kglobal))
                              ||
                              (mg_constrained_dofs->is_boundary_index(level1, jglobal) &&
                               mg_constrained_dofs->is_boundary_index(level2, kglobal) &&
                               jglobal == kglobal))
                            global.add(jglobal, kglobal, local(k,j));
                        }
                      else
                        global.add(jglobal, kglobal, local(k,j));
                    }
                }
            }
    }


    template <class MATRIX, typename number>
    template <class DOFINFO>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(const DOFINFO &info)
    {
      const unsigned int level = info.cell->level();

      for (unsigned int i=0; i<matrices->size(); ++i)
        {
          // Row and column index of
          // the block we are dealing with
          const unsigned int row = matrices->block(i)[level].row;
          const unsigned int col = matrices->block(i)[level].column;

          assemble(matrices->block(i)[level].matrix, info.matrix(i,false).matrix, row, col,
                   info.indices, info.indices, level, level);
          if (mg_constrained_dofs != 0)
            {
              if (interface_in != 0)
                assemble_in(interface_in->block(i)[level], info.matrix(i,false).matrix, row, col,
                            info.indices, info.indices, level, level);
              if (interface_out != 0)
                assemble_in(interface_out->block(i)[level], info.matrix(i,false).matrix, row, col,
                            info.indices, info.indices, level, level);

              assemble_in(matrices->block_in(i)[level], info.matrix(i,false).matrix, row, col,
                          info.indices, info.indices, level, level);
              assemble_out(matrices->block_out(i)[level], info.matrix(i,false).matrix, row, col,
                           info.indices, info.indices, level, level);
            }
        }
    }


    template <class MATRIX, typename number>
    template <class DOFINFO>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(
      const DOFINFO &info1,
      const DOFINFO &info2)
    {
      const unsigned int level1 = info1.cell->level();
      const unsigned int level2 = info2.cell->level();

      for (unsigned int i=0; i<matrices->size(); ++i)
        {
          MGLevelObject<MatrixBlock<MATRIX> > &o = matrices->block(i);

          // Row and column index of
          // the block we are dealing with
          const unsigned int row = o[level1].row;
          const unsigned int col = o[level1].column;

          if (level1 == level2)
            {
              if (mg_constrained_dofs == 0)
                {
                  assemble(o[level1].matrix, info1.matrix(i,false).matrix, row, col,
                           info1.indices, info1.indices, level1, level1);
                  assemble(o[level1].matrix, info1.matrix(i,true).matrix, row, col,
                           info1.indices, info2.indices, level1, level2);
                  assemble(o[level1].matrix, info2.matrix(i,false).matrix, row, col,
                           info2.indices, info2.indices, level2, level2);
                  assemble(o[level1].matrix, info2.matrix(i,true).matrix, row, col,
                           info2.indices, info1.indices, level2, level1);
                }
              else
                {
                  assemble_fluxes(o[level1].matrix, info1.matrix(i,false).matrix, row, col,
                                  info1.indices, info1.indices, level1, level1);
                  assemble_fluxes(o[level1].matrix, info1.matrix(i,true).matrix, row, col,
                                  info1.indices, info2.indices, level1, level2);
                  assemble_fluxes(o[level1].matrix, info2.matrix(i,false).matrix, row, col,
                                  info2.indices, info2.indices, level2, level2);
                  assemble_fluxes(o[level1].matrix, info2.matrix(i,true).matrix, row, col,
                                  info2.indices, info1.indices, level2, level1);
                }
            }
          else
            {
              Assert(level1 > level2, ExcNotImplemented());
              if (flux_up->size() != 0)
                {
                  // Do not add M22,
                  // which is done by
                  // the coarser cell
                  assemble_fluxes(o[level1].matrix, info1.matrix(i,false).matrix, row, col,
                                  info1.indices, info1.indices, level1, level1);
                  assemble_up(flux_up->block(i)[level1].matrix, info1.matrix(i,true).matrix, row, col,
                              info1.indices, info2.indices, level1, level2);
                  assemble_down(flux_down->block(i)[level1].matrix, info2.matrix(i,true).matrix, row, col,
                                info2.indices, info1.indices, level2, level1);
                }
            }
        }
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
