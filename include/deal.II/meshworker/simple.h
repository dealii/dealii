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


#ifndef __deal2__mesh_worker_simple_h
#define __deal2__mesh_worker_simple_h

#include <deal.II/base/named_data.h>
#include <deal.II/algorithms/any_data.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/functional.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>

/**
 * @file
 * @brief MeshWorker::Assember objects for systems without block structure
 *
 * The header containing the classes
 * MeshWorker::Assember::MatrixSimple,
 * MeshWorker::Assember::MGMatrixSimple,
 * MeshWorker::Assember::ResidualSimple, and
 * MeshWorker::Assember::SystemSimple.
 */

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  namespace Assembler
  {
    /**
     * Assemble residuals without block structure.
     *
     * The data structure for this Assembler class is a simple vector on
     * each cell with entries from zero to
     * FiniteElementData::dofs_per_cell and a simple global vector with
     * entries numbered from zero to DoFHandler::n_dofs(). No BlockInfo is
     * required and the global vector may be any type of vector having
     * element access through <tt>operator() (unsigned int)</tt>
     *
     * @ingroup MeshWorker
     * @author Guido Kanschat, 2009
     */
    template <class VECTOR>
    class ResidualSimple
    {
    public:
      /**
       * Initialize with an AnyData object holding the result of
       * assembling.
       *
       * Assembling currently writes into the first vector of <tt>results</tt>.
       */
      void initialize(AnyData &results);

      /**
       * @deprecated Use initialize(AnyData&) instead.
       */
      void initialize(NamedData<VECTOR *> &results);
      /**
       * Initialize the constraints.
       */
      void initialize(const ConstraintMatrix &constraints);

      /**
       * @deprecated This function is of no effect. Only the block info
       * structure in DoFInfo is being used.
       *
       * Store information on the local block structure. If the
       * assembler is inititialized with this function,
       * initialize_info() will generate one local matrix for each
       * block row and column, which will be numbered
       * lexicographically, row by row.
       *
       * In spite of using local block structure, all blocks will be
       * enteres into the same global matrix, disregarding any global
       * block structure.
       */

      void initialize_local_blocks(const BlockIndices &);

      /**
       * Initialize the local data in the DoFInfo object used later
       * for assembling.
       *
       * The info object refers to a cell if <code>!face</code>, or
       * else to an interior or boundary face.
       */
      template <class DOFINFO>
      void initialize_info(DOFINFO &info, bool face) const;

      /**
       * Assemble the local residuals into the global residuals.
       *
       * Values are added to the previous contents. If constraints are
       * active, ConstraintMatrix::distribute_local_to_global() is
       * used.
       */
      template <class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble both local residuals into the global residuals.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);
    private:
      /**
       * The global residal vectors filled by assemble().
       */
      AnyData residuals;
      /**
       * A pointer to the object containing constraints.
       */
      SmartPointer<const ConstraintMatrix,ResidualSimple<VECTOR> > constraints;
    };


    /**
     * Assemble local matrices into a single global matrix. If this global
     * matrix has a block structure, this structure is not used, but
     * rather the global numbering of degrees of freedom.
     *
     * After being initialized with a SparseMatrix object (or another
     * matrix offering the same functionality as SparseMatrix::add()),
     * this class can be used in a MeshWorker::loop() to assemble the cell
     * and face matrices into the global matrix.
     *
     * If a ConstraintMatrix has been provided during initialization, this
     * matrix will be used
     * (ConstraintMatrix::distribute_local_to_global(), to be precise) to
     * enter the local matrix into the global sparse matrix.
     *
     * The assembler can handle two different types of local data. First,
     * by default, the obvious choice of taking a single local matrix with
     * dimensions equal to the number of degrees of freedom of the
     * cell. Alternatively, a local block structure can be initialized
     * in DoFInfo. After this, the local data will be
     * arranged as an array of n by n FullMatrix blocks, which are
     * ordered lexicographically in DoFInfo.
     *
     * @ingroup MeshWorker
     * @author Guido Kanschat, 2009
     */
    template <class MATRIX>
    class MatrixSimple
    {
    public:
      /**
       * Constructor, initializing the #threshold, which limits how
       * small numbers may be to be entered into the matrix.
       */
      MatrixSimple(double threshold = 1.e-12);

      /**
       * Store the result matrix for later assembling.
       */
      void initialize(MATRIX &m);
      /**
       * Initialize the constraints. After this function has been
       * called with a valid ConstraintMatrix, the function
       * ConstraintMatrix::distribute_local_to_global() will be used
       * by assemble() to distribute the cell and face matrices into a
       * global sparse matrix.
       */
      void initialize(const ConstraintMatrix &constraints);

      /**
       * @deprecated This function is of no effect. Only the block info
       * structure in DoFInfo is being used.
       *
       * Store information on the local block structure. If the
       * assembler is inititialized with this function,
       * initialize_info() will generate one local matrix for each
       * block row and column, which will be numbered
       * lexicographically, row by row.
       *
       * In spite of using local block structure, all blocks will be
       * enteres into the same global matrix, disregarding any global
       * block structure.
       */

      void initialize_local_blocks(const BlockIndices &);

      /**
       * Initialize the local data in the DoFInfo object used later
       * for assembling.
       *
       * The info object refers to a cell if <code>!face</code>, or
       * else to an interior or boundary face.
       */
      template <class DOFINFO>
      void initialize_info(DOFINFO &info, bool face) const;

      /**
       * Assemble the matrix DoFInfo::M1[0] into the global matrix.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble both local matrices in the info objects into the
       * global matrix.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);
    private:
      /**
       * Assemble a single matrix
       * into #matrix.
       */
      void assemble(const FullMatrix<double> &M,
                    const std::vector<types::global_dof_index> &i1,
                    const std::vector<types::global_dof_index> &i2);

      /**
       * The global matrix being
       * assembled.
       */
      SmartPointer<MATRIX,MatrixSimple<MATRIX> > matrix;
      /**
       * A pointer to the object
       * containing constraints.
       */
      SmartPointer<const ConstraintMatrix,MatrixSimple<MATRIX> > constraints;

      /**
       * The smallest positive number that will be entered into the
       * global matrix. All smaller absolute values will be treated as
       * zero and will not be assembled.
       */
      const double threshold;

    };


    /**
     * Assemble local matrices into level matrices without using
     * block structure.
     *
     * @todo The matrix structures needed for assembling level matrices
     * with local refinement and continuous elements are missing.
     *
     * @ingroup MeshWorker
     * @author Guido Kanschat, 2009
     */
    template <class MATRIX>
    class MGMatrixSimple
    {
    public:
      /**
       * Constructor, initializing the #threshold, which limits how
       * small numbers may be to be entered into the matrix.
       */
      MGMatrixSimple(double threshold = 1.e-12);

      /**
       * Store the result matrix for later assembling.
       */
      void initialize(MGLevelObject<MATRIX> &m);

      /**
       * Initialize the multilevel constraints.
       */
      void initialize(const MGConstrainedDoFs &mg_constrained_dofs);

      /**
       * @deprecated This function is of no effect. Only the block
       * info structure in DoFInfo is being used.
       *
       * Store information on the local block structure. If the
       * assembler is inititialized with this function,
       * initialize_info() will generate one local matrix for each
       * block row and column, which will be numbered
       * lexicographically, row by row.
       *
       * In spite of using local block structure, all blocks will be
       * enteres into the same global matrix, disregarding any global
       * block structure.
       */
      void initialize_local_blocks(const BlockIndices &);

      /**
       * Initialize the matrices #flux_up and #flux_down used for
       * local refinement with discontinuous Galerkin methods.
       */
      void initialize_fluxes(MGLevelObject<MATRIX> &flux_up,
                             MGLevelObject<MATRIX> &flux_down);

      /**
       * Initialize the matrices #interface_in and #interface_out used
       * for local refinement with continuous Galerkin methods.
       */

      void initialize_interfaces(MGLevelObject<MATRIX> &interface_in,
                                 MGLevelObject<MATRIX> &interface_out);
      /**
       * Initialize the local data in the DoFInfo object used later
       * for assembling.
       *
       * The info object refers to a cell if <code>!face</code>, or
       * else to an interior or boundary face.
       */
      template <class DOFINFO>
      void initialize_info(DOFINFO &info, bool face) const;

      /**
       * Assemble the matrix DoFInfo::M1[0] into the global matrix.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble both local matrices in the info objects into the
       * global matrices.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);
    private:
      /**
       * Assemble a single matrix into a global matrix.
       */
      void assemble(MATRIX &G,
                    const FullMatrix<double> &M,
                    const std::vector<types::global_dof_index> &i1,
                    const std::vector<types::global_dof_index> &i2);

      /**
       * Assemble a single matrix into a global matrix.
       */
      void assemble(MATRIX &G,
                    const FullMatrix<double> &M,
                    const std::vector<types::global_dof_index> &i1,
                    const std::vector<types::global_dof_index> &i2,
                    const unsigned int level);

      /**
       * Assemble a single matrix into a global matrix.
       */

      void assemble_up(MATRIX &G,
                       const FullMatrix<double> &M,
                       const std::vector<types::global_dof_index> &i1,
                       const std::vector<types::global_dof_index> &i2,
                       const unsigned int level = numbers::invalid_unsigned_int);
      /**
       * Assemble a single matrix into a global matrix.
       */

      void assemble_down(MATRIX &G,
                         const FullMatrix<double> &M,
                         const std::vector<types::global_dof_index> &i1,
                         const std::vector<types::global_dof_index> &i2,
                         const unsigned int level = numbers::invalid_unsigned_int);

      /**
       * Assemble a single matrix into a global matrix.
       */

      void assemble_in(MATRIX &G,
                       const FullMatrix<double> &M,
                       const std::vector<types::global_dof_index> &i1,
                       const std::vector<types::global_dof_index> &i2,
                       const unsigned int level = numbers::invalid_unsigned_int);

      /**
       * Assemble a single matrix into a global matrix.
       */

      void assemble_out(MATRIX &G,
                        const FullMatrix<double> &M,
                        const std::vector<types::global_dof_index> &i1,
                        const std::vector<types::global_dof_index> &i2,
                        const unsigned int level = numbers::invalid_unsigned_int);

      /**
       * The global matrix being assembled.
       */
      SmartPointer<MGLevelObject<MATRIX>,MGMatrixSimple<MATRIX> > matrix;

      /**
       * The matrix used for face flux terms across the refinement
       * edge, coupling coarse to fine.
       */
      SmartPointer<MGLevelObject<MATRIX>,MGMatrixSimple<MATRIX> > flux_up;

      /**
       * The matrix used for face flux terms across the refinement
       * edge, coupling fine to coarse.
       */
      SmartPointer<MGLevelObject<MATRIX>,MGMatrixSimple<MATRIX> > flux_down;

      /**
       * The matrix used for face contributions for continuous
       * elements across the refinement edge, coupling coarse to fine.
       */
      SmartPointer<MGLevelObject<MATRIX>,MGMatrixSimple<MATRIX> > interface_in;

      /**
       * The matrix used for face contributions for continuous
       * elements across the refinement edge, coupling fine to coarse.
       */
      SmartPointer<MGLevelObject<MATRIX>,MGMatrixSimple<MATRIX> > interface_out;
      /**
       * A pointer to the object containing constraints.
       */
      SmartPointer<const MGConstrainedDoFs,MGMatrixSimple<MATRIX> > mg_constrained_dofs;

      /**
       * The smallest positive number that will be entered into the
       * global matrix. All smaller absolute values will be treated as
       * zero and will not be assembled.
       */
      const double threshold;

    };


    /**
     * Assemble a simple matrix and a simple right hand side at once. We
     * use a combination of MatrixSimple and ResidualSimple to achieve
     * this. Cell and face operators should fill the matrix and vector
     * objects in LocalResults and this class will assemble
     * them into matrix and vector objects.
     *
     * @ingroup MeshWorker
     * @author Guido Kanschat, 2009
     */
    template <class MATRIX, class VECTOR>
    class SystemSimple :
      private MatrixSimple<MATRIX>,
      private ResidualSimple<VECTOR>
    {
    public:
      /**
       * Constructor setting the
       * threshold value in
       * MatrixSimple.
       */
      SystemSimple(double threshold = 1.e-12);

      /**
       * Store the two objects data
       * is assembled into.
       */
      void initialize(MATRIX &m, VECTOR &rhs);

      /**
       * Initialize the constraints. After this function has been
       * called with a valid ConstraintMatrix, the function
       * ConstraintMatrix::distribute_local_to_global() will be used
       * by assemble() to distribute the cell and face matrices into a
       * global sparse matrix.
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
       * Assemble the matrix
       * DoFInfo::M1[0]
       * into the global matrix.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble both local
       * matrices in the info
       * objects into the global
       * matrix.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);
    };


//----------------------------------------------------------------------//

    template <class VECTOR>
    inline void
    ResidualSimple<VECTOR>::initialize(AnyData &results)
    {
      residuals = results;
    }

    template <class VECTOR>
    inline void
    ResidualSimple<VECTOR>::initialize(NamedData<VECTOR *> &results)
    {
      residuals = results;
    }

    template <class VECTOR>
    inline void
    ResidualSimple<VECTOR>::initialize(const ConstraintMatrix &c)
    {
      constraints = &c;
    }


    template <class MATRIX>
    inline void
    ResidualSimple<MATRIX>::initialize_local_blocks(const BlockIndices &)
    {}


    template <class VECTOR>
    template <class DOFINFO>
    inline void
    ResidualSimple<VECTOR>::initialize_info(DOFINFO &info, bool) const
    {
      info.initialize_vectors(residuals.size());
    }


    template <class VECTOR>
    template <class DOFINFO>
    inline void
    ResidualSimple<VECTOR>::assemble(const DOFINFO &info)
    {
      for (unsigned int k=0; k<residuals.size(); ++k)
        {
          VECTOR *v = residuals.entry<VECTOR *>(k);
          if (constraints == 0)
            {
              for (unsigned int i=0; i<info.vector(k).block(0).size(); ++i)
                (*v)(info.indices[i]) += info.vector(k).block(0)(i);
            }
          else
            {
              if (info.indices_by_block.size() == 0)
                constraints->distribute_local_to_global(info.vector(k).block(0), info.indices, *v);
              else
                for (unsigned int i=0; i != info.vector(k).n_blocks(); ++i)
                  constraints->distribute_local_to_global(info.vector(k).block(i), info.indices_by_block[i], *v);
            }
        }
    }

    template <class VECTOR>
    template <class DOFINFO>
    inline void
    ResidualSimple<VECTOR>::assemble(const DOFINFO &info1,
                                     const DOFINFO &info2)
    {
      for (unsigned int k=0; k<residuals.size(); ++k)
        {
          VECTOR *v = residuals.entry<VECTOR *>(k);
          if (constraints == 0)
            {
              for (unsigned int i=0; i<info1.vector(k).block(0).size(); ++i)
                (*v)(info1.indices[i]) += info1.vector(k).block(0)(i);
              for (unsigned int i=0; i<info2.vector(k).block(0).size(); ++i)
                (*v)(info2.indices[i]) += info2.vector(k).block(0)(i);
            }
          else
            {
              if (info1.indices_by_block.size() == 0 && info2.indices_by_block.size() == 0)
                {
                  constraints->distribute_local_to_global
                  (info1.vector(k).block(0), info1.indices, *v);
                  constraints->distribute_local_to_global
                  (info2.vector(k).block(0), info2.indices, *v);
                }
              else if (info1.indices_by_block.size() != 0 && info2.indices_by_block.size() != 0)
                {
                  for (unsigned int i=0; i<info1.vector(k).n_blocks(); ++i)
                    {
                      constraints->distribute_local_to_global
                      (info1.vector(k).block(i), info1.indices_by_block[i], *v);
                      constraints->distribute_local_to_global
                      (info2.vector(k).block(i), info2.indices_by_block[i], *v);
                    }
                }
              else
                {
                  Assert(false, ExcNotImplemented());
                }
            }
        }
    }


//----------------------------------------------------------------------//

    template <class MATRIX>
    inline
    MatrixSimple<MATRIX>::MatrixSimple(double threshold)
      :
      threshold(threshold)
    {}


    template <class MATRIX>
    inline void
    MatrixSimple<MATRIX>::initialize(MATRIX &m)
    {
      matrix = &m;
    }


    template <class MATRIX>
    inline void
    MatrixSimple<MATRIX>::initialize(const ConstraintMatrix &c)
    {
      constraints = &c;
    }


    template <class MATRIX>
    inline void
    MatrixSimple<MATRIX>::initialize_local_blocks(const BlockIndices &)
    {}


    template <class MATRIX >
    template <class DOFINFO>
    inline void
    MatrixSimple<MATRIX>::initialize_info(DOFINFO &info, bool face) const
    {
      const unsigned int n = info.indices_by_block.size();

      if (n == 0)
        info.initialize_matrices(1, face);
      else
        {
          info.initialize_matrices(n*n, face);
          unsigned int k=0;
          for (unsigned int i=0; i<n; ++i)
            for (unsigned int j=0; j<n; ++j,++k)
              {
                info.matrix(k,false).row = i;
                info.matrix(k,false).column = j;
                if (face)
                  {
                    info.matrix(k,true).row = i;
                    info.matrix(k,true).column = j;
                  }
              }
        }
    }



    template <class MATRIX>
    inline void
    MatrixSimple<MATRIX>::assemble(const FullMatrix<double> &M,
                                   const std::vector<types::global_dof_index> &i1,
                                   const std::vector<types::global_dof_index> &i2)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());

      if (constraints == 0)
        {
          for (unsigned int j=0; j<i1.size(); ++j)
            for (unsigned int k=0; k<i2.size(); ++k)
              if (std::fabs(M(j,k)) >= threshold)
                matrix->add(i1[j], i2[k], M(j,k));
        }
      else
        constraints->distribute_local_to_global(M, i1, i2, *matrix);
    }


    template <class MATRIX>
    template <class DOFINFO>
    inline void
    MatrixSimple<MATRIX>::assemble(const DOFINFO &info)
    {
      Assert(!info.level_cell, ExcMessage("Cell may not access level dofs"));

      if (info.indices_by_block.size() == 0)
        assemble(info.matrix(0,false).matrix, info.indices, info.indices);
      else
        {
          for (unsigned int k=0; k<info.n_matrices(); ++k)
            {
              assemble(info.matrix(k,false).matrix,
                       info.indices_by_block[info.matrix(k,false).row],
                       info.indices_by_block[info.matrix(k,false).column]);
            }
        }
    }


    template <class MATRIX>
    template <class DOFINFO>
    inline void
    MatrixSimple<MATRIX>::assemble(const DOFINFO &info1, const DOFINFO &info2)
    {
      Assert(!info1.level_cell, ExcMessage("Cell may not access level dofs"));
      Assert(!info2.level_cell, ExcMessage("Cell may not access level dofs"));

      if (info1.indices_by_block.size() == 0 && info2.indices_by_block.size() == 0)
        {
          assemble(info1.matrix(0,false).matrix, info1.indices, info1.indices);
          assemble(info1.matrix(0,true).matrix, info1.indices, info2.indices);
          assemble(info2.matrix(0,false).matrix, info2.indices, info2.indices);
          assemble(info2.matrix(0,true).matrix, info2.indices, info1.indices);
        }
      else if (info1.indices_by_block.size() != 0 && info2.indices_by_block.size() != 0)
        for (unsigned int k=0; k<info1.n_matrices(); ++k)
          {
            const unsigned int row = info1.matrix(k,false).row;
            const unsigned int column = info1.matrix(k,false).column;

            assemble(info1.matrix(k,false).matrix,
                     info1.indices_by_block[row], info1.indices_by_block[column]);
            assemble(info1.matrix(k,true).matrix,
                     info1.indices_by_block[row], info2.indices_by_block[column]);
            assemble(info2.matrix(k,false).matrix,
                     info2.indices_by_block[row], info2.indices_by_block[column]);
            assemble(info2.matrix(k,true).matrix,
                     info2.indices_by_block[row], info1.indices_by_block[column]);
          }
      else
        {
          Assert(false, ExcNotImplemented());
        }
    }


//----------------------------------------------------------------------//

    template <class MATRIX>
    inline
    MGMatrixSimple<MATRIX>::MGMatrixSimple(double threshold)
      :
      threshold(threshold)
    {}


    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::initialize(MGLevelObject<MATRIX> &m)
    {
      matrix = &m;
    }

    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::initialize(const MGConstrainedDoFs &c)
    {
      mg_constrained_dofs = &c;
    }

    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::initialize_local_blocks(const BlockIndices &)
    {}


    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::initialize_fluxes(
      MGLevelObject<MATRIX> &up, MGLevelObject<MATRIX> &down)
    {
      flux_up = &up;
      flux_down = &down;
    }


    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::initialize_interfaces(
      MGLevelObject<MATRIX> &in, MGLevelObject<MATRIX> &out)
    {
      interface_in = &in;
      interface_out = &out;
    }


    template <class MATRIX >
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MATRIX>::initialize_info(DOFINFO &info, bool face) const
    {
      const unsigned int n = info.indices_by_block.size();

      if (n == 0)
        info.initialize_matrices(1, face);
      else
        {
          info.initialize_matrices(n*n, face);
          unsigned int k=0;
          for (unsigned int i=0; i<n; ++i)
            for (unsigned int j=0; j<n; ++j,++k)
              {
                info.matrix(k,false).row = i;
                info.matrix(k,false).column = j;
                if (face)
                  {
                    info.matrix(k,true).row = i;
                    info.matrix(k,true).column = j;
                  }
              }
        }
    }


    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble(
      MATRIX &G,
      const FullMatrix<double> &M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());
      Assert(mg_constrained_dofs == 0, ExcInternalError());
//TODO: Possibly remove this function all together

      for (unsigned int j=0; j<i1.size(); ++j)
        for (unsigned int k=0; k<i2.size(); ++k)
          if (std::fabs(M(j,k)) >= threshold)
            G.add(i1[j], i2[k], M(j,k));
    }


    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble(
      MATRIX &G,
      const FullMatrix<double> &M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2,
      const unsigned int level)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());

      if (mg_constrained_dofs == 0)
        {
          for (unsigned int j=0; j<i1.size(); ++j)
            for (unsigned int k=0; k<i2.size(); ++k)
              if (std::fabs(M(j,k)) >= threshold)
                G.add(i1[j], i2[k], M(j,k));
        }
      else
        {
          for (unsigned int j=0; j<i1.size(); ++j)
            for (unsigned int k=0; k<i2.size(); ++k)
              {
                // Only enter the local values into the global matrix,
                //  if the value is larger than the threshold
                if (std::fabs(M(j,k)) < threshold)
                  continue;

                // Do not enter, if either the row or the column
                // corresponds to an index on the refinement edge. The
                // level problems are solved with homogeneous
                // Dirichlet boundary conditions, therefore we
                // eliminate these rows and columns. The corresponding
                // matrix entries are entered by assemble_in() and
                // assemble_out().
                if (mg_constrained_dofs->at_refinement_edge(level, i1[j]) ||
                    mg_constrained_dofs->at_refinement_edge(level, i2[k]))
                  continue;

                // At the boundary, only enter the term on the
                // diagonal, but not the coupling terms
                if ((mg_constrained_dofs->is_boundary_index(level, i1[j]) ||
                     mg_constrained_dofs->is_boundary_index(level, i2[k])) &&
                    (i1[j] != i2[k]))
                  continue;

                G.add(i1[j], i2[k], M(j,k));
              }
        }
    }


    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble_up(
      MATRIX &G,
      const FullMatrix<double> &M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2,
      const unsigned int level)
    {
      AssertDimension(M.n(), i1.size());
      AssertDimension(M.m(), i2.size());

      if (mg_constrained_dofs == 0)
        {
          for (unsigned int j=0; j<i1.size(); ++j)
            for (unsigned int k=0; k<i2.size(); ++k)
              if (std::fabs(M(k,j)) >= threshold)
                G.add(i1[j], i2[k], M(k,j));
        }
      else
        {
          for (unsigned int j=0; j<i1.size(); ++j)
            for (unsigned int k=0; k<i2.size(); ++k)
              if (std::fabs(M(k,j)) >= threshold)
                if (!mg_constrained_dofs->at_refinement_edge(level, i2[k]))
                  G.add(i1[j], i2[k], M(k,j));
        }
    }

    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble_down(
      MATRIX &G,
      const FullMatrix<double> &M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2,
      const unsigned int level)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());

      if (mg_constrained_dofs == 0)
        {
          for (unsigned int j=0; j<i1.size(); ++j)
            for (unsigned int k=0; k<i2.size(); ++k)
              if (std::fabs(M(j,k)) >= threshold)
                G.add(i1[j], i2[k], M(j,k));
        }
      else
        {
          for (unsigned int j=0; j<i1.size(); ++j)
            for (unsigned int k=0; k<i2.size(); ++k)
              if (std::fabs(M(j,k)) >= threshold)
                if (!mg_constrained_dofs->at_refinement_edge(level, i2[k]))
                  G.add(i1[j], i2[k], M(j,k));
        }
    }

    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble_in(
      MATRIX &G,
      const FullMatrix<double> &M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2,
      const unsigned int level)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());
      Assert(mg_constrained_dofs != 0, ExcInternalError());

      for (unsigned int j=0; j<i1.size(); ++j)
        for (unsigned int k=0; k<i2.size(); ++k)
          if (std::fabs(M(j,k)) >= threshold)
            // Enter values into matrix only if j corresponds to a
            // degree of freedom on the refinemenent edge, k does
            // not, and both are not on the boundary. This is part
            // the difference between the complete matrix with no
            // boundary condition at the refinement edge and and
            // the matrix assembled above by assemble().

            // Thus the logic is: enter the row if it is
            // constrained by hanging node constraints (actually,
            // the whole refinement edge), but not if it is
            // constrained by a boundary constraint.
            if (mg_constrained_dofs->at_refinement_edge(level, i1[j]) &&
                !mg_constrained_dofs->at_refinement_edge(level, i2[k]))
              {
                if ((!mg_constrained_dofs->is_boundary_index(level, i1[j]) &&
                     !mg_constrained_dofs->is_boundary_index(level, i2[k]))
                    ||
                    (mg_constrained_dofs->is_boundary_index(level, i1[j]) &&
                     mg_constrained_dofs->is_boundary_index(level, i2[k]) &&
                     i1[j] == i2[k]))
                  G.add(i1[j], i2[k], M(j,k));
              }
    }


    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble_out(
      MATRIX &G,
      const FullMatrix<double> &M,
      const std::vector<types::global_dof_index> &i1,
      const std::vector<types::global_dof_index> &i2,
      const unsigned int level)
    {
      AssertDimension(M.n(), i1.size());
      AssertDimension(M.m(), i2.size());
      Assert(mg_constrained_dofs != 0, ExcInternalError());

      for (unsigned int j=0; j<i1.size(); ++j)
        for (unsigned int k=0; k<i2.size(); ++k)
          if (std::fabs(M(k,j)) >= threshold)
            if (mg_constrained_dofs->at_refinement_edge(level, i1[j]) &&
                !mg_constrained_dofs->at_refinement_edge(level, i2[k]))
              {
                if ((!mg_constrained_dofs->is_boundary_index(level, i1[j]) &&
                     !mg_constrained_dofs->is_boundary_index(level, i2[k]))
                    ||
                    (mg_constrained_dofs->is_boundary_index(level, i1[j]) &&
                     mg_constrained_dofs->is_boundary_index(level, i2[k]) &&
                     i1[j] == i2[k]))
                  G.add(i1[j], i2[k], M(k,j));
              }
    }


    template <class MATRIX>
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MATRIX>::assemble(const DOFINFO &info)
    {
      Assert(info.level_cell, ExcMessage("Cell must access level dofs"));
      const unsigned int level = info.cell->level();

      if (info.indices_by_block.size() == 0)
        {
          assemble((*matrix)[level], info.matrix(0,false).matrix,
                   info.indices, info.indices, level);
          if (mg_constrained_dofs != 0)
            {
              assemble_in((*interface_in)[level], info.matrix(0,false).matrix,
                          info.indices, info.indices, level);
              assemble_out((*interface_out)[level],info.matrix(0,false).matrix,
                           info.indices, info.indices, level);
            }
        }
      else
        for (unsigned int k=0; k<info.n_matrices(); ++k)
          {
            const unsigned int row = info.matrix(k,false).row;
            const unsigned int column = info.matrix(k,false).column;

            assemble((*matrix)[level], info.matrix(k,false).matrix,
                     info.indices_by_block[row], info.indices_by_block[column], level);

            if (mg_constrained_dofs != 0)
              {
                assemble_in((*interface_in)[level], info.matrix(k,false).matrix,
                            info.indices_by_block[row], info.indices_by_block[column], level);
                assemble_out((*interface_out)[level],info.matrix(k,false).matrix,
                             info.indices_by_block[column], info.indices_by_block[row], level);
              }
          }
    }


    template <class MATRIX>
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MATRIX>::assemble(const DOFINFO &info1,
                                     const DOFINFO &info2)
    {
      Assert(info1.level_cell, ExcMessage("Cell must access level dofs"));
      Assert(info2.level_cell, ExcMessage("Cell must access level dofs"));
      const unsigned int level1 = info1.cell->level();
      const unsigned int level2 = info2.cell->level();

      if (info1.indices_by_block.size() == 0)
        {
          if (level1 == level2)
            {
              assemble((*matrix)[level1], info1.matrix(0,false).matrix, info1.indices, info1.indices, level1);
              assemble((*matrix)[level1], info1.matrix(0,true).matrix, info1.indices, info2.indices, level1);
              assemble((*matrix)[level1], info2.matrix(0,false).matrix, info2.indices, info2.indices, level1);
              assemble((*matrix)[level1], info2.matrix(0,true).matrix, info2.indices, info1.indices, level1);
            }
          else
            {
              Assert(level1 > level2, ExcInternalError());
              // Do not add info2.M1,
              // which is done by
              // the coarser cell
              assemble((*matrix)[level1], info1.matrix(0,false).matrix, info1.indices, info1.indices, level1);
              if (level1>0)
                {
                  assemble_up((*flux_up)[level1],info1.matrix(0,true).matrix, info2.indices, info1.indices, level1);
                  assemble_down((*flux_down)[level1], info2.matrix(0,true).matrix, info2.indices, info1.indices, level1);
                }
            }
        }
      else
        for (unsigned int k=0; k<info1.n_matrices(); ++k)
          {
            const unsigned int row = info1.matrix(k,false).row;
            const unsigned int column = info1.matrix(k,false).column;

            if (level1 == level2)
              {
                assemble((*matrix)[level1], info1.matrix(k,false).matrix, info1.indices_by_block[row], info1.indices_by_block[column], level1);
                assemble((*matrix)[level1], info1.matrix(k,true).matrix, info1.indices_by_block[row], info2.indices_by_block[column], level1);
                assemble((*matrix)[level1], info2.matrix(k,false).matrix, info2.indices_by_block[row], info2.indices_by_block[column], level1);
                assemble((*matrix)[level1], info2.matrix(k,true).matrix, info2.indices_by_block[row], info1.indices_by_block[column], level1);
              }
            else
              {
                Assert(level1 > level2, ExcInternalError());
                // Do not add info2.M1,
                // which is done by
                // the coarser cell
                assemble((*matrix)[level1], info1.matrix(k,false).matrix, info1.indices_by_block[row], info1.indices_by_block[column], level1);
                if (level1>0)
                  {
                    assemble_up((*flux_up)[level1],info1.matrix(k,true).matrix, info2.indices_by_block[row], info1.indices_by_block[column], level1);
                    assemble_down((*flux_down)[level1], info2.matrix(k,true).matrix, info2.indices_by_block[row], info1.indices_by_block[column], level1);
                  }
              }
          }
    }

//----------------------------------------------------------------------//

    template <class MATRIX, class VECTOR>
    SystemSimple<MATRIX,VECTOR>::SystemSimple(double t)
      :
      MatrixSimple<MATRIX>(t)
    {}


    template <class MATRIX, class VECTOR>
    inline void
    SystemSimple<MATRIX,VECTOR>::initialize(MATRIX &m, VECTOR &rhs)
    {
      AnyData data;
      VECTOR *p = &rhs;
      data.add(p, "right hand side");

      MatrixSimple<MATRIX>::initialize(m);
      ResidualSimple<VECTOR>::initialize(data);
    }

    template <class MATRIX, class VECTOR>
    inline void
    SystemSimple<MATRIX,VECTOR>::initialize(const ConstraintMatrix &c)
    {
      MatrixSimple<MATRIX>::initialize(c);
      ResidualSimple<VECTOR>::initialize(c);
    }


    template <class MATRIX, class VECTOR>
    template <class DOFINFO>
    inline void
    SystemSimple<MATRIX,VECTOR>::initialize_info(DOFINFO &info,
                                                 bool face) const
    {
      MatrixSimple<MATRIX>::initialize_info(info, face);
      ResidualSimple<VECTOR>::initialize_info(info, face);
    }


    template <class MATRIX, class VECTOR>
    template <class DOFINFO>
    inline void
    SystemSimple<MATRIX,VECTOR>::assemble(const DOFINFO &info)
    {
      MatrixSimple<MATRIX>::assemble(info);
      ResidualSimple<VECTOR>::assemble(info);
    }


    template <class MATRIX, class VECTOR>
    template <class DOFINFO>
    inline void
    SystemSimple<MATRIX,VECTOR>::assemble(const DOFINFO &info1,
                                          const DOFINFO &info2)
    {
      MatrixSimple<MATRIX>::assemble(info1, info2);
      ResidualSimple<VECTOR>::assemble(info1, info2);
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
