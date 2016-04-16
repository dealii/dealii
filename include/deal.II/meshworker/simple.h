// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2015 by the deal.II authors
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


#ifndef dealii__mesh_worker_simple_h
#define dealii__mesh_worker_simple_h

#include <deal.II/algorithms/any_data.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/functional.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>

/*
 * The header containing the classes MeshWorker::Assember::MatrixSimple,
 * MeshWorker::Assember::MGMatrixSimple, MeshWorker::Assember::ResidualSimple,
 * and MeshWorker::Assember::SystemSimple.
 */

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  namespace Assembler
  {
    /**
     * Assemble residuals without block structure.
     *
     * The data structure for this Assembler class is a simple vector on each
     * cell with entries from zero to FiniteElementData::dofs_per_cell and a
     * simple global vector with entries numbered from zero to
     * DoFHandler::n_dofs(). No BlockInfo is required and the global vector
     * may be any type of vector having element access through <tt>operator()
     * (unsigned int)</tt>
     *
     * @ingroup MeshWorker
     * @author Guido Kanschat, 2009
     */
    template <typename VectorType>
    class ResidualSimple
    {
    public:
      /**
       * Initialize with an AnyData object holding the result of assembling.
       *
       * Assembling currently writes into the first vector of
       * <tt>results</tt>.
       */
      void initialize(AnyData &results);

      /**
       * Initialize the constraints.
       */
      void initialize(const ConstraintMatrix &constraints);

      /**
       * @deprecated This function is of no effect. Only the block info
       * structure in DoFInfo is being used.
       *
       * Store information on the local block structure. If the assembler is
       * initialized with this function, initialize_info() will generate one
       * local matrix for each block row and column, which will be numbered
       * lexicographically, row by row.
       *
       * In spite of using local block structure, all blocks will be entered
       * into the same global matrix, disregarding any global block structure.
       */

      void initialize_local_blocks(const BlockIndices &);

      /**
       * Initialize the local data in the DoFInfo object used later for
       * assembling.
       *
       * The info object refers to a cell if <code>!face</code>, or else to an
       * interior or boundary face.
       */
      template <class DOFINFO>
      void initialize_info(DOFINFO &info, bool face) const;

      /**
       * Assemble the local residuals into the global residuals.
       *
       * Values are added to the previous contents. If constraints are active,
       * ConstraintMatrix::distribute_local_to_global() is used.
       */
      template <class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble both local residuals into the global residuals.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);
    protected:
      /**
       * The global residual vectors filled by assemble().
       */
      AnyData residuals;

      /**
       * A pointer to the object containing constraints.
       */
      SmartPointer<const ConstraintMatrix,ResidualSimple<VectorType> > constraints;
    };


    /**
     * Assemble local matrices into a single global matrix or several global
     * matrices associated with the same DoFHandler. If these global matrix
     * have a block structure, this structure is not used, but rather the
     * global numbering of degrees of freedom.
     *
     * After being initialized with a SparseMatrix object (or another matrix
     * offering the same functionality as SparseMatrix::add()) or a vector of
     * such, this class can be used in a MeshWorker::loop() to assemble the
     * cell and face matrices into the global matrix.
     *
     * If a ConstraintMatrix has been provided during initialization, this
     * matrix will be used (ConstraintMatrix::distribute_local_to_global(), to
     * be precise) to enter the local matrix into the global sparse matrix.
     *
     * The assembler can handle two different types of local data. First, by
     * default, the obvious choice of taking a single local matrix with
     * dimensions equal to the number of degrees of freedom of the cell.
     * Alternatively, a local block structure can be initialized in DoFInfo.
     * After this, the local data will be arranged as an array of n by n
     * FullMatrix blocks (n being the number of blocks in the FESystem used by
     * the DoFHandler in DoFInfo), which are ordered lexicographically with
     * column index fastest in DoFInfo. If the matrix was initialized with a
     * vector of several matrices and local block structure is used, then the
     * first n<sup>2</sup> matrices in LocalResults will be used for the first
     * matrix in this vector, the second set of n<sup>2</sup> for the second,
     * and so on.
     *
     * @ingroup MeshWorker
     * @author Guido Kanschat, 2009
     */
    template <typename MatrixType>
    class MatrixSimple
    {
    public:
      /**
       * Constructor, initializing the #threshold, which limits how small
       * numbers may be to be entered into the matrix.
       */
      MatrixSimple(double threshold = 1.e-12);

      /**
       * Store the result matrix for later assembling.
       */
      void initialize(MatrixType &m);

      /**
       * Store several result matrices for later assembling.
       */
      void initialize(std::vector<MatrixType> &m);

      /**
       * Initialize the constraints. After this function has been called with
       * a valid ConstraintMatrix, the function
       * ConstraintMatrix::distribute_local_to_global() will be used by
       * assemble() to distribute the cell and face matrices into a global
       * sparse matrix.
       */
      void initialize(const ConstraintMatrix &constraints);

      /**
       * Initialize the local data in the DoFInfo object used later for
       * assembling.
       *
       * The info object refers to a cell if <code>!face</code>, or else to an
       * interior or boundary face.
       */
      template <class DOFINFO>
      void initialize_info(DOFINFO &info, bool face) const;

      /**
       * Assemble the local matrices associated with a single cell into the
       * global matrix.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble all local matrices associated with an interior face in the
       * info objects into the global matrix.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);
    protected:
      /**
       * The vector of global matrices being assembled.
       */
      std::vector<SmartPointer<MatrixType,MatrixSimple<MatrixType> > > matrix;

      /**
        * The smallest positive number that will be entered into the global
        * matrix. All smaller absolute values will be treated as zero and will
        * not be assembled.
        */
      const double threshold;

    private:
      /**
       * Assemble a single matrix <code>M</code> into the element at
       * <code>index</code> in the vector #matrix.
       */
      void assemble(const FullMatrix<double> &M,
                    const unsigned int index,
                    const std::vector<types::global_dof_index> &i1,
                    const std::vector<types::global_dof_index> &i2);

      /**
       * A pointer to the object containing constraints.
       */
      SmartPointer<const ConstraintMatrix,MatrixSimple<MatrixType> > constraints;

    };


    /**
     * Assemble local matrices into level matrices without using block
     * structure.
     *
     * @todo The matrix structures needed for assembling level matrices with
     * local refinement and continuous elements are missing.
     *
     * @ingroup MeshWorker
     * @author Guido Kanschat, 2009
     */
    template <typename MatrixType>
    class MGMatrixSimple
    {
    public:
      /**
       * Constructor, initializing the #threshold, which limits how small
       * numbers may be to be entered into the matrix.
       */
      MGMatrixSimple(double threshold = 1.e-12);

      /**
       * Store the result matrix for later assembling.
       */
      void initialize(MGLevelObject<MatrixType> &m);

      /**
       * Initialize the multilevel constraints.
       */
      void initialize(const MGConstrainedDoFs &mg_constrained_dofs);

      /**
       * Initialize the matrices #flux_up and #flux_down used for local
       * refinement with discontinuous Galerkin methods.
       */
      void initialize_fluxes(MGLevelObject<MatrixType> &flux_up,
                             MGLevelObject<MatrixType> &flux_down);

      /**
       * Initialize the matrices #interface_in and #interface_out used for
       * local refinement with continuous Galerkin methods.
       */
      void initialize_interfaces(MGLevelObject<MatrixType> &interface_in,
                                 MGLevelObject<MatrixType> &interface_out);
      /**
       * Initialize the local data in the DoFInfo object used later for
       * assembling.
       *
       * The info object refers to a cell if <code>!face</code>, or else to an
       * interior or boundary face.
       */
      template <class DOFINFO>
      void initialize_info(DOFINFO &info, bool face) const;

      /**
       * Assemble the matrix DoFInfo::M1[0] into the global matrix.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble both local matrices in the info objects into the global
       * matrices.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);
    private:
      /**
       * Assemble a single matrix into a global matrix.
       */
      void assemble(MatrixType &G,
                    const FullMatrix<double> &M,
                    const std::vector<types::global_dof_index> &i1,
                    const std::vector<types::global_dof_index> &i2);

      /**
       * Assemble a single matrix into a global matrix.
       */
      void assemble(MatrixType &G,
                    const FullMatrix<double> &M,
                    const std::vector<types::global_dof_index> &i1,
                    const std::vector<types::global_dof_index> &i2,
                    const unsigned int level);

      /**
       * Assemble a single matrix into a global matrix.
       */

      void assemble_up(MatrixType &G,
                       const FullMatrix<double> &M,
                       const std::vector<types::global_dof_index> &i1,
                       const std::vector<types::global_dof_index> &i2,
                       const unsigned int level = numbers::invalid_unsigned_int);
      /**
       * Assemble a single matrix into a global matrix.
       */

      void assemble_down(MatrixType &G,
                         const FullMatrix<double> &M,
                         const std::vector<types::global_dof_index> &i1,
                         const std::vector<types::global_dof_index> &i2,
                         const unsigned int level = numbers::invalid_unsigned_int);

      /**
       * Assemble a single matrix into a global matrix.
       */

      void assemble_in(MatrixType &G,
                       const FullMatrix<double> &M,
                       const std::vector<types::global_dof_index> &i1,
                       const std::vector<types::global_dof_index> &i2,
                       const unsigned int level = numbers::invalid_unsigned_int);

      /**
       * Assemble a single matrix into a global matrix.
       */

      void assemble_out(MatrixType &G,
                        const FullMatrix<double> &M,
                        const std::vector<types::global_dof_index> &i1,
                        const std::vector<types::global_dof_index> &i2,
                        const unsigned int level = numbers::invalid_unsigned_int);

      /**
       * The global matrix being assembled.
       */
      SmartPointer<MGLevelObject<MatrixType>,MGMatrixSimple<MatrixType> > matrix;

      /**
       * The matrix used for face flux terms across the refinement edge,
       * coupling coarse to fine.
       */
      SmartPointer<MGLevelObject<MatrixType>,MGMatrixSimple<MatrixType> > flux_up;

      /**
       * The matrix used for face flux terms across the refinement edge,
       * coupling fine to coarse.
       */
      SmartPointer<MGLevelObject<MatrixType>,MGMatrixSimple<MatrixType> > flux_down;

      /**
       * The matrix used for face contributions for continuous elements across
       * the refinement edge, coupling coarse to fine.
       */
      SmartPointer<MGLevelObject<MatrixType>,MGMatrixSimple<MatrixType> > interface_in;

      /**
       * The matrix used for face contributions for continuous elements across
       * the refinement edge, coupling fine to coarse.
       */
      SmartPointer<MGLevelObject<MatrixType>,MGMatrixSimple<MatrixType> > interface_out;
      /**
       * A pointer to the object containing constraints.
       */
      SmartPointer<const MGConstrainedDoFs,MGMatrixSimple<MatrixType> > mg_constrained_dofs;

      /**
       * The smallest positive number that will be entered into the global
       * matrix. All smaller absolute values will be treated as zero and will
       * not be assembled.
       */
      const double threshold;

    };


    /**
     * Assemble a simple matrix and a simple right hand side at once. We use a
     * combination of MatrixSimple and ResidualSimple to achieve this. Cell
     * and face operators should fill the matrix and vector objects in
     * LocalResults and this class will assemble them into matrix and vector
     * objects.
     *
     * @ingroup MeshWorker
     * @author Guido Kanschat, 2009
     */
    template <typename MatrixType, typename VectorType>
    class SystemSimple :
      private MatrixSimple<MatrixType>,
      private ResidualSimple<VectorType>
    {
    public:
      /**
       * Constructor setting the threshold value in MatrixSimple.
       */
      SystemSimple(double threshold = 1.e-12);

      /**
       * Store the two objects data is assembled into.
       */
      void initialize(MatrixType &m, VectorType &rhs);

      /**
       * Initialize the constraints. After this function has been called with
       * a valid ConstraintMatrix, the function
       * ConstraintMatrix::distribute_local_to_global() will be used by
       * assemble() to distribute the cell and face matrices into a global
       * sparse matrix.
       */
      void initialize(const ConstraintMatrix &constraints);

      /**
       * Initialize the local data in the DoFInfo object used later for
       * assembling.
       *
       * The info object refers to a cell if <code>!face</code>, or else to an
       * interior or boundary face.
       */
      template <class DOFINFO>
      void initialize_info(DOFINFO &info, bool face) const;

      /**
       * Assemble the matrix DoFInfo::M1[0] into the global matrix.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info);

      /**
       * Assemble both local matrices in the info objects into the global
       * matrix.
       */
      template<class DOFINFO>
      void assemble(const DOFINFO &info1,
                    const DOFINFO &info2);

    private:
      /**
        * Assemble a single matrix <code>M</code> into the element at
        * <code>index</code> in the vector #matrix.
        */
      void assemble(const FullMatrix<double> &M,
                    const Vector<double> &vector,
                    const unsigned int index,
                    const std::vector<types::global_dof_index> &indices);

      void assemble(const FullMatrix<double> &M,
                    const Vector<double> &vector,
                    const unsigned int index,
                    const std::vector<types::global_dof_index> &i1,
                    const std::vector<types::global_dof_index> &i2);
    };


//----------------------------------------------------------------------//

    template <typename VectorType>
    inline void
    ResidualSimple<VectorType>::initialize(AnyData &results)
    {
      residuals = results;
    }

    template <typename VectorType>
    inline void
    ResidualSimple<VectorType>::initialize(const ConstraintMatrix &c)
    {
      constraints = &c;
    }


    template <typename MatrixType>
    inline void
    ResidualSimple<MatrixType>::initialize_local_blocks(const BlockIndices &)
    {}


    template <typename VectorType>
    template <class DOFINFO>
    inline void
    ResidualSimple<VectorType>::initialize_info(DOFINFO &info, bool) const
    {
      info.initialize_vectors(residuals.size());
    }


    template <typename VectorType>
    template <class DOFINFO>
    inline void
    ResidualSimple<VectorType>::assemble(const DOFINFO &info)
    {
      for (unsigned int k=0; k<residuals.size(); ++k)
        {
          VectorType *v = residuals.entry<VectorType *>(k);
          for (unsigned int i=0; i != info.vector(k).n_blocks(); ++i)
            {
              const std::vector<types::global_dof_index>  &ldi = info.vector(k).n_blocks()==1?
                                                                 info.indices:
                                                                 info.indices_by_block[i];

              if (constraints !=0)
                constraints->distribute_local_to_global(info.vector(k).block(i), ldi, *v);
              else
                v->add(ldi, info.vector(k).block(i));
            }
        }
    }

    template <typename VectorType>
    template <class DOFINFO>
    inline void
    ResidualSimple<VectorType>::assemble(const DOFINFO &info1,
                                         const DOFINFO &info2)
    {
      assemble(info1);
      assemble(info2);
    }


//----------------------------------------------------------------------//

    template <typename MatrixType>
    inline
    MatrixSimple<MatrixType>::MatrixSimple(double threshold)
      :
      threshold(threshold)
    {}


    template <typename MatrixType>
    inline void
    MatrixSimple<MatrixType>::initialize(MatrixType &m)
    {
      matrix.resize(1);
      matrix[0] = &m;
    }


    template <typename MatrixType>
    inline void
    MatrixSimple<MatrixType>::initialize(std::vector<MatrixType> &m)
    {
      matrix.resize(m.size());
      for (unsigned int i=0; i<m.size(); ++i)
        matrix[i] = &m[i];
    }


    template <typename MatrixType>
    inline void
    MatrixSimple<MatrixType>::initialize(const ConstraintMatrix &c)
    {
      constraints = &c;
    }


    template <typename MatrixType >
    template <class DOFINFO>
    inline void
    MatrixSimple<MatrixType>::initialize_info(DOFINFO &info, bool face) const
    {
      Assert(matrix.size() != 0, ExcNotInitialized());

      const unsigned int n = info.indices_by_block.size();

      if (n == 0)
        info.initialize_matrices(matrix.size(), face);
      else
        {
          info.initialize_matrices(matrix.size()*n*n, face);
          unsigned int k=0;
          for (unsigned int m=0; m<matrix.size(); ++m)
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



    template <typename MatrixType>
    inline void
    MatrixSimple<MatrixType>::assemble(const FullMatrix<double> &M,
                                       const unsigned int index,
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
                matrix[index]->add(i1[j], i2[k], M(j,k));
        }
      else
        constraints->distribute_local_to_global(M, i1, i2, *matrix[index]);
    }


    template <typename MatrixType>
    template <class DOFINFO>
    inline void
    MatrixSimple<MatrixType>::assemble(const DOFINFO &info)
    {
      Assert(!info.level_cell, ExcMessage("Cell may not access level dofs"));
      const unsigned int n = info.indices_by_block.size();

      if (n == 0)
        for (unsigned int m=0; m<matrix.size(); ++m)
          assemble(info.matrix(m,false).matrix, m, info.indices, info.indices);
      else
        {
          for (unsigned int m=0; m<matrix.size(); ++m)
            for (unsigned int k=0; k<n*n; ++k)
              {
                assemble(info.matrix(k+m*n*n,false).matrix, m,
                         info.indices_by_block[info.matrix(k+m*n*n,false).row],
                         info.indices_by_block[info.matrix(k+m*n*n,false).column]);
              }
        }
    }


    template <typename MatrixType>
    template <class DOFINFO>
    inline void
    MatrixSimple<MatrixType>::assemble(const DOFINFO &info1, const DOFINFO &info2)
    {
      Assert(!info1.level_cell, ExcMessage("Cell may not access level dofs"));
      Assert(!info2.level_cell, ExcMessage("Cell may not access level dofs"));
      AssertDimension(info1.indices_by_block.size(),info2.indices_by_block.size());

      const unsigned int n = info1.indices_by_block.size();

      if (n == 0)
        {
          for (unsigned int m=0; m<matrix.size(); ++m)
            {
              assemble(info1.matrix(m,false).matrix, m, info1.indices, info1.indices);
              assemble(info1.matrix(m,true).matrix, m, info1.indices, info2.indices);
              assemble(info2.matrix(m,false).matrix, m, info2.indices, info2.indices);
              assemble(info2.matrix(m,true).matrix, m, info2.indices, info1.indices);
            }
        }
      else
        {
          for (unsigned int m=0; m<matrix.size(); ++m)
            for (unsigned int k=0; k<n*n; ++k)
              {
                const unsigned int row = info1.matrix(k+m*n*n,false).row;
                const unsigned int column = info1.matrix(k+m*n*n,false).column;

                assemble(info1.matrix(k+m*n*n,false).matrix, m,
                         info1.indices_by_block[row], info1.indices_by_block[column]);
                assemble(info1.matrix(k+m*n*n,true).matrix, m,
                         info1.indices_by_block[row], info2.indices_by_block[column]);
                assemble(info2.matrix(k+m*n*n,false).matrix, m,
                         info2.indices_by_block[row], info2.indices_by_block[column]);
                assemble(info2.matrix(k+m*n*n,true).matrix, m,
                         info2.indices_by_block[row], info1.indices_by_block[column]);
              }
        }
    }


//----------------------------------------------------------------------//

    template <typename MatrixType>
    inline
    MGMatrixSimple<MatrixType>::MGMatrixSimple(double threshold)
      :
      threshold(threshold)
    {}


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::initialize(MGLevelObject<MatrixType> &m)
    {
      matrix = &m;
    }

    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::initialize(const MGConstrainedDoFs &c)
    {
      mg_constrained_dofs = &c;
    }


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::initialize_fluxes(MGLevelObject<MatrixType> &up,
                                                  MGLevelObject<MatrixType> &down)
    {
      flux_up = &up;
      flux_down = &down;
    }


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::initialize_interfaces
    (MGLevelObject<MatrixType> &in, MGLevelObject<MatrixType> &out)
    {
      interface_in = &in;
      interface_out = &out;
    }


    template <typename MatrixType >
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MatrixType>::initialize_info(DOFINFO &info, bool face) const
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


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble
    (MatrixType                                 &G,
     const FullMatrix<double>                   &M,
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


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble
    (MatrixType                                 &G,
     const FullMatrix<double>                   &M,
     const std::vector<types::global_dof_index> &i1,
     const std::vector<types::global_dof_index> &i2,
     const unsigned int                          level)
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


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble_up
    (MatrixType                                 &G,
     const FullMatrix<double>                   &M,
     const std::vector<types::global_dof_index> &i1,
     const std::vector<types::global_dof_index> &i2,
     const unsigned int                          level)
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

    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble_down
    (MatrixType                                 &G,
     const FullMatrix<double>                   &M,
     const std::vector<types::global_dof_index> &i1,
     const std::vector<types::global_dof_index> &i2,
     const unsigned int                          level)
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

    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble_in
    (MatrixType                                 &G,
     const FullMatrix<double>                   &M,
     const std::vector<types::global_dof_index> &i1,
     const std::vector<types::global_dof_index> &i2,
     const unsigned int                          level)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());
      Assert(mg_constrained_dofs != 0, ExcInternalError());

      for (unsigned int j=0; j<i1.size(); ++j)
        for (unsigned int k=0; k<i2.size(); ++k)
          if (std::fabs(M(j,k)) >= threshold)
            // Enter values into matrix only if j corresponds to a
            // degree of freedom on the refinement edge, k does
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


    template <typename MatrixType>
    inline void
    MGMatrixSimple<MatrixType>::assemble_out
    (MatrixType                                 &G,
     const FullMatrix<double>                   &M,
     const std::vector<types::global_dof_index> &i1,
     const std::vector<types::global_dof_index> &i2,
     const unsigned int                          level)
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


    template <typename MatrixType>
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MatrixType>::assemble(const DOFINFO &info)
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


    template <typename MatrixType>
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MatrixType>::assemble(const DOFINFO &info1,
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
                    assemble_up((*flux_up)[level1],info1.matrix(k,true).matrix, info2.indices_by_block[column], info1.indices_by_block[row], level1);
                    assemble_down((*flux_down)[level1], info2.matrix(k,true).matrix, info2.indices_by_block[row], info1.indices_by_block[column], level1);
                  }
              }
          }
    }

//----------------------------------------------------------------------//

    template <typename MatrixType, typename VectorType>
    SystemSimple<MatrixType,VectorType>::SystemSimple(double t)
      :
      MatrixSimple<MatrixType>(t)
    {}


    template <typename MatrixType, typename VectorType>
    inline void
    SystemSimple<MatrixType,VectorType>::initialize(MatrixType &m, VectorType &rhs)
    {
      AnyData data;
      VectorType *p = &rhs;
      data.add(p, "right hand side");

      MatrixSimple<MatrixType>::initialize(m);
      ResidualSimple<VectorType>::initialize(data);
    }

    template <typename MatrixType, typename VectorType>
    inline void
    SystemSimple<MatrixType,VectorType>::initialize(const ConstraintMatrix &c)
    {
      ResidualSimple<VectorType>::initialize(c);
    }


    template <typename MatrixType, typename VectorType>
    template <class DOFINFO>
    inline void
    SystemSimple<MatrixType,VectorType>::initialize_info(DOFINFO &info,
                                                         bool    face) const
    {
      MatrixSimple<MatrixType>::initialize_info(info, face);
      ResidualSimple<VectorType>::initialize_info(info, face);
    }

    template <typename MatrixType,typename VectorType>
    inline void
    SystemSimple<MatrixType,VectorType>::assemble(const FullMatrix<double> &M,
                                                  const Vector<double> &vector,
                                                  const unsigned int index,
                                                  const std::vector<types::global_dof_index> &indices)
    {
      AssertDimension(M.m(), indices.size());
      AssertDimension(M.n(), indices.size());

      AnyData residuals = ResidualSimple<VectorType>::residuals;
      VectorType *v = residuals.entry<VectorType *>(index);

      if (ResidualSimple<VectorType>::constraints == 0)
        {
          for (unsigned int i=0; i<indices.size(); ++i)
            (*v)(indices[i]) += vector(i);

          for (unsigned int j=0; j<indices.size(); ++j)
            for (unsigned int k=0; k<indices.size(); ++k)
              if (std::fabs(M(j,k)) >= MatrixSimple<MatrixType>::threshold)
                MatrixSimple<MatrixType>::matrix[index]->add(indices[j], indices[k], M(j,k));
        }
      else
        {
          ResidualSimple<VectorType>::constraints->distribute_local_to_global(M,vector,indices,*MatrixSimple<MatrixType>::matrix[index],*v, true);
        }
    }

    template <typename MatrixType,typename VectorType>
    inline void
    SystemSimple<MatrixType,VectorType>::assemble(const FullMatrix<double> &M,
                                                  const Vector<double> &vector,
                                                  const unsigned int index,
                                                  const std::vector<types::global_dof_index> &i1,
                                                  const std::vector<types::global_dof_index> &i2)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());

      AnyData residuals = ResidualSimple<VectorType>::residuals;
      VectorType *v = residuals.entry<VectorType *>(index);

      if (ResidualSimple<VectorType>::constraints == 0)
        {
          for (unsigned int j=0; j<i1.size(); ++j)
            for (unsigned int k=0; k<i2.size(); ++k)
              if (std::fabs(M(j,k)) >= MatrixSimple<MatrixType>::threshold)
                MatrixSimple<MatrixType>::matrix[index]->add(i1[j], i2[k], M(j,k));
        }
      else
        {
          ResidualSimple<VectorType>::constraints->distribute_local_to_global(vector, i1, i2, *v, M, false);
          ResidualSimple<VectorType>::constraints->distribute_local_to_global(M, i1, i2, *MatrixSimple<MatrixType>::matrix[index]);
        }
    }


    template <typename MatrixType, typename VectorType>
    template <class DOFINFO>
    inline void
    SystemSimple<MatrixType,VectorType>::assemble(const DOFINFO &info)
    {
      AssertDimension(MatrixSimple<MatrixType>::matrix.size(),ResidualSimple<VectorType>::residuals.size());
      Assert(!info.level_cell, ExcMessage("Cell may not access level dofs"));
      const unsigned int n = info.indices_by_block.size();

      if (n == 0)
        {
          for (unsigned int m=0; m<MatrixSimple<MatrixType>::matrix.size(); ++m)
            assemble(info.matrix(m,false).matrix,info.vector(m).block(0), m, info.indices);
        }
      else
        {
          for (unsigned int m=0; m<MatrixSimple<MatrixType>::matrix.size(); ++m)
            for (unsigned int k=0; k<n*n; ++k)
              {
                const unsigned int row = info.matrix(k+m*n*n,false).row;
                const unsigned int column = info.matrix(k+m*n*n,false).column;

                if (row == column)
                  assemble(info.matrix(k+m*n*n,false).matrix,
                           info.vector(m).block(row), m,
                           info.indices_by_block[row]);
                else
                  assemble(info.matrix(k+m*n*n,false).matrix,
                           info.vector(m).block(row), m,
                           info.indices_by_block[row],
                           info.indices_by_block[column]);
              }
        }

    }


    template <typename MatrixType, typename VectorType>
    template <class DOFINFO>
    inline void
    SystemSimple<MatrixType,VectorType>::assemble(const DOFINFO &info1,
                                                  const DOFINFO &info2)
    {
      Assert(!info1.level_cell, ExcMessage("Cell may not access level dofs"));
      Assert(!info2.level_cell, ExcMessage("Cell may not access level dofs"));
      AssertDimension(info1.indices_by_block.size(),info2.indices_by_block.size());

      const unsigned int n = info1.indices_by_block.size();

      if (n == 0)
        {
          for (unsigned int m=0; m<MatrixSimple<MatrixType>::matrix.size(); ++m)
            {
              assemble(info1.matrix(m,false).matrix, info1.vector(m).block(0), m, info1.indices);
              assemble(info1.matrix(m,true).matrix, info1.vector(m).block(0), m, info1.indices, info2.indices);
              assemble(info2.matrix(m,false).matrix, info2.vector(m).block(0), m, info2.indices);
              assemble(info2.matrix(m,true).matrix, info2.vector(m).block(0), m, info2.indices, info1.indices);
            }
        }
      else
        {
          for (unsigned int m=0; m<MatrixSimple<MatrixType>::matrix.size(); ++m)
            for (unsigned int k=0; k<n*n; ++k)
              {
                const unsigned int row = info1.matrix(k+m*n*n,false).row;
                const unsigned int column = info1.matrix(k+m*n*n,false).column;

                if (row == column)
                  {
                    assemble(info1.matrix(k+m*n*n,false).matrix, info1.vector(m).block(row), m,
                             info1.indices_by_block[row]);
                    assemble(info2.matrix(k+m*n*n,false).matrix, info2.vector(m).block(row), m,
                             info2.indices_by_block[row]);
                  }
                else
                  {
                    assemble(info1.matrix(k+m*n*n,false).matrix, info1.vector(m).block(row), m,
                             info1.indices_by_block[row], info1.indices_by_block[column]);
                    assemble(info2.matrix(k+m*n*n,false).matrix, info2.vector(m).block(row), m,
                             info2.indices_by_block[row], info2.indices_by_block[column]);
                  }
                assemble(info1.matrix(k+m*n*n,true).matrix, info1.vector(m).block(row), m,
                         info1.indices_by_block[row], info2.indices_by_block[column]);
                assemble(info2.matrix(k+m*n*n,true).matrix, info2.vector(m).block(row), m,
                         info2.indices_by_block[row], info1.indices_by_block[column]);
              }
        }
    }
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
