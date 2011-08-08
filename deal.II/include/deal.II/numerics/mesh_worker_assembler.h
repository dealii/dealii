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

#ifndef __deal2__mesh_worker_assembler_h
#define __deal2__mesh_worker_assembler_h

#include <deal.II/base/named_data.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/numerics/mesh_worker_info.h>
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
 * The class assembling local contributions to a functional into the
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
	void initialize(NamedData<VECTOR*>& results);
					 /**
					  * Initialize the constraints. 
					  */
        void initialize(const ConstraintMatrix& constraints);
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
	void initialize_info(DOFINFO& info, bool face) const;
	
					 /**
					  * Assemble the local residuals
					  * into the global residuals.
					  */
	template <class DOFINFO>
	void assemble(const DOFINFO& info);
	
					 /**
					  * Assemble both local residuals
					  * into the global residuals.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info1,
		      const DOFINFO& info2);
      private:
					 /**
					  * The global residal vectors
					  * filled by assemble().
					  */
	NamedData<SmartPointer<VECTOR,ResidualSimple<VECTOR> > > residuals;
      /**
       * A pointer to the object containing constraints.
       */
      SmartPointer<const ConstraintMatrix,ResidualSimple<VECTOR> > constraints;
    };
    
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
	void initialize(const BlockInfo* block_info,
			NamedData<VECTOR*>& residuals);
					 /**
					  * Initialize the constraints. 
					  */
        void initialize(const ConstraintMatrix& constraints);
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
					  * Assemble the local residuals
					  * into the global residuals.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info);
	
					 /**
					  * Assemble both local residuals
					  * into the global residuals.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info1,
		      const DOFINFO& info2);
      private:
					 /**
					  * Assemble a single local
					  * residual into the global.
					  */
	void assemble(VECTOR& global,
		      const BlockVector<double>& local,
		      const std::vector<unsigned int>& dof);
	
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
 * Assemble local matrices into a single global matrix without using
 * block structure.
 *
 * After being initialized with a SparseMatrix object (or another
 * matrix offering the same functionality as SparseMatrix::add()),
 * this class can be used in a MeshWorker::loop() to assemble the cell
 * and face matrices into the global matrix.
 *
 * @todo On locally refined meshes, a ConstraintMatrix should be used
 * to automatically eliminate hanging nodes.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
    template <class MATRIX>
    class MatrixSimple
    {
      public:
					 /**
					  * Constructor, initializing
					  * the #threshold, which
					  * limits how small numbers
					  * may be to be entered into
					  * the matrix.
					  */
	MatrixSimple(double threshold = 1.e-12);
	
					 /**
					  * Store the result matrix
					  * for later assembling.
					  */
	void initialize(MATRIX& m);
					 /**
					  * Initialize the constraints. 
					  */
        void initialize(const ConstraintMatrix& constraints);
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
					  * Assemble the matrix
					  * DoFInfo::M1[0]
					  * into the global matrix.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info);
	
					 /**
					  * Assemble both local
					  * matrices in the info
					  * objects into the global
					  * matrix.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info1,
		      const DOFINFO& info2);
      private:
					 /**
					  * Assemble a single matrix
					  * into #matrix.
					  */
	void assemble(const FullMatrix<double>& M,
		      const std::vector<unsigned int>& i1,
		      const std::vector<unsigned int>& i2);
	
					 /**
					  * The global matrix being
					  * assembled.
					  */
	SmartPointer<MATRIX,MatrixSimple<MATRIX> > matrix;
      /**
       * A pointer to the object containing constraints.
       */
        SmartPointer<const ConstraintMatrix,MatrixSimple<MATRIX> > constraints;
	
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
					  * Constructor, initializing
					  * the #threshold, which
					  * limits how small numbers
					  * may be to be entered into
					  * the matrix.
					  */
	MGMatrixSimple(double threshold = 1.e-12);
	
					 /**
					  * Store the result matrix
					  * for later assembling.
					  */
	void initialize(MGLevelObject<MATRIX>& m);

					 /**
					  * Initialize the multilevel
                                          * constraints. 
					  */
        void initialize(const MGConstrainedDoFs& mg_constrained_dofs);

					 /**
					  * Initialize the matrices
					  * #flux_up and #flux_down
					  * used for local refinement
					  * with discontinuous
					  * Galerkin methods.
					  */
	void initialize_fluxes(MGLevelObject<MATRIX>& flux_up,
			       MGLevelObject<MATRIX>& flux_down);
	
					 /**
					  * Initialize the matrices
					  * #interface_in and #interface_out
					  * used for local refinement
					  * with continuous
					  * Galerkin methods.
					  */

	void initialize_interfaces(MGLevelObject<MATRIX>& interface_in,
			       MGLevelObject<MATRIX>& interface_out);
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
					  * Assemble the matrix
					  * DoFInfo::M1[0]
					  * into the global matrix.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info);
	
					 /**
					  * Assemble both local
					  * matrices in the info
					  * objects into the global
					  * matrices.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info1,
		      const DOFINFO& info2);
      private:
					 /**
					  * Assemble a single matrix
					  * into a global matrix.
					  */
    void assemble(MATRIX& G,
		  const FullMatrix<double>& M,
		  const std::vector<unsigned int>& i1,
		  const std::vector<unsigned int>& i2);

					 /**
					  * Assemble a single matrix
					  * into a global matrix.
					  */
    void assemble(MATRIX& G,
		  const FullMatrix<double>& M,
		  const std::vector<unsigned int>& i1,
		  const std::vector<unsigned int>& i2,
                  const unsigned int level);
	
					 /**
					  * Assemble a single matrix
					  * into a global matrix.
					  */

    void assemble_up(MATRIX& G,
			    const FullMatrix<double>& M,
			    const std::vector<unsigned int>& i1,
			    const std::vector<unsigned int>& i2,
                            const unsigned int level = -1);
					 /**
					  * Assemble a single matrix
					  * into a global matrix.
					  */

    void assemble_down(MATRIX& G,
			    const FullMatrix<double>& M,
			    const std::vector<unsigned int>& i1,
			    const std::vector<unsigned int>& i2,
                            const unsigned int level = -1);

					 /**
					  * Assemble a single matrix
					  * into a global matrix.
					  */

    void assemble_in(MATRIX& G,
		  const FullMatrix<double>& M,
		  const std::vector<unsigned int>& i1,
		  const std::vector<unsigned int>& i2,
                  const unsigned int level = -1);

					 /**
					  * Assemble a single matrix
					  * into a global matrix.
					  */

    void assemble_out(MATRIX& G,
		  const FullMatrix<double>& M,
		  const std::vector<unsigned int>& i1,
		  const std::vector<unsigned int>& i2,
                  const unsigned int level = -1);
	
					 /**
					  * The global matrix being
					  * assembled.
					  */
	SmartPointer<MGLevelObject<MATRIX>,MGMatrixSimple<MATRIX> > matrix;
	
					 /**
					  * The matrix used for face
					  * flux terms across the
					  * refinement edge, coupling
					  * coarse to fine.
					  */
	SmartPointer<MGLevelObject<MATRIX>,MGMatrixSimple<MATRIX> > flux_up;
	
					 /**
					  * The matrix used for face
					  * flux terms across the
					  * refinement edge, coupling
					  * fine to coarse.
					  */
	SmartPointer<MGLevelObject<MATRIX>,MGMatrixSimple<MATRIX> > flux_down;

					 /**
					  * The matrix used for face
					  * contributions for continuous
                                          * elements across the
					  * refinement edge, coupling
					  * coarse to fine.
					  */
	SmartPointer<MGLevelObject<MATRIX>,MGMatrixSimple<MATRIX> > interface_in;
	
					 /**
					  * The matrix used for face
					  * contributions for continuous
                                          * elements across the
					  * refinement edge, coupling
					  * fine to coarse.
					  */
	SmartPointer<MGLevelObject<MATRIX>,MGMatrixSimple<MATRIX> > interface_out;
      /**
       * A pointer to the object containing constraints.
       */
        SmartPointer<const MGConstrainedDoFs,MGMatrixSimple<MATRIX> > mg_constrained_dofs;
	
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
      void initialize(const BlockInfo* block_info,
		      MatrixBlockVector<MATRIX>& matrices);

					 /**
					  * Initialize the constraints. 
					  */
      void initialize(const ConstraintMatrix& constraints);
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
					  * Assemble the local matrices
					  * into the global matrices.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info);
	
					 /**
					  * Assemble all local matrices
					  * into the global matrices.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info1,
		      const DOFINFO& info2);
	
      private:
					 /**
					  * Assemble a single local
					  * matrix into a global one.
					  */
	void assemble(
	  MatrixBlock<MATRIX>& global,
	  const FullMatrix<number>& local,
	  unsigned int block_row,
	  unsigned int block_col,
	  const std::vector<unsigned int>& dof1,
	  const std::vector<unsigned int>& dof2);
	
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
	void initialize(const BlockInfo* block_info,
			MatrixPtrVector& matrices);

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
	void initialize_edge_flux(MatrixPtrVector& up, MatrixPtrVector& down);
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
					  * Assemble the local matrices
					  * into the global matrices.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info);
	
					 /**
					  * Assemble all local matrices
					  * into the global matrices.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info1,
		      const DOFINFO& info2);
	
      private:
					 /**
					  * Assemble a single local
					  * matrix into a global one.
					  */
	void assemble(
	  MATRIX& global,
	  const FullMatrix<number>& local,
	  unsigned int block_row,
	  unsigned int block_col,
	  const std::vector<unsigned int>& dof1,
	  const std::vector<unsigned int>& dof2,
	  unsigned int level1,
	  unsigned int level2,
	  bool transpose = false);
	
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
	void initialize(MATRIX& m, VECTOR& rhs);
	
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
					  * Assemble the matrix
					  * DoFInfo::M1[0]
					  * into the global matrix.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info);
	
					 /**
					  * Assemble both local
					  * matrices in the info
					  * objects into the global
					  * matrix.
					  */
	template<class DOFINFO>
	void assemble(const DOFINFO& info1,
		      const DOFINFO& info2);
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
    


//----------------------------------------------------------------------//

    template <class VECTOR>
    inline void
    ResidualSimple<VECTOR>::initialize(NamedData<VECTOR*>& results)
    {
      residuals = results;
    }

    template <class VECTOR>
    inline void
    ResidualSimple<VECTOR>::initialize(const ConstraintMatrix& c)
    {
      constraints = &c;
    }

    
    template <class VECTOR>
    template <class DOFINFO>
    inline void
    ResidualSimple<VECTOR>::initialize_info(DOFINFO& info, bool) const
    {
      info.initialize_vectors(residuals.size());
    }

    
    template <class VECTOR>
    template <class DOFINFO>
    inline void
    ResidualSimple<VECTOR>::assemble(const DOFINFO& info)
    {
      for (unsigned int k=0;k<residuals.size();++k)
      {
        if(constraints == 0)
        {
          for (unsigned int i=0;i<info.vector(k).block(0).size();++i)
            (*residuals(k))(info.indices[i]) += info.vector(k).block(0)(i);
        }
        else
          constraints->distribute_local_to_global(
              info.vector(k).block(0), info.indices, (*residuals(k)));
      }
    }

    
    template <class VECTOR>
    template <class DOFINFO>
    inline void
    ResidualSimple<VECTOR>::assemble(const DOFINFO& info1,
					const DOFINFO& info2)
    {
      for (unsigned int k=0;k<residuals.size();++k)
	{
          if(constraints == 0)
          {
            for (unsigned int i=0;i<info1.vector(k).block(0).size();++i)
              (*residuals(k))(info1.indices[i]) += info1.vector(k).block(0)(i);
            for (unsigned int i=0;i<info2.vector(k).block(0).size();++i)
              (*residuals(k))(info2.indices[i]) += info2.vector(k).block(0)(i);
          }
          else
          {
            constraints->distribute_local_to_global(
              info1.vector(k).block(0), info1.indices, (*residuals(k)));
            constraints->distribute_local_to_global(
              info2.vector(k).block(0), info2.indices, (*residuals(k)));
          }
	}
    }

    
//----------------------------------------------------------------------//

    template <class VECTOR>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VECTOR>::initialize(const BlockInfo* b,
							  NamedData<VECTOR*>& m)
    {
      block_info = b;
      residuals = m;      
    }

    template <class VECTOR>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VECTOR>::initialize(
      const ConstraintMatrix& c)
    {
      constraints = &c;
    }


    template <class VECTOR>
    template <class DOFINFO>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VECTOR>::initialize_info(
      DOFINFO& info, bool) const
    {
      info.initialize_vectors(residuals.size());
    }

    template <class VECTOR>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VECTOR>::assemble(
      VECTOR& global,
      const BlockVector<double>& local,
      const std::vector<unsigned int>& dof)
    {
        if(constraints == 0)
        {
          for (unsigned int b=0;b<local.n_blocks();++b)
            for (unsigned int j=0;j<local.block(b).size();++j)
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
      const DOFINFO& info)
    {
      for (unsigned int i=0;i<residuals.size();++i)
	assemble(*residuals(i), info.vector(i), info.indices);
    }

    
    template <class VECTOR>
    template <class DOFINFO>
    inline void
    ResidualLocalBlocksToGlobalBlocks<VECTOR>::assemble(
      const DOFINFO& info1,
      const DOFINFO& info2)
    {
      for (unsigned int i=0;i<residuals.size();++i)
	{
	  assemble(*residuals(i), info1.vector(i), info1.indices);
	  assemble(*residuals(i), info2.vector(i), info2.indices);
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
    MatrixSimple<MATRIX>::initialize(MATRIX& m)
    {
      matrix = &m;
    }


    template <class MATRIX>
    inline void
    MatrixSimple<MATRIX>::initialize(const ConstraintMatrix& c)
    {
      constraints = &c;
    }
    

    template <class MATRIX >
    template <class DOFINFO>
    inline void
    MatrixSimple<MATRIX>::initialize_info(DOFINFO& info, bool face) const
    {
      info.initialize_matrices(1, face);
    }



    template <class MATRIX>
    inline void
    MatrixSimple<MATRIX>::assemble(const FullMatrix<double>& M,
				   const std::vector<unsigned int>& i1,
				   const std::vector<unsigned int>& i2)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());
     
     if(constraints == 0)
     { 
      for (unsigned int j=0; j<i1.size(); ++j)
	for (unsigned int k=0; k<i2.size(); ++k)
	  if (std::fabs(M(j,k)) >= threshold)
	    matrix->add(i1[j], i2[k], M(j,k));
     }
     else
       constraints->distribute_local_to_global(
           M, i1, i2, *matrix);
    }
    
    
    template <class MATRIX>
    template <class DOFINFO>
    inline void
    MatrixSimple<MATRIX>::assemble(const DOFINFO& info)
    {
      assemble(info.matrix(0,false).matrix, info.indices, info.indices);
    }
    

    template <class MATRIX>
    template <class DOFINFO>
    inline void
    MatrixSimple<MATRIX>::assemble(const DOFINFO& info1,
				   const DOFINFO& info2)
    {
      assemble(info1.matrix(0,false).matrix, info1.indices, info1.indices);
      assemble(info1.matrix(0,true).matrix, info1.indices, info2.indices);
      assemble(info2.matrix(0,false).matrix, info2.indices, info2.indices);
      assemble(info2.matrix(0,true).matrix, info2.indices, info1.indices);
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
    MGMatrixSimple<MATRIX>::initialize(MGLevelObject<MATRIX>& m)
    {
      matrix = &m;
    }
    
    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::initialize(const MGConstrainedDoFs& c)
    {
      mg_constrained_dofs = &c;
    }

    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::initialize_fluxes(
      MGLevelObject<MATRIX>& up, MGLevelObject<MATRIX>& down)
    {
      flux_up = &up;
      flux_down = &down;
    }
    

    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::initialize_interfaces(
      MGLevelObject<MATRIX>& in, MGLevelObject<MATRIX>& out)
    {
      interface_in = &in;
      interface_out = &out;
    }


    template <class MATRIX >
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MATRIX>::initialize_info(DOFINFO& info, bool face) const
    {
      info.initialize_matrices(1, face);
    }


    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble(
      MATRIX& G,
      const FullMatrix<double>& M,
      const std::vector<unsigned int>& i1,
      const std::vector<unsigned int>& i2)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());
      
      if(mg_constrained_dofs == 0)
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
            {
                if(!mg_constrained_dofs->continuity_across_refinement_edges())
                  G.add(i1[j], i2[k], M(j,k));
            }
      }
    }


    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble(
      MATRIX& G,
      const FullMatrix<double>& M,
      const std::vector<unsigned int>& i1,
      const std::vector<unsigned int>& i2,
      const unsigned int level)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());
      
      if(mg_constrained_dofs == 0)
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
              if (!mg_constrained_dofs->at_refinement_edge(level, i1[j]) &&
                  !mg_constrained_dofs->at_refinement_edge(level, i2[k]))
              {
                if (mg_constrained_dofs->set_boundary_values())
                {
                  if ((!mg_constrained_dofs->is_boundary_index(level, i1[j]) &&
                        !mg_constrained_dofs->is_boundary_index(level, i2[k]))
                      ||
                      (mg_constrained_dofs->is_boundary_index(level, i1[j]) &&
                       mg_constrained_dofs->is_boundary_index(level, i2[k]) &&
                       i1[j] == i2[k]))
                    G.add(i1[j], i2[k], M(j,k));
                }
                else
                    G.add(i1[j], i2[k], M(j,k));
            }
      }
    }
    
    
    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble_up(
      MATRIX& G,
      const FullMatrix<double>& M,
      const std::vector<unsigned int>& i1,
      const std::vector<unsigned int>& i2,
      const unsigned int level)
    {
      AssertDimension(M.n(), i1.size());
      AssertDimension(M.m(), i2.size());
      
      if(mg_constrained_dofs == 0)
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
              if(mg_constrained_dofs->at_refinement_edge(level, i1[j]) &&
                  !mg_constrained_dofs->at_refinement_edge(level, i2[k]))
                G.add(i1[j], i2[k], M(k,j));
      }
    }

    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble_down(
      MATRIX& G,
      const FullMatrix<double>& M,
      const std::vector<unsigned int>& i1,
      const std::vector<unsigned int>& i2,
      const unsigned int level)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());
      
      if(mg_constrained_dofs == 0)
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
              if(mg_constrained_dofs->at_refinement_edge(level, i1[j]) &&
                  !mg_constrained_dofs->at_refinement_edge(level, i2[k]))
                  G.add(i1[j], i2[k], M(j,k));
      }
    }

    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble_in(
      MATRIX& G,
      const FullMatrix<double>& M,
      const std::vector<unsigned int>& i1,
      const std::vector<unsigned int>& i2,
      const unsigned int level)
    {
      AssertDimension(M.m(), i1.size());
      AssertDimension(M.n(), i2.size());
      
      if(mg_constrained_dofs == 0)
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
              if(mg_constrained_dofs->at_refinement_edge(level, i1[j]) &&
                  !mg_constrained_dofs->at_refinement_edge(level, i2[k]))
              {
                if (mg_constrained_dofs->set_boundary_values())
                {
                  if((!mg_constrained_dofs->at_refinement_edge_boundary(level, i1[j]) &&
                        !mg_constrained_dofs->at_refinement_edge_boundary(level, i2[k]))
                      ||
                      (mg_constrained_dofs->at_refinement_edge_boundary(level, i1[j]) &&
                       mg_constrained_dofs->at_refinement_edge_boundary(level, i2[k]) &&
                       i1[j] == i2[k]))
                    G.add(i1[j], i2[k], M(j,k));
                }
                else
                  G.add(i1[j], i2[k], M(j,k));
              }
      }
    }
    
    template <class MATRIX>
    inline void
    MGMatrixSimple<MATRIX>::assemble_out(
      MATRIX& G,
      const FullMatrix<double>& M,
      const std::vector<unsigned int>& i1,
      const std::vector<unsigned int>& i2,
      const unsigned int level)
    {
      AssertDimension(M.n(), i1.size());
      AssertDimension(M.m(), i2.size());
      
      if(mg_constrained_dofs == 0)
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
            if (std::fabs(M(j,k)) >= threshold)
              if(mg_constrained_dofs->at_refinement_edge(level, i1[j]) &&
                  !mg_constrained_dofs->at_refinement_edge(level, i2[k]))
              {
                if (mg_constrained_dofs->set_boundary_values())
                {
                  if((!mg_constrained_dofs->at_refinement_edge_boundary(level, i1[j]) &&
                        !mg_constrained_dofs->at_refinement_edge_boundary(level, i2[k])) 
                      ||
                      (mg_constrained_dofs->at_refinement_edge_boundary(level, i1[j]) &&
                       mg_constrained_dofs->at_refinement_edge_boundary(level, i2[k]) &&
                       i1[j] == i2[k]))
                    G.add(i1[j], i2[k], M(k,j));
                }
                else
                    G.add(i1[j], i2[k], M(k,j));
              }
      }
    }
    
    
    template <class MATRIX>
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MATRIX>::assemble(const DOFINFO& info)
    {
      const unsigned int level = info.cell->level();
      assemble((*matrix)[level], info.matrix(0,false).matrix, info.indices, info.indices, level);

      if(mg_constrained_dofs != 0)
      {
        assemble_in((*interface_in)[level], info.matrix(0,false).matrix, info.indices, info.indices, level);
        assemble_out((*interface_out)[level],info.matrix(0,false).matrix, info.indices, info.indices, level);
      }
    }
    

    template <class MATRIX>
    template <class DOFINFO>
    inline void
    MGMatrixSimple<MATRIX>::assemble(const DOFINFO& info1,
				     const DOFINFO& info2)
    {
      const unsigned int level1 = info1.cell->level();
      const unsigned int level2 = info2.cell->level();
      
      if (level1 == level2)
	{
          if(mg_constrained_dofs == 0)
          {
            assemble((*matrix)[level1], info1.matrix(0,false).matrix, info1.indices, info1.indices);
            assemble((*matrix)[level1], info1.matrix(0,true).matrix, info1.indices, info2.indices);
            assemble((*matrix)[level1], info2.matrix(0,false).matrix, info2.indices, info2.indices);
            assemble((*matrix)[level1], info2.matrix(0,true).matrix, info2.indices, info1.indices);
          }
          else
          {
            assemble((*matrix)[level1], info1.matrix(0,false).matrix, info1.indices, info1.indices, level1);
            assemble((*matrix)[level1], info1.matrix(0,true).matrix, info1.indices, info2.indices, level1);
            assemble((*matrix)[level1], info2.matrix(0,false).matrix, info2.indices, info2.indices, level1);
            assemble((*matrix)[level1], info2.matrix(0,true).matrix, info2.indices, info1.indices, level1);
          }
	}
      else
	{
	  Assert(level1 > level2, ExcInternalError());
					   // Do not add info2.M1,
					   // which is done by
					   // the coarser cell
          assemble((*matrix)[level1], info1.matrix(0,false).matrix, info1.indices, info1.indices);
          //assemble_transpose((*flux_up)[level1],info1.matrix(0,true).matrix, info2.indices, info1.indices);
          //assemble((*flux_down)[level1], info2.matrix(0,true).matrix, info2.indices, info1.indices);
          if(level1>0)
          {
          assemble_up((*flux_up)[level1],info1.matrix(0,true).matrix, info2.indices, info1.indices, level1);
          assemble_down((*flux_down)[level1], info2.matrix(0,true).matrix, info2.indices, info1.indices, level1);
          }
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
      const BlockInfo* b,
      MatrixBlockVector<MATRIX>& m)
    {
      block_info = b;
      matrices = &m;
    }



    template <class MATRIX, typename number>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize(
      const ConstraintMatrix& c)
    {
      constraints = &c;
    }
    

    
    template <class MATRIX ,typename number>
    template <class DOFINFO>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize_info(
      DOFINFO& info,
      bool face) const
    {
      info.initialize_matrices(*matrices, face);
    }



    template <class MATRIX, typename number>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(
      MatrixBlock<MATRIX>& global,
      const FullMatrix<number>& local,
      unsigned int block_row,
      unsigned int block_col,
      const std::vector<unsigned int>& dof1,
      const std::vector<unsigned int>& dof2)
    {
      if(constraints == 0)
      {
        for (unsigned int j=0;j<local.n_rows();++j)
          for (unsigned int k=0;k<local.n_cols();++k)
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
        std::vector<unsigned int> sliced_row_indices (bi.block_size(block_row));
        for(unsigned int i=0; i<sliced_row_indices.size(); ++i)
          sliced_row_indices[i] = dof1[bi.block_start(block_row)+i];

        std::vector<unsigned int> sliced_col_indices (bi.block_size(block_col));
        for(unsigned int i=0; i<sliced_col_indices.size(); ++i)
          sliced_col_indices[i] = dof2[bi.block_start(block_col)+i];
        
        constraints->distribute_local_to_global(local,
            sliced_row_indices, sliced_col_indices, global);
      }
    }

    
    template <class MATRIX, typename number>
    template <class DOFINFO>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(
      const DOFINFO& info)
    {
      for (unsigned int i=0;i<matrices->size();++i)
	{
					   // Row and column index of
					   // the block we are dealing with
	  const unsigned int row = matrices->block(i).row;
	  const unsigned int col = matrices->block(i).column;

	  assemble(matrices->block(i), info.matrix(i,false).matrix, row, col, info.indices, info.indices);
	}
    }

    
    template <class MATRIX, typename number>
    template <class DOFINFO>
    inline void
    MatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(
      const DOFINFO& info1,
      const DOFINFO& info2)
    {
      for (unsigned int i=0;i<matrices->size();++i)
	{
					   // Row and column index of
					   // the block we are dealing with
	  const unsigned int row = matrices->block(i).row;
	  const unsigned int col = matrices->block(i).column;

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
      const BlockInfo* b,
      MatrixPtrVector& m)
    {
      block_info = b;
      AssertDimension(block_info->local().size(), block_info->global().size());
      matrices = &m;
    }
    

    template <class MATRIX ,typename number>
    template <class DOFINFO>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize_info(
      DOFINFO& info,
      bool face) const
    {
      info.initialize_matrices(*matrices, face);
    }



    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::initialize_edge_flux(
      MatrixPtrVector& up,
      MatrixPtrVector& down)
    {
      flux_up = up;
      flux_down = down;
    }
    
    
    template <class MATRIX, typename number>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(
      MATRIX& global,
      const FullMatrix<number>& local,
      unsigned int block_row,
      unsigned int block_col,
      const std::vector<unsigned int>& dof1,
      const std::vector<unsigned int>& dof2,
      unsigned int level1,
      unsigned int level2,
      bool transpose)
    {
      for (unsigned int j=0;j<local.n_rows();++j)
	for (unsigned int k=0;k<local.n_cols();++k)
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

	      if (transpose)
		global.add(kglobal, jglobal, local(j,k));
	      else
		global.add(jglobal, kglobal, local(j,k));
	    }
    }

    
    template <class MATRIX, typename number>
    template <class DOFINFO>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(const DOFINFO& info)
    {
      const unsigned int level = info.cell->level();
      
      for (unsigned int i=0;i<matrices->size();++i)
	{
					   // Row and column index of
					   // the block we are dealing with
	  const unsigned int row = matrices->block(i)[level].row;
	  const unsigned int col = matrices->block(i)[level].column;

	  assemble(matrices->block(i)[level].matrix, info.matrix(i,false).matrix, row, col,
		   info.indices, info.indices, level, level);
	}
    }

    
    template <class MATRIX, typename number>
    template <class DOFINFO>
    inline void
    MGMatrixLocalBlocksToGlobalBlocks<MATRIX, number>::assemble(
      const DOFINFO& info1,
      const DOFINFO& info2)
    {
      const unsigned int level1 = info1.cell->level();
      const unsigned int level2 = info2.cell->level();
      
      for (unsigned int i=0;i<matrices->size();++i)
	{
	  MGLevelObject<MatrixBlock<MATRIX> >& o = matrices->block(i);
	  
					   // Row and column index of
					   // the block we are dealing with
	  const unsigned int row = o[level1].row;
	  const unsigned int col = o[level1].column;

	  if (level1 == level2)
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
	      Assert(level1 > level2, ExcNotImplemented());
	      if (flux_up->size() != 0)
		{
						   // Do not add M22,
						   // which is done by
						   // the coarser cell
		  assemble(o[level1].matrix, info1.matrix(i,false).matrix, row, col,
			   info1.indices, info1.indices, level1, level1);	      
		  assemble(flux_up->block(i)[level1].matrix, info1.matrix(i,true).matrix, row, col,
			   info1.indices, info2.indices, level1, level2, true);
		  assemble(flux_down->block(i)[level1].matrix, info2.matrix(i,true).matrix, row, col,
			   info2.indices, info1.indices, level2, level1);
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
    SystemSimple<MATRIX,VECTOR>::initialize(MATRIX& m, VECTOR& rhs)
    {
      NamedData<VECTOR*> data;
      VECTOR* p = &rhs;
      data.add(p, "right hand side");
      
      MatrixSimple<MATRIX>::initialize(m);
      ResidualSimple<VECTOR>::initialize(data);
    }

    
    template <class MATRIX, class VECTOR>
    template <class DOFINFO>
    inline void
    SystemSimple<MATRIX,VECTOR>::initialize_info(DOFINFO& info,
						 bool face) const
    {
      MatrixSimple<MATRIX>::initialize_info(info, face);
      ResidualSimple<VECTOR>::initialize_info(info, face);
    }


    template <class MATRIX, class VECTOR>
    template <class DOFINFO>
    inline void
    SystemSimple<MATRIX,VECTOR>::assemble(const DOFINFO& info)
    {
      MatrixSimple<MATRIX>::assemble(info);
      ResidualSimple<VECTOR>::assemble(info);
    }


    template <class MATRIX, class VECTOR>
    template <class DOFINFO>
    inline void
    SystemSimple<MATRIX,VECTOR>::assemble(const DOFINFO& info1,
					  const DOFINFO& info2)
    {
      MatrixSimple<MATRIX>::assemble(info1, info2);
      ResidualSimple<VECTOR>::assemble(info1, info2);
    }    
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
