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

#ifndef __deal2__mesh_worker_simple_h
#define __deal2__mesh_worker_simple_h

#include <deal.II/base/named_data.h>
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
			 const unsigned int level = numbers::invalid_unsigned_int);
					 /**
					  * Assemble a single matrix
					  * into a global matrix.
					  */

	void assemble_down(MATRIX& G,
			   const FullMatrix<double>& M,
			   const std::vector<unsigned int>& i1,
			   const std::vector<unsigned int>& i2,
			   const unsigned int level = numbers::invalid_unsigned_int);

					 /**
					  * Assemble a single matrix
					  * into a global matrix.
					  */

	void assemble_in(MATRIX& G,
			 const FullMatrix<double>& M,
			 const std::vector<unsigned int>& i1,
			 const std::vector<unsigned int>& i2,
			 const unsigned int level = numbers::invalid_unsigned_int);

					 /**
					  * Assemble a single matrix
					  * into a global matrix.
					  */

	void assemble_out(MATRIX& G,
			  const FullMatrix<double>& M,
			  const std::vector<unsigned int>& i1,
			  const std::vector<unsigned int>& i2,
			  const unsigned int level = numbers::invalid_unsigned_int);

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
