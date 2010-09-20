//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#ifndef __deal2__mesh_worker_h
#define __deal2__mesh_worker_h

#include <base/config.h>
#include <base/std_cxx1x/function.h>
#include <base/geometry_info.h>
#include <lac/matrix_block.h>
#include <lac/block_vector.h>
#include <numerics/mesh_worker_vector_selector.h>

DEAL_II_NAMESPACE_OPEN

class BlockIndices;
template<int,int> class DoFHandler;
template<int,int> class MGDoFHandler;

/**
 * A collection of functions and classes for the mesh loops that are
 * an ubiquitous part of each finite element program.
 *
 * The workhorse of this namespace is the loop() function, which implements a
 * completely generic loop over all mesh cells.
 *
 * The loop() depends on certain objects handed to it as
 * arguments. These objects are of two types, info objects like
 * DoFInfo and IntegrationInfo and worker objects like LocalWorker and
 * IntegrationWorker.
 *
 * Worker objects usually do two different jobs: first, they compute
 * the local contribution of a cell or face to the global
 * operation. Second, they assemble this local contribution into the
 * global result, whether a functional, a form or a bilinear
 * form. While the first job is particular to the problem being
 * solved, the second is generic and only depends on the data
 * structures. Therefore, base classes for workers assembling into
 * global data are provided in the namespace Assembler.
 *
 * <h3>Template argument types</h3>
 *
 * The functions loop() and cell_action() take some arguments which
 * are template parameters. Let us list the minimum requirements for
 * these classes here and describe their properties.
 *
 * <h4>ITERATOR</h4>
 *
 * Any object that has an <tt>operator++()</tt> and points to a
 * TriaObjectAccessor.
 *
 * <h4>DOFINFO</h4>
 *
 * For an example implementation, refer to the class template DoFInfo.
 * In order to work with cell_action() and loop(), DOFINFO needs to
 * follow the following interface.
 * @code
 * class DOFINFO
 * {
 *   private:
 *     DOFINFO();
 *     DOFINFO(const DOFINFO&);
 *     DOFINFO& operator=(const DOFINFO&);
 *
 *   public:
 *     template <class CellIt>
 *     void reinit(const CellIt& c);
 *
 *     template <class CellIt, class FaceIt>
 *     void reinit(const CellIt& c, const FaceIt& f, unsigned int n);
 *
 *     template <class CellIt, class FaceIt>
 *     void reinit(const CellIt& c, const FaceIt& f, unsigned int n,
 *     unsigned int s);
 *
 *   friend template class DoFInfoBox<int dim, DOFINFO>;
 * };
 * @endcode
 *
 * The three private functions are called by DoFInfoBox and should not
 * be needed elsewhere. Obviously, they can be made public and then
 * the friend declaration at the end may be missing.
 *
 * Additionally, you will need at least one public constructor. Furthermore
 * DOFINFO is pretty useless yet: functions to interface with
 * INTEGRATIONINFO and ASSEMBLER are needed.
 *
 * DOFINFO objects are gathered in a DoFInfoBox. In those objects, we
 * store the results of local operations on each cel and its
 * faces. Once all this information has been gathered, an ASSEMBLER is
 * used to assemble it into golbal data.
 *
 * <h4>INFOBOX</h4>
 *
 * This type is exemplified in IntegrationInfoBox. It collects the
 * input data for actions on cells and faces in INFO objects (see
 * below). It provides the following interface to loop() and
 * cell_action():
 *
 * @code
 * class INFOBOX
 * {
 *   public:
 *     template <int dim, class DOFINFO>
 *     void post_cell(const DoFInfoBox<dim, DOFINFO>&);
 *
 *     template <int dim, class DOFINFO>
 *     void post_faces(const DoFInfoBox<dim, DOFINFO>&);
 *
 *     INFO cell;
 *     INFO boundary;
 *     INFO face;
 *     INFO subface;
 *     INFO neighbor;
 * };
 * @endcode
 *
 * The main purpose of this class is gathering the five INFO objects,
 * which contain the temporary data used on each cell or face. The
 * requirements on these objects are listed below. Here, we only note
 * that there need to be these 5 objects with the names listed above.
 *
 * The two function templates are call back functions called in
 * cell_action(). The first is called before the faces are worked on,
 * the second after the faces.
 *
 * <h4>INFO</h4>
 *
 * See IntegrationInfo for an example of these objects. They contain
 * the temporary data needed on each cell or face to compute the
 * result. The MeshWorker only uses the interface
 *
 * @code
 * class INFO
 * {
 *   public:
 *     void reinit(const DOFINFO& i);
 * };
 * @endcode
 *
 * <h3>Simplified interfaces</h3>
 *
 * Since the loop() is fairly general, a specialization
 * integration_loop() is available, which is a wrapper around loop()
 * with a simplified interface.
 *
 * The integration_loop() function loop takes most of the information
 * that it needs to pass to loop() from an IntegrationInfoBox
 * object. Its use is explained in step-12, but in
 * short it requires functions that do the local integration on a
 * cell, interior or boundary face, and it needs an object (called
 * "assembler") that copies these local contributions into the global
 * matrix and right hand side objects.
 *
 * Before we can run the integration loop, we have to initialize
 * several data structures in our IntegrationWorker and assembler
 * objects. For instance, we have to decide on the quadrature rule or
 * we may need more than the default update flags.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
namespace MeshWorker
{
/**
 * The class providing the scrapbook to fill with local integration
 * results. These can be values, local contributions to forms or cell
 * and face matrices.
 *
 * The local matrices initialized by reinit() of the info object and
 * then assembled into the global system by Assembler classes.
 *
 * @ingroup MeshWorker
 * @author Guido Kanschat, 2009
 */
  template <typename number>
  class LocalResults
  {
    private:
				       /**
					* Initialize a single local
					* matrix block. A helper
					* function for initialize()
					*/
      void initialize_local(MatrixBlock<FullMatrix<number> >& M,
			    const unsigned int row,
			    const unsigned int col);

    public:
      void initialize_numbers(const unsigned int n);
      void initialize_vectors(const unsigned int n);
				       /**
					* Allocate @p n local
					* matrices. Additionally,
					* set their block row and
					* column coordinates to
					* zero. The matrices
					* themselves are resized by
					* reinit().
					*
					* The template parameter @p
					* MatrixPtr should point to
					* a MatrixBlock
					* instantiation in order to
					* provide row and column info.
					*/
      void initialize_matrices(unsigned int n, bool both);

				       /**
					* Allocate a local matrix
					* for each of the global
					* ones in @p
					* matrices. Additionally,
					* set their block row and
					* column coordinates. The
					* matrices themselves are
					* resized by reinit().
					*
					* The template parameter @p
					* MatrixPtr should point to
					* a MatrixBlock
					* instantiation in order to
					* provide row and column info.
					*/
      template <class MATRIX>
      void initialize_matrices(const MatrixBlockVector<MATRIX>& matrices,
			       bool both);

				       /**
					* Initialize quadrature values
					* to <tt>nv</tt> values in
					* <tt>np</tt> quadrature points.
					*/
      void initialize_quadrature(unsigned int np, unsigned int nv);
      
				       /**
					* Reinitialize matrices for
					* new cell. Resizes the
					* matrices for hp and sets
					* them to zero.
					*/
      void reinit(const BlockIndices& local_sizes);

				       /**
					* The number of scalar values.
					*/
      unsigned int n_values () const;
      
				       /**
					* The number of vectors.
					*/
      unsigned int n_vectors () const;
      
				       /**
					* The number of matrices.
					*/
      unsigned int n_matrices () const;
      
				       /**
					* The number of quadrature
					* points in #quadrature_values.
					*/
      unsigned int n_quadrature_points() const;
      
				       /**
					* The number of values in each
					* quadrature point in
					* #quadrature_values.
					*/
      unsigned int n_quadrature_values() const;
      
				       /**
					* Access scalar value at index
					* @p i.
					*/
      number& value(unsigned int i);

				       /**
					* Read scalar value at index
					* @p i.
					*/
      number value(unsigned int i) const;

				       /**
					* Access vector at index @p i.
					*/
      BlockVector<number>& vector(unsigned int i);
      
				       /**
					* Read vector at index @p i.
					*/
      const BlockVector<number>& vector(unsigned int i) const;
      
				       /**
					* Access matrix at index @p
					* i. For results on internal
					* faces, a true value for @p
					* external refers to the flux
					* between cells, while false
					* refers to entries coupling
					* inside the cell.
					*/
      MatrixBlock<FullMatrix<number> >& matrix(unsigned int i, bool external = false);
      
				       /**
					* Read matrix at index @p
					* i. For results on internal
					* faces, a true value for @p
					* external refers to the flux
					* between cells, while false
					* refers to entries coupling
					* inside the cell.
					*/
      const MatrixBlock<FullMatrix<number> >& matrix(unsigned int i, bool external = false) const;
      
				       /**
					* Access the <i>i</i>th value
					* at quadrature point <i>k</i>
					*/
      number& quadrature_value(unsigned int k, unsigned int i);
      
				       /**
					* Read the <i>i</i>th value
					* at quadrature point <i>k</i>
					*/
      number quadrature_value(unsigned int k, unsigned int i) const;

				       /**
					* The memory used by this object.
					*/
      unsigned int memory_consumption () const;
      
    private:
				       /**
					* The local numbers,
					* computed on a cell or on a
					* face.
					*/
      std::vector<number> J;

				       /**
					* The local vectors. This
					* field is public, so that
					* local integrators can
					* write to it.
					*/
      std::vector<BlockVector<number> > R;

				       /**
					* The local matrices
					* coupling degrees of
					* freedom in the cell
					* itself or within the
					* first cell on a face.
					*/
      std::vector<MatrixBlock<FullMatrix<number> > > M1;

				       /**
					* The local matrices
					* coupling test functions on
					* the cell with trial
					* functions on the other
					* cell.
					*
					* Only used on interior
					* faces.
					*/
      std::vector<MatrixBlock<FullMatrix<number> > > M2;

				       /**
					* Values in quadrature points.
					*/
      std::vector<std::vector<number> > quadrature_values;
  };

//----------------------------------------------------------------------//

  template <typename number>
  inline void
  LocalResults<number>::initialize_numbers(unsigned int n)
  {
    J.resize(n);
  }


  template <typename number>
  inline void
  LocalResults<number>::initialize_vectors(const unsigned int n)
  {
    R.resize(n);
  }


  template <typename number>
  template <class MATRIX>
  inline void
  LocalResults<number>::initialize_matrices(
    const MatrixBlockVector<MATRIX>& matrices,
    bool both)
  {
    M1.resize(matrices.size());
    if (both)
      M2.resize(matrices.size());
    for (unsigned int i=0;i<matrices.size();++i)
      {
	const unsigned int row = matrices.block(i).row;
	const unsigned int col = matrices.block(i).column;

	M1[i].row = row;
	M1[i].column = col;
	if (both)
	  {
	    M2[i].row = row;
	    M2[i].column = col;
	  }
      }
  }


  template <typename number>
  inline void
  LocalResults<number>::initialize_matrices(const unsigned int n,
					    const bool both)
  {
    M1.resize(n);
    if (both)
      M2.resize(n);
    for (unsigned int i=0;i<n;++i)
      {
	M1[i].row = 0;
	M1[i].column = 0;
	if (both)
	  {
	    M2[i].row = 0;
	    M2[i].column = 0;
	  }
      }
  }
  

  template <typename number>
  inline void
  LocalResults<number>::initialize_quadrature(unsigned int np, unsigned int nv)
  {
    quadrature_values.resize(np, std::vector<number>(nv));
  }
  
  
  template <typename number>
  inline
  unsigned int
  LocalResults<number>::n_values() const
  {
    return J.size();
  }
  

  template <typename number>
  inline
  unsigned int
  LocalResults<number>::n_vectors() const
  {
    return R.size();
  }
  

  template <typename number>
  inline
  unsigned int
  LocalResults<number>::n_matrices() const
  {
    return M1.size();
  }
  

  template <typename number>
  inline
  unsigned int
  LocalResults<number>::n_quadrature_points() const
  {
    return quadrature_values.size();
  }

  
  template <typename number>
  inline
  unsigned int
  LocalResults<number>::n_quadrature_values() const
  {
    Assert(quadrature_values.size() != 0, ExcNotInitialized());
    return quadrature_values[0].size();
  }
  
  
  template <typename number>
  inline
  number&
  LocalResults<number>::value(unsigned int i)
  {
    AssertIndexRange(i,J.size());
    return J[i];
  }
  

  template <typename number>
  inline
  BlockVector<number>&
  LocalResults<number>::vector(unsigned int i)
  {
    AssertIndexRange(i,R.size());
    return R[i];
  }
  
  
  template <typename number>
  inline
  MatrixBlock<FullMatrix<number> >&
  LocalResults<number>::matrix(unsigned int i, bool external)
  {
    if (external)
      {
	AssertIndexRange(i,M2.size());
	return M2[i];
      }
    AssertIndexRange(i,M1.size());
    return M1[i];
  }
  
  
  template <typename number>
  inline
  number&
  LocalResults<number>::quadrature_value(unsigned int k, unsigned int i)
  {
    AssertIndexRange(k,quadrature_values.size());
    AssertIndexRange(i,quadrature_values[0].size());
    return quadrature_values[k][i];
  }
  
  
  template <typename number>
  inline
  number
  LocalResults<number>::value(unsigned int i) const
  {
    AssertIndexRange(i,J.size());
    return J[i];
  }
  

  template <typename number>
  inline
  const BlockVector<number>&
  LocalResults<number>::vector(unsigned int i) const
  {
    AssertIndexRange(i,R.size());
    return R[i];
  }
  
  
  template <typename number>
  inline
  const MatrixBlock<FullMatrix<number> >&
  LocalResults<number>::matrix(unsigned int i, bool external) const
  {
    if (external)
      {
	AssertIndexRange(i,M2.size());
	return M2[i];
      }
    AssertIndexRange(i,M1.size());
    return M1[i];
  }
  
  
  template <typename number>
  inline
  number
  LocalResults<number>::quadrature_value(unsigned int k, unsigned int i) const
  {
    AssertIndexRange(k,quadrature_values.size());
    AssertIndexRange(i,quadrature_values[0].size());
    return quadrature_values[k][i];
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
