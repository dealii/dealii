//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__trilinos_block_vector_h
#define __deal2__trilinos_block_vector_h


#include <base/config.h>

#ifdef DEAL_II_USE_TRILINOS

#  include <lac/trilinos_vector.h>
#  include <lac/block_indices.h>
#  include <lac/block_vector_base.h>
#  include <lac/exceptions.h>

DEAL_II_NAMESPACE_OPEN

                                   // forward declaration
template <typename Number> class BlockVector;

/*! @addtogroup TrilinosWrappers
 *@{
 */

namespace TrilinosWrappers
{
                                   // forward declaration
  namespace MPI
  {
    class BlockVector;
  }
  class BlockVector;
  class BlockSparseMatrix;


  namespace MPI
  {
/**
 * An implementation of block vectors based on the vector class
 * implemented in TrilinosWrappers. While the base class provides for
 * most of the interface, this class handles the actual allocation of
 * vectors and provides functions that are specific to the underlying
 * vector type.
 *
 * The model of distribution of data is such that each of the blocks
 * is distributed across all MPI processes named in the MPI
 * communicator. I.e. we don't just distribute the whole vector, but
 * each component. In the constructors and reinit() functions, one
 * therefore not only has to specify the sizes of the individual
 * blocks, but also the number of elements of each of these blocks to
 * be stored on the local process.
 *
 * @ingroup Vectors
 * @ingroup TrilinosWrappers
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Martin Kronbichler, Wolfgang Bangerth, 2008, 2009
 */
    class BlockVector : public BlockVectorBase<Vector>
    {
      public:
                                       /**
					* Typedef the base class for simpler
					* access to its own typedefs.
					*/
        typedef BlockVectorBase<Vector> BaseClass;

                                       /**
					* Typedef the type of the underlying
					* vector.
					*/
        typedef BaseClass::BlockType  BlockType;

                                       /**
					* Import the typedefs from the base
					* class.
					*/
        typedef BaseClass::value_type      value_type;
        typedef BaseClass::pointer         pointer;
        typedef BaseClass::const_pointer   const_pointer;
        typedef BaseClass::reference       reference;
        typedef BaseClass::const_reference const_reference;
        typedef BaseClass::size_type       size_type;
        typedef BaseClass::iterator        iterator;
        typedef BaseClass::const_iterator  const_iterator;

                                       /**
					* Default constructor. Generate an
					* empty vector without any blocks.
					*/
        BlockVector ();

                                       /**
					* Constructor. Generate a block
					* vector with as many blocks as
					* there are entries in @p
					* partitioning. Each Epetra_Map
					* contains the layout of the
					* distribution of data among the MPI
					* processes.
					*/
        BlockVector (const std::vector<Epetra_Map> &parallel_partitioning);

                                       /**
					* Constructor. Generate a block
					* vector with as many blocks as
					* there are entries in
					* @p partitioning.  Each IndexSet
					* together with the MPI communicator
					* contains the layout of the
					* distribution of data among the MPI
					* processes.
					*/
        BlockVector (const std::vector<IndexSet> &parallel_partitioning,
		     const MPI_Comm              &communicator = MPI_COMM_WORLD);

                                       /**
					* Copy-Constructor. Set all the
					* properties of the parallel vector
					* to those of the given argument and
					* copy the elements.
					*/
        BlockVector (const BlockVector  &V);

                                       /**
					* Creates a block vector
					* consisting of
					* <tt>num_blocks</tt>
					* components, but there is no
					* content in the individual
					* components and the user has to
					* fill appropriate data using a
					* reinit of the blocks.
					*/
        BlockVector (const unsigned int num_blocks);

                                       /**
					* Destructor. Clears memory
					*/
        ~BlockVector ();

                                       /**
					* Copy operator: fill all
					* components of the vector that
					* are locally stored with the
					* given scalar value.
					*/
        BlockVector &
	  operator = (const value_type s);

                                       /**
					* Copy operator for arguments of
					* the same type.
					*/
        BlockVector &
	  operator = (const BlockVector &V);

                                       /**
					* Copy operator for arguments of
					* the localized Trilinos vector
					* type.
					*/
        BlockVector &
	  operator = (const ::dealii::TrilinosWrappers::BlockVector &V);

                                       /**
					* Another copy function. This
					* one takes a deal.II block
					* vector and copies it into a
					* TrilinosWrappers block
					* vector. Note that the number
					* of blocks has to be the same
					* in the vector as in the input
					* vector. Use the reinit()
					* command for resizing the
					* BlockVector or for changing
					* the internal structure of the
					* block components.
					*
					* Since Trilinos only works on
					* doubles, this function is
					* limited to accept only one
					* possible number type in the
					* deal.II vector.
					*/
        template <typename Number>
	BlockVector &
	  operator = (const ::dealii::BlockVector<Number> &V);

                                         /**
                                          * Reinitialize the BlockVector to
                                          * contain as many blocks as there
					  * are Epetra_Maps given in the input
					  * argument, according to the
					  * parallel distribution of the
					  * individual components described
					  * in the maps.
                                          *
                                          * If <tt>fast==false</tt>, the vector
                                          * is filled with zeros.
                                          */
        void reinit (const std::vector<Epetra_Map> &parallel_partitioning,
		     const bool                     fast = false);

                                         /**
                                          * Reinitialize the BlockVector to
                                          * contain as many blocks as there
					  * are index sets given in the input
					  * argument, according to the
					  * parallel distribution of the
					  * individual components described
					  * in the maps.
                                          *
                                          * If <tt>fast==false</tt>, the vector
                                          * is filled with zeros.
                                          */
        void reinit (const std::vector<IndexSet> &parallel_partitioning,
		     const MPI_Comm              &communicator = MPI_COMM_WORLD,
		     const bool                   fast = false);

                                         /**
                                          * Change the dimension to that
                                          * of the vector <tt>V</tt>. The same
                                          * applies as for the other
                                          * reinit() function.
                                          *
                                          * The elements of <tt>V</tt> are not
                                          * copied, i.e.  this function is
                                          * the same as calling <tt>reinit
                                          * (V.size(), fast)</tt>.
                                          *
                                          * Note that you must call this
                                          * (or the other reinit()
                                          * functions) function, rather
                                          * than calling the reinit()
                                          * functions of an individual
                                          * block, to allow the block
                                          * vector to update its caches of
                                          * vector sizes. If you call
                                          * reinit() on one of the
                                          * blocks, then subsequent
                                          * actions on this object may
                                          * yield unpredictable results
                                          * since they may be routed to
                                          * the wrong block.
                                          */
        void reinit (const BlockVector &V,
		     const bool fast = false);

                                         /**
                                          * Change the number of blocks to
					  * <tt>num_blocks</tt>. The individual
					  * blocks will get initialized with
					  * zero size, so it is assumed that
					  * the user resizes the
					  * individual blocks by herself
					  * in an appropriate way, and
					  * calls <tt>collect_sizes</tt>
					  * afterwards.
                                          */
        void reinit (const unsigned int num_blocks);

				         /**
					  * This reinit function is meant to
					  * be used for parallel
					  * calculations where some
					  * non-local data has to be
					  * used. The typical situation
					  * where one needs this function is
					  * the call of the
					  * FEValues<dim>::get_function_values
					  * function (or of some
					  * derivatives) in parallel. Since
					  * it is usually faster to retrieve
					  * the data in advance, this
					  * function can be called before
					  * the assembly forks out to the
					  * different processors. What this
					  * function does is the following:
					  * It takes the information in the
					  * columns of the given matrix and
					  * looks which data couples between
					  * the different processors. That
					  * data is then queried from the
					  * input vector. Note that you
					  * should not write to the
					  * resulting vector any more, since
					  * the some data can be stored
					  * several times on different
					  * processors, leading to
					  * unpredictable results. In
					  * particular, such a vector cannot
					  * be used for matrix-vector
					  * products as for example done
					  * during the solution of linear
					  * systems.
					  */
	void import_nonlocal_data_for_fe (const TrilinosWrappers::BlockSparseMatrix &m,
					  const BlockVector                         &v);

					 /**
					  * Compress the underlying
					  * representation of the Trilinos
					  * object, i.e. flush the buffers
					  * of the vector object if it has
					  * any. This function is
					  * necessary after writing into a
					  * vector element-by-element and
					  * before anything else can be
					  * done on it.
					  *
					  * The (defaulted) argument can
					  * be used to specify the
					  * compress mode
					  * (<code>Add</code> or
					  * <code>Insert</code>) in case
					  * the vector has not been
					  * written to since the last
					  * time this function was
					  * called. The argument is
					  * ignored if the vector has
					  * been added or written to
					  * since the last time
					  * compress() was called.
					  *
					  * See @ref GlossCompress "Compressing distributed objects"
					  * for more information.
					  * more information.
					  */
	void compress (const Epetra_CombineMode last_action = Zero);

				         /**
					  * Returns the state of the
					  * vector, i.e., whether
					  * compress() needs to be
					  * called after an operation
					  * requiring data
					  * exchange. Does only return
					  * non-true values when used in
					  * <tt>debug</tt> mode, since
					  * it is quite expensive to
					  * keep track of all operations
					  * that lead to the need for
					  * compress().
					  */
	bool is_compressed () const;

                                         /**
                                          * Swap the contents of this
                                          * vector and the other vector
                                          * <tt>v</tt>. One could do this
                                          * operation with a temporary
                                          * variable and copying over the
                                          * data elements, but this
                                          * function is significantly more
                                          * efficient since it only swaps
                                          * the pointers to the data of
                                          * the two vectors and therefore
                                          * does not need to allocate
                                          * temporary storage and move
                                          * data around.
                                          *
                                          * Limitation: right now this
                                          * function only works if both
                                          * vectors have the same number
                                          * of blocks. If needed, the
                                          * numbers of blocks should be
                                          * exchanged, too.
                                          *
                                          * This function is analog to the
                                          * the swap() function of all C++
                                          * standard containers. Also,
                                          * there is a global function
                                          * swap(u,v) that simply calls
                                          * <tt>u.swap(v)</tt>, again in analogy
                                          * to standard functions.
                                          */
        void swap (BlockVector &v);

				     /**
				      * Print to a stream.
				      */
	void print (std::ostream       &out,
		    const unsigned int  precision = 3,
		    const bool          scientific = true,
		    const bool          across = true) const;

                                         /**
                                          * Exception
                                          */
        DeclException0 (ExcIteratorRangeDoesNotMatchVectorSize);

                                         /**
                                          * Exception
                                          */
        DeclException0 (ExcNonMatchingBlockVectors);
    };



/*----------------------- Inline functions ----------------------------------*/


    inline
    BlockVector::BlockVector ()
    {}



    inline
    BlockVector::BlockVector (const std::vector<Epetra_Map> &parallel_partitioning)
    {
      reinit (parallel_partitioning, false);
    }



    inline
    BlockVector::BlockVector (const std::vector<IndexSet> &parallel_partitioning,
			      const MPI_Comm              &communicator)
    {
      reinit (parallel_partitioning, communicator, false);
    }



    inline
    BlockVector::BlockVector (const unsigned int num_blocks)
    {
      reinit (num_blocks);
    }



    inline
    BlockVector::BlockVector (const BlockVector& v)
                    :
                    BlockVectorBase<Vector > ()
    {
      this->components.resize (v.n_blocks());
      this->block_indices = v.block_indices;

      for (unsigned int i=0; i<this->n_blocks(); ++i)
        this->components[i] = v.components[i];
    }



    inline
    bool
    BlockVector::is_compressed () const
    {
      bool compressed = true;
      for (unsigned int row=0; row<n_blocks(); ++row)
	if (block(row).is_compressed() == false)
	  {
	    compressed = false;
	    break;
	  }

      return compressed;
    }



    template <typename Number>
    BlockVector &
    BlockVector::operator = (const ::dealii::BlockVector<Number> &v)
    {
      if (n_blocks() != v.n_blocks())
	{
	  std::vector<unsigned int> block_sizes (v.n_blocks(), 0);
	  block_indices.reinit (block_sizes);
	  if (components.size() != n_blocks())
	    components.resize(n_blocks());
	}

      for (unsigned int i=0; i<this->n_blocks(); ++i)
	this->components[i] = v.block(i);

      collect_sizes();

      return *this;
    }


    inline
    void
    BlockVector::swap (BlockVector &v)
    {
      Assert (n_blocks() == v.n_blocks(),
	      ExcDimensionMismatch(n_blocks(),v.n_blocks()));

      for (unsigned int row=0; row<n_blocks(); ++row)
	block(row).swap (v.block(row));
    }



/**
 * Global function which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @relates TrilinosWrappers::MPI::BlockVector
 * @author Martin Kronbichler, Wolfgang Bangerth, 2008
 */
    inline
    void swap (BlockVector &u,
               BlockVector &v)
    {
      u.swap (v);
    }

  } /* end of namespace MPI */




/**
 * An implementation of block vectors based on the vector class
 * implemented in TrilinosWrappers. While the base class provides for
 * most of the interface, this class handles the actual allocation of
 * vectors and provides functions that are specific to the underlying
 * vector type.
 *
 * In contrast to the class MPI::BlockVector, this class is based on a
 * localized version of the vectors, which means that the whole vector
 * is stored on each processor. Note that matrix vector products with
 * this block vector class do only work in case the program is run on
 * only one processor, since the Trilinos matrices are inherently
 * parallel.
 *
 * @ingroup Vectors
 * @ingroup TrilinosWrappers
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Martin Kronbichler, 2008
 */
  class BlockVector : public BlockVectorBase<Vector>
  {
    public:
                                       /**
					* Typedef the base class for simpler
					* access to its own typedefs.
					*/
      typedef BlockVectorBase<Vector> BaseClass;

                                       /**
					* Typedef the type of the underlying
					* vector.
					*/
      typedef BaseClass::BlockType  BlockType;

                                       /**
					* Import the typedefs from the base
					* class.
					*/
      typedef BaseClass::value_type      value_type;
      typedef BaseClass::pointer         pointer;
      typedef BaseClass::const_pointer   const_pointer;
      typedef BaseClass::reference       reference;
      typedef BaseClass::const_reference const_reference;
      typedef BaseClass::size_type       size_type;
      typedef BaseClass::iterator        iterator;
      typedef BaseClass::const_iterator  const_iterator;

                                       /**
					* Default constructor. Generate an
					* empty vector without any blocks.
					*/
      BlockVector ();

                                       /**
					* Constructor. Generate a block
					* vector with as many blocks as
					* there are entries in Input_Maps.
					* For this non-distributed vector,
					* the %parallel partitioning is not
					* used, just the global size of the
					* partitioner.
					*/
      BlockVector (const std::vector<Epetra_Map> &partitioner);

                                       /**
					* Constructor. Generate a block
					* vector with as many blocks as
					* there are entries in Input_Maps.
					* For this non-distributed vector,
					* the %parallel partitioning is not
					* used, just the global size of the
					* partitioner.
					*/
      BlockVector (const std::vector<IndexSet> &partitioner,
		   const MPI_Comm              &communicator = MPI_COMM_WORLD);

                                       /**
					* Copy-Constructor. Set all the
					* properties of the non-%parallel
					* vector to those of the given
					* %parallel vector and import the
					* elements.
					*/
      BlockVector (const MPI::BlockVector &V);

                                       /**
					* Copy-Constructor. Set all the
					* properties of the vector to those
					* of the given input vector and copy
					* the elements.
					*/
      BlockVector (const BlockVector  &V);

                                       /**
					* Creates a block vector
					* consisting of
					* <tt>num_blocks</tt>
					* components, but there is no
					* content in the individual
					* components and the user has to
					* fill appropriate data using a
					* reinit of the blocks.
					*/
      BlockVector (const unsigned int num_blocks);

                                       /**
					* Constructor. Set the number of
					* blocks to <tt>n.size()</tt> and
					* initialize each block with
					* <tt>n[i]</tt> zero elements.
					*
					* References BlockVector.reinit().
					*/
      BlockVector (const std::vector<unsigned int> &N);

                                       /**
                                        * Constructor. Set the number of
                                        * blocks to
                                        * <tt>n.size()</tt>. Initialize the
                                        * vector with the elements
                                        * pointed to by the range of
                                        * iterators given as second and
                                        * third argument. Apart from the
                                        * first argument, this
                                        * constructor is in complete
                                        * analogy to the respective
                                        * constructor of the
                                        * <tt>std::vector</tt> class, but the
                                        * first argument is needed in
                                        * order to know how to subdivide
                                        * the block vector into
                                        * different blocks.
                                        */
      template <typename InputIterator>
      BlockVector (const std::vector<unsigned int> &n,
                   const InputIterator              first,
                   const InputIterator              end);

                                       /**
					* Destructor. Clears memory
					*/
      ~BlockVector ();

                                       /**
					* Copy operator: fill all
					* components of the vector that
					* are locally stored with the
					* given scalar value.
					*/
      BlockVector &
        operator = (const value_type s);

                                       /**
					* Copy operator for a
					* distributed Trilinos vector to
					* a localized one.
					*/
      BlockVector &
	operator = (const MPI::BlockVector &V);

                                       /**
					* Copy operator for arguments of
					* the same type.
					*/
      BlockVector &
	operator = (const BlockVector &V);

                                       /**
					* Another copy function. This
					* one takes a deal.II block
					* vector and copies it into a
					* TrilinosWrappers block
					* vector. Note that the number
					* of blocks has to be the same
					* in the vector as in the input
					* vector. Use the reinit()
					* command for resizing the
					* BlockVector or for changing
					* the internal structure of the
					* block components.
					*
					* Since Trilinos only works on
					* doubles, this function is
					* limited to accept only one
					* possible number type in the
					* deal.II vector.
					*/
      template <typename Number>
      BlockVector &
	operator = (const ::dealii::BlockVector<Number> &V);

                                         /**
                                          * Reinitialize the BlockVector to
                                          * contain as many blocks as there
                                          * are Epetra_Maps given in the
                                          * input argument, according to the
                                          * global size of the individual
                                          * components described in the
                                          * maps. Note that the resulting
                                          * vector will be stored completely
                                          * on each process. The Epetra_Map
                                          * is useful when data exchange
                                          * with a distributed vector based
                                          * on the same Epetra_map is
                                          * intended. In that case, the same
                                          * communicator is used for data
                                          * exchange.
                                          *
                                          * If <tt>fast==false</tt>, the vector
                                          * is filled with zeros.
                                          */
      void reinit (const std::vector<Epetra_Map> &partitioning,
		   const bool                     fast = false);

                                         /**
                                          * Reinitialize the BlockVector to
                                          * contain as many blocks as there
                                          * are index sets given in the
                                          * input argument, according to the
                                          * global size of the individual
                                          * components described in the
                                          * index set, and using a given MPI
                                          * communicator. The MPI
                                          * communicator is useful when data
                                          * exchange with a distributed
                                          * vector based on the same
                                          * initialization is intended. In
                                          * that case, the same communicator
                                          * is used for data exchange.
                                          *
                                          * If <tt>fast==false</tt>, the vector
                                          * is filled with zeros.
                                          */
      void reinit (const std::vector<IndexSet> &partitioning,
		   const MPI_Comm              &communicator = MPI_COMM_WORLD,
		   const bool                   fast = false);

                                         /**
                                          * Reinitialize the BlockVector to
                                          * contain as many blocks as there
                                          * are elements in the first
                                          * argument, and with the respective
                                          * sizes. Since no distribution map
                                          * is given, all vectors are local
                                          * vectors.
                                          *
                                          * If <tt>fast==false</tt>, the vector
                                          * is filled with zeros.
                                          */
      void reinit (const std::vector<unsigned int> &N,
		   const bool                       fast=false);

                                         /**
                                          * Reinit the function
                                          * according to a distributed
                                          * block vector. The elements
                                          * will be copied in this
                                          * process.
                                          */
      void reinit (const MPI::BlockVector &V);

                                         /**
                                          * Change the dimension to that
                                          * of the vector <tt>V</tt>. The same
                                          * applies as for the other
                                          * reinit() function.
                                          *
                                          * The elements of <tt>V</tt> are not
                                          * copied, i.e.  this function is
                                          * the same as calling <tt>reinit
                                          * (V.size(), fast)</tt>.
                                          *
                                          * Note that you must call this
                                          * (or the other reinit()
                                          * functions) function, rather
                                          * than calling the reinit()
                                          * functions of an individual
                                          * block, to allow the block
                                          * vector to update its caches of
                                          * vector sizes. If you call
                                          * reinit() on one of the
                                          * blocks, then subsequent
                                          * actions on this object may
                                          * yield unpredictable results
                                          * since they may be routed to
                                          * the wrong block.
                                          */
      void reinit (const BlockVector &V,
		   const bool fast = false);

                                         /**
                                          * Change the number of blocks to
					  * <tt>num_blocks</tt>. The individual
					  * blocks will get initialized with
					  * zero size, so it is assumed that
					  * the user resizes the
					  * individual blocks by herself
					  * in an appropriate way, and
					  * calls <tt>collect_sizes</tt>
					  * afterwards.
                                          */
      void reinit (const unsigned int num_blocks);

                                         /**
                                          * Swap the contents of this
                                          * vector and the other vector
                                          * <tt>v</tt>. One could do this
                                          * operation with a temporary
                                          * variable and copying over the
                                          * data elements, but this
                                          * function is significantly more
                                          * efficient since it only swaps
                                          * the pointers to the data of
                                          * the two vectors and therefore
                                          * does not need to allocate
                                          * temporary storage and move
                                          * data around.
                                          *
                                          * Limitation: right now this
                                          * function only works if both
                                          * vectors have the same number
                                          * of blocks. If needed, the
                                          * numbers of blocks should be
                                          * exchanged, too.
                                          *
                                          * This function is analog to the
                                          * the swap() function of all C++
                                          * standard containers. Also,
                                          * there is a global function
                                          * swap(u,v) that simply calls
                                          * <tt>u.swap(v)</tt>, again in analogy
                                          * to standard functions.
                                          */
      void swap (BlockVector &v);

				     /**
				      * Print to a stream.
				      */
      void print (std::ostream       &out,
		  const unsigned int  precision = 3,
		  const bool          scientific = true,
		  const bool          across = true) const;

                                         /**
                                          * Exception
                                          */
      DeclException0 (ExcIteratorRangeDoesNotMatchVectorSize);

                                         /**
                                          * Exception
                                          */
      DeclException0 (ExcNonMatchingBlockVectors);

                                       /**
					* Exception
					*/
      DeclException2 (ExcNonLocalizedMap,
      		      int, int,
      		      << "For the generation of a localized vector the map has "
		      << "to assign all elements to all vectors! "
		      << "local_size = global_size is a necessary condition, but"
		      << arg1 << " != " << arg2 << " was given!");

  };



/*----------------------- Inline functions ----------------------------------*/



  inline
  BlockVector::BlockVector ()
  {}



  inline
  BlockVector::BlockVector (const std::vector<Epetra_Map> &partitioning)
  {
    reinit (partitioning);
  }



  inline
  BlockVector::BlockVector (const std::vector<IndexSet> &partitioning,
			    const MPI_Comm              &communicator)
  {
    reinit (partitioning, communicator);
  }



  inline
  BlockVector::BlockVector (const std::vector<unsigned int> &N)
  {
    reinit (N);
  }



  template <typename InputIterator>
  BlockVector::BlockVector (const std::vector<unsigned int> &n,
                            const InputIterator              first,
                            const InputIterator              end)
  {
                                     // first set sizes of blocks, but
                                     // don't initialize them as we will
                                     // copy elements soon
    reinit (n, true);
    InputIterator start = first;
    for (unsigned int b=0; b<n.size(); ++b)
      {
        InputIterator end = start;
        std::advance (end, static_cast<signed int>(n[b]));

	for (unsigned int i=0; i<n[b]; ++i, ++start)
	  this->block(b)(i) = *start;
      }
    Assert (start == end, ExcIteratorRangeDoesNotMatchVectorSize());
  }



  inline
  BlockVector::BlockVector (const unsigned int num_blocks)
  {
    reinit (num_blocks);
  }



  inline
  BlockVector::~BlockVector()
  {}



  inline
  BlockVector::BlockVector (const MPI::BlockVector &v)
  {
    reinit (v);
  }



  inline
  BlockVector::BlockVector (const BlockVector& v)
                    :
                    BlockVectorBase<Vector > ()
  {
    this->components.resize (v.n_blocks());
    this->block_indices = v.block_indices;

    for (unsigned int i=0; i<this->n_blocks(); ++i)
      this->components[i] = v.components[i];
  }


  inline
  void
  BlockVector::swap (BlockVector &v)
  {
    Assert (n_blocks() == v.n_blocks(),
	    ExcDimensionMismatch(n_blocks(),v.n_blocks()));

    for (unsigned int row=0; row<n_blocks(); ++row)
      block(row).swap (v.block(row));
  }


  template <typename Number>
  BlockVector &
  BlockVector::operator = (const ::dealii::BlockVector<Number> &v)
  {
    if (n_blocks() != v.n_blocks())
      {
	std::vector<unsigned int> block_sizes (v.n_blocks(), 0);
	block_indices.reinit (block_sizes);
	if (components.size() != n_blocks())
	  components.resize(n_blocks());
      }

    for (unsigned int i=0; i<this->n_blocks(); ++i)
      this->components[i] = v.block(i);

    collect_sizes();

    return *this;
  }


/**
 * Global function which overloads the default implementation
 * of the C++ standard library which uses a temporary object. The
 * function simply exchanges the data of the two vectors.
 *
 * @relates TrilinosWrappers::BlockVector
 * @author Martin Kronbichler, 2008
 */
  inline
  void swap (BlockVector &u,
	     BlockVector &v)
  {
    u.swap (v);
  }

}

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif  // DEAL_II_USE_TRILINOS

#endif
