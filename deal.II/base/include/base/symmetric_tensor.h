//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__symmetric_tensor_h
#define __deal2__symmetric_tensor_h


#include <base/tensor.h>
#include <base/table_indices.h>


template <int rank, int dim> class SymmetricTensor;

template <int dim> SymmetricTensor<2,dim> unit_symmetric_tensor ();
template <int dim> SymmetricTensor<4,dim> deviator_tensor ();
template <int dim> SymmetricTensor<4,dim> identity_tensor ();
template <int dim> SymmetricTensor<4,dim> invert (const SymmetricTensor<4,dim> &);
template <int dim2> double trace (const SymmetricTensor<2,dim2> &);

template <int dim> SymmetricTensor<2,dim>
deviator (const SymmetricTensor<2,dim> &);


namespace internal
{
                                   /**
                                    * A namespace for classes that are
                                    * internal to how the SymmetricTensor
                                    * class works.
                                    */
  namespace SymmetricTensorAccessors
  {
				     /**
				      * Create a TableIndices<2>
				      * object where the first entries
				      * up to <tt>position-1</tt> are
				      * taken from previous_indices,
				      * and new_index is put at
				      * position
				      * <tt>position</tt>. The
				      * remaining indices remain in
				      * invalid state.
				      */
    inline
    TableIndices<2> merge (const TableIndices<2> &previous_indices,
			   const unsigned int     new_index,
			   const unsigned int     position)
    {
      Assert (position < 2, ExcIndexRange (position, 0, 2));
      
      if (position == 0)
	return TableIndices<2>(new_index);
      else
	return TableIndices<2>(previous_indices[0], new_index);
    }



				     /**
				      * Create a TableIndices<4>
				      * object where the first entries
				      * up to <tt>position-1</tt> are
				      * taken from previous_indices,
				      * and new_index is put at
				      * position
				      * <tt>position</tt>. The
				      * remaining indices remain in
				      * invalid state.
				      */
    inline
    TableIndices<4> merge (const TableIndices<4> &previous_indices,
			   const unsigned int     new_index,
			   const unsigned int     position)
    {
      Assert (position < 4, ExcIndexRange (position, 0, 4));

      switch (position)
	{
	  case 0:
		return TableIndices<4>(new_index);
	  case 1:
		return TableIndices<4>(previous_indices[0],
				       new_index);
	  case 2:
		return TableIndices<4>(previous_indices[0],
				       previous_indices[1],
				       new_index); 
	  case 3:
		return TableIndices<4>(previous_indices[0],
				       previous_indices[1],
				       previous_indices[2],
				       new_index);
	}
      Assert (false, ExcInternalError());
      return TableIndices<4>();
    }


				     /**
				      * Typedef template magic
				      * denoting the result of a
				      * double contraction between two
				      * tensors or ranks rank1 and
				      * rank2. In general, this is a
				      * tensor of rank
				      * <tt>rank1+rank2-4</tt>, but if
				      * this is zero it is a single
				      * scalar double. For this case,
				      * we have a specialization.
				      *
				      * @author Wolfgang Bangerth, 2005
				      */
    template <int rank1, int rank2, int dim>
    struct double_contraction_result 
    {
	typedef SymmetricTensor<rank1+rank2-4,dim> type;
    };
    

				     /**
				      * Typedef template magic
				      * denoting the result of a
				      * double contraction between two
				      * tensors or ranks rank1 and
				      * rank2. In general, this is a
				      * tensor of rank
				      * <tt>rank1+rank2-4</tt>, but if
				      * this is zero it is a single
				      * scalar double. For this case,
				      * we have a specialization.
				      *
				      * @author Wolfgang Bangerth, 2005
				      */
    template <int dim>
    struct double_contraction_result<2,2,dim> 
    {
	typedef double type;
    };
    
    
    
                                     /**
                                      * Declaration of typedefs for the type
                                      * of data structures which are used to
                                      * store symmetric tensors. For example,
                                      * for rank-2 symmetric tensors, we use a
                                      * flat vector to store all the
                                      * elements. On the other hand, symmetric
                                      * rank-4 tensors are mappings from
                                      * symmetric rank-2 tensors into
                                      * symmetric rank-2 tensors, so they can
                                      * be represented as matrices, etc.
                                      *
                                      * This information is probably of little
                                      * interest to all except the accessor
                                      * classes that need it. In particular,
                                      * you shouldn't make any assumptions
                                      * about the storage format in your
                                      * application programs.
                                      */
    template <int rank, int dim>
    struct StorageType;

                                     /**
                                      * Specialization of StorageType for
                                      * rank-2 tensors.
                                      */
    template <int dim>
    struct StorageType<2,dim> 
    {
                                         /**
                                          * Number of independent components of a
                                          * symmetric tensor of rank 2. We store
                                          * only the upper right half of it.
                                          */
        static const unsigned int
        n_independent_components = (dim*dim + dim)/2;

                                         /**
                                          * Declare the type in which we actually
                                          * store the data.
                                          */
        typedef Tensor<1,n_independent_components> base_tensor_type;
    };



                                     /**
                                      * Specialization of StorageType for
                                      * rank-4 tensors.
                                      */
    template <int dim>
    struct StorageType<4,dim> 
    {
                                         /**
                                          * Number of independent components
                                          * of a symmetric tensor of rank
                                          * 2. Since rank-4 tensors are
                                          * mappings between such objects, we
                                          * need this information.
                                          */
        static const unsigned int
        n_rank2_components = (dim*dim + dim)/2;

                                         /**
                                          * Declare the type in which we
                                          * actually store the data. Symmetric
                                          * rank-4 tensors are mappings
                                          * between symmetric rank-2 tensors,
                                          * so we can represent the data as a
                                          * matrix if we represent the rank-2
                                          * tensors as vectors.
                                          */
        typedef Tensor<2,n_rank2_components> base_tensor_type;
    };
    
    

                                     /**
                                      * Switch type to select a tensor of
                                      * rank 2 and dimension <tt>dim</tt>,
                                      * switching on whether the tensor
                                      * should be constant or not.
                                      */
    template <int rank, int dim, bool constness>
    struct AccessorTypes;

                                     /**
                                      * Switch type to select a tensor of
                                      * rank 2 and dimension <tt>dim</tt>,
                                      * switching on whether the tensor
                                      * should be constant or not.
                                      *
                                      * Specialization for constant tensors.
                                      */
    template <int rank, int dim>
    struct AccessorTypes<rank, dim,true>
    {
        typedef const SymmetricTensor<rank,dim> tensor_type;

        typedef double reference;
    };

                                     /**
                                      * Switch type to select a tensor of
                                      * rank 2 and dimension <tt>dim</tt>,
                                      * switching on whether the tensor
                                      * should be constant or not.
                                      *
                                      * Specialization for non-constant
                                      * tensors.
                                      */
    template <int rank, int dim>
    struct AccessorTypes<rank,dim,false>
    {
        typedef SymmetricTensor<rank,dim> tensor_type;

        typedef double &reference;
    };


/**
 * @internal
 *
 * Class that acts as accessor to elements of type
 * SymmetricTensor. The template parameter <tt>C</tt> may be either
 * true or false, and indicates whether the objects worked on are
 * constant or not (i.e. write access is only allowed if the value is
 * false).
 *
 * Since with <tt>N</tt> indices, the effect of applying
 * <tt>operator[]</tt> is getting access to something we <tt>N-1</tt>
 * indices, we have to implement these accessor classes recursively,
 * with stopping when we have only one index left. For the latter
 * case, a specialization of this class is declared below, where
 * calling <tt>operator[]</tt> gives you access to the objects
 * actually stored by the tensor; the tensor class also makes sure
 * that only those elements are actually accessed which we actually
 * store, i.e. it reorders indices if necessary. The template
 * parameter <tt>P</tt> indicates how many remaining indices there
 * are. For a rank-2 tensor, <tt>P</tt> may be two, and when using
 * <tt>operator[]</tt>, an object with <tt>P=1</tt> emerges.
 *
 * As stated for the entire namespace, you will not usually have to do
 * with these classes directly, and should not try to use their
 * interface directly as it may change without notice. In fact, since
 * the constructors are made private, you will not even be able to
 * generate objects of this class, as they are only thought as
 * temporaries for access to elements of the table class, not for
 * passing them around as arguments of functions, etc.
 *
 * This class is an adaptation of a similar class used for the Table
 * class.
 *
 * @author Wolfgang Bangerth, 2002, 2005
 */
    template <int rank, int dim, bool constness, int P>
    class Accessor
    {
      public:
                                         /**
                                          * Import two typedefs from the
                                          * switch class above.
                                          */
        typedef typename AccessorTypes<rank,dim,constness>::reference reference;
        typedef typename AccessorTypes<rank,dim,constness>::tensor_type tensor_type;

      private:
                                         /**
                                          * Constructor. Take a
                                          * reference to the tensor
                                          * object which we will
                                          * access.
                                          *
					  * The second argument
					  * denotes the values of
					  * previous indices into the
					  * tensor. For example, for a
					  * rank-4 tensor, if P=2,
					  * then we will already have
					  * had two successive element
					  * selections (e.g. through
					  * <tt>tensor[1][2]</tt>),
					  * and the two index values
					  * have to be stored
					  * somewhere. This class
					  * therefore only makes use
					  * of the first rank-P
					  * elements of this array,
					  * but passes it on to the
					  * next level with P-1 which
					  * fills the next entry, and
					  * so on.
					  *
                                          * The constructor is made
                                          * private in order to prevent
                                          * you having such objects
                                          * around. The only way to
                                          * create such objects is via
                                          * the <tt>Table</tt> class, which
                                          * only generates them as
                                          * temporary objects. This
                                          * guarantees that the accessor
                                          * objects go out of scope
                                          * earlier than the mother
                                          * object, avoid problems with
                                          * data consistency.
                                          */
        Accessor (tensor_type              &tensor,
                  const TableIndices<rank> &previous_indices);

                                         /**
                                          * Default constructor. Not
                                          * needed, and invisible, so
                                          * private.
                                          */
        Accessor ();
        
                                         /**
                                          * Copy constructor. Not
                                          * needed, and invisible, so
                                          * private.
                                          */
        Accessor (const Accessor &a);

      public:
      
                                         /**
                                          * Index operator.
                                          */
        Accessor<rank,dim,constness,P-1> operator [] (const unsigned int i);

      private:
                                         /**
                                          * Store the data given to the
                                          * constructor.
                                          */
        tensor_type             &tensor;
        const TableIndices<rank> previous_indices;

                                         // declare some other classes
                                         // as friends. make sure to
                                         // work around bugs in some
                                         // compilers
#ifndef DEAL_II_NAMESP_TEMPL_FRIEND_BUG
        template <int,int> friend class SymmetricTensor;
        template <int,int,bool,int>
        friend class Accessor;
#  ifndef DEAL_II_TEMPL_SPEC_FRIEND_BUG
        friend class SymmetricTensor<rank,dim>;
        friend class Accessor<rank,dim,constness,P+1>;
#  endif
#else
        friend class SymmetricTensor<rank,dim>;
        friend class Accessor<rank,dim,constness,P+1>;
#endif
    };


  
/**
 * @internal
 * Accessor class for SymmetricTensor. This is the specialization for the last
 * index, which actually allows access to the elements of the table,
 * rather than recursively returning access objects for further
 * subsets. The same holds for this specialization as for the general
 * template; see there for more information.
 *
 * @author Wolfgang Bangerth, 2002, 2005
 */
    template <int rank, int dim, bool constness>
    class Accessor<rank,dim,constness,1>
    {
      public:
                                         /**
                                          * Import two typedefs from the
                                          * switch class above.
                                          */
        typedef typename AccessorTypes<rank,dim,constness>::reference reference;
        typedef typename AccessorTypes<rank,dim,constness>::tensor_type tensor_type;

      private:
                                         /**
                                          * Constructor. Take a
                                          * reference to the tensor
                                          * object which we will
                                          * access.
                                          *
					  * The second argument
					  * denotes the values of
					  * previous indices into the
					  * tensor. For example, for a
					  * rank-4 tensor, if P=2,
					  * then we will already have
					  * had two successive element
					  * selections (e.g. through
					  * <tt>tensor[1][2]</tt>),
					  * and the two index values
					  * have to be stored
					  * somewhere. This class
					  * therefore only makes use
					  * of the first rank-P
					  * elements of this array,
					  * but passes it on to the
					  * next level with P-1 which
					  * fills the next entry, and
					  * so on.
					  *
					  * For this particular
					  * specialization, i.e. for
					  * P==1, all but the last
					  * index are already filled.
					  *
                                          * The constructor is made
                                          * private in order to prevent
                                          * you having such objects
                                          * around. The only way to
                                          * create such objects is via
                                          * the <tt>Table</tt> class, which
                                          * only generates them as
                                          * temporary objects. This
                                          * guarantees that the accessor
                                          * objects go out of scope
                                          * earlier than the mother
                                          * object, avoid problems with
                                          * data consistency.
                                          */
        Accessor (tensor_type              &tensor,
                  const TableIndices<rank> &previous_indices);      

                                         /**
                                          * Default constructor. Not
                                          * needed, and invisible, so
                                          * private.
                                          */
        Accessor ();

                                         /**
                                          * Copy constructor. Not
                                          * needed, and invisible, so
                                          * private.
                                          */
        Accessor (const Accessor &a);

      public:
      
                                         /**
                                          * Index operator.
                                          */
        reference operator [] (const unsigned int);
      
      private:
                                         /**
                                          * Store the data given to the
                                          * constructor.
                                          */
        tensor_type             &tensor;
        const TableIndices<rank> previous_indices;

                                         // declare some other classes
                                         // as friends. make sure to
                                         // work around bugs in some
                                         // compilers
#ifndef DEAL_II_NAMESP_TEMPL_FRIEND_BUG
        template <int,int> friend class SymmetricTensor;
        template <int,int,bool,int>
        friend class Accessor;
#  ifndef DEAL_II_TEMPL_SPEC_FRIEND_BUG
        friend class SymmetricTensor<rank,dim>;
        friend class Accessor<rank,dim,constness,2>;
#  endif
#else
        friend class SymmetricTensor<rank,dim>;
        friend class Accessor<rank,dim,constness,2>;
#endif
    };
  }
}

          


/**
 * Provide a class that stores symmetric tensors of rank 2,4,... efficiently,
 * i.e. only store those off-diagonal elements of the full tensor that are not
 * redundant. For example, for symmetric 2x2 tensors, this would be the
 * elements 11, 22, and 12, while the element 21 is equal to the 12 element.
 *
 * Using this class for symmetric tensors of rank 2 has advantages over
 * matrices in many cases since the dimension is known to the compiler as well
 * as the location of the data. It is therefore possible to produce far more
 * efficient code than for matrices with runtime-dependent dimension. It is
 * also more efficient than using the more general <tt>Tensor</tt> class,
 * since less elements are stored, and the class automatically makes sure that
 * the tensor represents a symmetric object.
 *
 * For tensors of higher rank, the savings in storage are even higher. For
 * example for the 3x3x3x3 tensors of rank 4, only 36 instead of the full 81
 * entries have to be stored.
 *
 * While the definition of a symmetric rank-2 tensor is obvious,
 * tensors of rank 4 are considered symmetric if they are operators
 * mapping symmetric rank-2 tensors onto symmetric rank-2
 * tensors. This entails certain symmetry properties on the elements
 * in their 4-dimensional index space, in particular that
 * <tt>C<sub>ijkl</sub>=C<sub>jikl</sub>=C<sub>ijlk</sub></tt>. However,
 * it does not imply the relation
 * <tt>C<sub>ijkl</sub>=C<sub>klij</sub></tt>. Consequently, symmetric
 * tensors of rank 4 as understood here are only tensors that map
 * symmetric tensors onto symmetric tensors, but they do not
 * necessarily induce a symmetric scalar product <tt>a:C:b=b:C:a</tt>
 * or even a positive (semi-)definite form <tt>a:C:a</tt>, where
 * <tt>a,b</tt> are symmetric rank-2 tensors and the colon indicates
 * the common double-index contraction that acts as a product for
 * symmetric tensors.
 *
 * Symmetric tensors are most often used in structural and fluid mechanics,
 * where strains and stresses are usually symmetric tensors, and the
 * stress-strain relationship is given by a symmetric rank-4 tensor.
 *
 * Note that symmetric tensors only exist with even numbers of indices. In
 * other words, the only objects that you can use are
 * <tt>SymmetricTensor<2,dim></tt>, <tt>SymmetricTensor<4,dim></tt>, etc, but
 * <tt>SymmetricTensor<1,dim></tt> and <tt>SymmetricTensor<3,dim></tt> do not
 * exist and their use will most likely lead to compiler errors.
 *
 *
 * <h3>Accessing elements</h3>
 *
 * The elements of a tensor <tt>t</tt> can be accessed using the
 * bracket operator, i.e. for a tensor of rank 4,
 * <tt>t[0][1][0][1]</tt> accesses the element
 * <tt>t<sub>0101</sub></tt>. This access can be used for both reading
 * and writing (if the tensor is non-constant at least). You may also
 * perform other operations on it, although that may lead to confusing
 * situations because several elements of the tensor are stored at the
 * same location. For example, for a rank-2 tensor that is assumed to
 * be zero at the beginning, writing <tt>t[0][1]+=1; t[1][0]+=1;</tt>
 * will lead to the same element being increased by one
 * <em>twice</em>, because even though the accesses use different
 * indices, the elements that are accessed are symmetric and therefore
 * stored at the same location. It may therefore be useful in
 * application programs to restrict operations on individual elements
 * to simple reads or writes.
 *
 * @author Wolfgang Bangerth, 2005
 */
template <int rank, int dim>
class SymmetricTensor
{
  public:
				     /**
				      * Provide a way to get the
				      * dimension of an object without
				      * explicit knowledge of it's
				      * data type. Implementation is
				      * this way instead of providing
				      * a function <tt>dimension()</tt>
				      * because now it is possible to
				      * get the dimension at compile
				      * time without the expansion and
				      * preevaluation of an inlined
				      * function; the compiler may
				      * therefore produce more
				      * efficient code and you may use
				      * this value to declare other
				      * data types.
				      */
    static const unsigned int dimension = dim;
    
                                     /**
                                      * Default constructor. Creates a zero
                                      * tensor.
                                      */
    SymmetricTensor ();

                                     /**
                                      * Constructor. Generate a symmetric
                                      * tensor from a general one. Assumes
                                      * that @p t is already symmetric, and in
                                      * debug mode this is in fact
                                      * checked. Note that no provision is
                                      * made to assure that the tensor is
                                      * symmetric only up to round-off error:
                                      * if the incoming tensor is not exactly
                                      * symmetric, then an exception is
                                      * thrown. If you know that incoming
                                      * tensor is symmetric only up to
                                      * round-off, then you may want to call
                                      * the <tt>symmetrize</tt> function
                                      * first. If you aren't sure, it is good
                                      * practice to check before calling
                                      * <tt>symmetrize</tt>.
                                      */
    SymmetricTensor (const Tensor<2,dim> &t);

				     /**
				      *  Assignment operator.
				      */
    SymmetricTensor & operator = (const SymmetricTensor &);

                                     /**
                                      * Convert the present symmetric tensor
                                      * into a full tensor with the same
                                      * elements, but using the different
                                      * storage scheme of full tensors.
                                      */
    operator Tensor<rank,dim> () const;

				     /**
				      *  Test for equality of two tensors.
				      */
    bool operator == (const SymmetricTensor &) const;

    				     /**
				      *  Test for inequality of two tensors.
				      */
    bool operator != (const SymmetricTensor &) const;

				     /**
				      *  Add another tensor.
				      */
    SymmetricTensor & operator += (const SymmetricTensor &);
    
				     /**
				      *  Subtract another tensor.
				      */
    SymmetricTensor & operator -= (const SymmetricTensor &);

				     /**
				      *  Scale the tensor by <tt>factor</tt>,
				      *  i.e. multiply all components by
				      *  <tt>factor</tt>.
				      */
    SymmetricTensor & operator *= (const double factor);

				     /**
				      *  Scale the vector by
				      *  <tt>1/factor</tt>.
				      */
    SymmetricTensor & operator /= (const double factor);

				     /**
				      *  Add two tensors. If possible, you
				      *  should use <tt>operator +=</tt>
				      *  instead since this does not need the
				      *  creation of a temporary.
				      */
    SymmetricTensor   operator + (const SymmetricTensor &s) const;

				     /**
				      *  Subtract two tensors. If possible,
				      *  you should use <tt>operator -=</tt>
				      *  instead since this does not need the
				      *  creation of a temporary.
				      */
    SymmetricTensor   operator - (const SymmetricTensor &s) const;

				     /**
				      * Unary minus operator. Negate all
				      * entries of a tensor.
				      */
    SymmetricTensor   operator - () const;

                                     /**
				      * Product between the present
				      * symmetric tensor and a tensor
				      * of rank 2. For example, if the
				      * present object is also a
				      * rank-2 tensor, then this is
				      * the scalar-product double
				      * contraction
				      * <tt>a<sub>ij</sub>b<sub>ij</sub></tt>
				      * over all indices
				      * <tt>i,j</tt>. In this case,
				      * the return value evaluates to
				      * a single scalar. While it is
				      * possible to define other
				      * scalar product (and associated
				      * induced norms), this one seems
				      * to be the most appropriate
				      * one.
				      *
				      * If the present object is a
				      * rank-4 tensor, the the result
				      * is a rank-2 tensor, the
				      * operation contracts over the
				      * last two indices of the
				      * present object and the indices
				      * of the argument, and the
				      * result is a tensor of rank 2.
				      *
				      * Note that the multiplication
				      * operator for symmetrict
				      * tensors is defined to be a
				      * double contraction over two
				      * indices, while it is defined
				      * as a single contraction over
				      * only one index for regular
				      * <tt>Tensor</tt> objects. For
				      * symmetric tensors it therefore
				      * acts in a way that is commonly
				      * denoted by a "colon
				      * multiplication" in the
				      * mathematica literature.
				      *
				      * There are global functions
				      * <tt>double_contract</tt> that
				      * do the same work as this
				      * operator, but rather than
				      * returning the result as a
				      * return value, they write it
				      * into the first argument to the
				      * function.
                                      */
    typename internal::SymmetricTensorAccessors::double_contraction_result<rank,2,dim>::type
    operator * (const SymmetricTensor<2,dim> &s) const;

				     /**
				      * Contraction over two indices
				      * of the present object with the
				      * rank-4 symmetric tensor given
				      * as argument.
				      */
    typename internal::SymmetricTensorAccessors::double_contraction_result<rank,4,dim>::type
    operator * (const SymmetricTensor<4,dim> &s) const;
    
                                     /**
                                      * Return a read-write reference
                                      * to the indicated element.
                                      */
    double & operator() (const TableIndices<rank> &indices);
  
                                     /**
                                      * Return the value of the
                                      * indicated element as a
                                      * read-only reference.
                                      *
                                      * We return the requested value
                                      * as a constant reference rather
                                      * than by value since this
                                      * object may hold data types
                                      * that may be large, and we
                                      * don't know here whether
                                      * copying is expensive or not.
                                      */
    double operator() (const TableIndices<rank> &indices) const;

                                     /**
                                      * Access the elements of a row of this
                                      * symmetric tensor. This function is
                                      * called for constant tensors.
                                      */
    internal::SymmetricTensorAccessors::Accessor<rank,dim,true,rank-1>
    operator [] (const unsigned int row) const;

                                     /**
                                      * Access the elements of a row of this
                                      * symmetric tensor. This function is
                                      * called for non-constant tensors.
                                      */
    internal::SymmetricTensorAccessors::Accessor<rank,dim,false,rank-1>
    operator [] (const unsigned int row);
    
                                     /**
                                      * Return the Frobenius-norm of a tensor,
                                      * i.e. the square root of the sum of
                                      * squares of all entries. This norm is
                                      * induced by the scalar product defined
                                      * above for two symmetric tensors. Note
                                      * that it includes <i>all</i> entries of
                                      * the tensor, counting symmetry, not
                                      * only the unique ones (for example, for
                                      * rank-2 tensors, this norm includes
                                      * adding up the squares of upper right
                                      * as well as lower left entries, not
                                      * just one of them, although they are
                                      * equal for symmetric tensors).
                                      */
    double norm () const;
    
    				     /**
				      * Reset all values to zero.
				      *
				      * Note that this is partly inconsistent
				      * with the semantics of the @p clear()
				      * member functions of the STL and of
				      * several other classes within deal.II
				      * which not only reset the values of
				      * stored elements to zero, but release
				      * all memory and return the object into
				      * a virginial state. However, since the
				      * size of objects of the present type is
				      * determined by its template parameters,
				      * resizing is not an option, and indeed
				      * the state where all elements have a
				      * zero value is the state right after
				      * construction of such an object.
				      */
    void clear ();

				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object.
				      */
    static unsigned int memory_consumption ();

    
  private:
				     /**
				      * A structure that describes
				      * properties of the base tensor.
				      */
    typedef 
    internal::SymmetricTensorAccessors::StorageType<rank,dim>
    base_tensor_descriptor;
    
                                     /**
                                      * Data storage type for a
                                      * symmetric tensor.
                                      */
    typedef typename base_tensor_descriptor::base_tensor_type base_tensor_type;
    
				     /**
				      * The place where we store the
				      * data of the tensor.
				      */
    base_tensor_type data;

				     /**
				      * Make all other symmetric tensors friends.
				      */
    template <int, int> friend class SymmetricTensor;

				     /**
				      * Make a few more functions friends.
				      */
    template <int dim2>
    friend double trace (const SymmetricTensor<2,dim2> &d);

    friend double determinant (const SymmetricTensor<2,1> &t);
    
    friend double determinant (const SymmetricTensor<2,2> &t);

    friend double determinant (const SymmetricTensor<2,3> &t);
    
    template <int dim2>
    friend SymmetricTensor<2,dim2>
    deviator (const SymmetricTensor<2,dim2> &t);
    
    template <int dim2>
    friend SymmetricTensor<2,dim2> unit_symmetric_tensor ();

    template <int dim2>
    friend SymmetricTensor<4,dim2> deviator_tensor ();

    template <int dim2>
    friend SymmetricTensor<4,dim2> identity_tensor ();

    template <int dim2>
    friend SymmetricTensor<4,dim2> invert (const SymmetricTensor<4,dim2> &);
};



// ------------------------- inline functions ------------------------

///@if NoDoc

namespace internal
{
  namespace SymmetricTensorAccessors
  {
    template <int rank, int dim, bool constness, int P>
    Accessor<rank,dim,constness,P>::
    Accessor (tensor_type              &tensor,
	      const TableIndices<rank> &previous_indices)
		    :
		    tensor (tensor),
		    previous_indices (previous_indices)
    {}



    template <int rank, int dim, bool constness, int P>
    Accessor<rank,dim,constness,P-1>
    Accessor<rank,dim,constness,P>::operator[] (const unsigned int i)
    {
      return Accessor<rank,dim,constness,P-1> (tensor,
					       merge (previous_indices, i, rank-P));
    }



    template <int rank, int dim, bool constness>
    Accessor<rank,dim,constness,1>::
    Accessor (tensor_type              &tensor,
	      const TableIndices<rank> &previous_indices)
		    :
		    tensor (tensor),
		    previous_indices (previous_indices)
    {}



    template <int rank, int dim, bool constness>
    typename Accessor<rank,dim,constness,1>::reference
    Accessor<rank,dim,constness,1>::operator[] (const unsigned int i)
    {
      return tensor(merge (previous_indices, i, rank-1));
    }
    
    
  }
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim>::SymmetricTensor ()
{}



template <>
inline
SymmetricTensor<2,2>::SymmetricTensor (const Tensor<2,2> &t)
{
  Assert (t[0][1] == t[1][0], ExcInternalError());

  data[0] = t[0][0];
  data[1] = t[1][1];
  data[2] = t[0][1];
}



template <>
inline
SymmetricTensor<2,3>::SymmetricTensor (const Tensor<2,3> &t)
{
  Assert (t[0][1] == t[1][0], ExcInternalError());
  Assert (t[0][2] == t[2][0], ExcInternalError());
  Assert (t[1][2] == t[2][1], ExcInternalError());
  
  data[0] = t[0][0];
  data[1] = t[1][1];
  data[2] = t[2][2];
  data[3] = t[0][1];
  data[4] = t[0][2];
  data[5] = t[1][2];
}


template <int rank, int dim>
inline
SymmetricTensor<rank,dim> &
SymmetricTensor<rank,dim>::operator = (const SymmetricTensor<rank,dim> &t)
{
  data = t.data;
  return *this;
}



template <>
inline
SymmetricTensor<2,1>::
operator Tensor<2,1> () const
{
  const double t[1][1] = {{data[0]}};
  return Tensor<2,1>(t);
}



template <>
inline
SymmetricTensor<2,2>::
operator Tensor<2,2> () const
{
  const double t[2][2] = {{data[0], data[2]},
                          {data[2], data[1]}};
  return Tensor<2,2>(t);
}



template <>
inline
SymmetricTensor<2,3>::
operator Tensor<2,3> () const
{
  const double t[3][3] = {{data[0], data[3], data[4]},
                          {data[3], data[1], data[5]},
                          {data[4], data[5], data[2]}};
  return Tensor<2,3>(t);
}



template <int rank, int dim>
inline
bool
SymmetricTensor<rank,dim>::operator == (const SymmetricTensor<rank,dim> &t) const
{
  return data == t.data;
}



template <int rank, int dim>
inline
bool
SymmetricTensor<rank,dim>::operator != (const SymmetricTensor<rank,dim> &t) const
{
  return data != t.data;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim> &
SymmetricTensor<rank,dim>::operator += (const SymmetricTensor<rank,dim> &t)
{
  data += t.data;
  return *this;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim> &
SymmetricTensor<rank,dim>::operator -= (const SymmetricTensor<rank,dim> &t)
{
  data -= t.data;
  return *this;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim> &
SymmetricTensor<rank,dim>::operator *= (const double d)
{
  data *= d;
  return *this;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim> &
SymmetricTensor<rank,dim>::operator /= (const double d)
{
  data /= d;
  return *this;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
SymmetricTensor<rank,dim>::operator + (const SymmetricTensor &t) const
{
  SymmetricTensor tmp = *this;
  tmp.data += t.data;
  return tmp;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
SymmetricTensor<rank,dim>::operator - (const SymmetricTensor &t) const
{
  SymmetricTensor tmp = *this;
  tmp.data -= t.data;
  return tmp;
}



template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
SymmetricTensor<rank,dim>::operator - () const
{
  SymmetricTensor tmp = *this;
  tmp.data = -tmp.data;
  return tmp;
}



template <int rank, int dim>
inline
void
SymmetricTensor<rank,dim>::clear ()
{
  data.clear ();
}



template <int rank, int dim>
inline
unsigned int
SymmetricTensor<rank,dim>::memory_consumption ()
{
  return
    internal::SymmetricTensorAccessors::StorageType<rank,dim>::memory_consumption ();
}



template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<2,2,1>::type
SymmetricTensor<2,1>::operator * (const SymmetricTensor<2,1> &s) const
{
  return data[0] * s.data[0];
}



template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<2,2,2>::type
SymmetricTensor<2,2>::operator * (const SymmetricTensor<2,2> &s) const
{
  return (data[0] * s.data[0] +
          data[1] * s.data[1] +
          2*data[2] * s.data[2]);
}



template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<2,2,3>::type
SymmetricTensor<2,3>::operator * (const SymmetricTensor<2,3> &s) const
{
  return (data[0] * s.data[0] +
          data[1] * s.data[1] +
          data[2] * s.data[2] +
          2*data[3] * s.data[3] +
          2*data[4] * s.data[4] +
          2*data[5] * s.data[5]);
}



template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<4,2,1>::type
SymmetricTensor<4,1>::
operator * (const SymmetricTensor<2,1> &s) const
{
  const unsigned int dim = 1;
  SymmetricTensor<2,dim> tmp;
  tmp.data[0] = data[0][0] * s.data[0];
  return tmp;
}



template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<2,4,1>::type
SymmetricTensor<2,1>::
operator * (const SymmetricTensor<4,1> &s) const
{
  const unsigned int dim = 1;
  SymmetricTensor<2,dim> tmp;
  tmp[0][0] = data[0] * s.data[0][0];
  return tmp;
}



template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<4,2,2>::type
SymmetricTensor<4,2>::
operator * (const SymmetricTensor<2,2> &s) const
{
  const unsigned int dim = 2;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<3; ++i)
    tmp.data[i] = data[i][0] * s.data[0] +
		  data[i][1] * s.data[1] +
		  2 * data[i][2] * s.data[2];

  return tmp;
}



template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<2,4,2>::type
SymmetricTensor<2,2>::
operator * (const SymmetricTensor<4,2> &s) const
{
  const unsigned int dim = 2;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<3; ++i)
    tmp.data[i] = data[0] * s.data[0][i] +
		  data[1] * s.data[1][i] +
		  2 * data[2] * s.data[2][i];

  return tmp;
}



template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<4,2,3>::type
SymmetricTensor<4,3>::
operator * (const SymmetricTensor<2,3> &s) const
{
  const unsigned int dim = 3;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<6; ++i)
    tmp.data[i] = data[i][0] * s.data[0] +
		  data[i][1] * s.data[1] +
		  data[i][2] * s.data[2] +
		  2 * data[i][3] * s.data[3] +
		  2 * data[i][4] * s.data[4] +
		  2 * data[i][5] * s.data[5];

  return tmp;
}



template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<2,4,3>::type
SymmetricTensor<2,3>::
operator * (const SymmetricTensor<4,3> &s) const
{
  const unsigned int dim = 3;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<6; ++i)
    tmp.data[i] = data[0] * s.data[0][i] +
		  data[1] * s.data[1][i] +
		  data[2] * s.data[2][i] +
		  2 * data[3] * s.data[3][i] +
		  2 * data[4] * s.data[4][i] +
		  2 * data[5] * s.data[5][i];

  return tmp;
}




template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<4,4,1>::type
SymmetricTensor<4,1>::
operator * (const SymmetricTensor<4,1> &s) const
{
  const unsigned int dim = 1;
  SymmetricTensor<4,dim> tmp;
  tmp.data[0][0] = data[0][0] * s.data[0][0];
  return tmp;
}



template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<4,4,2>::type
SymmetricTensor<4,2>::
operator * (const SymmetricTensor<4,2> &s) const
{
  const unsigned int dim = 2;
  SymmetricTensor<4,dim> tmp;
  for (unsigned int i=0; i<base_tensor_descriptor::n_rank2_components; ++i)
    for (unsigned int j=0; j<base_tensor_descriptor::n_rank2_components; ++j)
      tmp.data[i][j] = data[i][0] * s.data[0][j] +
		       data[i][1] * s.data[1][j] +
		       2*data[i][2] * s.data[2][j];
  return tmp;
}



template <>
inline
internal::SymmetricTensorAccessors::double_contraction_result<4,4,3>::type
SymmetricTensor<4,3>::
operator * (const SymmetricTensor<4,3> &s) const
{
  const unsigned int dim = 3;
  SymmetricTensor<4,dim> tmp;
  for (unsigned int i=0; i<base_tensor_descriptor::n_rank2_components; ++i)
    for (unsigned int j=0; j<base_tensor_descriptor::n_rank2_components; ++j)
      tmp.data[i][j] = data[i][0] * s.data[0][j] +
		       data[i][1] * s.data[1][j] +
		       data[i][2] * s.data[2][j] +
		       2*data[i][3] * s.data[3][j] +
  		       2*data[i][4] * s.data[4][j] +
		       2*data[i][5] * s.data[5][j];
  return tmp;
}



template <>
inline
double &
SymmetricTensor<2,1>::operator () (const TableIndices<2> &indices)
{
  const unsigned int rank = 2;
  const unsigned int dim  = 1;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

  return data[0];
}



template <>
inline
double
SymmetricTensor<2,1>::operator () (const TableIndices<2> &indices) const
{
  const unsigned int rank = 2;
  const unsigned int dim  = 1;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

  return data[0];
}



template <>
inline
double &
SymmetricTensor<2,2>::operator () (const TableIndices<2> &indices)
{
  const unsigned int rank = 2;
  const unsigned int dim  = 2;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

                                   // first treat the main diagonal
                                   // elements, which are stored
                                   // consecutively at the beginning
  if (indices[0] == indices[1])
    return data[indices[0]];

                                   // the rest is messier and requires a few
                                   // switches. at least for the 2x2 case it
                                   // is reasonably simple
  Assert (((indices[0]==1) && (indices[1]==0)) ||
          ((indices[0]==0) && (indices[1]==1)),
          ExcInternalError());  
  return data[2];
}



template <>
inline
double
SymmetricTensor<2,2>::operator () (const TableIndices<2> &indices) const
{
  const unsigned int rank = 2;
  const unsigned int dim  = 2;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

                                   // first treat the main diagonal
                                   // elements, which are stored
                                   // consecutively at the beginning
  if (indices[0] == indices[1])
    return data[indices[0]];

                                   // the rest is messier and requires a few
                                   // switches. at least for the 2x2 case it
                                   // is reasonably simple
  Assert (((indices[0]==1) && (indices[1]==0)) ||
          ((indices[0]==0) && (indices[1]==1)),
          ExcInternalError());  
  return data[2];
}



template <>
inline
double &
SymmetricTensor<2,3>::operator () (const TableIndices<2> &indices)
{
  const unsigned int rank = 2;
  const unsigned int dim  = 3;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

                                   // first treat the main diagonal
                                   // elements, which are stored
                                   // consecutively at the beginning
  if (indices[0] == indices[1])
    return data[indices[0]];

                                   // the rest is messier and requires a few
                                   // switches, but simpler if we just sort
                                   // our indices
  TableIndices<2> sorted_indices (indices);
  sorted_indices.sort ();
  
  if ((sorted_indices[0]==0) && (sorted_indices[1]==1))
    return data[3];
  else if ((sorted_indices[0]==0) && (sorted_indices[1]==2))
    return data[4];
  else if ((sorted_indices[0]==1) && (sorted_indices[1]==2))
    return data[5];
  else
    Assert (false, ExcInternalError());

  static double dummy_but_referenceable = 0;
  return dummy_but_referenceable;
}



template <>
inline
double
SymmetricTensor<2,3>::operator () (const TableIndices<2> &indices) const
{
  const unsigned int rank = 2;
  const unsigned int dim  = 3;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

                                   // first treat the main diagonal
                                   // elements, which are stored
                                   // consecutively at the beginning
  if (indices[0] == indices[1])
    return data[indices[0]];

                                   // the rest is messier and requires a few
                                   // switches, but simpler if we just sort
                                   // our indices
  TableIndices<2> sorted_indices (indices);
  sorted_indices.sort ();
  
  if ((sorted_indices[0]==0) && (sorted_indices[1]==1))
    return data[3];
  else if ((sorted_indices[0]==0) && (sorted_indices[1]==2))
    return data[4];
  else if ((sorted_indices[0]==1) && (sorted_indices[1]==2))
    return data[5];
  else
    Assert (false, ExcInternalError());

  static double dummy_but_referenceable = 0;
  return dummy_but_referenceable;
}



template <>
inline
double &
SymmetricTensor<4,1>::operator () (const TableIndices<4> &indices)
{
  const unsigned int rank = 4;
  const unsigned int dim  = 1;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

  return data[0][0];
}



template <>
inline
double
SymmetricTensor<4,1>::operator () (const TableIndices<4> &indices) const
{
  const unsigned int rank = 4;
  const unsigned int dim  = 1;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

  return data[0][0];
}



template <>
inline
double &
SymmetricTensor<4,2>::operator () (const TableIndices<4> &indices)
{
  const unsigned int rank = 4;
  const unsigned int dim  = 2;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

				   // each entry of the tensor can be
				   // thought of as an entry in a
				   // matrix that maps the rolled-out
				   // rank-2 tensors into rolled-out
				   // rank-2 tensors. this is the
				   // format in which we store rank-4
				   // tensors. determine which
				   // position the present entry is
				   // stored in
  unsigned int base_index[2] ;
  if ((indices[0] == 0) && (indices[1] == 0))
    base_index[0] = 0;
  else if ((indices[0] == 1) && (indices[1] == 1))
    base_index[0] = 1;
  else
    base_index[0] = 2;

  if ((indices[2] == 0) && (indices[3] == 0))
    base_index[1] = 0;
  else if ((indices[2] == 1) && (indices[3] == 1))
    base_index[1] = 1;
  else
    base_index[1] = 2;  

  return data[base_index[0]][base_index[1]];
}



template <>
inline
double
SymmetricTensor<4,2>::operator () (const TableIndices<4> &indices) const
{
  const unsigned int rank = 4;
  const unsigned int dim  = 2;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

				   // each entry of the tensor can be
				   // thought of as an entry in a
				   // matrix that maps the rolled-out
				   // rank-2 tensors into rolled-out
				   // rank-2 tensors. this is the
				   // format in which we store rank-4
				   // tensors. determine which
				   // position the present entry is
				   // stored in
  unsigned int base_index[2] ;
  if ((indices[0] == 0) && (indices[1] == 0))
    base_index[0] = 0;
  else if ((indices[0] == 1) && (indices[1] == 1))
    base_index[0] = 1;
  else
    base_index[0] = 2;

  if ((indices[2] == 0) && (indices[3] == 0))
    base_index[1] = 0;
  else if ((indices[2] == 1) && (indices[3] == 1))
    base_index[1] = 1;
  else
    base_index[1] = 2;  

  return data[base_index[0]][base_index[1]];
}



template <>
inline
double &
SymmetricTensor<4,3>::operator () (const TableIndices<4> &indices)
{
  const unsigned int rank = 4;
  const unsigned int dim  = 3;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

				   // each entry of the tensor can be
				   // thought of as an entry in a
				   // matrix that maps the rolled-out
				   // rank-2 tensors into rolled-out
				   // rank-2 tensors. this is the
				   // format in which we store rank-4
				   // tensors. determine which
				   // position the present entry is
				   // stored in
  unsigned int base_index[2] ;
  if ((indices[0] == 0) && (indices[1] == 0))
    base_index[0] = 0;
  else if ((indices[0] == 1) && (indices[1] == 1))
    base_index[0] = 1;
  else if ((indices[0] == 2) && (indices[1] == 2))
    base_index[0] = 2;
  else if (((indices[0] == 0) && (indices[1] == 1)) ||
	   ((indices[0] == 1) && (indices[1] == 0)))
    base_index[0] = 3;
  else if (((indices[0] == 0) && (indices[1] == 2)) ||
	   ((indices[0] == 2) && (indices[1] == 0)))
    base_index[0] = 4;
  else
    {
      Assert (((indices[0] == 1) && (indices[1] == 2)) ||
	      ((indices[0] == 2) && (indices[1] == 1)),
	      ExcInternalError());
      base_index[0] = 5;
    }

  if ((indices[2] == 0) && (indices[3] == 0))
    base_index[1] = 0;
  else if ((indices[2] == 1) && (indices[3] == 1))
    base_index[1] = 1;
  else if ((indices[2] == 2) && (indices[3] == 2))
    base_index[1] = 2;
  else if (((indices[2] == 0) && (indices[3] == 1)) ||
	   ((indices[2] == 1) && (indices[3] == 0)))
    base_index[1] = 3;
  else if (((indices[2] == 0) && (indices[3] == 2)) ||
	   ((indices[2] == 2) && (indices[3] == 0)))
    base_index[1] = 4;
  else
    {
      Assert (((indices[2] == 1) && (indices[3] == 2)) ||
	      ((indices[2] == 2) && (indices[3] == 1)),
	      ExcInternalError());
      base_index[1] = 5;
    }
  
  return data[base_index[0]][base_index[1]];
}



template <>
inline
double
SymmetricTensor<4,3>::operator () (const TableIndices<4> &indices) const
{
  const unsigned int rank = 4;
  const unsigned int dim  = 3;
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dim, ExcIndexRange (indices[r], 0, dim));

				   // each entry of the tensor can be
				   // thought of as an entry in a
				   // matrix that maps the rolled-out
				   // rank-2 tensors into rolled-out
				   // rank-2 tensors. this is the
				   // format in which we store rank-4
				   // tensors. determine which
				   // position the present entry is
				   // stored in
  unsigned int base_index[2] ;
  if ((indices[0] == 0) && (indices[1] == 0))
    base_index[0] = 0;
  else if ((indices[0] == 1) && (indices[1] == 1))
    base_index[0] = 1;
  else if ((indices[0] == 2) && (indices[1] == 2))
    base_index[0] = 2;
  else if (((indices[0] == 0) && (indices[1] == 1)) ||
	   ((indices[0] == 1) && (indices[1] == 0)))
    base_index[0] = 3;
  else if (((indices[0] == 0) && (indices[1] == 2)) ||
	   ((indices[0] == 2) && (indices[1] == 0)))
    base_index[0] = 4;
  else
    {
      Assert (((indices[0] == 1) && (indices[1] == 2)) ||
	      ((indices[0] == 2) && (indices[1] == 1)),
	      ExcInternalError());
      base_index[0] = 5;
    }

  if ((indices[2] == 0) && (indices[3] == 0))
    base_index[1] = 0;
  else if ((indices[2] == 1) && (indices[3] == 1))
    base_index[1] = 1;
  else if ((indices[2] == 2) && (indices[3] == 2))
    base_index[1] = 2;
  else if (((indices[2] == 0) && (indices[3] == 1)) ||
	   ((indices[2] == 1) && (indices[3] == 0)))
    base_index[1] = 3;
  else if (((indices[2] == 0) && (indices[3] == 2)) ||
	   ((indices[2] == 2) && (indices[3] == 0)))
    base_index[1] = 4;
  else
    {
      Assert (((indices[2] == 1) && (indices[3] == 2)) ||
	      ((indices[2] == 2) && (indices[3] == 1)),
	      ExcInternalError());
      base_index[1] = 5;
    }
  
  return data[base_index[0]][base_index[1]];
}



template <int rank, int dim>
internal::SymmetricTensorAccessors::Accessor<rank,dim,true,rank-1>
SymmetricTensor<rank,dim>::operator [] (const unsigned int row) const
{
  return
    internal::SymmetricTensorAccessors::
    Accessor<rank,dim,true,rank-1> (*this, TableIndices<rank> (row));
}



template <int rank, int dim>
internal::SymmetricTensorAccessors::Accessor<rank,dim,false,rank-1>
SymmetricTensor<rank,dim>::operator [] (const unsigned int row)
{
  return
    internal::SymmetricTensorAccessors::
    Accessor<rank,dim,false,rank-1> (*this, TableIndices<rank> (row));
}



template <>
inline
double
SymmetricTensor<2,1>::norm () const
{
  return std::fabs(data[0]);
}



template <>
inline
double
SymmetricTensor<2,2>::norm () const
{
  return std::sqrt(data[0]*data[0] + data[1]*data[1] + 2*data[2]*data[2]);
}



template <>
inline
double
SymmetricTensor<2,3>::norm () const
{
  return std::sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2] +
                   2*data[3]*data[3] + 2*data[4]*data[4] + 2*data[5]*data[5]);
}



template <>
inline
double
SymmetricTensor<4,1>::norm () const
{
  return std::fabs(data[0][0]);
}



template <>
inline
double
SymmetricTensor<4,2>::norm () const
{
  const unsigned int dim = 2;

				   // this is not really efficient and
				   // could be improved by counting
				   // how often each tensor entry is
				   // accessed, but this isn't a
				   // really frequent operation anyway
  double t = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  {
	    const double a = (*this)[i][j][k][l];
	    t += a * a;
	  }
  return std::sqrt(t);
}



template <>
inline
double
SymmetricTensor<4,3>::norm () const
{
  const unsigned int dim = 3;

				   // this is not really efficient and
				   // could be improved by counting
				   // how often each tensor entry is
				   // accessed, but this isn't a
				   // really frequent operation anyway
  double t = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  {
	    const double a = (*this)[i][j][k][l];
	    t += a * a;
	  }
  return std::sqrt(t);
}

///@endif

/* ----------------- Non-member functions operating on tensors. ------------ */

/**
 * Compute the determinant of a tensor or rank 2. The determinant is
 * also commonly referred to as the third invariant of rank-2 tensors.
 *
 * For the present case of one-dimensional tensors, the determinant
 * equals the only element and is therefore equivalent to the trace.
 *
 * For greater notational simplicity, there is also a <tt>third_invariant</tt>
 * function that returns the determinant of a tensor.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
inline
double determinant (const SymmetricTensor<2,1> &t)
{
  return t.data[0];
}



/**
 * Compute the determinant of a tensor or rank 2. The determinant is
 * also commonly referred to as the third invariant of rank-2 tensors.
 *
 * For greater notational simplicity, there is also a <tt>third_invariant</tt>
 * function that returns the determinant of a tensor.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
inline
double determinant (const SymmetricTensor<2,2> &t)
{
  return (t.data[0] * t.data[1] - t.data[2]*t.data[2]);
}




/**
 * Compute the determinant of a tensor or rank 2. The determinant is
 * also commonly referred to as the third invariant of rank-2 tensors.
 *
 * For greater notational simplicity, there is also a <tt>third_invariant</tt>
 * function that returns the determinant of a tensor.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
inline
double determinant (const SymmetricTensor<2,3> &t)
{
				   // in analogy to general tensors, but
				   // there's something to be simplified for
				   // the present case
  return ( t.data[0]*t.data[1]*t.data[2]
	   -t.data[0]*t.data[5]*t.data[5]
	   -t.data[1]*t.data[4]*t.data[4]
	   -t.data[2]*t.data[3]*t.data[3]
	   +2*t.data[3]*t.data[4]*t.data[5] );
}



/**
 * Compute the determinant of a tensor or rank 2. This function therefore
 * computes the same value as the <tt>determinant()</tt> functions and is only
 * provided for greater notational simplicity (since there are also functions
 * <tt>first_invariant</tt> and <tt>second_invariant</tt>).
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
double third_invariant (const SymmetricTensor<2,dim> &t)
{
  return determinant (t);
}



/**
 * Compute and return the trace of a tensor of rank 2, i.e. the sum of
 * its diagonal entries. The trace is the first invariant of a rank-2
 * tensor.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
double trace (const SymmetricTensor<2,dim> &d)
{
  double t=0;
  for (unsigned int i=0; i<dim; ++i)
    t += d.data[i];
  return t;
}


/**
 * Compute the trace of a tensor or rank 2. This function therefore
 * computes the same value as the <tt>trace()</tt> functions and is only
 * provided for greater notational simplicity (since there are also functions
 * <tt>second_invariant</tt> and <tt>third_invariant</tt>).
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
double first_invariant (const SymmetricTensor<2,dim> &t)
{
  return trace (t);
}


/**
 * Compute the second invariant a tensor or rank 2. The second invariant is
 * defined as <tt>I2 = 1/2[ (trace sigma)^2 - trace (sigma^2) ]</tt>.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
double second_invariant (const SymmetricTensor<2,dim> &t)
{
  return (trace(t)*trace(t) - t*t) / 2;
}



/**
 * Return the transpose of the given symmetric tensor. Since we are working
 * with symmetric objects, the transpose is of course the same as the original
 * tensor. This function mainly exists for compatibility with the Tensor
 * class.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
transpose (const SymmetricTensor<rank,dim> &t)
{
  return t;
}



/**
 * Compute the deviator of a symmetric tensor, which is defined as <tt>dev[s]
 * = s - 1/dim*tr[s]*I</tt>, where <tt>I</tt> is the identity operator. This
 * quantity equals the original tensor minus its contractive or dilative
 * component and refers to the shear in, for example, elasticity.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
SymmetricTensor<2,dim>
deviator (const SymmetricTensor<2,dim> &t)
{
  SymmetricTensor<2,dim> tmp = t;

                                   // subtract scaled trace from the diagonal
  const double tr = trace(t) / dim;
  for (unsigned int i=0; i<dim; ++i)
    tmp.data[i] -= tr;
  
  return tmp;
}



/**
 * Return a unit symmetric tensor of rank 2 and dimension 1.
 * 
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <>
inline
SymmetricTensor<2,1>
unit_symmetric_tensor<1> () 
{
  SymmetricTensor<2,1> tmp;
  tmp.data[0] = 1;
  return tmp;
}



/**
 * Return a unit symmetric tensor of rank 2 and dimension 2.
 * 
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <>
inline
SymmetricTensor<2,2>
unit_symmetric_tensor<2> () 
{
  SymmetricTensor<2,2> tmp;
  tmp.data[0] = tmp.data[1] = 1;
  return tmp;
}



/**
 * Return a unit symmetric tensor of rank 2 and dimension 3.
 * 
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <>
inline
SymmetricTensor<2,3>
unit_symmetric_tensor<3> () 
{
  SymmetricTensor<2,3> tmp;
  tmp.data[0] = tmp.data[1] = tmp.data[2] = 1;
  return tmp;
}




/**
 * Return the tensor of rank 4 that, when multiplied by a symmetric rank 2
 * tensor <tt>t</tt> returns the deviator <tt>dev t</tt>. It is the operator
 * representation of the linear deviator operator.
 *
 * For every tensor <tt>t</tt>, there holds the identity
 * <tt>deviator(t)==deviator_tensor&lt;dim&gt;()*t</tt>, up to numerical
 * round-off. The reason this operator representation is provided is that one
 * sometimes needs to invert operators like <tt>identity_tensor&lt;dim&gt;() +
 * delta_t*deviator_tensor&lt;dim&gt;()</tt> or similar.
 * 
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
SymmetricTensor<4,dim>
deviator_tensor () 
{
  SymmetricTensor<4,dim> tmp;
  
                                   // fill the elements treating the diagonal
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      tmp.data[i][j] = (i==j ? 1 : 0) - 1./dim;

                                   // then fill the ones that copy over the
                                   // non-diagonal elements. note that during
                                   // the double-contraction, we handle the
                                   // off-diagonal elements twice, so simply
                                   // copying requires a weight of 1/2
  for (unsigned int i=dim;
       i<internal::SymmetricTensorAccessors::StorageType<4,dim>::n_rank2_components;
       ++i)
    tmp.data[i][i] = 0.5;
  
  return tmp;
}



/**
 * Return the tensor of rank 4 that, when multiplied by a symmetric rank 2
 * tensor <tt>t</tt> returns the deviator <tt>dev t</tt>. It is the operator
 * representation of the linear deviator operator.
 *
 * Note that this tensor, even though it is the identity, has a somewhat funny
 * form, and in particular does not only consist of zeros and ones. For
 * example, for <tt>dim=2</tt>, the identity tensor has all zero entries
 * except for <tt>id[0][0][0][0]=id[1][1][1][1]=1</tt> and
 * <tt>id[0][1][0][1]=id[0][1][1][0]=id[1][0][0][1]=id[1][0][1][0]=1/2</tt>. To
 * see why this factor of 1/2 is necessary, consider computing <tt>A=Id
 * . B</tt>. For the element <tt>a_01</tt> we have <tt>a_01=id_0100 b_00 +
 * id_0111 b_11 + id_0101 b_01 + id_0110 b_10</tt>. On the other hand, we need
 * to have <tt>a_01=b_01</tt>, and symmetry implies <tt>b_01=b_10</tt>,
 * leading to <tt>a_01=(id_0101+id_0110) b_01</tt>, or, again by symmetry,
 * <tt>id_0101=id_0110=1/2</tt>. Similar considerations hold for the
 * three-dimensional case.
 * 
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
SymmetricTensor<4,dim>
identity_tensor () 
{
  SymmetricTensor<4,dim> tmp;
  
                                   // fill the elements treating the diagonal
  for (unsigned int i=0; i<dim; ++i)
    tmp.data[i][i] = 1;

                                   // then fill the ones that copy over the
                                   // non-diagonal elements. note that during
                                   // the double-contraction, we handle the
                                   // off-diagonal elements twice, so simply
                                   // copying requires a weight of 1/2
  for (unsigned int i=dim;
       i<internal::SymmetricTensorAccessors::StorageType<4,dim>::n_rank2_components;
       ++i)
    tmp.data[i][i] = 0.5;
  
  return tmp;
}



/**
 * Invert a symmetric rank-4 tensor. Since symmetric rank-4 tensors are
 * mappings from and to symmetric rank-2 tensors, they can have an
 * inverse. This function computes it, if it exists, for the case that the
 * dimension equals 1.
 *
 * If a tensor is not invertible, then the result is unspecified, but will
 * likely contain the results of a division by zero or a very small number at
 * the very least.
 * 
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <>
inline
SymmetricTensor<4,1>
invert (const SymmetricTensor<4,1> &t)
{
  SymmetricTensor<4,1> tmp;
  tmp.data[0][0] = 1./t.data[0][0];
  return tmp;
}



/**
 * Invert a symmetric rank-4 tensor. Since symmetric rank-4 tensors are
 * mappings from and to symmetric rank-2 tensors, they can have an
 * inverse. This function computes it, if it exists, for the case that the
 * dimension equals 2.
 *
 * If a tensor is not invertible, then the result is unspecified, but will
 * likely contain the results of a division by zero or a very small number at
 * the very least.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <>
inline
SymmetricTensor<4,2>
invert (const SymmetricTensor<4,2> &t)
{
  SymmetricTensor<4,2> tmp;

                                   // inverting this tensor is a little more
                                   // complicated than necessary, since we
                                   // store the data of 't' as a 3x3 matrix
                                   // t.data, but the product between a rank-4
                                   // and a rank-2 tensor is really not the
                                   // product between this matrix and the
                                   // 3-vector of a rhs, but rather
                                   //
                                   // B.vec = t.data * mult * A.vec
                                   //
                                   // where mult is a 3x3 matrix with
                                   // entries [[1,0,0],[0,1,0],[0,0,2]] to
                                   // capture the fact that we need to add up
                                   // both the c_ij12*a_12 and the c_ij21*a_21
                                   // terms
                                   //
                                   // in addition, in this scheme, the
                                   // identity tensor has the matrix
                                   // representation mult^-1.
                                   //
                                   // the inverse of 't' therefore has the
                                   // matrix representation
                                   //
                                   // inv.data = mult^-1 * t.data^-1 * mult^-1
                                   //
                                   // in order to compute it, let's first
                                   // compute the inverse of t.data and put it
                                   // into tmp.data; at the end of the
                                   // function we then scale the last row and
                                   // column of the inverse by 1/2,
                                   // corresponding to the left and right
                                   // multiplication with mult^-1
  const double t4 = t.data[0][0]*t.data[1][1],
               t6 = t.data[0][0]*t.data[1][2],
               t8 = t.data[0][1]*t.data[1][0],
              t00 = t.data[0][2]*t.data[1][0],
              t01 = t.data[0][1]*t.data[2][0],
              t04 = t.data[0][2]*t.data[2][0],
              t07 = 1.0/(t4*t.data[2][2]-t6*t.data[2][1]-t8*t.data[2][2]+
                         t00*t.data[2][1]+t01*t.data[1][2]-t04*t.data[1][1]);
  tmp.data[0][0] = (t.data[1][1]*t.data[2][2]-t.data[1][2]*t.data[2][1])*t07;
  tmp.data[0][1] = -(t.data[0][1]*t.data[2][2]-t.data[0][2]*t.data[2][1])*t07;
  tmp.data[0][2] = -(-t.data[0][1]*t.data[1][2]+t.data[0][2]*t.data[1][1])*t07;
  tmp.data[1][0] = -(t.data[1][0]*t.data[2][2]-t.data[1][2]*t.data[2][0])*t07;
  tmp.data[1][1] = (t.data[0][0]*t.data[2][2]-t04)*t07;
  tmp.data[1][2] = -(t6-t00)*t07;
  tmp.data[2][0] = -(-t.data[1][0]*t.data[2][1]+t.data[1][1]*t.data[2][0])*t07;
  tmp.data[2][1] = -(t.data[0][0]*t.data[2][1]-t01)*t07;
  tmp.data[2][2] = (t4-t8)*t07;

                                   // scale last row and column as mentioned
                                   // above
  tmp.data[2][0] /= 2;
  tmp.data[2][1] /= 2;
  tmp.data[0][2] /= 2;
  tmp.data[1][2] /= 2;
  tmp.data[2][2] /= 4;
  
  return tmp;
}



/**
 * Invert a symmetric rank-4 tensor. Since symmetric rank-4 tensors are
 * mappings from and to symmetric rank-2 tensors, they can have an
 * inverse. This function computes it, if it exists, for the case that the
 * dimension equals 3.
 *
 * If a tensor is not invertible, then the result is unspecified, but will
 * likely contain the results of a division by zero or a very small number at
 * the very least.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <>
SymmetricTensor<4,2>
invert (const SymmetricTensor<4,2> &t);
// this function is implemented in the .cc file


/**
 * Return the tensor of rank 4 that is the outer product of the two tensors
 * given as arguments, i.e. the result $T=t1 \otimes t2$ satisfies
 * <tt>T phi = t1 (t2 : phi)</tt> for all symmetric tensors <tt>phi</tt>.
 *
 * For example, the deviator tensor can be computed as
 * <tt>identity_tensor<dim>() -
 * 1/d*outer_product(unit_symmetric_tensor<dim>(),
 * unit_symmetric_tensor<dim>())</tt>, since the (double) contraction with the
 * unit tensor yields the trace of a symmetric tensor.
 * 
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
SymmetricTensor<4,dim>
outer_product (const SymmetricTensor<2,dim> &t1,
               const SymmetricTensor<2,dim> &t2) 
{
  SymmetricTensor<4,dim> tmp;

                                   // fill only the elements really needed
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=k; l<dim; ++l)
          tmp[i][j][k][l] = t1[i][j] * t2[k][l];

  return tmp;
}



/**
 * Return the symmetrized version of the full rank-2 tensor,
 * i.e. (t+transpose(t))/2, as a symmetric rank-2 tensor.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
SymmetricTensor<2,dim>
symmetrize (const Tensor<2,dim> &t)
{
  SymmetricTensor<2,dim> s;
  for (unsigned int i=0; i<dim; ++i)
    {
      s(TableIndices<2>(i,i)) = t[i][i];
  
      for (unsigned int j=i+1; j<dim; ++j)
	s(TableIndices<2>(i,j))
	  = (t[i][j] + t[j][i])/2;
    }
  return s;
}


/**
 * Multiplication of a symmetric tensor of general rank with a scalar double
 * from the right.
 *
 * @relates SymmetricTensor
 */
template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
operator * (const SymmetricTensor<rank,dim> &t,
	    const double            factor)
{
  SymmetricTensor<rank,dim> tt = t;
  tt *= factor;
  return tt;
}



/**
 * Multiplication of a symmetric tensor of general rank with a scalar double
 * from the left.
 *
 * @relates SymmetricTensor
 */
template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
operator * (const double            factor,
	    const SymmetricTensor<rank,dim> &t)
{
  SymmetricTensor<rank,dim> tt = t;
  tt *= factor;
  return tt;
}



/**
 * Division of a symmetric tensor of general rank by a scalar double.
 *
 * @relates SymmetricTensor
 */
template <int rank, int dim>
inline
SymmetricTensor<rank,dim>
operator / (const SymmetricTensor<rank,dim> &t,
	    const double            factor)
{
  SymmetricTensor<rank,dim> tt = t;
  tt /= factor;
  return tt;
}




/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in the symmetric tensor of rank 2 that is given as first
 * argument to this function. This operation is the symmetric tensor
 * analogon of a matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the
 * SymmetricTensor class. It should not be used, however, since the
 * member operator has knowledge of the actual data storage format and
 * is at least 2 orders of magnitude faster. This function mostly
 * exists for compatibility purposes with the general tensor class.
 *
 * @related SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
inline
void
double_contract (SymmetricTensor<2,1> &tmp,
		 const SymmetricTensor<4,1> &t,
		 const SymmetricTensor<2,1> &s)
{
  tmp[0][0] = t[0][0][0][0] * s[0][0];
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in the symmetric tensor of rank 2 that is given as first
 * argument to this function. This operation is the symmetric tensor
 * analogon of a matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the
 * SymmetricTensor class. It should not be used, however, since the
 * member operator has knowledge of the actual data storage format and
 * is at least 2 orders of magnitude faster. This function mostly
 * exists for compatibility purposes with the general tensor class.
 *
 * @related SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
inline
void
double_contract (SymmetricTensor<2,1> &tmp,
		 const SymmetricTensor<2,1> &s,
		 const SymmetricTensor<4,1> &t)
{
  tmp[0][0] = t[0][0][0][0] * s[0][0];
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in the symmetric tensor of rank 2 that is given as first
 * argument to this function. This operation is the symmetric tensor
 * analogon of a matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the
 * SymmetricTensor class. It should not be used, however, since the
 * member operator has knowledge of the actual data storage format and
 * is at least 2 orders of magnitude faster. This function mostly
 * exists for compatibility purposes with the general tensor class.
 *
 * @related SymmetricTensor @author Wolfgang Bangerth, 2005
 */
inline
void
double_contract (SymmetricTensor<2,2> &tmp,
		 const SymmetricTensor<4,2> &t,
		 const SymmetricTensor<2,2> &s)
{
  const unsigned int dim = 2;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = t[i][j][0][0] * s[0][0] +
		  t[i][j][1][1] * s[1][1] +
		  2 * t[i][j][0][1] * s[0][1];
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in the symmetric tensor of rank 2 that is given as first
 * argument to this function. This operation is the symmetric tensor
 * analogon of a matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the
 * SymmetricTensor class. It should not be used, however, since the
 * member operator has knowledge of the actual data storage format and
 * is at least 2 orders of magnitude faster. This function mostly
 * exists for compatibility purposes with the general tensor class.
 *
 * @related SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
inline
void
double_contract (SymmetricTensor<2,2> &tmp,
		 const SymmetricTensor<2,2> &s,
		 const SymmetricTensor<4,2> &t)
{
  const unsigned int dim = 2;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = s[0][0] * t[0][0][i][j] * +
		  s[1][1] * t[1][1][i][j] +
		  2 * s[0][1] * t[0][1][i][j];
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in the symmetric tensor of rank 2 that is given as first
 * argument to this function. This operation is the symmetric tensor
 * analogon of a matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the
 * SymmetricTensor class. It should not be used, however, since the
 * member operator has knowledge of the actual data storage format and
 * is at least 2 orders of magnitude faster. This function mostly
 * exists for compatibility purposes with the general tensor class.
 *
 * @related SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
inline
void
double_contract (SymmetricTensor<2,3> &tmp,
		 const SymmetricTensor<4,3> &t,
		 const SymmetricTensor<2,3> &s)
{
  const unsigned int dim = 3;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = t[i][j][0][0] * s[0][0] +
		  t[i][j][1][1] * s[1][1] +
		  t[i][j][2][2] * s[2][2] +
		  2 * t[i][j][0][1] * s[0][1] +
		  2 * t[i][j][0][2] * s[0][2] +
		  2 * t[i][j][1][2] * s[1][2];
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in the symmetric tensor of rank 2 that is given as first
 * argument to this function. This operation is the symmetric tensor
 * analogon of a matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the
 * SymmetricTensor class. It should not be used, however, since the
 * member operator has knowledge of the actual data storage format and
 * is at least 2 orders of magnitude faster. This function mostly
 * exists for compatibility purposes with the general tensor class.
 *
 * @related SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
inline
void
double_contract (SymmetricTensor<2,3> &tmp,
		 const SymmetricTensor<2,3> &s,
		 const SymmetricTensor<4,3> &t)
{
  const unsigned int dim = 3;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = s[0][0] * t[0][0][i][j] +
		  s[1][1] * t[1][1][i][j] +
		  s[2][2] * t[2][2][i][j] +
		  2 * s[0][1] * t[0][1][i][j] +
		  2 * s[0][2] * t[0][2][i][j] +
		  2 * s[1][2] * t[1][2][i][j];
}



#endif
