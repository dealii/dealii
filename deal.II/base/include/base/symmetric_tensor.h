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
 * Tensors of rank 4 are considered symmetric if they are operators mapping
 * symmetric rank-2 tensors onto symmetric rank-2 tensors. This entails
 * certain symmetry properties on the elements in their 4-dimensional index
 * space.
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
                                      * that @p t is already symmetric, but
                                      * this is not checked: we simply copy
                                      * only a subset of elements.
                                      */
    SymmetricTensor (const Tensor<2,dim> &t);

				     /**
				      *  Assignment operator.
				      */
    SymmetricTensor & operator = (const SymmetricTensor &);

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
                                      * Scalar product between two symmetric
                                      * tensors. It is the contraction
                                      * <tt>a<sub>ij</sub>b<sub>ij</sub></tt>
                                      * over all indices <tt>i,j</tt>. While
                                      * it is possible to define other scalar
                                      * products (and associated induced
                                      * norms), this one seems to be the most
                                      * appropriate one.
                                      */
    double operator * (const SymmetricTensor &s) const;
    
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
                                      * above for two symmetric tensors.
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
                                      * Data storage for a symmetric tensor.
                                      */
    typename internal::SymmetricTensorAccessors::StorageType<rank,dim>::base_tensor_type data;
};



// ------------------------- inline functions ------------------------

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
double
SymmetricTensor<2,1>::operator * (const SymmetricTensor<2,1> &s) const
{
  return data[0] * s.data[0];
}



template <>
double
SymmetricTensor<2,2>::operator * (const SymmetricTensor<2,2> &s) const
{
  return (data[0] * s.data[0] +
          data[1] * s.data[1] +
          2*data[2] * s.data[2]);
}



template <>
double
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
double
SymmetricTensor<4,1>::operator * (const SymmetricTensor<4,1> &s) const
{
  return data[0][0] * s.data[0][0];
}



template <>
double
SymmetricTensor<4,2>::operator * (const SymmetricTensor<4,2> &s) const
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
	  t += (*this)[i][j][k][l] * s[i][j][k][l];
  return t;
}



template <>
double
SymmetricTensor<4,3>::operator * (const SymmetricTensor<4,3> &s) const
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
	  t += (*this)[i][j][k][l] * s[i][j][k][l];
  return t;
}



template <>
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
double
SymmetricTensor<2,1>::norm () const
{
  return std::fabs(data[0]);
}



template <>
double
SymmetricTensor<2,2>::norm () const
{
  return std::sqrt(data[0]*data[0] + data[1]*data[1] + 2*data[2]*data[2]);
}



template <>
double
SymmetricTensor<2,3>::norm () const
{
  return std::sqrt(data[0]*data[0] + data[1]*data[1] + data[2]*data[2] +
                   2*data[3]*data[3] + 2*data[4]*data[4] + 2*data[5]*data[5]);
}



template <>
double
SymmetricTensor<4,1>::norm () const
{
  return std::fabs(data[0][0]);
}



template <>
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


/* ----------------- Non-member functions operating on tensors. ------------ */

/**
 * Compute the determinant of a tensor or rank 2, here for <tt>dim==2</tt>.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
inline
double determinant (const SymmetricTensor<2,2> &t)
{
  return (t[0][0] * t[1][1] - 2*t[0][1]*t[0][1]);
}




/**
 * Compute the determinant of a tensor or rank 2, here for <tt>dim==3</tt>.
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
  return ( t[0][0]*t[1][1]*t[2][2]
	   -t[0][0]*t[1][2]*t[1][2]
	   -t[1][1]*t[0][2]*t[0][2]
	   -t[2][2]*t[0][1]*t[0][1]
	   +2*t[0][1]*t[0][2]*t[1][2] );
}



/**
 * Compute and return the trace of a tensor of rank 2, i.e. the sum of
 * its diagonal entries.
 *
 * @relates SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int rank, int dim>
double trace (const SymmetricTensor<rank,dim> &d)
{
  double t=0;
  for (unsigned int i=0; i<dim; ++i)
    t += d[i][i];
  return t;
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
 * resulting in a symmetric tensor of rank 2. This operation is the
 * symmetric tensor analogon of a matrix-vector multiplication.
 *
 * @related SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
SymmetricTensor<2,1>
operator * (const SymmetricTensor<4,1> &t,
	    const SymmetricTensor<2,1> &s)
{
  const unsigned int dim = 1;
  SymmetricTensor<2,dim> tmp;
  tmp[0][0] = t[0][0][0][0] * s[0][0];
  return tmp;
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in a symmetric tensor of rank 2. This operation is the
 * symmetric tensor analogon of a matrix-vector multiplication.
 *
 * @related SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
SymmetricTensor<2,1>
operator * (const SymmetricTensor<2,1> &s,
	    const SymmetricTensor<4,1> &t)
{
  const unsigned int dim = 1;
  SymmetricTensor<2,dim> tmp;
  tmp[0][0] = t[0][0][0][0] * s[0][0];
  return tmp;
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in a symmetric tensor of rank 2. This operation is the
 * symmetric tensor analogon of a matrix-vector multiplication.
 *
 * @related SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
SymmetricTensor<2,2>
operator * (const SymmetricTensor<4,2> &t,
	    const SymmetricTensor<2,2> &s)
{
  const unsigned int dim = 2;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = t[i][j][0][0] * s[0][0] +
		  t[i][j][1][1] * s[1][1] +
		  2 * t[i][j][0][1] * s[0][1];

  return tmp;
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in a symmetric tensor of rank 2. This operation is the
 * symmetric tensor analogon of a matrix-vector multiplication.
 *
 * @related SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
SymmetricTensor<2,2>
operator * (const SymmetricTensor<2,2> &s,
	    const SymmetricTensor<4,2> &t)
{
  const unsigned int dim = 2;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = s[0][0] * t[0][0][i][j] * +
		  s[1][1] * t[1][1][i][j] +
		  2 * s[0][1] * t[0][1][i][j];

  return tmp;
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in a symmetric tensor of rank 2. This operation is the
 * symmetric tensor analogon of a matrix-vector multiplication.
 *
 * @related SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
SymmetricTensor<2,3>
operator * (const SymmetricTensor<4,3> &t,
	    const SymmetricTensor<2,3> &s)
{
  const unsigned int dim = 3;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = t[i][j][0][0] * s[0][0] +
		  t[i][j][1][1] * s[1][1] +
		  t[i][j][2][2] * s[2][2] +
		  2 * t[i][j][0][1] * s[0][1] +
		  2 * t[i][j][0][2] * s[0][2] +
		  2 * t[i][j][1][2] * s[1][2];

  return tmp;
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in a symmetric tensor of rank 2. This operation is the
 * symmetric tensor analogon of a matrix-vector multiplication.
 *
 * @related SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
SymmetricTensor<2,3>
operator * (const SymmetricTensor<2,3> &s,
	    const SymmetricTensor<4,3> &t)
{
  const unsigned int dim = 3;
  SymmetricTensor<2,dim> tmp;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      tmp[i][j] = s[0][0] * t[0][0][i][j] +
		  s[1][1] * t[1][1][i][j] +
		  s[2][2] * t[2][2][i][j] +
		  2 * s[0][1] * t[0][1][i][j] +
		  2 * s[0][2] * t[0][2][i][j] +
		  2 * s[1][2] * t[1][2][i][j];

  return tmp;
}




#endif
