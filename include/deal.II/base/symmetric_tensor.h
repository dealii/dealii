// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
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

#ifndef dealii_symmetric_tensor_h
#define dealii_symmetric_tensor_h


#include <deal.II/base/tensor.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/template_constraints.h>

#include <array>
#include <algorithm>
#include <functional>

DEAL_II_NAMESPACE_OPEN

template <int rank, int dim, typename Number=double> class SymmetricTensor;

template <int dim, typename Number>
SymmetricTensor<2,dim,Number> unit_symmetric_tensor ();

template <int dim, typename Number>
SymmetricTensor<4,dim,Number> deviator_tensor ();

template <int dim, typename Number>
SymmetricTensor<4,dim,Number> identity_tensor ();

template <int dim, typename Number>
SymmetricTensor<2,dim,Number> invert (const SymmetricTensor<2,dim,Number> &);

template <int dim, typename Number>
SymmetricTensor<4,dim,Number> invert (const SymmetricTensor<4,dim,Number> &);

template <int dim2, typename Number> Number
trace (const SymmetricTensor<2,dim2,Number> &);

template <int dim, typename Number>
SymmetricTensor<2,dim,Number> deviator (const SymmetricTensor<2,dim,Number> &);

template <int dim, typename Number>
Number determinant (const SymmetricTensor<2,dim,Number> &);



namespace internal
{
  /**
   * A namespace for functions and classes that are internal to how the
   * SymmetricTensor class (and its associate functions) works.
   */
  namespace SymmetricTensorImplementation
  {
    /**
     * Compute the inverse of a symmetric tensor of a
     * generic @p rank, @p dim and @p Number type.
     */
    template <int rank, int dim, typename Number>
    struct Inverse;
  }

  /**
   * A namespace for classes that are internal to how the SymmetricTensor
   * class works.
   */
  namespace SymmetricTensorAccessors
  {
    /**
     * Create a TableIndices<2> object where the first entries up to
     * <tt>position-1</tt> are taken from previous_indices, and new_index is
     * put at position <tt>position</tt>. The remaining indices remain in
     * invalid state.
     */
    inline
    TableIndices<2> merge (const TableIndices<2> &previous_indices,
                           const unsigned int     new_index,
                           const unsigned int     position)
    {
      Assert (position < 2, ExcIndexRange (position, 0, 2));

      if (position == 0)
        return TableIndices<2>(new_index, numbers::invalid_unsigned_int);
      else
        return TableIndices<2>(previous_indices[0], new_index);
    }



    /**
     * Create a TableIndices<4> object where the first entries up to
     * <tt>position-1</tt> are taken from previous_indices, and new_index is
     * put at position <tt>position</tt>. The remaining indices remain in
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
          return TableIndices<4>(new_index,
                                 numbers::invalid_unsigned_int,
                                 numbers::invalid_unsigned_int,
                                 numbers::invalid_unsigned_int);
        case 1:
          return TableIndices<4>(previous_indices[0],
                                 new_index,
                                 numbers::invalid_unsigned_int,
                                 numbers::invalid_unsigned_int);
        case 2:
          return TableIndices<4>(previous_indices[0],
                                 previous_indices[1],
                                 new_index,
                                 numbers::invalid_unsigned_int);
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
     * Typedef template magic denoting the result of a double contraction
     * between two tensors or ranks rank1 and rank2. In general, this is a
     * tensor of rank <tt>rank1+rank2-4</tt>, but if this is zero it is a
     * single scalar Number. For this case, we have a specialization.
     *
     * @author Wolfgang Bangerth, 2005, Jean-Paul Pelteret, 2017
     */
    template <int rank1, int rank2, int dim, typename Number, typename OtherNumber = Number>
    struct double_contraction_result
    {
      typedef typename ProductType<Number,OtherNumber>::type          value_type;
      typedef ::dealii::SymmetricTensor<rank1+rank2-4,dim,value_type> type;
    };


    /**
     * Typedef template magic denoting the result of a double contraction
     * between two tensors or ranks rank1 and rank2. In general, this is a
     * tensor of rank <tt>rank1+rank2-4</tt>, but if this is zero it is a
     * single scalar Number. For this case, we have a specialization.
     *
     * @author Wolfgang Bangerth, 2005, Jean-Paul Pelteret, 2017
     */
    template <int dim, typename Number, typename OtherNumber>
    struct double_contraction_result<2,2,dim,Number,OtherNumber>
    {
      typedef typename ProductType<Number,OtherNumber>::type type;
    };



    /**
     * Declaration of typedefs for the type of data structures which are used
     * to store symmetric tensors. For example, for rank-2 symmetric tensors,
     * we use a flat vector to store all the elements. On the other hand,
     * symmetric rank-4 tensors are mappings from symmetric rank-2 tensors
     * into symmetric rank-2 tensors, so they can be represented as matrices,
     * etc.
     *
     * This information is probably of little interest to all except the
     * accessor classes that need it. In particular, you shouldn't make any
     * assumptions about the storage format in your application programs.
     */
    template <int rank, int dim, typename Number>
    struct StorageType;

    /**
     * Specialization of StorageType for rank-2 tensors.
     */
    template <int dim, typename Number>
    struct StorageType<2,dim,Number>
    {
      /**
       * Number of independent components of a symmetric tensor of rank 2. We
       * store only the upper right half of it.
       */
      static const unsigned int
      n_independent_components = (dim*dim + dim)/2;

      /**
       * Declare the type in which we actually store the data.
       */
      typedef Tensor<1,n_independent_components,Number> base_tensor_type;
    };



    /**
     * Specialization of StorageType for rank-4 tensors.
     */
    template <int dim, typename Number>
    struct StorageType<4,dim,Number>
    {
      /**
       * Number of independent components of a symmetric tensor of rank 2.
       * Since rank-4 tensors are mappings between such objects, we need this
       * information.
       */
      static const unsigned int
      n_rank2_components = (dim*dim + dim)/2;

      /**
       * Number of independent components of a symmetric tensor of rank 4.
       */
      static const unsigned int
      n_independent_components = (n_rank2_components *
                                  StorageType<2,dim,Number>::n_independent_components);

      /**
       * Declare the type in which we actually store the data. Symmetric
       * rank-4 tensors are mappings between symmetric rank-2 tensors, so we
       * can represent the data as a matrix if we represent the rank-2 tensors
       * as vectors.
       */
      typedef Tensor<2,n_rank2_components,Number> base_tensor_type;
    };



    /**
     * Switch type to select a tensor of rank 2 and dimension <tt>dim</tt>,
     * switching on whether the tensor should be constant or not.
     */
    template <int rank, int dim, bool constness, typename Number>
    struct AccessorTypes;

    /**
     * Switch type to select a tensor of rank 2 and dimension <tt>dim</tt>,
     * switching on whether the tensor should be constant or not.
     *
     * Specialization for constant tensors.
     */
    template <int rank, int dim, typename Number>
    struct AccessorTypes<rank,dim,true,Number>
    {
      typedef const ::dealii::SymmetricTensor<rank,dim,Number> tensor_type;

      typedef Number reference;
    };

    /**
     * Switch type to select a tensor of rank 2 and dimension <tt>dim</tt>,
     * switching on whether the tensor should be constant or not.
     *
     * Specialization for non-constant tensors.
     */
    template <int rank, int dim, typename Number>
    struct AccessorTypes<rank,dim,false,Number>
    {
      typedef ::dealii::SymmetricTensor<rank,dim,Number> tensor_type;

      typedef Number &reference;
    };


    /**
     * @internal
     *
     * Class that acts as accessor to elements of type SymmetricTensor. The
     * template parameter <tt>C</tt> may be either true or false, and
     * indicates whether the objects worked on are constant or not (i.e. write
     * access is only allowed if the value is false).
     *
     * Since with <tt>N</tt> indices, the effect of applying
     * <tt>operator[]</tt> is getting access to something we <tt>N-1</tt>
     * indices, we have to implement these accessor classes recursively, with
     * stopping when we have only one index left. For the latter case, a
     * specialization of this class is declared below, where calling
     * <tt>operator[]</tt> gives you access to the objects actually stored by
     * the tensor; the tensor class also makes sure that only those elements
     * are actually accessed which we actually store, i.e. it reorders indices
     * if necessary. The template parameter <tt>P</tt> indicates how many
     * remaining indices there are. For a rank-2 tensor, <tt>P</tt> may be
     * two, and when using <tt>operator[]</tt>, an object with <tt>P=1</tt>
     * emerges.
     *
     * As stated for the entire namespace, you will not usually have to do
     * with these classes directly, and should not try to use their interface
     * directly as it may change without notice. In fact, since the
     * constructors are made private, you will not even be able to generate
     * objects of this class, as they are only thought as temporaries for
     * access to elements of the table class, not for passing them around as
     * arguments of functions, etc.
     *
     * This class is an adaptation of a similar class used for the Table
     * class.
     *
     * @author Wolfgang Bangerth, 2002, 2005
     */
    template <int rank, int dim, bool constness, int P, typename Number>
    class Accessor
    {
    public:
      /**
       * Import two typedefs from the switch class above.
       */
      typedef typename AccessorTypes<rank,dim,constness,Number>::reference reference;
      typedef typename AccessorTypes<rank,dim,constness,Number>::tensor_type tensor_type;

    private:
      /**
       * Constructor. Take a reference to the tensor object which we will
       * access.
       *
       * The second argument denotes the values of previous indices into the
       * tensor. For example, for a rank-4 tensor, if P=2, then we will
       * already have had two successive element selections (e.g. through
       * <tt>tensor[1][2]</tt>), and the two index values have to be stored
       * somewhere. This class therefore only makes use of the first rank-P
       * elements of this array, but passes it on to the next level with P-1
       * which fills the next entry, and so on.
       *
       * The constructor is made private in order to prevent you having such
       * objects around. The only way to create such objects is via the
       * <tt>Table</tt> class, which only generates them as temporary objects.
       * This guarantees that the accessor objects go out of scope earlier
       * than the mother object, avoid problems with data consistency.
       */
      Accessor (tensor_type              &tensor,
                const TableIndices<rank> &previous_indices);

      /**
       * Copy constructor.
       */
      Accessor (const Accessor &) = default;

    public:

      /**
       * Index operator.
       */
      Accessor<rank,dim,constness,P-1,Number> operator [] (const unsigned int i);

      /**
       * Index operator.
       */
      Accessor<rank,dim,constness,P-1,Number> operator [] (const unsigned int i) const;

    private:
      /**
       * Store the data given to the constructor.
       */
      tensor_type             &tensor;
      const TableIndices<rank> previous_indices;

      // declare some other classes
      // as friends. make sure to
      // work around bugs in some
      // compilers
      template <int,int,typename> friend class dealii::SymmetricTensor;
      template <int,int,bool,int,typename>
      friend class Accessor;
#  ifndef DEAL_II_TEMPL_SPEC_FRIEND_BUG
      friend class ::dealii::SymmetricTensor<rank,dim,Number>;
      friend class Accessor<rank,dim,constness,P+1,Number>;
#  endif
    };



    /**
     * @internal Accessor class for SymmetricTensor. This is the
     * specialization for the last index, which actually allows access to the
     * elements of the table, rather than recursively returning access objects
     * for further subsets. The same holds for this specialization as for the
     * general template; see there for more information.
     *
     * @author Wolfgang Bangerth, 2002, 2005
     */
    template <int rank, int dim, bool constness, typename Number>
    class Accessor<rank,dim,constness,1,Number>
    {
    public:
      /**
       * Import two typedefs from the switch class above.
       */
      typedef typename AccessorTypes<rank,dim,constness,Number>::reference reference;
      typedef typename AccessorTypes<rank,dim,constness,Number>::tensor_type tensor_type;

    private:
      /**
       * Constructor. Take a reference to the tensor object which we will
       * access.
       *
       * The second argument denotes the values of previous indices into the
       * tensor. For example, for a rank-4 tensor, if P=2, then we will
       * already have had two successive element selections (e.g. through
       * <tt>tensor[1][2]</tt>), and the two index values have to be stored
       * somewhere. This class therefore only makes use of the first rank-P
       * elements of this array, but passes it on to the next level with P-1
       * which fills the next entry, and so on.
       *
       * For this particular specialization, i.e. for P==1, all but the last
       * index are already filled.
       *
       * The constructor is made private in order to prevent you having such
       * objects around. The only way to create such objects is via the
       * <tt>Table</tt> class, which only generates them as temporary objects.
       * This guarantees that the accessor objects go out of scope earlier
       * than the mother object, avoid problems with data consistency.
       */
      Accessor (tensor_type              &tensor,
                const TableIndices<rank> &previous_indices);

      /**
       * Default constructor. Not needed, so deleted.
       */
      Accessor () = delete;

      /**
       * Copy constructor.
       */
      Accessor (const Accessor &) = default;

    public:

      /**
       * Index operator.
       */
      reference operator [] (const unsigned int);

      /**
       * Index operator.
       */
      reference operator [] (const unsigned int) const;

    private:
      /**
       * Store the data given to the constructor.
       */
      tensor_type             &tensor;
      const TableIndices<rank> previous_indices;

      // declare some other classes
      // as friends. make sure to
      // work around bugs in some
      // compilers
      template <int,int,typename> friend class dealii::SymmetricTensor;
      template <int,int,bool,int,typename>
      friend class SymmetricTensorAccessors::Accessor;
#  ifndef DEAL_II_TEMPL_SPEC_FRIEND_BUG
      friend class ::dealii::SymmetricTensor<rank,dim,Number>;
      friend class SymmetricTensorAccessors::Accessor<rank,dim,constness,2,Number>;
#  endif
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
 * since fewer elements are stored, and the class automatically makes sure that
 * the tensor represents a symmetric object.
 *
 * For tensors of higher rank, the savings in storage are even higher. For
 * example for the 3x3x3x3 tensors of rank 4, only 36 instead of the full 81
 * entries have to be stored.
 *
 * While the definition of a symmetric rank-2 tensor is obvious, tensors of
 * rank 4 are considered symmetric if they are operators mapping symmetric
 * rank-2 tensors onto symmetric rank-2 tensors. This entails certain symmetry
 * properties on the elements in their 4-dimensional index space, in
 * particular that
 * <tt>C<sub>ijkl</sub>=C<sub>jikl</sub>=C<sub>ijlk</sub></tt>. However, it
 * does not imply the relation <tt>C<sub>ijkl</sub>=C<sub>klij</sub></tt>.
 * Consequently, symmetric tensors of rank 4 as understood here are only
 * tensors that map symmetric tensors onto symmetric tensors, but they do not
 * necessarily induce a symmetric scalar product <tt>a:C:b=b:C:a</tt> or even
 * a positive (semi-)definite form <tt>a:C:a</tt>, where <tt>a,b</tt> are
 * symmetric rank-2 tensors and the colon indicates the common double-index
 * contraction that acts as a product for symmetric tensors.
 *
 * Symmetric tensors are most often used in structural and fluid
 * mechanics, where strains and stresses are usually symmetric
 * tensors, and the stress-strain relationship is given by a symmetric
 * rank-4 tensor.
 *
 * @note Symmetric tensors only exist with even numbers of indices. In
 * other words, the only objects that you can use are
 * <tt>SymmetricTensor<2,dim></tt>, <tt>SymmetricTensor<4,dim></tt>, etc, but
 * <tt>SymmetricTensor<1,dim></tt> and <tt>SymmetricTensor<3,dim></tt> do not
 * exist and their use will most likely lead to compiler errors.
 *
 *
 * <h3>Accessing elements</h3>
 *
 * The elements of a tensor <tt>t</tt> can be accessed using the bracket
 * operator, i.e. for a tensor of rank 4, <tt>t[0][1][0][1]</tt> accesses the
 * element <tt>t<sub>0101</sub></tt>. This access can be used for both reading
 * and writing (if the tensor is non-constant at least). You may also perform
 * other operations on it, although that may lead to confusing situations
 * because several elements of the tensor are stored at the same location. For
 * example, for a rank-2 tensor that is assumed to be zero at the beginning,
 * writing <tt>t[0][1]+=1; t[1][0]+=1;</tt> will lead to the same element
 * being increased by one <em>twice</em>, because even though the accesses use
 * different indices, the elements that are accessed are symmetric and
 * therefore stored at the same location. It may therefore be useful in
 * application programs to restrict operations on individual elements to
 * simple reads or writes.
 *
 * @ingroup geomprimitives
 * @author Wolfgang Bangerth, 2005
 */
template <int rank_, int dim, typename Number>
class SymmetricTensor
{
public:
  static_assert(rank_%2==0, "A SymmetricTensor must have even rank!");

  /**
   * Provide a way to get the dimension of an object without explicit
   * knowledge of it's data type. Implementation is this way instead of
   * providing a function <tt>dimension()</tt> because now it is possible to
   * get the dimension at compile time without the expansion and preevaluation
   * of an inlined function; the compiler may therefore produce more efficient
   * code and you may use this value to declare other data types.
   */
  static const unsigned int dimension = dim;

  /**
   * Publish the rank of this tensor to the outside world.
   */
  static const unsigned int rank = rank_;

  /**
   * An integer denoting the number of independent components that fully
   * describe a symmetric tensor. In $d$ space dimensions, this number equals
   * $\frac 12 (d^2+d)$ for symmetric tensors of rank 2.
   */
  static const unsigned int n_independent_components
    = internal::SymmetricTensorAccessors::StorageType<rank_,dim,Number>::
      n_independent_components;

  /**
   * Default constructor. Creates a tensor with all entries equal to zero.
   */
  SymmetricTensor ();

  /**
   * Constructor. Generate a symmetric tensor from a general one. Assumes that
   * @p t is already symmetric, and in debug mode this is in fact checked.
   * Note that no provision is made to assure that the tensor is symmetric
   * only up to round-off error: if the incoming tensor is not exactly
   * symmetric, then an exception is thrown. If you know that incoming tensor
   * is symmetric only up to round-off, then you may want to call the
   * <tt>symmetrize</tt> function first. If you aren't sure, it is good
   * practice to check before calling <tt>symmetrize</tt>.
   */
  template <typename OtherNumber>
  explicit
  SymmetricTensor (const Tensor<2,dim,OtherNumber> &t);

  /**
   * A constructor that creates a symmetric tensor from an array holding its
   * independent elements. Using this constructor assumes that the caller
   * knows the order in which elements are stored in symmetric tensors; its
   * use is therefore discouraged, but if you think you want to use it anyway
   * you can query the order of elements using the unrolled_index() function.
   *
   * This constructor is currently only implemented for symmetric tensors of
   * rank 2.
   *
   * The size of the array passed is equal to
   * SymmetricTensor<rank_,dim>::n_independent_component; the reason for using
   * the object from the internal namespace is to work around bugs in some
   * older compilers.
   */
  SymmetricTensor (const Number (&array) [n_independent_components]);

  /**
   * Copy constructor from tensors with different underlying scalar type. This
   * obviously requires that the @p OtherNumber type is convertible to @p
   * Number.
   */
  template <typename OtherNumber>
  explicit
  SymmetricTensor (const SymmetricTensor<rank_,dim,OtherNumber> &initializer);

  /**
   * Return a pointer to the first element of the underlying storage.
   */
  Number *
  begin_raw();

  /**
   * Return a const pointer to the first element of the underlying storage.
   */
  const Number *
  begin_raw() const;

  /**
   * Return a pointer to the element past the end of the underlying storage.
   */
  Number *
  end_raw();

  /**
   * Return a const pointer to the element past the end of the underlying
   * storage.
   */
  const Number *
  end_raw() const;

  /**
   * Assignment operator from symmetric tensors with different underlying scalar type.
   * This obviously requires that the @p OtherNumber type is convertible to
   * @p Number.
   */
  template <typename OtherNumber>
  SymmetricTensor &operator = (const SymmetricTensor<rank_,dim,OtherNumber> &rhs);

  /**
   * This operator assigns a scalar to a tensor. To avoid confusion with what
   * exactly it means to assign a scalar value to a tensor, zero is the only
   * value allowed for <tt>d</tt>, allowing the intuitive notation
   * <tt>t=0</tt> to reset all elements of the tensor to zero.
   */
  SymmetricTensor &operator = (const Number &d);

  /**
   * Convert the present symmetric tensor into a full tensor with the same
   * elements, but using the different storage scheme of full tensors.
   */
  operator Tensor<rank_,dim,Number> () const;

  /**
   * Test for equality of two tensors.
   */
  bool operator == (const SymmetricTensor &) const;

  /**
   * Test for inequality of two tensors.
   */
  bool operator != (const SymmetricTensor &) const;

  /**
   * Add another tensor.
   */
  template <typename OtherNumber>
  SymmetricTensor &operator += (const SymmetricTensor<rank_,dim,OtherNumber> &);

  /**
   * Subtract another tensor.
   */
  template <typename OtherNumber>
  SymmetricTensor &operator -= (const SymmetricTensor<rank_,dim,OtherNumber> &);

  /**
   * Scale the tensor by <tt>factor</tt>, i.e. multiply all components by
   * <tt>factor</tt>.
   */
  template <typename OtherNumber>
  SymmetricTensor &operator *= (const OtherNumber &factor);

  /**
   * Scale the tensor by <tt>1/factor</tt>.
   */
  template <typename OtherNumber>
  SymmetricTensor &operator /= (const OtherNumber &factor);

  /**
   * Unary minus operator. Negate all entries of a tensor.
   */
  SymmetricTensor   operator - () const;

  /**
   * Product between the present symmetric tensor and a tensor of rank 2. For
   * example, if the present object is also a rank-2 tensor, then this is the
   * scalar-product double contraction <tt>a<sub>ij</sub>b<sub>ij</sub></tt>
   * over all indices <tt>i,j</tt>. In this case, the return value evaluates
   * to a single scalar. While it is possible to define other scalar product
   * (and associated induced norms), this one seems to be the most appropriate
   * one.
   *
   * If the present object is a rank-4 tensor, then the result is a rank-2
   * tensor, i.e., the operation contracts over the last two indices of the
   * present object and the indices of the argument, and the result is a
   * tensor of rank 2.
   *
   * Note that the multiplication operator for symmetric tensors is defined to
   * be a double contraction over two indices, while it is defined as a single
   * contraction over only one index for regular <tt>Tensor</tt> objects. For
   * symmetric tensors it therefore acts in a way that is commonly denoted by
   * a "colon multiplication" in the mathematical literature.
   *
   * There are global functions <tt>double_contract</tt> that do the same work
   * as this operator, but rather than returning the result as a return value,
   * they write it into the first argument to the function.
   */
  template <typename OtherNumber>
  typename internal::SymmetricTensorAccessors::double_contraction_result<rank_,2,dim,Number,OtherNumber>::type
  operator * (const SymmetricTensor<2,dim,OtherNumber> &s) const;

  /**
   * Contraction over two indices of the present object with the rank-4
   * symmetric tensor given as argument.
   */
  template <typename OtherNumber>
  typename internal::SymmetricTensorAccessors::double_contraction_result<rank_,4,dim,Number,OtherNumber>::type
  operator * (const SymmetricTensor<4,dim,OtherNumber> &s) const;

  /**
   * Return a read-write reference to the indicated element.
   */
  Number &operator() (const TableIndices<rank_> &indices);

  /**
   * Return a @p const reference to the value referred to by the argument.
   */
  const Number &operator() (const TableIndices<rank_> &indices) const;

  /**
   * Access the elements of a row of this symmetric tensor. This function is
   * called for constant tensors.
   */
  internal::SymmetricTensorAccessors::Accessor<rank_,dim,true,rank_-1,Number>
  operator [] (const unsigned int row) const;

  /**
   * Access the elements of a row of this symmetric tensor. This function is
   * called for non-constant tensors.
   */
  internal::SymmetricTensorAccessors::Accessor<rank_,dim,false,rank_-1,Number>
  operator [] (const unsigned int row);

  /**
   * Return a @p const reference to the value referred to by the argument.
   *
   * Exactly the same as operator().
   */
  const Number &
  operator [] (const TableIndices<rank_> &indices) const;

  /**
   * Return a read-write reference to the indicated element.
   *
   * Exactly the same as operator().
   */
  Number &
  operator [] (const TableIndices<rank_> &indices);

  /**
   * Access to an element according to unrolled index. The function
   * <tt>s.access_raw_entry(unrolled_index)</tt> does the same as
   * <tt>s[s.unrolled_to_component_indices(i)]</tt>, but more efficiently.
   */
  const Number &
  access_raw_entry (const unsigned int unrolled_index) const;

  /**
   * Access to an element according to unrolled index. The function
   * <tt>s.access_raw_entry(unrolled_index)</tt> does the same as
   * <tt>s[s.unrolled_to_component_indices(i)]</tt>, but more efficiently.
   */
  Number &
  access_raw_entry (const unsigned int unrolled_index);

  /**
   * Return the Frobenius-norm of a tensor, i.e. the square root of the sum of
   * squares of all entries. This norm is induced by the scalar product
   * defined above for two symmetric tensors. Note that it includes <i>all</i>
   * entries of the tensor, counting symmetry, not only the unique ones (for
   * example, for rank-2 tensors, this norm includes adding up the squares of
   * upper right as well as lower left entries, not just one of them, although
   * they are equal for symmetric tensors).
   */
  typename numbers::NumberTraits<Number>::real_type
  norm () const;

  /**
   * Tensors can be unrolled by simply pasting all elements into one long
   * vector, but for this an order of elements has to be defined. For
   * symmetric tensors, this function returns which index within the range
   * <code>[0,n_independent_components)</code> the given entry in a symmetric
   * tensor has.
   */
  static
  unsigned int
  component_to_unrolled_index (const TableIndices<rank_> &indices);

  /**
   * The opposite of the previous function: given an index $i$ in the unrolled
   * form of the tensor, return what set of indices $(k,l)$ (for rank-2
   * tensors) or $(k,l,m,n)$ (for rank-4 tensors) corresponds to it.
   */
  static
  TableIndices<rank_>
  unrolled_to_component_indices (const unsigned int i);

  /**
   * Reset all values to zero.
   *
   * Note that this is partly inconsistent with the semantics of the @p
   * clear() member functions of the standard library containers and of
   * several other classes within deal.II, which not only reset the values of
   * stored elements to zero, but release all memory and return the object
   * into a virginial state. However, since the size of objects of the present
   * type is determined by its template parameters, resizing is not an option,
   * and indeed the state where all elements have a zero value is the state
   * right after construction of such an object.
   */
  void clear ();

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  static std::size_t memory_consumption ();

  /**
   * Read or write the data of this object to or from a stream for the purpose
   * of serialization
   */
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version);

private:
  /**
   * A structure that describes properties of the base tensor.
   */
  typedef
  internal::SymmetricTensorAccessors::StorageType<rank_,dim,Number>
  base_tensor_descriptor;

  /**
   * Data storage type for a symmetric tensor.
   */
  typedef typename base_tensor_descriptor::base_tensor_type base_tensor_type;

  /**
   * The place where we store the data of the tensor.
   */
  base_tensor_type data;

  /**
   * Make all other symmetric tensors friends.
   */
  template <int, int, typename> friend class SymmetricTensor;

  /**
   * Make a few more functions friends.
   */
  template <int dim2, typename Number2>
  friend Number2 trace (const SymmetricTensor<2,dim2,Number2> &d);

  template <int dim2, typename Number2>
  friend Number2 determinant (const SymmetricTensor<2,dim2,Number2> &t);

  template <int dim2, typename Number2>
  friend SymmetricTensor<2,dim2,Number2>
  deviator (const SymmetricTensor<2,dim2,Number2> &t);

  template <int dim2, typename Number2>
  friend SymmetricTensor<2,dim2,Number2> unit_symmetric_tensor ();

  template <int dim2, typename Number2>
  friend SymmetricTensor<4,dim2,Number2> deviator_tensor ();

  template <int dim2, typename Number2>
  friend SymmetricTensor<4,dim2,Number2> identity_tensor ();


  /**
   * Make a few helper classes friends as well.
   */
  friend struct internal::SymmetricTensorImplementation::Inverse<2,dim,Number>;

  friend struct internal::SymmetricTensorImplementation::Inverse<4,dim,Number>;
};



// ------------------------- inline functions ------------------------

#ifndef DOXYGEN

namespace internal
{
  namespace SymmetricTensorAccessors
  {
    template <int rank_, int dim, bool constness, int P, typename Number>
    Accessor<rank_,dim,constness,P,Number>::
    Accessor (tensor_type              &tensor,
              const TableIndices<rank_> &previous_indices)
      :
      tensor (tensor),
      previous_indices (previous_indices)
    {}



    template <int rank_, int dim, bool constness, int P, typename Number>
    Accessor<rank_,dim,constness,P-1,Number>
    Accessor<rank_,dim,constness,P,Number>::operator[] (const unsigned int i)
    {
      return Accessor<rank_,dim,constness,P-1,Number> (tensor,
                                                       merge (previous_indices, i, rank_-P));
    }



    template <int rank_, int dim, bool constness, int P, typename Number>
    Accessor<rank_,dim,constness,P-1,Number>
    Accessor<rank_,dim,constness,P,Number>::operator[] (const unsigned int i) const
    {
      return Accessor<rank_,dim,constness,P-1,Number> (tensor,
                                                       merge (previous_indices, i, rank_-P));
    }



    template <int rank_, int dim, bool constness, typename Number>
    Accessor<rank_,dim,constness,1,Number>::
    Accessor (tensor_type              &tensor,
              const TableIndices<rank_> &previous_indices)
      :
      tensor (tensor),
      previous_indices (previous_indices)
    {}



    template <int rank_, int dim, bool constness, typename Number>
    typename Accessor<rank_,dim,constness,1,Number>::reference
    Accessor<rank_,dim,constness,1,Number>::operator[] (const unsigned int i)
    {
      return tensor(merge (previous_indices, i, rank_-1));
    }


    template <int rank_, int dim, bool constness, typename Number>
    typename Accessor<rank_,dim,constness,1,Number>::reference
    Accessor<rank_,dim,constness,1,Number>::operator[] (const unsigned int i) const
    {
      return tensor(merge (previous_indices, i, rank_-1));
    }
  }
}



template <int rank_, int dim, typename Number>
inline
SymmetricTensor<rank_,dim,Number>::SymmetricTensor ()
{
  // Some auto-differentiable numbers need explicit
  // zero initialization.
  for (unsigned int i=0; i<base_tensor_type::dimension; ++i)
    data[i] = internal::NumberType<Number>::value(0.0);
}


template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
SymmetricTensor<rank_,dim,Number>::SymmetricTensor (const Tensor<2,dim,OtherNumber> &t)
{
  Assert (rank == 2, ExcNotImplemented());
  switch (dim)
    {
    case 2:
      Assert (t[0][1] == t[1][0], ExcInternalError());

      data[0] = t[0][0];
      data[1] = t[1][1];
      data[2] = t[0][1];

      break;
    case 3:
      Assert (t[0][1] == t[1][0], ExcInternalError());
      Assert (t[0][2] == t[2][0], ExcInternalError());
      Assert (t[1][2] == t[2][1], ExcInternalError());

      data[0] = t[0][0];
      data[1] = t[1][1];
      data[2] = t[2][2];
      data[3] = t[0][1];
      data[4] = t[0][2];
      data[5] = t[1][2];

      break;
    default:
      for (unsigned int d=0; d<dim; ++d)
        for (unsigned int e=0; e<d; ++e)
          Assert(t[d][e] == t[e][d], ExcInternalError());

      for (unsigned int d=0; d<dim; ++d)
        data[d] = t[d][d];

      for (unsigned int d=0, c=0; d<dim; ++d)
        for (unsigned int e=d+1; e<dim; ++e, ++c)
          data[dim+c] = t[d][e];
    }
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
SymmetricTensor<rank_,dim,Number>::
SymmetricTensor (const SymmetricTensor<rank_,dim,OtherNumber> &initializer)
{
  for (unsigned int i=0; i<base_tensor_type::dimension; ++i)
    data[i] = internal::NumberType<typename base_tensor_type::value_type>::value(initializer.data[i]);
}




template <int rank_, int dim, typename Number>
inline
SymmetricTensor<rank_,dim,Number>::SymmetricTensor (const Number (&array) [n_independent_components])
  :
  data (*reinterpret_cast<const typename base_tensor_type::array_type *>(array))
{
  // ensure that the reinterpret_cast above actually works
  Assert (sizeof(typename base_tensor_type::array_type)
          == sizeof(array),
          ExcInternalError());
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
SymmetricTensor<rank_,dim,Number> &
SymmetricTensor<rank_,dim,Number>::operator = (const SymmetricTensor<rank_,dim,OtherNumber> &t)
{
  for (unsigned int i=0; i<base_tensor_type::dimension; ++i)
    data[i] = t.data[i];
  return *this;
}



template <int rank_, int dim, typename Number>
inline
SymmetricTensor<rank_,dim,Number> &
SymmetricTensor<rank_,dim,Number>::operator = (const Number &d)
{
  Assert (d==internal::NumberType<Number>::value(0.0), ExcMessage ("Only assignment with zero is allowed"));
  (void) d;

  data = internal::NumberType<Number>::value(0.0);

  return *this;
}


namespace internal
{
  namespace SymmetricTensorImplementation
  {
    template <int dim, typename Number>
    inline DEAL_II_ALWAYS_INLINE
    dealii::Tensor<2,dim,Number>
    convert_to_tensor (const dealii::SymmetricTensor<2,dim,Number> &s)
    {
      dealii::Tensor<2,dim,Number> t;

      // diagonal entries are stored first
      for (unsigned int d=0; d<dim; ++d)
        t[d][d] = s.access_raw_entry(d);

      // off-diagonal entries come next, row by row
      for (unsigned int d=0, c=0; d<dim; ++d)
        for (unsigned int e=d+1; e<dim; ++e, ++c)
          {
            t[d][e] = s.access_raw_entry(dim+c);
            t[e][d] = s.access_raw_entry(dim+c);
          }
      return t;
    }


    template <int dim, typename Number>
    dealii::Tensor<4,dim,Number>
    convert_to_tensor (const dealii::SymmetricTensor<4,dim,Number> &st)
    {
      // utilize the symmetry properties of SymmetricTensor<4,dim>
      // discussed in the class documentation to avoid accessing all
      // independent elements of the input tensor more than once
      dealii::Tensor<4,dim,Number> t;

      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=i; j<dim; ++j)
          for (unsigned int k=0; k<dim; ++k)
            for (unsigned int l=k; l<dim; ++l)
              t[TableIndices<4>(i,j,k,l)]
                = t[TableIndices<4>(i,j,l,k)]
                  = t[TableIndices<4>(j,i,k,l)]
                    = t[TableIndices<4>(j,i,l,k)]
                      = st[TableIndices<4>(i,j,k,l)];

      return t;
    }


    template <typename Number>
    struct Inverse<2,1,Number>
    {
      static inline dealii::SymmetricTensor<2,1,Number>
      value (const dealii::SymmetricTensor<2,1,Number> &t)
      {
        dealii::SymmetricTensor<2,1,Number> tmp;

        tmp[0][0] = 1.0/t[0][0];

        return tmp;
      }
    };


    template <typename Number>
    struct Inverse<2,2,Number>
    {
      static inline dealii::SymmetricTensor<2,2,Number>
      value (const dealii::SymmetricTensor<2,2,Number> &t)
      {
        dealii::SymmetricTensor<2,2,Number> tmp;

        // Sympy result: ([
        // [ t11/(t00*t11 - t01**2), -t01/(t00*t11 - t01**2)],
        // [-t01/(t00*t11 - t01**2),  t00/(t00*t11 - t01**2)]  ])
        const TableIndices<2> idx_00 (0,0);
        const TableIndices<2> idx_01 (0,1);
        const TableIndices<2> idx_11 (1,1);
        const Number inv_det_t
          = 1.0/(t[idx_00]*t[idx_11]
                 - t[idx_01]*t[idx_01]);
        tmp[idx_00] = t[idx_11];
        tmp[idx_01] = -t[idx_01];
        tmp[idx_11] = t[idx_00];
        tmp *= inv_det_t;

        return tmp;
      }
    };


    template <typename Number>
    struct Inverse<2,3,Number>
    {
      static dealii::SymmetricTensor<2,3,Number>
      value (const dealii::SymmetricTensor<2,3,Number> &t)
      {
        dealii::SymmetricTensor<2,3,Number> tmp;

        // Sympy result: ([
        // [  (t11*t22 - t12**2)/(t00*t11*t22 - t00*t12**2 - t01**2*t22 + 2*t01*t02*t12 - t02**2*t11),
        //    (-t01*t22 + t02*t12)/(t00*t11*t22 - t00*t12**2 - t01**2*t22 + 2*t01*t02*t12 - t02**2*t11),
        //    (t01*t12 - t02*t11)/(t00*t11*t22 - t00*t12**2 - t01**2*t22 + 2*t01*t02*t12 - t02**2*t11)],
        // [  (-t01*t22 + t02*t12)/(t00*t11*t22 - t00*t12**2 - t01**2*t22 + 2*t01*t02*t12 - t02**2*t11),
        //    (t00*t22 - t02**2)/(t00*t11*t22 - t00*t12**2 - t01**2*t22 + 2*t01*t02*t12 - t02**2*t11),
        //    (t00*t12 - t01*t02)/(-t00*t11*t22 + t00*t12**2 + t01**2*t22 - 2*t01*t02*t12 + t02**2*t11)],
        // [  (t01*t12 - t02*t11)/(t00*t11*t22 - t00*t12**2 - t01**2*t22 + 2*t01*t02*t12 - t02**2*t11),
        //    (t00*t12 - t01*t02)/(-t00*t11*t22 + t00*t12**2 + t01**2*t22 - 2*t01*t02*t12 + t02**2*t11),
        //    (-t00*t11 + t01**2)/(-t00*t11*t22 + t00*t12**2 + t01**2*t22 - 2*t01*t02*t12 + t02**2*t11)]  ])
        const TableIndices<2> idx_00 (0,0);
        const TableIndices<2> idx_01 (0,1);
        const TableIndices<2> idx_02 (0,2);
        const TableIndices<2> idx_11 (1,1);
        const TableIndices<2> idx_12 (1,2);
        const TableIndices<2> idx_22 (2,2);
        const Number inv_det_t
          = 1.0/(t[idx_00]*t[idx_11]*t[idx_22]
                 - t[idx_00]*t[idx_12]*t[idx_12]
                 - t[idx_01]*t[idx_01]*t[idx_22]
                 + 2.0*t[idx_01]*t[idx_02]*t[idx_12]
                 - t[idx_02]*t[idx_02]*t[idx_11]);
        tmp[idx_00] = t[idx_11]*t[idx_22] - t[idx_12]*t[idx_12];
        tmp[idx_01] = -t[idx_01]*t[idx_22] + t[idx_02]*t[idx_12];
        tmp[idx_02] = t[idx_01]*t[idx_12] - t[idx_02]*t[idx_11];
        tmp[idx_11] = t[idx_00]*t[idx_22] - t[idx_02]*t[idx_02];
        tmp[idx_12] = -t[idx_00]*t[idx_12] + t[idx_01]*t[idx_02];
        tmp[idx_22] = t[idx_00]*t[idx_11] - t[idx_01]*t[idx_01];
        tmp *= inv_det_t;

        return tmp;
      }
    };


    template <typename Number>
    struct Inverse<4,1,Number>
    {
      static inline dealii::SymmetricTensor<4,1,Number>
      value (const dealii::SymmetricTensor<4,1,Number> &t)
      {
        dealii::SymmetricTensor<4,1,Number> tmp;
        tmp.data[0][0] = 1.0/t.data[0][0];
        return tmp;
      }
    };


    template <typename Number>
    struct Inverse<4,2,Number>
    {
      static inline dealii::SymmetricTensor<4,2,Number>
      value (const dealii::SymmetricTensor<4,2,Number> &t)
      {
        dealii::SymmetricTensor<4,2,Number> tmp;

        // Inverting this tensor is a little more complicated than necessary,
        // since we store the data of 't' as a 3x3 matrix t.data, but the
        // product between a rank-4 and a rank-2 tensor is really not the
        // product between this matrix and the 3-vector of a rhs, but rather
        //
        // B.vec = t.data * mult * A.vec
        //
        // where mult is a 3x3 matrix with entries [[1,0,0],[0,1,0],[0,0,2]] to
        // capture the fact that we need to add up both the c_ij12*a_12 and the
        // c_ij21*a_21 terms.
        //
        // In addition, in this scheme, the identity tensor has the matrix
        // representation mult^-1.
        //
        // The inverse of 't' therefore has the matrix representation
        //
        // inv.data = mult^-1 * t.data^-1 * mult^-1
        //
        // in order to compute it, let's first compute the inverse of t.data and
        // put it into tmp.data; at the end of the function we then scale the
        // last row and column of the inverse by 1/2, corresponding to the left
        // and right multiplication with mult^-1.
        const Number t4 = t.data[0][0]*t.data[1][1],
                     t6 = t.data[0][0]*t.data[1][2],
                     t8 = t.data[0][1]*t.data[1][0],
                     t00 = t.data[0][2]*t.data[1][0],
                     t01 = t.data[0][1]*t.data[2][0],
                     t04 = t.data[0][2]*t.data[2][0],
                     t07 = 1.0/(t4*t.data[2][2]-t6*t.data[2][1]-
                                t8*t.data[2][2]+t00*t.data[2][1]+
                                t01*t.data[1][2]-t04*t.data[1][1]);
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
    };


    template <typename Number>
    struct Inverse<4,3,Number>
    {
      static dealii::SymmetricTensor<4,3,Number>
      value (const dealii::SymmetricTensor<4,3,Number> &t)
      {
        dealii::SymmetricTensor<4,3,Number> tmp = t;

        // This function follows the exact same scheme as the 2d case, except
        // that hardcoding the inverse of a 6x6 matrix is pretty wasteful.
        // Instead, we use the Gauss-Jordan algorithm implemented for
        // FullMatrix. For historical reasons the following code is copied from
        // there, with the tangential benefit that we do not need to copy the
        // tensor entries to and from the FullMatrix.
        const unsigned int N = 6;

        // First get an estimate of the size of the elements of this matrix,
        // for later checks whether the pivot element is large enough, or
        // whether we have to fear that the matrix is not regular.
        Number diagonal_sum = internal::NumberType<Number>::value(0.0);
        for (unsigned int i=0; i<N; ++i)
          diagonal_sum += std::fabs(tmp.data[i][i]);
        const Number typical_diagonal_element = diagonal_sum/static_cast<double>(N);
        (void)typical_diagonal_element;

        unsigned int p[N];
        for (unsigned int i=0; i<N; ++i)
          p[i] = i;

        for (unsigned int j=0; j<N; ++j)
          {
            // Pivot search: search that part of the line on and right of the
            // diagonal for the largest element.
            Number       max = std::fabs(tmp.data[j][j]);
            unsigned int r   = j;
            for (unsigned int i=j+1; i<N; ++i)
              if (std::fabs(tmp.data[i][j]) > max)
                {
                  max = std::fabs(tmp.data[i][j]);
                  r = i;
                }

            // Check whether the pivot is too small
            Assert(max > 1.e-16*typical_diagonal_element,
                   ExcMessage("This tensor seems to be noninvertible"));

            // Row interchange
            if (r>j)
              {
                for (unsigned int k=0; k<N; ++k)
                  std::swap (tmp.data[j][k], tmp.data[r][k]);

                std::swap (p[j], p[r]);
              }

            // Transformation
            const Number hr = 1./tmp.data[j][j];
            tmp.data[j][j] = hr;
            for (unsigned int k=0; k<N; ++k)
              {
                if (k==j) continue;
                for (unsigned int i=0; i<N; ++i)
                  {
                    if (i==j) continue;
                    tmp.data[i][k] -= tmp.data[i][j]*tmp.data[j][k]*hr;
                  }
              }
            for (unsigned int i=0; i<N; ++i)
              {
                tmp.data[i][j] *= hr;
                tmp.data[j][i] *= -hr;
              }
            tmp.data[j][j] = hr;
          }

        // Column interchange
        Number hv[N];
        for (unsigned int i=0; i<N; ++i)
          {
            for (unsigned int k=0; k<N; ++k)
              hv[p[k]] = tmp.data[i][k];
            for (unsigned int k=0; k<N; ++k)
              tmp.data[i][k] = hv[k];
          }

        // Scale rows and columns. The mult matrix
        // here is diag[1, 1, 1, 1/2, 1/2, 1/2].
        for (unsigned int i=3; i<6; ++i)
          for (unsigned int j=0; j<3; ++j)
            tmp.data[i][j] /= 2;

        for (unsigned int i=0; i<3; ++i)
          for (unsigned int j=3; j<6; ++j)
            tmp.data[i][j] /= 2;

        for (unsigned int i=3; i<6; ++i)
          for (unsigned int j=3; j<6; ++j)
            tmp.data[i][j] /= 4;

        return tmp;
      }
    };

  }
}



template <int rank_, int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
SymmetricTensor<rank_,dim,Number>::
operator Tensor<rank_,dim,Number> () const
{
  return internal::SymmetricTensorImplementation::convert_to_tensor (*this);
}



template <int rank_, int dim, typename Number>
inline
bool
SymmetricTensor<rank_,dim,Number>::operator ==
(const SymmetricTensor<rank_,dim,Number> &t) const
{
  return data == t.data;
}



template <int rank_, int dim, typename Number>
inline
bool
SymmetricTensor<rank_,dim,Number>::operator !=
(const SymmetricTensor<rank_,dim,Number> &t) const
{
  return data != t.data;
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
SymmetricTensor<rank_,dim,Number> &
SymmetricTensor<rank_,dim,Number>::operator +=
(const SymmetricTensor<rank_,dim,OtherNumber> &t)
{
  data += t.data;
  return *this;
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
SymmetricTensor<rank_,dim,Number> &
SymmetricTensor<rank_,dim,Number>::operator -=
(const SymmetricTensor<rank_,dim,OtherNumber> &t)
{
  data -= t.data;
  return *this;
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
SymmetricTensor<rank_,dim,Number> &
SymmetricTensor<rank_,dim,Number>::operator *= (const OtherNumber &d)
{
  data *= d;
  return *this;
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
SymmetricTensor<rank_,dim,Number> &
SymmetricTensor<rank_,dim,Number>::operator /= (const OtherNumber &d)
{
  data /= d;
  return *this;
}



template <int rank_, int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
SymmetricTensor<rank_,dim,Number>
SymmetricTensor<rank_,dim,Number>::operator - () const
{
  SymmetricTensor tmp = *this;
  tmp.data = -tmp.data;
  return tmp;
}



template <int rank_, int dim, typename Number>
inline DEAL_II_ALWAYS_INLINE
void
SymmetricTensor<rank_,dim,Number>::clear ()
{
  data.clear ();
}



template <int rank_, int dim, typename Number>
inline
std::size_t
SymmetricTensor<rank_,dim,Number>::memory_consumption ()
{
  // all memory consists of statically allocated memory of the current
  // object, no pointers
  return sizeof(SymmetricTensor<rank_,dim,Number>);
}



namespace internal
{

  template <int dim, typename Number, typename OtherNumber = Number>
  inline DEAL_II_ALWAYS_INLINE
  typename SymmetricTensorAccessors::double_contraction_result<2,2,dim,Number,OtherNumber>::type
  perform_double_contraction (const typename SymmetricTensorAccessors::StorageType<2,dim,Number>::base_tensor_type      &data,
                              const typename SymmetricTensorAccessors::StorageType<2,dim,OtherNumber>::base_tensor_type &sdata)
  {
    typedef typename SymmetricTensorAccessors::double_contraction_result<2,2,dim,Number,OtherNumber>::type result_type;

    switch (dim)
      {
      case 1:
        return data[0] * sdata[0];
      default:
        // Start with the non-diagonal part to avoid some multiplications by
        // 2.

        result_type sum = data[dim] * sdata[dim];
        for (unsigned int d=dim+1; d<(dim*(dim+1)/2); ++d)
          sum += data[d] * sdata[d];
        sum += sum; // sum = sum * 2.;

        // Now add the contributions from the diagonal
        for (unsigned int d=0; d<dim; ++d)
          sum += data[d] * sdata[d];
        return sum;
      }
  }



  template <int dim, typename Number, typename OtherNumber = Number>
  inline
  typename SymmetricTensorAccessors::double_contraction_result<4,2,dim,Number,OtherNumber>::type
  perform_double_contraction (const typename SymmetricTensorAccessors::StorageType<4,dim,Number>::base_tensor_type      &data,
                              const typename SymmetricTensorAccessors::StorageType<2,dim,OtherNumber>::base_tensor_type &sdata)
  {
    typedef typename SymmetricTensorAccessors::double_contraction_result<4,2,dim,Number,OtherNumber>::type       result_type;
    typedef typename SymmetricTensorAccessors::double_contraction_result<4,2,dim,Number,OtherNumber>::value_type value_type;

    const unsigned int data_dim =
      SymmetricTensorAccessors::StorageType<2,dim,value_type>::n_independent_components;
    value_type tmp [data_dim];
    for (unsigned int i=0; i<data_dim; ++i)
      tmp[i] = perform_double_contraction<dim,Number,OtherNumber>(data[i], sdata);
    return result_type(tmp);
  }



  template <int dim, typename Number, typename OtherNumber = Number>
  inline
  typename SymmetricTensorAccessors::StorageType<2,dim,
           typename SymmetricTensorAccessors::double_contraction_result<2,4,dim,Number,OtherNumber>::value_type
           >::base_tensor_type
           perform_double_contraction (const typename SymmetricTensorAccessors::StorageType<2,dim,Number>::base_tensor_type      &data,
                                       const typename SymmetricTensorAccessors::StorageType<4,dim,OtherNumber>::base_tensor_type &sdata)
  {
    typedef typename SymmetricTensorAccessors::double_contraction_result<2,4,dim,Number,OtherNumber>::value_type value_type;
    typedef typename SymmetricTensorAccessors::StorageType<2,dim,value_type>::base_tensor_type base_tensor_type;

    base_tensor_type tmp;
    for (unsigned int i=0; i<tmp.dimension; ++i)
      {
        // Start with the non-diagonal part
        value_type sum = data[dim] * sdata[dim][i];
        for (unsigned int d=dim+1; d<(dim*(dim+1)/2); ++d)
          sum += data[d] * sdata[d][i];
        sum += sum; // sum = sum * 2.;

        // Now add the contributions from the diagonal
        for (unsigned int d=0; d<dim; ++d)
          sum += data[d] * sdata[d][i];
        tmp[i] = sum;
      }
    return tmp;
  }



  template <int dim, typename Number, typename OtherNumber = Number>
  inline
  typename SymmetricTensorAccessors::StorageType<4,dim,
           typename SymmetricTensorAccessors::double_contraction_result<4,4,dim,Number,OtherNumber>::value_type
           >::base_tensor_type
           perform_double_contraction (const typename SymmetricTensorAccessors::StorageType<4,dim,Number>::base_tensor_type      &data,
                                       const typename SymmetricTensorAccessors::StorageType<4,dim,OtherNumber>::base_tensor_type &sdata)
  {
    typedef typename SymmetricTensorAccessors::double_contraction_result<4,4,dim,Number,OtherNumber>::value_type value_type;
    typedef typename SymmetricTensorAccessors::StorageType<4,dim,value_type>::base_tensor_type base_tensor_type;

    const unsigned int data_dim =
      SymmetricTensorAccessors::StorageType<2,dim,value_type>::n_independent_components;
    base_tensor_type tmp;
    for (unsigned int i=0; i<data_dim; ++i)
      for (unsigned int j=0; j<data_dim; ++j)
        {
          // Start with the non-diagonal part
          for (unsigned int d=dim; d<(dim*(dim+1)/2); ++d)
            tmp[i][j] += data[i][d] * sdata[d][j];
          tmp[i][j] += tmp[i][j]; // tmp[i][j] = tmp[i][j] * 2;

          // Now add the contributions from the diagonal
          for (unsigned int d=0; d<dim; ++d)
            tmp[i][j] += data[i][d] * sdata[d][j];
        }
    return tmp;
  }

} // end of namespace internal



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
typename internal::SymmetricTensorAccessors::double_contraction_result<rank_,2,dim,Number,OtherNumber>::type
SymmetricTensor<rank_,dim,Number>::operator * (const SymmetricTensor<2,dim,OtherNumber> &s) const
{
  // need to have two different function calls
  // because a scalar and rank-2 tensor are not
  // the same data type (see internal function
  // above)
  return internal::perform_double_contraction<dim,Number,OtherNumber> (data, s.data);
}



template <int rank_, int dim, typename Number>
template <typename OtherNumber>
inline
typename internal::SymmetricTensorAccessors::double_contraction_result<rank_,4,dim,Number,OtherNumber>::type
SymmetricTensor<rank_,dim,Number>::operator * (const SymmetricTensor<4,dim,OtherNumber> &s) const
{
  typename internal::SymmetricTensorAccessors::
  double_contraction_result<rank_,4,dim,Number,OtherNumber>::type tmp;
  tmp.data = internal::perform_double_contraction<dim,Number,OtherNumber> (data,s.data);
  return tmp;
}



// internal namespace to switch between the
// access of different tensors. There used to
// be explicit instantiations before for
// different ranks and dimensions, but since
// we now allow for templates on the data
// type, and since we cannot partially
// specialize the implementation, this got
// into a separate namespace
namespace internal
{
  template <int dim, typename Number>
  inline
  Number &
  symmetric_tensor_access (const TableIndices<2> &indices,
                           typename SymmetricTensorAccessors::StorageType<2,dim,Number>::base_tensor_type &data)
  {
    // 1d is very simple and done first
    if (dim == 1)
      return data[0];

    // first treat the main diagonal elements, which are stored consecutively
    // at the beginning
    if (indices[0] == indices[1])
      return data[indices[0]];

    // the rest is messier and requires a few switches.
    switch (dim)
      {
      case 2:
        // at least for the 2x2 case it is reasonably simple
        Assert (((indices[0]==1) && (indices[1]==0)) ||
                ((indices[0]==0) && (indices[1]==1)),
                ExcInternalError());
        return data[2];

      default:
        // to do the rest, sort our indices before comparing
      {
        TableIndices<2> sorted_indices (indices);
        sorted_indices.sort ();

        for (unsigned int d=0, c=0; d<dim; ++d)
          for (unsigned int e=d+1; e<dim; ++e, ++c)
            if ((sorted_indices[0]==d) && (sorted_indices[1]==e))
              return data[dim+c];
        Assert (false, ExcInternalError());
      }
      }

    static Number dummy_but_referenceable = Number();
    return dummy_but_referenceable;
  }



  template <int dim, typename Number>
  inline
  const Number &
  symmetric_tensor_access (const TableIndices<2> &indices,
                           const typename SymmetricTensorAccessors::StorageType<2,dim,Number>::base_tensor_type &data)
  {
    // 1d is very simple and done first
    if (dim == 1)
      return data[0];

    // first treat the main diagonal elements, which are stored consecutively
    // at the beginning
    if (indices[0] == indices[1])
      return data[indices[0]];

    // the rest is messier and requires a few switches.
    switch (dim)
      {
      case 2:
        // at least for the 2x2 case it is reasonably simple
        Assert (((indices[0]==1) && (indices[1]==0)) ||
                ((indices[0]==0) && (indices[1]==1)),
                ExcInternalError());
        return data[2];

      default:
        // to do the rest, sort our indices before comparing
      {
        TableIndices<2> sorted_indices (indices);
        sorted_indices.sort ();

        for (unsigned int d=0, c=0; d<dim; ++d)
          for (unsigned int e=d+1; e<dim; ++e, ++c)
            if ((sorted_indices[0]==d) && (sorted_indices[1]==e))
              return data[dim+c];
        Assert (false, ExcInternalError());
      }
      }

    static Number dummy_but_referenceable = Number();
    return dummy_but_referenceable;
  }



  template <int dim, typename Number>
  inline
  Number &
  symmetric_tensor_access (const TableIndices<4> &indices,
                           typename SymmetricTensorAccessors::StorageType<4,dim,Number>::base_tensor_type &data)
  {
    switch (dim)
      {
      case 1:
        return data[0][0];

      case 2:
        // each entry of the tensor can be
        // thought of as an entry in a
        // matrix that maps the rolled-out
        // rank-2 tensors into rolled-out
        // rank-2 tensors. this is the
        // format in which we store rank-4
        // tensors. determine which
        // position the present entry is
        // stored in
      {
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

      case 3:
        // each entry of the tensor can be
        // thought of as an entry in a
        // matrix that maps the rolled-out
        // rank-2 tensors into rolled-out
        // rank-2 tensors. this is the
        // format in which we store rank-4
        // tensors. determine which
        // position the present entry is
        // stored in
      {
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

      default:
        Assert (false, ExcNotImplemented());
      }

    static Number dummy;
    return dummy;
  }


  template <int dim, typename Number>
  inline
  const Number &
  symmetric_tensor_access (const TableIndices<4> &indices,
                           const typename SymmetricTensorAccessors::StorageType<4,dim,Number>::base_tensor_type &data)
  {
    switch (dim)
      {
      case 1:
        return data[0][0];

      case 2:
        // each entry of the tensor can be
        // thought of as an entry in a
        // matrix that maps the rolled-out
        // rank-2 tensors into rolled-out
        // rank-2 tensors. this is the
        // format in which we store rank-4
        // tensors. determine which
        // position the present entry is
        // stored in
      {
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

      case 3:
        // each entry of the tensor can be
        // thought of as an entry in a
        // matrix that maps the rolled-out
        // rank-2 tensors into rolled-out
        // rank-2 tensors. this is the
        // format in which we store rank-4
        // tensors. determine which
        // position the present entry is
        // stored in
      {
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

      default:
        Assert (false, ExcNotImplemented());
      }

    static Number dummy;
    return dummy;
  }

} // end of namespace internal



template <int rank_, int dim, typename Number>
inline
Number &
SymmetricTensor<rank_,dim,Number>::operator () (const TableIndices<rank_> &indices)
{
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dimension, ExcIndexRange (indices[r], 0, dimension));
  return internal::symmetric_tensor_access<dim,Number> (indices, data);
}



template <int rank_, int dim, typename Number>
inline
const Number &
SymmetricTensor<rank_,dim,Number>::operator ()
(const TableIndices<rank_> &indices) const
{
  for (unsigned int r=0; r<rank; ++r)
    Assert (indices[r] < dimension, ExcIndexRange (indices[r], 0, dimension));
  return internal::symmetric_tensor_access<dim,Number> (indices, data);
}



namespace internal
{
  namespace SymmetricTensorImplementation
  {
    template <int rank_>
    TableIndices<rank_>
    get_partially_filled_indices (const unsigned int row,
                                  const std::integral_constant<int, 2> &)
    {
      return TableIndices<rank_> (row,
                                  numbers::invalid_unsigned_int);

    }


    template <int rank_>
    TableIndices<rank_>
    get_partially_filled_indices (const unsigned int row,
                                  const std::integral_constant<int, 4> &)
    {
      return TableIndices<rank_> (row,
                                  numbers::invalid_unsigned_int,
                                  numbers::invalid_unsigned_int,
                                  numbers::invalid_unsigned_int);

    }
  }
}


template <int rank_, int dim, typename Number>
internal::SymmetricTensorAccessors::Accessor<rank_,dim,true,rank_-1,Number>
SymmetricTensor<rank_,dim,Number>::operator [] (const unsigned int row) const
{
  return
    internal::SymmetricTensorAccessors::
    Accessor<rank_,dim,true,rank_-1,Number> (*this,
                                             internal::SymmetricTensorImplementation::get_partially_filled_indices<rank_> (row,
                                                 std::integral_constant<int, rank_>()));
}



template <int rank_, int dim, typename Number>
internal::SymmetricTensorAccessors::Accessor<rank_,dim,false,rank_-1,Number>
SymmetricTensor<rank_,dim,Number>::operator [] (const unsigned int row)
{
  return
    internal::SymmetricTensorAccessors::
    Accessor<rank_,dim,false,rank_-1,Number> (*this,
                                              internal::SymmetricTensorImplementation::get_partially_filled_indices<rank_> (row,
                                                  std::integral_constant<int, rank_>()));
}



template <int rank_, int dim, typename Number>
inline
const Number &
SymmetricTensor<rank_,dim,Number>::operator [] (const TableIndices<rank_> &indices) const
{
  return operator()(indices);
}



template <int rank_, int dim, typename Number>
inline
Number &
SymmetricTensor<rank_,dim,Number>::operator [] (const TableIndices<rank_> &indices)
{
  return operator()(indices);
}



template <int rank_, int dim, typename Number>
inline
Number *
SymmetricTensor<rank_,dim,Number>::begin_raw()
{
  return std::addressof(this->access_raw_entry(0));
}



template <int rank_, int dim, typename Number>
inline
const Number *
SymmetricTensor<rank_,dim,Number>::begin_raw() const
{
  return std::addressof(this->access_raw_entry(0));
}



template <int rank_, int dim, typename Number>
inline
Number *
SymmetricTensor<rank_,dim,Number>::end_raw()
{
  return begin_raw()+n_independent_components;
}



template <int rank_, int dim, typename Number>
inline
const Number *
SymmetricTensor<rank_,dim,Number>::end_raw() const
{
  return begin_raw()+n_independent_components;
}



namespace internal
{
  namespace SymmetricTensorImplementation
  {
    template <int dim, typename Number>
    unsigned int
    entry_to_indices (const dealii::SymmetricTensor<2,dim,Number> &,
                      const unsigned int index)
    {
      return index;
    }


    template <int dim, typename Number>
    dealii::TableIndices<2>
    entry_to_indices (const dealii::SymmetricTensor<4,dim,Number> &,
                      const unsigned int index)
    {
      return
        internal::SymmetricTensorAccessors::StorageType<4,dim,Number>::base_tensor_type::
        unrolled_to_component_indices(index);
    }

  }
}



template <int rank_, int dim, typename Number>
inline
const Number &
SymmetricTensor<rank_,dim,Number>::access_raw_entry (const unsigned int index) const
{
  AssertIndexRange (index, n_independent_components);
  return data[internal::SymmetricTensorImplementation::entry_to_indices(*this, index)];
}



template <int rank_, int dim, typename Number>
inline
Number &
SymmetricTensor<rank_,dim,Number>::access_raw_entry (const unsigned int index)
{
  AssertIndexRange (index, n_independent_components);
  return data[internal::SymmetricTensorImplementation::entry_to_indices(*this, index)];
}



namespace internal
{
  template <int dim, typename Number>
  inline
  typename numbers::NumberTraits<Number>::real_type
  compute_norm (const typename SymmetricTensorAccessors::StorageType<2,dim,Number>::base_tensor_type &data)
  {
    switch (dim)
      {
      case 1:
        return numbers::NumberTraits<Number>::abs(data[0]);

      case 2:
        return std::sqrt(numbers::NumberTraits<Number>::abs_square(data[0]) +
                         numbers::NumberTraits<Number>::abs_square(data[1]) +
                         2. * numbers::NumberTraits<Number>::abs_square(data[2]));

      case 3:
        return std::sqrt(numbers::NumberTraits<Number>::abs_square(data[0]) +
                         numbers::NumberTraits<Number>::abs_square(data[1]) +
                         numbers::NumberTraits<Number>::abs_square(data[2]) +
                         2. * numbers::NumberTraits<Number>::abs_square(data[3]) +
                         2. * numbers::NumberTraits<Number>::abs_square(data[4]) +
                         2. * numbers::NumberTraits<Number>::abs_square(data[5]));

      default:
      {
        typename numbers::NumberTraits<Number>::real_type return_value
          = typename numbers::NumberTraits<Number>::real_type();

        for (unsigned int d=0; d<dim; ++d)
          return_value += numbers::NumberTraits<Number>::abs_square(data[d]);
        for (unsigned int d=dim; d<(dim*dim+dim)/2; ++d)
          return_value += 2. * numbers::NumberTraits<Number>::abs_square(data[d]);

        return std::sqrt(return_value);
      }
      }
  }



  template <int dim, typename Number>
  inline
  typename numbers::NumberTraits<Number>::real_type
  compute_norm (const typename SymmetricTensorAccessors::StorageType<4,dim,Number>::base_tensor_type &data)
  {
    switch (dim)
      {
      case 1:
        return numbers::NumberTraits<Number>::abs (data[0][0]);

      default:
      {
        typename numbers::NumberTraits<Number>::real_type return_value
          = typename numbers::NumberTraits<Number>::real_type();

        const unsigned int n_independent_components = data.dimension;

        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=0; j<dim; ++j)
            return_value += numbers::NumberTraits<Number>::abs_square(data[i][j]);
        for (unsigned int i=0; i<dim; ++i)
          for (unsigned int j=dim; j<n_independent_components; ++j)
            return_value += 2. * numbers::NumberTraits<Number>::abs_square(data[i][j]);
        for (unsigned int i=dim; i<n_independent_components; ++i)
          for (unsigned int j=0; j<dim; ++j)
            return_value += 2. * numbers::NumberTraits<Number>::abs_square(data[i][j]);
        for (unsigned int i=dim; i<n_independent_components; ++i)
          for (unsigned int j=dim; j<n_independent_components; ++j)
            return_value += 4. * numbers::NumberTraits<Number>::abs_square(data[i][j]);

        return std::sqrt(return_value);
      }
      }
  }

} // end of namespace internal



template <int rank_, int dim, typename Number>
inline
typename numbers::NumberTraits<Number>::real_type
SymmetricTensor<rank_,dim,Number>::norm () const
{
  return internal::compute_norm<dim,Number> (data);
}



namespace internal
{
  namespace SymmetricTensorImplementation
  {
    namespace
    {
      // a function to do the unrolling from a set of indices to a
      // scalar index into the array in which we store the elements of
      // a symmetric tensor
      //
      // this function is for rank-2 tensors
      template <int dim>
      inline
      unsigned int
      component_to_unrolled_index
      (const TableIndices<2> &indices)
      {
        Assert (indices[0] < dim, ExcIndexRange(indices[0], 0, dim));
        Assert (indices[1] < dim, ExcIndexRange(indices[1], 0, dim));

        switch (dim)
          {
          case 1:
          {
            return 0;
          }

          case 2:
          {
            static const unsigned int table[2][2] = {{0, 2},
              {2, 1}
            };
            return table[indices[0]][indices[1]];
          }

          case 3:
          {
            static const unsigned int table[3][3] = {{0, 3, 4},
              {3, 1, 5},
              {4, 5, 2}
            };
            return table[indices[0]][indices[1]];
          }

          case 4:
          {
            static const unsigned int table[4][4] = {{0, 4, 5, 6},
              {4, 1, 7, 8},
              {5, 7, 2, 9},
              {6, 8, 9, 3}
            };
            return table[indices[0]][indices[1]];
          }

          default:
            // for the remainder, manually figure out the numbering
          {
            if (indices[0] == indices[1])
              return indices[0];

            TableIndices<2> sorted_indices (indices);
            sorted_indices.sort ();

            for (unsigned int d=0, c=0; d<dim; ++d)
              for (unsigned int e=d+1; e<dim; ++e, ++c)
                if ((sorted_indices[0]==d) && (sorted_indices[1]==e))
                  return dim+c;

            // should never get here:
            Assert(false, ExcInternalError());
            return 0;
          }
          }
      }

      // a function to do the unrolling from a set of indices to a
      // scalar index into the array in which we store the elements of
      // a symmetric tensor
      //
      // this function is for tensors of ranks not already handled
      // above
      template <int dim, int rank_>
      inline
      unsigned int
      component_to_unrolled_index
      (const TableIndices<rank_> &indices)
      {
        (void)indices;
        Assert (false, ExcNotImplemented());
        return numbers::invalid_unsigned_int;
      }
    }
  }
}


template <int rank_, int dim, typename Number>
inline
unsigned int
SymmetricTensor<rank_,dim,Number>::component_to_unrolled_index
(const TableIndices<rank_> &indices)
{
  return internal::SymmetricTensorImplementation::component_to_unrolled_index<dim> (indices);
}



namespace internal
{
  namespace SymmetricTensorImplementation
  {
    namespace
    {
      // a function to do the inverse of the unrolling from a set of
      // indices to a scalar index into the array in which we store
      // the elements of a symmetric tensor. in other words, it goes
      // from the scalar index into the array to a set of indices of
      // the tensor
      //
      // this function is for rank-2 tensors
      template <int dim>
      inline
      TableIndices<2>
      unrolled_to_component_indices
      (const unsigned int i,
       const std::integral_constant<int, 2> &)
      {
        Assert ((i < dealii::SymmetricTensor<2,dim,double>::n_independent_components),
                ExcIndexRange(i, 0, dealii::SymmetricTensor<2,dim,double>::n_independent_components));
        switch (dim)
          {
          case 1:
          {
            return TableIndices<2>(0,0);
          }

          case 2:
          {
            const TableIndices<2> table[3] =
            {
              TableIndices<2> (0,0),
              TableIndices<2> (1,1),
              TableIndices<2> (0,1)
            };
            return table[i];
          }

          case 3:
          {
            const TableIndices<2> table[6] =
            {
              TableIndices<2> (0,0),
              TableIndices<2> (1,1),
              TableIndices<2> (2,2),
              TableIndices<2> (0,1),
              TableIndices<2> (0,2),
              TableIndices<2> (1,2)
            };
            return table[i];
          }

          default:
            if (i<dim)
              return TableIndices<2> (i,i);

            for (unsigned int d=0, c=0; d<dim; ++d)
              for (unsigned int e=d+1; e<dim; ++e, ++c)
                if (c==i)
                  return TableIndices<2>(d,e);

            // should never get here:
            Assert(false, ExcInternalError());
            return TableIndices<2>(0, 0);
          }
      }

      // a function to do the inverse of the unrolling from a set of
      // indices to a scalar index into the array in which we store
      // the elements of a symmetric tensor. in other words, it goes
      // from the scalar index into the array to a set of indices of
      // the tensor
      //
      // this function is for tensors of a rank not already handled
      // above
      template <int dim, int rank_>
      inline
      TableIndices<rank_>
      unrolled_to_component_indices
      (const unsigned int i,
       const std::integral_constant<int, rank_> &)
      {
        (void)i;
        Assert ((i < dealii::SymmetricTensor<rank_,dim,double>::n_independent_components),
                ExcIndexRange(i, 0, dealii::SymmetricTensor<rank_,dim,double>::n_independent_components));
        Assert (false, ExcNotImplemented());
        return TableIndices<rank_>();
      }

    }
  }
}

template <int rank_, int dim, typename Number>
inline
TableIndices<rank_>
SymmetricTensor<rank_,dim,Number>::unrolled_to_component_indices
(const unsigned int i)
{
  return
    internal::SymmetricTensorImplementation::unrolled_to_component_indices<dim> (i,
        std::integral_constant<int, rank_>());
}



template <int rank_, int dim, typename Number>
template <class Archive>
inline
void
SymmetricTensor<rank_,dim,Number>::serialize(Archive &ar, const unsigned int)
{
  ar &data;
}


#endif // DOXYGEN

/* ----------------- Non-member functions operating on tensors. ------------ */


/**
 * Addition of two symmetric tensors of equal rank. The result is another
 * SymmetricTensor that has a number type that is compatible with the
 * operation.
 *
 * If possible (e.g. when @p Number and @p OtherNumber are of the same type,
 * or if the result of <code>Number() + OtherNumber()</code> is another @p Number),
 * you should use <tt>operator +=</tt> instead since this does not require the
 * creation of a temporary variable.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
inline
SymmetricTensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
operator+(const SymmetricTensor<rank_, dim, Number>      &left,
          const SymmetricTensor<rank_, dim, OtherNumber> &right)
{
  SymmetricTensor<rank_, dim, typename ProductType<Number, OtherNumber>::type> tmp = left;
  tmp += right;
  return tmp;
}


/**
 * Subtraction of two symmetric tensors of equal rank. The result is another
 * SymmetricTensor that has a number type that is compatible with the
 * operation.
 *
 * If possible (e.g. when @p Number and @p OtherNumber are of the same type,
 * or if the result of <code>Number() + OtherNumber()</code> is another @p Number),
 * you should use <tt>operator +=</tt> instead since this does not require the
 * creation of a temporary variable.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
inline
SymmetricTensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
operator-(const SymmetricTensor<rank_, dim, Number>      &left,
          const SymmetricTensor<rank_, dim, OtherNumber> &right)
{
  SymmetricTensor<rank_, dim, typename ProductType<Number, OtherNumber>::type> tmp = left;
  tmp -= right;
  return tmp;
}


/**
 * Addition of a SymmetricTensor and a general Tensor of equal rank. The
 * result is a general Tensor that has a number type that is compatible with the
 * operation.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
inline
Tensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
operator+(const SymmetricTensor<rank_, dim, Number> &left,
          const Tensor<rank_, dim, OtherNumber> &right)
{
  return Tensor<rank_, dim, Number>(left) + right;
}


/**
 * Addition of a general Tensor with a SymmetricTensor of equal rank. The
 * result is a general Tensor that has a number type that is compatible with the
 * operation.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
inline
Tensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
operator+(const Tensor<rank_, dim, Number> &left,
          const SymmetricTensor<rank_, dim, OtherNumber> &right)
{
  return left + Tensor<rank_, dim, OtherNumber>(right);
}


/**
 * Subtraction of a general Tensor from a SymmetricTensor of equal rank. The
 * result is a general Tensor that has a number type that is compatible with the
 * operation.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
inline
Tensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
operator-(const SymmetricTensor<rank_, dim, Number> &left,
          const Tensor<rank_, dim, OtherNumber> &right)
{
  return Tensor<rank_, dim, Number>(left) - right;
}


/**
 * Subtraction of a SymmetricTensor from a general Tensor of equal rank. The
 * result is a general Tensor that has a number type that is compatible with the
 * operation.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
inline
Tensor<rank_, dim, typename ProductType<Number, OtherNumber>::type>
operator-(const Tensor<rank_, dim, Number> &left,
          const SymmetricTensor<rank_, dim, OtherNumber> &right)
{
  return left - Tensor<rank_, dim, OtherNumber>(right);
}



/**
 * Compute the determinant of a tensor or rank 2. The determinant is also
 * commonly referred to as the third invariant of rank-2 tensors.
 *
 * For a one-dimensional tensor, the determinant equals the only element and
 * is therefore equivalent to the trace.
 *
 * For greater notational simplicity, there is also a <tt>third_invariant</tt>
 * function that returns the determinant of a tensor.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Number determinant (const SymmetricTensor<2,dim,Number> &t)
{
  switch (dim)
    {
    case 1:
      return t.data[0];
    case 2:
      return (t.data[0] * t.data[1] - t.data[2]*t.data[2]);
    case 3:
    {
      // in analogy to general tensors, but
      // there's something to be simplified for
      // the present case
      const Number tmp = t.data[3]*t.data[4]*t.data[5];
      return ( tmp + tmp
               +t.data[0]*t.data[1]*t.data[2]
               -t.data[0]*t.data[5]*t.data[5]
               -t.data[1]*t.data[4]*t.data[4]
               -t.data[2]*t.data[3]*t.data[3]);
    }
    default:
      Assert (false, ExcNotImplemented());
      return internal::NumberType<Number>::value(0.0);
    }
}



/**
 * Compute the determinant of a tensor or rank 2. This function therefore
 * computes the same value as the <tt>determinant()</tt> functions and is only
 * provided for greater notational simplicity (since there are also functions
 * first_invariant() and second_invariant()).
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Number third_invariant (const SymmetricTensor<2,dim,Number> &t)
{
  return determinant (t);
}



/**
 * Compute and return the trace of a tensor of rank 2, i.e. the sum of its
 * diagonal entries. The trace is the first invariant of a rank-2 tensor.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
Number trace (const SymmetricTensor<2,dim,Number> &d)
{
  Number t = d.data[0];
  for (unsigned int i=1; i<dim; ++i)
    t += d.data[i];
  return t;
}


/**
 * Compute the trace of a tensor or rank 2. This function therefore computes
 * the same value as the <tt>trace()</tt> functions and is only provided for
 * greater notational simplicity (since there are also functions
 * second_invariant() and third_invariant()).
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
Number first_invariant (const SymmetricTensor<2,dim,Number> &t)
{
  return trace (t);
}


/**
 * Compute the second invariant of a tensor of rank 2. The second invariant of
 * a tensor $\sigma$ is defined as
 * $II(\sigma) = \frac 12 \left[ (\text{trace}\ \sigma)^2
 *                               -\text{trace}\ (\sigma^2) \right]$.
 *
 * For the kind of arguments to this function, i.e., a rank-2 tensor of
 * size 1, the result is simply zero.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005, 2010
 */
template <typename Number>
inline
Number second_invariant (const SymmetricTensor<2,1,Number> &)
{
  return internal::NumberType<Number>::value(0.0);
}



/**
 * Compute the second invariant of a tensor of rank 2. The second invariant of
 * a tensor $\sigma$ is defined as
 * $II(\sigma) = \frac 12 \left[ (\text{trace}\ \sigma)^2
 *                               -\text{trace}\ (\sigma^2) \right]$.
 *
 * For the kind of arguments to this function, i.e., a rank-2 tensor of
 * size 2, the result is (counting indices starting at one)
 * $II(\sigma) = \frac 12 \left[ (\sigma_{11} + \sigma_{22})^2
 *                               -(\sigma_{11}^2+2\sigma_{12}^2+\sigma_{22}^2) \right]
 *             = \sigma_{11}\sigma_{22} - \sigma_{12}^2$.
 * As expected, for the $2\times 2$ symmetric tensors this function handles,
 * this equals the determinant of the tensor. (This is so because for
 * $2\times 2$ symmetric tensors, there really are only two invariants, so
 * the second and third invariant are the same; the determinant is the third
 * invariant.)
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005, 2010
 */
template <typename Number>
inline
Number second_invariant (const SymmetricTensor<2,2,Number> &t)
{
  return t[0][0]*t[1][1] - t[0][1]*t[0][1];
}



/**
 * Compute the second invariant of a tensor of rank 2. The second invariant of
 * a tensor $\sigma$ is defined as
 * $II(\sigma) = \frac 12 \left[ (\text{trace}\ \sigma)^2
 *                               -\text{trace}\ (\sigma^2) \right]$.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005, 2010
 */
template <typename Number>
inline
Number second_invariant (const SymmetricTensor<2,3,Number> &t)
{
  return (t[0][0]*t[1][1] + t[1][1]*t[2][2] + t[2][2]*t[0][0]
          - t[0][1]*t[0][1] - t[0][2]*t[0][2] - t[1][2]*t[1][2]);
}



/**
 * Return the eigenvalues of a symmetric 1x1 tensor of rank 2.
 * The (single) entry of the tensor is, of course, equal to the (single)
 * eigenvalue.
 *
 * @relatesalso SymmetricTensor
 * @author Jean-Paul Pelteret, 2017
 */
template <typename Number>
std::array<Number,1>
eigenvalues (const SymmetricTensor<2,1,Number> &T);



/**
 * Return the eigenvalues of a symmetric 2x2 tensor of rank 2.
 * The array of eigenvalues is sorted in descending order.
 *
 * For 2x2 tensors, the eigenvalues of tensor $T$ are the roots of
 * <a href="https://en.wikipedia.org/wiki/Eigenvalue_algorithm#2.C3.972_matrices">the characteristic polynomial</a>
 * $0 = \lambda^{2} - \lambda\textrm{tr}(T) + \textrm{det}(T)$
 * as given by
 * $\lambda = \frac{\textrm{tr}(T) \pm \sqrt{[\textrm{tr}(T)]^{2} - 4\textrm{det}(T)}}{2}$.
 *
 * @warning The algorithm employed here determines the eigenvalues by
 * computing the roots of the characteristic polynomial. In the case that there
 * exists a common root (the eigenvalues are equal), the computation is
 * <a href="https://scicomp.stackexchange.com/q/23686">subject to round-off errors</a>
 * of order $\sqrt{\epsilon}$.
 * As an alternative, the eigenvectors() function provides a more robust, but costly,
 * method to compute the eigenvalues of a symmetric tensor.
 *
 * @relatesalso SymmetricTensor
 * @author Jean-Paul Pelteret, 2017
 */
template <typename Number>
std::array<Number,2>
eigenvalues (const SymmetricTensor<2,2,Number> &T);



/**
 * Return the eigenvalues of a symmetric 3x3 tensor of rank 2.
 * The array of eigenvalues is sorted in descending order.
 *
 * For 3x3 tensors, the eigenvalues of tensor $T$ are the roots of
 * <a href="https://en.wikipedia.org/wiki/Eigenvalue_algorithm#3.C3.973_matrices">the characteristic polynomial</a>
 * $0 = \lambda^{3} - \lambda^{2}\textrm{tr}(T) - \frac{1}{2} \lambda (\textrm{tr}(T^{2}) - [\textrm{tr}(T)]^{2}) - \textrm{det}(T)$.
 *
 * @warning The algorithm employed here determines the eigenvalues by
 * computing the roots of the characteristic polynomial. In the case that there
 * exists a common root (the eigenvalues are equal), the computation is
 * <a href="https://scicomp.stackexchange.com/q/23686">subject to round-off errors</a>
 * of order $\sqrt{\epsilon}$.
 * As an alternative, the eigenvectors() function provides a more robust, but costly,
 * method to compute the eigenvalues of a symmetric tensor.
 *
 * @relatesalso SymmetricTensor
 * @author Jean-Paul Pelteret, 2017
 */
template <typename Number>
std::array<Number,3>
eigenvalues (const SymmetricTensor<2,3,Number> &T);



namespace internal
{
  namespace SymmetricTensorImplementation
  {
    /**
     * Tridiagonalize a rank-2 symmetric tensor using the Householder method.
     * The specialized algorithm implemented here is given in
     * @code{.bib}
     * @Article{Kopp2008,
     *   title        = {Efficient numerical diagonalization of hermitian 3x3 matrices},
     *   author       = {Kopp, J.},
     *   journal      = {International Journal of Modern Physics C},
     *   year         = {2008},
     *   volume       = {19},
     *   number       = {3},
     *   pages        = {523--548},
     *   doi          = {10.1142/S0129183108012303},
     *   eprinttype   = {arXiv},
     *   eprint       = {physics/0610206v3},
     *   eprintclass  = {physics.comp-ph},
     *   url          = {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
     * }
     * @endcode
     * and is based off of the generic algorithm presented in section 11.3.2 of
     * @code{.bib}
     * @Book{Press2007,
     *   title        = {Numerical recipes 3rd edition: The art of scientific computing},
     *   author       = {Press, W. H.},
     *   journal      = {Cambridge university press},
     *   year         = {2007}
     * }
     * @endcode
     *
     * @param[in]  A This tensor to be tridiagonalized
     * @param[out] Q The orthogonal matrix effecting the transformation
     * @param[out] d The diagonal elements of the tridiagonal matrix
     * @param[out] e The off-diagonal elements of the tridiagonal matrix
     *
     * @author Joachim Kopp, Jean-Paul Pelteret, 2017
     */
    template <int dim, typename Number>
    void
    tridiagonalize (const dealii::SymmetricTensor<2,dim,Number> &A,
                    dealii::Tensor<2,dim,Number>                &Q,
                    std::array<Number,dim>                      &d,
                    std::array<Number,dim-1>                    &e);



    /**
     * Compute the eigenvalues and eigenvectors of a real-valued rank-2
     * symmetric tensor using the QL algorithm with implicit shifts.
     * The specialized algorithm implemented here is given in
     * @code{.bib}
     * @Article{Kopp2008,
     *   title        = {Efficient numerical diagonalization of hermitian 3x3 matrices},
     *   author       = {Kopp, J.},
     *   journal      = {International Journal of Modern Physics C},
     *   year         = {2008},
     *   volume       = {19},
     *   number       = {3},
     *   pages        = {523--548},
     *   doi          = {10.1142/S0129183108012303},
     *   eprinttype   = {arXiv},
     *   eprint       = {physics/0610206v3},
     *   eprintclass  = {physics.comp-ph},
     *   url          = {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
     * }
     * @endcode
     * and is based off of the generic algorithm presented in section 11.4.3 of
     * @code{.bib}
     * @Book{Press2007,
     *   title        = {Numerical recipes 3rd edition: The art of scientific computing},
     *   author       = {Press, W. H.},
     *   journal      = {Cambridge university press},
     *   year         = {2007}
     * }
     * @endcode
     *
     * @param[in] A The tensor of which the eigenvectors and eigenvalues are
     * to be computed.
     *
     * @return An array containing the eigenvectors and the associated eigenvalues.
     * The array is not sorted in any particular order.
     *
     * @author Joachim Kopp, Jean-Paul Pelteret, 2017
     */
    template <int dim, typename Number>
    std::array<std::pair<Number, Tensor<1,dim,Number> >,dim>
    ql_implicit_shifts (const dealii::SymmetricTensor<2,dim,Number> &A);



    /**
     * Compute the eigenvalues and eigenvectors of a real-valued rank-2
     * symmetric tensor using the Jacobi algorithm.
     * The specialized algorithm implemented here is given in
     * @code{.bib}
     * @Article{Kopp2008,
     *   title        = {Efficient numerical diagonalization of hermitian 3x3 matrices},
     *   author       = {Kopp, J.},
     *   journal      = {International Journal of Modern Physics C},
     *   year         = {2008},
     *   volume       = {19},
     *   number       = {3},
     *   pages        = {523--548},
     *   doi          = {10.1142/S0129183108012303},
     *   eprinttype   = {arXiv},
     *   eprint       = {physics/0610206v3},
     *   eprintclass  = {physics.comp-ph},
     *   url          = {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
     * }
     * @endcode
     * and is based off of the generic algorithm presented in section 11.4.3 of
     * @code{.bib}
     * @Book{Press2007,
     *   title        = {Numerical recipes 3rd edition: The art of scientific computing},
     *   author       = {Press, W. H.},
     *   journal      = {Cambridge university press},
     *   year         = {2007}
     * }
     * @endcode
     *
     * @param[in] A The tensor of which the eigenvectors and eigenvalues are
     * to be computed.
     *
     * @return An array containing the eigenvectors and the associated eigenvalues.
     * The array is not sorted in any particular order.
     *
     * @author Joachim Kopp, Jean-Paul Pelteret, 2017
     */
    template <int dim, typename Number>
    std::array<std::pair<Number, Tensor<1,dim,Number> >,dim>
    jacobi (dealii::SymmetricTensor<2,dim,Number> A);



    /**
     * Compute the eigenvalues and eigenvectors of a real-valued rank-2
     * symmetric 2x2 tensor using the characteristic equation to compute eigenvalues
     * and an analytical approach based on the cross-product for the eigenvectors.
     * If the computations are deemed too inaccurate then the method falls
     * back to ql_implicit_shifts.
     *
     * @param[in] A The tensor of which the eigenvectors and eigenvalues are
     * to be computed.
     *
     * @return An array containing the eigenvectors and the associated eigenvalues.
     * The array is not sorted in any particular order.
     *
     * @author Joachim Kopp, Jean-Paul Pelteret, 2017
     */
    template <typename Number>
    std::array<std::pair<Number, Tensor<1,2,Number> >,2>
    hybrid (const dealii::SymmetricTensor<2,2,Number> &A);



    /**
     * Compute the eigenvalues and eigenvectors of a real-valued rank-2
     * symmetric 3x3 tensor using the characteristic equation to compute eigenvalues
     * and an analytical approach based on the cross-product for the eigenvectors.
     * If the computations are deemed too inaccurate then the method falls
     * back to ql_implicit_shifts.
     * The specialized algorithm implemented here is given in
     * @code{.bib}
     * @Article{Kopp2008,
     *   title        = {Efficient numerical diagonalization of hermitian 3x3 matrices},
     *   author       = {Kopp, J.},
     *   journal      = {International Journal of Modern Physics C},
     *   year         = {2008},
     *   volume       = {19},
     *   number       = {3},
     *   pages        = {523--548},
     *   doi          = {10.1142/S0129183108012303},
     *   eprinttype   = {arXiv},
     *   eprint       = {physics/0610206v3},
     *   eprintclass  = {physics.comp-ph},
     *   url          = {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
     * }
     * @endcode
     *
     * @param[in] A The tensor of which the eigenvectors and eigenvalues are
     * to be computed.
     *
     * @return An array containing the eigenvectors and the associated eigenvalues.
     * The array is not sorted in any particular order.
     *
     * @author Joachim Kopp, Jean-Paul Pelteret, 2017
     */
    template <typename Number>
    std::array<std::pair<Number, Tensor<1,3,Number> >,3>
    hybrid (const dealii::SymmetricTensor<2,3,Number> &A);

    namespace
    {

      /**
       * A struct that is used to sort arrays of pairs of eign=envalues and
       * eigenvectors. Sorting is performed in in descending order of eigenvalue.
       */
      template <int dim, typename Number>
      struct SortEigenValuesVectors
      {
        typedef std::pair<Number, Tensor<1,dim,Number> > EigValsVecs;
        bool operator() (const EigValsVecs &lhs,
                         const EigValsVecs &rhs)
        {
          return lhs.first > rhs.first;
        }
      };

    }

  } // namespace SymmetricTensor

} // namespace internal



// The line below is to ensure that doxygen puts the full description
// of this global enumeration into the documentation
// See https://stackoverflow.com/a/1717984
/** @file */
/**
 * An enumeration for the algorithm to be employed when performing
 * the computation of normalized eigenvectors and their corresponding
 * eigenvalues by the eigenvalues() and eigenvectors() methods operating on
 * SymmetricTensor objects.
 *
 * The specialized algorithms utilized in computing the eigenvectors are
 * presented in
 * @code{.bib}
 * @Article{Kopp2008,
 *   title        = {Efficient numerical diagonalization of hermitian 3x3 matrices},
 *   author       = {Kopp, J.},
 *   journal      = {International Journal of Modern Physics C},
 *   year         = {2008},
 *   volume       = {19},
 *   number       = {3},
 *   pages        = {523--548},
 *   doi          = {10.1142/S0129183108012303},
 *   eprinttype   = {arXiv},
 *   eprint       = {physics/0610206v3},
 *   eprintclass  = {physics.comp-ph},
 *   url          = {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
 * }
 * @endcode
 */
enum struct SymmetricTensorEigenvectorMethod
{
  /**
   * A hybrid approach that preferentially uses the characteristic equation to
   * compute eigenvalues and an analytical approach based on the cross-product
   * for the eigenvectors. If the computations are deemed too inaccurate then
   * the method falls back to ql_implicit_shifts.
   *
   * This method potentially offers the quickest computation if the pathological
   * case is not encountered.
   */
  hybrid,
  /**
   * The iterative QL algorithm with implicit shifts applied after
   * tridiagonalization of the tensor using the householder method.
   *
   * This method offers a compromise between speed of computation and its
   * robustness. This method is particularly useful when the elements
   * of $T$ have greatly varying magnitudes, which would typically lead to a
   * loss of accuracy when computing the smaller eigenvalues.
   */
  ql_implicit_shifts,
  /**
   * The iterative Jacobi algorithm.
   *
   * This method offers is the most robust of the available options, with
   * reliable results obtained for even the most pathological cases. It is,
   * however, the slowest algorithm of all of those implemented.
   */
  jacobi
};



/**
 * Return the eigenvalues and eigenvectors of a real-valued rank-2 symmetric
 * tensor $T$. The array of matched eigenvalue and eigenvector pairs is sorted
 * in descending order (determined by the eigenvalues).
 *
 * The specialized algorithms utilized in computing the eigenvectors are
 * presented in
 * @code{.bib}
 * @Article{Kopp2008,
 *   title        = {Efficient numerical diagonalization of hermitian 3x3 matrices},
 *   author       = {Kopp, J.},
 *   journal      = {International Journal of Modern Physics C},
 *   year         = {2008},
 *   volume       = {19},
 *   number       = {3},
 *   pages        = {523--548},
 *   doi          = {10.1142/S0129183108012303},
 *   eprinttype   = {arXiv},
 *   eprint       = {physics/0610206v3},
 *   eprintclass  = {physics.comp-ph},
 *   url          = {https://www.mpi-hd.mpg.de/personalhomes/globes/3x3/index.html}
 * }
 * @endcode
 *
 * @relatesalso SymmetricTensor
 * @author Joachim Kopp, Jean-Paul Pelteret, 2017
 */
template <int dim, typename Number>
std::array<std::pair<Number, Tensor<1,dim,Number> >,std::integral_constant<int, dim>::value>
eigenvectors (const SymmetricTensor<2,dim,Number>         &T,
              const SymmetricTensorEigenvectorMethod       method = SymmetricTensorEigenvectorMethod::ql_implicit_shifts);



/**
 * Return the transpose of the given symmetric tensor. Since we are working
 * with symmetric objects, the transpose is of course the same as the original
 * tensor. This function mainly exists for compatibility with the Tensor
 * class.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int rank_, int dim, typename Number>
inline
SymmetricTensor<rank_,dim,Number>
transpose (const SymmetricTensor<rank_,dim,Number> &t)
{
  return t;
}



/**
 * Compute the deviator of a symmetric tensor, which is defined as <tt>dev[s]
 * = s - 1/dim*tr[s]*I</tt>, where <tt>I</tt> is the identity operator. This
 * quantity equals the original tensor minus its contractive or dilative
 * component and refers to the shear in, for example, elasticity.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
SymmetricTensor<2,dim,Number>
deviator (const SymmetricTensor<2,dim,Number> &t)
{
  SymmetricTensor<2,dim,Number> tmp = t;

  // subtract scaled trace from the diagonal
  const Number tr = trace(t) / dim;
  for (unsigned int i=0; i<dim; ++i)
    tmp.data[i] -= tr;

  return tmp;
}



/**
 * Return a unit symmetric tensor of rank 2, i.e., the dim-by-dim identity
 * matrix.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
SymmetricTensor<2,dim,Number>
unit_symmetric_tensor ()
{
  // create a default constructed matrix filled with
  // zeros, then set the diagonal elements to one
  SymmetricTensor<2,dim,Number> tmp;
  switch (dim)
    {
    case 1:
      tmp.data[0] = 1;
      break;
    case 2:
      tmp.data[0] = tmp.data[1] = 1;
      break;
    case 3:
      tmp.data[0] = tmp.data[1] = tmp.data[2] = 1;
      break;
    default:
      for (unsigned int d=0; d<dim; ++d)
        tmp.data[d] = 1;
    }
  return tmp;
}



/**
 * Return a unit symmetric tensor of rank 2, i.e., the dim-by-dim identity
 * matrix. This specialization of the function uses <code>double</code> as the
 * data type for the elements.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
SymmetricTensor<2,dim>
unit_symmetric_tensor ()
{
  return unit_symmetric_tensor<dim,double>();
}



/**
 * Return the tensor of rank 4 that, when multiplied by a symmetric rank 2
 * tensor <tt>t</tt> returns the deviator $\textrm{dev}\ t$. It is the
 * operator representation of the linear deviator operator.
 *
 * For every tensor <tt>t</tt>, there holds the identity
 * <tt>deviator(t)==deviator_tensor&lt;dim&gt;()*t</tt>, up to numerical
 * round-off. The reason this operator representation is provided is that one
 * sometimes needs to invert operators like <tt>identity_tensor&lt;dim&gt;() +
 * delta_t*deviator_tensor&lt;dim&gt;()</tt> or similar.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
SymmetricTensor<4,dim,Number>
deviator_tensor ()
{
  SymmetricTensor<4,dim,Number> tmp;

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
       i<internal::SymmetricTensorAccessors::StorageType<4,dim,Number>::n_rank2_components;
       ++i)
    tmp.data[i][i] = 0.5;

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
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
SymmetricTensor<4,dim>
deviator_tensor ()
{
  return deviator_tensor<dim,double>();
}



/**
 * Return the fourth-order symmetric identity tensor which maps symmetric
 * second-order tensors to themselves.
 *
 * Note that this tensor, even though it is the identity, has a somewhat funny
 * form, and in particular does not only consist of zeros and ones. For
 * example, for <tt>dim=2</tt>, the identity tensor has all zero entries
 * except for <tt>id[0][0][0][0]=id[1][1][1][1]=1</tt> and
 * <tt>id[0][1][0][1]=id[0][1][1][0]=id[1][0][0][1]=id[1][0][1][0]=1/2</tt>.
 * To see why this factor of 1/2 is necessary, consider computing <tt>A=Id :
 * B</tt>. For the element <tt>a_01</tt> we have <tt>a_01=id_0100 b_00 +
 * id_0111 b_11 + id_0101 b_01 + id_0110 b_10</tt>. On the other hand, we need
 * to have <tt>a_01=b_01</tt>, and symmetry implies <tt>b_01=b_10</tt>,
 * leading to <tt>a_01=(id_0101+id_0110) b_01</tt>, or, again by symmetry,
 * <tt>id_0101=id_0110=1/2</tt>. Similar considerations hold for the three-
 * dimensional case.
 *
 * This issue is also explained in the introduction to step-44.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
SymmetricTensor<4,dim,Number>
identity_tensor ()
{
  SymmetricTensor<4,dim,Number> tmp;

  // fill the elements treating the diagonal
  for (unsigned int i=0; i<dim; ++i)
    tmp.data[i][i] = 1;

  // then fill the ones that copy over the
  // non-diagonal elements. note that during
  // the double-contraction, we handle the
  // off-diagonal elements twice, so simply
  // copying requires a weight of 1/2
  for (unsigned int i=dim;
       i<internal::SymmetricTensorAccessors::StorageType<4,dim,Number>::n_rank2_components;
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
 * <tt>id[0][1][0][1]=id[0][1][1][0]=id[1][0][0][1]=id[1][0][1][0]=1/2</tt>.
 * To see why this factor of 1/2 is necessary, consider computing <tt>A=Id .
 * B</tt>. For the element <tt>a_01</tt> we have <tt>a_01=id_0100 b_00 +
 * id_0111 b_11 + id_0101 b_01 + id_0110 b_10</tt>. On the other hand, we need
 * to have <tt>a_01=b_01</tt>, and symmetry implies <tt>b_01=b_10</tt>,
 * leading to <tt>a_01=(id_0101+id_0110) b_01</tt>, or, again by symmetry,
 * <tt>id_0101=id_0110=1/2</tt>. Similar considerations hold for the three-
 * dimensional case.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim>
inline
SymmetricTensor<4,dim>
identity_tensor ()
{
  return identity_tensor<dim,double>();
}



/**
 * Invert a symmetric rank-2 tensor.
 *
 * @note If a tensor is not invertible, then the result is unspecified, but will
 * likely contain the results of a division by zero or a very small number at
 * the very least.
 *
 * @relatesalso SymmetricTensor
 * @author Jean-Paul Pelteret, 2016
 */
template <int dim, typename Number>
inline
SymmetricTensor<2,dim,Number>
invert (const SymmetricTensor<2,dim,Number> &t)
{
  return internal::SymmetricTensorImplementation::Inverse<2,dim,Number>::value(t);
}



/**
 * Invert a symmetric rank-4 tensor. Since symmetric rank-4 tensors are
 * mappings from and to symmetric rank-2 tensors, they can have an inverse.
 *
 * If a tensor is not invertible, then the result is unspecified, but will
 * likely contain the results of a division by zero or a very small number at
 * the very least.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
SymmetricTensor<4,dim,Number>
invert (const SymmetricTensor<4,dim,Number> &t)
{
  return internal::SymmetricTensorImplementation::Inverse<4,dim,Number>::value(t);
}



/**
 * Return the tensor of rank 4 that is the outer product of the two tensors
 * given as arguments, i.e. the result $T=t1 \otimes t2$ satisfies <tt>T phi =
 * t1 (t2 : phi)</tt> for all symmetric tensors <tt>phi</tt>.
 *
 * For example, the deviator tensor can be computed as
 * <tt>identity_tensor<dim,Number>() -
 * 1/d*outer_product(unit_symmetric_tensor<dim,Number>(),
 * unit_symmetric_tensor<dim,Number>())</tt>, since the (double) contraction
 * with the unit tensor yields the trace of a symmetric tensor.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number>
inline
SymmetricTensor<4,dim,Number>
outer_product (const SymmetricTensor<2,dim,Number> &t1,
               const SymmetricTensor<2,dim,Number> &t2)
{
  SymmetricTensor<4,dim,Number> tmp;

  // fill only the elements really needed
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=i; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=k; l<dim; ++l)
          tmp[i][j][k][l] = t1[i][j] * t2[k][l];

  return tmp;
}



/**
 * Return the symmetrized version of a full rank-2 tensor, i.e.
 * (t+transpose(t))/2, as a symmetric rank-2 tensor. This is the version for
 * general dimensions.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim,typename Number>
inline
SymmetricTensor<2,dim,Number>
symmetrize (const Tensor<2,dim,Number> &t)
{
  Number array[(dim*dim+dim)/2];
  for (unsigned int d=0; d<dim; ++d)
    array[d] = t[d][d];
  for (unsigned int d=0, c=0; d<dim; ++d)
    for (unsigned int e=d+1; e<dim; ++e, ++c)
      array[dim+c] = (t[d][e]+t[e][d])*0.5;
  return SymmetricTensor<2,dim,Number>(array);
}



/**
 * Multiplication of a symmetric tensor of general rank with a scalar from the
 * right. This version of the operator is used if the scalar has the same data
 * type as is used to store the elements of the symmetric tensor.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim, typename Number>
inline
SymmetricTensor<rank_,dim,Number>
operator * (const SymmetricTensor<rank_,dim,Number> &t,
            const Number                           &factor)
{
  SymmetricTensor<rank_,dim,Number> tt = t;
  tt *= factor;
  return tt;
}



/**
 * Multiplication of a symmetric tensor of general rank with a scalar from the
 * left. This version of the operator is used if the scalar has the same data
 * type as is used to store the elements of the symmetric tensor.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim, typename Number>
inline
SymmetricTensor<rank_,dim,Number>
operator * (const Number                           &factor,
            const SymmetricTensor<rank_,dim,Number> &t)
{
  // simply forward to the other operator
  return t*factor;
}



/**
 * Multiplication of a symmetric tensor with a scalar number from the right.
 *
 * The purpose of this operator is to enable only multiplication of a tensor
 * by a scalar number (i.e., a floating point number, a complex floating point
 * number, etc.). The function is written in a way that only allows the
 * compiler to consider the function if the second argument is indeed a scalar
 * number -- in other words, @p OtherNumber will not match, for example
 * <code>std::vector@<double@></code> as the product of a tensor and a vector
 * clearly would make no sense. The mechanism by which the compiler is
 * prohibited of considering this operator for multiplication with non-scalar
 * types are explained in the documentation of the EnableIfScalar class.
 *
 * The return type of the function is chosen so that it matches the types of
 * both the tensor and the scalar argument. For example, if you multiply a
 * <code>SymmetricTensor@<2,dim,double@></code> by
 * <code>std::complex@<double@></code>, then the result will be a
 * <code>SymmetricTensor@<2,dim,std::complex@<double@>@></code>. In other
 * words, the type with which the returned tensor stores its components equals
 * the type you would get if you multiplied an individual component of the
 * input tensor by the scalar factor.
 *
 * @relatesalso SymmetricTensor
 * @relatesalso EnableIfScalar
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
inline
SymmetricTensor<rank_,dim,typename ProductType<Number,typename EnableIfScalar<OtherNumber>::type>::type>
operator * (const SymmetricTensor<rank_,dim,Number> &t,
            const OtherNumber                      &factor)
{
  // form the product. we have to convert the two factors into the final
  // type via explicit casts because, for awkward reasons, the C++
  // standard committee saw it fit to not define an
  //   operator*(float,std::complex<double>)
  // (as well as with switched arguments and double<->float).
  typedef typename ProductType<Number,OtherNumber>::type product_type;
  SymmetricTensor<rank_,dim,product_type> tt(t);
  // we used to shorten the following by 'tt *= product_type(factor);'
  // which requires that a converting constructor
  // 'product_type::product_type(const OtherNumber) is defined.
  // however, a user-defined constructor is not allowed for aggregates,
  // e.g. VectorizedArray. therefore, we work around this issue using a
  // copy-assignment operator 'product_type::operator=(const OtherNumber)'
  // which we assume to be defined.
  product_type new_factor;
  new_factor = factor;
  tt *= new_factor;
  return tt;
}



/**
 * Multiplication of a symmetric tensor with a scalar number from the left.
 * See the discussion with the operator with switched arguments for more
 * information about template arguments and the return type.
 *
 * @relatesalso SymmetricTensor
 * @relatesalso EnableIfScalar
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
inline
SymmetricTensor<rank_,dim,typename ProductType<OtherNumber,typename EnableIfScalar<Number>::type>::type>
operator * (const Number                                &factor,
            const SymmetricTensor<rank_,dim,OtherNumber> &t)
{
  // simply forward to the other operator with switched arguments
  return (t*factor);
}



/**
 * Division of a symmetric tensor of general rank by a scalar.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim, typename Number, typename OtherNumber>
inline
SymmetricTensor<rank_,dim,typename ProductType<Number,typename EnableIfScalar<OtherNumber>::type>::type>
operator / (const SymmetricTensor<rank_,dim,Number> &t,
            const OtherNumber                      &factor)
{
  SymmetricTensor<rank_,dim,typename ProductType<Number,OtherNumber>::type> tt = t;
  tt /= factor;
  return tt;
}



/**
 * Multiplication of a symmetric tensor of general rank with a scalar from the
 * right.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim>
inline
SymmetricTensor<rank_,dim>
operator * (const SymmetricTensor<rank_,dim> &t,
            const double                     factor)
{
  SymmetricTensor<rank_,dim> tt = t;
  tt *= factor;
  return tt;
}



/**
 * Multiplication of a symmetric tensor of general rank with a scalar from the
 * left.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim>
inline
SymmetricTensor<rank_,dim>
operator * (const double                     factor,
            const SymmetricTensor<rank_,dim> &t)
{
  SymmetricTensor<rank_,dim> tt = t;
  tt *= factor;
  return tt;
}



/**
 * Division of a symmetric tensor of general rank by a scalar.
 *
 * @relatesalso SymmetricTensor
 */
template <int rank_, int dim>
inline
SymmetricTensor<rank_,dim>
operator / (const SymmetricTensor<rank_,dim> &t,
            const double                     factor)
{
  SymmetricTensor<rank_,dim> tt = t;
  tt /= factor;
  return tt;
}

/**
 * Compute the scalar product $a:b=\sum_{i,j} a_{ij}b_{ij}$ between two
 * tensors $a,b$ of rank 2. In the current case where both arguments are
 * symmetric tensors, this is equivalent to calling the expression
 * <code>t1*t2</code> which uses the overloaded <code>operator*</code> between
 * two symmetric tensors of rank 2.
 *
 * @relatesalso SymmetricTensor
 */
template <int dim, typename Number, typename OtherNumber>
inline
typename ProductType<Number, OtherNumber>::type
scalar_product (const SymmetricTensor<2,dim,Number>      &t1,
                const SymmetricTensor<2,dim,OtherNumber> &t2)
{
  return (t1*t2);
}


/**
 * Compute the scalar product $a:b=\sum_{i,j} a_{ij}b_{ij}$ between two
 * tensors $a,b$ of rank 2. We don't use <code>operator*</code> for this
 * operation since the product between two tensors is usually assumed to be
 * the contraction over the last index of the first tensor and the first index
 * of the second tensor, for example $(a\cdot b)_{ij}=\sum_k a_{ik}b_{kj}$.
 *
 * @relatesalso Tensor @relatesalso SymmetricTensor
 */
template <int dim, typename Number, typename OtherNumber>
inline
typename ProductType<Number, OtherNumber>::type
scalar_product (const SymmetricTensor<2,dim,Number> &t1,
                const Tensor<2,dim,OtherNumber>     &t2)
{
  typename ProductType<Number, OtherNumber>::type s = internal::NumberType<typename ProductType<Number, OtherNumber>::type>::value(0.0);
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      s += t1[i][j] * t2[i][j];
  return s;
}


/**
 * Compute the scalar product $a:b=\sum_{i,j} a_{ij}b_{ij}$ between two
 * tensors $a,b$ of rank 2. We don't use <code>operator*</code> for this
 * operation since the product between two tensors is usually assumed to be
 * the contraction over the last index of the first tensor and the first index
 * of the second tensor, for example $(a\cdot b)_{ij}=\sum_k a_{ik}b_{kj}$.
 *
 * @relatesalso Tensor @relatesalso SymmetricTensor
 */
template <int dim, typename Number, typename OtherNumber>
inline
typename ProductType<Number, OtherNumber>::type
scalar_product (const Tensor<2,dim,Number>               &t1,
                const SymmetricTensor<2,dim,OtherNumber> &t2)
{
  return scalar_product(t2, t1);
}


/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in the symmetric tensor of rank 2 that is given as first argument
 * to this function. This operation is the symmetric tensor analogon of a
 * matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the SymmetricTensor
 * class. It should not be used, however, since the member operator has
 * knowledge of the actual data storage format and is at least 2 orders of
 * magnitude faster. This function mostly exists for compatibility purposes
 * with the general tensor class.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <typename Number, typename OtherNumber>
inline
void
double_contract (SymmetricTensor<2,1,typename ProductType<Number, OtherNumber>::type> &tmp,
                 const SymmetricTensor<4,1,Number>      &t,
                 const SymmetricTensor<2,1,OtherNumber> &s)
{
  tmp[0][0] = t[0][0][0][0] * s[0][0];
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in the symmetric tensor of rank 2 that is given as first argument
 * to this function. This operation is the symmetric tensor analogon of a
 * matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the SymmetricTensor
 * class. It should not be used, however, since the member operator has
 * knowledge of the actual data storage format and is at least 2 orders of
 * magnitude faster. This function mostly exists for compatibility purposes
 * with the general tensor class.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <typename Number, typename OtherNumber>
inline
void
double_contract (SymmetricTensor<2,1,typename ProductType<Number, OtherNumber>::type> &tmp,
                 const SymmetricTensor<2,1,Number>      &s,
                 const SymmetricTensor<4,1,OtherNumber> &t)
{
  tmp[0][0] = t[0][0][0][0] * s[0][0];
}



/**
 * Double contraction between a rank-4 and a rank-2 symmetric tensor,
 * resulting in the symmetric tensor of rank 2 that is given as first argument
 * to this function. This operation is the symmetric tensor analogon of a
 * matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the SymmetricTensor
 * class. It should not be used, however, since the member operator has
 * knowledge of the actual data storage format and is at least 2 orders of
 * magnitude faster. This function mostly exists for compatibility purposes
 * with the general tensor class.
 *
 * @relatesalso SymmetricTensor @author Wolfgang Bangerth, 2005
 */
template <typename Number, typename OtherNumber>
inline
void
double_contract (SymmetricTensor<2,2,typename ProductType<Number, OtherNumber>::type> &tmp,
                 const SymmetricTensor<4,2,Number>      &t,
                 const SymmetricTensor<2,2,OtherNumber> &s)
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
 * resulting in the symmetric tensor of rank 2 that is given as first argument
 * to this function. This operation is the symmetric tensor analogon of a
 * matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the SymmetricTensor
 * class. It should not be used, however, since the member operator has
 * knowledge of the actual data storage format and is at least 2 orders of
 * magnitude faster. This function mostly exists for compatibility purposes
 * with the general tensor class.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <typename Number, typename OtherNumber>
inline
void
double_contract (SymmetricTensor<2,2,typename ProductType<Number, OtherNumber>::type> &tmp,
                 const SymmetricTensor<2,2,Number>      &s,
                 const SymmetricTensor<4,2,OtherNumber> &t)
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
 * resulting in the symmetric tensor of rank 2 that is given as first argument
 * to this function. This operation is the symmetric tensor analogon of a
 * matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the SymmetricTensor
 * class. It should not be used, however, since the member operator has
 * knowledge of the actual data storage format and is at least 2 orders of
 * magnitude faster. This function mostly exists for compatibility purposes
 * with the general tensor class.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <typename Number, typename OtherNumber>
inline
void
double_contract (SymmetricTensor<2,3,typename ProductType<Number, OtherNumber>::type> &tmp,
                 const SymmetricTensor<4,3,Number>      &t,
                 const SymmetricTensor<2,3,OtherNumber> &s)
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
 * resulting in the symmetric tensor of rank 2 that is given as first argument
 * to this function. This operation is the symmetric tensor analogon of a
 * matrix-vector multiplication.
 *
 * This function does the same as the member operator* of the SymmetricTensor
 * class. It should not be used, however, since the member operator has
 * knowledge of the actual data storage format and is at least 2 orders of
 * magnitude faster. This function mostly exists for compatibility purposes
 * with the general tensor class.
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <typename Number, typename OtherNumber>
inline
void
double_contract (SymmetricTensor<2,3,typename ProductType<Number, OtherNumber>::type> &tmp,
                 const SymmetricTensor<2,3,Number>      &s,
                 const SymmetricTensor<4,3,OtherNumber> &t)
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



/**
 * Multiply a symmetric rank-2 tensor (i.e., a matrix) by a rank-1 tensor
 * (i.e., a vector). The result is a rank-1 tensor (i.e., a vector).
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number, typename OtherNumber>
Tensor<1,dim,typename ProductType<Number,OtherNumber>::type>
operator * (const SymmetricTensor<2,dim,Number> &src1,
            const Tensor<1,dim,OtherNumber>     &src2)
{
  Tensor<1,dim,typename ProductType<Number,OtherNumber>::type> dest;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      dest[i] += src1[i][j] * src2[j];
  return dest;
}


/**
 * Multiply a rank-1 tensor (i.e., a vector) by a symmetric rank-2 tensor
 * (i.e., a matrix). The result is a rank-1 tensor (i.e., a vector).
 *
 * @relatesalso SymmetricTensor
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, typename Number, typename OtherNumber>
Tensor<1,dim,typename ProductType<Number,OtherNumber>::type>
operator * (const Tensor<1,dim,Number>               &src1,
            const SymmetricTensor<2,dim,OtherNumber> &src2)
{
  // this is easy for symmetric tensors:
  return src2 * src1;
}



/**
 * The dot product (single contraction) for tensors: Return a tensor of rank
 * $(\text{rank}_1 + \text{rank}_2 - 2)$ that is the contraction of the last
 * index of a tensor @p src1 of rank @p rank_1 with the first index of a
 * tensor @p src2 of rank @p rank_2:
 * @f[
 *   \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 *   = \sum_{k}
 *     \text{left}_{i_1,\ldots,i_{r1}, k}
 *     \text{right}_{k, j_1,\ldots,j_{r2}}
 * @f]
 *
 * @note As one operand is a Tensor, the multiplication operator only performs a
 * contraction over a single pair of indices. This is in contrast to the
 * multiplication operator for SymmetricTensor, which does the double
 * contraction.
 *
 * @relatesalso SymmetricTensor
 * @author Matthias Maier, Jean-Paul Pelteret, 2017
 */
template <int rank_1, int rank_2, int dim,
          typename Number, typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
typename Tensor<rank_1 + rank_2 - 2, dim, typename ProductType<Number, OtherNumber>::type>::tensor_type
operator * (const Tensor<rank_1, dim, Number>               &src1,
            const SymmetricTensor<rank_2, dim, OtherNumber> &src2s)
{
  typename Tensor<rank_1 + rank_2 - 2, dim, typename ProductType<Number, OtherNumber>::type>::tensor_type result;
  const Tensor<rank_2, dim, OtherNumber> src2 (src2s);
  return src1*src2;
}



/**
 * The dot product (single contraction) for tensors: Return a tensor of rank
 * $(\text{rank}_1 + \text{rank}_2 - 2)$ that is the contraction of the last
 * index of a tensor @p src1 of rank @p rank_1 with the first index of a
 * tensor @p src2 of rank @p rank_2:
 * @f[
 *   \text{result}_{i_1,\ldots,i_{r1},j_1,\ldots,j_{r2}}
 *   = \sum_{k}
 *     \text{left}_{i_1,\ldots,i_{r1}, k}
 *     \text{right}_{k, j_1,\ldots,j_{r2}}
 * @f]
 *
 * @note As one operand is a Tensor, the multiplication operator only performs a
 * contraction over a single pair of indices. This is in contrast to the
 * multiplication operator for SymmetricTensor, which does the double
 * contraction.
 *
 * @relatesalso SymmetricTensor
 * @author Matthias Maier, Jean-Paul Pelteret, 2017
 */
template <int rank_1, int rank_2, int dim,
          typename Number, typename OtherNumber>
inline DEAL_II_ALWAYS_INLINE
typename Tensor<rank_1 + rank_2 - 2, dim, typename ProductType<Number, OtherNumber>::type>::tensor_type
operator * (const SymmetricTensor<rank_1, dim, Number> &src1s,
            const Tensor<rank_2, dim, OtherNumber>     &src2)
{
  typename Tensor<rank_1 + rank_2 - 2, dim, typename ProductType<Number, OtherNumber>::type>::tensor_type result;
  const Tensor<rank_2, dim, OtherNumber> src1 (src1s);
  return src1*src2;
}



/**
 * Output operator for symmetric tensors of rank 2. Print the elements
 * consecutively, with a space in between, two spaces between rank 1
 * subtensors, three between rank 2 and so on. No special amends are made to
 * represents the symmetry in the output, for example by outputting only the
 * unique entries.
 *
 * @relatesalso SymmetricTensor
 */
template <int dim, typename Number>
inline
std::ostream &operator << (std::ostream &out,
                           const SymmetricTensor<2,dim,Number> &t)
{
  //make our lives a bit simpler by outputting
  //the tensor through the operator for the
  //general Tensor class
  Tensor<2,dim,Number> tt;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      tt[i][j] = t[i][j];

  return out << tt;
}



/**
 * Output operator for symmetric tensors of rank 4. Print the elements
 * consecutively, with a space in between, two spaces between rank 1
 * subtensors, three between rank 2 and so on. No special amends are made to
 * represents the symmetry in the output, for example by outputting only the
 * unique entries.
 *
 * @relatesalso SymmetricTensor
 */
template <int dim, typename Number>
inline
std::ostream &operator << (std::ostream &out,
                           const SymmetricTensor<4,dim,Number> &t)
{
  //make our lives a bit simpler by outputting
  //the tensor through the operator for the
  //general Tensor class
  Tensor<4,dim,Number> tt;

  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
        for (unsigned int l=0; l<dim; ++l)
          tt[i][j][k][l] = t[i][j][k][l];

  return out << tt;
}


DEAL_II_NAMESPACE_CLOSE

#endif
