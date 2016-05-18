// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2016 by the deal.II authors
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

#ifndef dealii__table_h
#define dealii__table_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/table_indices.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/aligned_vector.h>

#include <cstddef>
#include <algorithm>

DEAL_II_NAMESPACE_OPEN

// forward declaration
template <int N, typename T> class TableBase;
template <int N, typename T> class Table;
template <typename T> class Table<1,T>;
template <typename T> class Table<2,T>;
template <typename T> class Table<3,T>;
template <typename T> class Table<4,T>;
template <typename T> class Table<5,T>;
template <typename T> class Table<6,T>;



namespace internal
{

  /**
   * @internal Have a namespace in which we declare some classes that are used
   * to access the elements of tables using the <tt>operator[]</tt>. These are
   * quite technical, since they have to do their work recursively (due to the
   * fact that the number of indices is not known, we have to return an
   * iterator into the next lower dimension object if we access one object,
   * until we are on the lowest level and can actually return a reference to
   * the stored data type itself).  This is so technical that you will not
   * usually want to look at these classes at all, except possibly for
   * educational reasons.  None of the classes herein has a interface that you
   * should use explicitly in your programs (except, of course, through access
   * to the elements of tables with <tt>operator[]</tt>, which generates
   * temporary objects of the types of this namespace).
   *
   * @author Wolfgang Bangerth, 2002
   */
  namespace TableBaseAccessors
  {
    /**
     * @internal Have a class which declares some nested typedefs, depending
     * on its template parameters. The general template declares nothing, but
     * there are more useful specializations regarding the last parameter
     * indicating constness of the table for which accessor objects are to be
     * generated in this namespace.
     */
    template <int N, typename T, bool Constness>
    struct Types
    {};

    /**
     * @internal Have a class which declares some nested typedefs, depending
     * on its template parameters. Specialization for accessors to constant
     * objects.
     */
    template <int N, typename T> struct Types<N,T,true>
    {
      typedef const T value_type;
      typedef const TableBase<N,T> TableType;

      typedef typename AlignedVector<T>::const_iterator iterator;
      typedef typename AlignedVector<T>::const_iterator const_iterator;

      typedef typename AlignedVector<T>::const_reference reference;
      typedef typename AlignedVector<T>::const_reference const_reference;
    };

    /**
     * @internal Have a class which declares some nested typedefs, depending
     * on its template parameters. Specialization for accessors to non-
     * constant objects.
     */
    template <int N, typename T> struct Types<N,T,false>
    {
      typedef T value_type;
      typedef TableBase<N,T> TableType;

      typedef typename AlignedVector<T>::iterator iterator;
      typedef typename AlignedVector<T>::const_iterator const_iterator;

      typedef typename AlignedVector<T>::reference reference;
      typedef typename AlignedVector<T>::const_reference const_reference;
    };


    /**
     * @internal Class that acts as accessor to subobjects of tables of type
     * <tt>Table<N,T></tt>. The template parameter <tt>C</tt> may be either
     * true or false, and indicates whether the objects worked on are constant
     * or not (i.e. write access is only allowed if the value is false).
     *
     * Since with <tt>N</tt> indices, the effect of applying
     * <tt>operator[]</tt> is getting access to something we <tt>N-1</tt>
     * indices, we have to implement these accessor classes recursively, with
     * stopping when we have only one index left. For the latter case, a
     * specialization of this class is declared below, where calling
     * <tt>operator[]</tt> gives you access to the objects actually stored by
     * the table. In the value given to the index operator needs to be checked
     * whether it is inside its bounds, for which we need to know which index
     * of the table we are actually accessing presently. This is done through
     * the template parameter <tt>P</tt>: it indicates, how many remaining
     * indices there are. For a vector, <tt>P</tt> may only be one (and then
     * the specialization below is used). For a table this value may be two,
     * and when using <tt>operator[]</tt>, an object with <tt>P=1</tt>
     * emerges.
     *
     * The value of <tt>P</tt> is also used to determine the stride: this
     * object stores a pointer indicating the beginning of the range of
     * objects that it may access. When we apply <tt>operator[]</tt> on this
     * object, the resulting new accessor may only access a subset of these
     * elements, and to know which subset we need to know the dimensions of
     * the table and the present index, which is indicated by <tt>P</tt>.
     *
     * As stated for the entire namespace, you will not usually have to do
     * with these classes directly, and should not try to use their interface
     * directly as it may change without notice. In fact, since the
     * constructors are made private, you will not even be able to generate
     * objects of this class, as they are only thought as temporaries for
     * access to elements of the table class, not for passing them around as
     * arguments of functions, etc.
     *
     * @author Wolfgang Bangerth, 2002
     */
    template <int N, typename T, bool C, unsigned int P>
    class Accessor
    {
    public:
      typedef typename Types<N,T,C>::TableType TableType;

      typedef typename Types<N,T,C>::iterator iterator;
      typedef typename Types<N,T,C>::const_iterator const_iterator;

      typedef size_t size_type;
      typedef ptrdiff_t difference_type;
    private:
      /**
       * Constructor. Take a pointer to the table object to know about the
       * sizes of the various dimensions, and a pointer to the subset of data
       * we may access.
       */
      Accessor (const TableType &table,
                const iterator    data);

      /**
       * Default constructor. Not needed, and invisible, so private.
       */
      Accessor ();

    public:

      /**
       * Copy constructor. This constructor is public so that one can pass
       * sub-tables to functions as arguments, as in <code>f(table[i])</code>.
       *
       * Using this constructor is risky if accessors are stored longer than
       * the table it points to. Don't do this.
       */
      Accessor (const Accessor &a);

      /**
       * Index operator. Performs a range check.
       */
      Accessor<N,T,C,P-1> operator [] (const unsigned int i) const;

      /**
       * Exception for range check. Do not use global exception since this way
       * we can output which index is the wrong one.
       */
      DeclException3 (ExcIndexRange, int, int, int,
                      << "Index " << N-P+1 << "has a value of "
                      << arg1 << " but needs to be in the range ["
                      << arg2 << "," << arg3 << "[.");
    private:
      /**
       * Store the data given to the constructor. There are no non-const
       * member functions of this class, so there is no reason not to make
       * these elements constant.
       */
      const TableType &table;
      const iterator   data;

      // declare some other classes
      // as friends. make sure to
      // work around bugs in some
      // compilers
      template <int N1, typename T1> friend class dealii::Table;
      template <int N1, typename T1, bool C1, unsigned int P1>
      friend class Accessor;
#  ifndef DEAL_II_TEMPL_SPEC_FRIEND_BUG
      friend class dealii::Table<N,T>;
      friend class Accessor<N,T,C,P+1>;
#  endif
    };



    /**
     * @internal Accessor class for tables. This is the specialization for the
     * last index, which actually allows access to the elements of the table,
     * rather than recursively returning access objects for further subsets.
     * The same holds for this specialization as for the general template; see
     * there for more information.
     *
     * @author Wolfgang Bangerth, 2002
     */
    template <int N, typename T, bool C>
    class Accessor<N,T,C,1>
    {
    public:
      /**
       * Typedef constant and non-constant iterator types to the elements of
       * this row, as well as all the other types usually required for the
       * standard library algorithms.
       */
      typedef typename Types<N,T,C>::value_type value_type;

      typedef typename Types<N,T,C>::iterator iterator;
      typedef typename Types<N,T,C>::const_iterator const_iterator;

      typedef typename Types<N,T,C>::reference reference;
      typedef typename Types<N,T,C>::const_reference const_reference;

      typedef size_t size_type;
      typedef ptrdiff_t difference_type;

      /**
       * Import a typedef from the switch class above.
       */
      typedef typename Types<N,T,C>::TableType    TableType;

    private:

      /**
       * Constructor. Take an iterator to the table object to know about the
       * sizes of the various dimensions, and a iterator to the subset of data
       * we may access (which in this particular case is only one row).
       *
       * The constructor is made private in order to prevent you having such
       * objects around. The only way to create such objects is via the
       * <tt>Table</tt> class, which only generates them as temporary objects.
       * This guarantees that the accessor objects go out of scope earlier
       * than the mother object, avoid problems with data consistency.
       */
      Accessor (const TableType &table,
                const iterator    data);

      /**
       * Default constructor. Not needed, so private.
       */
      Accessor ();

    public:
      /**
       * Copy constructor. This constructor is public so that one can pass
       * sub-tables to functions as arguments, as in <code>f(table[i])</code>.
       *
       * Using this constructor is risky if accessors are stored longer than
       * the table it points to. Don't do this.
       */
      Accessor (const Accessor &a);


      /**
       * Index operator. Performs a range check.
       */
      reference operator [] (const unsigned int) const;

      /**
       * Return the length of one row, i.e. the number of elements
       * corresponding to the last index of the table object.
       */
      unsigned int size () const;

      /**
       * Return an iterator to the first element of this row.
       */
      iterator begin () const;

      /**
       * Return an iterator to the element past the end of this row.
       */
      iterator end () const;

    private:
      /**
       * Store the data given to the constructor. There are no non-const
       * member functions of this class, so there is no reason not to make
       * these elements constant.
       */
      const TableType &table;
      const iterator   data;

      // declare some other classes
      // as friends. make sure to
      // work around bugs in some
      // compilers
      template <int N1, typename T1> friend class dealii::Table;
      template <int N1, typename T1, bool C1, unsigned int P1>
      friend class Accessor;
#  ifndef DEAL_II_TEMPL_SPEC_FRIEND_BUG
      friend class dealii::Table<2,T>;
      friend class Accessor<N,T,C,2>;
#  endif
    };
  }

} // namespace internal




/**
 * General class holding an array of objects of templated type in multiple
 * dimensions. If the template parameter indicating the number of dimensions
 * is one, then this is more or less a vector, if it is two then it is a
 * matrix, and so on.
 *
 * Previously, this data type was emulated in this library by constructs like
 * <tt>std::vector<std::vector<T>></tt>, or even higher nested constructs.
 * However, this has the disadvantage that it is hard to initialize, and most
 * importantly that it is very inefficient if all rows have the same size
 * (which is the usual case), since then the memory for each row is allocated
 * independently, both wasting time and memory. This can be made more
 * efficient by allocating only one chunk of memory for the entire object.
 *
 * Therefore, this data type was invented. Its implementation is rather
 * straightforward, with two exceptions. The first thing to think about is how
 * to pass the size in each of the coordinate directions to the object; this
 * is done using the TableIndices class. Second, how to access the individual
 * elements. The basic problem here is that we would like to make the number
 * of arguments to be passed to the constructor as well as the access
 * functions dependent on the template parameter <tt>N</tt> indicating the
 * number of dimensions. Of course, this is not possible.
 *
 * The way out of the first problem (and partly the second one as well) is to
 * have a common base class TableBase and a derived class for each value of
 * <tt>N</tt>.  This derived class has a constructor with the correct number
 * of arguments, namely <tt>N</tt>. These then transform their arguments into
 * the data type the base class (this class in fact) uses in the constructor
 * as well as in element access through operator() functions.
 *
 * The second problem is that we would like to allow access through a sequence
 * of <tt>operator[]</tt> calls. This mostly because, as said, this class is a
 * replacement for previous use of nested <tt>std::vector</tt> objects, where
 * we had to use the <tt>operator[]</tt> access function recursively until we
 * were at the innermost object. Emulating this behavior without losing the
 * ability to do index checks, and in particular without losing performance is
 * possible but nontrivial, and done in the TableBaseAccessors namespace.
 *
 *
 * <h3>Comparison with the Tensor class</h3>
 *
 * In some way, this class is similar to the Tensor class, in that it
 * templatizes on the number of dimensions. However, there are two major
 * differences. The first is that the Tensor class stores only numeric values
 * (as <tt>double</tt>s), while the Table class stores arbitrary objects. The
 * second is that the Tensor class has fixed dimensions, also given as a
 * template argument, while this class can handle arbitrary dimensions, which
 * may also be different between different indices.
 *
 * This has two consequences. First, since the size is not known at compile
 * time, it has to do explicit memory allocating. Second, the layout of
 * individual elements is not known at compile time, so access is slower than
 * for the Tensor class where the number of elements are their location is
 * known at compile time and the compiler can optimize with this knowledge
 * (for example when unrolling loops). On the other hand, this class is of
 * course more flexible, for example when you want a two-dimensional table
 * with the number of rows equal to the number of degrees of freedom on a
 * cell, and the number of columns equal to the number of quadrature points.
 * Both numbers may only be known at run-time, so a flexible table is needed
 * here. Furthermore, you may want to store, say, the gradients of shape
 * functions, so the data type is not a single scalar value, but a tensor
 * itself.
 *
 * @ingroup data
 * @author Wolfgang Bangerth, 2002.
 */
template <int N, typename T>
class TableBase : public Subscriptor
{
public:
  typedef T value_type;

  /**
   * Integer type used to count the number of elements in this container.
   */
  typedef typename AlignedVector<T>::size_type size_type;


  /**
   * Default constructor. Set all dimensions to zero.
   */
  TableBase ();

  /**
   * Constructor. Initialize the array with the given dimensions in each index
   * component.
   */
  TableBase (const TableIndices<N> &sizes);

  /**
   * Constructor. Initialize the array with the given dimensions in each index
   * component, and then initialize the elements of the table using the second
   * and third argument by calling fill(entries,C_style_indexing).
   */
  template <typename InputIterator>
  TableBase (const TableIndices<N> &sizes,
             InputIterator entries,
             const bool      C_style_indexing = true);

  /**
   * Copy constructor. Performs a deep copy.
   */
  TableBase (const TableBase<N,T> &src);

  /**
   * Copy constructor. Performs a deep copy from a table object storing some
   * other data type.
   */
  template <typename T2>
  TableBase (const TableBase<N,T2> &src);

#ifdef DEAL_II_WITH_CXX11
  /**
   * Move constructor. Transfers the contents of another Table.
   *
   * @note This constructor is only available if deal.II is built with C++11.
   */
  TableBase (TableBase<N,T> &&src);
#endif

  /**
   * Destructor. Free allocated memory.
   */
  ~TableBase ();

  /**
   * Assignment operator. Copy all elements of <tt>src</tt> into the matrix.
   * The size is adjusted if needed.
   *
   * We can't use the other, templatized version since if we don't declare
   * this one, the compiler will happily generate a predefined copy operator
   * which is not what we want.
   */
  TableBase<N,T> &operator = (const TableBase<N,T> &src);

  /**
   * Copy operator. Copy all elements of <tt>src</tt> into the array. The size
   * is adjusted if needed.
   *
   * This function requires that the type <tt>T2</tt> is convertible to
   * <tt>T</tt>.
   */
  template<typename T2>
  TableBase<N,T> &operator = (const TableBase<N,T2> &src);

#ifdef DEAL_II_WITH_CXX11
  /**
   * Move assignment operator. Transfer all elements of <tt>src</tt> into the
   * table.
   *
   * @note This operator is only available if deal.II is built with C++11.
   */
  TableBase<N,T> &operator = (TableBase<N,T> &&src);
#endif

  /**
   * Test for equality of two tables.
   */
  bool operator == (const TableBase<N,T> &T2)  const;

  /**
   * Set all entries to their default value (i.e. copy them over with default
   * constructed objects). Do not change the size of the table, though.
   */
  void reset_values ();

  /**
   * Set the dimensions of this object to the sizes given in the argument, and
   * newly allocate the required memory. If
   * <tt>omit_default_initialization</tt> is set to <tt>false</tt>, all
   * elements of the table are set to a default constructed object for the
   * element type. Otherwise the memory is left in an uninitialized or
   * otherwise undefined state.
   */
  void reinit (const TableIndices<N> &new_size,
               const bool             omit_default_initialization = false);

  /**
   * Size of the table in direction <tt>i</tt>.
   */
  unsigned int size (const unsigned int i) const;

  /**
   * Return the sizes of this object in each direction.
   */
  const TableIndices<N> &size () const;

  /**
   * Return the number of elements stored in this object, which is the product
   * of the extensions in each dimension.
   */
  size_type n_elements () const;

  /**
   * Return whether the object is empty, i.e. one of the directions is zero.
   * This is equivalent to <tt>n_elements()==0</tt>.
   */
  bool empty () const;

  /**
   * Fill this table (which is assumed to already have the correct size) from
   * a source given by dereferencing the given forward iterator (which could,
   * for example, be a pointer to the first element of an array, or an
   * inserting std::istream_iterator). The second argument denotes whether the
   * elements pointed to are arranged in a way that corresponds to the last
   * index running fastest or slowest. The default is to use C-style indexing
   * where the last index runs fastest (as opposed to Fortran-style where the
   * first index runs fastest when traversing multidimensional arrays. For
   * example, if you try to fill an object of type Table<2,T>, then calling
   * this function with the default value for the second argument will result
   * in the equivalent of doing
   * @code
   *   Table<2,T> t;
   *   for (unsigned int i=0; i<t.sizes()[0]; ++i)
   *     for (unsigned int j=0; j<t.sizes()[1]; ++j)
   *       t[i][j] = *entries++;
   * @endcode
   * On the other hand, if the second argument to this function is false, then
   * this would result in code of the following form:
   * @code
   *   Table<2,T> t;
   *   for (unsigned int j=0; j<t.sizes()[1]; ++j)
   *     for (unsigned int i=0; i<t.sizes()[0]; ++i)
   *       t[i][j] = *entries++;
   * @endcode
   * Note the switched order in which we fill the table elements by traversing
   * the given set of iterators.
   *
   * @param entries An iterator to a set of elements from which to initialize
   * this table. It is assumed that iterator can be incremented and
   * dereferenced a sufficient number of times to fill this table.
   * @param C_style_indexing If true, run over elements of the table with the
   * last index changing fastest as we dereference subsequent elements of the
   * input range. If false, change the first index fastest.
   */
  template <typename InputIterator>
  void fill (InputIterator entries,
             const bool      C_style_indexing = true);

  /**
   * Fill all table entries with the same value.
   */
  void fill (const T &value);

  /**
   * Return a read-write reference to the indicated element.
   */
  typename AlignedVector<T>::reference
  operator () (const TableIndices<N> &indices);

  /**
   * Return the value of the indicated element as a read-only reference.
   *
   * We return the requested value as a constant reference rather than by
   * value since this object may hold data types that may be large, and we
   * don't know here whether copying is expensive or not.
   */
  typename AlignedVector<T>::const_reference
  operator () (const TableIndices<N> &indices) const;

  /**
   * Swap the contents of this table and the other table @p v. One could do
   * this operation with a temporary variable and copying over the data
   * elements, but this function is significantly more efficient since it only
   * swaps the pointers to the data of the two vectors and therefore does not
   * need to allocate temporary storage and move data around.
   *
   * This function is analog to the the @p swap function of all C++ standard
   * containers. Also, there is a global function <tt>swap(u,v)</tt> that
   * simply calls <tt>u.swap(v)</tt>, again in analogy to standard functions.
   */
  void swap (TableBase<N,T> &v);

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t memory_consumption () const;

  /**
   * Write or read the data of this object to or from a stream for the purpose
   * of serialization.
   */
  template <class Archive>
  void serialize (Archive &ar, const unsigned int version);

protected:
  /**
   * Return the position of the indicated element within the array of elements
   * stored one after the other. This function does no index checking.
   */
  size_type position (const TableIndices<N> &indices) const;

  /**
   * Return a read-write reference to the indicated element.
   *
   * This function does no bounds checking and is only to be used internally
   * and in functions already checked.
   */
  typename AlignedVector<T>::reference el (const TableIndices<N> &indices);

  /**
   * Return the value of the indicated element as a read-only reference.
   *
   * This function does no bounds checking and is only to be used internally
   * and in functions already checked.
   *
   * We return the requested value as a constant reference rather than by
   * value since this object may hold data types that may be large, and we
   * don't know here whether copying is expensive or not.
   */
  typename AlignedVector<T>::const_reference el (const TableIndices<N> &indices) const;

protected:
  /**
   * Component-array.
   */
  AlignedVector<T> values;

  /**
   * Size in each direction of the table.
   */
  TableIndices<N> table_size;

  /**
   * Make all other table classes friends.
   */
  template <int, typename> friend class TableBase;
};


/**
 * A class representing a table with arbitrary but fixed number of indices.
 * This general template implements some additional functions over those
 * provided by the TableBase class, such as indexing functions taking the
 * correct number of arguments, etc.
 *
 * Rather than this general template, these functions are implemented in
 * partial specializations of this class, with fixed numbers of dimensions.
 * See there, and in the documentation of the base class for more information.
 *
 * @ingroup data
 * @author Wolfgang Bangerth, 2002
 */
template <int N,typename T>
class Table : public TableBase<N,T>
{
};


/**
 * A class representing a one-dimensional table, i.e. a vector-like class.
 * Since the C++ library has a vector class, there is probably not much need
 * for this particular class, but since it is so simple to implement on top of
 * the template base class, we provide it anyway.
 *
 * For the rationale of this class, and a description of the interface, see
 * the base class.
 *
 * @ingroup data
 * @author Wolfgang Bangerth, 2002
 */
template <typename T>
class Table<1,T> : public TableBase<1,T>
{
public:
  /**
   * Integer type used to count the number of elements in this container.
   */
  typedef typename TableBase<1,T>::size_type size_type;

  /**
   * Default constructor. Set all dimensions to zero.
   */
  Table ();

  /**
   * Constructor. Pass down the given dimension to the base class.
   */
  Table (const unsigned int size);

  /**
   * Constructor. Create a table with a given size and initialize it from a
   * set of iterators.
   *
   * This function is entirely equivalent to creating a table <code>t</code>
   * of the given size and then calling
   * @code
   *   t.fill (entries, C_style_indexing);
   * @endcode
   * on it, using the TableBase::fill() function where the arguments are
   * explained in more detail. The point, however, is that that is only
   * possible if the table can be changed after running the constructor,
   * whereas calling the current constructor allows sizing and initializing an
   * object right away so that it can be marked const.
   *
   * Using this constructor, you can do things like this:
   * @code
   *   const double values[] = { 1, 2, 3 };
   *   const Table<1,double> t(3, entries, true);
   * @endcode
   * You can also initialize a table right from a file, using input iterators:
   * @code
   *   std::ifstream input ("myfile");
   *   const Table<1,double> t(3,
   *                           std::istream_iterator<double>(input),
   *                           true);
   * @endcode
   *
   *
   * @param size The size of this one-dimensional table.
   * @param entries An iterator to a set of elements from which to initialize
   * this table. It is assumed that iterator can be incremented and
   * dereferenced a sufficient number of times to fill this table.
   * @param C_style_indexing If true, run over elements of the table with the
   * last index changing fastest as we dereference subsequent elements of the
   * input range. If false, change the first index fastest.
   */
  template <typename InputIterator>
  Table (const unsigned int size,
         InputIterator entries,
         const bool      C_style_indexing = true);

  /**
   * Access operator. Since this is a one-dimensional object, this simply
   * accesses the requested data element. Returns a read-only reference.
   */
  typename AlignedVector<T>::const_reference
  operator [] (const unsigned int i) const;

  /**
   * Access operator. Since this is a one-dimensional object, this simply
   * accesses the requested data element. Returns a read-write reference.
   */
  typename AlignedVector<T>::reference
  operator [] (const unsigned int i);

  /**
   * Access operator. Since this is a one-dimensional object, this simply
   * accesses the requested data element. Returns a read-only reference.
   */
  typename AlignedVector<T>::const_reference
  operator () (const unsigned int i) const;

  /**
   * Access operator. Since this is a one-dimensional object, this simply
   * accesses the requested data element. Returns a read-write reference.
   */
  typename AlignedVector<T>::reference
  operator () (const unsigned int i);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::reference
  operator () (const TableIndices<1> &indices);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::const_reference
  operator () (const TableIndices<1> &indices) const;
};



/**
 * A class representing a two-dimensional table, i.e. a matrix of objects (not
 * necessarily only numbers).
 *
 * For the rationale of this class, and a description of the interface, see
 * the base class. Since this serves as the base class of the full matrix
 * classes in this library, and to keep a minimal compatibility with a
 * predecessor class (<tt>vector2d</tt>), some additional functions are
 * provided.
 *
 * @ingroup data
 * @author Wolfgang Bangerth, 2002
 */
template <typename T>
class Table<2,T> : public TableBase<2,T>
{
public:
  /**
   * Integer type used to count the number of elements in this container.
   */
  typedef typename TableBase<2,T>::size_type size_type;

  /**
   * Default constructor. Set all dimensions to zero.
   */
  Table ();

  /**
   * Constructor. Pass down the given dimensions to the base class.
   */
  Table (const unsigned int size1,
         const unsigned int size2);

  /**
   * Constructor. Create a table with a given size and initialize it from a
   * set of iterators.
   *
   * This function is entirely equivalent to creating a table <code>t</code>
   * of the given size and then calling
   * @code
   *   t.fill (entries, C_style_indexing);
   * @endcode
   * on it, using the TableBase::fill() function where the arguments are
   * explained in more detail. The point, however, is that that is only
   * possible if the table can be changed after running the constructor,
   * whereas calling the current constructor allows sizing and initializing an
   * object right away so that it can be marked const.
   *
   * Using this constructor, you can do things like this:
   * @code
   *   const double values[] = { 1, 2, 3, 4, 5, 6 };
   *   const Table<2,double> t(2, 3, entries, true);
   * @endcode
   * You can also initialize a table right from a file, using input iterators:
   * @code
   *   std::ifstream input ("myfile");
   *   const Table<2,double> t(2, 3,
   *                           std::istream_iterator<double>(input),
   *                           true);
   * @endcode
   *
   *
   * @param size1 The size of this table in the first dimension.
   * @param size2 The size of this table in the second dimension.
   * @param entries An iterator to a set of elements from which to initialize
   * this table. It is assumed that iterator can be incremented and
   * dereferenced a sufficient number of times to fill this table.
   * @param C_style_indexing If true, run over elements of the table with the
   * last index changing fastest as we dereference subsequent elements of the
   * input range. If false, change the first index fastest.
   */
  template <typename InputIterator>
  Table (const unsigned int size1,
         const unsigned int size2,
         InputIterator entries,
         const bool      C_style_indexing = true);

  /**
   * Reinitialize the object. This function is mostly here for compatibility
   * with the earlier <tt>vector2d</tt> class. Passes down to the base class
   * by converting the arguments to the data type requested by the base class.
   */
  void reinit (const unsigned int size1,
               const unsigned int size2,
               const bool         omit_default_initialization = false);

  using TableBase<2,T>::reinit;

  /**
   * Access operator. Generate an object that accesses the requested row of
   * this two-dimensional table. Range checks are performed.
   *
   * This version of the function only allows read access.
   */
  dealii::internal::TableBaseAccessors::Accessor<2,T,true,1>
  operator [] (const unsigned int i) const;

  /**
   * Access operator. Generate an object that accesses the requested row of
   * this two-dimensional table. Range checks are performed.
   *
   * This version of the function allows read-write access.
   */
  dealii::internal::TableBaseAccessors::Accessor<2,T,false,1>
  operator [] (const unsigned int i);

  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function only allows read access.
   */
  typename AlignedVector<T>::const_reference
  operator () (const unsigned int i,
               const unsigned int j) const;


  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function allows read-write access.
   */
  typename AlignedVector<T>::reference
  operator () (const unsigned int i,
               const unsigned int j);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::reference
  operator () (const TableIndices<2> &indices);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::const_reference
  operator () (const TableIndices<2> &indices) const;


  /**
   * Number of rows. This function really makes only sense since we have a
   * two-dimensional object here.
   */
  unsigned int n_rows () const;

  /**
   * Number of columns. This function really makes only sense since we have a
   * two-dimensional object here.
   */
  unsigned int n_cols () const;

protected:
  /**
   * Return a read-write reference to the element <tt>(i,j)</tt>.
   *
   * This function does no bounds checking and is only to be used internally
   * and in functions already checked.
   *
   * These functions are mainly here for compatibility with a former
   * implementation of these table classes for 2d arrays, then called
   * <tt>vector2d</tt>.
   */
  typename AlignedVector<T>::reference el (const unsigned int i,
                                           const unsigned int j);

  /**
   * Return the value of the element <tt>(i,j)</tt> as a read-only reference.
   *
   * This function does no bounds checking and is only to be used internally
   * and in functions already checked.
   *
   * We return the requested value as a constant reference rather than by
   * value since this object may hold data types that may be large, and we
   * don't know here whether copying is expensive or not.
   *
   * These functions are mainly here for compatibility with a former
   * implementation of these table classes for 2d arrays, then called
   * <tt>vector2d</tt>.
   */
  typename AlignedVector<T>::const_reference el (const unsigned int i,
                                                 const unsigned int j) const;
};



/**
 * A class representing a three-dimensional table of objects (not necessarily
 * only numbers).
 *
 * For the rationale of this class, and a description of the interface, see
 * the base class.
 *
 * @ingroup data
 * @author Wolfgang Bangerth, 2002
 */
template <typename T>
class Table<3,T> : public TableBase<3,T>
{
public:
  /**
   * Integer type used to count the number of elements in this container.
   */
  typedef typename TableBase<3,T>::size_type size_type;

  /**
   * Default constructor. Set all dimensions to zero.
   */
  Table ();

  /**
   * Constructor. Pass down the given dimensions to the base class.
   */
  Table (const unsigned int size1,
         const unsigned int size2,
         const unsigned int size3);

  /**
   * Constructor. Create a table with a given size and initialize it from a
   * set of iterators.
   *
   * This function is entirely equivalent to creating a table <code>t</code>
   * of the given size and then calling
   * @code
   *   t.fill (entries, C_style_indexing);
   * @endcode
   * on it, using the TableBase::fill() function where the arguments are
   * explained in more detail. The point, however, is that that is only
   * possible if the table can be changed after running the constructor,
   * whereas calling the current constructor allows sizing and initializing an
   * object right away so that it can be marked const.
   *
   * Using this constructor, you can do things like this (shown here for a
   * two-dimensional table, but the same works for the current class):
   * @code
   *   const double values[] = { 1, 2, 3, 4, 5, 6 };
   *   const Table<2,double> t(2, 3, entries, true);
   * @endcode
   * You can also initialize a table right from a file, using input iterators:
   * @code
   *   std::ifstream input ("myfile");
   *   const Table<2,double> t(2, 3,
   *                           std::istream_iterator<double>(input),
   *                           true);
   * @endcode
   *
   *
   * @param size1 The size of this table in the first dimension.
   * @param size2 The size of this table in the second dimension.
   * @param size3 The size of this table in the third dimension.
   * @param entries An iterator to a set of elements from which to initialize
   * this table. It is assumed that iterator can be incremented and
   * dereferenced a sufficient number of times to fill this table.
   * @param C_style_indexing If true, run over elements of the table with the
   * last index changing fastest as we dereference subsequent elements of the
   * input range. If false, change the first index fastest.
   */
  template <typename InputIterator>
  Table (const unsigned int size1,
         const unsigned int size2,
         const unsigned int size3,
         InputIterator entries,
         const bool      C_style_indexing = true);

  /**
   * Access operator. Generate an object that accesses the requested two-
   * dimensional subobject of this three-dimensional table. Range checks are
   * performed.
   *
   * This version of the function only allows read access.
   */
  dealii::internal::TableBaseAccessors::Accessor<3,T,true,2>
  operator [] (const unsigned int i) const;

  /**
   * Access operator. Generate an object that accesses the requested two-
   * dimensional subobject of this three-dimensional table. Range checks are
   * performed.
   *
   * This version of the function allows read-write access.
   */
  dealii::internal::TableBaseAccessors::Accessor<3,T,false,2>
  operator [] (const unsigned int i);

  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function only allows read access.
   */
  typename AlignedVector<T>::const_reference operator () (const unsigned int i,
                                                          const unsigned int j,
                                                          const unsigned int k) const;


  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function allows read-write access.
   */
  typename AlignedVector<T>::reference operator () (const unsigned int i,
                                                    const unsigned int j,
                                                    const unsigned int k);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::reference operator () (const TableIndices<3> &indices);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::const_reference operator () (const TableIndices<3> &indices) const;
};



/**
 * A class representing a four-dimensional table of objects (not necessarily
 * only numbers).
 *
 * For the rationale of this class, and a description of the interface, see
 * the base class.
 *
 * @ingroup data
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <typename T>
class Table<4,T> : public TableBase<4,T>
{
public:
  /**
   * Integer type used to count the number of elements in this container.
   */
  typedef typename TableBase<4,T>::size_type size_type;

  /**
   * Default constructor. Set all dimensions to zero.
   */
  Table ();

  /**
   * Constructor. Pass down the given dimensions to the base class.
   */
  Table (const unsigned int size1,
         const unsigned int size2,
         const unsigned int size3,
         const unsigned int size4);

  /**
   * Access operator. Generate an object that accesses the requested three-
   * dimensional subobject of this four-dimensional table. Range checks are
   * performed.
   *
   * This version of the function only allows read access.
   */
  dealii::internal::TableBaseAccessors::Accessor<4,T,true,3>
  operator [] (const unsigned int i) const;

  /**
   * Access operator. Generate an object that accesses the requested three-
   * dimensional subobject of this four-dimensional table. Range checks are
   * performed.
   *
   * This version of the function allows read-write access.
   */
  dealii::internal::TableBaseAccessors::Accessor<4,T,false,3>
  operator [] (const unsigned int i);

  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function only allows read access.
   */
  typename AlignedVector<T>::const_reference operator () (const unsigned int i,
                                                          const unsigned int j,
                                                          const unsigned int k,
                                                          const unsigned int l) const;


  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function allows read-write access.
   */
  typename AlignedVector<T>::reference operator () (const unsigned int i,
                                                    const unsigned int j,
                                                    const unsigned int k,
                                                    const unsigned int l);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::reference
  operator () (const TableIndices<4> &indices);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::const_reference
  operator () (const TableIndices<4> &indices) const;
};



/**
 * A class representing a five-dimensional table of objects (not necessarily
 * only numbers).
 *
 * For the rationale of this class, and a description of the interface, see
 * the base class.
 *
 * @ingroup data
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <typename T>
class Table<5,T> : public TableBase<5,T>
{
public:
  /**
   * Integer type used to count the number of elements in this container.
   */
  typedef typename TableBase<5,T>::size_type size_type;


  /**
   * Default constructor. Set all dimensions to zero.
   */
  Table ();

  /**
   * Constructor. Pass down the given dimensions to the base class.
   */
  Table (const unsigned int size1,
         const unsigned int size2,
         const unsigned int size3,
         const unsigned int size4,
         const unsigned int size5);

  /**
   * Access operator. Generate an object that accesses the requested four-
   * dimensional subobject of this five-dimensional table. Range checks are
   * performed.
   *
   * This version of the function only allows read access.
   */
  dealii::internal::TableBaseAccessors::Accessor<5,T,true,4>
  operator [] (const unsigned int i) const;

  /**
   * Access operator. Generate an object that accesses the requested four-
   * dimensional subobject of this five-dimensional table. Range checks are
   * performed.
   *
   * This version of the function allows read-write access.
   */
  dealii::internal::TableBaseAccessors::Accessor<5,T,false,4>
  operator [] (const unsigned int i);

  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function only allows read access.
   */
  typename AlignedVector<T>::const_reference operator () (const unsigned int i,
                                                          const unsigned int j,
                                                          const unsigned int k,
                                                          const unsigned int l,
                                                          const unsigned int m) const;

  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function allows read-write access.
   */
  typename AlignedVector<T>::reference operator () (const unsigned int i,
                                                    const unsigned int j,
                                                    const unsigned int k,
                                                    const unsigned int l,
                                                    const unsigned int m);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::reference
  operator () (const TableIndices<5> &indices);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::const_reference
  operator () (const TableIndices<5> &indices) const;
};



/**
 * A class representing a six-dimensional table of objects (not necessarily
 * only numbers).
 *
 * For the rationale of this class, and a description of the interface, see
 * the base class.
 *
 * @ingroup data
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <typename T>
class Table<6,T> : public TableBase<6,T>
{
public:
  /**
   * Integer type used to count the number of elements in this container.
   */
  typedef typename TableBase<6,T>::size_type size_type;

  /**
   * Default constructor. Set all dimensions to zero.
   */
  Table ();

  /**
   * Constructor. Pass down the given dimensions to the base class.
   */
  Table (const unsigned int size1,
         const unsigned int size2,
         const unsigned int size3,
         const unsigned int size4,
         const unsigned int size5,
         const unsigned int size6);

  /**
   * Access operator. Generate an object that accesses the requested five-
   * dimensional subobject of this six-dimensional table. Range checks are
   * performed.
   *
   * This version of the function only allows read access.
   */
  dealii::internal::TableBaseAccessors::Accessor<6,T,true,5>
  operator [] (const unsigned int i) const;

  /**
   * Access operator. Generate an object that accesses the requested five-
   * dimensional subobject of this six-dimensional table. Range checks are
   * performed.
   *
   * This version of the function allows read-write access.
   */
  dealii::internal::TableBaseAccessors::Accessor<6,T,false,5>
  operator [] (const unsigned int i);

  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function only allows read access.
   */
  typename AlignedVector<T>::const_reference operator () (const unsigned int i,
                                                          const unsigned int j,
                                                          const unsigned int k,
                                                          const unsigned int l,
                                                          const unsigned int m,
                                                          const unsigned int n) const;

  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function allows read-write access.
   */
  typename AlignedVector<T>::reference operator () (const unsigned int i,
                                                    const unsigned int j,
                                                    const unsigned int k,
                                                    const unsigned int l,
                                                    const unsigned int m,
                                                    const unsigned int n);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::reference
  operator () (const TableIndices<6> &indices);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::const_reference
  operator () (const TableIndices<6> &indices) const;
};


/**
 * A class representing a seven-dimensional table of objects (not necessarily
 * only numbers).
 *
 * For the rationale of this class, and a description of the interface, see
 * the base class.
 *
 * @ingroup data
 * @author Wolfgang Bangerth, 2002, Ralf Hartmann 2004
 */
template <typename T>
class Table<7,T> : public TableBase<7,T>
{
public:
  /**
   * Integer type used to count the number of elements in this container.
   */
  typedef typename TableBase<7,T>::size_type size_type;

  /**
   * Default constructor. Set all dimensions to zero.
   */
  Table ();

  /**
   * Constructor. Pass down the given dimensions to the base class.
   */
  Table (const unsigned int size1,
         const unsigned int size2,
         const unsigned int size3,
         const unsigned int size4,
         const unsigned int size5,
         const unsigned int size6,
         const unsigned int size7);

  /**
   * Access operator. Generate an object that accesses the requested six-
   * dimensional subobject of this seven-dimensional table. Range checks are
   * performed.
   *
   * This version of the function only allows read access.
   */
  dealii::internal::TableBaseAccessors::Accessor<7,T,true,6>
  operator [] (const unsigned int i) const;

  /**
   * Access operator. Generate an object that accesses the requested six-
   * dimensional subobject of this seven-dimensional table. Range checks are
   * performed.
   *
   * This version of the function allows read-write access.
   */
  dealii::internal::TableBaseAccessors::Accessor<7,T,false,6>
  operator [] (const unsigned int i);

  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function only allows read access.
   */
  typename AlignedVector<T>::const_reference operator () (const unsigned int i,
                                                          const unsigned int j,
                                                          const unsigned int k,
                                                          const unsigned int l,
                                                          const unsigned int m,
                                                          const unsigned int n,
                                                          const unsigned int o) const;

  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function allows read-write access.
   */
  typename AlignedVector<T>::reference operator () (const unsigned int i,
                                                    const unsigned int j,
                                                    const unsigned int k,
                                                    const unsigned int l,
                                                    const unsigned int m,
                                                    const unsigned int n,
                                                    const unsigned int o);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::reference
  operator () (const TableIndices<7> &indices);

  /**
   * Make the corresponding operator () from the TableBase base class
   * available also in this class.
   */
  typename AlignedVector<T>::const_reference
  operator () (const TableIndices<7> &indices) const;
};



/**
 * A class representing a transpose two-dimensional table, i.e. a matrix of
 * objects (not necessarily only numbers) in column first numbering (FORTRAN
 * convention). The only real difference is therefore really in the storage
 * format.
 *
 * This class copies the functions of Table<2,T>, but the element access and
 * the dimensions will be for the transpose ordering of the data field in
 * TableBase.
 *
 * @ingroup data
 * @author Guido Kanschat, 2005
 */
template <typename T>
class TransposeTable : public TableBase<2,T>
{
public:
  /**
   * Integer type used to count the number of elements in this container.
   */
  typedef typename TableBase<2,T>::size_type size_type;

  /**
   * Default constructor. Set all dimensions to zero.
   */
  TransposeTable ();

  /**
   * Constructor. Pass down the given dimensions to the base class.
   */
  TransposeTable (const unsigned int size1,
                  const unsigned int size2);

  /**
   * Reinitialize the object. This function is mostly here for compatibility
   * with the earlier <tt>vector2d</tt> class. Passes down to the base class
   * by converting the arguments to the data type requested by the base class.
   */
  void reinit (const unsigned int size1,
               const unsigned int size2,
               const bool         omit_default_initialization = false);

  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function only allows read access.
   */
  typename AlignedVector<T>::const_reference operator () (const unsigned int i,
                                                          const unsigned int j) const;

  /**
   * Direct access to one element of the table by specifying all indices at
   * the same time. Range checks are performed.
   *
   * This version of the function allows read-write access.
   */
  typename AlignedVector<T>::reference operator () (const unsigned int i,
                                                    const unsigned int j);

  /**
   * Number of rows. This function really makes only sense since we have a
   * two-dimensional object here.
   */
  unsigned int n_rows () const;

  /**
   * Number of columns. This function really makes only sense since we have a
   * two-dimensional object here.
   */
  unsigned int n_cols () const;

protected:
  /**
   * Return a read-write reference to the element <tt>(i,j)</tt>.
   *
   * This function does no bounds checking and is only to be used internally
   * and in functions already checked.
   *
   * These functions are mainly here for compatibility with a former
   * implementation of these table classes for 2d arrays, then called
   * <tt>vector2d</tt>.
   */
  typename AlignedVector<T>::reference el (const unsigned int i,
                                           const unsigned int j);

  /**
   * Return the value of the element <tt>(i,j)</tt> as a read-only reference.
   *
   * This function does no bounds checking and is only to be used internally
   * and in functions already checked.
   *
   * We return the requested value as a constant reference rather than by
   * value since this object may hold data types that may be large, and we
   * don't know here whether copying is expensive or not.
   *
   * These functions are mainly here for compatibility with a former
   * implementation of these table classes for 2d arrays, then called
   * <tt>vector2d</tt>.
   */
  typename AlignedVector<T>::const_reference el (const unsigned int i,
                                                 const unsigned int j) const;
};




/* --------------------- Template and inline functions ---------------- */

#ifndef DOXYGEN

template <int N, typename T>
TableBase<N,T>::TableBase ()
{}



template <int N, typename T>
TableBase<N,T>::TableBase (const TableIndices<N> &sizes)
{
  reinit (sizes);
}



template <int N, typename T>
template <typename InputIterator>
TableBase<N,T>::
TableBase (const TableIndices<N> &sizes,
           InputIterator entries,
           const bool      C_style_indexing)
{
  reinit (sizes);
  fill (entries, C_style_indexing);
}




template <int N, typename T>
TableBase<N,T>::TableBase (const TableBase<N,T> &src)
  :
  Subscriptor ()
{
  reinit (src.table_size, true);
  values = src.values;
}



template <int N, typename T>
template <typename T2>
TableBase<N,T>::TableBase (const TableBase<N,T2> &src)
{
  reinit (src.table_size);
  if (src.n_elements() != 0)
    std::copy (src.values.begin(), src.values.end(), values.begin());
}



#ifdef DEAL_II_WITH_CXX11

template <int N, typename T>
TableBase<N,T>::TableBase (TableBase<N,T> &&src)
  :
  Subscriptor (std::move(src)),
  values (std::move(src.values)),
  table_size (src.table_size)
{
  src.table_size = TableIndices<N>();
}

#endif



template <int N, typename T>
template <class Archive>
inline
void
TableBase<N,T>::serialize (Archive &ar, const unsigned int)
{
  ar &static_cast<Subscriptor &>(*this);

  ar &values &table_size;
}



namespace internal
{
  namespace TableBaseAccessors
  {
    template <int N, typename T, bool C, unsigned int P>
    inline
    Accessor<N,T,C,P>::Accessor (const TableType &table,
                                 const iterator    data)
      :
      table (table),
      data (data)
    {}



    template <int N, typename T, bool C, unsigned int P>
    inline
    Accessor<N,T,C,P>::Accessor (const Accessor &a)
      :
      table (a.table),
      data (a.data)
    {}



    template <int N, typename T, bool C, unsigned int P>
    inline
    Accessor<N,T,C,P>::Accessor ()
      :
      table (*static_cast<const TableType *>(0)),
      data (0)
    {
      // accessor objects are only
      // temporary objects, so should
      // not need to be copied around
      Assert (false, ExcInternalError());
    }



    template <int N, typename T, bool C, unsigned int P>
    inline
    Accessor<N,T,C,P-1>
    Accessor<N,T,C,P>::operator [] (const unsigned int i) const
    {
      Assert (i < table.size()[N-P],
              ExcIndexRange (i, 0, table.size()[N-P]));

      // access i-th
      // subobject. optimize on the
      // case i==0
      if (i==0)
        return Accessor<N,T,C,P-1> (table, data);
      else
        {
          // note: P>1, otherwise the
          // specialization would have
          // been taken!
          size_type subobject_size = table.size()[N-1];
          for (int p=P-1; p>1; --p)
            subobject_size *= table.size()[N-p];
          const iterator new_data = data + i*subobject_size;
          return Accessor<N,T,C,P-1> (table, new_data);
        }
    }



    template <int N, typename T, bool C>
    inline
    Accessor<N,T,C,1>::Accessor (const TableType &table,
                                 const iterator    data)
      :
      table (table),
      data (data)
    {}



    template <int N, typename T, bool C>
    inline
    Accessor<N,T,C,1>::Accessor ()
      :
      table (*static_cast<const TableType *>(0)),
      data (0)
    {
      // accessor objects are only
      // temporary objects, so should
      // not need to be copied around
      Assert (false, ExcInternalError());
    }



    template <int N, typename T, bool C>
    inline
    Accessor<N,T,C,1>::Accessor (const Accessor &a)
      :
      table (a.table),
      data (a.data)
    {}



    template <int N, typename T, bool C>
    inline
    typename Accessor<N,T,C,1>::reference
    Accessor<N,T,C,1>::operator [] (const unsigned int i) const
    {
      Assert (i < table.size()[N-1],
              ExcIndexRange (i, 0, table.size()[N-1]));
      return *(data+i);
    }



    template <int N, typename T, bool C>
    inline
    unsigned int
    Accessor<N,T,C,1>::size () const
    {
      return table.size()[N-1];
    }



    template <int N, typename T, bool C>
    inline
    typename Accessor<N,T,C,1>::iterator
    Accessor<N,T,C,1>::begin () const
    {
      return data;
    }



    template <int N, typename T, bool C>
    inline
    typename Accessor<N,T,C,1>::iterator
    Accessor<N,T,C,1>::end () const
    {
      return data+table.size()[N-1];
    }
  }
}



template <int N, typename T>
inline
TableBase<N,T>::~TableBase ()
{}



template <int N, typename T>
inline
TableBase<N,T> &
TableBase<N,T>::operator = (const TableBase<N,T> &m)
{
  if (!m.empty())
    values = m.values;
  reinit (m.size(), true);

  return *this;
}



template <int N, typename T>
template <typename T2>
inline
TableBase<N,T> &
TableBase<N,T>::operator = (const TableBase<N,T2> &m)
{
  reinit (m.size(), true);
  if (!empty())
    std::copy (m.values.begin(), m.values.begin() + n_elements(),
               values.begin());

  return *this;
}



#ifdef DEAL_II_WITH_CXX11

template <int N, typename T>
inline
TableBase<N,T> &
TableBase<N,T>::operator = (TableBase<N,T> &&m)
{
  static_cast<Subscriptor &>(*this) = std::move(m);
  values = std::move(m.values);
  table_size = m.table_size;
  m.table_size = TableIndices<N>();

  return *this;
}

#endif



template <int N, typename T>
inline
bool
TableBase<N,T>::operator == (const TableBase<N,T> &T2)  const
{
  return (values == T2.values);
}



template <int N, typename T>
inline
void
TableBase<N,T>::reset_values ()
{
  // use parallel set operation
  if (n_elements() != 0)
    values.fill(T());
}



template <int N, typename T>
inline
void
TableBase<N,T>::fill (const T &value)
{
  if (n_elements() != 0)
    values.fill(value);
}




template <int N, typename T>
inline
void
TableBase<N,T>::reinit (const TableIndices<N> &new_sizes,
                        const bool             omit_default_initialization)
{
  table_size = new_sizes;

  const size_type new_size = n_elements();

  // if zero size was given: free all memory
  if (new_size == 0)
    {
      values.resize (0);
      // set all sizes to zero, even
      // if one was previously
      // nonzero. This simplifies
      // some assertions.
      table_size = TableIndices<N>();

      return;
    }

  // adjust values field. If it was empty before, we can simply call resize(),
  // which can set all the data fields. Otherwise, select the fast resize and
  // manually fill in all the elements as required by the design of this
  // class. (Selecting another code for the empty case ensures that we touch
  // the memory only once for non-trivial classes that need to initialize the
  // memory also in resize_fast.)
  if (!omit_default_initialization)
    {
      if (values.empty())
        values.resize(new_size, T());
      else
        {
          values.resize_fast(new_size);
          values.fill(T());
        }
    }
  else
    values.resize_fast (new_size);
}



template <int N, typename T>
inline
const TableIndices<N> &
TableBase<N,T>::size () const
{
  return table_size;
}



template <int N, typename T>
inline
unsigned int
TableBase<N,T>::size (const unsigned int i) const
{
  Assert (i<N, ExcIndexRange(i,0,N));
  return table_size[i];
}



template <int N, typename T>
inline
typename TableBase<N,T>::size_type
TableBase<N,T>::n_elements () const
{
  size_type s = 1;
  for (unsigned int n=0; n<N; ++n)
    s *= table_size[n];
  return s;
}



template <int N, typename T>
inline
bool
TableBase<N,T>::empty () const
{
  return (n_elements() == 0);
}



namespace internal
{
  namespace Table
  {
    template <typename InputIterator, typename T>
    void fill_Fortran_style (InputIterator  entries,
                             TableBase<1,T>  &table)
    {
      for (unsigned int i=0; i<table.size()[0]; ++i)
        table(TableIndices<1>(i)) = *entries++;
    }


    template <typename InputIterator, typename T>
    void fill_Fortran_style (InputIterator  entries,
                             TableBase<2,T>  &table)
    {
      for (unsigned int j=0; j<table.size()[1]; ++j)
        for (unsigned int i=0; i<table.size()[0]; ++i)
          table(TableIndices<2>(i,j)) = *entries++;
    }


    template <typename InputIterator, typename T>
    void fill_Fortran_style (InputIterator  entries,
                             TableBase<3,T>  &table)
    {
      for (unsigned int k=0; k<table.size()[2]; ++k)
        for (unsigned int j=0; j<table.size()[1]; ++j)
          for (unsigned int i=0; i<table.size()[0]; ++i)
            table(TableIndices<3>(i,j,k)) = *entries++;
    }


    template <typename InputIterator, typename T, int N>
    void fill_Fortran_style (InputIterator,
                             TableBase<N,T> &)
    {
      Assert (false, ExcNotImplemented());
    }
  }
}


template <int N, typename T>
template <typename InputIterator>
inline
void
TableBase<N,T>::fill (InputIterator entries,
                      const bool C_style_indexing)
{
  Assert (n_elements() != 0,
          ExcMessage("Trying to fill an empty matrix."));

  if (C_style_indexing)
    for (typename AlignedVector<T>::iterator p = values.begin();
         p != values.end(); ++p)
      *p = *entries++;
  else
    internal::Table::fill_Fortran_style (entries, *this);
}



template <int N, typename T>
inline
void
TableBase<N,T>::swap (TableBase<N,T> &v)
{
  values.swap(v.values);
  std::swap (table_size, v.table_size);
}



template <int N, typename T>
inline
std::size_t
TableBase<N,T>::memory_consumption () const
{
  return sizeof(*this) + MemoryConsumption::memory_consumption(values);
}



template <int N, typename T>
inline
typename TableBase<N,T>::size_type
TableBase<N,T>::position (const TableIndices<N> &indices) const
{
  // specialize this for the
  // different numbers of dimensions,
  // to make the job somewhat easier
  // for the compiler. have the
  // general formula nevertheless:
  switch (N)
    {
    case 1:
      return indices[0];
    case 2:
      return size_type(indices[0])*table_size[1] + indices[1];
    case 3:
      return ((size_type(indices[0])*table_size[1] + indices[1])*table_size[2]
              + indices[2]);
    default:
    {
      size_type s = indices[0];
      for (unsigned int n=1; n<N; ++n)
        s = s*table_size[n] + indices[n];
      return s;
    }
    }
}



template <int N, typename T>
inline
typename AlignedVector<T>::const_reference
TableBase<N,T>::operator () (const TableIndices<N> &indices) const
{
  for (unsigned int n=0; n<N; ++n)
    Assert (indices[n] < table_size[n],
            ExcIndexRange (indices[n], 0, table_size[n]));
  return el(indices);
}



template <int N, typename T>
inline
typename AlignedVector<T>::reference
TableBase<N,T>::operator () (const TableIndices<N> &indices)
{
  for (unsigned int n=0; n<N; ++n)
    Assert (indices[n] < table_size[n],
            ExcIndexRange (indices[n], 0, table_size[n]));
  return el(indices);
}



template <int N, typename T>
inline
typename AlignedVector<T>::const_reference
TableBase<N,T>::el (const TableIndices<N> &indices) const
{
  return values[position(indices)];
}



template <int N, typename T>
inline
typename AlignedVector<T>::reference
TableBase<N,T>::el (const TableIndices<N> &indices)
{
  Assert (position(indices) < values.size(),
          ExcIndexRange (position(indices), 0, values.size()));
  return values[position(indices)];
}



template <typename T>
inline
Table<1,T>::Table ()
{}



template <typename T>
inline
Table<1,T>::Table (const unsigned int size)
  :
  TableBase<1,T> (TableIndices<1> (size))
{}



template <typename T>
template <typename InputIterator>
inline
Table<1,T>::Table (const unsigned int size,
                   InputIterator entries,
                   const bool C_style_indexing)
  :
  TableBase<1,T> (TableIndices<1> (size),
                  entries,
                  C_style_indexing)
{}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<1,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return this->values[i];
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<1,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return this->values[i];
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<1,T>::operator () (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return this->values[i];
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<1,T>::operator () (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return this->values[i];
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<1,T>::operator () (const TableIndices<1> &indices) const
{
  return TableBase<1,T>::operator () (indices);
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<1,T>::operator () (const TableIndices<1> &indices)
{
  return TableBase<1,T>::operator () (indices);
}


//---------------------------------------------------------------------------

template <typename T>
inline
Table<2,T>::Table ()
{}



template <typename T>
inline
Table<2,T>::Table (const unsigned int size1,
                   const unsigned int size2)
  :
  TableBase<2,T> (TableIndices<2> (size1, size2))
{}



template <typename T>
template <typename InputIterator>
inline
Table<2,T>::Table (const unsigned int size1,
                   const unsigned int size2,
                   InputIterator entries,
                   const bool C_style_indexing)
  :
  TableBase<2,T> (TableIndices<2> (size1, size2),
                  entries,
                  C_style_indexing)
{}



template <typename T>
inline
void
Table<2,T>::reinit (const unsigned int size1,
                    const unsigned int size2,
                    const bool         omit_default_initialization)
{
  this->TableBase<2,T>::reinit (TableIndices<2> (size1, size2),omit_default_initialization);
}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<2,T,true,1>
Table<2,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return dealii::internal::TableBaseAccessors::Accessor<2,T,true,1>(*this,
         this->values.begin()+size_type(i)*n_cols());
}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<2,T,false,1>
Table<2,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return dealii::internal::TableBaseAccessors::Accessor<2,T,false,1>(*this,
         this->values.begin()+size_type(i)*n_cols());
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<2,T>::operator () (const unsigned int i,
                         const unsigned int j) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  return this->values[size_type(i)*this->table_size[1]+j];
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<2,T>::operator () (const unsigned int i,
                         const unsigned int j)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  return this->values[size_type(i)*this->table_size[1]+j];
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<2,T>::operator () (const TableIndices<2> &indices) const
{
  return TableBase<2,T>::operator () (indices);
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<2,T>::operator () (const TableIndices<2> &indices)
{
  return TableBase<2,T>::operator () (indices);
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<2,T>::el (const unsigned int i,
                const unsigned int j) const
{
  return this->values[size_type(i)*this->table_size[1]+j];
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<2,T>::el (const unsigned int i,
                const unsigned int j)
{
  return this->values[size_type(i)*this->table_size[1]+j];
}



template <typename T>
inline
unsigned int
Table<2,T>::n_rows () const
{
  return this->table_size[0];
}



template <typename T>
inline
unsigned int
Table<2,T>::n_cols () const
{
  return this->table_size[1];
}



//---------------------------------------------------------------------------

template <typename T>
inline
TransposeTable<T>::TransposeTable ()
{}



template <typename T>
inline
TransposeTable<T>::TransposeTable (const unsigned int size1,
                                   const unsigned int size2)
  :
  TableBase<2,T> (TableIndices<2> (size2, size1))
{}



template <typename T>
inline
void
TransposeTable<T>::reinit (const unsigned int size1,
                           const unsigned int size2,
                           const bool         omit_default_initialization)
{
  this->TableBase<2,T>::reinit (TableIndices<2> (size2, size1), omit_default_initialization);
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
TransposeTable<T>::operator () (const unsigned int i,
                                const unsigned int j) const
{
  Assert (i < this->table_size[1],
          ExcIndexRange (i, 0, this->table_size[1]));
  Assert (j < this->table_size[0],
          ExcIndexRange (j, 0, this->table_size[0]));
  return this->values[size_type(j)*this->table_size[1]+i];
}



template <typename T>
inline
typename AlignedVector<T>::reference
TransposeTable<T>::operator () (const unsigned int i,
                                const unsigned int j)
{
  Assert (i < this->table_size[1],
          ExcIndexRange (i, 0, this->table_size[1]));
  Assert (j < this->table_size[0],
          ExcIndexRange (j, 0, this->table_size[0]));
  return this->values[size_type(j)*this->table_size[1]+i];
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
TransposeTable<T>::el (const unsigned int i,
                       const unsigned int j) const
{
  return this->values[size_type(j)*this->table_size[1]+i];
}



template <typename T>
inline
typename AlignedVector<T>::reference
TransposeTable<T>::el (const unsigned int i,
                       const unsigned int j)
{
  return this->values[size_type(j)*this->table_size[1]+i];
}



template <typename T>
inline
unsigned int
TransposeTable<T>::n_rows () const
{
  return this->table_size[1];
}



template <typename T>
inline
unsigned int
TransposeTable<T>::n_cols () const
{
  return this->table_size[0];
}



//---------------------------------------------------------------------------


template <typename T>
inline
Table<3,T>::Table ()
{}



template <typename T>
inline
Table<3,T>::Table (const unsigned int size1,
                   const unsigned int size2,
                   const unsigned int size3)
  :
  TableBase<3,T> (TableIndices<3> (size1, size2, size3))
{}



template <typename T>
template <typename InputIterator>
inline
Table<3,T>::Table (const unsigned int size1,
                   const unsigned int size2,
                   const unsigned int size3,
                   InputIterator entries,
                   const bool C_style_indexing)
  :
  TableBase<3,T> (TableIndices<3> (size1, size2, size3),
                  entries,
                  C_style_indexing)
{}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<3,T,true,2>
Table<3,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2];
  return (dealii::internal::TableBaseAccessors::Accessor<3,T,true,2>
          (*this,
           this->values.begin() + i*subobject_size));
}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<3,T,false,2>
Table<3,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2];
  return (dealii::internal::TableBaseAccessors::Accessor<3,T,false,2>
          (*this,
           this->values.begin() + i*subobject_size));
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<3,T>::operator () (const unsigned int i,
                         const unsigned int j,
                         const unsigned int k) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  Assert (k < this->table_size[2],
          ExcIndexRange (k, 0, this->table_size[2]));
  return this->values[(size_type(i)*this->table_size[1]+j)
                      *this->table_size[2] + k];
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<3,T>::operator () (const unsigned int i,
                         const unsigned int j,
                         const unsigned int k)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  Assert (k < this->table_size[2],
          ExcIndexRange (k, 0, this->table_size[2]));
  return this->values[(size_type(i)*this->table_size[1]+j)
                      *this->table_size[2] + k];
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<3,T>::operator () (const TableIndices<3> &indices) const
{
  return TableBase<3,T>::operator () (indices);
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<3,T>::operator () (const TableIndices<3> &indices)
{
  return TableBase<3,T>::operator () (indices);
}



template <typename T>
inline
Table<4,T>::Table ()
{}



template <typename T>
inline
Table<4,T>::Table (const unsigned int size1,
                   const unsigned int size2,
                   const unsigned int size3,
                   const unsigned int size4)
  :
  TableBase<4,T> (TableIndices<4> (size1, size2, size3, size4))
{}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<4,T,true,3>
Table<4,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] *
                                   this->table_size[3];
  return (dealii::internal::TableBaseAccessors::Accessor<4,T,true,3>
          (*this,
           this->values.begin() + i*subobject_size));
}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<4,T,false,3>
Table<4,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] *
                                   this->table_size[3];
  return (dealii::internal::TableBaseAccessors::Accessor<4,T,false,3>
          (*this,
           this->values.begin() + i*subobject_size));
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<4,T>::operator () (const unsigned int i,
                         const unsigned int j,
                         const unsigned int k,
                         const unsigned int l) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  Assert (k < this->table_size[2],
          ExcIndexRange (k, 0, this->table_size[2]));
  Assert (l < this->table_size[3],
          ExcIndexRange (l, 0, this->table_size[3]));
  return this->values[((size_type(i)*this->table_size[1]+j)
                       *this->table_size[2] + k)
                      *this->table_size[3] + l];
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<4,T>::operator () (const unsigned int i,
                         const unsigned int j,
                         const unsigned int k,
                         const unsigned int l)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  Assert (k < this->table_size[2],
          ExcIndexRange (k, 0, this->table_size[2]));
  Assert (l < this->table_size[3],
          ExcIndexRange (l, 0, this->table_size[3]));
  return this->values[((size_type(i)*this->table_size[1]+j)
                       *this->table_size[2] + k)
                      *this->table_size[3] + l];
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<4,T>::operator () (const TableIndices<4> &indices) const
{
  return TableBase<4,T>::operator () (indices);
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<4,T>::operator () (const TableIndices<4> &indices)
{
  return TableBase<4,T>::operator () (indices);
}



template <typename T>
inline
Table<5,T>::Table ()
{}



template <typename T>
inline
Table<5,T>::Table (const unsigned int size1,
                   const unsigned int size2,
                   const unsigned int size3,
                   const unsigned int size4,
                   const unsigned int size5)
  :
  TableBase<5,T> (TableIndices<5> (size1, size2, size3, size4, size5))
{}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<5,T,true,4>
Table<5,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] *
                                   this->table_size[3] *
                                   this->table_size[4];
  return (dealii::internal::TableBaseAccessors::Accessor<5,T,true,4>
          (*this,
           this->values.begin() + i*subobject_size));
}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<5,T,false,4>
Table<5,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] *
                                   this->table_size[3] *
                                   this->table_size[4];
  return (dealii::internal::TableBaseAccessors::Accessor<5,T,false,4>
          (*this,
           this->values.begin() + i*subobject_size));
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<5,T>::operator () (const unsigned int i,
                         const unsigned int j,
                         const unsigned int k,
                         const unsigned int l,
                         const unsigned int m) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  Assert (k < this->table_size[2],
          ExcIndexRange (k, 0, this->table_size[2]));
  Assert (l < this->table_size[3],
          ExcIndexRange (l, 0, this->table_size[3]));
  Assert (m < this->table_size[4],
          ExcIndexRange (m, 0, this->table_size[4]));
  return this->values[(((size_type(i)*this->table_size[1]+j)
                        *this->table_size[2] + k)
                       *this->table_size[3] + l)
                      *this->table_size[4] + m];
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<5,T>::operator () (const unsigned int i,
                         const unsigned int j,
                         const unsigned int k,
                         const unsigned int l,
                         const unsigned int m)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  Assert (k < this->table_size[2],
          ExcIndexRange (k, 0, this->table_size[2]));
  Assert (l < this->table_size[3],
          ExcIndexRange (l, 0, this->table_size[3]));
  Assert (m < this->table_size[4],
          ExcIndexRange (m, 0, this->table_size[4]));
  return this->values[(((size_type(i)*this->table_size[1]+j)
                        *this->table_size[2] + k)
                       *this->table_size[3] + l)
                      *this->table_size[4] + m];
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<5,T>::operator () (const TableIndices<5> &indices) const
{
  return TableBase<5,T>::operator () (indices);
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<5,T>::operator () (const TableIndices<5> &indices)
{
  return TableBase<5,T>::operator () (indices);
}



template <typename T>
inline
Table<6,T>::Table ()
{}



template <typename T>
inline
Table<6,T>::Table (const unsigned int size1,
                   const unsigned int size2,
                   const unsigned int size3,
                   const unsigned int size4,
                   const unsigned int size5,
                   const unsigned int size6)
  :
  TableBase<6,T> (TableIndices<6> (size1, size2, size3, size4, size5, size6))
{}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<6,T,true,5>
Table<6,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] *
                                   this->table_size[3] *
                                   this->table_size[4] *
                                   this->table_size[5];
  return (dealii::internal::TableBaseAccessors::Accessor<6,T,true,5>
          (*this,
           this->values.begin() + i*subobject_size));
}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<6,T,false,5>
Table<6,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] *
                                   this->table_size[3] *
                                   this->table_size[4] *
                                   this->table_size[5];
  return (dealii::internal::TableBaseAccessors::Accessor<6,T,false,5>
          (*this,
           this->values.begin() + i*subobject_size));
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<6,T>::operator () (const unsigned int i,
                         const unsigned int j,
                         const unsigned int k,
                         const unsigned int l,
                         const unsigned int m,
                         const unsigned int n) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  Assert (k < this->table_size[2],
          ExcIndexRange (k, 0, this->table_size[2]));
  Assert (l < this->table_size[3],
          ExcIndexRange (l, 0, this->table_size[3]));
  Assert (m < this->table_size[4],
          ExcIndexRange (m, 0, this->table_size[4]));
  Assert (n < this->table_size[5],
          ExcIndexRange (n, 0, this->table_size[5]));
  return this->values[((((size_type(i)*this->table_size[1]+j)
                         *this->table_size[2] + k)
                        *this->table_size[3] + l)
                       *this->table_size[4] + m)
                      *this->table_size[5] + n];
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<6,T>::operator () (const unsigned int i,
                         const unsigned int j,
                         const unsigned int k,
                         const unsigned int l,
                         const unsigned int m,
                         const unsigned int n)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  Assert (k < this->table_size[2],
          ExcIndexRange (k, 0, this->table_size[2]));
  Assert (l < this->table_size[3],
          ExcIndexRange (l, 0, this->table_size[3]));
  Assert (m < this->table_size[4],
          ExcIndexRange (m, 0, this->table_size[4]));
  Assert (n < this->table_size[5],
          ExcIndexRange (n, 0, this->table_size[5]));
  return this->values[((((size_type(i)*this->table_size[1]+j)
                         *this->table_size[2] + k)
                        *this->table_size[3] + l)
                       *this->table_size[4] + m)
                      *this->table_size[5] + n];
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<6,T>::operator () (const TableIndices<6> &indices) const
{
  return TableBase<6,T>::operator () (indices);
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<6,T>::operator () (const TableIndices<6> &indices)
{
  return TableBase<6,T>::operator () (indices);
}



template <typename T>
inline
Table<7,T>::Table ()
{}



template <typename T>
inline
Table<7,T>::Table (const unsigned int size1,
                   const unsigned int size2,
                   const unsigned int size3,
                   const unsigned int size4,
                   const unsigned int size5,
                   const unsigned int size6,
                   const unsigned int size7)
  :
  TableBase<7,T> (TableIndices<7> (size1, size2, size3, size4, size5, size6, size7))
{}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<7,T,true,6>
Table<7,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] *
                                   this->table_size[3] *
                                   this->table_size[4] *
                                   this->table_size[5] *
                                   this->table_size[6];
  return (dealii::internal::TableBaseAccessors::Accessor<7,T,true,6>
          (*this,
           this->values.begin() + i*subobject_size));
}



template <typename T>
inline
dealii::internal::TableBaseAccessors::Accessor<7,T,false,6>
Table<7,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const size_type subobject_size = size_type(this->table_size[1]) *
                                   this->table_size[2] *
                                   this->table_size[3] *
                                   this->table_size[4] *
                                   this->table_size[5] *
                                   this->table_size[6];
  return (dealii::internal::TableBaseAccessors::Accessor<7,T,false,6>
          (*this,
           this->values.begin() + i*subobject_size));
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<7,T>::operator () (const unsigned int i,
                         const unsigned int j,
                         const unsigned int k,
                         const unsigned int l,
                         const unsigned int m,
                         const unsigned int n,
                         const unsigned int o) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  Assert (k < this->table_size[2],
          ExcIndexRange (k, 0, this->table_size[2]));
  Assert (l < this->table_size[3],
          ExcIndexRange (l, 0, this->table_size[3]));
  Assert (m < this->table_size[4],
          ExcIndexRange (m, 0, this->table_size[4]));
  Assert (n < this->table_size[5],
          ExcIndexRange (n, 0, this->table_size[5]));
  Assert (o < this->table_size[6],
          ExcIndexRange (o, 0, this->table_size[6]));
  return this->values[(((((size_type(i)*this->table_size[1]+j)
                          *this->table_size[2] + k)
                         *this->table_size[3] + l)
                        *this->table_size[4] + m)
                       *this->table_size[5] + n)
                      *this->table_size[6] + o];
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<7,T>::operator () (const unsigned int i,
                         const unsigned int j,
                         const unsigned int k,
                         const unsigned int l,
                         const unsigned int m,
                         const unsigned int n,
                         const unsigned int o)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  Assert (k < this->table_size[2],
          ExcIndexRange (k, 0, this->table_size[2]));
  Assert (l < this->table_size[3],
          ExcIndexRange (l, 0, this->table_size[3]));
  Assert (m < this->table_size[4],
          ExcIndexRange (m, 0, this->table_size[4]));
  Assert (n < this->table_size[5],
          ExcIndexRange (n, 0, this->table_size[5]));
  Assert (o < this->table_size[5],
          ExcIndexRange (o, 0, this->table_size[6]));
  return this->values[(((((size_type(i)*this->table_size[1]+j)
                          *this->table_size[2] + k)
                         *this->table_size[3] + l)
                        *this->table_size[4] + m)
                       *this->table_size[5] + n)
                      *this->table_size[6] + o];
}



template <typename T>
inline
typename AlignedVector<T>::const_reference
Table<7,T>::operator () (const TableIndices<7> &indices) const
{
  return TableBase<7,T>::operator () (indices);
}



template <typename T>
inline
typename AlignedVector<T>::reference
Table<7,T>::operator () (const TableIndices<7> &indices)
{
  return TableBase<7,T>::operator () (indices);
}


#endif // DOXYGEN



/**
 * Global function @p swap which overloads the default implementation of the
 * C++ standard library which uses a temporary object. The function simply
 * exchanges the data of the two tables.
 *
 * @author Martin Kronbichler, 2013
 */
template <int N, typename T>
inline
void swap (TableBase<N,T> &u, TableBase<N,T> &v)
{
  u.swap (v);
}

DEAL_II_NAMESPACE_CLOSE

#endif
