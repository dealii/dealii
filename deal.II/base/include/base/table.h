//-----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------
#ifndef __deal2__table_h
#define __deal2__table_h

#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>

#include <cstddef>
#include <algorithm>


// forward declaration
template <int N, typename T> class TableBase;
template <int N, typename T> class Table;
template <typename T> class Table<1,T>;
template <typename T> class Table<2,T>;
template <typename T> class Table<3,T>;
template <typename T> class Table<4,T>;
template <typename T> class Table<5,T>;
template <typename T> class Table<6,T>;




/**
 * Base class for an array of indices of fixed size used for the
 * @ref{TableBase} class. Actually, this class serves a dual purpose,
 * as it not only stores indices into said class, but also the sizes
 * of the table in its various coordinates.
 *
 * @author Wolfgang Bangerth, 2002
 */
template <int N>
class TableIndicesBase
{
  public:
                                     /**
                                      * Access the value of the
                                      * @p{i}th index.
                                      */
    unsigned int operator[] (const unsigned int i) const;
    
  protected:
                                     /**
                                      * Store the indices in an array.
                                      */
    unsigned indices[N];
};


/**
 * Array of indices of fixed size used for the @ref{TableBase}
 * class.
 *
 * This is the general template, and has no implementation. There are
 * a number of specializations that are actually implemented (one for
 * each used value of @p{N}), which only differ in the way they
 * implement their constructors (they take @p{N} arguments, something
 * that cannot be represented by a general template). Actual storage
 * of and access to data is done by the @ref{TableIndicesBase} base
 * class of a specializations.
 *
 * @author Wolfgang Bangerth, 2002
 */
template <int N>
class TableIndices
{};



/**
 * Array of indices of fixed size used for the @ref{TableBase}
 * class.
 *
 * This is the specialization for a one-dimensional table, i.e. a
 * vector. This class only differs in the non-default constructors
 * from the other specializations. Actual storage of and access to
 * data is done by the @ref{TableIndicesBase} base class of a
 * specializations.
 *
 * @author Wolfgang Bangerth, 2002
 */
template <>
class TableIndices<1> : public TableIndicesBase<1>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1);
};



/**
 * Array of indices of fixed size used for the @ref{TableBase}
 * class.
 *
 * This is the specialization for a two-dimensional table. This class
 * only differs in the non-default constructors from the other
 * specializations. Actual storage of and access to data is done by
 * the @ref{TableIndicesBase} base class of a specializations.
 *
 * @author Wolfgang Bangerth, 2002
 */
template <>
class TableIndices<2> : public TableIndicesBase<2>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1,
                  const unsigned int index2);
};



/**
 * Array of indices of fixed size used for the @ref{TableBase}
 * class.
 *
 * This is the specialization for a three-dimensional table. This class
 * only differs in the non-default constructors from the other
 * specializations. Actual storage of and access to data is done by
 * the @ref{TableIndicesBase} base class of a specializations.
 *
 * @author Wolfgang Bangerth, 2002
 */
template <>
class TableIndices<3> : public TableIndicesBase<3>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1,
                  const unsigned int index2,
                  const unsigned int index3);
};


/**
 * Array of indices of fixed size used for the @ref{TableBase}
 * class.
 *
 * This is the specialization for a four-dimensional table. This class
 * only differs in the non-default constructors from the other
 * specializations. Actual storage of and access to data is done by
 * the @ref{TableIndicesBase} base class of a specializations.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <>
class TableIndices<4> : public TableIndicesBase<4>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1,
                  const unsigned int index2,
                  const unsigned int index3,
		  const unsigned int index4);
};


/**
 * Array of indices of fixed size used for the @ref{TableBase}
 * class.
 *
 * This is the specialization for a five-dimensional table. This class
 * only differs in the non-default constructors from the other
 * specializations. Actual storage of and access to data is done by
 * the @ref{TableIndicesBase} base class of a specializations.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <>
class TableIndices<5> : public TableIndicesBase<5>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1,
                  const unsigned int index2,
                  const unsigned int index3,
		  const unsigned int index4,
		  const unsigned int index5);
};


/**
 * Array of indices of fixed size used for the @ref{TableBase}
 * class.
 *
 * This is the specialization for a four-dimensional table. This class
 * only differs in the non-default constructors from the other
 * specializations. Actual storage of and access to data is done by
 * the @ref{TableIndicesBase} base class of a specializations.
 *
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <>
class TableIndices<6> : public TableIndicesBase<6>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * indices to zero.
                                      */
    TableIndices ();

                                     /**
                                      * Constructor. Set indices to
                                      * the given values.
                                      */
    TableIndices (const unsigned int index1,
                  const unsigned int index2,
                  const unsigned int index3,
		  const unsigned int index4,
		  const unsigned int index5,
		  const unsigned int index6);
};


/**
 * Have a namespace in which we declare some classes that are used to
 * access the elements of tables using the @p{operator[]}. These are
 * quite technical, since they have to do their work recursively (due
 * to the fact that the number of indices is not known, we have to
 * return an iterator into the next lower dimension object if we
 * access one object, until we are on the lowest level and can
 * actually return a reference to the stored data type itself).  This
 * is so technical that you will not usually want to look at these
 * classes at all, except possibly for educational reasons.  None of
 * the classes herein has a interface that you should use explicitly
 * in your programs (except, of course, through access to the elements
 * of tables with @p{operator[]}, which generates temporary objects of
 * the types of this namespace).
 *
 * @author Wolfgang Bangerth, 2002
 */
namespace TableBaseAccessors
{
/**
 * Have a class which declares some nested typedefs, depending on its
 * template parameters. The general template declares nothing, but
 * there are more useful specializations regaring the last parameter
 * indicating constness of the table for which accessor objects are to
 * be generated in this namespace.
 */
  template <int N, typename T, bool Constness>
  class Types
  {};

/**
 * Have a class which declares some nested typedefs, depending on its
 * template parameters. Specialization for accessors to constant
 * objects.
 */
  template <int N, typename T> struct Types<N,T,true> 
  {
      typedef const T value_type;
      typedef const TableBase<N,T> TableType;
  };

/**
 * Have a class which declares some nested typedefs, depending on its
 * template parameters. Specialization for accessors to non-constant
 * objects.
 */
  template <int N, typename T> struct Types<N,T,false> 
  {
      typedef T value_type;
      typedef TableBase<N,T> TableType;
  };
  

/**
 * Class that acts as accessor to subobjects of tables of type
 * @p{Table<N,T>}. The template parameter @p{C} may be either true or
 * false, and indicates whether the objects worked on are constant or
 * not (i.e. write access is only allowed if the value is false).
 *
 * Since with @p{N} indices, the effect of applying @p{operator[]} is
 * getting access to something we @p{N-1} indices, we have to
 * implement these accessor classes recursively, with stopping when we
 * have only one index left. For the latter case, a specialization of
 * this class is declared below, where calling @p{operator[]} gives
 * you access to the objects actually stored by the table. In the
 * value given to the index operator needs to be checked whether it is
 * inside its bounds, for which we need to know which index of the
 * table we are actually accessing presently. This is done through the
 * template parameter @p{P}: it indicates, how many remaining indices
 * there are. For a vector, @p{P} may only be one (and then the
 * specialization below is used). For a table this value may be two,
 * and when using @p{operator[]}, an object with @p{P=1} emerges.
 *
 * The value of @p{P} is also used to determine the stride: this
 * object stores a pointer indicating the beginning of the range of
 * objects that it may access. When we apply @p{operator[]} on this
 * object, the resulting new accessor may only access a subset of
 * these elements, and to know which subset we need to know the
 * dimensions of the table and the present index, which is indicated
 * by @p{P}.
 *
 * As stated for the entire namespace, you will not usually have to do
 * with these classes directly, and should not try to use their
 * interface directly as it may change without notice. In fact, since
 * the constructors are made private, you will not even be able to
 * generate objects of this class, as they are only thought as
 * temporaries for access to elements of the table class, not for
 * passing them around as arguments of functions, etc.
 * 
 * @author Wolfgang Bangerth, 2002
 */
  template <int N, typename T, bool C, unsigned int P>
  class Accessor
  {
    public:
                                       /**
                                        * Import two typedefs from the
                                        * switch class above.
                                        */
      typedef typename Types<N,T,C>::value_type * pointer;      
      typedef typename Types<N,T,C>::TableType    TableType;

    private:
                                       /**
                                        * Constructor. Take a pointer
                                        * to the table object to know
                                        * about the sizes of the
                                        * various dimensions, and a
                                        * pointer to the subset of
                                        * data we may access.
                                        *
                                        * The constructor is made
                                        * private in order to prevent
                                        * you having such objects
                                        * around. The only way to
                                        * create such objects is via
                                        * the @p{Table} class, which
                                        * only generates them as
                                        * temporary objects. This
                                        * guarantees that the accessor
                                        * objects go out of scope
                                        * earlier than the mother
                                        * object, avoid problems with
                                        * data consistency.
                                        */
      Accessor (const TableType &table,
                const pointer    data);

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
                                        * Index operator. Performs a
                                        * range check.
                                        */
      Accessor<N,T,C,P-1> operator [] (const unsigned int i) const;

                                       /**
                                        * Exception for range
                                        * check. Do not use global
                                        * exception since this way we
                                        * can output which index is
                                        * the wrong one.
                                        */
      DeclException3 (ExcIndexRange, int, int, int,
                      << "The " << N-P+1 << "th index has a value of "
                      << arg1 << " but needs to be in the range ["
                      << arg2 << "," << arg3 << "[");
    private:
                                       /**
                                        * Store the data given to the
                                        * constructor. There are no
                                        * non-const member functions
                                        * of this class, so there is
                                        * no reason not to make these
                                        * elements constant.
                                        */
      const TableType &table;
      const pointer   data;

                                       // declare some other classes
                                       // as friends. make sure to
                                       // work around bugs in some
                                       // compilers
#ifndef DEAL_II_NAMESP_TEMPL_FRIEND_BUG      
      template <int N1, typename T1> friend class Table;
      template <int N1, typename T1, bool C1, unsigned int P1>
      friend class Accessor;
#  ifndef DEAL_II_TEMPL_SPEC_FRIEND_BUG
      friend class Table<N,T>;
      friend class Accessor<N,T,C,P+1>;
#  endif
#else
      friend class Table<N,T>;
      friend class Accessor<N,T,C,P+1>;
#endif
  };


  
/**
 * Accessor class for tables. This is the specialization for the last
 * index, which actually allows access to the elements of the table,
 * rather than recursively returning access objects for further
 * subsets. The same holds for this specialization as for the general
 * template; see there for more information.
 *
 * @author Wolfgang Bangerth, 2002
 */
  template <int N, typename T, bool C>
  class Accessor<N,T,C,1>
  {
    public:
                                       /**
                                        * Typedef constant and
                                        * non-constant iterator
                                        * types to the elements of
                                        * this row, as well as all
                                        * the other types usually
                                        * required for the standard
                                        * library algorithms.
                                        */
      typedef typename Types<N,T,C>::value_type value_type;
      typedef value_type* pointer;
      typedef const value_type* const_pointer;
      typedef value_type* iterator;
      typedef const value_type* const_iterator;
      typedef value_type& reference;
      typedef const value_type& const_reference;
      typedef size_t size_type;
      typedef ptrdiff_t difference_type;

                                       /**
                                        * Import a typedef from the
                                        * switch class above.
                                        */
      typedef typename Types<N,T,C>::TableType    TableType;

    private:
      
                                       /**
                                        * Constructor. Take a pointer
                                        * to the table object to know
                                        * about the sizes of the
                                        * various dimensions, and a
                                        * pointer to the subset of
                                        * data we may access (which in
                                        * this particular case is only
                                        * one row).
                                        *
                                        * The constructor is made
                                        * private in order to prevent
                                        * you having such objects
                                        * around. The only way to
                                        * create such objects is via
                                        * the @p{Table} class, which
                                        * only generates them as
                                        * temporary objects. This
                                        * guarantees that the accessor
                                        * objects go out of scope
                                        * earlier than the mother
                                        * object, avoid problems with
                                        * data consistency.
                                        */
      Accessor (const TableType &table,
                const pointer    data);

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
                                        * Index operator. Performs a
                                        * range check.
                                        */
      reference operator [] (const unsigned int) const;
      
                                       /**
                                        * Return the length of one row,
                                        * i.e. the number of elements
                                        * corresponding to the last
                                        * index of the table object.
                                        */
      unsigned int size () const;
        
                                       /**
                                        * Return an iterator to the
                                        * first element of this
                                        * row.
                                        */
      iterator begin () const;
      
                                       /**
                                        * Return an interator to the
                                        * element past the end of
                                        * this row.
                                        */
      iterator end () const;
      
                                       /**
                                        * Exception for range
                                        * check. Do not use global
                                        * exception since this way we
                                        * can output which index is
                                        * the wrong one.
                                        */
      DeclException3 (ExcIndexRange, int, int, int,
                      << "The " << N << "th index has a value of "
                      << arg1 << " but needs to be in the range ["
                      << arg2 << "," << arg3 << "[");
    private:
                                       /**
                                        * Store the data given to the
                                        * constructor. There are no
                                        * non-const member functions
                                        * of this class, so there is
                                        * no reason not to make these
                                        * elements constant.
                                        */
      const TableType &table;
      const pointer   data;

                                       // declare some other classes
                                       // as friends. make sure to
                                       // work around bugs in some
                                       // compilers
#ifndef DEAL_II_NAMESP_TEMPL_FRIEND_BUG
      template <int N1, typename T1> friend class Table;
      template <int N1, typename T1, bool C1, unsigned int P1>
      friend class Accessor;
#  ifndef DEAL_II_TEMPL_SPEC_FRIEND_BUG
      friend class Table<2,T>;
      friend class Accessor<N,T,C,2>;
#  endif
#else
      friend class Table<2,T>;
      friend class Accessor<N,T,C,2>;
#endif
  };
};
  





/**
 * General class holding an array of objects of templated type in
 * multiple dimensions. If the template parameter indicating the
 * number of dimensions is one, then this is more or less a vector, if
 * it is two then it is a matrix, and so on.
 *
 * Previously, this data type was emulated in this library by
 * constructs like @p{std::vector<std::vector<T>>}, or even higher
 * nested constructs.  However, this has the disadvantage that it is
 * hard to initialize, and most importantly that it is very
 * inefficient if all rows have the same size (which is the usual
 * case), since then the memory for each row is allocated
 * independently, both wasting time and memory. This can be made more
 * efficient by allocating only one chunk of memory for the entire
 * object.
 *
 * Therefore, this data type was invented. Its implementation is
 * rather straightforward, with two exceptions. The first thing to
 * think about is how to pass the size in each of the coordinate
 * directions to the object; this is done using the @ref{TableIndices}
 * class. Second, how to access the individual elements. The basic
 * problem here is that we would like to make the number of arguments
 * to be passed to the constructor as well as the access functions
 * dependent on the template parameter @p{N} indicating the number of
 * dimensions. Of course, this is not possible.
 *
 * The way out of the first problem (and partly the second one as
 * well) is to have derived class for each value of @p{N} that have a
 * constructor with the right number of arguments, one for each
 * dimension. These then transform their arguments into the data type
 * this class wants to see, both for construction as well as access
 * through the @p{operator()} function.
 *
 * The second problem is that we would like to allow access through a
 * sequence of @p{operator[]} calls. This mostly because, as said,
 * this class is a replacement for previous use of nested
 * @p{std::vector} objects, where we had to use the @p{operator[]}
 * access function recursively until we were at the innermost
 * object. Emulating this behavior without losing the ability to do
 * index checks, and in particular without losing performance is
 * possible but nontrivial, and done in the @ref{TableBaseAccessors}
 * namespace.
 *
 *
 * @sect3{Comparison with the Tensor class}
 *
 * In some way, this class is similar to the @ref{Tensor} class, in
 * that it templatizes on the number of dimensions. However, there are
 * two major differences. The first is that the @ref{Tensor} class
 * stores only numeric values (as @p{double}s), while the @p{Table}
 * class stores arbitrary objects. The second is that the @ref{Tensor}
 * class has fixed dimensions, also give as a template argument, while
 * this class can handle arbitrary dimensions, which may also be
 * different between different indices.
 *
 * This has two consequences. First, since the size is not known at
 * compile time, it has to do explicit memory allocating. Second, the
 * layout of individual elements is not known at compile time, so
 * access is slower than for the @ref{Tensor} class where the number
 * of elements are their location is known at compile time and the
 * compiler can optimize with this knowledge (for example when
 * unrolling loops). On the other hand, this class is of course more
 * flexible, for example when you want a two-dimensional table with
 * the number of rows equal to the number of degrees of freedom on a
 * cell, and the number of columns equal to the number of quadrature
 * points. Both numbers may only be known at run-time, so a flexible
 * table is needed here. Furthermore, you may want to store, say, the
 * gradients of shape functions, so the data type is not a single
 * scalar value, but a tensor itself.
 * 
 * @author Wolfgang Bangerth, 2002.
 */
template <int N, typename T>
class TableBase : public Subscriptor
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * dimensions to zero.
                                      */
    TableBase ();
    
                                     /**
                                      * Constructor. Initialize the
                                      * array with the given
                                      * dimensions in each index
                                      * component.
                                      */
    TableBase (const TableIndices<N> &sizes);

                                     /**
                                      * Copy constructor. Performs a
                                      * deep copy.
                                      */
    TableBase (const TableBase<N,T> &src);

                                     /**
                                      * Copy constructor. Performs a
                                      * deep copy from a table object
                                      * storing some other data type.
                                      */
    template <typename T2>
    TableBase (const TableBase<N,T2> &src);
    
                                     /**
                                      * Destructor. Free allocated memory.
                                      */
    ~TableBase ();
    
                                     /**
                                      * Assignment operator.
                                      * Copy all elements of @p{src}
                                      * into the matrix. The size is
                                      * adjusted if needed.
                                      *
                                      * We can't use the other, templatized
                                      * version since if we don't declare
                                      * this one, the compiler will happily
                                      * generate a predefined copy
                                      * operator which is not what we want.
                                      */
    TableBase<N,T>& operator = (const TableBase<N,T>& src);
    
                                     /**
                                      * Copy operator.
                                      * Copy all elements of @p{src}
                                      * into the array. The size is
                                      * adjusted if needed.
                                      *
                                      * This function requires that the
                                      * type @p{T2} is convertible to
                                      * @p{T}.
                                      */
    template<typename T2>
    TableBase<N,T>& operator = (const TableBase<N,T2> &src);
    
                                     /**
                                      * Set all entries to their
                                      * default value (i.e. copy them
                                      * over with default constructed
                                      * objects).
                                      */
    void clear ();
    
                                     /**
                                      * Set the dimensions of this
                                      * object to the sizes given in
                                      * the argument, and newly
                                      * allocate the required
                                      * memory. Forget the previous
                                      * content of the array.
                                      */
    void reinit (const TableIndices<N> &new_size);

                                     /**
                                      * Return the sizes of this
                                      * object in each direction.
                                      */
    const TableIndices<N> & size () const;

                                     /**
                                      * Return the number of elements
                                      * stored in this object, which
                                      * is the product of the
                                      * extensions in each dimension.
                                      */
    unsigned int n_elements () const;

                                     /**
                                      * Return whether the object is
                                      * empty, i.e. one of the
                                      * directions is zero. This is
                                      * equivalent to
                                      * @p{n_elements()==0}.
                                      */
    bool empty () const;
    
                                     /**
                                      * Fill array with an array of
                                      * elements. The input array must
                                      * be arranged in usual C style,
                                      * i.e. with the last index
                                      * running fastest. For
                                      * two-dimensional tables, this
                                      * means line by line. No range
                                      * checking is performed, i.e.,
                                      * it is assumed that the input
                                      * array @p{entries} contains
                                      * @p{n_rows()*n_cols()}
                                      * elements, and that the layout
                                      * refers to the desired shape of
                                      * this table. The only check we
                                      * do is that the present array
                                      * is non-empty.
                                      *
                                      * Note also that the type of the
                                      * objects of the input array,
                                      * @p{T2}, must be convertible to
                                      * the type of the objects of
                                      * this array.
                                      */
    template<typename T2>
    void fill (const T2 *entries);
    
                                     /**
                                      * Return a read-write reference
                                      * to the indicated element.
                                      */
    T & operator() (const TableIndices<N> &indices);
  
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
    const T & operator() (const TableIndices<N> &indices) const;

                                     /**
                                      * Determine an estimate for the
                                      * memory consumption (in bytes)
                                      * of this object.
                                      */
    unsigned int memory_consumption () const;

  protected:
                                     /**
                                      * Return the position of the
                                      * indicated element within the
                                      * array of elements stored one
                                      * after the other. This function
                                      * does no index checking.
                                      */
    unsigned int position (const TableIndices<N> &indices) const;
    
                                     /**
                                      * Return a read-write reference
                                      * to the indicated element.
                                      *
                                      * This function does no bounds
                                      * checking and is only to be
                                      * used internally and in
                                      * functions already checked.
                                      */
    T & el (const TableIndices<N> &indices);
  
                                     /**
                                      * Return the value of the
                                      * indicated element as a
                                      * read-only reference.
                                      *
                                      * This function does no bounds
                                      * checking and is only to be
                                      * used internally and in
                                      * functions already checked.
                                      *
                                      * We return the requested value
                                      * as a constant reference rather
                                      * than by value since this
                                      * object may hold data types
                                      * that may be large, and we
                                      * don't know here whether
                                      * copying is expensive or not.
                                      */
    const T & el (const TableIndices<N> &indices) const;    
  
                                     /**
                                      * Direct read-only access to
                                      * data field. Used by
                                      * @ref{FullMatrix} (there even
                                      * with a cast from const),
                                      * otherwise, keep away!
                                      */
    const T* data () const;
    
  protected:
                                     /**
                                      * Component-array.
                                      */
    T* val;
    
                                     /**
                                      * Size of array. This may be
                                      * larger than the number of
                                      * actually used elements, since
                                      * we don't shrink the array upon
                                      * calls to @p{reinit} unless the
                                      * new size is zero.
                                      */
    unsigned int val_size;    

                                     /**
                                      * Size in each direction of the
                                      * table.
                                      */
    TableIndices<N> table_size;
};


/**
 * A class representing a table with arbitrary but fixed number of
 * indices. This general template implements some additional functions
 * over those provided bythe @ref{TableBase} class, such as indexing
 * functions taking the correct number of arguments, etc.
 *
 * Rather than this general template, these functions are implemented
 * in partial specializations of this class, with fixed numbers of
 * dimensions. See there, and in the documentation of the base class
 * for more information.
 * 
 * @author Wolfgang Bangerth, 2002
 */
template <int N,typename T>
class Table : public TableBase<N,T>
{};


/**
 * A class representing a one-dimensional table, i.e. a vector-like
 * class. Since the C++ library has a vector class, there is probably
 * not much need for this particular class, but since it is so simple
 * to implement on top of the template base class, we provide it
 * anyway.
 *
 * For the rationale of this class, and a description of the
 * interface, see the base class.
 * 
 * @author Wolfgang Bangerth, 2002
 */
template <typename T>
class Table<1,T> : public TableBase<1,T>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * dimensions to zero.
                                      */
    Table ();
    
                                     /**
                                      * Constructor. Pass down the
                                      * given dimension to the base
                                      * class.
                                      */
    Table (const unsigned int size);

                                     /**
                                      * Access operator. Since this is
                                      * a one-dimensional object, this
                                      * simply accesses the requested
                                      * data element. Returns a
                                      * read-only reference.
                                      */
    const T &
    operator [] (const unsigned int i) const;

                                     /**
                                      * Access operator. Since this is
                                      * a one-dimensional object, this
                                      * simply accesses the requested
                                      * data element. Returns a
                                      * read-write reference.
                                      */
    T &
    operator [] (const unsigned int i);

                                     /**
                                      * Access operator. Since this is
                                      * a one-dimensional object, this
                                      * simply accesses the requested
                                      * data element. Returns a
                                      * read-only reference.
                                      */
    const T &
    operator () (const unsigned int i) const;

                                     /**
                                      * Access operator. Since this is
                                      * a one-dimensional object, this
                                      * simply accesses the requested
                                      * data element. Returns a
                                      * read-write reference.
                                      */
    T &
    operator () (const unsigned int i);
};



/**
 * A class representing a two-dimensional table, i.e. a matrix of
 * objects (not necessarily only numbers).
 *
 * For the rationale of this class, and a description of the
 * interface, see the base class. Since this serves as the base class
 * of the full matrix classes in this library, and to keep a minimal
 * compatibility with a predecessor class (@p{vector2d}), some
 * additional functions are provided.
 * 
 * @author Wolfgang Bangerth, 2002
 */
template <typename T>
class Table<2,T> : public TableBase<2,T>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * dimensions to zero.
                                      */
    Table ();

                                     /**
                                      * Constructor. Pass down the
                                      * given dimensions to the base
                                      * class.
                                      */
    Table (const unsigned int size1,
           const unsigned int size2);

                                     /**
                                      * Reinitialize the object. This
                                      * function is mostly here for
                                      * compatibility with the earlier
                                      * @p{vector2d} class. Passes
                                      * down to the base class by
                                      * converting the arguments to
                                      * the data type requested by the
                                      * base class.
                                      */
    void reinit (const unsigned int size1,
                 const unsigned int size2);

                                     /**
                                      * Access operator. Generate an
                                      * object that accesses the
                                      * requested row of this
                                      * two-dimensional table. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * only allows read access.
                                      */
    TableBaseAccessors::Accessor<2,T,true,1>
    operator [] (const unsigned int i) const;

                                     /**
                                      * Access operator. Generate an
                                      * object that accesses the
                                      * requested row of this
                                      * two-dimensional table. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * allows read-write access.
                                      */
    TableBaseAccessors::Accessor<2,T,false,1>
    operator [] (const unsigned int i);

                                     /**
                                      * Direct access to one element
                                      * of the table by specifying all
                                      * indices at the same time. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * only allows read access.
                                      */
    const T & operator () (const unsigned int i,
                           const unsigned int j) const;
    

                                     /**
                                      * Direct access to one element
                                      * of the table by specifying all
                                      * indices at the same time. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * allows read-write access.
                                      */
    T & operator () (const unsigned int i,
                     const unsigned int j);

    
                                     /**
                                      * Number of rows. This function
                                      * really makes only sense since
                                      * we have a two-dimensional
                                      * object here.
                                      */
    unsigned int n_rows () const;
    
                                     /**
                                      * Number of columns. This function
                                      * really makes only sense since
                                      * we have a two-dimensional
                                      * object here.
                                      */
    unsigned int n_cols () const;

  protected:
                                     /**
                                      * Return a read-write reference
                                      * to the element @p{(i,j)}.
                                      *
                                      * This function does no bounds
                                      * checking and is only to be
                                      * used internally and in
                                      * functions already checked.
                                      *
                                      * These functions are mainly
                                      * here for compatibility with a
                                      * former implementation of these
                                      * table classes for 2d arrays,
                                      * then called @p{vector2d}.
                                      */
    T & el (const unsigned int i,
            const unsigned int j);
  
                                     /**
                                      * Return the value of the
                                      * element @p{(i,j)} as a
                                      * read-only reference.
                                      *
                                      * This function does no bounds
                                      * checking and is only to be
                                      * used internally and in
                                      * functions already checked.
                                      *
                                      * We return the requested value
                                      * as a constant reference rather
                                      * than by value since this
                                      * object may hold data types
                                      * that may be large, and we
                                      * don't know here whether
                                      * copying is expensive or not.
                                      *
                                      * These functions are mainly
                                      * here for compatibility with a
                                      * former implementation of these
                                      * table classes for 2d arrays,
                                      * then called @p{vector2d}.
                                      */
    const T & el (const unsigned int i,
                  const unsigned int j) const;
};



/**
 * A class representing a three-dimensional table of objects (not
 * necessarily only numbers).
 *
 * For the rationale of this class, and a description of the
 * interface, see the base class.
 * 
 * @author Wolfgang Bangerth, 2002
 */
template <typename T>
class Table<3,T> : public TableBase<3,T>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * dimensions to zero.
                                      */
    Table ();

                                     /**
                                      * Constructor. Pass down the
                                      * given dimensions to the base
                                      * class.
                                      */
    Table (const unsigned int size1,
           const unsigned int size2,
           const unsigned int size3);

                                     /**
                                      * Access operator. Generate an
                                      * object that accesses the
                                      * requested two-dimensional
                                      * subobject of this
                                      * three-dimensional table. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * only allows read access.
                                      */
    TableBaseAccessors::Accessor<3,T,true,2>
    operator [] (const unsigned int i) const;

                                     /**
                                      * Access operator. Generate an
                                      * object that accesses the
                                      * requested two-dimensional
                                      * subobject of this
                                      * three-dimensional table. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * allows read-write access.
                                      */
    TableBaseAccessors::Accessor<3,T,false,2>
    operator [] (const unsigned int i);

                                     /**
                                      * Direct access to one element
                                      * of the table by specifying all
                                      * indices at the same time. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * only allows read access.
                                      */
    const T & operator () (const unsigned int i,
                           const unsigned int j,
                           const unsigned int k) const;
    

                                     /**
                                      * Direct access to one element
                                      * of the table by specifying all
                                      * indices at the same time. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * allows read-write access.
                                      */
    T & operator () (const unsigned int i,
                     const unsigned int j,
                     const unsigned int k);
};



/**
 * A class representing a four-dimensional table of objects (not
 * necessarily only numbers).
 *
 * For the rationale of this class, and a description of the
 * interface, see the base class.
 * 
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <typename T>
class Table<4,T> : public TableBase<4,T>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * dimensions to zero.
                                      */
    Table ();

                                     /**
                                      * Constructor. Pass down the
                                      * given dimensions to the base
                                      * class.
                                      */
    Table (const unsigned int size1,
           const unsigned int size2,
           const unsigned int size3,
	   const unsigned int size4);
    
                                     /**
                                      * Access operator. Generate an
                                      * object that accesses the
                                      * requested two-dimensional
                                      * subobject of this
                                      * three-dimensional table. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * only allows read access.
                                      */
    TableBaseAccessors::Accessor<4,T,true,3>
    operator [] (const unsigned int i) const;

                                     /**
                                      * Access operator. Generate an
                                      * object that accesses the
                                      * requested two-dimensional
                                      * subobject of this
                                      * three-dimensional table. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * allows read-write access.
                                      */
    TableBaseAccessors::Accessor<4,T,false,3>
    operator [] (const unsigned int i);

                                     /**
                                      * Direct access to one element
                                      * of the table by specifying all
                                      * indices at the same time. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * only allows read access.
                                      */
    const T & operator () (const unsigned int i,
                           const unsigned int j,
                           const unsigned int k,
			   const unsigned int l) const;
    

                                     /**
                                      * Direct access to one element
                                      * of the table by specifying all
                                      * indices at the same time. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * allows read-write access.
                                      */
    T & operator () (const unsigned int i,
                     const unsigned int j,
                     const unsigned int k,
		     const unsigned int l);
};



/**
 * A class representing a five-dimensional table of objects (not
 * necessarily only numbers).
 *
 * For the rationale of this class, and a description of the
 * interface, see the base class.
 * 
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <typename T>
class Table<5,T> : public TableBase<5,T>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * dimensions to zero.
                                      */
    Table ();

                                     /**
                                      * Constructor. Pass down the
                                      * given dimensions to the base
                                      * class.
                                      */
    Table (const unsigned int size1,
           const unsigned int size2,
           const unsigned int size3,
	   const unsigned int size4,
	   const unsigned int size5);

                                     /**
                                      * Access operator. Generate an
                                      * object that accesses the
                                      * requested two-dimensional
                                      * subobject of this
                                      * three-dimensional table. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * only allows read access.
                                      */
    TableBaseAccessors::Accessor<5,T,true,4>
    operator [] (const unsigned int i) const;

                                     /**
                                      * Access operator. Generate an
                                      * object that accesses the
                                      * requested two-dimensional
                                      * subobject of this
                                      * three-dimensional table. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * allows read-write access.
                                      */
    TableBaseAccessors::Accessor<5,T,false,4>
    operator [] (const unsigned int i);

                                     /**
                                      * Direct access to one element
                                      * of the table by specifying all
                                      * indices at the same time. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * only allows read access.
                                      */
    const T & operator () (const unsigned int i,
                           const unsigned int j,
                           const unsigned int k,
			   const unsigned int l,
			   const unsigned int m) const;
    

                                     /**
                                      * Direct access to one element
                                      * of the table by specifying all
                                      * indices at the same time. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * allows read-write access.
                                      */
    T & operator () (const unsigned int i,
                     const unsigned int j,
                     const unsigned int k,
		     const unsigned int l,
		     const unsigned int m);
};



/**
 * A class representing a six-dimensional table of objects (not
 * necessarily only numbers).
 *
 * For the rationale of this class, and a description of the
 * interface, see the base class.
 * 
 * @author Wolfgang Bangerth, Ralf Hartmann 2002
 */
template <typename T>
class Table<6,T> : public TableBase<6,T>
{
  public:
                                     /**
                                      * Default constructor. Set all
                                      * dimensions to zero.
                                      */
    Table ();

                                     /**
                                      * Constructor. Pass down the
                                      * given dimensions to the base
                                      * class.
                                      */
    Table (const unsigned int size1,
           const unsigned int size2,
           const unsigned int size3,
	   const unsigned int size4,
	   const unsigned int size5,
	   const unsigned int size6);

                                     /**
                                      * Access operator. Generate an
                                      * object that accesses the
                                      * requested two-dimensional
                                      * subobject of this
                                      * three-dimensional table. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * only allows read access.
                                      */
    TableBaseAccessors::Accessor<6,T,true,5>
    operator [] (const unsigned int i) const;

                                     /**
                                      * Access operator. Generate an
                                      * object that accesses the
                                      * requested two-dimensional
                                      * subobject of this
                                      * three-dimensional table. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * allows read-write access.
                                      */
    TableBaseAccessors::Accessor<6,T,false,5>
    operator [] (const unsigned int i);

                                     /**
                                      * Direct access to one element
                                      * of the table by specifying all
                                      * indices at the same time. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * only allows read access.
                                      */
    const T & operator () (const unsigned int i,
                           const unsigned int j,
                           const unsigned int k,
			   const unsigned int l,
			   const unsigned int m,
			   const unsigned int n) const;
    

                                     /**
                                      * Direct access to one element
                                      * of the table by specifying all
                                      * indices at the same time. Range
                                      * checks are performed.
                                      *
                                      * This version of the function
                                      * allows read-write access.
                                      */
    T & operator () (const unsigned int i,
                     const unsigned int j,
                     const unsigned int k,
		     const unsigned int l,
		     const unsigned int m,
		     const unsigned int n);
};





/* --------------------- Template and inline functions ---------------- */


template <int N>
inline
unsigned int
TableIndicesBase<N>::operator [] (const unsigned int i) const 
{
  Assert (i < N, ExcIndexRange (i, 0, N));
  return indices[i];
};



inline
TableIndices<1>::TableIndices () 
{
  this->indices[0] = 0;
};



inline
TableIndices<1>::TableIndices (const unsigned int index1) 
{
  this->indices[0] = index1;
};



inline
TableIndices<2>::TableIndices () 
{
  this->indices[0] = this->indices[1] = 0;
};



inline
TableIndices<2>::TableIndices (const unsigned int index1,
                               const unsigned int index2)
{
  this->indices[0] = index1;
  this->indices[1] = index2;
};



inline
TableIndices<3>::TableIndices () 
{
  this->indices[0] = this->indices[1] = this->indices[2] = 0;
};



inline
TableIndices<3>::TableIndices (const unsigned int index1,
                               const unsigned int index2,
                               const unsigned int index3)
{
  this->indices[0] = index1;
  this->indices[1] = index2;
  this->indices[2] = index3;
};


inline
TableIndices<4>::TableIndices () 
{
  this->indices[0] = this->indices[1] = this->indices[2] = this->indices[3] = 0;
};



inline
TableIndices<4>::TableIndices (const unsigned int index1,
                               const unsigned int index2,
                               const unsigned int index3,
			       const unsigned int index4)
{
  this->indices[0] = index1;
  this->indices[1] = index2;
  this->indices[2] = index3;
  this->indices[3] = index4;
};



inline
TableIndices<5>::TableIndices () 
{
  this->indices[0] = this->indices[1] = this->indices[2] = this->indices[3] = this->indices[4] = 0;
};



inline
TableIndices<5>::TableIndices (const unsigned int index1,
                               const unsigned int index2,
                               const unsigned int index3,
			       const unsigned int index4,
			       const unsigned int index5)
{
  this->indices[0] = index1;
  this->indices[1] = index2;
  this->indices[2] = index3;
  this->indices[3] = index4;
  this->indices[4] = index5;
};



inline
TableIndices<6>::TableIndices () 
{
  this->indices[0] = this->indices[1] = this->indices[2]
		   = this->indices[3] = this->indices[4] = 0;
};



inline
TableIndices<6>::TableIndices (const unsigned int index1,
                               const unsigned int index2,
                               const unsigned int index3,
			       const unsigned int index4,
			       const unsigned int index5,
			       const unsigned int index6)
{
  this->indices[0] = index1;
  this->indices[1] = index2;
  this->indices[2] = index3;
  this->indices[3] = index4;
  this->indices[4] = index5;
  this->indices[5] = index6;
};



template <int N, typename T>
TableBase<N,T>::TableBase ()
                :
                val (0),
                val_size (0)
{};



template <int N, typename T>
TableBase<N,T>::TableBase (const TableIndices<N> &sizes)
                :
                val (0),
                val_size (0)
{
  reinit (sizes);
};



template <int N, typename T>
TableBase<N,T>::TableBase (const TableBase<N,T> &src)
                :
                Subscriptor (),
                val (0),
                val_size (0)
{
  reinit (src.table_size);
  if (src.n_elements() != 0)
    fill (src.data());
};



template <int N, typename T>
template <typename T2>
TableBase<N,T>::TableBase (const TableBase<N,T2> &src)
                :
                val (0),
                val_size (0)
{
  reinit (src.table_size);
  fill (src.data());
};



namespace TableBaseAccessors 
{
  template <int N, typename T, bool C, unsigned int P>
  inline
  Accessor<N,T,C,P>::Accessor (const TableType &table,
                               const pointer    data)
                  :
                  table (table),
                  data (data)
  {};



  template <int N, typename T, bool C, unsigned int P>
  inline
  Accessor<N,T,C,P>::Accessor (const Accessor &)
                  :
                  table (*static_cast<const TableType*>(0)),
                  data (0)
  {
                                     // accessor objects are only
                                     // temporary objects, so should
                                     // not need to be copied around
    Assert (false, ExcInternalError());
  };
  

  
  template <int N, typename T, bool C, unsigned int P>
  inline
  Accessor<N,T,C,P>::Accessor ()
                  :
                  table (*static_cast<const TableType*>(0)),
                  data (0)
  {
                                     // accessor objects are only
                                     // temporary objects, so should
                                     // not need to be copied around
    Assert (false, ExcInternalError());
  };
  

  
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
        unsigned int subobject_size = table.size()[N-1];
        for (int p=P-1; p>1; --p)
          subobject_size *= table.size()[N-p];
        const pointer new_data = data + i*subobject_size;
        return Accessor<N,T,C,P-1> (table, new_data);
      };
  };



  template <int N, typename T, bool C>
  inline
  Accessor<N,T,C,1>::Accessor (const TableType &table,
                               const pointer    data)
                  :
                  table (table),
                  data (data)
  {};



  template <int N, typename T, bool C>
  inline
  Accessor<N,T,C,1>::Accessor ()
                  :
                  table (*static_cast<const TableType*>(0)),
                  data (0)
  {
                                     // accessor objects are only
                                     // temporary objects, so should
                                     // not need to be copied around
    Assert (false, ExcInternalError());
  };
  


  template <int N, typename T, bool C>
  inline
  Accessor<N,T,C,1>::Accessor (const Accessor &)
                  :
                  table (*static_cast<const TableType*>(0)),
                  data (0)
  {
                                     // accessor objects are only
                                     // temporary objects, so should
                                     // not need to be copied around
    Assert (false, ExcInternalError());
  };  


  
  template <int N, typename T, bool C>
  inline
  typename Accessor<N,T,C,1>::reference
  Accessor<N,T,C,1>::operator [] (const unsigned int i) const 
  {
    Assert (i < table.size()[N-1],
            ExcIndexRange (i, 0, table.size()[N-1]));
    return data[i];
  };


  
  template <int N, typename T, bool C>
  inline
  unsigned int
  Accessor<N,T,C,1>::size () const
  {
    return table.size()[N-1];
  };



  template <int N, typename T, bool C>
  inline
  typename Accessor<N,T,C,1>::iterator
  Accessor<N,T,C,1>::begin () const
  {
    return data;
  };



  template <int N, typename T, bool C>
  inline
  typename Accessor<N,T,C,1>::iterator
  Accessor<N,T,C,1>::end () const
  {
    return data+table.size()[N-1];
  };
};



template <int N, typename T>
inline
TableBase<N,T>::~TableBase ()
{
  if (val != 0)
    delete[] val;
};



template <int N, typename T>
TableBase<N,T>&
TableBase<N,T>::operator = (const TableBase<N,T>& m) 
{
  reinit (m.size());
  if (!empty())
    std::copy (&m.val[0], &m.val[n_elements()], &val[0]);
  
  return *this;
}



template <int N, typename T>
template <typename T2>
TableBase<N,T>&
TableBase<N,T>::operator = (const TableBase<N,T2>& m) 
{
  reinit (m.size());
  if (!empty())
    std::copy (&m.val[0], &m.val[n_elements()], &val[0]);
  
  return *this;
}



template <int N, typename T>
inline
void
TableBase<N,T>::clear ()
{
  if (n_elements() != 0)
    std::fill_n (val, n_elements(), T());
};



template <int N, typename T>
void
TableBase<N,T>::reinit (const TableIndices<N> &new_sizes)
{
  table_size = new_sizes;
  
  const unsigned int new_size = n_elements();
  
                                   // if zero size was given: free all
                                   // memory
  if (new_size == 0)
    {
      if (val != 0)
        delete[] val;

      val      = 0;
      val_size = 0;

                                       // set all sizes to zero, even
                                       // if one was previously
                                       // nonzero. This simplifies
                                       // some assertions.
      table_size = TableIndices<N>();

      return;
    };
  
                                   // if new size is nonzero:
                                   // if necessary allocate
                                   // additional memory
  if (val_size<new_size)
    {
      if (val != 0)
        delete[] val;

      val_size = new_size;
      val      = new T[val_size];
    };

                                   // reinitialize contents of old or
                                   // new memory.
  clear ();
};



template <int N, typename T>
const TableIndices<N> &
TableBase<N,T>::size () const
{
  return table_size;
};



template <int N, typename T>
unsigned int
TableBase<N,T>::n_elements () const
{
  unsigned s = 1;
  for (unsigned int n=0; n<N; ++n)
    s *= table_size[n];
  return s;
};



template <int N, typename T>
bool
TableBase<N,T>::empty () const
{
  return (n_elements() == 0);
};



template <int N, typename T>
template <typename T2>
inline
void
TableBase<N,T>::fill (const T2* entries)
{
  Assert (n_elements() != 0,
          ExcMessage("Trying to fill an empty matrix."));

  std::copy (entries, entries+n_elements(), val);
}



template <int N, typename T>
unsigned int
TableBase<N,T>::memory_consumption () const
{
  return sizeof(*this) + val_size*sizeof(T);
}


template <int N, typename T>
inline
unsigned int
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
            return indices[0]*table_size[1] + indices[1];
      case 3:
            return ((indices[0]*table_size[1] + indices[1])*table_size[2]
                    + indices[2]);
      default:
      {
        unsigned int s = indices[0];
        for (unsigned int n=1; n<N; ++n)
          s = s*table_size[n] + indices[n];
        return s;
      };
    };
};




template <int N, typename T>
inline const T &
TableBase<N,T>::operator() (const TableIndices<N> &indices) const
{
  for (unsigned int n=0; n<N; ++n)
    Assert (indices[n] < table_size[n],
            ExcIndexRange (indices[n], 0, table_size[n]));
  return el(indices);
};



template <int N, typename T>
inline T &
TableBase<N,T>::operator() (const TableIndices<N> &indices)
{
  for (unsigned int n=0; n<N; ++n)
    Assert (indices[n] < table_size[n],
            ExcIndexRange (indices[n], 0, table_size[n]));
  return el(indices);
};



template <int N, typename T>
inline const T &
TableBase<N,T>::el (const TableIndices<N> &indices) const
{  
  return val[position(indices)];
};



template <int N, typename T>
inline T &
TableBase<N,T>::el (const TableIndices<N> &indices)
{
  return val[position(indices)];
};



template <int N, typename T>
inline
const T *
TableBase<N,T>::data () const
{
  return val;
}




template <typename T>
Table<1,T>::Table () 
{};



template <typename T>
Table<1,T>::Table (const unsigned int size)
                :
                TableBase<1,T> (TableIndices<1> (size))
{};



template <typename T>
const T &
Table<1,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return this->val[i];
};



template <typename T>
T &
Table<1,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return this->val[i];
};



template <typename T>
inline
const T &
Table<1,T>::operator () (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return this->val[i];
};



template <typename T>
inline
T &
Table<1,T>::operator () (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return this->val[i];
};




template <typename T>
Table<2,T>::Table ()
{};



template <typename T>
Table<2,T>::Table (const unsigned int size1,
                   const unsigned int size2)
                :
                TableBase<2,T> (TableIndices<2> (size1, size2))
{};



template <typename T>
void
Table<2,T>::reinit (const unsigned int size1,
                    const unsigned int size2)
{
  this->TableBase<2,T>::reinit (TableIndices<2> (size1, size2));
};



template <typename T>
inline
TableBaseAccessors::Accessor<2,T,true,1>
Table<2,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return TableBaseAccessors::Accessor<2,T,true,1>(*this,
                                                  this->val+i*n_cols());
};



template <typename T>
inline
TableBaseAccessors::Accessor<2,T,false,1>
Table<2,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  return TableBaseAccessors::Accessor<2,T,false,1>(*this,
                                                   this->val+i*n_cols());
};



template <typename T>
inline
const T &
Table<2,T>::operator () (const unsigned int i,
                         const unsigned int j) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  return this->val[i*this->table_size[1]+j];
};



template <typename T>
inline
T &
Table<2,T>::operator () (const unsigned int i,
                         const unsigned int j)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  Assert (j < this->table_size[1],
          ExcIndexRange (j, 0, this->table_size[1]));
  return this->val[i*this->table_size[1]+j];
};



template <typename T>
inline
const T &
Table<2,T>::el (const unsigned int i,
                const unsigned int j) const
{
  return this->val[i*this->table_size[1]+j];
};



template <typename T>
inline
T &
Table<2,T>::el (const unsigned int i,
                const unsigned int j)
{
  return this->val[i*this->table_size[1]+j];
};



template <typename T>
inline
unsigned int
Table<2,T>::n_rows () const
{
  return this->table_size[0];
};



template <typename T>
inline
unsigned int
Table<2,T>::n_cols () const
{
  return this->table_size[1];
};





template <typename T>
Table<3,T>::Table () 
{};



template <typename T>
Table<3,T>::Table (const unsigned int size1,
                   const unsigned int size2,
                   const unsigned int size3)
                :
                TableBase<3,T> (TableIndices<3> (size1, size2, size3))
{};



template <typename T>
inline
TableBaseAccessors::Accessor<3,T,true,2>
Table<3,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const unsigned int subobject_size = this->table_size[1] *
                                      this->table_size[2];
  return (TableBaseAccessors::Accessor<3,T,true,2>
          (*this,
           this->val+i*subobject_size));
};



template <typename T>
inline
TableBaseAccessors::Accessor<3,T,false,2>
Table<3,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const unsigned int subobject_size = this->table_size[1] *
                                      this->table_size[2];
  return (TableBaseAccessors::Accessor<3,T,false,2>
          (*this,
           this->val+i*subobject_size));
};



template <typename T>
inline
const T &
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
  return this->val[(i*this->table_size[1]+j)
		  *this->table_size[2] + k];
};



template <typename T>
inline
T &
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
  return this->val[(i*this->table_size[1]+j)
		  *this->table_size[2] + k];
};



template <typename T>
Table<4,T>::Table () 
{};



template <typename T>
Table<4,T>::Table (const unsigned int size1,
                   const unsigned int size2,
                   const unsigned int size3,
		   const unsigned int size4)
                :
                TableBase<4,T> (TableIndices<4> (size1, size2, size3, size4))
{};



template <typename T>
inline
TableBaseAccessors::Accessor<4,T,true,3>
Table<4,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const unsigned int subobject_size = this->table_size[1] *
                                      this->table_size[2] *
				      this->table_size[3];
  return (TableBaseAccessors::Accessor<4,T,true,3>
          (*this,
           this->val+i*subobject_size));
};



template <typename T>
inline
TableBaseAccessors::Accessor<4,T,false,3>
Table<4,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const unsigned int subobject_size = this->table_size[1] *
                                      this->table_size[2] *
				      this->table_size[3];
  return (TableBaseAccessors::Accessor<4,T,false,3>
          (*this,
           this->val+i*subobject_size));
};



template <typename T>
inline
const T &
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
  return this->val[((i*this->table_size[1]+j)
		    *this->table_size[2] + k)
		  *this->table_size[3] + l];
};



template <typename T>
inline
T &
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
  return this->val[((i*this->table_size[1]+j)
		    *this->table_size[2] + k)
		  *this->table_size[3] + l];
};




template <typename T>
Table<5,T>::Table () 
{};



template <typename T>
Table<5,T>::Table (const unsigned int size1,
                   const unsigned int size2,
                   const unsigned int size3,
		   const unsigned int size4,
		   const unsigned int size5)
                :
                TableBase<5,T> (TableIndices<5> (size1, size2, size3, size4, size5))
{};



template <typename T>
inline
TableBaseAccessors::Accessor<5,T,true,4>
Table<5,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const unsigned int subobject_size = this->table_size[1] *
                                      this->table_size[2] *
				      this->table_size[3] *
				      this->table_size[4];
  return (TableBaseAccessors::Accessor<5,T,true,4>
          (*this,
           this->val+i*subobject_size));
};



template <typename T>
inline
TableBaseAccessors::Accessor<5,T,false,4>
Table<5,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const unsigned int subobject_size = this->table_size[1] *
                                      this->table_size[2] *
				      this->table_size[3] *
				      this->table_size[4];
  return (TableBaseAccessors::Accessor<5,T,false,4>
          (*this,
           this->val+i*subobject_size));
};



template <typename T>
inline
const T &
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
  return this->val[(((i*this->table_size[1]+j)
		     *this->table_size[2] + k)
		    *this->table_size[3] + l)
		  *this->table_size[4] + m];
};



template <typename T>
inline
T &
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
  return this->val[(((i*this->table_size[1]+j)
		     *this->table_size[2] + k)
		    *this->table_size[3] + l)
		  *this->table_size[4] + m];
};



template <typename T>
Table<6,T>::Table () 
{};



template <typename T>
Table<6,T>::Table (const unsigned int size1,
                   const unsigned int size2,
                   const unsigned int size3,
		   const unsigned int size4,
		   const unsigned int size5,
		   const unsigned int size6)
                :
                TableBase<6,T> (TableIndices<6> (size1, size2, size3, size4, size5, size6))
{};



template <typename T>
inline
TableBaseAccessors::Accessor<6,T,true,5>
Table<6,T>::operator [] (const unsigned int i) const
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const unsigned int subobject_size = this->table_size[1] *
                                      this->table_size[2] *
				      this->table_size[3] *
				      this->table_size[4] *
				      this->table_size[5];
  return (TableBaseAccessors::Accessor<6,T,true,5>
          (*this,
           this->val+i*subobject_size));
};



template <typename T>
inline
TableBaseAccessors::Accessor<6,T,false,5>
Table<6,T>::operator [] (const unsigned int i)
{
  Assert (i < this->table_size[0],
          ExcIndexRange (i, 0, this->table_size[0]));
  const unsigned int subobject_size = this->table_size[1] *
                                      this->table_size[2] *
				      this->table_size[3] *
				      this->table_size[4] *
				      this->table_size[5];
  return (TableBaseAccessors::Accessor<6,T,false,5>
          (*this,
           this->val+i*subobject_size));
};



template <typename T>
inline
const T &
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
  return this->val[((((i*this->table_size[1]+j)
		      *this->table_size[2] + k)
		     *this->table_size[3] + l)
		    *this->table_size[4] + m)
		  *this->table_size[5] + n];
};



template <typename T>
inline
T &
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
  return this->val[((((i*this->table_size[1]+j)
		      *this->table_size[2] + k)
		     *this->table_size[3] + l)
		    *this->table_size[4] + m)
		  *this->table_size[5] + n];
};



#endif
