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
#ifndef __deal2__vector2d_h
#define __deal2__vector2d_h

#include <base/config.h>
#include <base/subscriptor.h>


template <typename T> class vector2d;
template <typename T> class FullMatrix;


/**
 * Two-dimensional array of arbitrary data.
 *
 * This is an implementation of a two-dimensional array with access by
 * pairs of coordinates. Access is either by @p{x(i,j)} to keep with
 * matrix notation, or @p{x[i][j]} in accordance to C-style arrays.
 *
 * The name is chosen to be @p{vector2d} to be conformant with the
 * standard C++ classes. This, although array would have been
 * preferrable, since a vector should be element of a vector space.
 *
 * @author Guido Kanschat, 2001. Parts by Wolfgang Bangerth and others.
 */
template<typename T>
class vector2d : public Subscriptor
{
  public:
				     /**
				      * Object representing one row of
				      * a @ref{vector2d} object. It
				      * allows to access the elements
				      * through the
				      * @p{operator[]}. Since objects
				      * of this type are also
				      * generated through the
				      * @p{vector2d::operator[]}, this
				      * allows accessing the elements
				      * of @p{vector2d} objects just
				      * as those of two-dimensional
				      * C-style arrays.
				      *
				      * This class is used for
				      * constant @p{vector2d}
				      * objects. It only allows
				      * read-only access to the
				      * elements.
				      */
    class ConstRowAccessor
    {
      public:
					 /**
					  * Constructor.
					  */
	ConstRowAccessor (const vector2d<T>  &parent,
			  const unsigned int  row);

					 /**
					  * Access operator.
					  */
	T  operator [] (const unsigned int column) const;
	
      protected:
					 /**
					  * Pointer to the parent
					  * object. Used only to check
					  * access indices for
					  * validity.
					  */
	const vector2d<T>  &parent;

					 /**
					  * Pointer to the start of
					  * the row pointed to by this
					  * object.
					  */
	const T * const     row_start;
    };    


				     /**
				      * Object representing one row of
				      * a @ref{vector2d} object. It
				      * allows to access the elements
				      * through the
				      * @p{operator[]}. Since objects
				      * of this type are also
				      * generated through the
				      * @p{vector2d::operator[]}, this
				      * allows accessing the elements
				      * of @p{vector2d} objects just
				      * as those of two-dimensional
				      * C-style arrays.
				      *
				      * This class is used for
				      * non-constant @p{vector2d}
				      * objects. It allows read and
				      * write access to the elements.
				      */
    class NonConstRowAccessor
    {
      public:
					 /**
					  * Constructor.
					  */
	NonConstRowAccessor (vector2d<T>        &parent,
			     const unsigned int  row);


					 /**
					  * Access operator.
					  */
	T & operator [] (const unsigned int column) const;
	
      private:
					 /**
					  * Pointer to the parent
					  * object. Used only to check
					  * access indices for
					  * validity.
					  */
	const vector2d<T>  &parent;

					 /**
					  * Pointer to the start of
					  * the row pointed to by this
					  * object.
					  */
	T * const           row_start;
    };    
				     /**
				      * Constructor for a quadratic
				      * @p{rows} by @p{rows} array. The
				      * standard constructor creates an
				      * empty object.
				      */
    vector2d (const unsigned int rows = 0);
    
				     /**
				      * Constructor for a @p{rows} by
				      * @p{cols} array.
				      */
    vector2d (const unsigned int rows,
	      const unsigned int cols);
    
				     /**
				      * Copy-constructor for deep copy.
				      */
    vector2d (const vector2d&);
 				     /**
				      * Constructor initializing from
				      * an array of data elements. The array
				      * is arranged line by line. No
				      * range checking is performed.
				      */
    vector2d (const unsigned int rows,
	      const unsigned int cols,
	      const T* entries);
    
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
    vector2d<T>& operator = (const vector2d<T>& src);
    
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
    vector2d<T>& operator = (const vector2d<T2> &src);
				      
				     /**
				      * Destructor. Free allocated memory.
				      */
    ~vector2d ();
    
				     /**
				      * Set dimension to $m\times n$
				      * and allocate memory if
				      * necessary. Forget the previous
				      * content of the array.
				      */
    void reinit (const unsigned int m,
		 const unsigned int n);
    
				     /**
				      * Set dimension to $n\times n$
				      * and allocate memory if
				      * necessary. Forget the previous
				      * content of the array.
				      */
    void reinit (const unsigned int n);

				     /**
				      * Set the dimensions to the same
				      * as another array. The other
				      * array will not be copied,
				      * though. The entries of this
				      * array will be zero.
				      */
    template <class T2>
    void reinit (const vector2d<T2> &shape);
    
				     /**
				      * Number of rows.
				      */
    unsigned int n_rows () const;
    
				     /**
				      * Number of columns.
				      */
    unsigned int n_cols () const;
    
				     /**
				      * Return the value of the
				      * element at position @p{(i,j)}.
				      */
    T operator() (const unsigned int i,
		  const unsigned int j) const;
    
				     /**
				      * Return a read-write reference to
				      * the element at position @p{(i,j)}.
				      */
    T & operator() (const unsigned int i,
		    const unsigned int j);

				     /**
				      * Return an object representing
				      * a certain row of this
				      * array. This object in turn has
				      * an @p{operator[]}, so that the
				      * elements of this array can be
				      * accessed through @p{x[i][j]}
				      * as well as through @p{x(i,j)}.
				      *
				      * This function is for constant
				      * objects. The returned row
				      * representing object only
				      * allows read access to the
				      * elements of the row pointed
				      * to.
				      */
    ConstRowAccessor    operator [] (const unsigned int row) const;

				     /**
				      * Return an object representing
				      * a certain row of this
				      * array. This object in turn has
				      * an @p{operator[]}, so that the
				      * elements of this array can be
				      * accessed through @p{x[i][j]}
				      * as well as through @p{x(i,j)}.
				      *
				      * This function is for
				      * non-constant objects. The
				      * returned row representing
				      * object allows read and write
				      * access to the elements of the
				      * row pointed to.
				      */
    NonConstRowAccessor operator [] (const unsigned int row);
    

				     /**
				      * Set all entries to their
				      * default value (zero).
				      */
    void clear ();
    
                                     /**
				      * Fill array with an array of
				      * elements.  The array is
				      * arranged line by line. No
				      * range checking is performed.
				      */
    template<typename T2>
    void fill (const T2 *entries);
    
				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

  protected:
				     /**
				      * Return a read-write reference
				      * to the element @p{(i,j)}.
				      *
				      * This function does no bounds
				      * checking and is only to be
				      * used internally and in
				      * functions already checked.
				      */
    T & el (const unsigned int i, const unsigned int j);
  
				     /**
				      * Return the value of the
				      * element @p{(i,j)}.
				      *
				      * This function does no bounds
				      * checking and is only to be
				      * used internally and in
				      * functions already checked.
				      */
    T el (const unsigned int i, const unsigned int j) const;    
  
				     /**
				      * Direct read-only access to
				      * data field. Used by
				      * @ref{FullMatrix} (there even
				      * with a cast from const),
				      * otherwise, keep away!
				      */
    const T* data () const;
    
  private:
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
				      * Number of Columns
				      */
    unsigned int num_cols;
    
				     /**
				      * Number of Rows
				      */
    unsigned int num_rows;
				     /**
				      * Friend declaration needed for
				      * inter-type copy operator.
				      */
    template <typename T2> friend class vector2d;
    
				     /**
				      * This is unfortunately needed.
				      */
    template <typename T2> friend class FullMatrix;

    friend class ConstRowAccessor;
    friend class NonConstRowAccessor;    
};


/* --------------------- Template and inline functions ---------------- */

template <typename T>
inline
vector2d<T>::ConstRowAccessor::
ConstRowAccessor (const vector2d<T>  &parent,
		  const unsigned int  row)
			:
		parent (parent),
		row_start(&parent.val[row*parent.n_cols()])
{
  Assert (row < parent.n_rows(), ExcInternalError());
};


template <typename T>
inline
T
vector2d<T>::ConstRowAccessor::
operator [] (const unsigned int column) const
{
  Assert (column < parent.n_cols(), ExcInternalError());
  return *(row_start+column);
};



template <typename T>
inline
vector2d<T>::NonConstRowAccessor::
NonConstRowAccessor (vector2d<T>        &parent,
		     const unsigned int  row)
			:
		parent (parent),
		row_start(&parent.val[row*parent.n_cols()])
{
  Assert (row < parent.n_rows(), ExcInternalError());
};


template <typename T>
inline
T &
vector2d<T>::NonConstRowAccessor::
operator [] (const unsigned int column) const
{
  Assert (column < parent.n_cols(), ExcInternalError());
  return *(row_start+column);
};



template <typename T>
inline
vector2d<T>::~vector2d ()
{
  if (val != 0)
    delete[] val;
};



template <typename T>
inline
void
vector2d<T>::clear ()
{
  if (val != 0)
    std::fill_n (val, num_rows*num_cols, T());
};



template <typename T>
template <typename T2>
inline
void
vector2d<T>::fill (const T2* entries)
{
  if (val_size != 0)
    std::copy (entries, entries+(num_rows*num_cols), val);
}



template <typename T>
void
vector2d<T>::reinit (const unsigned int mm,
		     const unsigned int nn)
{
  num_cols = nn;
  num_rows = mm;

				   // if zero size was given: free all
				   // memory
  if ((num_cols==0) || (num_rows == 0))
    {
      if (val != 0)
	delete[] val;

      val      = 0;
      val_size = 0;

				       // set both sizes to zero, even
				       // if one was previously
				       // nonzero. This simplifies
				       // some Assertions.
      num_cols = num_rows = 0;

      return;
    };
  
				   // if new size is nonzero:
				   // if necessary: allocate
				   // additional memory
  if (val_size<nn*mm)
    {
      if (val != 0)
	delete[] val;

      val_size = num_cols * num_rows;
      val      = new T[val_size];
    };

				   // Clear contents of old or new
				   // memory.
  clear ();
};



template <typename T>
void
vector2d<T>::reinit (const unsigned int n)
{
  reinit (n, n);
};



template <typename T>
template <typename T2>
void
vector2d<T>::reinit (const vector2d<T2> &B)
{
  reinit (B.n_rows(), B.n_cols());
};



template <typename T>
vector2d<T>::vector2d (const unsigned int m,
		       const unsigned int n) :
		val (0),
                val_size (0),
                num_cols (0),
                num_rows (0)
{
  reinit (m,n);
};



template <typename T>
vector2d<T>::vector2d (const unsigned int m) :
		val (0),
                val_size (0),
                num_cols (0),
                num_rows (0)
{
  reinit (m,m);
};



template <typename T>
vector2d<T>::vector2d (const unsigned int m,
				const unsigned int n,
				const T* entries) :
		val (0),
                val_size (0),
		num_cols (0),
		num_rows (0)
{
  reinit (m,n);

  if (num_cols*num_rows != 0)
    std::copy (entries, entries+num_rows*num_cols, val);
};



template <typename T>
vector2d<T>::vector2d (const vector2d &m) :
		Subscriptor (),
		val (0),
		val_size (0),
		num_cols (0),
		num_rows (0)
{
  reinit (m.num_rows, m.num_cols);
  if (num_cols*num_rows != 0)
    std::copy (&m.val[0], &m.val[num_rows*num_cols],
	       &val[0]);
};



template <typename T>
vector2d<T>&
vector2d<T>::operator = (const vector2d<T>& m) 
{
  reinit (m);
  if (num_cols*num_rows != 0)
    std::copy (&m.val[0], &m.val[num_rows*num_cols],
	       &val[0]);
  
  return *this;
}



template <typename T>
template <typename T2>
vector2d<T>&
vector2d<T>::operator = (const vector2d<T2>& m) 
{
  reinit(m);
  if (num_cols*num_rows != 0)
    copy (&m.val[0], &m.val[num_rows*num_cols],
	  &val[0]);

  return *this;
}



template <class T>
inline 
unsigned int
vector2d<T>::n_rows () const
{
  return num_rows;
  
}



template <class T>
inline 
unsigned int
vector2d<T>::n_cols () const
{
  return num_cols;
  
}



template <typename T>
inline T &
vector2d<T>::el (const unsigned int i, const unsigned int j)
{
  return val[i*n_cols()+j];
};



template <typename T>
inline T
vector2d<T>::el (const unsigned int i, const unsigned int j) const
{
  return val[i*n_cols()+j];
};



template <typename T>
inline T
vector2d<T>::operator() (const unsigned int i, const unsigned int j) const
{  
  Assert (i<num_rows, ExcIndexRange (i, 0, num_rows));
  Assert (j<num_cols, ExcIndexRange (j, 0, num_cols));
  return el(i,j);
};



template <typename T>
inline
T &
vector2d<T>::operator() (const unsigned int i, const unsigned int j)
{
  Assert (i<num_rows, ExcIndexRange (i, 0, num_rows));
  Assert (j<num_cols, ExcIndexRange (j, 0, num_cols));
  return el(i,j);
};



template <typename T>
inline
typename vector2d<T>::ConstRowAccessor 
vector2d<T>::operator [] (const unsigned int row) const 
{
				   // Note: check for validity of row
				   // number is done in the
				   // constructor of the created
				   // object
  return ConstRowAccessor(*this, row);
};



template <typename T>
inline
typename vector2d<T>::NonConstRowAccessor 
vector2d<T>::operator [] (const unsigned int row)
{
				   // Note: check for validity of row
				   // number is done in the
				   // constructor of the created
				   // object
  return NonConstRowAccessor(*this, row);
};





template <typename T>
inline
const T *
vector2d<T>::data () const
{
  return val;
}



template <typename T>
inline
unsigned int
vector2d<T>::memory_consumption () const
{
  return sizeof (*this) + val_size * sizeof (T);
}

#endif
