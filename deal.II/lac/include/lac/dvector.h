// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier


#ifndef __lac_dvector_h
#define __lac_dvector_h
#include <stdio.h>
#ifndef __base_types_h
#include <base/types.h>
#endif
//#ifndef __base_errors_h
//#include <base/errors.h>
//#endif
#ifndef __lac_vectorbase_h
#include <lac/vectorbase.h>
#endif

#include <base/exceptions.h>



/**
 *  Double precision Vector.
 *  Memory for Components is supplied explicitly <p>
 *  ( ! Amount of memory needs not to comply with actual dimension due to reinitializations ! ) <p>
 *  - all necessary methods for Vectors are supplied <p>
 *  - operators available are `=` , `*` and `( )` <p>
 *  CONVENTIONS for used `equations` : <p>
 *  - THIS vector is always named `U` <p>
 *  - vectors are always uppercase , scalars are lowercase
 *
 * @author Roland Becker, Guido Kanschat, Franz-Theo Suttmeier, revised and extended by Wolfgang Bangerth, documented by Klaus Mampel and Wolfgang Bangerth
 */
class dVector : public VectorBase
{
  friend class dFMatrix;

  protected:

				     /**
				      * Dimension. Actual number of components
				      * contained in the vector.
				      * Get this number by calling #size()#.
				      */
    unsigned int dim;

				     /**
				      * Amount of memory actually reserved for
				      * this vector. This number may be greater
				      * than #dim# if a #reinit# was called with
				      * less memory requirements than the vector
				      * needed last time. At present #reinit#
				      * does not free memory when the number of
				      * needed elements is reduced.
				      */
    unsigned int maxdim;

				     /**
				      * Pointer to the array of components.
				      */
    double *val;

  public:

				     /**
				      * Declare iterator types just like those
				      * for the C++ standard library:
				      *
				      * Data type stored by this container.
				      */
    typedef double value_type;

				     /**
				      * Declare standard types used in all
				      * containers.
				      */
    typedef value_type* pointer;
    typedef const value_type* const_pointer;
    typedef value_type* iterator;
    typedef const value_type* const_iterator;
    typedef value_type& reference;
    typedef const value_type& const_reference;
    typedef size_t size_type;
    

				     /**
				      * @name 1: Basic Object-handling 
				      */
				     //@{
				     /**
				      *  Dummy-Constructor. Dimension=0
				      */
    dVector ();
    
				     /**
				      *   Copy-Constructor. Dimension set to that of V , <p>
				      *                     all components are copied from V
				      */
    dVector (const dVector& V);
    
				     /**
				      *        Constructor. Dimension = N (>0)
				      */
    dVector (const unsigned int N);
    
				     /**
				      *         Destructor. Clears memory
				      */
    ~dVector ();

				     /**
				      * Set all entries to zero. Equivalent to
				      * #v = 0#, but more obvious and faster.
				      */
    void clear ();
    
				     /**
				      *  U(0-N) = s       . Fill all components
				      */
    dVector& operator= (const double s);
    
				     /**
				      *  U = V            . Copy all components
				      */
    dVector& operator= (const dVector& V);
    
				     /**
				      *  U = U * V        . Scalar Produkt
				      */
    double operator* (const dVector& V) const;

				     /**
				      * Return square of the l2-norm.
				      */
    double norm_sqr () const;

				     /**
				      * Return the mean value of the elements of
				      * this vector.
				      */
    double mean_value () const;

				     /**
				      * Return the l1-norm of the vector, i.e.
				      * the sum of the absolute values.
				      */
    double l1_norm () const;

				     /**
				      * Return the l2-norm of the vector, i.e.
				      * the square root of the sum of the
				      * squares of the elements.
				      */
    double l2_norm () const;

				     /**
				      * Return the maximum absolute value of the
				      * elements of this vector.
				      */
    double linfty_norm () const;
    
    
				     /**
				      * Change the dimension of the vector to
				      * #N#. The reserved memory for this vector
				      * remains unchanged, however, to make
				      * things faster, but this may waste some
				      * memory, so take this in the back of your
				      * head.
				      *
				      * On #fast==false#, the vector is filled by
				      * zeros.
				      */ 
    void reinit (const unsigned int N, const bool fast=false);
    
				     /**
				      * Adjust  Dimension. <p>
				      * Set dimension to n(V) <p>
				      * ! reserved memory for This remains unchanged ! <p>
				      * ! components of V are not copied in any case ! <p>
				      * on fast==false vector is filled by 0.
				      */
    void reinit (const dVector& V, const bool fast=false);
    
				     /**
				      *  Inquire Dimension. returns Dimension , 
				      *             INLINE
				      */
    unsigned int size () const;

				     /**
				      * Make the #dVector# class a bit like the
				      * #vector<># class of the C++ standard
				      * library by returning iterators to
				      * the start and end of the elements of this
				      * vector.
				      */
    iterator begin ();

				     /**
				      * Return constant iterator to the start of
				      * the vectors.
				      */
    const_iterator begin () const;

				     /**
				      * Return an iterator pointing to the
				      * element past the end of the array.
				      */
    iterator end ();

    				     /**
				      * Return a constant iterator pointing to
				      * the element past the end of the array.
				      */
    const_iterator end () const;  
				     //@}
    
    
				     /**
				      * @name 2: Data-Access
				      */
				     //@{
				     /**
				      *  Access Components. returns U(i) , 
				      *             INLINE
				      */
    double operator() (const unsigned int i) const;
    
				     /**
				      *  Access Components. returns U(i) , 
				      *             INLINE
				      */
    double& operator() (const unsigned int i);
				     //@}
    
    
				     /**
				      * @name 3: Modification of vectors
				      */
				     //@{
				     /**
				      * Fast equivalent to #U.add(1, V)#.
				      */
    dVector & operator += (const dVector &V);

    				     /**
				      * Fast equivalent to #U.add(-1, V)#.
				      */
    dVector & operator -= (const dVector &V);

				     /**
				      *  U(0-DIM)+=s      . Addition of S to all components
				      */
    void add (const double s);
    
				     /**
				      *  U+=V             . Simple addition
				      */
    void add (const dVector& V);
    
				     /**
				      *  U+=a*V           . Simple addition
				      */
    void add (const double a, const dVector& V);
    
				     /**
				      *  U+=a*V+b*W       . Multiple addition
				      */
    void add (const double a, const dVector& V,
	      const double b, const dVector& W);
    
				     /**
				      *  U=s*U+V          . Scaling + simple addition
				      */
    void sadd (const double s, const dVector& V);
    
				     /**
				      *  U=s*U+a*V        . Scaling + simple addition
				      */
    void sadd (const double s, const double a, const dVector& V);
    
				     /**
				      *  U=s*U+a*V+b*W    . Scaling + multiple addition
				      */
    void sadd (const double s, const double a,
	       const dVector& V, const double b, const dVector& W);
    
				     /**
				      *  U=s*U+a*V+b*W+c*X. Scaling + multiple addition
				      */
    void sadd (const double s, const double a,
	       const dVector& V, const double b, const dVector& W, 
	       const double c, const dVector& X);
    
				     /**
				      * Scale each element of the vector by the
				      * given factor. This function was
				      * previously called #equ(double)#, which
				      * in my eyes is an extremely unintuitive
				      * naming and was thus replaced.
				      */
    void scale (const double factor);
    
				     /**
				      *  U=a*V            . Replacing
				      */
    void equ (const double a, const dVector& V);
    
				     /**
				      *  U=a*V+b*W        . Replacing by sum
				      */
    void equ (const double a, const dVector& V,
	      const double b, const dVector& W);

				     /**
				      * Compute the elementwise ratio of the
				      * two given vectors, that is let
				      * #this[i] = a[i]/b[i]#. This is useful
				      * for example if you want to compute
				      * the cellwise ratio of true to estimated
				      * error.
				      *
				      * This vector is appropriately scaled to
				      * hold the result.
				      *
				      * If any of the #b[i]# is zero, the result
				      * is undefined. No attempt is made to
				      * catch such situations.
				      */
    void ratio (const dVector &a, const dVector &b);
				     //@}
    
    
				     /**
				      * @name 4: Modification of components
				      */
				     //@{
				     /**
				      *  U(i)=0             . ONE Component only
				      */
    void czero (const unsigned int i);
    
				     /**
				      *  U(i)=a*V(j)        . Replacing
				      */
    void cequ(unsigned int i, const VectorBase& V,
	      double a, unsigned int j);
    
				     /**
				      *  U(i)=a*V(j)+b*V(k) . Replacing by sum
				      */
    void cequ (const unsigned int i, const VectorBase& V,
	       const double a, const unsigned int j,
	       const double b, const unsigned int k);
    
				     /**
				      *  U(i)=a*V(j)+b*V(k)+c*V(l)+d*V(m) . Replacing by sum
				      */
    void cequ (const unsigned int i, const VectorBase& V,
	       const double a, const unsigned int j,
	       const double b, const unsigned int k,
	       const double c, const unsigned int l,
	       const double d, const unsigned int m);
    
				     /**
				      *  U(i)+=a*V(j)       . Simple addition
				      */
    void cadd (const unsigned int i, const VectorBase& V,
	       const double a, const unsigned int j);
    
				     /**
				      *  U(i)+=a*V(j)+b*V(k). Multiple addition
				      */ 
    void cadd (const unsigned int i, const VectorBase& V,
	       const double a, const unsigned int j,
	       const double b, const unsigned int k);
    
				     /**
				      *  U(i)+=a*V(j)+b*V(k)+c*V(l)+d*V(m) . Multiple addition
				      */
    void cadd (const unsigned int i, const VectorBase& V,
	       const double a, const unsigned int j,
	       const double b, const unsigned int k,
	       const double c, const unsigned int l,
	       const double d, const unsigned int m);
				     //@}
    
    
				     /**
				      * @name 5: Mixed stuff
				      */
				     //@{
				     ///
    virtual const char* name () const;
    
				     /**
				      *  Output of vector in user-defined format.
				      */
    void print (FILE* fp, const char* format = 0) const;
    
				     /**
				      *  Output of vector in user-defined format.
				      */
    void print (const char* format = 0) const;

				     /**
				      * Print to given stream, one element per line.
				      */
    void print (ostream &) const;
				     //@}

				     /**
				      * Exception
				      */
    DeclException2 (ExcDimensionsDontMatch,
		    int, int,
		    << "The dimensions " << arg1 << " and " << arg2
		    << " do not match here.");
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidIndex,
		    int, int,
		    << "The given index " << arg1
		    << " should be less than " << arg2 << ".");
				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidNumber,
		    int,
		    << "The provided number is invalid here: " << arg1);
				     /**
				      * Exception
				      */
    DeclException0 (ExcOutOfMemory);
};





inline unsigned int dVector::size () const
{
  return dim;
}



inline
dVector::iterator dVector::begin () {
  return &val[0];
};



inline
dVector::const_iterator dVector::begin () const {
  return &val[0];
};



inline
dVector::iterator dVector::end () {
  return &val[dim];
};



inline
dVector::const_iterator dVector::end () const {
  return &val[dim];
};



inline double dVector::operator() (const unsigned int i) const
{
  Assert (i<dim, ExcInvalidIndex(i,dim));
  return val[i];
}



inline double& dVector::operator() (const unsigned int i)
{
  Assert (i<dim, ExcInvalidIndex(i,dim));
  return val[i];
}

#endif

