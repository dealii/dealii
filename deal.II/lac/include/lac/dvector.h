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
 */
class dVector : public VectorBase
{
  friend class dFMatrix;

  protected:

				     /// Dimension. Actual number of components
    unsigned int dim;

				     /// Dimension. Determines amount of reserved memory , evtl. >DIM !
    unsigned int maxdim;

				     /// Component-array. 
    double *val;

  public:

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
				      * Return square of the norm.
				      */
    double norm_sqr () const;
    
				     /**
				      * Change  Dimension. <p>
				      * Set dimension to N <p>
				      * ! reserved memory for This remains unchanged ! <p>
				      * on fast==false vector is filled by 0.
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
    unsigned int n () const;
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
				      *  U=s*U            . Scaling
				      */
    void equ (const double s);
    
				     /**
				      *  U=a*V            . Replacing
				      */
    void equ (const double a, const dVector& V);
    
				     /**
				      *  U=a*V+b*W        . Replacing by sum
				      */
    void equ (const double a, const dVector& V,
	      const double b, const dVector& W);
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





inline unsigned int dVector::n () const
{
  return dim;
}

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
