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
    int dim;

				     /// Dimension. Determines amount of reserved memory , evtl. >DIM !
    int maxdim;

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
    dVector();
    
				     /**
				      *   Copy-Constructor. Dimension set to that of V , <p>
				      *                     all components are copied from V
				      */
    dVector(const dVector& V);
    
				     /**
				      *        Constructor. Dimension = N (>0)
				      */
    dVector(int N);
    
				     /**
				      *         Destructor. Clears memory
				      */
    ~dVector();
    
				     /**
				      *  U(0-N) = s       . Fill all components
				      */
    dVector& operator=(double s);
    
				     /**
				      *  U = V            . Copy all components
				      */
    dVector& operator=(const dVector& V);
    
				     /**
				      *  U = U * V        . Scalar Produkt
				      */
    double operator*(const dVector& V) const;
    
				     /**
				      * Change  Dimension. <p>
				      * Set dimension to N <p>
				      * ! reserved memory for This remains unchanged ! <p>
				      * on fast=0 vector is filled by 0.
				      */ 
    void reinit(int N, int fast = 0);
    
				     /**
				      * Adjust  Dimension. <p>
				      * Set dimension to n(V) <p>
				      * ! reserved memory for This remains unchanged ! <p>
				      * ! components of V are not copied in any case ! <p>
				      * on fast=0 vector is filled by 0.
				      */
    void reinit(const dVector& V, int fast = 0);
    
				     /**
				      *  Inquire Dimension. returns Dimension , 
				      *             INLINE
				      */
    int n() const;
				     //@}
    
    
				     /**
				      * @name 2: Data-Access
				      */
				     //@{
				     /**
				      *  Access Components. returns U(i) , 
				      *             INLINE
				      */
    double operator()(int i) const;
    
				     /**
				      *  Access Components. returns U(i) , 
				      *             INLINE
				      */
    double& operator()(int i);
				     //@}
    
    
				     /**
				      * @name 3: Modification of vectors
				      */
				     //@{
				     /**
				      *  U(0-DIM)+=s      . Addition of S to all components
				      */
    void add(const double s);
    
				     /**
				      *  U+=V             . Simple addition
				      */
    void add(const dVector& V);
    
				     /**
				      *  U+=a*V           . Simple addition
				      */
    void add(double a, const dVector& V);
    
				     /**
				      *  U+=a*V+b*W       . Multiple addition
				      */
    void add(double a, const dVector& V, double b, const dVector& W);
    
				     /**
				      *  U=s*U+V          . Scaling + simple addition
				      */
    void sadd(double s, const dVector& V);
    
				     /**
				      *  U=s*U+a*V        . Scaling + simple addition
				      */
    void sadd(double s, double a, const dVector& V);
    
				     /**
				      *  U=s*U+a*V+b*W    . Scaling + multiple addition
				      */
    void sadd(double s, double a, const dVector& V, double b, const dVector& W);
    
				     /**
				      *  U=s*U+a*V+b*W+c*X. Scaling + multiple addition
				      */
    void sadd(double s, double a, const dVector& V, double b, const dVector& W, 
	      double c, const dVector& X);
    
				     /**
				      *  U=s*U            . Scaling
				      */
    void equ(double s);
    
				     /**
				      *  U=a*V            . Replacing
				      */
    void equ(double a, const dVector& V);
    
				     /**
				      *  U=a*V+b*W        . Replacing by sum
				      */
    void equ(double a, const dVector& V, double b, const dVector& W);
				     //@}
    
    
				     /**
				      * @name 4: Modification of components
				      */
				     //@{
				     /**
				      *  U(i)=0             . ONE Component only
				      */
    void czero(int i);
    
				     /**
				      *  U(i)=a*V(j)        . Replacing
				      */
    void cequ(int i, const VectorBase& V, double a, int j);
    
				     /**
				      *  U(i)=a*V(j)+b*V(k) . Replacing by sum
				      */
    void cequ(int i, const VectorBase& V, double a, int j, double b, int k);
    
				     /**
				      *  U(i)=a*V(j)+b*V(k)+c*V(l)+d*V(m) . Replacing by sum
				      */
    void cequ(int i, const VectorBase& V, double a, int j, double b,
	      int k, double c, int l, double d, int m);
    
				     /**
				      *  U(i)+=a*V(j)       . Simple addition
				      */
    void cadd(int i, const VectorBase& V, double a, int j);
    
				     /**
				      *  U(i)+=a*V(j)+b*V(k). Multiple addition
				      */ 
    void cadd(int i, const VectorBase& V, double a, int j, double b, int k);
    
				     /**
				      *  U(i)+=a*V(j)+b*V(k)+c*V(l)+d*V(m) . Multiple addition
				      */
    void cadd(int i, const VectorBase& V, double a, int j, double b,
	      int k, double c, int l, double d, int m);
				     //@}
    
    
				     /**
				      * @name 5: Mixed stuff
				      */
				     //@{
				     ///
    virtual const char* name() const;
    
				     /**
				      *  Output of vector in user-defined format.
				      */
    void print(FILE* fp, const char* format = 0) const;
    
				     /**
				      *  Output of vector in user-defined format.
				      */
    void print(const char* format = 0) const;
				     //@}
};





inline int dVector::n() const
{
  return dim;
}

inline double dVector::operator() (int i) const
{
				    //THROW2( (i<0) || (i>=dim), IntError(IntError::Range,i));
  return val[i];
}

inline double& dVector::operator() (int i)
{
				    //THROW2( (i<0) || (i>=dim), IntError(IntError::Range,i));
  return val[i];
}

#endif
