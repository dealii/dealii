// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier


#ifndef __lac_ivector_h
#define __lac_ivector_h
#ifndef __base_types_h
#include <base/types.h>
#endif



/**
 *  Integer Vector.
 *  Memory for Components is supplied explicitly <p>
 *  ( ! Amount of memory needs not to comply with actual dimension due to reinitializations ! ) <p>
 *  - all defined methods for iVectors are supplied <p>
 *  - operators available are `=` and `( )` <p>
 *  CONVENTIONS for used `equations` : <p>
 *  - THIS vector is always named `U` <p>
 *  - vectors are always uppercase , scalars are lowercase
 */
class iVector
{
  friend class dFMatrix;

protected:

				     /// Dimension. Actual number of components
    int dim;

				     /// Dimension. Determines amount of reserved memory , evtl. >DIM !
    int maxdim;
    
				     /// Component-array.
    int *val;
    
  public:
    
				     /**@name 1: Basic Object-handling */
				     //@{
				     /**
				      *  Dummy-Constructor. Dimension=1
				      */
    iVector();
    
				     /**
				      *   Copy-Constructor. Dimension set to that of V , <p>
				      *                     all components are copied from V
				      */
    iVector(const iVector& V);
    
				     /**
				      *        Constructor. Dimension = N (>0)
				      */
    iVector(int N);
    
				     /**
				      *         Destructor. Clears memory
				      */
    ~iVector();
    
				     /**
				      *  U(0-N) = s       . Fill all components
				      */
    iVector& operator=(int s);
    
				     /**
				      *  U = V            . Copy all components
				      */
    iVector& operator=(const iVector& V);
    
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
    void reinit(const iVector& V, int fast = 0);
    
				     /**
				      *  Inquire Dimension. returns Dimension ,
				      *             INLINE
				      */
    int n() const;
				     //@}
    
    
				     /**@name 2: Data-Access
				      */
				     //@{
				     /**
				      *  Access Components. returns U(i) ,
				      *             INLINE
				      */
    int operator()(int i) const;
    
				     /**
				      *  Access Components. returns U(i) ,
				      *             INLINE
				      */
    int& operator()(int i);
				     //@}
    
    
				     /**@name 3: Modification of vectors
				      */
				     //@{
				     /**
				      *  U+=V             . Simple addition
				      */
    void add(const iVector& V);
    
				     /**
				      *  U+=a*V           . Simple addition
				      */
    void add(int a, const iVector& V);
    
				     /**
				      *  U=a*V            . Replacing
				      */
    void equ(int a, const iVector& V);
				     //@}
};










inline int iVector::n() const
{
  return dim;
}

inline int iVector::operator() (int i) const
{
//THROW2( (i<0) || (i>=dim), IntError(IntError::Range,i));
  return val[i];
}

inline int& iVector::operator() (int i)
{
//  THROW2( (i<0) || (i>=dim), IntError(IntError::Range,i));
  return val[i];
}
#endif
