// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier


#ifndef __lac_ivector_h
#define __lac_ivector_h

#include <base/exceptions.h>


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
    unsigned int dim;

				     /// Dimension. Determines amount of reserved memory , evtl. >DIM !
    unsigned int maxdim;
    
				     /// Component-array.
    int *val;
    
  public:
    
				     /**@name 1: Basic Object-handling */
				     //@{
				     /**
				      *  Dummy-Constructor. Dimension=1
				      */
    iVector ();
    
				     /**
				      *   Copy-Constructor. Dimension set to that of V , <p>
				      *                     all components are copied from V
				      */
    iVector (const iVector& V);
    
				     /**
				      *        Constructor. Dimension = N (>0)
				      */
    iVector (const unsigned int N);
    
				     /**
				      *         Destructor. Clears memory
				      */
    ~iVector ();
    
				     /**
				      *  U(0-N) = s       . Fill all components
				      */
    iVector& operator = (const int s);
    
				     /**
				      *  U = V            . Copy all components
				      */
    iVector& operator = (const iVector& V);
    
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
    void reinit (const iVector& V, const bool fast=false);

				     /**
				      * Set all elements to zero.
				      */
    void clear ();
    
				     /**
				      *  Inquire Dimension. returns Dimension ,
				      *             INLINE
				      */
    unsigned int n () const;
				     //@}
    
    
				     /**@name 2: Data-Access
				      */
				     //@{
				     /**
				      *  Access Components. returns U(i) ,
				      *             INLINE
				      */
    int operator() (const unsigned int i) const;
    
				     /**
				      *  Access Components. returns U(i) ,
				      *             INLINE
				      */
    int& operator() (const unsigned int i);
				     //@}
    
    
				     /**@name 3: Modification of vectors
				      */
				     //@{
				     /**
				      *  U+=V             . Simple addition
				      */
    void add (const iVector& V);
    
				     /**
				      *  U+=a*V           . Simple addition
				      */
    void add (const unsigned int a, const iVector& V);
    
				     /**
				      *  U=a*V            . Replacing
				      */
    void equ (const unsigned int a, const iVector& V);
				     //@}

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










inline unsigned int iVector::n() const
{
  return dim;
}



inline int iVector::operator() (unsigned int i) const
{
  Assert (i<dim, ExcInvalidIndex(i,dim));
  return val[i];
}



inline int& iVector::operator() (unsigned int i)
{
  Assert (i<dim, ExcInvalidIndex(i,dim));
  return val[i];
}


#endif
