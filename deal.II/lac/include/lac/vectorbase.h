// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier


#ifndef __lac_vectorbase_h
#define __lac_vectorbase_h



/**
 * Vector Baseclass (abstract).
 *  CONVENTIONS for used `equations` : <p>
 *  - THIS vector is always named `U` <p>
 *  - vectors are always uppercase , scalars are lowercase
 */
class VectorBase
{
  public:
    
				     /**@name 1: Basic Object-handling
				      */
				     //@{
				     /**
				      *           Destructor. Should clear memory
				      */
    virtual ~VectorBase() {}
				     //@}
    
    
				     /**@name 2: Modification of components
				      */
				     //@{
				     /**
				      *  U(i)=0             . ONE Component only
				      */
    virtual void czero(int i) = 0;
    
				     /**
				      *  U(i)=a*V(j)        . Replacing
				      */
    virtual void cequ(int i, const VectorBase& V, double a, int j) = 0;
    
				     /**
				      *  U(i)=a*V(j)+b*V(k) . Replacing by sum
				      */
    virtual void cequ(int i, const VectorBase& V, double a, int j, double b, int k) = 0;
    
				     /**
				      *  U(i)=a*V(j)+b*V(k)+c*V(l)+d*V(m) . Replacing by sum
				      */
    virtual void cequ(int i, const VectorBase& V, double a, int j, double b,
		      int k, double c, int l, double d, int m) = 0;
    
				     /**
				      *  U(i)+=a*V(j)       . Simple addition
				      */
    virtual void cadd(int i, const VectorBase& V, double a, int j) = 0;
    
				     /**
				      *  U(i)+=a*V(j)+b*V(k). Multiple addition
				      */
    virtual void cadd(int i, const VectorBase& V, double a, int j, double b, int k) = 0;
    
				     /**
				      *  U(i)+=a*V(j)+b*V(k)+c*V(l)+d*V(m) . Multiple addition
				      */
    virtual void cadd(int i, const VectorBase& V, double a, int j, double b,
		      int k, double c, int l, double d, int m) = 0;
				     //@}
    
    
				     /**@name 3: Mixed stuff
				      */
				     //@{
				     ///
    virtual const char* name() const = 0;
				     //@}
};

#endif
