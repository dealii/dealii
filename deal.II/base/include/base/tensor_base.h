/*----------------------------   tensor_base.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tensor_base_H
#define __tensor_base_H
/*----------------------------   tensor_base.h     ---------------------------*/


// single this file out from tensor.h, since we want to derive Point<dim>
// from Tensor<1,dim>. However, the point class will not need all the
// tensor stuff, so we don't want the whole tensor package to be included
// everytime we use a point.

#include <base/exceptions.h>




// general template; specialized for rank==1; the general template is in
// tensor.h
template <int rank, int dim> class Tensor;



template <int dim>
class Tensor<1,dim> {
      public:
				     /**
				      * Provide a way to get the dimension of
				      * an object without explicit knowledge
				      * of it's data type. Implementation is
				      * this way instead of providing a function
				      * #dimension()# because now it is possible
				      * to get the dimension at compile time
				      * without the expansion and preevaluation
				      * of an inlined function; the compiler may
				      * therefore produce more efficient code
				      * and you may use this value to declare
				      * other data types.
				      */
    static const unsigned int dimension = dim;

				     /**
				      * Publish the rank of this tensor to
				      * the outside world.
				      */
    static const unsigned int rank      = 1;

				     /**
				      * Constructor. Initialize all entries
				      * to zero.
				      */
    explicit Tensor ();
    
    				     /**
				      *  Copy constructor.
				      */
    Tensor (const Tensor<1,dim> &);

				     /**
				      *  Read access to the #index#th coordinate.
				      *
				      * Note that the derived #Point# class also
				      * provides access through the #()#
				      * operator for backcompatibility.
				      */
    double   operator [] (const unsigned int index) const;

    				     /**
				      *  Read and write access to the #index#th
				      *  coordinate.
				      *
				      * Note that the derived #Point# class also
				      * provides access through the #()#
				      * operator for backcompatibility.
				      */
    double & operator [] (const unsigned int index);

				     /**
				      * Reset all values to zero.
				      */
    void clear ();
    
				     /**
				      *  Exception
				      */
    DeclException1 (ExcDimTooSmall,
		    int,
		    << "Given dimensions must be >= 1, but was " << arg1);

				     /**
				      *  Exception
				      */
    DeclException1 (ExcInvalidIndex,
		    int,
		    << "Invalid index " << arg1);
  protected:
				     /**
				      *  Stores the values in a simple array.
				      */
    double values[dim];
};




/*------------------------------- Inline functions: Tensor ---------------------------*/


template <int dim>
inline
Tensor<1,dim>::Tensor () {
  Assert (dim>0, ExcDimTooSmall(dim));

  for (unsigned int i=0; i<dim; ++i)
    values[i] = 0;
};



template <int dim>
inline
Tensor<1,dim>::Tensor (const Tensor<1,dim> &p) {
  for (unsigned int i=0; i<dim; ++i)
    values[i] = p.values[i];
};



template <int dim>
inline
double Tensor<1,dim>::operator [] (const unsigned int index) const {
  Assert (index<dim, ExcInvalidIndex (index));
  return values[index];
};



template <int dim>
inline
double & Tensor<1,dim>::operator [] (const unsigned int index) {
  Assert (index<dim, ExcInvalidIndex (index));
  return values[index];
};



template <int dim>
inline
void Tensor<1,dim>::clear () {
  for (unsigned int i=0; i<dim; ++i)
    values[i] = 0;
};




/*----------------------------   tensor_base.h     ---------------------------*/
/* end of #ifndef __tensor_base_H */
#endif
/*----------------------------   tensor_base.h     ---------------------------*/
