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
#include <iostream>




// general template; specialized for rank==1; the general template is in
// tensor.h
template <int rank, int dim> class Tensor;


/**
 * This class is a specialized version of the #Tensor<rank,dim># class.
 * It handles tensors with one index, i.e. vectors, of fixed dimension
 * and offers the functionality needed for tensors of higher rank.
 *
 * In many cases, you will want to use the more specialized #Point# class
 * which acts as a tensor of rank one but has more functionality.
 */
template <int dim>
class Tensor<1,dim> {
      public:
				     /**
				      * Provide a way to get the
				      * dimension of an object without
				      * explicit knowledge of it's
				      * data type. Implementation is
				      * this way instead of providing
				      * a function #dimension()#
				      * because now it is possible to
				      * get the dimension at compile
				      * time without the expansion and
				      * preevaluation of an inlined
				      * function; the compiler may
				      * therefore produce more
				      * efficient code and you may use
				      * this value to declare other
				      * data types.
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
				      *  Assignment operator.
				      */
    Tensor<1,dim> & operator = (const Tensor<1,dim> &);

				     /**
				      *  Test for equality of two points.
				      */
    bool operator == (const Tensor<1,dim> &) const;

    				     /**
				      *  Test for inequality of two points.
				      */
    bool operator != (const Tensor<1,dim> &) const;

				     /**
				      *  Add another vector, i.e. move this point by
				      *  the given offset.
				      */
    Tensor<1,dim> & operator += (const Tensor<1,dim> &);
				     /**
				      *  Subtract another vector.
				      */
    Tensor<1,dim> & operator -= (const Tensor<1,dim> &);

				     /**
				      *  Scale the vector by #factor#, i.e. multiply
				      *  all coordinates by #factor#.
				      */
    Tensor<1,dim> & operator *= (const double &factor);

				     /**
				      *  Scale the vector by #1/factor#.
				      */
    Tensor<1,dim> & operator /= (const double &factor);

				     /**
				      *  Returns the scalar product of two vectors.
				      */
    double              operator * (const Tensor<1,dim> &) const;

				     /**
				      * Reset all values to zero.
				      */
    void clear ();
    
  protected:
				     /**
				      *  Stores the values in a simple array.
				      */
    double values[dim];
};




DeclException1 (ExcDimTooSmall,
		int,
		<< "Given dimensions must be >= 1, but was " << arg1);
DeclException1 (ExcInvalidIndex,
		int,
		<< "Invalid index " << arg1);




				 /**
				  *  Prints the values of this point in the
				  *  form #x1 x2 x3 etc#.
				  */
template <int dim>
ostream & operator << (ostream &out, const Tensor<1,dim> &p);






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
Tensor<1,dim> & Tensor<1,dim>::operator = (const Tensor<1,dim> &p) {
  for (unsigned int i=0; i<dim; ++i)
    values[i] = p.values[i];
  return *this;
};



template <int dim>
inline
bool Tensor<1,dim>::operator == (const Tensor<1,dim> &p) const {
  for (unsigned int i=0; i<dim; ++i)
    if (values[i] != p.values[i]) return false;
  return true;
};



template <int dim>
inline
bool Tensor<1,dim>::operator != (const Tensor<1,dim> &p) const {
  return !((*this) == p);
};



template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator += (const Tensor<1,dim> &p) {
  for (unsigned int i=0; i<dim; ++i)
    values[i] += p.values[i];
  return *this;
};



template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator -= (const Tensor<1,dim> &p) {
  for (unsigned int i=0; i<dim; ++i)
    values[i] -= p.values[i];
  return *this;
};



template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator *= (const double &s) {
  for (unsigned int i=0; i<dim; ++i)
    values[i] *= s;
  return *this;
};



template <int dim>
inline
Tensor<1,dim> & Tensor<1,dim>::operator /= (const double &s) {
  for (unsigned int i=0; i<dim; ++i)
    values[i] /= s;
  return *this;
};



template <int dim>
inline
double Tensor<1,dim>::operator * (const Tensor<1,dim> &p) const {
  double q=0;
  for (unsigned int i=0; i<dim; ++i)
    q += values[i] * p.values[i];
  return q;
};



template <int dim>
inline
void Tensor<1,dim>::clear () {
  for (unsigned int i=0; i<dim; ++i)
    values[i] = 0;
};






template <int dim>
inline
ostream & operator << (ostream &out, const Tensor<1,dim> &p) {
  for (unsigned int i=0; i<dim-1; ++i)
    out << p[i] << ' ';
  out << p[dim-1];

  if (!out)
    throw GlobalExcIO ();

  return out;
};


template <>
inline
ostream & operator << (ostream &out, const Tensor<1,1> &p) {
  out << p[0];

  if (!out)
    throw GlobalExcIO ();

  return out;
};
  





/*----------------------------   tensor_base.h     ---------------------------*/
/* end of #ifndef __tensor_base_H */
#endif
/*----------------------------   tensor_base.h     ---------------------------*/
