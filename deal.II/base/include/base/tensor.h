/*----------------------------   tensor.h     ---------------------------*/
/*      $Id$                 */
#ifndef __tensor_H
#define __tensor_H
/*----------------------------   tensor.h     ---------------------------*/


#include <base/tensor_base.h>



/**
 * Provide a general tensor class with an arbitrary rank, i.e. with
 * an arbitrary number of indices. The Tensor class provides an
 * indexing operator and a bit of infrastructure, but most
 * functionality is recursively handed down to tensors of rank 1 or
 * put into external templated functions, e.g. the #contract# family.
 *
 * Using this tensor class for objects of rank to has advantages over
 * matrices in many cases since the dimension is known to the compiler
 * as well as the location of the data. It is therefore possible to
 * produce far more efficient code than for matrices with
 * runtime-dependant dimension.
 */
template <int rank_, int dim>
class Tensor {
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
    static const unsigned int rank      = rank_;
    
				     /**
				      * Read-Write access operator.
				      */
    Tensor<rank_-1,dim> &operator [] (const unsigned int i);

				     /**
				      * Read-only access operator.
				      */
    const Tensor<rank_-1,dim> &operator [] (const unsigned int i) const;

				     /**
				      *  Assignment operator.
				      */
    Tensor & operator = (const Tensor<rank_,dim> &);

				     /**
				      *  Test for equality of two points.
				      */
    bool operator == (const Tensor<rank_,dim> &) const;

    				     /**
				      *  Test for inequality of two points.
				      */
    bool operator != (const Tensor<rank_,dim> &) const;

				     /**
				      * Reset all values to zero.
				      */
    void clear ();

    
				     /**
				      *  Exception
				      */
    DeclException1 (ExcInvalidIndex,
		    int,
		    << "Invalid index " << arg1);
  private:
				     /**
				      * Array of tensors holding the
				      * subelements.
				      */
    Tensor<rank_-1,dim> subtensor[dim];
};




/*--------------------------- Inline functions -----------------------------*/

template <int rank_, int dim>
inline
Tensor<rank_-1,dim> &
Tensor<rank_,dim>::operator[] (const unsigned int i) {
  Assert (i<dim, ExcInvalidIndex(i));
  
  return subtensor[i];
};



template <int rank_, int dim>
inline
const Tensor<rank_-1,dim> &
Tensor<rank_,dim>::operator[] (const unsigned int i) const {
  Assert (i<dim, ExcInvalidIndex(i));
  
  return subtensor[i];
};



template <int rank_, int dim>
inline
Tensor<rank_,dim> & Tensor<rank_,dim>::operator = (const Tensor<rank_,dim> &t) {
  for (unsigned int i=0; i<dim; ++i)
    values[i] = t.values[i];
  return *this;
};



template <int rank_, int dim>
inline
bool Tensor<rank_,dim>::operator == (const Tensor<rank_,dim> &p) const {
  for (unsigned int i=0; i<dim; ++i)
    if (values[i] != p.values[i]) return false;
  return true;
};



template <int rank_, int dim>
inline
bool Tensor<rank_,dim>::operator != (const Tensor<rank_,dim> &p) const {
  return !((*this) == p);
};



template <int rank_, int dim>
inline
void Tensor<rank_,dim>::clear () {
  for (unsigned int i=0; i<dim; ++i)
    subtensor[i].clear();
};





/* ----------------- Non-member functions operating on tensors. ------------ */




template <int dim>
inline
void contract (Tensor<2,dim>       &dest,
	       const Tensor<2,dim> &src1,
	       const Tensor<2,dim> &src2) {
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	dest[i][j] += src1[i][k] * src2[k][j];
};



template <int dim>
inline
void contract (Tensor<3,dim>       &dest,
	       const Tensor<3,dim> &src1,
	       const Tensor<2,dim> &src2) {
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  dest[i][j][k] += src1[i][j][l] * src2[l][k];
};



template <int dim>
inline
void contract (Tensor<3,dim>       &dest,
	       const Tensor<2,dim> &src1,
	       const Tensor<3,dim> &src2) {
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  dest[i][j][k] += src1[i][l] * src2[l][j][k];
};



template <int dim>
inline
void contract (Tensor<4,dim>       &dest,
	       const Tensor<3,dim> &src1,
	       const Tensor<3,dim> &src2) {
  dest.clear ();
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      for (unsigned int k=0; k<dim; ++k)
	for (unsigned int l=0; l<dim; ++l)
	  for (unsigned int m=0; m<dim; ++m)
	    dest[i][j][k][l] += src1[i][j][m] * src2[m][k][l];
};

    
    
    
    



/*----------------------------   tensor.h     ---------------------------*/
/* end of #ifndef __tensor_H */
#endif
/*----------------------------   tensor.h     ---------------------------*/
