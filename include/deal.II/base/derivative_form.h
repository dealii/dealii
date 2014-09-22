// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__derivative_form_h
#define __deal2__derivative_form_h

#include <deal.II/base/tensor.h>

DEAL_II_NAMESPACE_OPEN

/**
   This class represents the (tangential) derivatives of a function
   $ f: {\mathbb R}^{\text{dim}} \rightarrow {\mathbb R}^{\text{spacedim}}$.
   Such functions are always used to map the reference dim-dimensional
   cell into spacedim-dimensional space.
   For such objects, the first  derivative of the function is a linear map from
   ${\mathbb R}^{\text{dim}}$  to ${\mathbb R}^{\text{spacedim}}$,
   the second derivative a bilinear map
   from  ${\mathbb R}^{\text{dim}} \times  {\mathbb R}^{\text{dim}}$
   to ${\mathbb R}^{\text{spacedim}}$ and so on.

   In deal.II we represent these derivaties using objects of
   type DerivativeForm<1,dim,spacedim>, DerivativeForm<2,dim,spacedim> and so on.

   @author Sebastian Pauletti, 2011
*/
template <int order, int dim, int spacedim>
class DerivativeForm
{
public:
  /**
   * Constructor. Initialize all entries
   * to zero.
   */
  DerivativeForm ();

  /**
     Constructor from a second order tensor.
  */
  DerivativeForm (const Tensor<2,dim> &);

  /**
   * Read-Write access operator.
   */
  Tensor<order,dim> &operator [] (const unsigned int i);

  /**
   * Read-only access operator.
   */
  const Tensor<order,dim> &operator [] (const unsigned int i) const;

  /**
   *  Assignment operator.
   */
  DerivativeForm   &operator = (const DerivativeForm <order, dim, spacedim> &);


  /**
   *  Assignment operator.
   */
  DerivativeForm   &operator = (const Tensor<2,dim> &);

  /**
   *  Assignment operator.
   */
  DerivativeForm   &operator = (const Tensor<1,dim> &);

  /**
     Converts a DerivativeForm <1,dim, dim>
     to Tensor<2,dim>.
     If the derivative is the Jacobian of F,
     then Tensor[i] = grad(F^i).
  */
  operator Tensor<2,dim>() const;

  /**
     Converts a DerivativeForm <1, dim, 1>
     to Tensor<1,dim>.
  */
  operator Tensor<1,dim>() const;

  /**
     Return the transpose of a rectangular DerivativeForm,
     that is to say viewed as a two dimensional matrix.
  */
  DerivativeForm<1, spacedim, dim> transpose () const;


  /**
     Computes the volume element associated with the
     jacobian of the tranformation F.
     That is to say if $DF$ is square, it computes
     $\det(DF)$, in case DF is not square returns
     $\sqrt(\det(DF^{t} * DF))$.
  */
  double determinant () const;


  /**
     Assuming (*this) stores the jacobian of
     the mapping F, it computes its covariant
     matrix, namely $DF*G^{-1}$, where
     $G = DF^{t}*DF$.
     If $DF$ is square, covariant from
     gives $DF^{-t}$.
  */
  DerivativeForm<1, dim, spacedim> covariant_form() const;




  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  static std::size_t memory_consumption ();

  /**
   * Exception.
   */
  DeclException1 (ExcInvalidTensorIndex,
                  int,
                  << "Invalid DerivativeForm index " << arg1);


private:
  /** Auxiliary function that computes
      (*this) * T^{t} */
  DerivativeForm<1, dim, spacedim> times_T_t (Tensor<2,dim> T) const;


private:
  /**
   * Array of tensors holding the
   * subelements.
   */
  Tensor<order,dim> tensor[spacedim];


};


/*--------------------------- Inline functions -----------------------------*/

#ifndef DOXYGEN




template <int order, int dim, int spacedim>
inline
DerivativeForm<order, dim, spacedim>::DerivativeForm  ()
{
// default constructor. not specifying an initializer list calls
// the default constructor of the subobjects, which initialize them
// selves. therefore, the tensor array  is set to zero this way
}



template <int order, int dim, int spacedim>
inline
DerivativeForm<order, dim, spacedim>::DerivativeForm(const Tensor<2,dim> &T)
{
  Assert( (dim == spacedim) && (order==1),
          ExcMessage("Only allowed for square tensors."));
  if ((dim == spacedim) && (order==1))
    for (unsigned int j=0; j<dim; ++j)
      (*this)[j] = T[j];
}



template <int order, int dim, int spacedim>
inline
DerivativeForm<order, dim, spacedim> &
DerivativeForm<order, dim, spacedim>::
operator = (const DerivativeForm<order, dim, spacedim> &ta)
{
  for (unsigned int j=0; j<spacedim; ++j)
    (*this)[j] = ta[j];
  return *this;
}



template <int order, int dim, int spacedim>
inline
DerivativeForm<order, dim, spacedim> &DerivativeForm<order, dim, spacedim>::
operator = (const Tensor<2,dim> &ta)
{
  Assert( (dim == spacedim) && (order==1),
          ExcMessage("Only allowed for square tensors."));

  if ((dim == spacedim) && (order==1))
    for (unsigned int j=0; j<dim; ++j)
      (*this)[j] = ta[j];
  return *this;

}



template <int order, int dim, int spacedim>
inline
DerivativeForm<order, dim, spacedim> &DerivativeForm<order, dim, spacedim>::
operator = (const Tensor<1,dim> &T)
{
  Assert( (1 == spacedim) && (order==1),
          ExcMessage("Only allowed for spacedim==1 and order==1."));

  (*this)[0] = T;

  return *this;

}



template <int order, int dim, int spacedim>
inline
Tensor<order,dim> &DerivativeForm<order, dim, spacedim>::
operator[] (const unsigned int i)
{
  Assert (i<spacedim, ExcIndexRange(i, 0, spacedim));

  return tensor[i];
}



template <int order, int dim, int spacedim>
inline
const Tensor<order,dim> &DerivativeForm<order, dim, spacedim>::
operator[] (const unsigned int i) const
{
  Assert (i<spacedim, ExcIndexRange(i, 0, spacedim));

  return tensor[i];
}



template <int order, int dim, int spacedim>
inline
DerivativeForm<order, dim, spacedim>::operator Tensor<1,dim>() const
{
  Assert( (1 == spacedim) && (order==1),
          ExcMessage("Only allowed for spacedim==1."));

  return (*this)[0];

}



template <int order, int dim, int spacedim>
inline
DerivativeForm<order, dim, spacedim>::operator Tensor<2,dim>() const
{
  Assert( (dim == spacedim) && (order==1),
          ExcMessage("Only allowed for square tensors."));



  Tensor<2,dim> t;

  if ((dim == spacedim) && (order==1))
    for (unsigned int j=0; j<dim; ++j)
      t[j] = (*this)[j];

  return t;

}



template <int order, int dim, int spacedim>
inline
DerivativeForm<1,spacedim,dim>
DerivativeForm<order,dim,spacedim>::
transpose () const
{
  Assert(order==1, ExcMessage("Only for rectangular DerivativeForm."));
  DerivativeForm<1,spacedim,dim> tt;

  for (unsigned int i=0; i<spacedim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      tt[j][i] = (*this)[i][j];

  return tt;
}



template <int order, int dim, int spacedim>
inline
DerivativeForm<1, dim, spacedim>
DerivativeForm<order,dim,spacedim>::times_T_t (Tensor<2,dim> T) const
{
  Assert( order==1, ExcMessage("Only for order == 1."));
  DerivativeForm<1,dim, spacedim> dest;
  for (unsigned int i=0; i<spacedim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      dest[i][j] = (*this)[i] * T[j];

  return dest;
}


template <int order, int dim, int spacedim>
inline
double
DerivativeForm<order,dim,spacedim>::determinant () const
{
  Assert( order==1, ExcMessage("Only for order == 1."));
  if (dim == spacedim)
    {
      Tensor<2,dim> T = (Tensor<2,dim>)( (*this) );
      return dealii::determinant(T);
    }
  else
    {
      Assert( spacedim>dim, ExcMessage("Only for spacedim>dim."));
      DerivativeForm<1,spacedim,dim> DF_t = this->transpose();
      Tensor<2, dim> G; //First fundamental form
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          G[i][j] = DF_t[i] * DF_t[j];

      return ( sqrt(dealii::determinant(G)) );

    }

}



template <int order, int dim, int spacedim>
inline
DerivativeForm<1,dim,spacedim>
DerivativeForm<order,dim,spacedim>::covariant_form() const
{

  if (dim == spacedim)
    {

      Tensor<2,dim> DF_t (dealii::transpose(invert(  (Tensor<2,dim>)(*this)   )));
      DerivativeForm<1,dim, spacedim> result = DF_t;
      return (result);
    }
  else
    {

      DerivativeForm<1,spacedim,dim> DF_t = this->transpose();
      Tensor<2, dim> G; //First fundamental form
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          G[i][j] = DF_t[i] * DF_t[j];

      return (this->times_T_t(invert(G)));

    }

}


template <int order, int dim, int spacedim>
inline
std::size_t
DerivativeForm<order, dim, spacedim>::memory_consumption ()
{
  return sizeof(DerivativeForm<order, dim, spacedim>);
}

#endif // DOXYGEN





/**
   One of the uses of DerivativeForm is to apply it as a transformation.
   This is what this function does.

   If @p T is DerivativeForm<1,dim,1> it computes $DF * T$,
   if @p T is DerivativeForm<1,dim,rank> it computes $T*DF^{t}$.

   @relates DerivativeForm
   @author Sebastian Pauletti, 2011
*/
template <int spacedim, int dim>
inline
Tensor<1, spacedim>
apply_transformation (const DerivativeForm<1,dim,spacedim> &DF,
                      const Tensor<1,dim>               &T)
{
  Tensor<1, spacedim> dest;
  for (unsigned int i=0; i<spacedim; ++i)
    dest[i] = DF[i] * T;
  return dest;
}



/**
   Similar to previous apply_transformation.
   It computes $T*DF^{t}$
   @relates DerivativeForm
   @author Sebastian Pauletti, 2011
*/
//rank=2
template <int spacedim, int dim>
inline
DerivativeForm<1, spacedim, dim>
apply_transformation (const DerivativeForm<1,dim,spacedim> &DF,
                      const Tensor<2,dim>               &T)
{

  DerivativeForm<1, spacedim, dim> dest;
  for (unsigned int i=0; i<dim; ++i)
    dest[i] = apply_transformation(DF, T[i]);

  return dest;
}

/**
   Similar to previous apply_transformation.
   It computes $DF2*DF1^{t}$
   @relates DerivativeForm
   @author Sebastian Pauletti, 2011
*/
template <int spacedim, int dim>
inline
Tensor<2, spacedim>
apply_transformation (const DerivativeForm<1,dim,spacedim> &DF1,
                      const DerivativeForm<1,dim,spacedim> &DF2)
{
  Tensor<2, spacedim> dest;

  for (unsigned int i=0; i<spacedim; ++i)
    dest[i] = apply_transformation(DF1, DF2[i]);

  return dest;
}


/**
   Transpose of a rectangular DerivativeForm DF,
   mostly for compatibility reasons.
   @relates DerivativeForm
   @author Sebastian Pauletti, 2011
*/
template <int dim, int spacedim>
inline
DerivativeForm<1,spacedim,dim>
transpose (const DerivativeForm<1,dim,spacedim> &DF)
{
  DerivativeForm<1,spacedim,dim> tt;
  tt = DF.transpose();
  return tt;
}


DEAL_II_NAMESPACE_CLOSE

#endif
