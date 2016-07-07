// ---------------------------------------------------------------------
//
// Copyright (C) 2015 by the deal.II authors
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

#ifndef dealii__fe_p1nc_h
#define dealii__fe_p1nc_h

#include <deal.II/base/config.h>
#include <deal.II/fe/fe.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>


DEAL_II_NAMESPACE_OPEN


/*!@addtogroup fe */
/*@{*/


// parametric version of P1 Nonconforming FE
class FE_P1NCParametric : public FiniteElement<2,2>
{
private:


  static
  std::vector<ComponentMask>
  get_nonzero_component();


  /**
   * Return the vector consists of the numbers of degrees of freedom per objects.
   */
  static
  std::vector<unsigned int>
  get_dpo_vector ();


  class InternalData : public FiniteElement<2,2>::InternalDataBase
  {
  public:

    Table<2,double > shape_values;
    Table<2,Tensor<1,2> > shape_gradients;
    mutable std::vector<Tensor<1,2> > transformed_shape_gradients ;
    Table<2,Tensor<2,2> > shape_hessians;
  };

  /**
   * Compute the data on the reference domain.
   */
virtual FiniteElement<2,2>::InternalDataBase *
  get_data (const UpdateFlags update_flags,
            const Mapping<2,2> &,
            const Quadrature<2> &,
            dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const ;




/**
 * Compute the data on the current cell.
 */
  virtual
  void
  fill_fe_values (const Triangulation<2,2>::cell_iterator           &cell,
                  const CellSimilarity::Similarity                   ,
                  const Quadrature<2>                               &quadrature,
                  const Mapping<2,2>                                &mapping,
                  const Mapping<2,2>::InternalDataBase              &,
                  const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                  const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                  internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const;




  /**
   * Compute the data on the face of the current cell.
   */
  virtual
  void
  fill_fe_face_values (const Triangulation<2,2>::cell_iterator           &cell,
                       const unsigned int                                                   face_no,
                       const Quadrature<1>                                             &quadrature,
                       const Mapping<2,2>                                         &mapping,
                       const Mapping<2,2>::InternalDataBase              &mapping_internal,
                       const dealii::internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                       const InternalDataBase                                              &fe_internal,
                       dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const;




  /**
   * Compute the data on the subface of the current cell.
   */
  virtual
  void
  fill_fe_subface_values (const Triangulation<2,2>::cell_iterator           &cell,
                          const unsigned int                                                   face_no,
                          const unsigned int                                                   sub_no,
                          const Quadrature<1>                                             &quadrature,
                          const Mapping<2,2>                                         &mapping,
                          const Mapping<2,2>::InternalDataBase              &mapping_internal,
                          const dealii::internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                          const InternalDataBase                                              &fe_internal,
                          dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const;



public:
  /**
   * Constructor for parametric version of P1 nonconforming element.
   */
  FE_P1NCParametric() ;


  /**
   * Return the name of the class for the element.
   */
  virtual std::string get_name () const ;

  /**
   * Return the update flags which are needed.
   */
  virtual UpdateFlags     requires_update_flags (const UpdateFlags flags) const ;

  /**
   * Return the unit face support points.
   */
  const std::vector<Point<1> > &get_unit_face_support_points () const ;

  /**
   * Return the unit face support point with given index.
   */
  Point<1> unit_face_support_point (const unsigned index) const ;

  bool has_face_support_points () const ;

  /**
   * Copy constructor.
   */
  virtual FiniteElement<2,2> *clone () const ;

  /**
   * Destructor.
   */
  virtual ~FE_P1NCParametric ();




  /**
   * The value of the basis function with index i at given point.
   */
  virtual double shape_value (const unsigned int i,
                              const Point<2>    &p) const ;

  /**
   * The gradient of the basis function with index i at given point.
   */
  virtual Tensor<1,2> shape_grad (const unsigned int i,
                                  const Point<2>    &p) const ;

  /**
   * The Hessian matrix of the basis function with index i at given point.
   * Notice that it always returns the trivial matrix because the element space contains polynomials upto the degree 1.
   */
  virtual Tensor<2,2> shape_grad_grad (const unsigned int i,
				       const Point<2>    &p) const ;


protected:

  /**
   * Create the constraints matrix for hanging edges.
   */
  void initialize_constraints () ;

};










// Nonparametric version of P1 Nonconforming FE
class FE_P1NCNonparametric : public FiniteElement<2,2>
{
private:


  static
  std::vector<ComponentMask>
  get_nonzero_component();


  /**
   * Return the vector consists of the numbers of degrees of freedom per objects.
   */
  static
  std::vector<unsigned int>
  get_dpo_vector ();



  /**
   * Do the work which is needed before cellwise data computation.
   * Since the basis functions are constructed independently on each cell,
   * the data on the reference cell is not necessary.
   */
  virtual FiniteElement<2,2>::InternalDataBase *
  get_data (const UpdateFlags update_flags,
            const Mapping<2,2> &,
            const Quadrature<2> &,
            dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const ;




  /**
   * Compute the data on the current cell.
   */
   virtual
  void
  fill_fe_values (const Triangulation<2,2>::cell_iterator           &cell,
                  const CellSimilarity::Similarity                   ,
                  const Quadrature<2>                               &quadrature,
                  const Mapping<2,2>                                &mapping,
                  const Mapping<2,2>::InternalDataBase              &,
                  const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                  const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                  internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const;



  /**
   * Compute the data on the face of the current cell.
   */
  virtual
  void
  fill_fe_face_values (const Triangulation<2,2>::cell_iterator           &cell,
                       const unsigned int                                                   face_no,
                       const Quadrature<1>                                             &quadrature,
                       const Mapping<2,2>                                         &mapping,
                       const Mapping<2,2>::InternalDataBase              &mapping_internal,
                       const dealii::internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                       const InternalDataBase                                              &fe_internal,
                       dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const;




  /**
   * Compute the data on the subface of the current cell.
   */
  virtual
  void
  fill_fe_subface_values (const Triangulation<2,2>::cell_iterator           &cell,
                          const unsigned int                                                   face_no,
                          const unsigned int                                                   sub_no,
                          const Quadrature<1>                                             &quadrature,
                          const Mapping<2,2>                                         &mapping,
                          const Mapping<2,2>::InternalDataBase              &mapping_internal,
                          const dealii::internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                          const InternalDataBase                                              &fe_internal,
                          dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const;



public:
  /**
   * Constructor for nonparametric version of P1 nonconforming element.
   */
  FE_P1NCNonparametric() ;

  /**
   * Return the name of the class for the element.
   */
  virtual std::string get_name () const ;

  /**
   * Return the update flags which are needed.
   */
  virtual UpdateFlags     requires_update_flags (const UpdateFlags flags) const ;

  /**
   * Copy constructor.
   */
  virtual FiniteElement<2,2> *clone () const ;

  /**
   * Destructor.
   */
  virtual ~FE_P1NCNonparametric ();


protected:

  /**
   * Create the constraints matrix for hanging edges.
   */
  void initialize_constraints () ;

};




/** @}*/

DEAL_II_NAMESPACE_CLOSE

#endif
