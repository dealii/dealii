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


  static
  std::vector<unsigned int>
  get_dpo_vector ();


  class InternalData : public FiniteElement<2,2>::InternalDataBase
  {
  public:

    Table<2,double > shape_values;
    Table<2,Tensor<1,2> > shape_gradients;
    mutable std::vector<Tensor<1,2> > transformed_shape_gradients ;
  };


virtual FiniteElement<2,2>::InternalDataBase *
  get_data (const UpdateFlags update_flags,
            const Mapping<2,2> &,
            const Quadrature<2> &,
            dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const ;




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
  FE_P1NCParametric() ;



  virtual std::string get_name () const ;

  virtual UpdateFlags     requires_update_flags (const UpdateFlags flags) const ;

  const std::vector<Point<1> > &get_unit_face_support_points () const ;

  Point<1> unit_face_support_point (const unsigned index) const ;

  bool has_face_support_points () const ;

  virtual FiniteElement<2,2> *clone () const ;

  virtual ~FE_P1NCParametric ();





  virtual double shape_value (const unsigned int i,
                              const Point<2>    &p) const ;

  virtual Tensor<1,2> shape_grad (const unsigned int i,
                                  const Point<2>    &p) const ;



protected:

  void initialize_constraints () ;

};










// Nonparametric version of P1 Nonconforming FE
class FE_P1NCNonparametric : public FiniteElement<2,2>
{
private:


  static
  std::vector<ComponentMask>
  get_nonzero_component();


  static
  std::vector<unsigned int>
  get_dpo_vector ();



  virtual FiniteElement<2,2>::InternalDataBase *
  get_data (const UpdateFlags update_flags,
            const Mapping<2,2> &,
            const Quadrature<2> &,
            dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const ;




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
  FE_P1NCNonparametric() ;

  virtual std::string get_name () const ;

  virtual UpdateFlags     requires_update_flags (const UpdateFlags flags) const ;

  virtual FiniteElement<2,2> *clone () const ;

  virtual ~FE_P1NCNonparametric ();


protected:

  void initialize_constraints () ;

};




/** @}*/

DEAL_II_NAMESPACE_CLOSE

#endif
