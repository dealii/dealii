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


class FE_P1NC : public FiniteElement<2,2>
{

public:
  /**
   * Constructor for nonparametric version of P1 nonconforming element.
   */
  FE_P1NC() ;

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
  virtual ~FE_P1NC ();



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
   * Compute the linear shape functions phi(x,y) = ax + by + c
   * such that each midpoint value on two connecting edges is a half,
   * and two other midpoint values are all zero.
   */
  static
  void
  get_linear_shape (const Triangulation<2,2>::cell_iterator &cell,
                    std::vector<double> &a,
                    std::vector<double> &b,
                    std::vector<double> &c);



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
                  const Mapping<2,2>::InternalDataBase &,
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



  /**
   * Create the constraints matrix for hanging edges.
   */
  void initialize_constraints () ;

};




/** @}*/

DEAL_II_NAMESPACE_CLOSE

#endif
