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


#include <deal.II/fe/fe_p1nc.h>

DEAL_II_NAMESPACE_OPEN


std::vector<ComponentMask>
FE_P1NC::get_nonzero_component()
{
  const unsigned int dofs_per_cell = 4;
  std::vector<ComponentMask> masks (dofs_per_cell);
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    masks[i] = ComponentMask(1,true);
  return masks;
}


std::vector<unsigned int>
FE_P1NC::get_dpo_vector ()
{
  std::vector<unsigned int> dpo(3);
  dpo[0] = 1; // dofs per object: vertex
  dpo[1] = 0; // line
  dpo[2] = 0; // quad
  return dpo;
}


FiniteElement<2,2>::InternalDataBase *
FE_P1NC::get_data (const UpdateFlags update_flags,
                                const Mapping<2,2> &,
                                const Quadrature<2> &,
				dealii::internal::FEValues::FiniteElementRelatedData<2,2> &) const
{
  FiniteElement<2,2>::InternalDataBase *data = new FiniteElement<2,2>::InternalDataBase;

  data->update_each = requires_update_flags(update_flags);

  return data;




}



void
FE_P1NC::fill_fe_values (const Triangulation<2,2>::cell_iterator           &cell,
                                      const CellSimilarity::Similarity                   ,
                                      const Quadrature<2>                               &quadrature,
                                      const Mapping<2,2>                                &mapping,
                                      const Mapping<2,2>::InternalDataBase              &,
                                      const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                                      const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                                      internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const
{
  const UpdateFlags flags(fe_internal.update_each) ;

  const unsigned int n_q_points = mapping_data.quadrature_points.size();

  std::vector<double> values(flags & update_values ? this->dofs_per_cell : 0);
  std::vector<Tensor<1,2> > grads(flags & update_gradients ? this->dofs_per_cell : 0);



  // edge midpoints
  std::vector<Point<2> > mpt(4) ;

  mpt[0](0) = (cell->vertex(0)(0) + cell->vertex(2)(0))/2.0 ;
  mpt[0](1) = (cell->vertex(0)(1) + cell->vertex(2)(1))/2.0 ;

  mpt[1](0) = (cell->vertex(1)(0) + cell->vertex(3)(0))/2.0 ;
  mpt[1](1) = (cell->vertex(1)(1) + cell->vertex(3)(1))/2.0 ;

  mpt[2](0) = (cell->vertex(0)(0) + cell->vertex(1)(0))/2.0 ;
  mpt[2](1) = (cell->vertex(0)(1) + cell->vertex(1)(1))/2.0 ;

  mpt[3](0) = (cell->vertex(2)(0) + cell->vertex(3)(0))/2.0 ;
  mpt[3](1) = (cell->vertex(2)(1) + cell->vertex(3)(1))/2.0 ;

  // center point
  Point<2> cpt ;
  cpt(0) = (mpt[0](0) + mpt[1](0) + mpt[2](0) + mpt[3](0))/4.0 ;
  cpt(1) = (mpt[0](1) + mpt[1](1) + mpt[2](1) + mpt[3](1))/4.0 ;


  // basis functions with a half value: phi(x,y) = ax + by + c
  std::vector<double> a(4), b(4), c(4) ;
  double det ;

  det = (mpt[0](0)-mpt[1](0))*(mpt[2](1)-mpt[3](1)) - (mpt[2](0)-mpt[3](0))*(mpt[0](1)-mpt[1](1)) ;

  a[0] = ((mpt[2](1)-mpt[3](1))*(0.5) -(mpt[0](1)-mpt[1](1))*(0.5))/det ;
  a[1] = ((mpt[2](1)-mpt[3](1))*(-0.5) -(mpt[0](1)-mpt[1](1))*(0.5))/det ;
  a[2] = ((mpt[2](1)-mpt[3](1))*(0.5) -(mpt[0](1)-mpt[1](1))*(-0.5))/det ;
  a[3] = ((mpt[2](1)-mpt[3](1))*(-0.5) -(mpt[0](1)-mpt[1](1))*(-0.5))/det ;

  b[0] = (-(mpt[2](0)-mpt[3](0))*(0.5) +(mpt[0](0)-mpt[1](0))*(0.5))/det ;
  b[1] = (-(mpt[2](0)-mpt[3](0))*(-0.5) +(mpt[0](0)-mpt[1](0))*(0.5))/det ;
  b[2] = (-(mpt[2](0)-mpt[3](0))*(0.5) +(mpt[0](0)-mpt[1](0))*(-0.5))/det ;
  b[3] = (-(mpt[2](0)-mpt[3](0))*(-0.5) +(mpt[0](0)-mpt[1](0))*(-0.5))/det ;

  c[0] = 0.25 - cpt(0)*a[0] - cpt(1)*b[0] ;
  c[1] = 0.25 - cpt(0)*a[1] - cpt(1)*b[1] ;
  c[2] = 0.25 - cpt(0)*a[2] - cpt(1)*b[2] ;
  c[3] = 0.25 - cpt(0)*a[3] - cpt(1)*b[3] ;



  // compute basis functions
  if (flags & (update_values | update_gradients))
    for (unsigned int i=0; i<n_q_points; ++i)
      {
        for (unsigned int k=0; k<this->dofs_per_cell; ++k)
          {
            if (flags & update_values)
              {
                values[k] = a[k]*mapping_data.quadrature_points[i](0) + b[k]*mapping_data.quadrature_points[i](1) + c[k] ;
                output_data.shape_values[k][i] = values[k];
              }

            if (flags & update_gradients)
              {
                grads[k][0] = a[k] ;
                grads[k][1] = b[k] ;
                output_data.shape_gradients[k][i] = grads[k];
              }
          }

      }

      // hessian
      std::vector<Tensor<2,2> > hessians(flags & update_hessians ? this->dofs_per_cell : 0);
      if (flags & update_hessians)
      {
        for (unsigned int i=0; i<n_q_points; ++i)
        {
            for (unsigned int k=0; k<this->dofs_per_cell; ++k)
            {
                hessians[k][0][0] = 0.0 ;
                hessians[k][0][1] = 0.0 ;
                hessians[k][1][0] = 0.0 ;
                hessians[k][1][1] = 0.0 ;
                output_data.shape_hessians[k][i] = hessians[k];
            }
        }
      }

  // When this function is called for the graphic out,
  // MappingRelatedData does not work properly.
  // In this case, the quadrature points on the real cell is computed in manual sense, using 'mapping' and 'quadrature'.
  // This is a temporary solution. It needs to be fixed fundamentally.
  if (n_q_points==0)
    {
      for (unsigned int i=0; i<quadrature.size(); ++i)
        for (unsigned int k=0; k<this->dofs_per_cell; ++k)
          {
            Point<2> realquadrature ;

            realquadrature = mapping.transform_unit_to_real_cell(cell, quadrature.point(i)) ;
            values[k] = a[k]*realquadrature(0) + b[k]*realquadrature(1) + c[k] ;
            output_data.shape_values[k][i] = values[k];

          }
    }

}


void
FE_P1NC::fill_fe_face_values (const Triangulation<2,2>::cell_iterator           &cell,
                                           const unsigned int                                 face_no,
                                           const Quadrature<1>                                             &quadrature,
                                           const Mapping<2,2>                                         &mapping,
                                           const Mapping<2,2>::InternalDataBase              &mapping_internal,
                                           const dealii::internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                                           const InternalDataBase                                              &fe_internal,
                                           dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const

{
  const UpdateFlags flags(fe_internal.update_each) ;

  const unsigned int n_q_points = mapping_data.quadrature_points.size();

  std::vector<double> values(flags & update_values ? this->dofs_per_cell : 0);
  std::vector<Tensor<1,2> > grads(flags & update_gradients ? this->dofs_per_cell : 0);



  // edge midpoints
  std::vector<Point<2> > mpt(4) ;

  mpt[0](0) = (cell->vertex(0)(0) + cell->vertex(2)(0))/2.0 ;
  mpt[0](1) = (cell->vertex(0)(1) + cell->vertex(2)(1))/2.0 ;

  mpt[1](0) = (cell->vertex(1)(0) + cell->vertex(3)(0))/2.0 ;
  mpt[1](1) = (cell->vertex(1)(1) + cell->vertex(3)(1))/2.0 ;

  mpt[2](0) = (cell->vertex(0)(0) + cell->vertex(1)(0))/2.0 ;
  mpt[2](1) = (cell->vertex(0)(1) + cell->vertex(1)(1))/2.0 ;

  mpt[3](0) = (cell->vertex(2)(0) + cell->vertex(3)(0))/2.0 ;
  mpt[3](1) = (cell->vertex(2)(1) + cell->vertex(3)(1))/2.0 ;

  // center point
  Point<2> cpt ;
  cpt(0) = (mpt[0](0) + mpt[1](0) + mpt[2](0) + mpt[3](0))/4.0 ;
  cpt(1) = (mpt[0](1) + mpt[1](1) + mpt[2](1) + mpt[3](1))/4.0 ;


  // basis functions with a half value: phi(x,y) = ax + by + c
  std::vector<double> a(4), b(4), c(4) ;
  double det ;

  det = (mpt[0](0)-mpt[1](0))*(mpt[2](1)-mpt[3](1)) - (mpt[2](0)-mpt[3](0))*(mpt[0](1)-mpt[1](1)) ;

  a[0] = ((mpt[2](1)-mpt[3](1))*(0.5) -(mpt[0](1)-mpt[1](1))*(0.5))/det ;
  a[1] = ((mpt[2](1)-mpt[3](1))*(-0.5) -(mpt[0](1)-mpt[1](1))*(0.5))/det ;
  a[2] = ((mpt[2](1)-mpt[3](1))*(0.5) -(mpt[0](1)-mpt[1](1))*(-0.5))/det ;
  a[3] = ((mpt[2](1)-mpt[3](1))*(-0.5) -(mpt[0](1)-mpt[1](1))*(-0.5))/det ;

  b[0] = (-(mpt[2](0)-mpt[3](0))*(0.5) +(mpt[0](0)-mpt[1](0))*(0.5))/det ;
  b[1] = (-(mpt[2](0)-mpt[3](0))*(-0.5) +(mpt[0](0)-mpt[1](0))*(0.5))/det ;
  b[2] = (-(mpt[2](0)-mpt[3](0))*(0.5) +(mpt[0](0)-mpt[1](0))*(-0.5))/det ;
  b[3] = (-(mpt[2](0)-mpt[3](0))*(-0.5) +(mpt[0](0)-mpt[1](0))*(-0.5))/det ;

  c[0] = 0.25 - cpt(0)*a[0] - cpt(1)*b[0] ;
  c[1] = 0.25 - cpt(0)*a[1] - cpt(1)*b[1] ;
  c[2] = 0.25 - cpt(0)*a[2] - cpt(1)*b[2] ;
  c[3] = 0.25 - cpt(0)*a[3] - cpt(1)*b[3] ;



  // compute basis functions
  if (flags & (update_values | update_gradients))
    for (unsigned int i=0; i<n_q_points; ++i)
      {
        for (unsigned int k=0; k<this->dofs_per_cell; ++k)
          {
            if (flags & update_values)
              {
                values[k] = a[k]*mapping_data.quadrature_points[i](0) + b[k]*mapping_data.quadrature_points[i](1) + c[k] ;
                output_data.shape_values[k][i] = values[k];
              }

            if (flags & update_gradients)
              {
                grads[k][0] = a[k] ;
                grads[k][1] = b[k] ;
                output_data.shape_gradients[k][i] = grads[k];
              }
          }

      }

  // When this function is called for computation of facial jump residual,
  // MappingRelatedData does not work properly.
  // In this case, the quadrature points on the real cell is computed in manual sense, using 'mapping' and 'quadrature'.
  // This is a temporary solution. It needs to be fixed fundamentally.
  if (n_q_points==0)
    {
      Quadrature<2> cellquadrature = QProjector<2>::project_to_face(quadrature, face_no) ;
      for(unsigned int i=0;i<cellquadrature.size();++i)
               for (unsigned int k=0; k<this->dofs_per_cell; ++k)
               {
                   if (flags & update_values)
                   {
                   Point<2> realquadrature ;

                   realquadrature = mapping.transform_unit_to_real_cell(cell, cellquadrature.point(i)) ;
                   values[k] = a[k]*realquadrature(0) + b[k]*realquadrature(1) + c[k] ;
                   output_data.shape_values[k][i] = values[k];
                   }

                   if (flags & update_gradients)
                   {
                   grads[k][0] = a[k] ;
                   grads[k][1] = b[k] ;
                   output_data.shape_gradients[k][i] = grads[k];
                   }

               }

    }

}






void
FE_P1NC::fill_fe_subface_values (const Triangulation<2,2>::cell_iterator           &cell,
                                              const unsigned int                                                  face_no,
                                              const unsigned int                                                   sub_no,
                                              const Quadrature<1>                                             &quadrature,
                                              const Mapping<2,2>                                         &mapping,
                                              const Mapping<2,2>::InternalDataBase              &mapping_internal,
                                              const dealii::internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                                              const InternalDataBase                                              &fe_internal,
                                              dealii::internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const

{
  const UpdateFlags flags(fe_internal.update_each) ;

  const unsigned int n_q_points = mapping_data.quadrature_points.size();

  std::vector<double> values(flags & update_values ? this->dofs_per_cell : 0);
  std::vector<Tensor<1,2> > grads(flags & update_gradients ? this->dofs_per_cell : 0);



  // edge midpoints
  std::vector<Point<2> > mpt(4) ;

  mpt[0](0) = (cell->vertex(0)(0) + cell->vertex(2)(0))/2.0 ;
  mpt[0](1) = (cell->vertex(0)(1) + cell->vertex(2)(1))/2.0 ;

  mpt[1](0) = (cell->vertex(1)(0) + cell->vertex(3)(0))/2.0 ;
  mpt[1](1) = (cell->vertex(1)(1) + cell->vertex(3)(1))/2.0 ;

  mpt[2](0) = (cell->vertex(0)(0) + cell->vertex(1)(0))/2.0 ;
  mpt[2](1) = (cell->vertex(0)(1) + cell->vertex(1)(1))/2.0 ;

  mpt[3](0) = (cell->vertex(2)(0) + cell->vertex(3)(0))/2.0 ;
  mpt[3](1) = (cell->vertex(2)(1) + cell->vertex(3)(1))/2.0 ;

  // center point
  Point<2> cpt ;
  cpt(0) = (mpt[0](0) + mpt[1](0) + mpt[2](0) + mpt[3](0))/4.0 ;
  cpt(1) = (mpt[0](1) + mpt[1](1) + mpt[2](1) + mpt[3](1))/4.0 ;


  // basis functions with a half value: phi(x,y) = ax + by + c
  std::vector<double> a(4), b(4), c(4) ;
  double det ;

  det = (mpt[0](0)-mpt[1](0))*(mpt[2](1)-mpt[3](1)) - (mpt[2](0)-mpt[3](0))*(mpt[0](1)-mpt[1](1)) ;

  a[0] = ((mpt[2](1)-mpt[3](1))*(0.5) -(mpt[0](1)-mpt[1](1))*(0.5))/det ;
  a[1] = ((mpt[2](1)-mpt[3](1))*(-0.5) -(mpt[0](1)-mpt[1](1))*(0.5))/det ;
  a[2] = ((mpt[2](1)-mpt[3](1))*(0.5) -(mpt[0](1)-mpt[1](1))*(-0.5))/det ;
  a[3] = ((mpt[2](1)-mpt[3](1))*(-0.5) -(mpt[0](1)-mpt[1](1))*(-0.5))/det ;

  b[0] = (-(mpt[2](0)-mpt[3](0))*(0.5) +(mpt[0](0)-mpt[1](0))*(0.5))/det ;
  b[1] = (-(mpt[2](0)-mpt[3](0))*(-0.5) +(mpt[0](0)-mpt[1](0))*(0.5))/det ;
  b[2] = (-(mpt[2](0)-mpt[3](0))*(0.5) +(mpt[0](0)-mpt[1](0))*(-0.5))/det ;
  b[3] = (-(mpt[2](0)-mpt[3](0))*(-0.5) +(mpt[0](0)-mpt[1](0))*(-0.5))/det ;

  c[0] = 0.25 - cpt(0)*a[0] - cpt(1)*b[0] ;
  c[1] = 0.25 - cpt(0)*a[1] - cpt(1)*b[1] ;
  c[2] = 0.25 - cpt(0)*a[2] - cpt(1)*b[2] ;
  c[3] = 0.25 - cpt(0)*a[3] - cpt(1)*b[3] ;



  // compute basis functions
  if (flags & (update_values | update_gradients))
    for (unsigned int i=0; i<n_q_points; ++i)
      {
        for (unsigned int k=0; k<this->dofs_per_cell; ++k)
          {
            if (flags & update_values)
              {
                values[k] = a[k]*mapping_data.quadrature_points[i](0) + b[k]*mapping_data.quadrature_points[i](1) + c[k] ;
                output_data.shape_values[k][i] = values[k];
              }

            if (flags & update_gradients)
              {
                grads[k][0] = a[k] ;
                grads[k][1] = b[k] ;
                output_data.shape_gradients[k][i] = grads[k];
              }
          }

      }

  // When this function is called for computation of facial jump residual,
  // MappingRelatedData does not work properly.
  // In this case, the quadrature points on the real cell is computed in manual sense, using 'mapping' and 'quadrature'.
  // This is a temporary solution. It needs to be fixed fundamentally.
  if (n_q_points==0)
    {
      Quadrature<2> cellquadrature = QProjector<2>::project_to_subface(quadrature, face_no, sub_no) ;
      for(unsigned int i=0;i<cellquadrature.size();++i)
               for (unsigned int k=0; k<this->dofs_per_cell; ++k)
               {
                   if (flags & update_values)
                   {
                   Point<2> realquadrature ;

                   realquadrature = mapping.transform_unit_to_real_cell(cell, cellquadrature.point(i)) ;
                   values[k] = a[k]*realquadrature(0) + b[k]*realquadrature(1) + c[k] ;
                   output_data.shape_values[k][i] = values[k];
                   }

                   if (flags & update_gradients)
                   {
                   grads[k][0] = a[k] ;
                   grads[k][1] = b[k] ;
                   output_data.shape_gradients[k][i] = grads[k];
                   }

               }

    }
}



FE_P1NC::FE_P1NC()
  :
  FiniteElement<2,2>(FiniteElementData<2>(get_dpo_vector(),
                                          1,
                                          1,
                                          FiniteElementData<2>::L2),
                     std::vector<bool> (1, false),
                     get_nonzero_component())
{

  // face support points: 2 end vertices
  unit_face_support_points.resize(2);
  unit_face_support_points[0][0] = 0.0 ;
  unit_face_support_points[1][0] = 1.0 ;


  // initialize constraints matrix
  initialize_constraints () ;
}



std::string FE_P1NC::get_name () const
{
  return "FE_P1NC";
}


UpdateFlags     FE_P1NC::requires_update_flags (const UpdateFlags flags) const
{
    UpdateFlags out = update_default;

    if (flags & update_values)
      out |= update_values | update_quadrature_points ;
    if (flags & update_gradients)
      out |= update_gradients ;
    if (flags & update_cell_normal_vectors)
      out |= update_cell_normal_vectors | update_JxW_values;
    if (flags & update_hessians)
      out |= update_hessians;
    if (flags & update_hessians)
      out |= update_hessians;
 
  return out;
}


FiniteElement<2,2> *FE_P1NC::clone () const
{
  return new FE_P1NC(*this);
}

FE_P1NC::~FE_P1NC () {}




void FE_P1NC::initialize_constraints ()
{

  std::vector<Point<1> > constraint_points;

  // Add midpoint
  constraint_points.push_back (Point<1> (0.5));

  // coefficient relation between children and mother
  interface_constraints
  .TableBase<2,double>::reinit (interface_constraints_size());


  interface_constraints(0,0) = 0.5 ;
  interface_constraints(0,1) = 0.5 ;

}





DEAL_II_NAMESPACE_CLOSE
