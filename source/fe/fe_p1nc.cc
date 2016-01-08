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


// FE_P1NCParametric

std::vector<ComponentMask>
FE_P1NCParametric::get_nonzero_component()
{
  const unsigned int dofs_per_cell = 4;
  std::vector<ComponentMask> masks (dofs_per_cell);
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    masks[i] = ComponentMask(1,true);
  return masks;
}


std::vector<unsigned int>
FE_P1NCParametric::get_dpo_vector ()
{
  std::vector<unsigned int> dpo(3);
  dpo[0] = 1; //dofs per object: vertex
  dpo[1] = 0; // line
  dpo[2] = 0; // quad
  return dpo;
}


FiniteElement<2,2>::InternalDataBase *
FE_P1NCParametric::get_data (const UpdateFlags update_flags,
                             const Mapping<2,2> &mapping,
                             const Quadrature<2> &quadrature) const
{

  InternalData *data = new InternalData;


  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  const UpdateFlags flags(data->update_flags);
  const unsigned int n_q_points = quadrature.size();


  // initialize shape value
  if (flags & update_values)
    data->shape_values.resize (this->dofs_per_cell,
                               std::vector<double> (n_q_points));

  // initialize shape gradient
  if (flags & update_gradients)
    data->shape_gradients.resize (this->dofs_per_cell,
                                  std::vector<Tensor<1,2> > (n_q_points));


  // update
  if (flags & (update_values | update_gradients))
    for (unsigned int i=0; i<n_q_points; ++i)
      {
        if (flags & update_values)
          for (unsigned int k=0; k<this->dofs_per_cell; ++k)
            {
              data->shape_values[k][i] = shape_value(k,quadrature.point(i));
            }


        if (flags & update_gradients)
          for (unsigned int k=0; k<this->dofs_per_cell; ++k)
            {
              data->shape_gradients[k][i] = shape_grad(k,quadrature.point(i));
            }

      }
  return data;

}


// get fe value in real cell from fedata in reference cell
void
FE_P1NCParametric::fill_fe_values (const Mapping<2,2>                                &mapping,
                                   const Triangulation<2,2>::cell_iterator           &cell,
                                   const Quadrature<2>                               &quadrature,
                                   const Mapping<2,2>::InternalDataBase              &mapping_internal,
                                   const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                                   const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                                   internal::FEValues::FiniteElementRelatedData<2,2> &output_data,
                                   const CellSimilarity::Similarity                   cell_similarity) const

{

  Assert (dynamic_cast<const InternalData *> (&fe_internal) != 0, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fe_internal);

  const UpdateFlags flags(fe_data.current_update_flags());

  // gets shape function values from its input argument fedata
  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      // shape value
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          output_data.shape_values(k,i) = fe_data.shape_values[k][i];


      // shape gradient
      if (flags & update_gradients && cell_similarity != CellSimilarity::translation)
        mapping.transform(fe_data.shape_gradients[k],
                          mapping_covariant,
                          mapping_internal,
                          output_data.shape_gradients[k]);

    }

}


void
FE_P1NCParametric::fill_fe_face_values (const Mapping<2,2>                                &mapping,
                                        const Triangulation<2,2>::cell_iterator           &cell,
                                        const unsigned int                                 face_no,
                                        const Quadrature<1>                               &quadrature,
                                        const Mapping<2,2>::InternalDataBase              &mapping_internal,
                                        const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                                        const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                                        internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const


{

  Assert (dynamic_cast<const InternalData *> (&fe_internal) != 0, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fe_internal);


  const QProjector<2>::DataSetDescriptor offset
    = QProjector<2>::DataSetDescriptor::face (face_no,
                                              cell->face_orientation(face_no),
                                              cell->face_flip(face_no),
                                              cell->face_rotation(face_no),
                                              quadrature.size());

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          output_data.shape_values(k,i) = fe_data.shape_values[k][i+offset];

      if (flags & update_gradients)
        mapping.transform(make_slice(fe_data.shape_gradients[k], offset, quadrature.size()),
                          mapping_covariant,
                          mapping_internal,
                          output_data.shape_gradients[k]);
    }


}






void
FE_P1NCParametric::fill_fe_subface_values (const Mapping<2,2>                                &mapping,
                                           const Triangulation<2,2>::cell_iterator           &cell,
                                           const unsigned int                                 face_no,
                                           const unsigned int                                 sub_no,
                                           const Quadrature<1>                               &quadrature,
                                           const Mapping<2,2>::InternalDataBase              &mapping_internal,
                                           const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                                           const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                                           internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const

{

  Assert (dynamic_cast<const InternalData *> (&fe_internal) != 0, ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &> (fe_internal);



  const QProjector<2>::DataSetDescriptor offset
    = QProjector<2>::DataSetDescriptor::subface (face_no, sub_no,
                                                 cell->face_orientation(face_no),
                                                 cell->face_flip(face_no),
                                                 cell->face_rotation(face_no),
                                                 quadrature.size(),
                                                 cell->subface_case(face_no));

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.size(); ++i)
          output_data.shape_values(k,i) = fe_data.shape_values[k][i+offset];

      if (flags & update_gradients)
        mapping.transform(make_slice(fe_data.shape_gradients[k], offset, quadrature.size()),
                          mapping_covariant,
                          mapping_internal,
                          output_data.shape_gradients[k]);
    }

}



FE_P1NCParametric::FE_P1NCParametric()
  :
  FiniteElement<2,2>(FiniteElementData<2>(get_dpo_vector(),
                                          1,
                                          1,
                                          FiniteElementData<2>::L2,
                                          numbers::invalid_unsigned_int),
                     std::vector<bool> (1, false),
                     get_nonzero_component())
{

  // support points: 4 vertices
  unit_support_points.resize(4) ;
  unit_support_points[0][0] = 0.0 ;
  unit_support_points[0][1] = 0.0 ;
  unit_support_points[1][0] = 1.0 ;
  unit_support_points[1][1] = 0.0 ;
  unit_support_points[2][0] = 0.0 ;
  unit_support_points[2][1] = 1.0 ;
  unit_support_points[3][0] = 1.0 ;
  unit_support_points[3][1] = 1.0 ;

  // face support points: 2 end vertices
  unit_face_support_points.resize(2);
  unit_face_support_points[0][0] = 0.0 ;
  unit_face_support_points[1][0] = 1.0 ;


  // CONSTRAINTS MATRIX INITIALIZATION
  initialize_constraints () ;
}



std::string FE_P1NCParametric::get_name () const
{
  return "FE_P1NCParametric";
}



UpdateFlags     FE_P1NCParametric::update_once (const UpdateFlags flags) const
{
  return (update_default | (flags & update_values));
}

UpdateFlags     FE_P1NCParametric::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;
  if (flags & update_cell_normal_vectors)
    out |= update_cell_normal_vectors | update_JxW_values;

  return out;
}


const std::vector<Point<1> > &FE_P1NCParametric::get_unit_face_support_points () const
{

  Assert ((unit_face_support_points.size() == 0) ||
          (unit_face_support_points.size() == this->dofs_per_face),
          ExcInternalError());
  return unit_face_support_points;
}

Point<1> FE_P1NCParametric::unit_face_support_point (const unsigned index) const
{
  Assert (index < this->dofs_per_face,
          ExcIndexRange (index, 0, this->dofs_per_face));
  Assert (unit_face_support_points.size() == this->dofs_per_face,
          ExcFEHasNoSupportPoints ());
  return unit_face_support_points[index];
}

bool FE_P1NCParametric::has_face_support_points () const
{
  return (unit_face_support_points.size() != 0);
}

FiniteElement<2,2> *FE_P1NCParametric::clone () const
{
  return new FE_P1NCParametric(*this);
}

FE_P1NCParametric::~FE_P1NCParametric () {}

// Define the exact value of the shape function for the P1NC scalar
// function in the unit cell = (0,1)^2.
//
//  2-----------------|------------------3
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  -                                    -
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  |                                    |
//  0-----------------|------------------1
//


double FE_P1NCParametric::shape_value (const unsigned int i,
                                       const Point<2>    &p) const
{


  double value = 0.0 ;

  // P1NC with HALF COEFFICIENT
  switch (i)
    {
    case 0:
      value = -p[0]-p[1]+1.5 ;
      break ;
    case 1:
      value =  p[0]-p[1]+0.5 ;
      break ;
    case 2:
      value = -p[0]+p[1]+0.5 ;
      break ;
    case 3:
      value =  p[0]+p[1]-0.5 ;
      break ;
    default:
      value = 0.0 ;
    };

  return value*0.5 ;


  // // Q1 for test
  // switch (i)
  //  {
  //  case 0:
  //    return  (1.0-p[0])*(1.0-p[1]) ;
  //  case 1:
  //    return  p[0]*(1.0-p[1]) ;
  //  case 2:
  //    return  (1.0-p[0])*p[1] ;
  //  case 3:
  //    return  p[0]*p[1] ;
  //  default:
  //    return 0 ;
  //  };


}




Tensor<1,2> FE_P1NCParametric::shape_grad (const unsigned int i,
                                           const Point<2>    &p) const
{
  Tensor<1,2> grad ;


  // P1NC with HALF COEFFICIENT
  switch (i)
    {
    case 0:
      grad[0] = -1.0 ;
      grad[1] = -1.0 ;
      break ;
    case 1:
      grad[0] =  1.0 ;
      grad[1] = -1.0 ;
      break ;
    case 2:
      grad[0] = -1.0 ;
      grad[1] =  1.0 ;
      break ;
    case 3:
      grad[0] =  1.0 ;
      grad[1] =  1.0 ;
      break ;
    default:
      grad = 0.0 ;
    };

  return grad*0.5 ;


  // // Q1 for test
  // switch (i)
  //  {
  //  case 0:
  //    grad[0] = p[1]-1.0 ;
  //    grad[1] = p[0]-1.0 ;
  //    return grad ;
  //  case 1:
  //    grad[0] = 1.0-p[1] ;
  //    grad[1] = -p[0] ;
  //    return grad ;
  //  case 2:
  //    grad[0] = -p[1] ;
  //    grad[1] = 1.0-p[0] ;
  //    return grad ;
  //  case 3:
  //    grad[0] =  p[1] ;
  //    grad[1] =  p[0] ;
  //    return grad ;
  //  default:
  //    return grad ;
  //  };


}


// 2015 05 13 CONSTRAINTS MATRIX FOR HANGING NODE
// BASED ON 'fe_q_base.cc'

void FE_P1NCParametric::initialize_constraints ()
{

  std::vector<Point<1> > constraint_points;
  // Add midpoint
  constraint_points.push_back (Point<1> (0.5));

  // Now construct relation between destination (child) and source (mother)
  // dofs.

  interface_constraints
  .TableBase<2,double>::reinit (interface_constraints_size());


  interface_constraints(0,0) = 0.5 ;
  interface_constraints(0,1) = 0.5 ;

}








// FE_P1NCNonparametric

std::vector<ComponentMask>
FE_P1NCNonparametric::get_nonzero_component()
{
  const unsigned int dofs_per_cell = 4;
  std::vector<ComponentMask> masks (dofs_per_cell);
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    masks[i] = ComponentMask(1,true);
  return masks;
}


std::vector<unsigned int>
FE_P1NCNonparametric::get_dpo_vector ()
{
  std::vector<unsigned int> dpo(3);
  dpo[0] = 1; // dofs per object: vertex
  dpo[1] = 0; // line
  dpo[2] = 0; // quad
  return dpo;
}


FiniteElement<2,2>::InternalDataBase *
FE_P1NCNonparametric::get_data (const UpdateFlags update_flags,
                                const Mapping<2,2> &,
                                const Quadrature<2> &) const
{
  FiniteElement<2,2>::InternalDataBase *data = new FiniteElement<2,2>::InternalDataBase;

  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  return data;




}



void
FE_P1NCNonparametric::fill_fe_values (const Mapping<2,2>                                &mapping,
                                      const Triangulation<2,2>::cell_iterator           &cell,
                                      const Quadrature<2>                               &quadrature,
                                      const Mapping<2,2>::InternalDataBase              &mapping_internal,
                                      const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                                      const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                                      internal::FEValues::FiniteElementRelatedData<2,2> &output_data,
                                      const CellSimilarity::Similarity                   cell_similarity) const

{
  const UpdateFlags flags(fe_internal.current_update_flags()) ;

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
FE_P1NCNonparametric::fill_fe_face_values (const Mapping<2,2>                                &mapping,
                                           const Triangulation<2,2>::cell_iterator           &cell,
                                           const unsigned int                                 face_no,
                                           const Quadrature<1>                               &quadrature,
                                           const Mapping<2,2>::InternalDataBase              &mapping_internal,
                                           const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                                           const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                                           internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const


{}






void
FE_P1NCNonparametric::fill_fe_subface_values (const Mapping<2,2>                                &mapping,
                                              const Triangulation<2,2>::cell_iterator           &cell,
                                              const unsigned int                                 face_no,
                                              const unsigned int                                 sub_no,
                                              const Quadrature<1>                               &quadrature,
                                              const Mapping<2,2>::InternalDataBase              &mapping_internal,
                                              const FiniteElement<2,2>::InternalDataBase        &fe_internal,
                                              const internal::FEValues::MappingRelatedData<2,2> &mapping_data,
                                              internal::FEValues::FiniteElementRelatedData<2,2> &output_data) const

{}



FE_P1NCNonparametric::FE_P1NCNonparametric()
  :
  FiniteElement<2,2>(FiniteElementData<2>(get_dpo_vector(),
                                          1,
                                          1,
                                          FiniteElementData<2>::L2,
                                          numbers::invalid_unsigned_int),
                     std::vector<bool> (1, false),
                     get_nonzero_component())
{

  // face support points: 2 end vertices
  unit_face_support_points.resize(2);
  unit_face_support_points[0][0] = 0.0 ;
  unit_face_support_points[1][0] = 1.0 ;


  // CONSTRAINTS MATRIX INITIALIZATION
  initialize_constraints () ;
}



std::string FE_P1NCNonparametric::get_name () const
{
  return "FE_P1NCNonparametric";
}



UpdateFlags     FE_P1NCNonparametric::update_once (const UpdateFlags flags) const
{
  return (update_default );

}

UpdateFlags     FE_P1NCNonparametric::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = (update_default | (flags & update_values) | (flags & update_quadrature_points)) ;

  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation | update_JxW_values;
  if (flags & update_cell_normal_vectors)
    out |= update_cell_normal_vectors | update_JxW_values;

  return out;
}


FiniteElement<2,2> *FE_P1NCNonparametric::clone () const
{
  return new FE_P1NCNonparametric(*this);
}

FE_P1NCNonparametric::~FE_P1NCNonparametric () {}




void FE_P1NCNonparametric::initialize_constraints ()
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
