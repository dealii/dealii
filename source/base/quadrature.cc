// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/utilities.h>

#include <cmath>
#include <cstdlib>
#include <iterator>

DEAL_II_NAMESPACE_OPEN


template <>
Quadrature<0>::Quadrature (const unsigned int n_q)
  :
  quadrature_points (n_q),
  weights (n_q, 0)
{}



template <>
Quadrature<0>::~Quadrature ()
{}



template <int dim>
Quadrature<dim>::Quadrature (const unsigned int n_q)
  :
  quadrature_points (n_q, Point<dim>()),
  weights (n_q, 0)
{}



template <int dim>
void
Quadrature<dim>::initialize (const std::vector<Point<dim> > &p,
                             const std::vector<double>      &w)
{
  AssertDimension (w.size(), p.size());
  quadrature_points = p;
  weights = w;
}



template <int dim>
Quadrature<dim>::Quadrature (const std::vector<Point<dim> > &points,
                             const std::vector<double>      &weights)
  :
  quadrature_points(points),
  weights(weights)
{
  Assert (weights.size() == points.size(),
          ExcDimensionMismatch(weights.size(), points.size()));
}



template <int dim>
Quadrature<dim>::Quadrature (const std::vector<Point<dim> > &points)
  :
  quadrature_points(points),
  weights(points.size(), std::atof("Inf"))
{
  Assert(weights.size() == points.size(),
         ExcDimensionMismatch(weights.size(), points.size()));
}



template <int dim>
Quadrature<dim>::Quadrature (const Point<dim> &point)
  :
  quadrature_points(std::vector<Point<dim> > (1, point)),
  weights(std::vector<double> (1, 1.))
{}


template <>
Quadrature<0>::Quadrature (const SubQuadrature &,
                           const Quadrature<1> &)
{
  Assert(false, ExcImpossibleInDim(0));
}



template <int dim>
Quadrature<dim>::Quadrature (const SubQuadrature &q1,
                             const Quadrature<1> &q2)
  :
  quadrature_points (q1.size() * q2.size()),
  weights (q1.size() * q2.size())
{
  unsigned int present_index = 0;
  for (unsigned int i2=0; i2<q2.size(); ++i2)
    for (unsigned int i1=0; i1<q1.size(); ++i1)
      {
        // compose coordinates of
        // new quadrature point by tensor
        // product in the last component
        for (unsigned int d=0; d<dim-1; ++d)
          quadrature_points[present_index](d)
            = q1.point(i1)(d);
        quadrature_points[present_index](dim-1)
          = q2.point(i2)(0);

        weights[present_index] = q1.weight(i1) * q2.weight(i2);

        ++present_index;
      };

#ifdef DEBUG
  if (size() > 0)
    {
      double sum = 0;
      for (unsigned int i=0; i<size(); ++i)
        sum += weights[i];
      // we cannot guarantee the sum of weights
      // to be exactly one, but it should be
      // near that.
      Assert ((sum>0.999999) && (sum<1.000001), ExcInternalError());
    }
#endif
}



template <>
Quadrature<1>::Quadrature (const SubQuadrature &,
                           const Quadrature<1> &q2)
  :
  quadrature_points (q2.size()),
  weights (q2.size())
{
  unsigned int present_index = 0;
  for (unsigned int i2=0; i2<q2.size(); ++i2)
    {
      // compose coordinates of
      // new quadrature point by tensor
      // product in the last component
      quadrature_points[present_index](0)
        = q2.point(i2)(0);

      weights[present_index] = q2.weight(i2);

      ++present_index;
    }

#ifdef DEBUG
  if (size() > 0)
    {
      double sum = 0;
      for (unsigned int i=0; i<size(); ++i)
        sum += weights[i];
      // we cannot guarantee the sum of weights
      // to be exactly one, but it should be
      // near that.
      Assert ((sum>0.999999) && (sum<1.000001), ExcInternalError());
    }
#endif
}



template <>
Quadrature<0>::Quadrature (const Quadrature<1> &)
  :
  Subscriptor(),
//              quadrature_points(1),
  weights(1,1.)
{}


template <>
Quadrature<1>::Quadrature (const Quadrature<0> &)
  :
  Subscriptor()
{
  // this function should never be
  // called -- this should be the
  // copy constructor in 1d...
  Assert (false, ExcImpossibleInDim(1));
}



template <int dim>
Quadrature<dim>::Quadrature (const Quadrature<dim != 1 ? 1 : 0> &q)
:
Subscriptor(),
            quadrature_points (Utilities::fixed_power<dim>(q.size())),
            weights (Utilities::fixed_power<dim>(q.size()))
{
  Assert (dim <= 3, ExcNotImplemented());

  const unsigned int n0 = q.size();
  const unsigned int n1 = (dim>1) ? n0 : 1;
  const unsigned int n2 = (dim>2) ? n0 : 1;

  unsigned int k=0;
  for (unsigned int i2=0; i2<n2; ++i2)
    for (unsigned int i1=0; i1<n1; ++i1)
      for (unsigned int i0=0; i0<n0; ++i0)
        {
          quadrature_points[k](0) = q.point(i0)(0);
          if (dim>1)
            quadrature_points[k](1) = q.point(i1)(0);
          if (dim>2)
            quadrature_points[k](2) = q.point(i2)(0);
          weights[k] = q.weight(i0);
          if (dim>1)
            weights[k] *= q.weight(i1);
          if (dim>2)
            weights[k] *= q.weight(i2);
          ++k;
        }
}



template <int dim>
Quadrature<dim>::Quadrature (const Quadrature<dim> &q)
  :
  Subscriptor(),
  quadrature_points (q.quadrature_points),
  weights (q.weights)
{}


template <int dim>
Quadrature<dim> &
Quadrature<dim>::operator= (const Quadrature<dim> &q)
{
  weights = q.weights;
  quadrature_points = q.quadrature_points;
  return *this;
}



template <int dim>
bool
Quadrature<dim>::operator == (const Quadrature<dim> &q) const
{
  return ((quadrature_points == q.quadrature_points)
          &&
          (weights == q.weights));
}



template <int dim>
Quadrature<dim>::~Quadrature ()
{}



template <int dim>
std::size_t
Quadrature<dim>::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (quadrature_points) +
          MemoryConsumption::memory_consumption (weights));
}


//---------------------------------------------------------------------------
template<int dim>
QAnisotropic<dim>::QAnisotropic(const Quadrature<1> &qx)
  : Quadrature<dim>(qx.size())
{
  Assert (dim==1, ExcImpossibleInDim(dim));
  unsigned int k=0;
  for (unsigned int k1=0; k1<qx.size(); ++k1)
    {
      this->quadrature_points[k](0) = qx.point(k1)(0);
      this->weights[k++] = qx.weight(k1);
    }
  Assert (k==this->size(), ExcInternalError());
}



template<int dim>
QAnisotropic<dim>::QAnisotropic(const Quadrature<1> &qx,
                                const Quadrature<1> &qy)
  : Quadrature<dim>(qx.size()
                    *qy.size())
{
  Assert (dim==2, ExcImpossibleInDim(dim));
  unsigned int k=0;
  for (unsigned int k2=0; k2<qy.size(); ++k2)
    for (unsigned int k1=0; k1<qx.size(); ++k1)
      {
        this->quadrature_points[k](0) = qx.point(k1)(0);
        this->quadrature_points[k](1) = qy.point(k2)(0);
        this->weights[k++] = qx.weight(k1) * qy.weight(k2);
      }
  Assert (k==this->size(), ExcInternalError());
}



template<int dim>
QAnisotropic<dim>::QAnisotropic(const Quadrature<1> &qx,
                                const Quadrature<1> &qy,
                                const Quadrature<1> &qz)
  : Quadrature<dim>(qx.size()
                    *qy.size()
                    *qz.size())
{
  Assert (dim==3, ExcImpossibleInDim(dim));
  unsigned int k=0;
  for (unsigned int k3=0; k3<qz.size(); ++k3)
    for (unsigned int k2=0; k2<qy.size(); ++k2)
      for (unsigned int k1=0; k1<qx.size(); ++k1)
        {
          this->quadrature_points[k](0) = qx.point(k1)(0);
          this->quadrature_points[k](1) = qy.point(k2)(0);
          this->quadrature_points[k](2) = qz.point(k3)(0);
          this->weights[k++] = qx.weight(k1) * qy.weight(k2) * qz.weight(k3);
        }
  Assert (k==this->size(), ExcInternalError());
}



//---------------------------------------------------------------------------



template <int dim>
Quadrature<2>
QProjector<dim>::reflect (const Quadrature<2> &q)
{
  std::vector<Point<2> > q_points (q.size());
  std::vector<double>    weights (q.size());
  for (unsigned int i=0; i<q.size(); ++i)
    {
      q_points[i][0] = q.point(i)[1];
      q_points[i][1] = q.point(i)[0];

      weights[i] = q.weight(i);
    }

  return Quadrature<2> (q_points, weights);
}


template <int dim>
Quadrature<2>
QProjector<dim>::rotate (const Quadrature<2> &q,
                         const unsigned int   n_times)
{
  std::vector<Point<2> > q_points (q.size());
  std::vector<double>    weights (q.size());
  for (unsigned int i=0; i<q.size(); ++i)
    {
      switch (n_times%4)
        {
        case 0:
          // 0 degree
          q_points[i][0] = q.point(i)[0];
          q_points[i][1] = q.point(i)[1];
          break;
        case 1:
          // 90 degree counterclockwise
          q_points[i][0] = 1.0 - q.point(i)[1];
          q_points[i][1] = q.point(i)[0];
          break;
        case 2:
          // 180 degree counterclockwise
          q_points[i][0] = 1.0 - q.point(i)[0];
          q_points[i][1] = 1.0 - q.point(i)[1];
          break;
        case 3:
          // 270 degree counterclockwise
          q_points[i][0] = q.point(i)[1];
          q_points[i][1] = 1.0 - q.point(i)[0];
          break;
        }

      weights[i] = q.weight(i);
    }

  return Quadrature<2> (q_points, weights);
}


template <>
void
QProjector<1>::project_to_face (const Quadrature<0> &,
                                const unsigned int face_no,
                                std::vector<Point<1> > &q_points)
{
  const unsigned int dim=1;
  AssertIndexRange (face_no, GeometryInfo<dim>::faces_per_cell);
  AssertDimension (q_points.size(), 1);

  q_points[0] = Point<dim>((double) face_no);
}



template <>
void
QProjector<2>::project_to_face (const Quadrature<1>      &quadrature,
                                const unsigned int        face_no,
                                std::vector<Point<2> >   &q_points)
{
  const unsigned int dim=2;
  AssertIndexRange (face_no, GeometryInfo<dim>::faces_per_cell);
  Assert (q_points.size() == quadrature.size(),
          ExcDimensionMismatch (q_points.size(), quadrature.size()));

  for (unsigned int p=0; p<quadrature.size(); ++p)
    switch (face_no)
      {
      case 0:
        q_points[p] = Point<dim>(0,quadrature.point(p)(0));
        break;
      case 1:
        q_points[p] = Point<dim>(1,quadrature.point(p)(0));
        break;
      case 2:
        q_points[p] = Point<dim>(quadrature.point(p)(0),0);
        break;
      case 3:
        q_points[p] = Point<dim>(quadrature.point(p)(0),1);
        break;
      default:
        Assert (false, ExcInternalError());
      };
}



template <>
void
QProjector<3>::project_to_face (const Quadrature<2>    &quadrature,
                                const unsigned int      face_no,
                                std::vector<Point<3> > &q_points)
{
  const unsigned int dim=3;
  AssertIndexRange (face_no, GeometryInfo<dim>::faces_per_cell);
  Assert (q_points.size() == quadrature.size(),
          ExcDimensionMismatch (q_points.size(), quadrature.size()));

  for (unsigned int p=0; p<quadrature.size(); ++p)
    switch (face_no)
      {
      case 0:
        q_points[p] = Point<dim>(0,
                                 quadrature.point(p)(0),
                                 quadrature.point(p)(1));
        break;
      case 1:
        q_points[p] = Point<dim>(1,
                                 quadrature.point(p)(0),
                                 quadrature.point(p)(1));
        break;
      case 2:
        q_points[p] = Point<dim>(quadrature.point(p)(1),
                                 0,
                                 quadrature.point(p)(0));
        break;
      case 3:
        q_points[p] = Point<dim>(quadrature.point(p)(1),
                                 1,
                                 quadrature.point(p)(0));
        break;
      case 4:
        q_points[p] = Point<dim>(quadrature.point(p)(0),
                                 quadrature.point(p)(1),
                                 0);
        break;
      case 5:
        q_points[p] = Point<dim>(quadrature.point(p)(0),
                                 quadrature.point(p)(1),
                                 1);
        break;

      default:
        Assert (false, ExcInternalError());
      };
}



template <>
void
QProjector<1>::project_to_subface (const Quadrature<0> &,
                                   const unsigned int face_no,
                                   const unsigned int,
                                   std::vector<Point<1> > &q_points,
                                   const RefinementCase<0> &)
{
  const unsigned int dim=1;
  AssertIndexRange (face_no, GeometryInfo<dim>::faces_per_cell);
  AssertDimension (q_points.size(), 1);

  q_points[0] = Point<dim>((double) face_no);
}



template <>
void
QProjector<2>::project_to_subface (const Quadrature<1>    &quadrature,
                                   const unsigned int      face_no,
                                   const unsigned int      subface_no,
                                   std::vector<Point<2> > &q_points,
                                   const RefinementCase<1> &)
{
  const unsigned int dim=2;
  AssertIndexRange (face_no, GeometryInfo<dim>::faces_per_cell);
  AssertIndexRange (subface_no, GeometryInfo<dim>::max_children_per_face);

  Assert (q_points.size() == quadrature.size(),
          ExcDimensionMismatch (q_points.size(), quadrature.size()));

  for (unsigned int p=0; p<quadrature.size(); ++p)
    switch (face_no)
      {
      case 0:
        switch (subface_no)
          {
          case 0:
            q_points[p] = Point<dim>(0,quadrature.point(p)(0)/2);
            break;
          case 1:
            q_points[p] = Point<dim>(0,quadrature.point(p)(0)/2+0.5);
            break;
          default:
            Assert (false, ExcInternalError());
          };
        break;
      case 1:
        switch (subface_no)
          {
          case 0:
            q_points[p] = Point<dim>(1,quadrature.point(p)(0)/2);
            break;
          case 1:
            q_points[p] = Point<dim>(1,quadrature.point(p)(0)/2+0.5);
            break;
          default:
            Assert (false, ExcInternalError());
          };
        break;
      case 2:
        switch (subface_no)
          {
          case 0:
            q_points[p]
              = Point<dim>(quadrature.point(p)(0)/2,0);
            break;
          case 1:
            q_points[p]
              = Point<dim>(quadrature.point(p)(0)/2+0.5,0);
            break;
          default:
            Assert (false, ExcInternalError());
          };
        break;
      case 3:
        switch (subface_no)
          {
          case 0:
            q_points[p] = Point<dim>(quadrature.point(p)(0)/2,1);
            break;
          case 1:
            q_points[p] = Point<dim>(quadrature.point(p)(0)/2+0.5,1);
            break;
          default:
            Assert (false, ExcInternalError());
          };
        break;

      default:
        Assert (false, ExcInternalError());
      };
}



template <>
void
QProjector<3>::project_to_subface (const Quadrature<2>    &quadrature,
                                   const unsigned int      face_no,
                                   const unsigned int      subface_no,
                                   std::vector<Point<3> > &q_points,
                                   const RefinementCase<2> &ref_case)
{
  const unsigned int dim=3;
  AssertIndexRange (face_no, GeometryInfo<dim>::faces_per_cell);
  AssertIndexRange (subface_no, GeometryInfo<dim>::max_children_per_face);
  Assert (q_points.size() == quadrature.size(),
          ExcDimensionMismatch (q_points.size(), quadrature.size()));

  // one coordinate is at a const value. for
  // faces 0, 2 and 4 this value is 0.0, for
  // faces 1, 3 and 5 it is 1.0
  double const_value=face_no%2;
  // local 2d coordinates are xi and eta,
  // global 3d coordinates are x, y and
  // z. those have to be mapped. the following
  // indices tell, which global coordinate
  // (0->x, 1->y, 2->z) corresponds to which
  // local one
  unsigned int xi_index   = deal_II_numbers::invalid_unsigned_int,
               eta_index   = deal_II_numbers::invalid_unsigned_int,
               const_index = face_no/2;
  // the xi and eta values have to be scaled
  // (by factor 0.5 or factor 1.0) depending on
  // the refinement case and translated (by 0.0
  // or 0.5) depending on the refinement case
  // and subface_no.
  double xi_scale=1.0,
         eta_scale=1.0,
         xi_translation=0.0,
         eta_translation=0.0;
  // set the index mapping between local and
  // global coordinates
  switch (face_no/2)
    {
    case 0:
      xi_index=1;
      eta_index=2;
      break;
    case 1:
      xi_index=2;
      eta_index=0;
      break;
    case 2:
      xi_index=0;
      eta_index=1;
      break;
    }
  // set the scale and translation parameter
  // for individual subfaces
  switch ((unsigned char)ref_case)
    {
    case RefinementCase<dim-1>::cut_x:
      xi_scale=0.5;
      xi_translation=subface_no%2 * 0.5;
      break;
    case RefinementCase<dim-1>::cut_y:
      eta_scale=0.5;
      eta_translation=subface_no%2 * 0.5;
      break;
    case RefinementCase<dim-1>::cut_xy:
      xi_scale= 0.5;
      eta_scale=0.5;
      xi_translation =subface_no%2 * 0.5;
      eta_translation=subface_no/2 * 0.5;
      break;
    default:
      Assert(false,ExcInternalError());
      break;
    }
  // finally, compute the scaled, translated,
  // projected quadrature points
  for (unsigned int p=0; p<quadrature.size(); ++p)
    {
      q_points[p][xi_index]    = xi_scale  * quadrature.point(p)(0) + xi_translation;
      q_points[p][eta_index]   = eta_scale * quadrature.point(p)(1) + eta_translation;
      q_points[p][const_index] = const_value;
    }
}


template <>
Quadrature<1>
QProjector<1>::project_to_all_faces (const Quadrature<0> &quadrature)
{
  const unsigned int dim = 1;

  const unsigned int n_points = 1,
                     n_faces  = GeometryInfo<dim>::faces_per_cell;

  // first fix quadrature points
  std::vector<Point<dim> > q_points;
  q_points.reserve(n_points * n_faces);
  std::vector <Point<dim> > help(n_points);


  // project to each face and append
  // results
  for (unsigned int face=0; face<n_faces; ++face)
    {
      project_to_face(quadrature, face, help);
      std::copy (help.begin(), help.end(),
                 std::back_inserter (q_points));
    }

  // next copy over weights
  std::vector<double> weights;
  weights.reserve (n_points * n_faces);
  for (unsigned int face=0; face<n_faces; ++face)
    std::copy (quadrature.get_weights().begin(),
               quadrature.get_weights().end(),
               std::back_inserter (weights));

  Assert (q_points.size() == n_points * n_faces,
          ExcInternalError());
  Assert (weights.size() == n_points * n_faces,
          ExcInternalError());

  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<2>
QProjector<2>::project_to_all_faces (const SubQuadrature &quadrature)
{
  const unsigned int dim = 2;

  const unsigned int n_points = quadrature.size(),
                     n_faces  = GeometryInfo<dim>::faces_per_cell;

  // first fix quadrature points
  std::vector<Point<dim> > q_points;
  q_points.reserve(n_points * n_faces);
  std::vector <Point<dim> > help(n_points);

  // project to each face and append
  // results
  for (unsigned int face=0; face<n_faces; ++face)
    {
      project_to_face(quadrature, face, help);
      std::copy (help.begin(), help.end(),
                 std::back_inserter (q_points));
    }

  // next copy over weights
  std::vector<double> weights;
  weights.reserve (n_points * n_faces);
  for (unsigned int face=0; face<n_faces; ++face)
    std::copy (quadrature.get_weights().begin(),
               quadrature.get_weights().end(),
               std::back_inserter (weights));

  Assert (q_points.size() == n_points * n_faces,
          ExcInternalError());
  Assert (weights.size() == n_points * n_faces,
          ExcInternalError());

  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<3>
QProjector<3>::project_to_all_faces (const SubQuadrature &quadrature)
{
  const unsigned int dim = 3;

  SubQuadrature q_reflected=reflect (quadrature);
  SubQuadrature q[8]=
  {
    quadrature,
    rotate (quadrature,1),
    rotate (quadrature,2),
    rotate (quadrature,3),
    q_reflected,
    rotate (q_reflected,3),
    rotate (q_reflected,2),
    rotate (q_reflected,1)
  };



  const unsigned int n_points = quadrature.size(),
                     n_faces  = GeometryInfo<dim>::faces_per_cell;

  // first fix quadrature points
  std::vector<Point<dim> > q_points;
  q_points.reserve(n_points * n_faces * 8);
  std::vector <Point<dim> > help(n_points);

  std::vector<double> weights;
  weights.reserve (n_points * n_faces * 8);

  // do the following for all possible
  // mutations of a face (mutation==0
  // corresponds to a face with standard
  // orientation, no flip and no rotation)
  for (unsigned int mutation=0; mutation<8; ++mutation)
    {
      // project to each face and append
      // results
      for (unsigned int face=0; face<n_faces; ++face)
        {
          project_to_face(q[mutation], face, help);
          std::copy (help.begin(), help.end(),
                     std::back_inserter (q_points));
        }

      // next copy over weights
      for (unsigned int face=0; face<n_faces; ++face)
        std::copy (q[mutation].get_weights().begin(),
                   q[mutation].get_weights().end(),
                   std::back_inserter (weights));
    }


  Assert (q_points.size() == n_points * n_faces * 8,
          ExcInternalError());
  Assert (weights.size() == n_points * n_faces * 8,
          ExcInternalError());

  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<1>
QProjector<1>::project_to_all_subfaces (const Quadrature<0> &quadrature)
{
  const unsigned int dim = 1;

  const unsigned int n_points          = 1,
                     n_faces           = GeometryInfo<dim>::faces_per_cell,
                     subfaces_per_face = GeometryInfo<dim>::max_children_per_face;

  // first fix quadrature points
  std::vector<Point<dim> > q_points;
  q_points.reserve (n_points * n_faces * subfaces_per_face);
  std::vector <Point<dim> > help(n_points);

  // project to each face and copy
  // results
  for (unsigned int face=0; face<n_faces; ++face)
    for (unsigned int subface=0; subface<subfaces_per_face; ++subface)
      {
        project_to_subface(quadrature, face, subface, help);
        std::copy (help.begin(), help.end(),
                   std::back_inserter (q_points));
      };

  // next copy over weights
  std::vector<double> weights;
  weights.reserve (n_points * n_faces * subfaces_per_face);
  for (unsigned int face=0; face<n_faces; ++face)
    for (unsigned int subface=0; subface<subfaces_per_face; ++subface)
      std::copy (quadrature.get_weights().begin(),
                 quadrature.get_weights().end(),
                 std::back_inserter (weights));

  Assert (q_points.size() == n_points * n_faces * subfaces_per_face,
          ExcInternalError());
  Assert (weights.size() == n_points * n_faces * subfaces_per_face,
          ExcInternalError());

  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<2>
QProjector<2>::project_to_all_subfaces (const SubQuadrature &quadrature)
{
  const unsigned int dim = 2;

  const unsigned int n_points          = quadrature.size(),
                     n_faces           = GeometryInfo<dim>::faces_per_cell,
                     subfaces_per_face = GeometryInfo<dim>::max_children_per_face;

  // first fix quadrature points
  std::vector<Point<dim> > q_points;
  q_points.reserve (n_points * n_faces * subfaces_per_face);
  std::vector <Point<dim> > help(n_points);

  // project to each face and copy
  // results
  for (unsigned int face=0; face<n_faces; ++face)
    for (unsigned int subface=0; subface<subfaces_per_face; ++subface)
      {
        project_to_subface(quadrature, face, subface, help);
        std::copy (help.begin(), help.end(),
                   std::back_inserter (q_points));
      };

  // next copy over weights
  std::vector<double> weights;
  weights.reserve (n_points * n_faces * subfaces_per_face);
  for (unsigned int face=0; face<n_faces; ++face)
    for (unsigned int subface=0; subface<subfaces_per_face; ++subface)
      std::copy (quadrature.get_weights().begin(),
                 quadrature.get_weights().end(),
                 std::back_inserter (weights));

  Assert (q_points.size() == n_points * n_faces * subfaces_per_face,
          ExcInternalError());
  Assert (weights.size() == n_points * n_faces * subfaces_per_face,
          ExcInternalError());

  return Quadrature<dim>(q_points, weights);
}



template <>
Quadrature<3>
QProjector<3>::project_to_all_subfaces (const SubQuadrature &quadrature)
{
  const unsigned int dim = 3;
  SubQuadrature q_reflected=reflect (quadrature);
  SubQuadrature q[8]=
  {
    quadrature,
    rotate (quadrature,1),
    rotate (quadrature,2),
    rotate (quadrature,3),
    q_reflected,
    rotate (q_reflected,3),
    rotate (q_reflected,2),
    rotate (q_reflected,1)
  };

  const unsigned int n_points          = quadrature.size(),
                     n_faces           = GeometryInfo<dim>::faces_per_cell,
                     total_subfaces_per_face = 2 + 2 + 4;

  // first fix quadrature points
  std::vector<Point<dim> > q_points;
  q_points.reserve (n_points * n_faces * total_subfaces_per_face * 8);
  std::vector <Point<dim> > help(n_points);

  std::vector<double> weights;
  weights.reserve (n_points * n_faces * total_subfaces_per_face * 8);

  // do the following for all possible
  // mutations of a face (mutation==0
  // corresponds to a face with standard
  // orientation, no flip and no rotation)
  for (unsigned int mutation=0; mutation<8; ++mutation)
    {
      // project to each face and copy
      // results
      for (unsigned int face=0; face<n_faces; ++face)
        for (unsigned int ref_case=RefinementCase<dim-1>::cut_xy;
             ref_case>=RefinementCase<dim-1>::cut_x;
             --ref_case)
          for (unsigned int subface=0; subface<GeometryInfo<dim-1>::n_children(RefinementCase<dim-1>(ref_case)); ++subface)
            {
              project_to_subface(q[mutation], face, subface, help,
                                 RefinementCase<dim-1>(ref_case));
              std::copy (help.begin(), help.end(),
                         std::back_inserter (q_points));
            }

      // next copy over weights
      for (unsigned int face=0; face<n_faces; ++face)
        for (unsigned int ref_case=RefinementCase<dim-1>::cut_xy;
             ref_case>=RefinementCase<dim-1>::cut_x;
             --ref_case)
          for (unsigned int subface=0; subface<GeometryInfo<dim-1>::n_children(RefinementCase<dim-1>(ref_case)); ++subface)
            std::copy (q[mutation].get_weights().begin(),
                       q[mutation].get_weights().end(),
                       std::back_inserter (weights));
    }

  Assert (q_points.size() == n_points * n_faces * total_subfaces_per_face * 8,
          ExcInternalError());
  Assert (weights.size() == n_points * n_faces * total_subfaces_per_face * 8,
          ExcInternalError());

  return Quadrature<dim>(q_points, weights);
}



// This function is not used in the library
template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_child (const Quadrature<dim>    &quadrature,
                                   const unsigned int        child_no)
{
  Assert (child_no < GeometryInfo<dim>::max_children_per_cell,
          ExcIndexRange (child_no, 0, GeometryInfo<dim>::max_children_per_cell));

  const unsigned int n_q_points = quadrature.size();

  std::vector<Point<dim> > q_points(n_q_points);
  for (unsigned int i=0; i<n_q_points; ++i)
    q_points[i]=GeometryInfo<dim>::child_to_cell_coordinates(
                  quadrature.point(i), child_no);

  // for the weights, things are
  // equally simple: copy them and
  // scale them
  std::vector<double> weights = quadrature.get_weights ();
  for (unsigned int i=0; i<n_q_points; ++i)
    weights[i] *= (1./GeometryInfo<dim>::max_children_per_cell);

  return Quadrature<dim> (q_points, weights);
}


template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_all_children (const Quadrature<dim> &quadrature)
{
  const unsigned int n_points = quadrature.size(),
                     n_children  = GeometryInfo<dim>::max_children_per_cell;

  std::vector<Point<dim> > q_points(n_points * n_children);
  std::vector<double> weights(n_points * n_children);

  // project to each child and copy
  // results
  for (unsigned int child=0; child<n_children; ++child)
    {
      Quadrature<dim> help = project_to_child(quadrature, child);
      for (unsigned int i=0; i<n_points; ++i)
        {
          q_points[child*n_points+i] = help.point(i);
          weights[child*n_points+i] = help.weight(i);
        }
    }
  return Quadrature<dim>(q_points, weights);
}




template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_line(
  const Quadrature<1> &quadrature,
  const Point<dim> &p1,
  const Point<dim> &p2)
{
  const unsigned int n = quadrature.size();
  std::vector<Point<dim> > points(n);
  std::vector<double> weights(n);
  const double length = p1.distance(p2);

  for (unsigned int k=0; k<n; ++k)
    {
      const double alpha = quadrature.point(k)(0);
      points[k] = alpha * p2;
      points[k] += (1.-alpha) * p1;
      weights[k] = length * quadrature.weight(k);
    }
  return Quadrature<dim> (points, weights);
}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::
face (const unsigned int face_no,
      const bool         face_orientation,
      const bool         face_flip,
      const bool         face_rotation,
      const unsigned int n_quadrature_points)
{
  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
          ExcInternalError());

  switch (dim)
    {
    case 1:
    case 2:
      return face_no * n_quadrature_points;


    case 3:
    {
      // in 3d, we have to account for faces that
      // have non-standard face orientation, flip
      // and rotation. thus, we have to store
      // _eight_ data sets per face or subface

      // set up a table with the according offsets
      // for non-standard orientation, first index:
      // face_orientation (standard true=1), second
      // index: face_flip (standard false=0), third
      // index: face_rotation (standard false=0)
      //
      // note, that normally we should use the
      // obvious offsets 0,1,2,3,4,5,6,7. However,
      // prior to the changes enabling flipped and
      // rotated faces, in many places of the
      // library the convention was used, that the
      // first dataset with offset 0 corresponds to
      // a face in standard orientation. therefore
      // we use the offsets 4,5,6,7,0,1,2,3 here to
      // stick to that (implicit) convention
      static const unsigned int offset[2][2][2]=
      {
        { {4*GeometryInfo<dim>::faces_per_cell, 5*GeometryInfo<dim>::faces_per_cell},    // face_orientation=false; face_flip=false; face_rotation=false and true
          {6*GeometryInfo<dim>::faces_per_cell, 7*GeometryInfo<dim>::faces_per_cell}
        },   // face_orientation=false; face_flip=true;  face_rotation=false and true
        { {0*GeometryInfo<dim>::faces_per_cell, 1*GeometryInfo<dim>::faces_per_cell},    // face_orientation=true;  face_flip=false; face_rotation=false and true
          {2*GeometryInfo<dim>::faces_per_cell, 3*GeometryInfo<dim>::faces_per_cell}
        }
      };  // face_orientation=true;  face_flip=true;  face_rotation=false and true

      return ((face_no
               + offset[face_orientation][face_flip][face_rotation])
              * n_quadrature_points);
    }

    default:
      Assert (false, ExcInternalError());
    }
  return numbers::invalid_unsigned_int;
}



template <>
QProjector<1>::DataSetDescriptor
QProjector<1>::DataSetDescriptor::
subface (const unsigned int face_no,
         const unsigned int subface_no,
         const bool,
         const bool,
         const bool,
         const unsigned int n_quadrature_points,
         const internal::SubfaceCase<1>)
{
  Assert (face_no < GeometryInfo<1>::faces_per_cell,
          ExcInternalError());
  Assert (subface_no < GeometryInfo<1>::max_children_per_face,
          ExcInternalError());

  return ((face_no * GeometryInfo<1>::max_children_per_face +
           subface_no)
          * n_quadrature_points);
}



template <>
QProjector<2>::DataSetDescriptor
QProjector<2>::DataSetDescriptor::
subface (const unsigned int face_no,
         const unsigned int subface_no,
         const bool,
         const bool,
         const bool,
         const unsigned int n_quadrature_points,
         const internal::SubfaceCase<2>)
{
  Assert (face_no < GeometryInfo<2>::faces_per_cell,
          ExcInternalError());
  Assert (subface_no < GeometryInfo<2>::max_children_per_face,
          ExcInternalError());

  return ((face_no * GeometryInfo<2>::max_children_per_face +
           subface_no)
          * n_quadrature_points);
}


template <>
QProjector<3>::DataSetDescriptor
QProjector<3>::DataSetDescriptor::
subface (const unsigned int face_no,
         const unsigned int subface_no,
         const bool         face_orientation,
         const bool         face_flip,
         const bool         face_rotation,
         const unsigned int n_quadrature_points,
         const internal::SubfaceCase<3> ref_case)
{
  const unsigned int dim = 3;

  Assert (face_no < GeometryInfo<dim>::faces_per_cell,
          ExcInternalError());
  Assert (subface_no < GeometryInfo<dim>::max_children_per_face,
          ExcInternalError());

  // As the quadrature points created by
  // QProjector are on subfaces in their
  // "standard location" we have to use a
  // permutation of the equivalent subface
  // number in order to respect face
  // orientation, flip and rotation. The
  // information we need here is exactly the
  // same as the
  // GeometryInfo<3>::child_cell_on_face info
  // for the bottom face (face 4) of a hex, as
  // on this the RefineCase of the cell matches
  // that of the face and the subfaces are
  // numbered in the same way as the child
  // cells.

  // in 3d, we have to account for faces that
  // have non-standard face orientation, flip
  // and rotation. thus, we have to store
  // _eight_ data sets per face or subface
  // already for the isotropic
  // case. Additionally, we have three
  // different refinement cases, resulting in
  // <tt>4 + 2 + 2 = 8</tt> different subfaces
  // for each face.
  const unsigned int total_subfaces_per_face=8;

  // set up a table with the according offsets
  // for non-standard orientation, first index:
  // face_orientation (standard true=1), second
  // index: face_flip (standard false=0), third
  // index: face_rotation (standard false=0)
  //
  // note, that normally we should use the
  // obvious offsets 0,1,2,3,4,5,6,7. However,
  // prior to the changes enabling flipped and
  // rotated faces, in many places of the
  // library the convention was used, that the
  // first dataset with offset 0 corresponds to
  // a face in standard orientation. therefore
  // we use the offsets 4,5,6,7,0,1,2,3 here to
  // stick to that (implicit) convention
  static const unsigned int orientation_offset[2][2][2]=
  {
    {
      // face_orientation=false; face_flip=false; face_rotation=false and true
      {
        4*GeometryInfo<dim>::faces_per_cell*total_subfaces_per_face,
        5*GeometryInfo<dim>::faces_per_cell *total_subfaces_per_face
      },
      // face_orientation=false; face_flip=true;  face_rotation=false and true
      {
        6*GeometryInfo<dim>::faces_per_cell*total_subfaces_per_face,
        7*GeometryInfo<dim>::faces_per_cell *total_subfaces_per_face
      }
    },
    {
      // face_orientation=true;  face_flip=false; face_rotation=false and true
      {
        0*GeometryInfo<dim>::faces_per_cell*total_subfaces_per_face,
        1*GeometryInfo<dim>::faces_per_cell *total_subfaces_per_face
      },
      // face_orientation=true;  face_flip=true;  face_rotation=false and true
      {
        2*GeometryInfo<dim>::faces_per_cell*total_subfaces_per_face,
        3*GeometryInfo<dim>::faces_per_cell *total_subfaces_per_face
      }
    }
  };

  // set up a table with the offsets for a
  // given refinement case respecting the
  // corresponding number of subfaces. the
  // index corresponds to (RefineCase::Type - 1)

  // note, that normally we should use the
  // obvious offsets 0,2,6. However, prior to
  // the implementation of anisotropic
  // refinement, in many places of the library
  // the convention was used, that the first
  // dataset with offset 0 corresponds to a
  // standard (isotropic) face
  // refinement. therefore we use the offsets
  // 6,4,0 here to stick to that (implicit)
  // convention
  static const unsigned int ref_case_offset[3]=
  {
    6,  //cut_x
    4,  //cut_y
    0   //cut_xy
  };


  // for each subface of a given FaceRefineCase
  // there is a corresponding equivalent
  // subface number of one of the "standard"
  // RefineCases (cut_x, cut_y, cut_xy). Map
  // the given values to those equivalent
  // ones.

  // first, define an invalid number
  static const unsigned int e = deal_II_numbers::invalid_unsigned_int;

  static const RefinementCase<dim-1>
  equivalent_refine_case[internal::SubfaceCase<dim>::case_isotropic+1][GeometryInfo<3>::max_children_per_face]
  =
  {
    // case_none. there should be only
    // invalid values here. However, as
    // this function is also called (in
    // tests) for cells which have no
    // refined faces, use isotropic
    // refinement instead
    {
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy
    },
    // case_x
    {
      RefinementCase<dim-1>::cut_x,
      RefinementCase<dim-1>::cut_x,
      RefinementCase<dim-1>::no_refinement,
      RefinementCase<dim-1>::no_refinement
    },
    // case_x1y
    {
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_x,
      RefinementCase<dim-1>::no_refinement
    },
    // case_x2y
    {
      RefinementCase<dim-1>::cut_x,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::no_refinement
    },
    // case_x1y2y
    {
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy
    },
    // case_y
    {
      RefinementCase<dim-1>::cut_y,
      RefinementCase<dim-1>::cut_y,
      RefinementCase<dim-1>::no_refinement,
      RefinementCase<dim-1>::no_refinement
    },
    // case_y1x
    {
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_y,
      RefinementCase<dim-1>::no_refinement
    },
    // case_y2x
    {
      RefinementCase<dim-1>::cut_y,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::no_refinement
    },
    // case_y1x2x
    {
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy
    },
    // case_xy (case_isotropic)
    {
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy,
      RefinementCase<dim-1>::cut_xy
    }
  };

  static const unsigned int
  equivalent_subface_number[internal::SubfaceCase<dim>::case_isotropic+1][GeometryInfo<3>::max_children_per_face]
  =
  {
    // case_none, see above
    {0,1,2,3},
    // case_x
    {0,1,e,e},
    // case_x1y
    {0,2,1,e},
    // case_x2y
    {0,1,3,e},
    // case_x1y2y
    {0,2,1,3},
    // case_y
    {0,1,e,e},
    // case_y1x
    {0,1,1,e},
    // case_y2x
    {0,2,3,e},
    // case_y1x2x
    {0,1,2,3},
    // case_xy (case_isotropic)
    {0,1,2,3}
  };

  // If face-orientation or face_rotation are
  // non-standard, cut_x and cut_y have to be
  // exchanged.
  static const RefinementCase<dim-1> ref_case_permutation[4]
    = {RefinementCase<dim-1>::no_refinement,
       RefinementCase<dim-1>::cut_y,
       RefinementCase<dim-1>::cut_x,
       RefinementCase<dim-1>::cut_xy
      };

  // set a corresponding (equivalent)
  // RefineCase and subface number
  const RefinementCase<dim-1> equ_ref_case=equivalent_refine_case[ref_case][subface_no];
  const unsigned int equ_subface_no=equivalent_subface_number[ref_case][subface_no];
  // make sure, that we got a valid subface and RefineCase
  Assert(equ_ref_case!=RefinementCase<dim>::no_refinement, ExcInternalError());
  Assert(equ_subface_no!=e, ExcInternalError());
  // now, finally respect non-standard faces
  const RefinementCase<dim-1>
  final_ref_case = (face_orientation==face_rotation
                    ?
                    ref_case_permutation[equ_ref_case]
                    :
                    equ_ref_case);

  // what we have now is the number of
  // the subface in the natural
  // orientation of the *face*. what we
  // need to know is the number of the
  // subface concerning the standard face
  // orientation as seen from the *cell*.

  // this mapping is not trivial, but we
  // have done exactly this stuff in the
  // child_cell_on_face function. in
  // order to reduce the amount of code
  // as well as to make maintaining the
  // functionality easier we want to
  // reuse that information. So we note
  // that on the bottom face (face 4) of
  // a hex cell the local x and y
  // coordinates of the face and the cell
  // coincide, thus also the refinement
  // case of the face corresponds to the
  // refinement case of the cell
  // (ignoring cell refinement along the
  // z direction). Using this knowledge
  // we can (ab)use the
  // child_cell_on_face function to do
  // exactly the transformation we are in
  // need of now
  const unsigned int
  final_subface_no = GeometryInfo<dim>::child_cell_on_face(RefinementCase<dim>(final_ref_case),
                                                           4,
                                                           equ_subface_no,
                                                           face_orientation,
                                                           face_flip,
                                                           face_rotation,
                                                           equ_ref_case);

  return (((face_no * total_subfaces_per_face
            + ref_case_offset[final_ref_case-1]
            + final_subface_no)
           + orientation_offset[face_orientation][face_flip][face_rotation]
          )
          * n_quadrature_points);
}



template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_face(const SubQuadrature &quadrature,
                                 const unsigned int face_no)
{
  std::vector<Point<dim> > points(quadrature.size());
  project_to_face(quadrature, face_no, points);
  return Quadrature<dim>(points, quadrature.get_weights());
}


template <int dim>
Quadrature<dim>
QProjector<dim>::project_to_subface(const SubQuadrature       &quadrature,
                                    const unsigned int         face_no,
                                    const unsigned int         subface_no,
                                    const RefinementCase<dim-1> &ref_case)
{
  std::vector<Point<dim> > points(quadrature.size());
  project_to_subface(quadrature, face_no, subface_no, points, ref_case);
  return Quadrature<dim>(points, quadrature.get_weights());
}


// ------------------------------------------------------------ //


template <>
bool
QIterated<1>::uses_both_endpoints (const Quadrature<1> &base_quadrature)
{
  bool at_left = false,
       at_right = false;
  for (unsigned int i=0; i<base_quadrature.size(); ++i)
    {
      if (base_quadrature.point(i) == Point<1>(0.0))
        at_left = true;
      if (base_quadrature.point(i) == Point<1>(1.0))
        at_right = true;
    };

  return (at_left && at_right);
}


// template<>
// void
// QIterated<1>::fill(Quadrature<1>& dst,
//                 const Quadrature<1> &base_quadrature,
//                 const unsigned int   n_copies)
// {
//   Assert (n_copies > 0, ExcZero());
//   Assert (base_quadrature.size() > 0, ExcZero());

//   const unsigned int np =
//     uses_both_endpoints(base_quadrature)
//     ? (base_quadrature.size()-1) * n_copies + 1
//     : base_quadrature.size() * n_copies;

//   dst.quadrature_points.resize(np);
//   dst.weights.resize(np);

//   if (!uses_both_endpoints(base_quadrature))
//                                   // we don't have to skip some
//                                   // points in order to get a
//                                   // reasonable quadrature formula
//     {
//       unsigned int next_point = 0;
//       for (unsigned int copy=0; copy<n_copies; ++copy)
//      for (unsigned int q_point=0; q_point<base_quadrature.size(); ++q_point)
//        {
//          dst.quadrature_points[next_point](0)
//            = (copy + base_quadrature.point(q_point)(0)) / n_copies;
//          dst.weights[next_point]
//            = base_quadrature.weight(q_point) / n_copies;
//          ++next_point;
//        };
//     }
//   else
//                                   // skip doubly available points
//     {
//       unsigned int next_point = 0;

//                                     // first find out the weights of
//                                     // the left and the right boundary
//                                     // points. note that these usually
//                                     // are but need not necessarily be
//                                     // the same
//       double double_point_weight = 0;
//       unsigned int n_end_points = 0;
//       for (unsigned int i=0; i<base_quadrature.size(); ++i)
//                                       // add up the weight if this
//                                       // is an endpoint
//      if ((base_quadrature.point(i)(0) == 0.) ||
//          (base_quadrature.point(i)(0) == 1.))
//        {
//          double_point_weight += base_quadrature.weight(i);
//          ++n_end_points;
//        };
//                                     // scale the weight correctly
//       double_point_weight /= n_copies;

//                                     // make sure the base quadrature formula
//                                     // has only one quadrature point
//                                     // per end point
//       Assert (n_end_points == 2, ExcInvalidQuadratureFormula());


//       for (unsigned int copy=0; copy<n_copies; ++copy)
//      for (unsigned int q_point=0; q_point<base_quadrature.size(); ++q_point)
//        {
//                                           // skip the left point of
//                                           // this copy since we
//                                           // have already entered
//                                           // it the last time
//          if ((copy > 0) &&
//              (base_quadrature.point(q_point)(0) == 0.))
//            continue;

//          dst.quadrature_points[next_point](0)
//            = (copy+base_quadrature.point(q_point)(0)) / n_copies;

//                                           // if this is the
//                                           // rightmost point of one
//                                           // of the non-last
//                                           // copies: give it the
//                                           // double weight
//          if ((copy != n_copies-1) &&
//              (base_quadrature.point(q_point)(0) == 1.))
//            dst.weights[next_point] = double_point_weight;
//          else
//            dst.weights[next_point] = base_quadrature.weight(q_point) /
//                                        n_copies;

//          ++next_point;
//        };
//     };

// #if DEBUG
//   double sum_of_weights = 0;
//   for (unsigned int i=0; i<dst.size(); ++i)
//     sum_of_weights += dst.weight(i);
//   Assert (std::fabs(sum_of_weights-1) < 1e-15,
//        ExcInternalError());
// #endif

// }


template <>
QIterated<0>::QIterated (const Quadrature<1> &,
                         const unsigned int   )
  :
  Quadrature<0>()
{
  Assert (false, ExcNotImplemented());
}



template <>
QIterated<1>::QIterated (const Quadrature<1> &base_quadrature,
                         const unsigned int   n_copies)
  :
  Quadrature<1> (uses_both_endpoints(base_quadrature) ?
                 (base_quadrature.size()-1) * n_copies + 1 :
                 base_quadrature.size() * n_copies)
{
//  fill(*this, base_quadrature, n_copies);
  Assert (base_quadrature.size() > 0, ExcNotInitialized());
  Assert (n_copies > 0, ExcZero());

  if (!uses_both_endpoints(base_quadrature))
    // we don't have to skip some
    // points in order to get a
    // reasonable quadrature formula
    {
      unsigned int next_point = 0;
      for (unsigned int copy=0; copy<n_copies; ++copy)
        for (unsigned int q_point=0; q_point<base_quadrature.size(); ++q_point)
          {
            this->quadrature_points[next_point]
              = Point<1>(base_quadrature.point(q_point)(0) / n_copies
                         +
                         (1.0*copy)/n_copies);
            this->weights[next_point]
              = base_quadrature.weight(q_point) / n_copies;

            ++next_point;
          };
    }
  else
    // skip doubly available points
    {
      unsigned int next_point = 0;

      // first find out the weights of
      // the left and the right boundary
      // points. note that these usually
      // are but need not necessarily be
      // the same
      double double_point_weight = 0;
      unsigned int n_end_points = 0;
      for (unsigned int i=0; i<base_quadrature.size(); ++i)
        // add up the weight if this
        // is an endpoint
        if ((base_quadrature.point(i) == Point<1>(0.0)) ||
            (base_quadrature.point(i) == Point<1>(1.0)))
          {
            double_point_weight += base_quadrature.weight(i);
            ++n_end_points;
          };
      // scale the weight correctly
      double_point_weight /= n_copies;

      // make sure the base quadrature formula
      // has only one quadrature point
      // per end point
      Assert (n_end_points == 2, ExcInvalidQuadratureFormula());


      for (unsigned int copy=0; copy<n_copies; ++copy)
        for (unsigned int q_point=0; q_point<base_quadrature.size(); ++q_point)
          {
            // skip the left point of
            // this copy since we
            // have already entered
            // it the last time
            if ((copy > 0) &&
                (base_quadrature.point(q_point) == Point<1>(0.0)))
              continue;

            this->quadrature_points[next_point]
              = Point<1>(base_quadrature.point(q_point)(0) / n_copies
                         +
                         (1.0*copy)/n_copies);

            // if this is the
            // rightmost point of one
            // of the non-last
            // copies: give it the
            // double weight
            if ((copy != n_copies-1) &&
                (base_quadrature.point(q_point) == Point<1>(1.0)))
              this->weights[next_point] = double_point_weight;
            else
              this->weights[next_point] = base_quadrature.weight(q_point) /
                                          n_copies;

            ++next_point;
          };
    };

#if DEBUG
  double sum_of_weights = 0;
  for (unsigned int i=0; i<this->size(); ++i)
    sum_of_weights += this->weight(i);
  Assert (std::fabs(sum_of_weights-1) < 1e-13,
          ExcInternalError());
#endif
}


// template <int dim>
// void
// QIterated<dim>::fill(Quadrature<dim>&, const Quadrature<1>&, unsigned int)
// {
//   Assert(false, ExcNotImplemented());
// }


// construct higher dimensional quadrature formula by tensor product
// of lower dimensional iterated quadrature formulae
template <int dim>
QIterated<dim>::QIterated (const Quadrature<1> &base_quadrature,
                           const unsigned int   N)
  :
  Quadrature<dim> (QIterated<dim-1>(base_quadrature, N),
                   QIterated<1>(base_quadrature, N))
{}



// explicit instantiations; note: we need them all for all dimensions
template class Quadrature<0>;
template class Quadrature<1>;
template class Quadrature<2>;
template class Quadrature<3>;
template class QAnisotropic<1>;
template class QAnisotropic<2>;
template class QAnisotropic<3>;
template class QIterated<1>;
template class QIterated<2>;
template class QIterated<3>;
template class QProjector<1>;
template class QProjector<2>;
template class QProjector<3>;

DEAL_II_NAMESPACE_CLOSE
