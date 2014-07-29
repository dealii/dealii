/* ---------------------------------------------------------------------
 * $Id$
 *
 * Copyright (C) 2014 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2014
 */

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/base/std_cxx1x/array.h>

#include <fstream>


using namespace dealii;


template <int dim>
class SphereGeometry : public Boundary<dim>
{
public:
    SphereGeometry (const Point<dim> &center);
  virtual
  Point<dim>
  get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const;

  virtual
  Point<dim>
  get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const;
private:

  std_cxx1x::array<double,dim> pull_back (const Point<dim> &p) const;
  Point<dim>                   push_forward (const std_cxx1x::array<double,dim> &preimage) const;

  template <int N>
  static std_cxx1x::array<double,dim> average (const std_cxx1x::array<double,dim> (&array)[N]);

  const Point<dim> center;
};


template <int dim>
SphereGeometry<dim>::SphereGeometry (const Point<dim> &center)
:
center (center)
{}


template <>
std_cxx1x::array<double,2>
SphereGeometry<2>::pull_back (const Point<2> &p) const
{
  const Point<2> relative_point = p - center;

  const double r = relative_point.norm();
  const double phi = std::atan2 (relative_point[1], relative_point[0]);

  std_cxx1x::array<double,2> result;
  result[0] = r;
  result[1] = phi;

  return result;
}


template <>
Point<2>
SphereGeometry<2>::push_forward (const std_cxx1x::array<double,2> &preimage) const
{
  const Point<2> relative_point = preimage[0] * Point<2>(std::cos(preimage[1]), std::sin(preimage[1]));

  return relative_point + center;
}



template <>
std_cxx1x::array<double,3>
SphereGeometry<3>::pull_back (const Point<3> &p) const
{
  const Point<3> relative_point = p - center;

  const double r     = relative_point.norm();
  const double phi   = std::atan2 (relative_point[1], relative_point[0]);
  const double theta = std::atan2 (relative_point[2], std::sqrt (relative_point[0]*relative_point[0] +
                                                                 relative_point[1]*relative_point[1]));

  std_cxx1x::array<double,3> result;
  result[0] = r;
  result[1] = phi;
  result[2] = theta;

  return result;
}


template <>
Point<3>
SphereGeometry<3>::push_forward (const std_cxx1x::array<double,3> &preimage) const
{
  const Point<3> relative_point = (preimage[0] *
                                   Point<3>(std::cos(preimage[1]),
                                            std::sin(preimage[1]),
                                            std::cos(preimage[2])));

  return relative_point + center;
}



template <>
template <int N>
std_cxx1x::array<double,2>
SphereGeometry<2>::average (const std_cxx1x::array<double,2> (&array)[N])
{
  std_cxx1x::array<double,2> result;

  // average the radii first. this is uncritical
  {
    result[0] = 0;
    for (unsigned int i=0; i<N; ++i)
      result[0] += array[i][0];
    result[0] /= N;
  }

  // now do the angle. there, we need to
  // be more careful because the average of 0.9*pi and -0.9*pi should not
  // be zero but pi. to this end, bring everything that is farther than
  // pi away from the first angle we want to average with, back to within pi
  // by adding/subtracting 2*pi
  //
  // we also want to make sure that we exclude from averaging all points
  // that lie at the origin since they have no angle at all
  {
    bool origin_is_one_point = false;

    result[1] = 0;
    for (unsigned int i=0; i<N; ++i)
      if (array[i][0] > 1e-10)
        {
          const double angle = array[i][1];
          if (angle - array[0][1] > numbers::PI)
            result[1] += angle-2*numbers::PI;
          else if (angle - array[0][1] < -numbers::PI)
            result[1] += angle+2*numbers::PI;
          else
            result[1] += angle;
        }
      else
        origin_is_one_point = true;

    if (origin_is_one_point == false)
      result[1] /= N;
    else
      result[1] /= (N-1);
  }

  return result;
}



template <>
template <int N>
std_cxx1x::array<double,3>
SphereGeometry<3>::average (const std_cxx1x::array<double,3> (&array)[N])
{
  std_cxx1x::array<double,3> result;

  // average the radii first. this is uncritical
  {
    result[0] = 0;
    for (unsigned int i=0; i<N; ++i)
      result[0] += array[i][0];
    result[0] /= N;
  }

  // now do the angle along the equatorial direction. do the same as we did
  // in the 2d case
  {
    bool origin_is_one_point = false;

    result[1] = 0;
    for (unsigned int i=0; i<N; ++i)
      if (array[i][0] > 1e-10)
        {
          const double angle = array[i][1];
          if (angle - array[0][1] > numbers::PI)
            result[1] += angle-2*numbers::PI;
          else if (angle - array[0][1] < -numbers::PI)
            result[1] += angle+2*numbers::PI;
          else
            result[1] += angle;
        }
      else
        origin_is_one_point = true;

    if (origin_is_one_point == false)
      result[1] /= N;
    else
      result[1] /= (N-1);
  }

  // finally for the polar angle. the difficulty here is that we have, for
  // example, two points at $(r,\phi,\theta)$ and $(r,\phi+\pi,\pm \pi/2)$, then we want
  // to average these to $(r,\ast,\pi)$ where the equatorial angle does not matter
  {
    //??? not sure how exactly to do this
  }


  return result;
}




template <int dim>
Point<dim>
SphereGeometry<dim>::
get_new_point_on_line (const typename Triangulation<dim>::line_iterator &line) const
{
  std_cxx1x::array<double,dim> preimages[2] = { pull_back (line->vertex(0)),
                                                pull_back (line->vertex(1)) };

  return push_forward(average (preimages));
}


template <int dim>
Point<dim>
SphereGeometry<dim>::
get_new_point_on_quad (const typename Triangulation<dim>::quad_iterator &quad) const
{
  std_cxx1x::array<double,dim> preimages[4] = { pull_back (quad->vertex(0)),
                                                pull_back (quad->vertex(1)),
                                                pull_back (quad->vertex(2)),
                                                pull_back (quad->vertex(3)) };

  return push_forward(average(preimages));
}




template <int dim>
void make_grid ()
{
  Point<dim> center;
  for (unsigned int i=0; i<dim; ++i)
    center[i] = .25;
  const double radius=center.norm();

  SphereGeometry<dim> geometry(center);
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube (triangulation);

  triangulation.refine_global(1);

  for (typename Triangulation<dim>::active_cell_iterator cell=triangulation.begin_active();
      cell!=triangulation.end(); ++cell)
    {
      if (cell->center().distance(center)< radius)
        cell->set_manifold_id(1);

      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->face(f)->center().distance(center)< radius)
          cell->face(f)->set_all_manifold_ids(1);
    }

  triangulation.set_manifold(1,geometry);
  triangulation.refine_global(1);

  const std::string filename = "mesh-" + Utilities::int_to_string(dim) + "d.vtk";
  std::ofstream out (filename.c_str());
  GridOut grid_out;
  grid_out.write_vtk (triangulation, out);
}


  

// @sect3{The main function}

// Finally, the main function. There isn't much to do here, only to call the
// subfunctions.
int main ()
{
  make_grid<2> ();
//  make_grid<3> ();
}
