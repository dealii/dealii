// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
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


// test for GridTools::guess_point_owner

#include "../tests.h"

#include <deal.II/base/point.h>
#include <deal.II/base/bounding_box.h>
#include <deal.II/grid/grid_tools.h>

#include <tuple>

template <int spacedim>
void
test_point_owner(unsigned int n_procs)
{
  // Step 1: creating the vector of bounding boxes
  // We want to simulate n_procs bounding boxes:
  // each process of rank rk has rk + 1 bounding boxes
  // the bounding boxes describe a mesh which is a
  // segment (in dim 1), a rectangle (in dim 2),
  // a parallelogram (in dim 3)
  // and using the center point of each bounding box

  std::vector< std::vector< BoundingBox<spacedim> > >
  global_bboxes(n_procs);

  unsigned int tot_bbox = 0;

  for (unsigned int rk=0; rk<n_procs; ++rk)
    for (unsigned int box=0; box<rk; ++box)
      {
        std::pair<Point<spacedim>,Point<spacedim>> boundaries;
        boundaries.first[0] = tot_bbox;
        boundaries.second[0] = tot_bbox + 1 ;
        for (int i=1; i<spacedim; i++)
          boundaries.second[i] = 1;

        BoundingBox<spacedim> new_box(boundaries);
        global_bboxes[rk].push_back(new_box);
        ++tot_bbox;
      }

  // Creating the vector of center points
  std::vector< Point< spacedim > > points(tot_bbox);
  for (unsigned int i=0; i<tot_bbox; ++i)
    {
      points[i][0] = i + 0.5;
      for (unsigned int j=1; j<spacedim; j++)
        points[i][j] = 0.5;
    }

  // Step 2: testing the function
  auto output_tp = GridTools::guess_point_owner(global_bboxes,points);

  bool test_passed = true;

  if (std::get<2>(output_tp).size()!=0)
    {
      test_passed = false;
      deallog << " Test 1 failed: multiple owners" << std::endl;
    }
  else
    {
      // Check the points are all in the correct rank
      unsigned int tot_pt = 0;
      for (unsigned int rk=0; rk<n_procs; ++rk)
        {
          const auto &rk_points = std::get<0>(output_tp)[rk];
          for (unsigned int box=0; box<rk; ++box)
            {
              if ( std::find( rk_points.begin(), rk_points.end(),tot_pt) == rk_points.end())
                {
                  deallog << "Point " << tot_pt << " not found in rank " << rk << std::endl;
                  test_passed = false;
                }
              else
                ++tot_pt;
            }
        }
    }

  // Step 3: creating random points and adding bounding boxes

  // Creating the random points
  points.clear();

  unsigned int n_points = 2*tot_bbox;
  for (size_t i=0; i<n_points; ++i)
    {
      Point<spacedim> p;
      p[0] = tot_bbox*double(Testing::rand())/RAND_MAX;
      for (unsigned int d=1; d<spacedim; ++d)
        p[d] = double(Testing::rand())/RAND_MAX;
      points.push_back(p);
    }

  // Adding a bounding box to some processes which covers
  // part of the domain already covered by other processes

  for (unsigned int rk=1; rk<n_procs; rk+=2)
    {
      std::pair<Point<spacedim>,Point<spacedim>> boundaries
                                              = global_bboxes[rk][0].get_boundary_points();

      boundaries.first[0] -= 0.5*double(Testing::rand())/RAND_MAX;
      boundaries.second[0] -= 0.5*double(Testing::rand())/RAND_MAX;

      BoundingBox<spacedim> new_box(boundaries);
      global_bboxes[rk].push_back(new_box);
    }

  // Running the function again and checking the output
  output_tp = GridTools::guess_point_owner(global_bboxes,points);

  for (unsigned int rk=0; rk<n_procs ; ++rk)
    {
      for (const auto &pt: std::get<0>(output_tp)[rk])
        {
          bool found = false;
          for (const auto &bbox: global_bboxes[rk])
            if (bbox.point_inside(points[pt]))
              {
                found = true;
                break;
              }

          if ( ! found )
            {
              test_passed = false;
              deallog << "Couldn't find point " << pt << " in process " << rk << std::endl;
            }
          else
            {
              // Check if the point is contained also in one of the maps
              if ( std::get<1>(output_tp).find(pt) == std::get<1>(output_tp).end() &&
                   std::get<2>(output_tp).find(pt) == std::get<2>(output_tp).end())
                {
                  test_passed = false;
                  deallog << "Output maps missing point " << pt << " in process " << rk << std::endl;
                }
            }
        }
    }

  if (test_passed)
    deallog << "TEST PASSED" << std::endl;
  else
    deallog << "TEST FAILED" << std::endl;
}

int
main()
{
  initlog();

  deallog << "Test for GridTools::guess_point_owner " << std::endl;
  deallog << std::endl << "Test for dimension 1" << std::endl;
  for (unsigned int d=1; d<4; ++d)
    {
      deallog << "Simulating " << d << " processes" << std::endl;
      test_point_owner<1>(d);
    }

  deallog << std::endl << "Test for dimension 2" << std::endl;
  for (unsigned int d=2; d<10; ++d)
    {
      deallog << "Simulating " << d << " processes" << std::endl;
      test_point_owner<2>(d);
    }


  deallog << std::endl << "Test for dimension 3" << std::endl;
  for (unsigned int d=1; d<17; d+=2)
    {
      deallog << "Simulating " << d << " processes" << std::endl;
      test_point_owner<3>(d);
    }

}
