// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

#include <deal.II/python/triangulation_wrapper.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

namespace PyDealII
{
  TriangulationWrapper::TriangulationWrapper(const std::string &dimension)
    :
    triangulation(nullptr)
  {
    if ((dimension.compare("2D")==0) || (dimension.compare("2d")==0))
      {
        dim = 2;
        triangulation = new dealii::Triangulation<2>();
      }
    else if ((dimension.compare("3D")==0) || (dimension.compare("3d")==0))
      {
        dim = 3;
        triangulation = new dealii::Triangulation<3>();
      }
    else
      AssertThrow(false, dealii::ExcMessage("Dimension needs to be 2D or 3D"));
  }



  TriangulationWrapper::~TriangulationWrapper()
  {
    if (triangulation != nullptr)
      {
        if (dim == 2)
          {
            // We cannot call delete on a void pointer so cast the void pointer back
            // first.
            dealii::Triangulation<2> *tmp = static_cast<dealii::Triangulation<2>*>(triangulation);
            delete tmp;
          }
        else
          {
            dealii::Triangulation<3> *tmp = static_cast<dealii::Triangulation<3>*>(triangulation);
            delete tmp;
          }
        triangulation = nullptr;
      }
    dim = -1;
  }



  unsigned int TriangulationWrapper::n_active_cells() const
  {
    if (dim == 2)
      {
        return (*static_cast<dealii::Triangulation<2>*>(triangulation)).n_active_cells();
      }
    else
      return (*static_cast<dealii::Triangulation<3>*>(triangulation)).n_active_cells();
  }



  void TriangulationWrapper::generate_hyper_cube(const double left,
                                                 const double right,
                                                 const bool   colorize)
  {
    if (dim == 2)
      {
        dealii::Triangulation<2> *tria =
          static_cast<dealii::Triangulation<2>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::hyper_cube(*tria, left, right, colorize);
      }
    else
      {
        dealii::Triangulation<3> *tria =
          static_cast<dealii::Triangulation<3>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::hyper_cube(*tria, left, right, colorize);
      }
  }



  void TriangulationWrapper::generate_simplex(boost::python::list &vertices)
  {
    // Extract the PointWrapper object from the python list
    std::vector<PointWrapper> wrapped_points(dim);
    for (int i=0; i<dim; ++i)
      wrapped_points[i] = boost::python::extract<PointWrapper>(vertices[i]);

    if (dim == 2)
      {
        // Cast the PointWrapper objects to Point<dim>
        std::vector<dealii::Point<2>> points(dim);
        for (int i=0; i<dim; ++i)
          points[i] = *(static_cast<dealii::Point<2>*>((wrapped_points[i]).get_point()));

        dealii::Triangulation<2> *tria =
          static_cast<dealii::Triangulation<2>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::simplex(*tria, points);
      }
    else
      {
        std::vector<dealii::Point<3>> points(dim);
        for (int i=0; i<dim; ++i)
          points[i] = *(static_cast<dealii::Point<3>*>((wrapped_points[i]).get_point()));

        dealii::Triangulation<3> *tria =
          static_cast<dealii::Triangulation<3>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::simplex(*tria, points);
      }
  }



  void TriangulationWrapper::generate_subdivided_hyper_cube(const unsigned int repetitions,
                                                            const double       left,
                                                            const double       right)
  {
    if (dim == 2)
      {
        dealii::Triangulation<2> *tria =
          static_cast<dealii::Triangulation<2>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::subdivided_hyper_cube(*tria, repetitions, left, right);
      }
    else
      {
        dealii::Triangulation<3> *tria =
          static_cast<dealii::Triangulation<3>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::subdivided_hyper_cube(*tria, repetitions, left, right);
      }
  }



  void TriangulationWrapper::generate_hyper_rectangle(PointWrapper &p1,
                                                      PointWrapper &p2,
                                                      const bool    colorize)
  {
    if (dim == 2)
      {
        // Cast the PointWrapper object to Point<dim>
        dealii::Point<2> point_1 = *(static_cast<dealii::Point<2>*>(p1.get_point()));
        dealii::Point<2> point_2 = *(static_cast<dealii::Point<2>*>(p2.get_point()));

        dealii::Triangulation<2> *tria =
          static_cast<dealii::Triangulation<2>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::hyper_rectangle(*tria, point_1, point_2, colorize);
      }
    else
      {
        dealii::Point<3> point_1 = *(static_cast<dealii::Point<3>*>(p1.get_point()));
        dealii::Point<3> point_2 = *(static_cast<dealii::Point<3>*>(p2.get_point()));

        dealii::Triangulation<3> *tria =
          static_cast<dealii::Triangulation<3>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::hyper_rectangle(*tria, point_1, point_2, colorize);
      }

  }



  void TriangulationWrapper::generate_subdivided_hyper_rectangle(boost::python::list &repetition_list,
      PointWrapper        &p1,
      PointWrapper        &p2,
      const bool           colorize)
  {
    // Extract the repetitions from the python list
    std::vector<unsigned int> repetitions(dim);
    for (int i=0; i<dim; ++i)
      repetitions[i] = boost::python::extract<unsigned int>(repetition_list[i]);

    if (dim == 2)
      {
        // Cast the PointWrapper object to Point<dim>
        dealii::Point<2> point_1 = *(static_cast<dealii::Point<2>*>(p1.get_point()));
        dealii::Point<2> point_2 = *(static_cast<dealii::Point<2>*>(p2.get_point()));

        dealii::Triangulation<2> *tria =
          static_cast<dealii::Triangulation<2>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::subdivided_hyper_rectangle(*tria, repetitions, point_1,
                                                          point_2, colorize);
      }
    else
      {
        dealii::Point<3> point_1 = *(static_cast<dealii::Point<3>*>(p1.get_point()));
        dealii::Point<3> point_2 = *(static_cast<dealii::Point<3>*>(p2.get_point()));

        dealii::Triangulation<3> *tria =
          static_cast<dealii::Triangulation<3>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::subdivided_hyper_rectangle(*tria, repetitions, point_1,
                                                          point_2, colorize);
      }
  }



  void TriangulationWrapper::generate_hyper_ball(PointWrapper &center,
                                                 const double  radius)
  {
    if (dim == 2)
      {
        // Cast the PointWrapper object to Point<dim>
        dealii::Point<2> center_point = *(static_cast<dealii::Point<2>*>(
                                            center.get_point()));

        dealii::Triangulation<2> *tria =
          static_cast<dealii::Triangulation<2>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::hyper_ball(*tria, center_point, radius);
      }
    else
      {
        dealii::Point<3> center_point = *(static_cast<dealii::Point<3>*>(
                                            center.get_point()));

        dealii::Triangulation<3> *tria =
          static_cast<dealii::Triangulation<3>*>(triangulation);
        tria->clear();
        dealii::GridGenerator::hyper_ball(*tria, center_point, radius);
      }
  }



  void TriangulationWrapper::shift(boost::python::list &shift_list)
  {
    if (dim == 2)
      {
        // Extract the shift vector from the python list
        dealii::Tensor<1,2> shift_vector;
        for (int i=0; i<dim; ++i)
          shift_vector[i] = boost::python::extract<double>(shift_list[i]);

        dealii::Triangulation<2> *tria =
          static_cast<dealii::Triangulation<2>*>(triangulation);
        dealii::GridTools::shift(shift_vector, *tria);
      }
    else
      {
        dealii::Tensor<1,3> shift_vector;
        for (int i=0; i<dim; ++i)
          shift_vector[i] = boost::python::extract<double>(shift_list[i]);

        dealii::Triangulation<3> *tria =
          static_cast<dealii::Triangulation<3>*>(triangulation);
        dealii::GridTools::shift(shift_vector, *tria);
      }
  }



  void TriangulationWrapper::merge_triangulations(TriangulationWrapper &triangulation_1,
                                                  TriangulationWrapper &triangulation_2)
  {
    AssertThrow(triangulation_1.get_dim() == triangulation_2.get_dim(),
                dealii::ExcMessage("Triangulation_1 and Triangulation_2 should have the same dimension."));

    if (triangulation_1.get_dim() == 2)
      {
        dealii::Triangulation<2> *tria =
          static_cast<dealii::Triangulation<2>*>(triangulation);
        tria->clear();
        dealii::Triangulation<2> *tria_1 =
          static_cast<dealii::Triangulation<2>*>(triangulation_1.get_triangulation());
        dealii::Triangulation<2> *tria_2 =
          static_cast<dealii::Triangulation<2>*>(triangulation_2.get_triangulation());
        dealii::GridGenerator::merge_triangulations(*tria_1, *tria_2, *tria);
        // We need to reassign tria to triangulation because tria was cleared
        // inside merge_triangulations.
        triangulation = tria;
      }
    else
      {
        dealii::Triangulation<3> *tria =
          static_cast<dealii::Triangulation<3>*>(triangulation);
        tria->clear();
        dealii::Triangulation<3> *tria_1 =
          static_cast<dealii::Triangulation<3>*>(triangulation_1.get_triangulation());
        dealii::Triangulation<3> *tria_2 =
          static_cast<dealii::Triangulation<3>*>(triangulation_2.get_triangulation());
        dealii::GridGenerator::merge_triangulations(*tria_1, *tria_2, *tria);
        triangulation = tria;
      }
  }



  void TriangulationWrapper::save(const std::string &filename) const
  {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);

    if (dim == 2)
      {
        dealii::Triangulation<2> *tria =
          static_cast<dealii::Triangulation<2>*>(triangulation);

        oa << *tria;
      }
    else
      {
        dealii::Triangulation<3> *tria =
          static_cast<dealii::Triangulation<3>*>(triangulation);

        oa << *tria;
      }
  }



  void TriangulationWrapper::load(const std::string &filename)
  {
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);

    if (dim == 2)
      {
        dealii::Triangulation<2> *tria =
          static_cast<dealii::Triangulation<2>*>(triangulation);

        ia >> *tria;
      }
    else
      {
        dealii::Triangulation<3> *tria =
          static_cast<dealii::Triangulation<3>*>(triangulation);

        ia >> *tria;
      }
  }



  void TriangulationWrapper::refine_global(const unsigned int n)
  {
    if (dim == 2)
      {
        dealii::Triangulation<2> *tria =
          static_cast<dealii::Triangulation<2>*>(triangulation);
        tria->refine_global(n);
      }
    else
      {
        dealii::Triangulation<3> *tria =
          static_cast<dealii::Triangulation<3>*>(triangulation);
        tria->refine_global(n);
      }
  }
}
