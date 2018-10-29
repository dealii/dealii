// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <triangulation_wrapper.h>
#include <cell_accessor_wrapper.h>
#include <deal.II/base/types.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <fstream>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  namespace internal
  {
    template <int dim, int spacedim>
    void generate_hyper_cube(const double left,
                             const double right,
                             const bool   colorize,
                             void        *triangulation)
    {
      Triangulation<dim,spacedim> *tria =
        static_cast<Triangulation<dim,spacedim>*>(triangulation);
      tria->clear();
      GridGenerator::hyper_cube(*tria, left, right, colorize);
    }



    template <int dim>
    void generate_simplex(std::vector<PointWrapper> &wrapped_points,
                          void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      std::vector<Point<dim>> points(dim+1);
      for (int i=0; i<dim+1; ++i)
        points[i] = *(static_cast<Point<dim>*>((wrapped_points[i]).get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::simplex(*tria, points);
    }



    template <int dim, int spacedim>
    void generate_subdivided_hyper_cube(const unsigned int repetitions,
                                        const double       left,
                                        const double       right,
                                        void              *triangulation)
    {
      Triangulation<dim,spacedim> *tria =
        static_cast<Triangulation<dim,spacedim>*>(triangulation);
      tria->clear();
      GridGenerator::subdivided_hyper_cube(*tria, repetitions, left, right);
    }



    template <int dim, int spacedim>
    void generate_hyper_rectangle(PointWrapper &p1,
                                  PointWrapper &p2,
                                  const bool    colorize,
                                  void         *triangulation)
    {
      AssertThrow(p1.get_dim() == dim,
                  ExcMessage("Dimension of p1 is not the same as the dimension of the Triangulation."));
      AssertThrow(p2.get_dim() == dim,
                  ExcMessage("Dimension of p2 is not the same as the dimension of the Triangulation."));
      // Cast the PointWrapper object to Point<dim>
      Point<dim> point_1 = *(static_cast<Point<dim>*>(p1.get_point()));
      Point<dim> point_2 = *(static_cast<Point<dim>*>(p2.get_point()));

      Triangulation<dim,spacedim> *tria =
        static_cast<Triangulation<dim,spacedim>*>(triangulation);
      tria->clear();
      GridGenerator::hyper_rectangle(*tria, point_1, point_2, colorize);
    }



    template <int dim, int spacedim>
    void generate_subdivided_hyper_rectangle(const std::vector<unsigned int> &repetitions,
                                             PointWrapper                    &p1,
                                             PointWrapper                    &p2,
                                             const bool                       colorize,
                                             void                            *triangulation)
    {
      AssertThrow(p1.get_dim() == dim,
                  ExcMessage("Dimension of p1 is not the same as the dimension of the Triangulation."));
      AssertThrow(p2.get_dim() == dim,
                  ExcMessage("Dimension of p2 is not the same as the dimension of the Triangulation."));
      // Cast the PointWrapper object to Point<dim>
      Point<dim> point_1 = *(static_cast<Point<dim>*>(p1.get_point()));
      Point<dim> point_2 = *(static_cast<Point<dim>*>(p2.get_point()));

      Triangulation<dim,spacedim> *tria =
        static_cast<Triangulation<dim,spacedim>*>(triangulation);
      tria->clear();
      GridGenerator::subdivided_hyper_rectangle(*tria, repetitions, point_1,
                                                point_2, colorize);
    }



    template <int dim>
    void generate_subdivided_steps_hyper_rectangle(const std::vector<std::vector<double>> &step_sizes,
                                                   PointWrapper                           &p1,
                                                   PointWrapper                           &p2,
                                                   const bool                              colorize,
                                                   void                                   *triangulation)
    {
      AssertThrow(p1.get_dim() == dim,
                  ExcMessage("Dimension of p1 is not the same as the dimension of the Triangulation."));
      AssertThrow(p2.get_dim() == dim,
                  ExcMessage("Dimension of p2 is not the same as the dimension of the Triangulation."));
      // Cast the PointWrapper object to Point<dim>
      Point<dim> point_1 = *(static_cast<Point<dim>*>(p1.get_point()));
      Point<dim> point_2 = *(static_cast<Point<dim>*>(p2.get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::subdivided_hyper_rectangle(*tria, step_sizes, point_1,
                                                point_2, colorize);
    }



    template <int dim>
    void generate_subdivided_material_hyper_rectangle(const std::vector<std::vector<double>> &spacing,
                                        PointWrapper                          &p,
                                        const Table<dim,types::material_id>   &material_ids,
                                        const bool                             colorize,
                                        void                                  *triangulation)
    {
      AssertThrow(p.get_dim() == dim,
                  ExcMessage("Dimension of p is not the same as the dimension of the Triangulation."));
      // Cast the PointWrapper object to Point<dim>
      Point<dim> point = *(static_cast<Point<dim>*>(p.get_point()));
      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::subdivided_hyper_rectangle(*tria, spacing, point,
                                                material_ids, colorize);
    }



    template <int dim, int spacedim>
    void generate_cheese(const std::vector<unsigned int> &holes, 
                         void                            *triangulation)
    {
      Triangulation<dim,spacedim> *tria =
        static_cast<Triangulation<dim,spacedim>*>(triangulation);
      tria->clear();
      GridGenerator::cheese(*tria, holes);
    }



    template <int dim>
    void generate_general_cell(std::vector<PointWrapper> &wrapped_points,
                       const bool                 colorize,
                       void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      const unsigned int size = wrapped_points.size();
      std::vector<Point<dim>> points(size);
      for (unsigned int i=0; i<size; ++i)
        points[i] = *(static_cast<Point<dim>*>((wrapped_points[i]).get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::general_cell(*tria, points, colorize);
    }



    template <int dim>
    void generate_parallelogram(std::vector<PointWrapper> &wrapped_points,
                                const bool                 colorize,
                                void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      Point<dim> points[dim];
      for (unsigned int i=0; i<dim; ++i)
        points[i] = *(static_cast<Point<dim>*>((wrapped_points[i]).get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::parallelogram(*tria, points, colorize);
    }



    template <int dim>
    void generate_parallelepiped(std::vector<PointWrapper> &wrapped_points,
                                const bool                 colorize,
                                void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      Point<dim> points[dim];
      for (unsigned int i=0; i<dim; ++i)
        points[i] = *(static_cast<Point<dim>*>((wrapped_points[i]).get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::parallelepiped(*tria, points, colorize);
    }



    template <int dim>
    void generate_fixed_subdivided_parallelepiped(unsigned int               n_subdivisions,
                                                  std::vector<PointWrapper> &wrapped_points, 
                                                  const bool                 colorize,
                                                  void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      Point<dim> points[dim];
      for (unsigned int i=0; i<dim; ++i)
        points[i] = *(static_cast<Point<dim>*>((wrapped_points[i]).get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::subdivided_parallelepiped(*tria, n_subdivisions, points, colorize);
    }



    template <int dim>
    void generate_varying_subdivided_parallelepiped(std::vector<unsigned int> &n_subdivisions,
                                                  std::vector<PointWrapper> &wrapped_points, 
                                                  const bool                 colorize,
                                                  void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      Point<dim> points[dim];
      unsigned int subdivisions[dim];
      for (unsigned int i=0; i<dim; ++i)
      {
        points[i] = *(static_cast<Point<dim>*>((wrapped_points[i]).get_point()));
        subdivisions[i] = n_subdivisions[i];
      }

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::subdivided_parallelepiped(*tria, subdivisions, points, colorize);
    }


    
    template <int dim>
    void generate_enclosed_hyper_cube(const double  left, 
                                      const double  right, 
                                      const double  thickness, 
                                      const double  colorize,
                                      void         *triangulation)
    {
      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::enclosed_hyper_cube(*tria, left, right, thickness, colorize);
    }
     


    template <int dim>
    void generate_hyper_ball(PointWrapper &center,
                             const double  radius,
                             void         *triangulation)
    {
      // Cast the PointWrapper object to Point<dim>
      Point<dim> center_point = *(static_cast<Point<dim>*>(
                                    center.get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::hyper_ball(*tria, center_point, radius);
    }



    template <int dim, int spacedim>
    void generate_hyper_sphere(PointWrapper &center, 
                               const double  radius, 
                               void         *triangulation)
    {
      // Cast the PointWrapper object to Point<dim>
      Point<spacedim> center_point = *(static_cast<Point<spacedim>*>(
                                    center.get_point()));

      Triangulation<dim,spacedim> *tria =
        static_cast<Triangulation<dim,spacedim>*>(triangulation);
      tria->clear();
      GridGenerator::hyper_sphere(*tria, center_point, radius);
    }


    
    template <int dim>
    void generate_quarter_hyper_ball(PointWrapper &center, 
                                     const double  radius, 
                                     void         *triangulation)
    {
      // Cast the PointWrapper object to Point<dim>
      Point<dim> center_point = *(static_cast<Point<dim>*>(
                                    center.get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::quarter_hyper_ball(*tria, center_point, radius);
    }


    
    template <int dim>
    void generate_half_hyper_ball(PointWrapper &center, 
                                     const double  radius, 
                                     void         *triangulation)
    {
      // Cast the PointWrapper object to Point<dim>
      Point<dim> center_point = *(static_cast<Point<dim>*>(
                                    center.get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim>*>(triangulation);
      tria->clear();
      GridGenerator::half_hyper_ball(*tria, center_point, radius);
    }



    template <int dim, int spacedim>
    void shift(boost::python::list &shift_list,
               void                *triangulation)
    {
      // Extract the shift vector from the python list
      Tensor<1,spacedim> shift_vector;
      for (int i=0; i<spacedim; ++i)
        shift_vector[i] = boost::python::extract<double>(shift_list[i]);

      Triangulation<dim,spacedim> *tria =
        static_cast<Triangulation<dim,spacedim>*>(triangulation);
      GridTools::shift(shift_vector, *tria);
    }



    template <int dim, int spacedim>
    void merge_triangulations(TriangulationWrapper &triangulation_1,
                              TriangulationWrapper &triangulation_2,
                              void                 *triangulation)
    {
      Triangulation<dim,spacedim> *tria =
        static_cast<Triangulation<dim,spacedim>*>(triangulation);
      tria->clear();
      Triangulation<dim,spacedim> *tria_1 =
        static_cast<Triangulation<dim,spacedim>*>(triangulation_1.get_triangulation());
      Triangulation<dim,spacedim> *tria_2 =
        static_cast<Triangulation<dim,spacedim>*>(triangulation_2.get_triangulation());
      GridGenerator::merge_triangulations(*tria_1, *tria_2, *tria);
      // We need to reassign tria to triangulation because tria was cleared
      // inside merge_triangulations.
      triangulation = tria;
    }



    template <int dim, int spacedim_1, int spacedim_2>
    void flatten_triangulation(void *triangulation, TriangulationWrapper &tria_out)
    {
      Triangulation<dim,spacedim_1> *tria =
        static_cast<Triangulation<dim,spacedim_1>*>(triangulation);
      Triangulation<dim,spacedim_2> *tria_2 =
        static_cast<Triangulation<dim,spacedim_2>*>(tria_out.get_triangulation());
      GridGenerator::flatten_triangulation(*tria, *tria_2);
    }



    template <int dim, int spacedim>
    boost::python::list active_cells(TriangulationWrapper &triangulation_wrapper)
    {
      Triangulation<dim,spacedim> *tria =
        static_cast<Triangulation<dim,spacedim>*>(
          triangulation_wrapper.get_triangulation());
      boost::python::list cells;
      for (auto cell : tria->active_cell_iterators())
        cells.append(CellAccessorWrapper(triangulation_wrapper, cell->level(),
                                         cell->index()));

      return cells;
    }



    template <int dim, int spacedim>
    void write(const std::string &filename,
               const std::string &format,
               const void        *triangulation)
    {
      const Triangulation<dim,spacedim> *tria =
        static_cast<const Triangulation<dim,spacedim>*>(triangulation);

      GridOut::OutputFormat output_format;
      if (format.compare("dx") == 0)
        output_format = GridOut::OutputFormat::dx;
      else if (format.compare("gnuplot") == 0)
        output_format = GridOut::OutputFormat::gnuplot;
      else if (format.compare("eps") == 0)
        output_format = GridOut::OutputFormat::eps;
      else if (format.compare("ucd") == 0)
        output_format = GridOut::OutputFormat::ucd;
      else if (format.compare("xfig") == 0)
        output_format = GridOut::OutputFormat::xfig;
      else if (format.compare("msh") == 0)
        output_format = GridOut::OutputFormat::msh;
      else if (format.compare("svg") == 0)
        output_format = GridOut::OutputFormat::svg;
      else if (format.compare("mathgl") == 0)
        output_format = GridOut::OutputFormat::mathgl;
      else if (format.compare("vtk") == 0)
        output_format = GridOut::OutputFormat::vtk;
      else if (format.compare("vtu") == 0)
        output_format = GridOut::OutputFormat::vtu;
      else
        output_format = GridOut::OutputFormat::none;

      GridOut mesh_writer;
      std::ofstream ofs(filename);
      mesh_writer.write(*tria, ofs, output_format);
      ofs.close();
    }
  }



  TriangulationWrapper::TriangulationWrapper(const std::string &dimension)
  {
    if ((dimension.compare("2D")==0) || (dimension.compare("2d")==0))
      setup("2D", "2D");
    else if ((dimension.compare("3D")==0) || (dimension.compare("3d")==0))
      setup("3D", "3D");
    else
      AssertThrow(false, ExcMessage("Dimension needs to be 2D or 3D"));
  }



  TriangulationWrapper::TriangulationWrapper(const std::string &dimension,
                                             const std::string &spacedimension)
  {
    setup(dimension, spacedimension);
  }



  TriangulationWrapper::~TriangulationWrapper()
  {
    if (triangulation != nullptr)
      {
        if (dim == 2)
          {
            if (spacedim == 2)
              {
                // We cannot call delete on a void pointer so cast the void pointer back
                // first.
                Triangulation<2,2> *tmp =
                  static_cast<Triangulation<2,2>*>(triangulation);
                delete tmp;
              }
            else
              {
                Triangulation<2,3> *tmp =
                  static_cast<Triangulation<2,3>*>(triangulation);
                delete tmp;
              }
          }
        else
          {
            Triangulation<3,3> *tmp = static_cast<Triangulation<3,3>*>(triangulation);
            delete tmp;
          }
        triangulation = nullptr;
      }
    dim = -1;
  }



  unsigned int TriangulationWrapper::n_active_cells() const
  {
    if ((dim == 2) && (spacedim == 2))
      return (*static_cast<Triangulation<2,2>*>(triangulation)).n_active_cells();
    else if ((dim == 2) && (spacedim == 3))
      return (*static_cast<Triangulation<2,3>*>(triangulation)).n_active_cells();
    else
      return (*static_cast<Triangulation<3,3>*>(triangulation)).n_active_cells();
  }



  void TriangulationWrapper::generate_hyper_cube(const double left,
                                                 const double right,
                                                 const bool   colorize)
  {
    if ((dim == 2) && (spacedim == 2))
      internal::generate_hyper_cube<2,2>(left, right, colorize, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::generate_hyper_cube<2,3>(left, right, colorize, triangulation);
    else
      internal::generate_hyper_cube<3,3>(left, right, colorize, triangulation);
  }



  void TriangulationWrapper::generate_simplex(boost::python::list &vertices)
  {
    AssertThrow(boost::python::len(vertices) == dim+1,
                ExcMessage("The number of vertices should be equal to dim+1."));
    AssertThrow(dim == spacedim,
                ExcMessage("This function is only implemented for dim equal to spacedim."));
    // Extract the PointWrapper object from the python list
    std::vector<PointWrapper> wrapped_points(dim+1);
    for (int i=0; i<dim+1; ++i)
      {
        wrapped_points[i] = boost::python::extract<PointWrapper>(vertices[i]);
        AssertThrow(wrapped_points[i].get_dim() == dim,
                    ExcMessage("Point of wrong dimension."));
      }

    if (dim == 2)
      internal::generate_simplex<2>(wrapped_points, triangulation);
    else
      internal::generate_simplex<3>(wrapped_points, triangulation);
  }



  void TriangulationWrapper::generate_subdivided_hyper_cube(const unsigned int repetitions,
                                                            const double       left,
                                                            const double       right)
  {
    if ((dim == 2) && (spacedim == 2))
      internal::generate_subdivided_hyper_cube<2,2>(repetitions, left, right, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::generate_subdivided_hyper_cube<2,3>(repetitions, left, right, triangulation);
    else
      internal::generate_subdivided_hyper_cube<3,3>(repetitions, left, right, triangulation);
  }



  void TriangulationWrapper::generate_hyper_rectangle(PointWrapper &p1,
                                                      PointWrapper &p2,
                                                      const bool    colorize)
  {
    if ((dim == 2) && (spacedim == 2))
      internal::generate_hyper_rectangle<2,2>(p1, p2, colorize, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::generate_hyper_rectangle<2,3>(p1, p2, colorize, triangulation);
    else
      internal::generate_hyper_rectangle<3,3>(p1, p2, colorize, triangulation);
  }



  void TriangulationWrapper::generate_subdivided_hyper_rectangle(boost::python::list &repetition_list,
      PointWrapper        &p1,
      PointWrapper        &p2,
      const bool           colorize)
  {
    AssertThrow(boost::python::len(repetition_list) == dim,
                ExcMessage("The list of repetitions must have the same length as the number of dimension."));

    // Extract the repetitions from the python list
    std::vector<unsigned int> repetitions(dim);
    for (int i=0; i<dim; ++i)
      repetitions[i] = boost::python::extract<unsigned int>(repetition_list[i]);

    if ((dim == 2) && (spacedim == 2))
      internal::generate_subdivided_hyper_rectangle<2,2>(repetitions, p1, p2, colorize,
                                                         triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::generate_subdivided_hyper_rectangle<2,3>(repetitions, p1, p2, colorize,
                                                         triangulation);
    else
      internal::generate_subdivided_hyper_rectangle<3,3>(repetitions, p1, p2, colorize,
                                                         triangulation);
  }



  void TriangulationWrapper::generate_subdivided_steps_hyper_rectangle(boost::python::list &step_sizes_list, 
                                                 PointWrapper        &p1,
                                                 PointWrapper        &p2,
                                                 const bool           colorize)
  {
    AssertThrow(spacedim == dim,
                ExcMessage("This function is only implemented for dim equal to spacedim."));
    AssertThrow(boost::python::len(step_sizes_list) == dim,
                ExcMessage("The list of step_sizes must have the same length as the number of dimension."));

    // Extract the step sizes from the python list
    std::vector<std::vector<double>> step_sizes(dim);
    for (int i=0; i<dim; ++i)
    {
      step_sizes[i].resize(boost::python::len(step_sizes_list[i]));
      for (unsigned int j=0; j<step_sizes[i].size(); ++j)
        step_sizes[i][j] = boost::python::extract<double>(step_sizes_list[i][j]);
    }

    if (dim == 2)
      internal::generate_subdivided_steps_hyper_rectangle<2>(step_sizes,
          p1, p2, colorize, triangulation);
    else
      internal::generate_subdivided_steps_hyper_rectangle<3>(step_sizes,
          p1, p2, colorize, triangulation);
  }



  void TriangulationWrapper::generate_subdivided_material_hyper_rectangle(boost::python::list &spacing_list,
                                           PointWrapper        &p,
                                           boost::python::list &material_id_list,
                                           const bool           colorize)
  {
    AssertThrow(spacedim == dim,
                ExcMessage("This function is only implemented for dim equal to spacedim."));
    AssertThrow(boost::python::len(spacing_list) == dim,
                ExcMessage("The list of spacing must have the same length as the number of dimension."));

    // Extract the spacing and the material ID from the python list
    std::vector<std::vector<double>> spacing(dim);
    for (int i=0; i<dim; ++i)
    {
      spacing[i].resize(boost::python::len(spacing_list[i]));
      for (unsigned int j=0; j<spacing[i].size(); ++j)
        spacing[i][j] = boost::python::extract<double>(spacing_list[i][j]);
    }
    if (dim == 2)
    {
      const unsigned int index_0 = boost::python::len(material_id_list);
      const unsigned int index_1 = boost::python::len(material_id_list[0]);
      Table<2, types::material_id> material_ids(index_0, index_1);
      for (unsigned int i=0; i<index_0; ++i)
        for (unsigned int j=0; j<index_1; ++j)
          // We cannot use extract<types::material_id> because boost will throw
          // an exception if we try to extract -1
          material_ids[i][j] = boost::python::extract<int>(material_id_list[i][j]);

        internal::generate_subdivided_material_hyper_rectangle<2>(spacing, p, 
            material_ids, colorize, triangulation);
    }
    else 
    {
      const unsigned int index_0 = boost::python::len(material_id_list);
      const unsigned int index_1 = boost::python::len(material_id_list[0]);
      const unsigned int index_2 = boost::python::len(material_id_list[0][0]);
      Table<3, types::material_id> material_ids(index_0, index_1, index_2);
      for (unsigned int i=0; i<index_0; ++i)
        for (unsigned int j=0; j<index_1; ++j)
        for (unsigned int k=0; k<index_2; ++k)
          material_ids[i][j][k] = boost::python::extract<int>(material_id_list[i][j][k]);
      internal::generate_subdivided_material_hyper_rectangle<3>(spacing, p, 
          material_ids, colorize, triangulation);
    }
  }



  void TriangulationWrapper::generate_cheese(boost::python::list &holes_list)
  {
    const unsigned int size = boost::python::len(holes_list);
    std::vector<unsigned int> holes(size);
    for (unsigned int i=0; i<size; ++i)
      holes[i] = boost::python::extract<unsigned int>(holes_list[i]);

    if ((dim == 2) && (spacedim == 2))
      internal::generate_cheese<2,2>(holes, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::generate_cheese<2,3>(holes, triangulation);
    else
      internal::generate_cheese<3,3>(holes, triangulation);
  }



  void TriangulationWrapper::generate_general_cell(boost::python::list &vertices, 
                                                   const bool colorize)
  {
    AssertThrow(spacedim == dim,
                ExcMessage("This function is only implementd for dim equal to spacedim."));
    // Extract the PointWrapper object from the python list
    const int size = boost::python::len(vertices);
    AssertThrow(size > 0, ExcMessage("The vertices list is empty."));
    std::vector<PointWrapper> wrapped_points(size);
    for (int i=0; i<size; ++i)
        wrapped_points[i] = boost::python::extract<PointWrapper>(vertices[i]);
    if (dim == 2)
      internal::generate_general_cell<2>(wrapped_points, colorize, triangulation);
    else
      internal::generate_general_cell<3>(wrapped_points, colorize, triangulation);
  }



  void TriangulationWrapper::generate_parallelogram(boost::python::list &corners,
                                                     const bool colorize)
  {
    AssertThrow(spacedim == dim,
                ExcMessage("This function is only implemented for dim equal to spacedim."));
    // Extract the PointWrapper object from the python list
    AssertThrow(boost::python::len(corners) == dim,
                ExcMessage("The list of corners must have the same length as the number of dimension."));
    std::vector<PointWrapper> wrapped_points(dim);
    for (int i=0; i<dim; ++i)
        wrapped_points[i] = boost::python::extract<PointWrapper>(corners[i]);
    if (dim == 2)
      internal::generate_parallelogram<2>(wrapped_points, colorize, triangulation);
    else
      internal::generate_parallelogram<3>(wrapped_points, colorize, triangulation);
  }



  void TriangulationWrapper::generate_parallelepiped(boost::python::list &corners,
                                                     const bool colorize)
  {
    AssertThrow(spacedim == dim,
                ExcMessage("This function is only implemented for dim equal to spacedim."));
    // Extract the PointWrapper object from the python list
    AssertThrow(boost::python::len(corners) == dim,
                ExcMessage("The list of corners must have the same length as the number of dimension."));
    std::vector<PointWrapper> wrapped_points(dim);
    for (int i=0; i<dim; ++i)
        wrapped_points[i] = boost::python::extract<PointWrapper>(corners[i]);
    if (dim == 2)
      internal::generate_parallelepiped<2>(wrapped_points, colorize, triangulation);
    else
      internal::generate_parallelepiped<3>(wrapped_points, colorize, triangulation);
  }



  void TriangulationWrapper::generate_fixed_subdivided_parallelepiped(
    const unsigned int   n_subdivisions, 
    boost::python::list &corners,
    const bool           colorize)
  {
    AssertThrow(spacedim == dim,
                ExcMessage("This function is only implemented for dim equal to spacedim."));
    // Extract the PointWrapper object from the python list
    AssertThrow(boost::python::len(corners) == dim,
                ExcMessage("The list of corners must have the same length as the number of dimension."));
    std::vector<PointWrapper> wrapped_points(dim);
    for (int i=0; i<dim; ++i)
        wrapped_points[i] = boost::python::extract<PointWrapper>(corners[i]);
    if ((dim == 2) && (spacedim == 2))
      internal::generate_fixed_subdivided_parallelepiped<2>(n_subdivisions, 
                                          wrapped_points, colorize, triangulation);
    else
      internal::generate_fixed_subdivided_parallelepiped<3>(n_subdivisions, 
                                          wrapped_points, colorize, triangulation);
  }



  void TriangulationWrapper::generate_varying_subdivided_parallelepiped(
    boost::python::list &n_subdivisions,
    boost::python::list &corners,
    const bool           colorize)
  {
    AssertThrow(spacedim == dim,
                ExcMessage("This function is only implemented for dim equal to spacedim."));
    // Extract the subdivisions from the python list
    AssertThrow(boost::python::len(n_subdivisions) == dim,
                ExcMessage("The list of subdivisions must have the same length as the number of dimension."));
    std::vector<unsigned int> subdivisions(dim);
    for (int i=0; i<dim; ++i)
        subdivisions[i] = boost::python::extract<unsigned int>(n_subdivisions[i]);
    // Extract the PointWrapper object from the python list
    AssertThrow(boost::python::len(corners) == dim,
                ExcMessage("The list of corners must have the same length as the number of dimension."));
    std::vector<PointWrapper> wrapped_points(dim);
    for (int i=0; i<dim; ++i)
        wrapped_points[i] = boost::python::extract<PointWrapper>(corners[i]);
    if (dim == 2)
      internal::generate_varying_subdivided_parallelepiped<2>(subdivisions, 
                                          wrapped_points, colorize, triangulation);
    else
      internal::generate_varying_subdivided_parallelepiped<3>(subdivisions, 
                                          wrapped_points, colorize, triangulation);
  }



  void TriangulationWrapper::generate_enclosed_hyper_cube(const double left,
                                                          const double right,
                                                          const double thickness,
                                                          const bool colorize)
  {
    AssertThrow(spacedim == dim,
                ExcMessage("This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_enclosed_hyper_cube<2>(left, right, thickness, colorize, triangulation);
    else
      internal::generate_enclosed_hyper_cube<3>(left, right, thickness, colorize, triangulation);
  }


    
  void TriangulationWrapper::generate_hyper_ball(PointWrapper &center,
                                                 const double  radius)
  {
    AssertThrow(dim == spacedim,
                ExcMessage("This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_hyper_ball<2>(center, radius, triangulation);
    else
      internal::generate_hyper_ball<3>(center, radius, triangulation);
  }


  void TriangulationWrapper::generate_hyper_sphere(PointWrapper &center,
                                                   const double  radius)
  {
    AssertThrow(spacedim == dim+1,
                ExcMessage("This function is only implemented for spacedim equal to dim+1."));
    internal::generate_hyper_sphere<2,3>(center, radius, triangulation);
  }



  void TriangulationWrapper::generate_quarter_hyper_ball(PointWrapper &center,
                                                         const double  radius)
  {
    AssertThrow(dim == spacedim,
                ExcMessage("This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_quarter_hyper_ball<2>(center, radius, triangulation);
    else
      internal::generate_quarter_hyper_ball<3>(center, radius, triangulation);
  }


  void TriangulationWrapper::generate_half_hyper_ball(PointWrapper &center,
                                                 const double  radius)
  {
    AssertThrow(dim == spacedim,
                ExcMessage("This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_half_hyper_ball<2>(center, radius, triangulation);
    else
      internal::generate_half_hyper_ball<3>(center, radius, triangulation);
  }



  void TriangulationWrapper::shift(boost::python::list &shift_list)
  {
    AssertThrow(boost::python::len(shift_list) == spacedim,
                ExcMessage("Size of the shift vector is not equal to spacedim."));
    if ((dim == 2) && (spacedim == 2))
      internal::shift<2,2>(shift_list, triangulation);
    else if ((dim ==2) && (spacedim == 3))
      internal::shift<2,3>(shift_list, triangulation);
    else
      internal::shift<3,3>(shift_list, triangulation);
  }



  void TriangulationWrapper::merge_triangulations(TriangulationWrapper &triangulation_1,
                                                  TriangulationWrapper &triangulation_2)
  {
    AssertThrow(triangulation_1.get_dim() == triangulation_2.get_dim(),
                ExcMessage("Triangulation_1 and Triangulation_2 should have the same dimension."));
    AssertThrow(dim == triangulation_2.get_dim(),
                ExcMessage("Triangulation and Triangulation_2 should have the same dimension."));
    AssertThrow(triangulation_1.get_spacedim() == triangulation_2.get_spacedim(),
                ExcMessage("Triangulation_1 and Triangulation_2 should have the same space dimension."));
    AssertThrow(spacedim == triangulation_2.get_spacedim(),
                ExcMessage("Triangulation and Triangulation_2 should have the same space dimension."));

    if ((triangulation_1.get_dim() == 2) && (triangulation_1.get_spacedim() == 2))
      internal::merge_triangulations<2,2>(triangulation_1, triangulation_2, triangulation);
    else if ((triangulation_1.get_dim() == 2) && (triangulation_1.get_spacedim() == 3))
      internal::merge_triangulations<2,3>(triangulation_1, triangulation_2, triangulation);
    else
      internal::merge_triangulations<3,3>(triangulation_1, triangulation_2, triangulation);
  }



  void TriangulationWrapper::flatten_triangulation(TriangulationWrapper &tria_out)
  {
    AssertThrow(dim == tria_out.get_dim(),
                ExcMessage("The Triangulation and tria_out should have the same dimension."));
    AssertThrow(spacedim >= tria_out.get_spacedim(),
                ExcMessage("The Triangulation should have a spacedim greater or equal "
                           "to the spacedim of tria_out."));
    int spacedim_out = tria_out.get_spacedim();
    if ((dim == 2) && (spacedim == 2) && (spacedim_out == 2))
      internal::flatten_triangulation<2,2,2>(triangulation, tria_out);
    else if ((dim == 2) && (spacedim == 3) && (spacedim_out == 2))
      internal::flatten_triangulation<2,3,2>(triangulation, tria_out);
    else
      internal::flatten_triangulation<3,3,3>(triangulation, tria_out);
  }



  void TriangulationWrapper::refine_global(const unsigned int n)
  {
    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2,2> *tria =
          static_cast<Triangulation<2,2>*>(triangulation);
        tria->refine_global(n);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        Triangulation<2,3> *tria =
          static_cast<Triangulation<2,3>*>(triangulation);
        tria->refine_global(n);
      }
    else
      {
        Triangulation<3,3> *tria =
          static_cast<Triangulation<3,3>*>(triangulation);
        tria->refine_global(n);
      }
  }



  void TriangulationWrapper::execute_coarsening_and_refinement()
  {
    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2,2> *tria =
          static_cast<Triangulation<2,2>*>(triangulation);
        tria->execute_coarsening_and_refinement();
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        Triangulation<2,3> *tria =
          static_cast<Triangulation<2,3>*>(triangulation);
        tria->execute_coarsening_and_refinement();
      }
    else
      {
        Triangulation<3,3> *tria =
          static_cast<Triangulation<3,3>*>(triangulation);
        tria->execute_coarsening_and_refinement();
      }
  }



  void TriangulationWrapper::write(const std::string &filename,
                                   const std::string  format) const
  {
    if ((dim == 2) && (spacedim == 2))
      internal::write<2,2>(filename, format, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::write<2,3>(filename, format, triangulation);
    else
      internal::write<3,3>(filename, format, triangulation);
  }



  void TriangulationWrapper::save(const std::string &filename) const
  {
    std::ofstream ofs(filename);
    boost::archive::text_oarchive oa(ofs);

    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2,2> *tria =
          static_cast<Triangulation<2,2>*>(triangulation);

        oa << *tria;
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        {
          Triangulation<2,3> *tria =
            static_cast<Triangulation<2,3>*>(triangulation);

          oa << *tria;
        }
      }
    else
      {
        Triangulation<3,3> *tria =
          static_cast<Triangulation<3,3>*>(triangulation);

        oa << *tria;
      }
  }



  void TriangulationWrapper::load(const std::string &filename)
  {
    std::ifstream ifs(filename);
    boost::archive::text_iarchive ia(ifs);

    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2,2> *tria =
          static_cast<Triangulation<2,2>*>(triangulation);

        ia >> *tria;
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        Triangulation<2,3> *tria =
          static_cast<Triangulation<2,3>*>(triangulation);

        ia >> *tria;
      }
    else
      {
        Triangulation<3> *tria =
          static_cast<Triangulation<3,3>*>(triangulation);

        ia >> *tria;
      }
  }



  boost::python::list TriangulationWrapper::active_cells()
  {
    if ((dim == 2) && (spacedim == 2))
      return internal::active_cells<2,2>(*this);
    else if ((dim == 2) && (spacedim == 3))
      return internal::active_cells<2,3>(*this);
    else
      return internal::active_cells<3,3>(*this);
  }


  void TriangulationWrapper::setup(const std::string &dimension,
                                   const std::string &spacedimension)
  {
    if ((dimension.compare("2D") == 0) || (dimension.compare("2d") == 0))
      {
        dim = 2;

        if ((spacedimension.compare("2D") == 0) ||
            (spacedimension.compare("2d") == 0))
          {
            spacedim = 2;
            triangulation = new Triangulation<2,2>();
          }
        else if ((spacedimension.compare("3D") == 0) ||
                 (spacedimension.compare("3d") == 0))
          {
            spacedim = 3;
            triangulation = new Triangulation<2,3>();
          }
        else
          AssertThrow(false, ExcMessage("Spacedimension needs to be 2D or 3D."));
      }
    else if ((dimension.compare("3D") == 0) || (dimension.compare("3d") == 0))
      {
        if ((spacedimension.compare("3D") != 0) &&
            (spacedimension.compare("3d") != 0))
          AssertThrow(false, ExcMessage("Spacedimension needs to be 3D."));
        dim = 3;
        spacedim = 3;
        triangulation = new Triangulation<3,3>();
      }
    else
      AssertThrow(false, ExcMessage("Dimension needs to be 2D or 3D."));
  }
}

DEAL_II_NAMESPACE_CLOSE
