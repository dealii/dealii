// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/types.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <cell_accessor_wrapper.h>
#include <triangulation_wrapper.h>

#include <fstream>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  namespace internal
  {
    template <int dim, int spacedim>
    void
    create_triangulation(const boost::python::list &vertices_list,
                         const boost::python::list &cells_vertices,
                         void                      *triangulation)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);
      tria->clear();

      const size_t n_vertices = boost::python::len(vertices_list);
      std::vector<Point<spacedim>> vertices(n_vertices);
      for (size_t i = 0; i < n_vertices; ++i)
        {
          boost::python::list vertex =
            boost::python::extract<boost::python::list>(vertices_list[i]);
          for (int d = 0; d < spacedim; ++d)
            vertices[i][d] = boost::python::extract<double>(vertex[d]);
        }

      const size_t               n_cells = boost::python::len(cells_vertices);
      std::vector<CellData<dim>> cell_data(n_cells);
      for (size_t i = 0; i < n_cells; ++i)
        {
          boost::python::list vertex_indices =
            boost::python::extract<boost::python::list>(cells_vertices[i]);

          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
               ++v)
            cell_data[i].vertices[v] =
              boost::python::extract<unsigned int>(vertex_indices[v]);
        }

      tria->create_triangulation(vertices, cell_data, SubCellData());
    }


    template <int dim, int spacedim>
    void
    generate_hyper_cube(const double left,
                        const double right,
                        const bool   colorize,
                        void        *triangulation)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);
      tria->clear();
      GridGenerator::hyper_cube(*tria, left, right, colorize);
    }



    template <int dim>
    void
    generate_simplex(std::vector<PointWrapper> &wrapped_points,
                     void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      std::vector<Point<dim>> points(dim + 1);
      for (int i = 0; i < dim + 1; ++i)
        points[i] =
          *(static_cast<Point<dim> *>((wrapped_points[i]).get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::simplex(*tria, points);
    }



    template <int dim, int spacedim>
    void
    generate_subdivided_hyper_cube(const unsigned int repetitions,
                                   const double       left,
                                   const double       right,
                                   void              *triangulation)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);
      tria->clear();
      GridGenerator::subdivided_hyper_cube(*tria, repetitions, left, right);
    }



    template <int dim, int spacedim>
    void
    generate_hyper_rectangle(PointWrapper &p1,
                             PointWrapper &p2,
                             const bool    colorize,
                             void         *triangulation)
    {
      AssertThrow(
        p1.get_dim() == dim,
        ExcMessage(
          "Dimension of p1 is not the same as the dimension of the Triangulation."));
      AssertThrow(
        p2.get_dim() == dim,
        ExcMessage(
          "Dimension of p2 is not the same as the dimension of the Triangulation."));
      // Cast the PointWrapper object to Point<dim>
      Point<dim> point_1 = *(static_cast<Point<dim> *>(p1.get_point()));
      Point<dim> point_2 = *(static_cast<Point<dim> *>(p2.get_point()));

      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);
      tria->clear();
      GridGenerator::hyper_rectangle(*tria, point_1, point_2, colorize);
    }


    template <int dim>
    void
    generate_hyper_cube_with_cylindrical_hole(const double       inner_radius,
                                              const double       outer_radius,
                                              const double       L,
                                              const unsigned int repetitions,
                                              const bool         colorize,
                                              void              *triangulation)
    {
      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::hyper_cube_with_cylindrical_hole(
        *tria, inner_radius, outer_radius, L, repetitions, colorize);
    }



    template <int dim, int spacedim>
    void
    generate_subdivided_hyper_rectangle(
      const std::vector<unsigned int> &repetitions,
      PointWrapper                    &p1,
      PointWrapper                    &p2,
      const bool                       colorize,
      void                            *triangulation)
    {
      AssertThrow(
        p1.get_dim() == dim,
        ExcMessage(
          "Dimension of p1 is not the same as the dimension of the Triangulation."));
      AssertThrow(
        p2.get_dim() == dim,
        ExcMessage(
          "Dimension of p2 is not the same as the dimension of the Triangulation."));
      // Cast the PointWrapper object to Point<dim>
      Point<dim> point_1 = *(static_cast<Point<dim> *>(p1.get_point()));
      Point<dim> point_2 = *(static_cast<Point<dim> *>(p2.get_point()));

      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);
      tria->clear();
      GridGenerator::subdivided_hyper_rectangle(
        *tria, repetitions, point_1, point_2, colorize);
    }



    template <int dim>
    void
    generate_subdivided_steps_hyper_rectangle(
      const std::vector<std::vector<double>> &step_sizes,
      PointWrapper                           &p1,
      PointWrapper                           &p2,
      const bool                              colorize,
      void                                   *triangulation)
    {
      AssertThrow(
        p1.get_dim() == dim,
        ExcMessage(
          "Dimension of p1 is not the same as the dimension of the Triangulation."));
      AssertThrow(
        p2.get_dim() == dim,
        ExcMessage(
          "Dimension of p2 is not the same as the dimension of the Triangulation."));
      // Cast the PointWrapper object to Point<dim>
      Point<dim> point_1 = *(static_cast<Point<dim> *>(p1.get_point()));
      Point<dim> point_2 = *(static_cast<Point<dim> *>(p2.get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::subdivided_hyper_rectangle(
        *tria, step_sizes, point_1, point_2, colorize);
    }



    template <int dim>
    void
    generate_subdivided_material_hyper_rectangle(
      const std::vector<std::vector<double>> &spacing,
      PointWrapper                           &p,
      const Table<dim, types::material_id>   &material_ids,
      const bool                              colorize,
      void                                   *triangulation)
    {
      AssertThrow(
        p.get_dim() == dim,
        ExcMessage(
          "Dimension of p is not the same as the dimension of the Triangulation."));
      // Cast the PointWrapper object to Point<dim>
      Point<dim>          point = *(static_cast<Point<dim> *>(p.get_point()));
      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::subdivided_hyper_rectangle(
        *tria, spacing, point, material_ids, colorize);
    }



    template <int dim, int spacedim>
    void
    generate_cheese(const std::vector<unsigned int> &holes, void *triangulation)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);
      tria->clear();
      GridGenerator::cheese(*tria, holes);
    }



    template <int dim>
    void
    generate_plate_with_a_hole(const double        inner_radius,
                               const double        outer_radius,
                               const double        pad_bottom,
                               const double        pad_top,
                               const double        pad_left,
                               const double        pad_right,
                               const PointWrapper &center,
                               const int           polar_manifold_id,
                               const int           tfi_manifold_id,
                               const double        L,
                               const unsigned int  n_slices,
                               const bool          colorize,
                               void               *triangulation)
    {
      // Cast the PointWrapper object to Point<dim>
      Point<dim> point =
        center.get_point() ?
          *(static_cast<const Point<dim> *>(center.get_point())) :
          Point<dim>();

      Triangulation<dim, dim> *tria =
        static_cast<Triangulation<dim, dim> *>(triangulation);
      tria->clear();
      GridGenerator::plate_with_a_hole(*tria,
                                       inner_radius,
                                       outer_radius,
                                       pad_bottom,
                                       pad_top,
                                       pad_left,
                                       pad_right,
                                       point,
                                       polar_manifold_id,
                                       tfi_manifold_id,
                                       L,
                                       n_slices,
                                       colorize);
    }



    template <int dim>
    void
    generate_channel_with_cylinder(const double       shell_region_width,
                                   const unsigned int n_shells,
                                   const double       skewness,
                                   const bool         colorize,
                                   void              *triangulation)
    {
      Triangulation<dim, dim> *tria =
        static_cast<Triangulation<dim, dim> *>(triangulation);
      tria->clear();
      GridGenerator::channel_with_cylinder(
        *tria, shell_region_width, n_shells, skewness, colorize);
    }



    template <int dim>
    void
    generate_general_cell(std::vector<PointWrapper> &wrapped_points,
                          const bool                 colorize,
                          void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      const unsigned int      size = wrapped_points.size();
      std::vector<Point<dim>> points(size);
      for (unsigned int i = 0; i < size; ++i)
        points[i] =
          *(static_cast<Point<dim> *>((wrapped_points[i]).get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::general_cell(*tria, points, colorize);
    }



    template <int dim>
    void
    generate_parallelogram(std::vector<PointWrapper> &wrapped_points,
                           const bool                 colorize,
                           void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      Point<dim> points[dim];
      for (unsigned int i = 0; i < dim; ++i)
        points[i] =
          *(static_cast<Point<dim> *>((wrapped_points[i]).get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::parallelogram(*tria, points, colorize);
    }



    template <int dim>
    void
    generate_parallelepiped(std::vector<PointWrapper> &wrapped_points,
                            const bool                 colorize,
                            void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      Point<dim> points[dim];
      for (unsigned int i = 0; i < dim; ++i)
        points[i] =
          *(static_cast<Point<dim> *>((wrapped_points[i]).get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::parallelepiped(*tria, points, colorize);
    }



    template <int dim>
    void
    generate_fixed_subdivided_parallelepiped(
      unsigned int               n_subdivisions,
      std::vector<PointWrapper> &wrapped_points,
      const bool                 colorize,
      void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      Point<dim> points[dim];
      for (unsigned int i = 0; i < dim; ++i)
        points[i] =
          *(static_cast<Point<dim> *>((wrapped_points[i]).get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::subdivided_parallelepiped(*tria,
                                               n_subdivisions,
                                               points,
                                               colorize);
    }



    template <int dim>
    void
    generate_varying_subdivided_parallelepiped(
      std::vector<unsigned int> &n_subdivisions,
      std::vector<PointWrapper> &wrapped_points,
      const bool                 colorize,
      void                      *triangulation)
    {
      // Cast the PointWrapper objects to Point<dim>
      Point<dim>   points[dim];
      unsigned int subdivisions[dim];
      for (unsigned int i = 0; i < dim; ++i)
        {
          points[i] =
            *(static_cast<Point<dim> *>((wrapped_points[i]).get_point()));
          subdivisions[i] = n_subdivisions[i];
        }

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::subdivided_parallelepiped(*tria,
                                               subdivisions,
                                               points,
                                               colorize);
    }



    template <int dim>
    void
    generate_enclosed_hyper_cube(const double left,
                                 const double right,
                                 const double thickness,
                                 const bool   colorize,
                                 void        *triangulation)
    {
      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::enclosed_hyper_cube(
        *tria, left, right, thickness, colorize);
    }



    template <int dim>
    void
    generate_hyper_ball(PointWrapper &center,
                        const double  radius,
                        void         *triangulation)
    {
      // Cast the PointWrapper object to Point<dim>
      Point<dim> center_point =
        *(static_cast<Point<dim> *>(center.get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::hyper_ball(*tria, center_point, radius);
    }



    template <int dim>
    void
    generate_hyper_ball_balanced(const PointWrapper &center,
                                 const double        radius,
                                 void               *triangulation)
    {
      // Cast the PointWrapper object to Point<dim>
      Point<dim> center_point =
        center.get_point() ?
          *(static_cast<const Point<dim> *>(center.get_point())) :
          Point<dim>();

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::hyper_ball_balanced(*tria, center_point, radius);
    }



    template <int dim, int spacedim>
    void
    generate_hyper_sphere(PointWrapper &center,
                          const double  radius,
                          void         *triangulation)
    {
      // Cast the PointWrapper object to Point<dim>
      Point<spacedim> center_point =
        *(static_cast<Point<spacedim> *>(center.get_point()));

      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);
      tria->clear();
      GridGenerator::hyper_sphere(*tria, center_point, radius);
    }


    template <int dim>
    void
    generate_hyper_shell(PointWrapper  &center,
                         const double   inner_radius,
                         const double   outer_radius,
                         const unsigned n_cells,
                         bool           colorize,
                         void          *triangulation)
    {
      // Cast the PointWrapper object to Point<dim>
      Point<dim> center_point =
        *(static_cast<Point<dim> *>(center.get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::hyper_shell(
        *tria, center_point, inner_radius, outer_radius, n_cells, colorize);
    }


    template <int dim>
    void
    generate_quarter_hyper_ball(PointWrapper &center,
                                const double  radius,
                                void         *triangulation)
    {
      // Cast the PointWrapper object to Point<dim>
      Point<dim> center_point =
        *(static_cast<Point<dim> *>(center.get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::quarter_hyper_ball(*tria, center_point, radius);
    }



    template <int dim>
    void
    generate_half_hyper_ball(PointWrapper &center,
                             const double  radius,
                             void         *triangulation)
    {
      // Cast the PointWrapper object to Point<dim>
      Point<dim> center_point =
        *(static_cast<Point<dim> *>(center.get_point()));

      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::half_hyper_ball(*tria, center_point, radius);
    }



    template <int dim>
    void
    generate_cylinder(const double radius,
                      const double half_length,
                      void        *triangulation)
    {
      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::cylinder(*tria, radius, half_length);
    }



    template <int dim>
    void
    generate_subdivided_cylinder(const unsigned int x_subdivisions,
                                 const double       radius,
                                 const double       half_length,
                                 void              *triangulation)
    {
      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::subdivided_cylinder(*tria,
                                         x_subdivisions,
                                         radius,
                                         half_length);
    }



    template <int dim>
    void
    generate_truncated_cone(const double radius_0,
                            const double radius_1,
                            const double half_length,
                            void        *triangulation)
    {
      Triangulation<dim> *tria =
        static_cast<Triangulation<dim> *>(triangulation);
      tria->clear();
      GridGenerator::truncated_cone(*tria, radius_0, radius_1, half_length);
    }



    template <int dim, int spacedim>
    void
    scale(const double scaling_factor, void *triangulation)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);
      GridTools::scale(scaling_factor, *tria);
    }



    template <int dim, int spacedim>
    void
    shift(boost::python::list &shift_list, void *triangulation)
    {
      // Extract the shift vector from the python list
      Tensor<1, spacedim> shift_vector;
      for (int i = 0; i < spacedim; ++i)
        shift_vector[i] = boost::python::extract<double>(shift_list[i]);

      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);
      GridTools::shift(shift_vector, *tria);
    }



    template <int dim, int spacedim>
    void
    merge_triangulations(boost::python::list &triangulations,
                         void                *triangulation,
                         const double         duplicated_vertex_tolerance,
                         const bool           copy_manifold_ids)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);
      tria->clear();

      std::vector<const Triangulation<dim, spacedim> *> tria_ptrs;
      size_t n_triangulations = boost::python::len(triangulations);
      for (size_t i = 0; i < n_triangulations; ++i)
        {
          TriangulationWrapper &tria_wrapper =
            boost::python::extract<TriangulationWrapper &>(triangulations[i]);
          tria_ptrs.push_back(static_cast<const Triangulation<dim, spacedim> *>(
            tria_wrapper.get_triangulation()));
        }

      GridGenerator::merge_triangulations(tria_ptrs,
                                          *tria,
                                          duplicated_vertex_tolerance,
                                          copy_manifold_ids);
    }



    template <int dim, int spacedim_1, int spacedim_2>
    void
    flatten_triangulation(void *triangulation, TriangulationWrapper &tria_out)
    {
      Triangulation<dim, spacedim_1> *tria =
        static_cast<Triangulation<dim, spacedim_1> *>(triangulation);
      Triangulation<dim, spacedim_2> *tria_2 =
        static_cast<Triangulation<dim, spacedim_2> *>(
          tria_out.get_triangulation());
      GridGenerator::flatten_triangulation(*tria, *tria_2);
    }



    template <int dim, int spacedim>
    void
    convert_hypercube_to_simplex_mesh(void                 *triangulation,
                                      TriangulationWrapper &tria_out)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);
      Triangulation<dim, spacedim> *tria_2 =
        static_cast<Triangulation<dim, spacedim> *>(
          tria_out.get_triangulation());
      GridGenerator::convert_hypercube_to_simplex_mesh(*tria, *tria_2);
    }



    template <int dim, int spacedim>
    void
    replicate_triangulation(TriangulationWrapper      &tria_in,
                            const boost::python::list &extents,
                            void                      *tria_out)
    {
      Triangulation<dim, spacedim> *replicated_triangulation =
        static_cast<Triangulation<dim, spacedim> *>(tria_out);
      const Triangulation<dim, spacedim> *triangulation =
        static_cast<Triangulation<dim, spacedim> *>(
          tria_in.get_triangulation());
      replicated_triangulation->clear();

      std::vector<unsigned int> extents_vector(dim);
      for (int d = 0; d < dim; ++d)
        extents_vector[d] = boost::python::extract<unsigned>(extents[d]);

      GridGenerator::replicate_triangulation(*triangulation,
                                             extents_vector,
                                             *replicated_triangulation);
    }



    template <int dim, int spacedim>
    void
    distort_random(const double factor,
                   const bool   keep_boundary,
                   void        *triangulation)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);

      GridTools::distort_random(factor, *tria, keep_boundary);
    }



    template <int dim, int spacedim>
    void
    transform(boost::python::object &transformation, void *triangulation)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);

      auto func = [&transformation](const Point<spacedim> &p_in) {
        boost::python::list p_list_in;

        for (int d = 0; d < spacedim; ++d)
          p_list_in.append(p_in[d]);

        boost::python::list p_list_out =
          boost::python::extract<boost::python::list>(
            transformation(p_list_in));

        Point<spacedim> p_out;

        for (int d = 0; d < spacedim; ++d)
          p_out[d] = boost::python::extract<double>(p_list_out[d]);

        return p_out;
      };

      GridTools::transform(func, *tria);
    }



    template <int dim, int spacedim>
    std::pair<int, int>
    find_active_cell_around_point(PointWrapper    &p,
                                  MappingQWrapper &mapping_wrapper,
                                  void            *triangulation)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);

      Point<spacedim> point = *(static_cast<Point<spacedim> *>(p.get_point()));

      if (mapping_wrapper.get_mapping() != nullptr)
        {
          const MappingQ<dim, spacedim> *mapping =
            static_cast<const MappingQ<dim, spacedim> *>(
              mapping_wrapper.get_mapping());

          auto cell_pair =
            GridTools::find_active_cell_around_point(*mapping, *tria, point);
          return std::make_pair(cell_pair.first->level(),
                                cell_pair.first->index());
        }
      else
        {
          auto cell = GridTools::find_active_cell_around_point(*tria, point);
          return std::make_pair(cell->level(), cell->index());
        }
    }


    template <int dim, int spacedim>
    boost::python::list
    compute_aspect_ratio_of_cells(
      const MappingQWrapper      &mapping_wrapper,
      const QuadratureWrapper    &quadrature_wrapper,
      const TriangulationWrapper &triangulation_wrapper)
    {
      const Triangulation<dim, spacedim> *tria =
        static_cast<const Triangulation<dim, spacedim> *>(
          triangulation_wrapper.get_triangulation());

      const Quadrature<dim> *quad = static_cast<const Quadrature<dim> *>(
        quadrature_wrapper.get_quadrature());

      const MappingQ<dim, spacedim> *mapping =
        static_cast<const MappingQ<dim, spacedim> *>(
          mapping_wrapper.get_mapping());

      auto aspect_ratios =
        GridTools::compute_aspect_ratio_of_cells(*mapping, *tria, *quad);

      boost::python::list ratios;
      for (size_t i = 0; i < aspect_ratios.size(); ++i)
        ratios.append(aspect_ratios[i]);

      return ratios;
    }


    template <int dim, int spacedim>
    boost::python::list
    find_cells_adjacent_to_vertex(const unsigned int    vertex_index,
                                  TriangulationWrapper &triangulation_wrapper)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(
          triangulation_wrapper.get_triangulation());

      std::vector<typename Triangulation<dim, spacedim>::active_cell_iterator>
        adjacent_cells =
          GridTools::find_cells_adjacent_to_vertex(*tria, vertex_index);

      boost::python::list cells;
      for (auto &cell : adjacent_cells)
        cells.append(CellAccessorWrapper(triangulation_wrapper,
                                         cell->level(),
                                         cell->index()));

      return cells;
    }



    template <int dim, int spacedim>
    boost::python::list
    active_cells(TriangulationWrapper &triangulation_wrapper)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(
          triangulation_wrapper.get_triangulation());
      boost::python::list cells;
      for (auto &cell : tria->active_cell_iterators())
        cells.append(CellAccessorWrapper(triangulation_wrapper,
                                         cell->level(),
                                         cell->index()));

      return cells;
    }



    template <int dim, int spacedim>
    boost::python::list
    cells(TriangulationWrapper &triangulation_wrapper)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(
          triangulation_wrapper.get_triangulation());
      boost::python::list cells_list;
      for (auto &cell : tria->cell_iterators())
        cells_list.append(CellAccessorWrapper(triangulation_wrapper,
                                              cell->level(),
                                              cell->index()));

      return cells_list;
    }



    template <int dim, int spacedim>
    double
    maximal_cell_diameter(const void *triangulation)
    {
      const Triangulation<dim, spacedim> *tria =
        static_cast<const Triangulation<dim, spacedim> *>(triangulation);
      return GridTools::maximal_cell_diameter(*tria);
    }



    template <int dim, int spacedim>
    double
    minimal_cell_diameter(const void *triangulation)
    {
      const Triangulation<dim, spacedim> *tria =
        static_cast<const Triangulation<dim, spacedim> *>(triangulation);
      return GridTools::minimal_cell_diameter(*tria);
    }



    template <int dim, int spacedim>
    void
    write(const std::string &filename,
          const std::string &format,
          const void        *triangulation)
    {
      const Triangulation<dim, spacedim> *tria =
        static_cast<const Triangulation<dim, spacedim> *>(triangulation);

      GridOut mesh_writer;

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
        {
          output_format = GridOut::OutputFormat::vtu;
          GridOutFlags::Vtu flags(/* serialize_triangulation = */ true);
          mesh_writer.set_flags(flags);
        }
      else
        output_format = GridOut::OutputFormat::none;

      std::ofstream ofs(filename);
      mesh_writer.write(*tria, ofs, output_format);
      ofs.close();
    }



    template <int dim, int spacedim>
    void
    read(const std::string &filename,
         const std::string &format,
         void              *triangulation)
    {
      Triangulation<dim, spacedim> *tria =
        static_cast<Triangulation<dim, spacedim> *>(triangulation);

      tria->clear();

      typename GridIn<dim, spacedim>::Format input_format =
        GridIn<dim, spacedim>::Format::Default;
      if (format.compare("msh") == 0)
        input_format = GridIn<dim, spacedim>::Format::msh;
      else if (format.compare("vtk") == 0)
        input_format = GridIn<dim, spacedim>::Format::vtk;
      else
        Assert(false,
               ExcMessage("Cannot read triangulation of the given format."));

      GridIn<dim, spacedim> mesh_reader;
      mesh_reader.attach_triangulation(*tria);
      std::ifstream ifs(filename);
      AssertThrow(ifs, ExcIO());
      mesh_reader.read(ifs, input_format);
      ifs.close();
    }
  } // namespace internal



  TriangulationWrapper::TriangulationWrapper(
    const std::string &dimension,
    const int          mesh_smoothing,
    const bool         check_for_distorted_cells)
  {
    if ((dimension.compare("2D") == 0) || (dimension.compare("2d") == 0))
      setup("2D", "2D", mesh_smoothing, check_for_distorted_cells);
    else if ((dimension.compare("3D") == 0) || (dimension.compare("3d") == 0))
      setup("3D", "3D", mesh_smoothing, check_for_distorted_cells);
    else
      AssertThrow(false, ExcMessage("Dimension needs to be 2D or 3D"));
  }



  TriangulationWrapper::TriangulationWrapper(
    const std::string &dimension,
    const std::string &spacedimension,
    const int          mesh_smoothing,
    const bool         check_for_distorted_cells)
  {
    setup(dimension, spacedimension, mesh_smoothing, check_for_distorted_cells);
  }



  TriangulationWrapper::~TriangulationWrapper()
  {
    if (triangulation != nullptr)
      {
        if (dim == 2)
          {
            if (spacedim == 2)
              {
                // We cannot call delete on a void pointer so cast the void
                // pointer back first.
                Triangulation<2, 2> *tmp =
                  static_cast<Triangulation<2, 2> *>(triangulation);
                delete tmp;
              }
            else
              {
                Triangulation<2, 3> *tmp =
                  static_cast<Triangulation<2, 3> *>(triangulation);
                delete tmp;
              }
          }
        else
          {
            Triangulation<3, 3> *tmp =
              static_cast<Triangulation<3, 3> *>(triangulation);
            delete tmp;
          }
        triangulation = nullptr;
      }
    dim = -1;
  }



  unsigned int
  TriangulationWrapper::n_active_cells() const
  {
    if ((dim == 2) && (spacedim == 2))
      return (*static_cast<Triangulation<2, 2> *>(triangulation))
        .n_active_cells();
    else if ((dim == 2) && (spacedim == 3))
      return (*static_cast<Triangulation<2, 3> *>(triangulation))
        .n_active_cells();
    else
      return (*static_cast<Triangulation<3, 3> *>(triangulation))
        .n_active_cells();
  }



  unsigned int
  TriangulationWrapper::n_cells() const
  {
    if ((dim == 2) && (spacedim == 2))
      return (*static_cast<Triangulation<2, 2> *>(triangulation)).n_cells();
    else if ((dim == 2) && (spacedim == 3))
      return (*static_cast<Triangulation<2, 3> *>(triangulation)).n_cells();
    else
      return (*static_cast<Triangulation<3, 3> *>(triangulation)).n_cells();
  }



  void
  TriangulationWrapper::create_triangulation(
    const boost::python::list &vertices,
    const boost::python::list &cells_vertices)
  {
    if ((dim == 2) && (spacedim == 2))
      internal::create_triangulation<2, 2>(vertices,
                                           cells_vertices,
                                           triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::create_triangulation<2, 3>(vertices,
                                           cells_vertices,
                                           triangulation);
    else
      internal::create_triangulation<3, 3>(vertices,
                                           cells_vertices,
                                           triangulation);
  }



  void
  TriangulationWrapper::generate_hyper_cube(const double left,
                                            const double right,
                                            const bool   colorize)
  {
    if ((dim == 2) && (spacedim == 2))
      internal::generate_hyper_cube<2, 2>(left, right, colorize, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::generate_hyper_cube<2, 3>(left, right, colorize, triangulation);
    else
      internal::generate_hyper_cube<3, 3>(left, right, colorize, triangulation);
  }



  void
  TriangulationWrapper::generate_simplex(boost::python::list &vertices)
  {
    AssertThrow(boost::python::len(vertices) == dim + 1,
                ExcMessage("The number of vertices should be equal to dim+1."));
    AssertThrow(
      dim == spacedim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    // Extract the PointWrapper object from the python list
    std::vector<PointWrapper> wrapped_points(dim + 1);
    for (int i = 0; i < dim + 1; ++i)
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



  void
  TriangulationWrapper::generate_subdivided_hyper_cube(
    const unsigned int repetitions,
    const double       left,
    const double       right)
  {
    if ((dim == 2) && (spacedim == 2))
      internal::generate_subdivided_hyper_cube<2, 2>(repetitions,
                                                     left,
                                                     right,
                                                     triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::generate_subdivided_hyper_cube<2, 3>(repetitions,
                                                     left,
                                                     right,
                                                     triangulation);
    else
      internal::generate_subdivided_hyper_cube<3, 3>(repetitions,
                                                     left,
                                                     right,
                                                     triangulation);
  }



  void
  TriangulationWrapper::generate_hyper_rectangle(PointWrapper &p1,
                                                 PointWrapper &p2,
                                                 const bool    colorize)
  {
    if ((dim == 2) && (spacedim == 2))
      internal::generate_hyper_rectangle<2, 2>(p1, p2, colorize, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::generate_hyper_rectangle<2, 3>(p1, p2, colorize, triangulation);
    else
      internal::generate_hyper_rectangle<3, 3>(p1, p2, colorize, triangulation);
  }



  void
  TriangulationWrapper::generate_subdivided_hyper_rectangle(
    boost::python::list &repetition_list,
    PointWrapper        &p1,
    PointWrapper        &p2,
    const bool           colorize)
  {
    AssertThrow(
      boost::python::len(repetition_list) == dim,
      ExcMessage(
        "The list of repetitions must have the same length as the number of dimension."));

    // Extract the repetitions from the python list
    std::vector<unsigned int> repetitions(dim);
    for (int i = 0; i < dim; ++i)
      repetitions[i] = boost::python::extract<unsigned int>(repetition_list[i]);

    if ((dim == 2) && (spacedim == 2))
      internal::generate_subdivided_hyper_rectangle<2, 2>(
        repetitions, p1, p2, colorize, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::generate_subdivided_hyper_rectangle<2, 3>(
        repetitions, p1, p2, colorize, triangulation);
    else
      internal::generate_subdivided_hyper_rectangle<3, 3>(
        repetitions, p1, p2, colorize, triangulation);
  }



  void
  TriangulationWrapper::generate_subdivided_steps_hyper_rectangle(
    boost::python::list &step_sizes_list,
    PointWrapper        &p1,
    PointWrapper        &p2,
    const bool           colorize)
  {
    AssertThrow(
      spacedim == dim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    AssertThrow(
      boost::python::len(step_sizes_list) == dim,
      ExcMessage(
        "The list of step_sizes must have the same length as the number of dimension."));

    // Extract the step sizes from the python list
    std::vector<std::vector<double>> step_sizes(dim);
    for (int i = 0; i < dim; ++i)
      {
        step_sizes[i].resize(boost::python::len(step_sizes_list[i]));
        for (unsigned int j = 0; j < step_sizes[i].size(); ++j)
          step_sizes[i][j] =
            boost::python::extract<double>(step_sizes_list[i][j]);
      }

    if (dim == 2)
      internal::generate_subdivided_steps_hyper_rectangle<2>(
        step_sizes, p1, p2, colorize, triangulation);
    else
      internal::generate_subdivided_steps_hyper_rectangle<3>(
        step_sizes, p1, p2, colorize, triangulation);
  }



  void
  TriangulationWrapper::generate_subdivided_material_hyper_rectangle(
    boost::python::list &spacing_list,
    PointWrapper        &p,
    boost::python::list &material_id_list,
    const bool           colorize)
  {
    AssertThrow(
      spacedim == dim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    AssertThrow(
      boost::python::len(spacing_list) == dim,
      ExcMessage(
        "The list of spacing must have the same length as the number of dimension."));

    // Extract the spacing and the material ID from the python list
    std::vector<std::vector<double>> spacing(dim);
    for (int i = 0; i < dim; ++i)
      {
        spacing[i].resize(boost::python::len(spacing_list[i]));
        for (unsigned int j = 0; j < spacing[i].size(); ++j)
          spacing[i][j] = boost::python::extract<double>(spacing_list[i][j]);
      }
    if (dim == 2)
      {
        const unsigned int index_0 = boost::python::len(material_id_list);
        const unsigned int index_1 = boost::python::len(material_id_list[0]);
        Table<2, types::material_id> material_ids(index_0, index_1);
        for (unsigned int i = 0; i < index_0; ++i)
          for (unsigned int j = 0; j < index_1; ++j)
            // We cannot use extract<types::material_id> because boost will
            // throw an exception if we try to extract -1
            material_ids[i][j] =
              boost::python::extract<int>(material_id_list[i][j]);

        internal::generate_subdivided_material_hyper_rectangle<2>(
          spacing, p, material_ids, colorize, triangulation);
      }
    else
      {
        const unsigned int index_0 = boost::python::len(material_id_list);
        const unsigned int index_1 = boost::python::len(material_id_list[0]);
        const unsigned int index_2 = boost::python::len(material_id_list[0][0]);
        Table<3, types::material_id> material_ids(index_0, index_1, index_2);
        for (unsigned int i = 0; i < index_0; ++i)
          for (unsigned int j = 0; j < index_1; ++j)
            for (unsigned int k = 0; k < index_2; ++k)
              material_ids[i][j][k] =
                boost::python::extract<int>(material_id_list[i][j][k]);
        internal::generate_subdivided_material_hyper_rectangle<3>(
          spacing, p, material_ids, colorize, triangulation);
      }
  }



  void
  TriangulationWrapper::generate_cheese(boost::python::list &holes_list)
  {
    const unsigned int        size = boost::python::len(holes_list);
    std::vector<unsigned int> holes(size);
    for (unsigned int i = 0; i < size; ++i)
      holes[i] = boost::python::extract<unsigned int>(holes_list[i]);

    if ((dim == 2) && (spacedim == 2))
      internal::generate_cheese<2, 2>(holes, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::generate_cheese<2, 3>(holes, triangulation);
    else
      internal::generate_cheese<3, 3>(holes, triangulation);
  }



  void
  TriangulationWrapper::generate_plate_with_a_hole(const double inner_radius,
                                                   const double outer_radius,
                                                   const double pad_bottom,
                                                   const double pad_top,
                                                   const double pad_left,
                                                   const double pad_right,
                                                   const PointWrapper &center,
                                                   const int polar_manifold_id,
                                                   const int tfi_manifold_id,
                                                   const double       L,
                                                   const unsigned int n_slices,
                                                   const bool         colorize)
  {
    if (dim == 2)
      internal::generate_plate_with_a_hole<2>(inner_radius,
                                              outer_radius,
                                              pad_bottom,
                                              pad_top,
                                              pad_left,
                                              pad_right,
                                              center,
                                              polar_manifold_id,
                                              tfi_manifold_id,
                                              L,
                                              n_slices,
                                              colorize,
                                              triangulation);
    else
      internal::generate_plate_with_a_hole<3>(inner_radius,
                                              outer_radius,
                                              pad_bottom,
                                              pad_top,
                                              pad_left,
                                              pad_right,
                                              center,
                                              polar_manifold_id,
                                              tfi_manifold_id,
                                              L,
                                              n_slices,
                                              colorize,
                                              triangulation);
  }



  void
  TriangulationWrapper::generate_channel_with_cylinder(
    const double       shell_region_width,
    const unsigned int n_shells,
    const double       skewness,
    const bool         colorize)
  {
    if (dim == 2)
      internal::generate_channel_with_cylinder<2>(
        shell_region_width, n_shells, skewness, colorize, triangulation);
    else
      internal::generate_channel_with_cylinder<3>(
        shell_region_width, n_shells, skewness, colorize, triangulation);
  }



  void
  TriangulationWrapper::generate_general_cell(boost::python::list &vertices,
                                              const bool           colorize)
  {
    AssertThrow(
      spacedim == dim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    // Extract the PointWrapper object from the python list
    const int size = boost::python::len(vertices);
    AssertThrow(size > 0, ExcMessage("The vertices list is empty."));
    std::vector<PointWrapper> wrapped_points(size);
    for (int i = 0; i < size; ++i)
      wrapped_points[i] = boost::python::extract<PointWrapper>(vertices[i]);
    if (dim == 2)
      internal::generate_general_cell<2>(wrapped_points,
                                         colorize,
                                         triangulation);
    else
      internal::generate_general_cell<3>(wrapped_points,
                                         colorize,
                                         triangulation);
  }



  void
  TriangulationWrapper::generate_parallelogram(boost::python::list &corners,
                                               const bool           colorize)
  {
    AssertThrow(
      spacedim == dim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    // Extract the PointWrapper object from the python list
    AssertThrow(
      boost::python::len(corners) == dim,
      ExcMessage(
        "The list of corners must have the same length as the number of dimension."));
    std::vector<PointWrapper> wrapped_points(dim);
    for (int i = 0; i < dim; ++i)
      wrapped_points[i] = boost::python::extract<PointWrapper>(corners[i]);
    if (dim == 2)
      internal::generate_parallelogram<2>(wrapped_points,
                                          colorize,
                                          triangulation);
    else
      internal::generate_parallelogram<3>(wrapped_points,
                                          colorize,
                                          triangulation);
  }



  void
  TriangulationWrapper::generate_parallelepiped(boost::python::list &corners,
                                                const bool           colorize)
  {
    AssertThrow(
      spacedim == dim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    // Extract the PointWrapper object from the python list
    AssertThrow(
      boost::python::len(corners) == dim,
      ExcMessage(
        "The list of corners must have the same length as the number of dimension."));
    std::vector<PointWrapper> wrapped_points(dim);
    for (int i = 0; i < dim; ++i)
      wrapped_points[i] = boost::python::extract<PointWrapper>(corners[i]);
    if (dim == 2)
      internal::generate_parallelepiped<2>(wrapped_points,
                                           colorize,
                                           triangulation);
    else
      internal::generate_parallelepiped<3>(wrapped_points,
                                           colorize,
                                           triangulation);
  }



  void
  TriangulationWrapper::generate_fixed_subdivided_parallelepiped(
    const unsigned int   n_subdivisions,
    boost::python::list &corners,
    const bool           colorize)
  {
    AssertThrow(
      spacedim == dim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    // Extract the PointWrapper object from the python list
    AssertThrow(
      boost::python::len(corners) == dim,
      ExcMessage(
        "The list of corners must have the same length as the number of dimension."));
    std::vector<PointWrapper> wrapped_points(dim);
    for (int i = 0; i < dim; ++i)
      wrapped_points[i] = boost::python::extract<PointWrapper>(corners[i]);
    if ((dim == 2) && (spacedim == 2))
      internal::generate_fixed_subdivided_parallelepiped<2>(n_subdivisions,
                                                            wrapped_points,
                                                            colorize,
                                                            triangulation);
    else
      internal::generate_fixed_subdivided_parallelepiped<3>(n_subdivisions,
                                                            wrapped_points,
                                                            colorize,
                                                            triangulation);
  }



  void
  TriangulationWrapper::generate_varying_subdivided_parallelepiped(
    boost::python::list &n_subdivisions,
    boost::python::list &corners,
    const bool           colorize)
  {
    AssertThrow(
      spacedim == dim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    // Extract the subdivisions from the python list
    AssertThrow(
      boost::python::len(n_subdivisions) == dim,
      ExcMessage(
        "The list of subdivisions must have the same length as the number of dimension."));
    std::vector<unsigned int> subdivisions(dim);
    for (int i = 0; i < dim; ++i)
      subdivisions[i] = boost::python::extract<unsigned int>(n_subdivisions[i]);
    // Extract the PointWrapper object from the python list
    AssertThrow(
      boost::python::len(corners) == dim,
      ExcMessage(
        "The list of corners must have the same length as the number of dimension."));
    std::vector<PointWrapper> wrapped_points(dim);
    for (int i = 0; i < dim; ++i)
      wrapped_points[i] = boost::python::extract<PointWrapper>(corners[i]);
    if (dim == 2)
      internal::generate_varying_subdivided_parallelepiped<2>(subdivisions,
                                                              wrapped_points,
                                                              colorize,
                                                              triangulation);
    else
      internal::generate_varying_subdivided_parallelepiped<3>(subdivisions,
                                                              wrapped_points,
                                                              colorize,
                                                              triangulation);
  }



  void
  TriangulationWrapper::generate_enclosed_hyper_cube(const double left,
                                                     const double right,
                                                     const double thickness,
                                                     const bool   colorize)
  {
    AssertThrow(
      spacedim == dim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_enclosed_hyper_cube<2>(
        left, right, thickness, colorize, triangulation);
    else
      internal::generate_enclosed_hyper_cube<3>(
        left, right, thickness, colorize, triangulation);
  }



  void
  TriangulationWrapper::generate_hyper_ball(PointWrapper &center,
                                            const double  radius)
  {
    AssertThrow(
      dim == spacedim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_hyper_ball<2>(center, radius, triangulation);
    else
      internal::generate_hyper_ball<3>(center, radius, triangulation);
  }



  void
  TriangulationWrapper::generate_hyper_ball_balanced(const PointWrapper &center,
                                                     const double        radius)
  {
    AssertThrow(
      dim == spacedim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_hyper_ball_balanced<2>(center, radius, triangulation);
    else
      internal::generate_hyper_ball_balanced<3>(center, radius, triangulation);
  }



  void
  TriangulationWrapper::generate_hyper_shell(PointWrapper  &center,
                                             const double   inner_radius,
                                             const double   outer_radius,
                                             const unsigned n_cells,
                                             bool           colorize)
  {
    AssertThrow(
      dim == spacedim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_hyper_shell<2>(
        center, inner_radius, outer_radius, n_cells, colorize, triangulation);
    else
      internal::generate_hyper_shell<3>(
        center, inner_radius, outer_radius, n_cells, colorize, triangulation);
  }



  void
  TriangulationWrapper::generate_hyper_sphere(PointWrapper &center,
                                              const double  radius)
  {
    AssertThrow(
      spacedim == dim + 1,
      ExcMessage(
        "This function is only implemented for spacedim equal to dim+1."));
    internal::generate_hyper_sphere<2, 3>(center, radius, triangulation);
  }



  void
  TriangulationWrapper::generate_quarter_hyper_ball(PointWrapper &center,
                                                    const double  radius)
  {
    AssertThrow(
      dim == spacedim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_quarter_hyper_ball<2>(center, radius, triangulation);
    else
      internal::generate_quarter_hyper_ball<3>(center, radius, triangulation);
  }


  void
  TriangulationWrapper::generate_half_hyper_ball(PointWrapper &center,
                                                 const double  radius)
  {
    AssertThrow(
      dim == spacedim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_half_hyper_ball<2>(center, radius, triangulation);
    else
      internal::generate_half_hyper_ball<3>(center, radius, triangulation);
  }



  void
  TriangulationWrapper::generate_cylinder(const double radius,
                                          const double half_length)
  {
    AssertThrow(
      dim == spacedim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_cylinder<2>(radius, half_length, triangulation);
    else
      internal::generate_cylinder<3>(radius, half_length, triangulation);
  }



  void
  TriangulationWrapper::generate_subdivided_cylinder(
    const unsigned int x_subdivisions,
    const double       radius,
    const double       half_length)
  {
    AssertThrow(
      dim == spacedim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_subdivided_cylinder<2>(radius,
                                                x_subdivisions,
                                                half_length,
                                                triangulation);
    else
      internal::generate_subdivided_cylinder<3>(radius,
                                                x_subdivisions,
                                                half_length,
                                                triangulation);
  }



  void
  TriangulationWrapper::generate_truncated_cone(const double radius_0,
                                                const double radius_1,
                                                const double half_length)
  {
    AssertThrow(
      dim == spacedim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));
    if (dim == 2)
      internal::generate_truncated_cone<2>(radius_0,
                                           radius_1,
                                           half_length,
                                           triangulation);
    else
      internal::generate_truncated_cone<3>(radius_0,
                                           radius_1,
                                           half_length,
                                           triangulation);
  }



  void
  TriangulationWrapper::generate_hyper_cube_with_cylindrical_hole(
    const double       inner_radius,
    const double       outer_radius,
    const double       L,
    const unsigned int repetitions,
    const bool         colorize)
  {
    AssertThrow(
      dim == spacedim,
      ExcMessage(
        "This function is only implemented for dim equal to spacedim."));

    if (dim == 2)
      internal::generate_hyper_cube_with_cylindrical_hole<2>(
        inner_radius, outer_radius, L, repetitions, colorize, triangulation);
    else
      internal::generate_hyper_cube_with_cylindrical_hole<3>(
        inner_radius, outer_radius, L, repetitions, colorize, triangulation);
  }



  void
  TriangulationWrapper::shift(boost::python::list &shift_list)
  {
    AssertThrow(boost::python::len(shift_list) == spacedim,
                ExcMessage(
                  "Size of the shift vector is not equal to spacedim."));
    if ((dim == 2) && (spacedim == 2))
      internal::shift<2, 2>(shift_list, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::shift<2, 3>(shift_list, triangulation);
    else
      internal::shift<3, 3>(shift_list, triangulation);
  }



  void
  TriangulationWrapper::scale(const double scaling_factor)
  {
    if ((dim == 2) && (spacedim == 2))
      internal::scale<2, 2>(scaling_factor, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::scale<2, 3>(scaling_factor, triangulation);
    else
      internal::scale<3, 3>(scaling_factor, triangulation);
  }



  void
  TriangulationWrapper::merge_triangulations(
    boost::python::list &triangulations,
    const double         duplicated_vertex_tolerance,
    const bool           copy_manifold_ids)
  {
    AssertThrow(boost::python::len(triangulations) >= 2,
                ExcMessage(
                  "You must provide at least two triangulations to merge."));

    size_t n_triangulations = boost::python::len(triangulations);
    for (size_t i = 0; i < n_triangulations; ++i)
      {
        TriangulationWrapper &tria_wrapper =
          boost::python::extract<TriangulationWrapper &>(triangulations[i]);

        AssertThrow(dim == tria_wrapper.get_dim(),
                    ExcMessage(
                      "All triangulations must have the same dimension."));
        AssertThrow(
          spacedim == tria_wrapper.get_spacedim(),
          ExcMessage("All triangulations must have the same space dimension."));
      }

    if ((dim == 2) && (spacedim == 2))
      internal::merge_triangulations<2, 2>(triangulations,
                                           triangulation,
                                           duplicated_vertex_tolerance,
                                           copy_manifold_ids);
    else if ((dim == 2) && (spacedim == 3))
      internal::merge_triangulations<2, 3>(triangulations,
                                           triangulation,
                                           duplicated_vertex_tolerance,
                                           copy_manifold_ids);
    else
      internal::merge_triangulations<3, 3>(triangulations,
                                           triangulation,
                                           duplicated_vertex_tolerance,
                                           copy_manifold_ids);
  }



  void
  TriangulationWrapper::flatten_triangulation(TriangulationWrapper &tria_out)
  {
    AssertThrow(
      dim == tria_out.get_dim(),
      ExcMessage(
        "The Triangulation and tria_out should have the same dimension."));
    AssertThrow(spacedim >= tria_out.get_spacedim(),
                ExcMessage(
                  "The Triangulation should have a spacedim greater or equal "
                  "to the spacedim of tria_out."));
    int spacedim_out = tria_out.get_spacedim();
    if ((dim == 2) && (spacedim == 2) && (spacedim_out == 2))
      internal::flatten_triangulation<2, 2, 2>(triangulation, tria_out);
    else if ((dim == 2) && (spacedim == 3) && (spacedim_out == 2))
      internal::flatten_triangulation<2, 3, 2>(triangulation, tria_out);
    else
      internal::flatten_triangulation<3, 3, 3>(triangulation, tria_out);
  }



  void
  TriangulationWrapper::convert_hypercube_to_simplex_mesh(
    TriangulationWrapper &tria_out)
  {
    AssertThrow(
      (tria_out.get_dim() == dim) && (tria_out.get_spacedim() == spacedim),
      ExcMessage("The output Triangulation must be of the same dimension"));

    if ((dim == 2) && (spacedim == 2))
      internal::convert_hypercube_to_simplex_mesh<2, 2>(triangulation,
                                                        tria_out);
    else if ((dim == 2) && (spacedim == 3))
      internal::convert_hypercube_to_simplex_mesh<2, 3>(triangulation,
                                                        tria_out);
    else
      internal::convert_hypercube_to_simplex_mesh<3, 3>(triangulation,
                                                        tria_out);
  }



  void
  TriangulationWrapper::extrude_triangulation(
    const unsigned int    n_slices,
    const double          height,
    TriangulationWrapper &triangulation_out)
  {
    AssertThrow((dim == 2) && (spacedim == 2),
                ExcMessage(
                  "Extrude can only be applied to the dim and spacedim two."));
    AssertThrow((triangulation_out.get_dim() == 3) &&
                  (triangulation_out.get_spacedim() == 3),
                ExcMessage(
                  "The output Triangulation must be of dimension three"));


    Triangulation<2, 2> *tria =
      static_cast<Triangulation<2, 2> *>(triangulation);
    Triangulation<3, 3> *tria_out =
      static_cast<Triangulation<3, 3> *>(triangulation_out.get_triangulation());
    GridGenerator::extrude_triangulation(*tria, n_slices, height, *tria_out);
  }



  void
  TriangulationWrapper::replicate_triangulation(
    TriangulationWrapper &triangulation_in,
    boost::python::list  &extents)
  {
    AssertThrow(
      boost::python::len(extents) == dim,
      ExcMessage(
        "Size of the extents list must be equal to the triangulation dimension."));

    AssertThrow((triangulation_in.get_dim() == dim) &&
                  (triangulation_in.get_spacedim() == spacedim),
                ExcMessage(
                  "The input Triangulation must be of the same dimension"));


    if ((dim == 2) && (spacedim == 2))
      internal::replicate_triangulation<2, 2>(triangulation_in,
                                              extents,
                                              triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::replicate_triangulation<2, 3>(triangulation_in,
                                              extents,
                                              triangulation);
    else
      internal::replicate_triangulation<3, 3>(triangulation_in,
                                              extents,
                                              triangulation);
  }



  void
  TriangulationWrapper::distort_random(const double factor,
                                       const bool   keep_boundary)
  {
    if ((dim == 2) && (spacedim == 2))
      internal::distort_random<2, 2>(factor, keep_boundary, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::distort_random<2, 3>(factor, keep_boundary, triangulation);
    else
      internal::distort_random<3, 3>(factor, keep_boundary, triangulation);
  }



  void
  TriangulationWrapper::transform(boost::python::object &transformation)
  {
    if ((dim == 2) && (spacedim == 2))
      internal::transform<2, 2>(transformation, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::transform<2, 3>(transformation, triangulation);
    else
      internal::transform<3, 3>(transformation, triangulation);
  }



  CellAccessorWrapper
  TriangulationWrapper::find_active_cell_around_point(PointWrapper   &p,
                                                      MappingQWrapper mapping)
  {
    std::pair<int, int> level_index_pair;
    if ((dim == 2) && (spacedim == 2))
      level_index_pair =
        internal::find_active_cell_around_point<2, 2>(p,
                                                      mapping,
                                                      triangulation);
    else if ((dim == 2) && (spacedim == 3))
      level_index_pair =
        internal::find_active_cell_around_point<2, 3>(p,
                                                      mapping,
                                                      triangulation);
    else
      level_index_pair =
        internal::find_active_cell_around_point<3, 3>(p,
                                                      mapping,
                                                      triangulation);

    return CellAccessorWrapper(*this,
                               level_index_pair.first,
                               level_index_pair.second);
  }



  boost::python::list
  TriangulationWrapper::find_cells_adjacent_to_vertex(
    const unsigned int vertex_index)
  {
    if ((dim == 2) && (spacedim == 2))
      return internal::find_cells_adjacent_to_vertex<2, 2>(vertex_index, *this);
    else if ((dim == 2) && (spacedim == 3))
      return internal::find_cells_adjacent_to_vertex<2, 3>(vertex_index, *this);
    else
      return internal::find_cells_adjacent_to_vertex<3, 3>(vertex_index, *this);
  }



  boost::python::list
  TriangulationWrapper::compute_aspect_ratio_of_cells(
    const MappingQWrapper   &mapping,
    const QuadratureWrapper &quadrature)
  {
    if ((dim == 2) && (spacedim == 2))
      return internal::compute_aspect_ratio_of_cells<2, 2>(mapping,
                                                           quadrature,
                                                           *this);
    else if ((dim == 3) && (spacedim == 3))
      return internal::compute_aspect_ratio_of_cells<3, 3>(mapping,
                                                           quadrature,
                                                           *this);
    else
      AssertThrow(false,
                  ExcMessage(
                    "This combination of dim-spacedim is not supported."));
  }



  void
  TriangulationWrapper::refine_global(const unsigned int n)
  {
    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2, 2> *tria =
          static_cast<Triangulation<2, 2> *>(triangulation);
        tria->refine_global(n);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        Triangulation<2, 3> *tria =
          static_cast<Triangulation<2, 3> *>(triangulation);
        tria->refine_global(n);
      }
    else
      {
        Triangulation<3, 3> *tria =
          static_cast<Triangulation<3, 3> *>(triangulation);
        tria->refine_global(n);
      }
  }



  void
  TriangulationWrapper::execute_coarsening_and_refinement()
  {
    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2, 2> *tria =
          static_cast<Triangulation<2, 2> *>(triangulation);
        tria->execute_coarsening_and_refinement();
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        Triangulation<2, 3> *tria =
          static_cast<Triangulation<2, 3> *>(triangulation);
        tria->execute_coarsening_and_refinement();
      }
    else
      {
        Triangulation<3, 3> *tria =
          static_cast<Triangulation<3, 3> *>(triangulation);
        tria->execute_coarsening_and_refinement();
      }
  }



  double
  TriangulationWrapper::minimal_cell_diameter() const
  {
    if ((dim == 2) && (spacedim == 2))
      return internal::minimal_cell_diameter<2, 2>(triangulation);
    else if ((dim == 2) && (spacedim == 3))
      return internal::minimal_cell_diameter<2, 3>(triangulation);
    else
      return internal::minimal_cell_diameter<3, 3>(triangulation);
  }



  double
  TriangulationWrapper::maximal_cell_diameter() const
  {
    if ((dim == 2) && (spacedim == 2))
      return internal::maximal_cell_diameter<2, 2>(triangulation);
    else if ((dim == 2) && (spacedim == 3))
      return internal::maximal_cell_diameter<2, 3>(triangulation);
    else
      return internal::maximal_cell_diameter<3, 3>(triangulation);
  }



  void
  TriangulationWrapper::write(const std::string &filename,
                              const std::string  format) const
  {
    if ((dim == 2) && (spacedim == 2))
      internal::write<2, 2>(filename, format, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::write<2, 3>(filename, format, triangulation);
    else
      internal::write<3, 3>(filename, format, triangulation);
  }



  void
  TriangulationWrapper::read(const std::string &filename,
                             const std::string  format) const
  {
    if ((dim == 2) && (spacedim == 2))
      internal::read<2, 2>(filename, format, triangulation);
    else if ((dim == 2) && (spacedim == 3))
      internal::read<2, 3>(filename, format, triangulation);
    else
      internal::read<3, 3>(filename, format, triangulation);
  }



  void
  TriangulationWrapper::save(const std::string &filename) const
  {
    std::ofstream                   ofs(filename);
    boost::archive::binary_oarchive oa(ofs);

    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2, 2> *tria =
          static_cast<Triangulation<2, 2> *>(triangulation);

        tria->save(oa, 0);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        {
          Triangulation<2, 3> *tria =
            static_cast<Triangulation<2, 3> *>(triangulation);

          tria->save(oa, 0);
        }
      }
    else
      {
        Triangulation<3, 3> *tria =
          static_cast<Triangulation<3, 3> *>(triangulation);

        tria->save(oa, 0);
      }
  }



  void
  TriangulationWrapper::load(const std::string &filename)
  {
    std::ifstream                   ifs(filename);
    boost::archive::binary_iarchive ia(ifs);

    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2, 2> *tria =
          static_cast<Triangulation<2, 2> *>(triangulation);

        tria->load(ia, 0);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        Triangulation<2, 3> *tria =
          static_cast<Triangulation<2, 3> *>(triangulation);

        tria->load(ia, 0);
      }
    else
      {
        Triangulation<3> *tria =
          static_cast<Triangulation<3, 3> *>(triangulation);

        tria->load(ia, 0);
      }
  }



  boost::python::list
  TriangulationWrapper::active_cells()
  {
    if ((dim == 2) && (spacedim == 2))
      return internal::active_cells<2, 2>(*this);
    else if ((dim == 2) && (spacedim == 3))
      return internal::active_cells<2, 3>(*this);
    else
      return internal::active_cells<3, 3>(*this);
  }



  boost::python::list
  TriangulationWrapper::cells()
  {
    if ((dim == 2) && (spacedim == 2))
      return internal::cells<2, 2>(*this);
    else if ((dim == 2) && (spacedim == 3))
      return internal::cells<2, 3>(*this);
    else
      return internal::cells<3, 3>(*this);
  }



  void
  TriangulationWrapper::set_manifold(const int        number,
                                     ManifoldWrapper &manifold)
  {
    AssertThrow(
      dim == manifold.get_dim(),
      ExcMessage(
        "The Triangulation and Manifold should have the same dimension."));
    AssertThrow(
      spacedim == manifold.get_spacedim(),
      ExcMessage(
        "The Triangulation and Manifold should have the same space dimension."));

    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2, 2> *tria =
          static_cast<Triangulation<2, 2> *>(triangulation);
        Manifold<2, 2> *m =
          static_cast<Manifold<2, 2> *>(manifold.get_manifold());
        tria->set_manifold(number, *m);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        Triangulation<2, 3> *tria =
          static_cast<Triangulation<2, 3> *>(triangulation);
        Manifold<2, 3> *m =
          static_cast<Manifold<2, 3> *>(manifold.get_manifold());
        tria->set_manifold(number, *m);
      }
    else
      {
        Triangulation<3, 3> *tria =
          static_cast<Triangulation<3, 3> *>(triangulation);
        Manifold<3, 3> *m =
          static_cast<Manifold<3, 3> *>(manifold.get_manifold());
        tria->set_manifold(number, *m);
      }
  }



  void
  TriangulationWrapper::reset_manifold(const int number)
  {
    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2, 2> *tria =
          static_cast<Triangulation<2, 2> *>(triangulation);
        tria->reset_manifold(number);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        Triangulation<2, 3> *tria =
          static_cast<Triangulation<2, 3> *>(triangulation);
        tria->reset_manifold(number);
      }
    else
      {
        Triangulation<3, 3> *tria =
          static_cast<Triangulation<3, 3> *>(triangulation);
        tria->reset_manifold(number);
      }
  }



  int
  TriangulationWrapper::get_mesh_smoothing()
  {
    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2, 2> *tria =
          static_cast<Triangulation<2, 2> *>(triangulation);
        return (int)tria->get_mesh_smoothing();
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        Triangulation<2, 3> *tria =
          static_cast<Triangulation<2, 3> *>(triangulation);
        return (int)tria->get_mesh_smoothing();
      }
    else
      {
        Triangulation<3, 3> *tria =
          static_cast<Triangulation<3, 3> *>(triangulation);
        return (int)tria->get_mesh_smoothing();
      }
  }



  void
  TriangulationWrapper::set_mesh_smoothing(const int mesh_smoothing)
  {
    if ((dim == 2) && (spacedim == 2))
      {
        Triangulation<2, 2> *tria =
          static_cast<Triangulation<2, 2> *>(triangulation);
        tria->set_mesh_smoothing(
          (Triangulation<2, 2>::MeshSmoothing)mesh_smoothing);
      }
    else if ((dim == 2) && (spacedim == 3))
      {
        Triangulation<2, 3> *tria =
          static_cast<Triangulation<2, 3> *>(triangulation);
        tria->set_mesh_smoothing(
          (Triangulation<2, 3>::MeshSmoothing)mesh_smoothing);
      }
    else
      {
        Triangulation<3, 3> *tria =
          static_cast<Triangulation<3, 3> *>(triangulation);
        tria->set_mesh_smoothing(
          (Triangulation<3, 3>::MeshSmoothing)mesh_smoothing);
      }
  }



  void
  TriangulationWrapper::setup(const std::string &dimension,
                              const std::string &spacedimension,
                              const int          mesh_smoothing,
                              const bool         check_for_distorted_cells)
  {
    if ((dimension.compare("2D") == 0) || (dimension.compare("2d") == 0))
      {
        dim = 2;

        if ((spacedimension.compare("2D") == 0) ||
            (spacedimension.compare("2d") == 0))
          {
            spacedim      = 2;
            triangulation = new Triangulation<2, 2>(
              (Triangulation<2, 2>::MeshSmoothing)mesh_smoothing,
              check_for_distorted_cells);
          }
        else if ((spacedimension.compare("3D") == 0) ||
                 (spacedimension.compare("3d") == 0))
          {
            spacedim      = 3;
            triangulation = new Triangulation<2, 3>(
              (Triangulation<2, 3>::MeshSmoothing)mesh_smoothing,
              check_for_distorted_cells);
          }
        else
          AssertThrow(false,
                      ExcMessage("Spacedimension needs to be 2D or 3D."));
      }
    else if ((dimension.compare("3D") == 0) || (dimension.compare("3d") == 0))
      {
        if ((spacedimension.compare("3D") != 0) &&
            (spacedimension.compare("3d") != 0))
          AssertThrow(false, ExcMessage("Spacedimension needs to be 3D."));
        dim           = 3;
        spacedim      = 3;
        triangulation = new Triangulation<3, 3>(
          (Triangulation<3, 3>::MeshSmoothing)mesh_smoothing,
          check_for_distorted_cells);
      }
    else
      AssertThrow(false, ExcMessage("Dimension needs to be 2D or 3D."));
  }

} // namespace python

DEAL_II_NAMESPACE_CLOSE
