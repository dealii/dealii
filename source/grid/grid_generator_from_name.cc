// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/mutable_bind.h>

#include <deal.II/grid/grid_generator.h>

DEAL_II_NAMESPACE_OPEN

namespace GridGenerator
{
  namespace
  {
    /**
     * Given a GridGenerator function pointer, a  string containing a text
     * version of the function arguments and an empty Triangulation, calls the
     * corresponding function, after parsing the arguments from the given
     * string.
     */
    template <int dim, int spacedim, class... Arguments>
    void
    parse_and_create(void (*generator)(Triangulation<dim, spacedim> &,
                                       Arguments...),
                     const std::string            &arguments,
                     Triangulation<dim, spacedim> &tria)
    {
      std::function<void(Arguments...)> wrapper =
        [&tria, &generator](Arguments... args) { generator(tria, args...); };
      auto bound_function = Utilities::mutable_bind(wrapper);
      bound_function.parse_arguments(arguments);
      bound_function();
    }


    /**
     * Create grids that only exist for codimension zero combinations.
     *
     * Return true if a grid was actually generated, false otherwise.
     */
    template <int dim, int spacedim>
    std::enable_if_t<dim != spacedim, bool>
    generate_codimension_zero_grid(const std::string &,
                                   const std::string &,
                                   Triangulation<dim, spacedim> &)
    {
      return false;
    }

    /**
     * Create grids that only exist for codimension zero combinations.
     *
     * Return true if a grid was actually generated, false otherwise.
     */
    template <int dim>
    bool
    generate_codimension_zero_grid(const std::string  &name,
                                   const std::string  &arguments,
                                   Triangulation<dim> &tria)
    {
      if (name == "simplex")
        parse_and_create<dim, dim, const std::vector<Point<dim>> &>(simplex,
                                                                    arguments,
                                                                    tria);
      else if (name == "subdivided_hyper_rectangle")
        {
          // subdivided_hyper_rectangle is polymorphic, and can be called with
          // different sets of arguments. We support two of these function
          // calls. We try the first one, and if parsing fails, we go to the
          // second one.
          try
            {
              parse_and_create<dim,
                               dim,
                               const std::vector<unsigned int> &,
                               const Point<dim> &,
                               const Point<dim> &,
                               bool>(subdivided_hyper_rectangle,
                                     arguments,
                                     tria);
            }
          catch (Patterns::Tools::ExcNoMatch &)
            {
              parse_and_create<dim,
                               dim,
                               const std::vector<std::vector<double>> &,
                               const Point<dim> &,
                               const Point<dim> &,
                               bool>(subdivided_hyper_rectangle,
                                     arguments,
                                     tria);
            }
        }
      else if (name == "plate_with_a_hole")
        parse_and_create<dim,
                         dim,
                         double,
                         double,
                         double,
                         double,
                         double,
                         double,
                         const Point<dim> &,
                         types::manifold_id,
                         types::manifold_id,
                         double,
                         unsigned int,
                         bool>(plate_with_a_hole, arguments, tria);
      else if (name == "channel_with_cylinder")
        parse_and_create<dim, dim, double, unsigned int, double, bool>(
          channel_with_cylinder, arguments, tria);
      else if (name == "uniform_channel_with_cylinder")
        parse_and_create<dim,
                         dim,
                         const std::vector<unsigned int> &,
                         double,
                         unsigned int,
                         double,
                         unsigned int,
                         double,
                         bool,
                         bool>(uniform_channel_with_cylinder, arguments, tria);
      else if (name == "enclosed_hyper_cube")
        parse_and_create<dim, dim, double, double, double, bool>(
          enclosed_hyper_cube, arguments, tria);

      else if (name == "hyper_ball")
        parse_and_create<dim, dim, const Point<dim> &, double, bool>(hyper_ball,
                                                                     arguments,
                                                                     tria);
      else if (name == "hyper_ball_balanced")
        parse_and_create<dim, dim, const Point<dim> &, double>(
          hyper_ball_balanced, arguments, tria);

      else if (name == "quarter_hyper_ball")
        parse_and_create<dim, dim, const Point<dim> &, double>(
          quarter_hyper_ball, arguments, tria);

      else if (name == "half_hyper_ball")
        parse_and_create<dim, dim, const Point<dim> &, double>(half_hyper_ball,
                                                               arguments,
                                                               tria);

      else if (name == "cylinder")
        parse_and_create<dim, dim, double, double>(cylinder, arguments, tria);

      else if (name == "subdivided_cylinder")
        parse_and_create<dim, dim, unsigned int, double, double>(
          subdivided_cylinder, arguments, tria);

      else if (name == "truncated_cone")
        parse_and_create<dim, dim, double, double, double>(truncated_cone,
                                                           arguments,
                                                           tria);

      else if (name == "pipe_junction")
        parse_and_create<dim,
                         dim,
                         const std::vector<std::pair<Point<dim>, double>> &,
                         const std::pair<Point<dim>, double> &,
                         double>(pipe_junction, arguments, tria);

      else if (name == "hyper_L")
        parse_and_create<dim, dim, double, double, bool>(hyper_L,
                                                         arguments,
                                                         tria);

      else if (name == "hyper_cube_slit")
        parse_and_create<dim, dim, double, double, bool>(hyper_cube_slit,
                                                         arguments,
                                                         tria);

      else if (name == "hyper_shell")
        parse_and_create<dim,
                         dim,
                         const Point<dim> &,
                         double,
                         double,
                         unsigned int,
                         bool>(hyper_shell, arguments, tria);

      else if (name == "half_hyper_shell")
        parse_and_create<dim,
                         dim,
                         const Point<dim> &,
                         double,
                         double,
                         unsigned int,
                         bool>(half_hyper_shell, arguments, tria);

      else if (name == "quarter_hyper_shell")
        parse_and_create<dim,
                         dim,
                         const Point<dim> &,
                         double,
                         double,
                         unsigned int,
                         bool>(quarter_hyper_shell, arguments, tria);

      else if (name == "eccentric_hyper_shell")
        parse_and_create<dim,
                         dim,
                         const Point<dim> &,
                         const Point<dim> &,
                         double,
                         double,
                         unsigned int>(eccentric_hyper_shell, arguments, tria);

      else if (name == "cylinder_shell")
        parse_and_create<dim,
                         dim,
                         double,
                         double,
                         double,
                         unsigned int,
                         unsigned int,
                         bool>(cylinder_shell, arguments, tria);

      else if (name == "hyper_cube_with_cylindrical_hole")
        parse_and_create<dim, dim, double, double, double, unsigned int, bool>(
          hyper_cube_with_cylindrical_hole, arguments, tria);

      else if (name == "concentric_hyper_shells")
        parse_and_create<dim,
                         dim,
                         const Point<dim> &,
                         double,
                         double,
                         unsigned int,
                         double,
                         unsigned int,
                         bool>(concentric_hyper_shells, arguments, tria);

      else if (name == "subdivided_hyper_cube_with_simplices")
        parse_and_create<dim, dim, unsigned int, double, double, bool>(
          subdivided_hyper_cube_with_simplices, arguments, tria);

      else if (name == "subdivided_hyper_rectangle_with_simplices")
        parse_and_create<dim,
                         dim,
                         const std::vector<unsigned int> &,
                         const Point<dim> &,
                         const Point<dim> &,
                         bool>(subdivided_hyper_rectangle_with_simplices,
                               arguments,
                               tria);

      else if (name == "subdivided_hyper_L")
        parse_and_create<dim,
                         dim,
                         const std::vector<unsigned int> &,
                         const Point<dim> &,
                         const Point<dim> &,
                         const std::vector<int> &>(subdivided_hyper_L,
                                                   arguments,
                                                   tria);

      else
        return false;

      return true;
    }


    /**
     * Create grids that only exist for codimension one combinations.
     *
     * Return true if a grid was actually generated, false otherwise.
     */
    template <int dim, int spacedim>
    std::enable_if_t<dim != spacedim - 1, bool>
    generate_codimension_one_grid(const std::string &,
                                  const std::string &,
                                  Triangulation<dim, spacedim> &)
    {
      return false;
    }

    /**
     * Create grids that only exist for codimension one combinations.
     *
     * Return true if a grid was actually generated, false otherwise.
     */
    template <int dim>
    bool
    generate_codimension_one_grid(const std::string           &name,
                                  const std::string           &arguments,
                                  Triangulation<dim, dim + 1> &tria)
    {
      if (name == "hyper_sphere")
        parse_and_create<dim, dim + 1, const Point<dim + 1> &, double>(
          hyper_sphere, arguments, tria);
      else
        return false;
      return true;
    }

    /**
     * Methods only implemented for special combinations of dim and spacedim.
     *
     * Return true if a grid was actually generated, false otherwise.
     */
    template <int dim, int spacedim>
    bool
    generate_special(const std::string &,
                     const std::string &,
                     Triangulation<dim, spacedim> &)
    {
      return false;
    }

    /**
     * Methods implemented only for Triangulation<3,3>.
     *
     * Return true if a grid was actually generated, false otherwise.
     */
    bool
    generate_special(const std::string   &name,
                     const std::string   &arguments,
                     Triangulation<3, 3> &tria)
    {
      if (name == "moebius")
        parse_and_create<3, 3, unsigned int, unsigned int, double, double>(
          moebius, arguments, tria);
      else if (name == "torus")
        parse_and_create<3, 3, double, double, unsigned int, double>(torus,
                                                                     arguments,
                                                                     tria);
      else
        {
          return false;
        }
      return true;
    }

    /**
     * Methods implemented only for Triangulation<2,3>.
     *
     * Return true if a grid was actually generated, false otherwise.
     */
    bool
    generate_special(const std::string   &name,
                     const std::string   &arguments,
                     Triangulation<2, 3> &tria)
    {
      if (name == "torus")
        parse_and_create<2, 3, double, double, unsigned int, double>(torus,
                                                                     arguments,
                                                                     tria);
      else
        {
          return false;
        }
      return true;
    }
  } // namespace



  template <int dim, int spacedim>
  void
  generate_from_name_and_arguments(Triangulation<dim, spacedim> &tria,
                                   const std::string            &name,
                                   const std::string            &arguments)
  {
    // We begin with all function calls that are implemented for all
    // combinations of dim and spacedim.
    if (name == "hyper_cube")
      parse_and_create<dim, spacedim, double, double, bool>(hyper_cube,
                                                            arguments,
                                                            tria);
    else if (name == "subdivided_hyper_cube")
      parse_and_create<dim, spacedim, unsigned int, double, double, bool>(
        subdivided_hyper_cube, arguments, tria);
    else if (name == "hyper_rectangle")
      parse_and_create<dim,
                       spacedim,
                       const Point<dim> &,
                       const Point<dim> &,
                       bool>(hyper_rectangle, arguments, tria);
    else if (name == "cheese")
      parse_and_create<dim, spacedim, const std::vector<unsigned int> &>(
        cheese, arguments, tria);
    else if (name == "general_cell")
      parse_and_create<dim,
                       spacedim,
                       const std::vector<Point<spacedim>> &,
                       bool>(general_cell, arguments, tria);
    else if (name == "hyper_cross")
      parse_and_create<dim, spacedim, const std::vector<unsigned int> &, bool>(
        hyper_cross, arguments, tria);
    // If none of the above worked, than we try with more specific function
    // calls. First we try to call functions that are only implemented when
    // dim == spacedim, then when dim == spacedim-1, and lastly, we try to see
    // if the name, dim, and spacedim match some of the very special grid
    // generator functions, like torus, moebius, etc.
    //
    // If one of the function call succeeds, we skip the rest and return.
    else if (generate_codimension_zero_grid(name, arguments, tria))
      {
      }
    else if (generate_codimension_one_grid(name, arguments, tria))
      {
      }
    else if (generate_special(name, arguments, tria))
      {
      }
    else
      // If we got here, we really have no idea what grid the user wants to
      // generate.
      AssertThrow(false,
                  ExcMessage(name + "(" + arguments + ") not implemented"));
  }
} // namespace GridGenerator

#include "grid/grid_generator_from_name.inst"

DEAL_II_NAMESPACE_CLOSE
