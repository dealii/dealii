// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2020 by the deal.II authors
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


#ifndef dealii_gmsh_parameters_h
#define dealii_gmsh_parameters_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_GMSH

#  ifdef DEAL_II_WITH_OPENCASCADE
#    include <TopoDS_CompSolid.hxx>
#    include <TopoDS_Compound.hxx>
#    include <TopoDS_Edge.hxx>
#    include <TopoDS_Face.hxx>
#    include <TopoDS_Shape.hxx>
#    include <TopoDS_Vertex.hxx>
#  endif

#  include <deal.II/base/parameter_handler.h>

#  include <deal.II/grid/tria.h>

DEAL_II_NAMESPACE_OPEN

/**
 * A collection of %Gmsh related utilities and classes.
 *
 * @author Luca Heltai, Dirk Peschka, 2018
 */
namespace Gmsh
{
  /**
   * A parameter class used to pass options to the %Gmsh executable.
   *
   * @author Luca Heltai, 2018
   */
  class AdditionalParameters
  {
  public:
    /**
     * Set all additional parameters to their default values.
     */
    AdditionalParameters(const double       characteristic_length = 1.0,
                         const std::string &output_base_name      = "");

    /**
     * Call prm.add_parameter for each member of the AdditionalParameters class.
     */
    void
    add_parameters(ParameterHandler &prm);

    /**
     * The characteristic length used for the definition of the %Gmsh grid.
     *
     * %Gmsh will try to make sure that the size of each edge is as close as
     * possible to this value.
     */
    double characteristic_length = 1.0;

    /**
     * Basename for the output files.
     *
     * If this is left empty, then temporary files are used, and removed when
     * not needed any more.
     */
    std::string output_base_name = "";
  };

#  ifdef DEAL_II_WITH_OPENCASCADE
  /**
   * Given a smooth closed curve, create a triangulation from it using
   * %Gmsh.
   *
   * The input curve @p boundary should be closed.
   *
   * @authors Luca Heltai, Dirk Peschka, 2018
   */
  template <int spacedim>
  void
  create_triangulation_from_boundary_curve(
    const TopoDS_Edge &         boundary,
    Triangulation<2, spacedim> &tria,
    const AdditionalParameters &prm = AdditionalParameters());
#  endif
} // namespace Gmsh

DEAL_II_NAMESPACE_CLOSE

#endif
#endif
