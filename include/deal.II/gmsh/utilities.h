// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


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

#  include <deal.II/grid/tria.h>

DEAL_II_NAMESPACE_OPEN

#  ifndef DOXYGEN
class ParameterHandler;
#  endif

/**
 * A collection of %Gmsh related utilities and classes.
 */
namespace Gmsh
{
  /**
   * A parameter class used to pass options to the %Gmsh executable.
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
     * Base name for the output files. The base name may contain a directory
     * (followed by a slash), and must contain the base of the names of files
     * to be created.
     *
     * If this variable is left empty, then a temporary directory will be used,
     * and both the files written into it as well as the temporary directory
     * will be removed when not needed any more.
     */
    std::string output_base_name = {};
  };

#  ifdef DEAL_II_WITH_OPENCASCADE
  /**
   * Given a smooth closed curve, create a triangulation from it using
   * %Gmsh.
   *
   * The input curve @p boundary should be closed.
   */
  template <int spacedim>
  void
  create_triangulation_from_boundary_curve(
    const TopoDS_Edge          &boundary,
    Triangulation<2, spacedim> &tria,
    const AdditionalParameters &prm = AdditionalParameters());
#  endif
} // namespace Gmsh

DEAL_II_NAMESPACE_CLOSE

#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif
#endif
