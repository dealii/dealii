// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <deal.II/gmsh/utilities.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/opencascade/utilities.h>

#include <cstdio>

#ifdef DEAL_II_WITH_GMSH

#  ifdef DEAL_II_WITH_OPENCASCADE
#    include <TopoDS_Edge.hxx>
#  endif

DEAL_II_NAMESPACE_OPEN

namespace Gmsh
{
  AdditionalParameters::AdditionalParameters(
    const double       characteristic_length,
    const std::string &output_base_name)
    : characteristic_length(characteristic_length)
    , output_base_name(output_base_name)
  {}



  void
  AdditionalParameters::add_parameters(ParameterHandler &prm)
  {
    prm.add_parameter("Characteristic length", characteristic_length);
    prm.add_parameter("Intermediate file name base",
                      output_base_name,
                      "Keep empty, if you want the program to generate "
                      "temporary files, and then remove them when they "
                      "are no longer used.");
  }



#  ifdef DEAL_II_WITH_OPENCASCADE
  template <int spacedim>
  void
  create_triangulation_from_boundary_curve(const TopoDS_Edge          &boundary,
                                           Triangulation<2, spacedim> &tria,
                                           const AdditionalParameters &prm)
  {
    std::string base_name = prm.output_base_name;

    // If necessary, create a temp directory to put files into. The
    // following variable will hold the name of the tmp dir; it is
    // initialized with a template, and 'mkdtemp' then overwrites it
    // with the name of the directory it creates.
    char tmp_dir_name[] = "ctfbc-XXXXXX";
    if (prm.output_base_name.empty())
      {
        const char *temp = mkdtemp(tmp_dir_name);
        AssertThrow(temp != nullptr,
                    ExcMessage("Creating temporary directory failed!"));
        base_name = temp;
        base_name += "tmp";
      }

    const std::string iges_file_name     = base_name + ".iges";
    const std::string geo_file_name      = base_name + ".geo";
    const std::string msh_file_name      = base_name + ".msh";
    const std::string log_file_name      = base_name + ".log";
    const std::string warnings_file_name = base_name + "_warn.log";

    dealii::OpenCASCADE::write_IGES(boundary, iges_file_name);

    std::ofstream geofile;
    geofile.open(geo_file_name);
    geofile << "Merge \"" << iges_file_name << "\";" << std::endl
            << "Line Loop (2) = {1};" << std::endl
            << "Plane Surface (3) = {2};" << std::endl
            << "Characteristic Length { 1 } = " << prm.characteristic_length
            << ";" << std::endl
            << "Mesh.RecombineAll = 1;" << std::endl
            << "Mesh.SubdivisionAlgorithm = 1;" << std::endl;
    geofile.close();

    std::stringstream command;
    command << DEAL_II_GMSH_EXECUTABLE_PATH << " -2 " << geo_file_name << " 1> "
            << log_file_name << " 2> " << warnings_file_name;

    const auto ret_value = std::system(command.str().c_str());
    AssertThrow(ret_value == 0,
                ExcMessage("Gmsh failed to run. Check the " + log_file_name +
                           " file."));

    std::ifstream grid_file(msh_file_name);
    Assert(grid_file.fail() == false, ExcIO());

    GridIn<2, spacedim> gridin;
    gridin.attach_triangulation(tria);
    gridin.read_msh(grid_file);

    // Clean up files if a tmp directory was used:
    if (prm.output_base_name.empty())
      {
        // declaring the list without a type, i.e.,
        //
        //     auto filenames = {{iges_file_name, geo_file_name, ...}})
        //
        // causes internal compiler errors with GCC's concepts implementation,
        // so give it an explicit type:
        const std::array<const std::string *, 5> filenames{
          {&iges_file_name,
           &geo_file_name,
           &msh_file_name,
           &log_file_name,
           &warnings_file_name}};
        for (const std::string *filename : filenames)
          {
            const auto ret_value = std::remove(filename->c_str());
            AssertThrow(ret_value == 0,
                        ExcMessage("Failed to remove " + *filename));
          }

        // Finally also remove the tmp directory:
        const auto ret_value = std::remove(tmp_dir_name);
        AssertThrow(ret_value == 0,
                    ExcMessage("Failed to remove " +
                               std::string(tmp_dir_name)));
      }
  }
#  endif

  // explicit instantiations
#  ifdef DEAL_II_WITH_OPENCASCADE
// We don't build the utilities.inst file if deal.II isn't configured
// with GMSH, but doxygen doesn't know that and tries to find that
// file anyway for parsing -- which then of course it fails on. So
// exclude the following from doxygen consideration.
#    ifndef DOXYGEN
#      include "gmsh/utilities.inst"
#    endif
#  endif
} // namespace Gmsh

DEAL_II_NAMESPACE_CLOSE

#endif
