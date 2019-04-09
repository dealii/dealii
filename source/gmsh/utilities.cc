// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

#include <deal.II/gmsh/utilities.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/opencascade/utilities.h>

#include <cstdio>

#ifdef DEAL_II_WITH_GMSH

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
  create_triangulation_from_boundary_curve(const TopoDS_Edge &         boundary,
                                           Triangulation<2, spacedim> &tria,
                                           const AdditionalParameters &prm)
  {
    std::string base_name      = prm.output_base_name;
    char        dir_template[] = "ctfbc-XXXXXX";
    if (base_name.empty())
      {
        const char *temp = mkdtemp(dir_template);
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

    ofstream geofile;
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
    Assert(grid_file, ExcIO());

    GridIn<2, spacedim> gridin;
    gridin.attach_triangulation(tria);
    gridin.read_msh(grid_file);

    if (base_name != prm.output_base_name)
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
        const auto ret_value = std::remove(dir_template);
        AssertThrow(ret_value == 0,
                    ExcMessage("Failed to remove " +
                               std::string(dir_template)));
      }
  }
#  endif

  // explicit instantiations
#  ifdef DEAL_II_WITH_OPENCASCADE
#    include "utilities.inst"
#  endif
} // namespace Gmsh

DEAL_II_NAMESPACE_CLOSE

#endif
