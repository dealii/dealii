#include <deal.II/gmsh/utilities.h>
#include <deal.II/opencascade/utilities.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <cstdio>

#ifdef DEAL_II_WITH_GMSH

DEAL_II_NAMESPACE_OPEN

namespace Gmsh
{
  AdditionalParameters::AdditionalParameters(const double characteristic_length,
                                             const std::string &output_base_name):
    characteristic_length(characteristic_length),
    output_base_name(output_base_name)
  {}



  void AdditionalParameters::add_parameters(ParameterHandler &prm)
  {
    prm.add_parameter("Characteristic length", characteristic_length);
    prm.add_parameter("Intermediate file name base", output_base_name,
                      "Keep empty, if you want the program to generate "
                      "temporary files, and then remove them when they "
                      "are no longer used.");
  }



#ifdef DEAL_II_WITH_OPENCASCADE
  template <int spacedim>
  void
  create_triangulation_from_boundary_curve(const TopoDS_Edge &boundary,
                                           Triangulation<2,spacedim> &tria,
                                           const AdditionalParameters &prm)
  {
    std::string base_name = prm.output_base_name;
    if (base_name == "")
      base_name = std::tmpnam(nullptr);

    dealii::OpenCASCADE::write_IGES(boundary, base_name+".iges");

    ofstream geofile;
    geofile.open( base_name+".geo");
    geofile << "Merge \"" << base_name << ".iges\";" << std::endl
            << "Line Loop (2) = {1};" << std::endl
            << "Plane Surface (3) = {2};" << std::endl
            << "Characteristic Length { 1 } = " << prm.characteristic_length << ";" << std::endl
            << "Mesh.RecombineAll = 1;" << std::endl
            << "Mesh.SubdivisionAlgorithm=1;" << std::endl;

    geofile.close();

    // std::system("gmsh -2 -algo front3d temp_model.geo 1>temp_out.log 2>temp_warn.log");
    std::stringstream command;
    command << DEAL_II_GMSH_EXECUTABLE_PATH << " -2 "
            << base_name << ".geo 1> "
            << base_name << ".log 2> "
            << base_name << "_warn.log";

    auto ret_value = std::system(command.str().c_str());
    AssertThrow(ret_value == 0,
                ExcMessage("Gmsh failed to run. Check the "+base_name+".log file."));

    std::ifstream grid_file(base_name+".msh");
    Assert(grid_file, ExcIO());

    GridIn<2,spacedim> gridin;
    gridin.attach_triangulation(tria);
    gridin.read_msh(grid_file);

    if (base_name != prm.output_base_name)
      {
        std::remove((base_name + ".geo").c_str());
        std::remove((base_name + ".log").c_str());
        std::remove((base_name + "_warn.log").c_str());
        std::remove((base_name + ".msh").c_str());
        std::remove((base_name + ".iges").c_str());
      }
  }
#endif

#include "utilities.inst"
}

DEAL_II_NAMESPACE_CLOSE

#endif
