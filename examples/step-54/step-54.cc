/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2009 - 2019 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 *  Authors: Andrea Mola, Luca Heltai, 2014
 */


// @sect3{Include files}

// We start with including a bunch of files that we will use in the
// various parts of the program. Most of them have been discussed in
// previous tutorials already:
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// These are the headers of the opencascade support classes and
// functions. Notice that these will contain sensible data only if you
// compiled your deal.II library with support for OpenCASCADE, i.e.,
// specifying <code>-DDEAL_II_WITH_OPENCASCADE=ON</code> and
// <code>-DOPENCASCADE_DIR=/path/to/your/opencascade/installation</code>
// when calling <code>cmake</code> during deal.II configuration.
#include <deal.II/opencascade/manifold_lib.h>
#include <deal.II/opencascade/utilities.h>


// Finally, a few C++ standard header files
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

// We isolate the rest of the program in its own namespace
namespace Step54
{
  using namespace dealii;



  // @sect3{The TriangulationOnCAD class}

  // This is the main class. All it really does is store names for
  // input and output files, and a triangulation. It then provides
  // a function that generates such a triangulation from a coarse
  // mesh, using one of the strategies discussed in the introduction
  // and listed in the enumeration type at the top of the class.
  //
  // The member functions of this class are similar to what you can
  // find in most of the other tutorial programs in the setup stage of
  // the grid for the simulations.

  class TriangulationOnCAD
  {
  public:
    enum ProjectionType
    {
      NormalProjection       = 0,
      DirectionalProjection  = 1,
      NormalToMeshProjection = 2
    };


    TriangulationOnCAD(
      const std::string &  initial_mesh_filename,
      const std::string &  cad_file_name,
      const std::string &  output_filename,
      const ProjectionType surface_projection_kind = NormalProjection);

    void run();

  private:
    void read_domain();

    void refine_mesh();

    void output_results(const unsigned int cycle);

    Triangulation<2, 3> tria;

    const std::string initial_mesh_filename;
    const std::string cad_file_name;
    const std::string output_filename;

    const ProjectionType surface_projection_kind;
  };


  // @sect4{TriangulationOnCAD::TriangulationOnCAD}

  // The constructor of the TriangulationOnCAD class is very simple.
  // The input arguments are strings for the input and output file
  // names, and the enumeration type that determines which kind of
  // surface projector is used in the mesh refinement cycles (see
  // below for details).

  TriangulationOnCAD::TriangulationOnCAD(
    const std::string &  initial_mesh_filename,
    const std::string &  cad_file_name,
    const std::string &  output_filename,
    const ProjectionType surface_projection_kind)
    : initial_mesh_filename(initial_mesh_filename)
    , cad_file_name(cad_file_name)
    , output_filename(output_filename)
    , surface_projection_kind(surface_projection_kind)
  {}


  // @sect4{TriangulationOnCAD::read_domain}


  // The following function represents the core of this program.  In
  // this function we import the CAD shape upon which we want to
  // generate and refine our triangulation. We assume that the CAD
  // surface is contained in the @p cad_file_name file (we provide an
  // example IGES file in the input directory called
  // "input/DTMB-5415_bulbous_bow.iges" that represents the bulbous bow of a
  // ship). The presence of several convex and concave high curvature
  // regions makes the geometry we provided a particularly meaningful
  // example.
  //
  // After importing the hull bow surface, we extract some of the
  // curves and surfaces composing it, and use them to generate a set
  // of projectors. Such projectors define the rules the Triangulation
  // has to follow to position each new node during cell refinement.
  //
  // To initialize the Triangulation, as done in previous tutorial
  // programs, we import a pre-existing grid saved in VTK format. We
  // assume here that the user has generated a coarse mesh
  // externally, which matches the IGES geometry. At the moment of
  // writing this tutorial, the
  // deal.II library does not automatically support generation of such
  // meshes, but there are several tools which can provide you with
  // reasonable initial meshes starting from CAD files.
  // In our example, the imported mesh is composed of a single
  // quadrilateral cell whose vertices have been placed on the CAD
  // shape.
  //
  // After importing both the IGES geometry and the initial mesh, we
  // assign the projectors previously discussed to each of the edges
  // and cells which will have to be refined on the CAD surface.
  //
  // In this tutorial, we will test the three different CAD surface
  // projectors described in the introduction, and will analyze the
  // results obtained with each of them.  As mentioned, each of these
  // projection strategies has been implemented in a different class,
  // and objects of these types can be assigned to a triangulation
  // using the Triangulation::set_manifold method.
  //
  // The following function then first imports the given CAD file.
  // The function arguments are a string containing the desired file
  // name, and a scale factor. In this example, the scale factor is
  // set to 1e-3, as the original geometry is written in millimeters
  // (which is the typical unit of measure for most IGES files),
  // while we prefer to work in meters.  The output of the function
  // is an object of OpenCASCADE generic topological shape class,
  // namely a @p TopoDS_Shape.
  void TriangulationOnCAD::read_domain()
  {
    TopoDS_Shape bow_surface = OpenCASCADE::read_IGES(cad_file_name, 1e-3);

    // Each CAD geometrical object is defined along with a tolerance,
    // which indicates possible inaccuracy of its placement. For
    // instance, the tolerance @p tol of a vertex indicates that it can
    // be located in any point contained in a sphere centered in the
    // nominal position and having radius @p tol. While projecting a
    // point onto a surface (which will in turn have its tolerance) we
    // must keep in mind that the precision of the projection will be
    // limited by the tolerance with which the surface is built.

    // The following method extracts the tolerance of the given shape and
    // makes it a bit bigger to stay our of trouble:
    const double tolerance = OpenCASCADE::get_shape_tolerance(bow_surface) * 5;

    // We now want to extract a set of composite sub-shapes from the
    // generic shape. In particular, each face of the CAD file
    // is composed of a trimming curve of type @p TopoDS_Wire, which is
    // the collection of @p TopoDS_Edges that compose the boundary of a
    // surface, and a NURBS description of the surface itself. We will
    // use a line projector to associate the boundary of our
    // Triangulation to the wire delimiting the surface.  To extract
    // all compound sub-shapes, like wires, shells, or solids, we
    // resort to a method of the OpenCASCADE namespace.  The input of
    // OpenCASCADE::extract_compound_shapes is a shape and a set of empty
    // std::vectors of subshapes, which will be filled with all
    // compound shapes found in the given topological shape:
    std::vector<TopoDS_Compound>  compounds;
    std::vector<TopoDS_CompSolid> compsolids;
    std::vector<TopoDS_Solid>     solids;
    std::vector<TopoDS_Shell>     shells;
    std::vector<TopoDS_Wire>      wires;

    OpenCASCADE::extract_compound_shapes(
      bow_surface, compounds, compsolids, solids, shells, wires);

    // The next few steps are more familiar, and allow us to import an existing
    // mesh from an external VTK file, and convert it to a deal triangulation.
    std::ifstream in;

    in.open(initial_mesh_filename);

    GridIn<2, 3> gi;
    gi.attach_triangulation(tria);
    gi.read_vtk(in);

    // We output this initial mesh saving it as the refinement step 0.
    output_results(0);

    // The mesh imported has a single, two-dimensional cell located in
    // three-dimensional space. We now want to ensure that it is refined
    // according to the CAD geometry imported above. This this end, we get an
    // iterator to that cell and assign to it the manifold_id 1 (see
    // @ref GlossManifoldIndicator "this glossary entry").
    // We also get an iterator to its four faces, and assign each of them
    // the manifold_id 2:
    Triangulation<2, 3>::active_cell_iterator cell = tria.begin_active();
    cell->set_manifold_id(1);

    for (const auto &face : cell->face_iterators())
      face->set_manifold_id(2);

    // Once both the CAD geometry and the initial mesh have been
    // imported and digested, we use the CAD surfaces and curves to
    // define the projectors and assign them to the manifold ids just
    // specified.

    // A first projector is defined using the single wire contained in
    // our CAD file.  The ArclengthProjectionLineManifold will make
    // sure that every mesh edge located on the wire is refined with a
    // point that lies on the wire and splits it into two equal arcs
    // lying between the edge vertices. We first check
    // that the wires vector contains at least one element and then
    // create a Manifold object for it.
    //
    // Once the projector is created, we then assign it to all the parts of
    // the triangulation with manifold_id = 2:
    Assert(
      wires.size() > 0,
      ExcMessage(
        "I could not find any wire in the CAD file you gave me. Bailing out."));

    OpenCASCADE::ArclengthProjectionLineManifold<2, 3> line_projector(
      wires[0], tolerance);

    tria.set_manifold(2, line_projector);

    // The surface projector is created according to what is specified
    // with the @p surface_projection_kind option of the constructor. In particular,
    // if the surface_projection_kind value equals @p NormalProjection, we select the
    // OpenCASCADE::NormalProjectionManifold. The new mesh points will
    // then initially be generated at the barycenter of the cell/edge
    // considered, and then projected on the CAD surface along its
    // normal direction.  The NormalProjectionManifold constructor
    // only needs a shape and a tolerance, and we then assign it to
    // the triangulation for use with all parts that manifold having id 1:
    switch (surface_projection_kind)
      {
        case NormalProjection:
          {
            OpenCASCADE::NormalProjectionManifold<2, 3> normal_projector(
              bow_surface, tolerance);
            tria.set_manifold(1, normal_projector);

            break;
          }

        // @p If surface_projection_kind value is @p DirectionalProjection, we select the
        // OpenCASCADE::DirectionalProjectionManifold class. The new mesh points
        // will then initially be generated at the barycenter of the cell/edge
        // considered, and then projected on the CAD surface along a
        // direction that is specified to the
        // OpenCASCADE::DirectionalProjectionManifold constructor. In this case,
        // the projection is done along the y-axis.
        case DirectionalProjection:
          {
            OpenCASCADE::DirectionalProjectionManifold<2, 3>
              directional_projector(bow_surface,
                                    Point<3>(0.0, 1.0, 0.0),
                                    tolerance);
            tria.set_manifold(1, directional_projector);

            break;
          }

        // As a third option, if @p surface_projection_kind value
        // is @p NormalToMeshProjection, we select the
        // OpenCASCADE::NormalToMeshProjectionManifold. The new mesh points will
        // again initially be generated at the barycenter of the cell/edge
        // considered, and then projected on the CAD surface along a
        // direction that is an estimate of the mesh normal direction.
        // The OpenCASCADE::NormalToMeshProjectionManifold constructor only
        // requires a shape (containing at least a face) and a
        // tolerance.
        case NormalToMeshProjection:
          {
            OpenCASCADE::NormalToMeshProjectionManifold<2, 3>
              normal_to_mesh_projector(bow_surface, tolerance);
            tria.set_manifold(1, normal_to_mesh_projector);

            break;
          }

        // Finally, we use good software cleanliness by ensuring that this
        // really covers all possible options of the @p case statement. If we
        // get any other value, we simply abort the program:
        default:
          AssertThrow(false, ExcInternalError());
      }
  }


  // @sect4{TriangulationOnCAD::refine_mesh}

  // This function globally refines the mesh. In other tutorials, it
  // would typically also distribute degrees of freedom, and resize
  // matrices and vectors. These tasks are not carried out here, since
  // we are not running any simulation on the Triangulation produced.
  //
  // While the function looks innocent, this is where most of the work we are
  // interested in for this tutorial program actually happens. In particular,
  // when refining the quads and lines that define the surface of the ship's
  // hull, the Triangulation class will ask the various objects we have
  // assigned to handle individual manifold ids for where the new vertices
  // should lie.
  void TriangulationOnCAD::refine_mesh()
  {
    tria.refine_global(1);
  }



  // @sect4{TriangulationOnCAD::output_results}

  // Outputting the results of our computations is a rather mechanical
  // task. All the components of this function have been discussed
  // before:
  void TriangulationOnCAD::output_results(const unsigned int cycle)
  {
    const std::string filename =
      (output_filename + "_" + Utilities::int_to_string(cycle) + ".vtk");
    std::ofstream logfile(filename);
    GridOut       grid_out;
    grid_out.write_vtk(tria, logfile);
  }


  // @sect4{TriangulationOnCAD::run}

  // This is the main function. It should be self explanatory in its
  // briefness:
  void TriangulationOnCAD::run()
  {
    read_domain();

    const unsigned int n_cycles = 5;
    for (unsigned int cycle = 0; cycle < n_cycles; ++cycle)
      {
        refine_mesh();
        output_results(cycle + 1);
      }
  }
} // namespace Step54


// @sect3{The main() function}

// This is the main function of this program. It is in its basic structure
// like all previous tutorial programs, but runs the main class through the
// three possibilities of new vertex placement:
int main()
{
  try
    {
      using namespace Step54;

      const std::string in_mesh_filename = "input/initial_mesh_3d.vtk";
      const std::string cad_file_name    = "input/DTMB-5415_bulbous_bow.iges";

      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << "Testing projection in direction normal to CAD surface"
                << std::endl;
      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::string        out_mesh_filename = ("3d_mesh_normal_projection");
      TriangulationOnCAD tria_on_cad_norm(in_mesh_filename,
                                          cad_file_name,
                                          out_mesh_filename,
                                          TriangulationOnCAD::NormalProjection);
      tria_on_cad_norm.run();
      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << "Testing projection in y-axis direction" << std::endl;
      std::cout << "----------------------------------------------------------"
                << std::endl;
      out_mesh_filename = ("3d_mesh_directional_projection");
      TriangulationOnCAD tria_on_cad_dir(
        in_mesh_filename,
        cad_file_name,
        out_mesh_filename,
        TriangulationOnCAD::DirectionalProjection);
      tria_on_cad_dir.run();
      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;

      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << "Testing projection in direction normal to mesh elements"
                << std::endl;
      std::cout << "----------------------------------------------------------"
                << std::endl;
      out_mesh_filename = ("3d_mesh_normal_to_mesh_projection");
      TriangulationOnCAD tria_on_cad_norm_to_mesh(
        in_mesh_filename,
        cad_file_name,
        out_mesh_filename,
        TriangulationOnCAD::NormalToMeshProjection);
      tria_on_cad_norm_to_mesh.run();
      std::cout << "----------------------------------------------------------"
                << std::endl;
      std::cout << std::endl;
      std::cout << std::endl;
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
