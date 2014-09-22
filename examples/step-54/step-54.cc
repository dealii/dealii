/* ---------------------------------------------------------------------
 * $Id$
 *
 * Copyright (C) 2009 - 2013 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 *  Authors: Andrea Mola, Luca Heltai, 2014
 */


// @sect3{Include files}

// The program starts with including a bunch of include files that we will use
// in the various parts of the program. Most of them have been discussed in
// previous tutorials already:
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/convergence_table.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
#include <deal.II/base/parsed_function.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/solver_gmres.h>
#include <deal.II/lac/precondition.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>

#include <deal.II/opencascade/boundary_lib.h>
#include <deal.II/opencascade/utilities.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

// And here are a few C++ standard header files that we will need:
#include <cmath>
#include <iostream>
#include <fstream>
#include <string>

// The last part of this preamble is to import everything in the dealii
// namespace into the one into which everything in this program will go:
namespace Step54
{
  using namespace dealii;



  // @sect3{The TriangulationOnCAD class}

  // The structure of a boundary element method code is very similar to the
  // structure of a finite element code, and so the member functions of this
  // class are like those of most of the other tutorial programs. In
  // particular, by now you should be familiar with reading parameters from an
  // external file, and with the splitting of the different tasks into
  // different modules. The same applies to boundary element methods, and we
  // won't comment too much on them, except on the differences.

  class TriangulationOnCAD
  {
  public:
    TriangulationOnCAD(const unsigned int fe_degree = 1,
               const unsigned int mapping_degree = 1);

    
    ~TriangulationOnCAD();

    void run();

  private:

    void read_parameters (const std::string &filename);

    void read_domain();

    void refine_and_resize();

    void output_results(const unsigned int cycle);

    Triangulation<2, 3>   tria;
    FE_Q<2,3>             fe;
    DoFHandler<2,3>       dh;
    MappingQ<2,3>      mapping;

    std_cxx11::shared_ptr<Quadrature<2> > quadrature;

    SolverControl solver_control;

    unsigned int n_cycles;

    OpenCASCADE::ArclengthProjectionLineManifold<2,3> *line_projector;
    OpenCASCADE::NormalProjectionBoundary<2,3> *normal_projector;
    OpenCASCADE::DirectionalProjectionBoundary<2,3> *directional_projector;

  };


  // @sect4{TriangulationOnCAD::TriangulationOnCAD and TriangulationOnCAD::read_parameters}

  // The constructor initializes the various object in much the same way as
  // done in the finite element programs such as step-4 or step-6. The only
  // new ingredient here is the ParsedFunction object, which needs, at
  // construction time, the specification of the number of components.
  //
  // For the exact solution the number of vector components is one, and no
  // action is required since one is the default value for a ParsedFunction
  // object. The wind, however, requires dim components to be
  // specified. Notice that when declaring entries in a parameter file for the
  // expression of the Functions::ParsedFunction, we need to specify the
  // number of components explicitly, since the function
  // Functions::ParsedFunction::declare_parameters is static, and has no
  // knowledge of the number of components.

  TriangulationOnCAD::TriangulationOnCAD(const unsigned int fe_degree,
                              const unsigned int mapping_degree)
    :
    fe(fe_degree),
    dh(tria),
    mapping(mapping_degree, true)
  {}

  TriangulationOnCAD::~TriangulationOnCAD()
  {
  tria.set_manifold(1);
  tria.set_manifold(2);
  delete line_projector;
  delete normal_projector;
  delete directional_projector;
  }

  void TriangulationOnCAD::read_parameters (const std::string &filename)
  {
    deallog << std::endl << "Parsing parameter file " << filename << std::endl
            << "for a three dimensional geometry. " << std::endl;

    ParameterHandler prm;

    prm.declare_entry("Number of cycles", "4",
                      Patterns::Integer());

    prm.enter_subsection("Quadrature rules");
    {
      prm.declare_entry("Quadrature type", "gauss",
                        Patterns::Selection(QuadratureSelector<(2)>::get_quadrature_names()));
      prm.declare_entry("Quadrature order", "4", Patterns::Integer());
    }
    prm.leave_subsection();

    // After declaring all these parameters to the ParameterHandler object,
    // let's read an input file that will give the parameters their values. We
    // then proceed to extract these values from the ParameterHandler object:
    prm.read_input(filename);

    n_cycles = prm.get_integer("Number of cycles");

    prm.enter_subsection("Quadrature rules");
    {
      quadrature =
        std_cxx11::shared_ptr<Quadrature<2> >
        (new QuadratureSelector<2> (prm.get("Quadrature type"),
                                        prm.get_integer("Quadrature order")));
    }
    prm.leave_subsection();

  }


  // @sect4{TriangulationOnCAD::read_domain}


  // Some of the mesh formats supported in deal.II use by default three
  // dimensional points to describe meshes. These are the formats which are
  // compatible with the boundary element method capabilities of deal.II. In
  // particular we can use either UCD or GMSH formats. In both cases, we have
  // to be particularly careful with the orientation of the mesh, because,
  // unlike in the standard finite element case, no reordering or
  // compatibility check is performed here.  All meshes are considered as
  // oriented, because they are embedded in a higher dimensional space. (See
  // the documentation of the GridIn and of the Triangulation for further
  // details on orientation of cells in a triangulation.) In our case, the
  // normals to the mesh are external to both the circle in 2d or the sphere
  // in 3d.
  //
  // The other detail that is required for appropriate refinement of the
  // boundary element mesh, is an accurate description of the manifold that
  // the mesh is approximating. We already saw this several times for the
  // boundary of standard finite element meshes (for example in step-5 and
  // step-6), and here the principle and usage is the same, except that the
  // HyperBallBoundary class takes an additional template parameter that
  // specifies the embedding space dimension. The function object still has to
  // be static to live at least as long as the triangulation object to which
  // it is attached.


  void TriangulationOnCAD::read_domain()
  {


    std::ifstream in;

    in.open ("initial_mesh_3d.inp");

    GridIn<2,3> gi;
    gi.attach_triangulation (tria);
    gi.read_ucd (in);


    Triangulation<2,3>::active_cell_iterator cell = tria.begin_active();
    cell->set_manifold_id(1);

    for (unsigned int f=0; f<GeometryInfo<2>::faces_per_cell; ++f)
        cell->face(f)->set_manifold_id(2);

    TopoDS_Shape bow_surface = OpenCASCADE::read_IGES("DTMB-5415_bulbous_bow.iges",1e-3);

    std::vector<TopoDS_Compound> compounds;
    std::vector<TopoDS_CompSolid> compsolids;
    std::vector<TopoDS_Solid> solids;
    std::vector<TopoDS_Shell> shells;
    std::vector<TopoDS_Wire> wires;    

    OpenCASCADE::extract_compound_shapes(bow_surface,
                                         compounds,
                                         compsolids,
                                         solids,
                                         shells,
                                         wires);


    line_projector = new OpenCASCADE::ArclengthProjectionLineManifold<2,3>(wires[0]);
    normal_projector = new OpenCASCADE::NormalProjectionBoundary<2,3>(bow_surface);
    directional_projector = new OpenCASCADE::DirectionalProjectionBoundary<2,3>(bow_surface, Point<3>(0.0,1.0,0.0));

    tria.set_manifold(2, *line_projector);
    tria.set_manifold(1,*directional_projector);
  }


  // @sect4{TriangulationOnCAD::refine_and_resize}

  // This function globally refines the mesh, distributes degrees of freedom,
  // and resizes matrices and vectors.


  void TriangulationOnCAD::refine_and_resize()
  {
    tria.refine_global(1);

    dh.distribute_dofs(fe);

    const unsigned int n_dofs =  dh.n_dofs();

  }




  // @sect4{TriangulationOnCAD::output_results}

  // Outputting the results of our computations is a rather mechanical
  // tasks. All the components of this function have been discussed before.

  void TriangulationOnCAD::output_results(const unsigned int cycle)
  {

    std::string filename = ( Utilities::int_to_string(3) +
                             "d_mesh_" +
                             Utilities::int_to_string(cycle) +
                             ".inp" );
  std::ofstream logfile(filename.c_str());
  GridOut grid_out;
  grid_out.write_ucd(tria, logfile);


  }


  // @sect4{TriangulationOnCAD::run}

  // This is the main function. It should be self explanatory in its
  // briefness:

  void TriangulationOnCAD::run()
  {

    read_parameters("parameters.prm");

    read_domain();

    for (unsigned int cycle=0; cycle<n_cycles; ++cycle)
      {
        refine_and_resize();
        output_results(cycle);
      }

  }
}


// @sect3{The main() function}

// This is the main function of this program. It is exactly like all previous
// tutorial programs:
int main ()
{
  try
    {
      using namespace dealii;
      using namespace Step54;

      const unsigned int degree = 1;
      const unsigned int mapping_degree = 1;

      deallog.depth_console (3);
      TriangulationOnCAD triangulation_on_cad(degree, mapping_degree);
      triangulation_on_cad.run();

    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
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
      std::cerr << std::endl << std::endl
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

