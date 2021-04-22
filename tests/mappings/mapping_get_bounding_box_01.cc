// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2020 by the deal.II authors
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

// Check that bounding boxes are created correctly with high order mappings

#include <deal.II/base/bounding_box_data_out.h>
#include <deal.II/base/utilities.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools_cache.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <boost/core/demangle.hpp>

#include "../tests.h"


template <int dim, int spacedim>
void
test_bounding_boxes(const Mapping<dim, spacedim> &mapping,
                    const unsigned int            degree)
{
  deallog << "Testing " << boost::core::demangle(typeid(mapping).name()) << "("
          << degree << ")" << std::endl;

  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_ball(triangulation);
  GridTools::Cache<dim, spacedim> cache(triangulation, mapping);

  const auto boxes = cache.get_cell_bounding_boxes_rtree();

  {
    std::string fname = "boxes_" + std::to_string(dim) + "_" +
                        std::to_string(spacedim) + "_" +
                        std::to_string(degree) + ".vtk";
    std::ofstream           ofile(fname);
    BoundingBoxDataOut<dim> data_out;
    DataOutBase::VtkFlags   flags;
    flags.print_date_and_time = false;
    data_out.set_flags(flags);
    data_out.build_patches(boxes);
    data_out.write_vtk(ofile);
    cat_file(fname.c_str());
  }
}


int
main()
{
  initlog();
  {
    unsigned int degree = 2;
    MappingQ<2>  mapping(degree);
    test_bounding_boxes(mapping, degree);
  }
  {
    unsigned int       degree = 2;
    MappingQGeneric<2> mapping(degree);
    test_bounding_boxes(mapping, degree);
  }
}
