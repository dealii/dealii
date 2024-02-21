// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// tests DataOut with HDF5 and zero-dimensional output (i.e. points)

#include <deal.II/base/data_out_base.h>

#include <string>
#include <vector>

#include "../tests.h"

double cell_coordinates[3][8] = {{0, 1, 0, 1, 0, 1, 0, 1},
                                 {0, 0, 1, 1, 0, 0, 1, 1},
                                 {0, 0, 0, 0, 1, 1, 1, 1}};


// This function is a copy from tests/base/patches.h, included here
// to not introduce dependencies between different test targets
template <int dim, int spacedim>
void
create_patches(std::vector<DataOutBase::Patch<dim, spacedim>> &patches)
{
  for (unsigned int p = 0; p < patches.size(); ++p)
    {
      DataOutBase::Patch<dim, spacedim> &patch = patches[p];

      const unsigned int nsub  = p + 1;
      const unsigned int nsubp = nsub + 1;

      patch.n_subdivisions = nsub;
      patch.reference_cell = ReferenceCells::get_hypercube<dim>();
      for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
        for (unsigned int d = 0; d < spacedim; ++d)
          patch.vertices[v](d) =
            p + cell_coordinates[d][v] + ((d >= dim) ? v : 0);

      unsigned int n1 = (dim > 0) ? nsubp : 1;
      unsigned int n2 = (dim > 1) ? nsubp : 1;
      unsigned int n3 = (dim > 2) ? nsubp : 1;
      unsigned int n4 = (dim > 3) ? nsubp : 1;
      patch.data.reinit(5, n1 * n2 * n3 * n4);

      for (unsigned int i4 = 0; i4 < n4; ++i4)
        for (unsigned int i3 = 0; i3 < n3; ++i3)
          for (unsigned int i2 = 0; i2 < n2; ++i2)
            for (unsigned int i1 = 0; i1 < n1; ++i1)
              {
                const unsigned int i =
                  i1 + nsubp * (i2 + nsubp * (i3 + nsubp * i4));
                const float x1 = 1. * i1 / nsub;
                const float x2 = 1. * i2 / nsub;
                const float x3 = 1. * i3 / nsub;
                const float x4 = 1. * i4 / nsub;

                patch.data(0, i) = p + x1;
                patch.data(1, i) = p + x2;
                patch.data(2, i) = p + x3;
                patch.data(3, i) = p + x4;
                patch.data(4, i) = i;
              }
      patch.patch_index = p;
    }
}



template <int dim, int spacedim>
void
check()
{
  const unsigned int np = 4;

  std::vector<DataOutBase::Patch<dim, spacedim>> patches(np);

  create_patches(patches);

  std::vector<std::string> names(5);
  names[0] = "x1";
  names[1] = "x2";
  names[2] = "x3";
  names[3] = "x4";
  names[4] = "i";
  std::vector<std::tuple<unsigned int, unsigned int, std::string>> vectors;

  DataOutBase::DataOutFilter data_filter(
    DataOutBase::DataOutFilterFlags(false, false));

  DataOutBase::write_filtered_data(patches, names, vectors, data_filter);

  std::string output_basename = std::to_string(dim) + std::to_string(spacedim);

  DataOutBase::write_hdf5_parallel(patches,
                                   data_filter,
                                   output_basename + ".h5",
                                   MPI_COMM_SELF);

  const double current_time = 0.0;
  XDMFEntry    entry(output_basename + ".h5",
                  output_basename + ".h5",
                  current_time,
                  data_filter.n_nodes(),
                  data_filter.n_cells(),
                  dim,
                  spacedim);
  unsigned int n_data_sets = data_filter.n_data_sets();

  // The vector names generated here must match those generated in the HDF5 file
  for (unsigned int i = 0; i < n_data_sets; ++i)
    {
      entry.add_attribute(data_filter.get_data_set_name(i),
                          data_filter.get_data_set_dim(i));
    }

  std::ofstream xdmf_file(output_basename + ".xdmf");

  xdmf_file << "<?xml version=\"1.0\" ?>\n";
  xdmf_file << "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n";
  xdmf_file << "<Xdmf Version=\"2.0\">\n";
  xdmf_file << "  <Domain>\n";
  xdmf_file
    << "    <Grid Name=\"CellTime\" GridType=\"Collection\" CollectionType=\"Temporal\">\n";

  // Write out the entry
  xdmf_file << entry.get_xdmf_content(3);

  xdmf_file << "    </Grid>\n";
  xdmf_file << "  </Domain>\n";
  xdmf_file << "</Xdmf>\n";

  xdmf_file.close();

  // Sadly hdf5 is binary and we can not use hd5dump because it might
  // not be in the path. At least we can look at the xdmf
  // and make sure that the h5 file is created:
  deallog << "\n==============================\n"
          << output_basename + ".xdmf"
          << "\n==============================" << std::endl;

  cat_file((output_basename + ".xdmf").c_str());
  std::ifstream f(output_basename + ".h5");
  AssertThrow(f.good(), ExcIO());

  deallog << std::flush;
}



int
main(int argc, char *argv[])
{
  initlog();

  try
    {
      check<0, 2>();
      check<1, 2>();
      check<2, 2>();
      check<0, 3>();
      check<1, 3>();
      check<2, 3>();
      check<3, 3>();

      return 0;
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
