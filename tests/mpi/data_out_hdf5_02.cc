// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2024 by the deal.II authors
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

#include "../data_out/patches.h"

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
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vectors;

  DataOutBase::DataOutFilter data_filter(
    DataOutBase::DataOutFilterFlags(false, false));

  DataOutBase::write_filtered_data(patches, names, vectors, data_filter);

  std::string output_basename = std::to_string(dim) + std::to_string(spacedim);

  DataOutBase::Hdf5Flags hdf5Flags;
  hdf5Flags.compression_level = DataOutBase::CompressionLevel::no_compression;

  DataOutBase::write_hdf5_parallel(
    patches, data_filter, hdf5Flags, output_basename + ".h5", MPI_COMM_WORLD);

  const double current_time = 0.0;
  XDMFEntry    entry(output_basename + ".h5",
                  output_basename + ".h5",
                  current_time,
                  data_filter.n_nodes(),
                  data_filter.n_cells(),
                  dim,
                  spacedim,
                  ReferenceCell());
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

  if (0 == Utilities::MPI::this_mpi_process(MPI_COMM_WORLD))
    {
      cat_file((output_basename + ".xdmf").c_str());
      std::ifstream f(output_basename + ".h5");
      AssertThrow(f.good(), ExcIO());
    }

  deallog << std::flush;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

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
