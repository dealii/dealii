// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Test HDF5. Copy-paste from
// https://support.hdfgroup.org/HDF5/Tutor/crtfile.html

#include <deal.II/base/exceptions.h>

#include <hdf5.h>

#include <cstdio>
int
main()
{
  hid_t  file_id;
  herr_t status;

  /* Create a new file using default properties. */
  file_id = H5Fcreate("file.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

  /* Terminate access to the file. */
  status = H5Fclose(file_id);
  AssertThrow(status >= 0, dealii::ExcInternalError());
  std::remove("file.h5");
}
