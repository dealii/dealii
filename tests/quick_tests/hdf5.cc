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
