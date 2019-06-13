/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2019 by the deal.II authors
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

 */

// @sect3{Include files}

#include <deal.II/base/mpi.h>


namespace step67
{
  using namespace dealii;
}

// @sect3{The main function}

// Finally, the main function. There isn't much to do here, only to call the
// two subfunctions, which produce the two grids.
int main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_init_finalize(argc, argv, 1);
}
