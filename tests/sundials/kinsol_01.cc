//-----------------------------------------------------------
//
//    Copyright (C) 2017 - 2020 by the deal.II authors
//
//    This file is part of the deal.II library.
//
//    The deal.II library is free software; you can use it, redistribute
//    it, and/or modify it under the terms of the GNU Lesser General
//    Public License as published by the Free Software Foundation; either
//    version 2.1 of the License, or (at your option) any later version.
//    The full text of the license can be found in the file LICENSE.md at
//    the top level directory of deal.II.
//
//-----------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/sundials/kinsol.h>

#include "../tests.h"

// provide only residual function, use internal solver.

/**
 * Solve the non linear problem
 *
 * F(u) = 0 , where f_i(u) = u_i^2 - i^2,  0 <= i < N
 *
 */
int
main(int argc, char **argv)
{
  initlog();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, numbers::invalid_unsigned_int);

  typedef Vector<double> VectorType;

  SUNDIALS::KINSOL<VectorType>::AdditionalData data;
  ParameterHandler                             prm;
  data.add_parameters(prm);

  if (false)
    {
      std::ofstream ofile(SOURCE_DIR "/kinsol_01.prm");
      prm.print_parameters(ofile, ParameterHandler::ShortText);
      ofile.close();
    }

  std::ifstream ifile(SOURCE_DIR "/kinsol_01.prm");
  prm.parse_input(ifile);

  // Size of the problem
  unsigned int N = 10;

  SUNDIALS::KINSOL<VectorType> kinsol(data);

  kinsol.reinit_vector = [N](VectorType &v) { v.reinit(N); };

  kinsol.residual = [](const VectorType &u, VectorType &F) -> int {
    for (unsigned int i = 0; i < u.size(); ++i)
      F[i] = u[i] * u[i] - (i + 1) * (i + 1);
    return 0;
  };


  kinsol.iteration_function = [](const VectorType &u, VectorType &F) -> int {
    for (unsigned int i = 0; i < u.size(); ++i)
      F[i] = u[i] * u[i] - i * i - u[i];
    return 0;
  };

  VectorType v(N);
  v          = 1.0;
  auto niter = kinsol.solve(v);
  v.print(deallog.get_file_stream());
  deallog << "Converged in " << niter << " iterations." << std::endl;
}
