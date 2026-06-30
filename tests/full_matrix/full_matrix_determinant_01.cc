// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------


// Check that the determinant for a FullMatrix of
// size greater than 3 is correct.

#include <sstream>

#include "../tests.h"

#include "full_matrix_common.h"



template <typename number>
void
check()
{
  // 4x4
  {
    FullMatrix<number> A(4, 4);

    // Wolfram Alpha: determinant {{3,0,2,-1},{1,2,0,-2},{4,0,6,-3},{5,0,2,0}}
    // det(A) = 20
    A(0, 0) = 3;
    A(0, 1) = 0;
    A(0, 2) = 2;
    A(0, 3) = -1;
    A(1, 0) = 1;
    A(1, 1) = 2;
    A(1, 2) = 0;
    A(1, 3) = -2;
    A(2, 0) = 4;
    A(2, 1) = 0;
    A(2, 2) = 6;
    A(2, 3) = -3;
    A(3, 0) = 5;
    A(3, 1) = 0;
    A(3, 2) = 2;
    A(3, 3) = 0;

    const number soln_det_A = 20.0;
    const number rel_tol    = 1e-6;
    if (std::abs(A.determinant() - soln_det_A) < std::abs(rel_tol * soln_det_A))
      deallog << "OK" << std::endl;
    else
      {
        std::ostringstream ss;
        ss << A.determinant();
        deallog << "Incorrect 4x4 determinant. Computed " << ss.str()
                << std::endl;
      }
  }

  // 5x5
  {
    // https://math.stackexchange.com/questions/1955784/how-to-find-the-determinant-of-a-5x5-matrix
    // Wolfram Alpha: determinant
    // {{0,6,-2,-1,5},{0,0,0,-9,-7},{0,15,35,0,0},{0,-1,-11,-2,1},{-2,-2,3,0,2}}
    FullMatrix<number> A(5, 5);

    // det(A) = 2480
    A(0, 0) = 0;
    A(0, 1) = 6;
    A(0, 2) = -2;
    A(0, 3) = -1;
    A(0, 4) = 5;
    A(1, 0) = 0;
    A(1, 1) = 0;
    A(1, 2) = 0;
    A(1, 3) = -9;
    A(1, 4) = -7;
    A(2, 0) = 0;
    A(2, 1) = 15;
    A(2, 2) = 35;
    A(2, 3) = 0;
    A(2, 4) = 0;
    A(3, 0) = 0;
    A(3, 1) = -1;
    A(3, 2) = -11;
    A(3, 3) = -2;
    A(3, 4) = 1;
    A(4, 0) = -2;
    A(4, 1) = -2;
    A(4, 2) = 3;
    A(4, 3) = 0;
    A(4, 4) = 2;

    const number soln_det_A = 2480.0;
    const number rel_tol    = 1e-6;
    if (std::abs(A.determinant() - soln_det_A) < std::abs(rel_tol * soln_det_A))
      deallog << "OK" << std::endl;
    else
      {
        std::ostringstream ss;
        ss << A.determinant();
        deallog << "Incorrect 5x5 determinant. Computed " << ss.str()
                << std::endl;
      }
  }

  // 8x8
  {
    // A matrix generation through Wolfram Alpha:
    // RandomReal[{0,1},{n_rows,n_cols}]
    std::vector<number> entries = {
      0.275566,  0.384198, 0.0601487, 0.631826, 0.657166, 0.567488, 0.815103,
      0.522585,  0.721087, 0.869243,  0.256089, 0.276967, 0.346444, 0.850311,
      0.0443296, 0.417804, 0.673882,  0.854942, 0.686326, 0.23442,  0.157204,
      0.109822,  0.660962, 0.49385,   0.242054, 0.678298, 0.152465, 0.0110467,
      0.915869,  0.292891, 0.519057,  0.292487, 0.226375, 0.306915, 0.421131,
      0.0935529, 0.47106,  0.314059,  0.901549, 0.770425, 0.827285, 0.230283,
      0.408039,  0.662393, 0.103886,  0.169664, 0.247615, 0.960181, 0.175275,
      0.915451,  0.80152,  0.665238,  0.623794, 0.440641, 0.582839, 0.629464,
      0.0614058, 0.94639,  0.198588,  0.85945,  0.169552, 0.325892, 0.882038,
      0.0442052};
    FullMatrix<number> A(8, 8, entries.data());

    const number soln_det_A = 0.0626735;
    const number rel_tol    = 1e-6;
    if (std::abs(A.determinant() - soln_det_A) < std::abs(rel_tol * soln_det_A))
      deallog << "OK" << std::endl;
    else
      {
        std::ostringstream ss;
        ss << A.determinant();
        deallog << "Incorrect 8x8 determinant. Computed " << ss.str()
                << std::endl;
      }
  }

  // 7x7
  {
    std::vector<number> entries = {
      0.44029686, 0.4804838,  0.48551407, 0.7560161,  0.29350683, 0.9115595,
      0.02037634, 0.6845042,  0.30150527, 0.5702821,  0.23861961, 0.16337413,
      0.13380605, 0.92997324, 0.21977824, 0.8540596,  0.5130981,  0.165114,
      0.57020974, 0.9423338,  0.7913783,  0.3021855,  0.34198377, 0.49957728,
      0.5278706,  0.6399972,  0.08217729, 0.4763273,  0.9345214,  0.56911814,
      0.58146656, 0.94194585, 0.14227125, 0.2101821,  0.5399678,  0.8991025,
      0.08847044, 0.31052095, 0.9138506,  0.08069261, 0.94184846, 0.9575666,
      0.595712,   0.9861003,  0.754456,   0.17039633, 0.5623334,  0.80614275,
      0.7259835};

    FullMatrix<number> A(7, 7, entries.data());

    const number soln_det_A = 0.036075283;
    const number rel_tol    = 1e-5;
    if (std::abs(A.determinant() - soln_det_A) < std::abs(rel_tol * soln_det_A))
      deallog << "OK" << std::endl;
    else
      {
        std::ostringstream ss;
        ss << A.determinant();
        deallog << "Incorrect 7x7 determinant. Computed " << ss.str()
                << std::endl;
      }
  }
}
