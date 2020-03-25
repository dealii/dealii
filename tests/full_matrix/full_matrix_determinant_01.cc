// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

  // 13x13
  {
    std::vector<number> entries = {
      0.302024,  0.80384,   0.693137,  0.748042,   0.404792,  0.910457,
      0.852546,  0.307047,  0.0610559, 0.446627,   0.213923,  0.884411,
      0.965113,  0.835556,  0.136569,  0.76459,    0.122895,  0.353153,
      0.725082,  0.0164721, 0.732707,  0.708028,   0.326752,  0.1112,
      0.344864,  0.964341,  0.411462,  0.158079,   0.933206,  0.629767,
      0.177095,  0.543207,  0.651351,  0.330828,   0.95767,   0.0437675,
      0.406998,  0.0614156, 0.229926,  0.210028,   0.123443,  0.653746,
      0.977856,  0.70126,   0.98382,   0.420717,   0.142192,  0.0750294,
      0.447152,  0.684233,  0.777946,  0.866352,   0.46331,   0.147653,
      0.748947,  0.581605,  0.369571,  0.00815871, 0.587899,  0.328085,
      0.433468,  0.05921,   0.226504,  0.551432,   0.473922,  0.961679,
      0.564591,  0.764575,  0.748242,  0.358666,   0.709481,  0.858687,
      0.644768,  0.230691,  0.309604,  0.925409,   0.732361,  0.199902,
      0.164676,  0.818629,  0.16853,   0.581234,   0.736271,  0.591147,
      0.32312,   0.818551,  0.541567,  0.724956,   0.984736,  0.357938,
      0.446607,  0.541354,  0.609206,  0.576774,   0.910311,  0.614716,
      0.60952,   0.979587,  0.838219,  0.568433,   0.881019,  0.916143,
      0.34111,   0.0607764, 0.88064,   0.0948878,  0.233047,  0.394719,
      0.971708,  0.757256,  0.146067,  0.840419,   0.926348,  0.184664,
      0.0336569, 0.440521,  0.825535,  0.730938,   0.772993,  0.92274,
      0.106755,  0.138245,  0.566582,  0.96567,    0.949705,  0.306954,
      0.538506,  0.188428,  0.282917,  0.617393,   0.638839,  0.901067,
      0.188109,  0.599985,  0.860141,  0.320569,   0.852939,  0.049655,
      0.34805,   0.0983851, 0.906911,  0.136246,   0.0789367, 0.548858,
      0.157769,  0.444907,  0.476137,  0.27633,    0.480897,  0.13293,
      0.822298,  0.467239,  0.0368684, 0.47657,    0.401819,  0.0374826,
      0.962836,  0.855119,  0.702701,  0.809599,   0.346854,  0.658641,
      0.349523,  0.900854,  0.497165,  0.237687,   0.0472559, 0.491044,
      0.486307};

    FullMatrix<number> A(13, 13, entries.data());

    const number soln_det_A = -0.0136107;
    const number rel_tol    = 1e-5;
    if (std::abs(A.determinant() - soln_det_A) < std::abs(rel_tol * soln_det_A))
      deallog << "OK" << std::endl;
    else
      {
        std::ostringstream ss;
        ss << A.determinant();
        deallog << "Incorrect 13x13 determinant. Computed " << ss.str()
                << std::endl;
      }
  }
}
