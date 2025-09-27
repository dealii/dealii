// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Tests FE_Nedelec finite elements of orders 0...4 on a simple 2D mesh
// with a shared face orientations 0 and 1 (the are only two possible face
// orientations in 2D). See the notes in the header file for more detail.

#include "fe_nedelec_face_edge_orientation.h"

int
main()
{
  unsigned char combined_face_orientation = 0;

  initlog();

  run<2, 0>(combined_face_orientation);
  run<2, 1>(combined_face_orientation);
  run<2, 2>(combined_face_orientation);
  run<2, 3>(combined_face_orientation);
  run<2, 4>(combined_face_orientation);

  combined_face_orientation = 1;

  run<2, 0>(combined_face_orientation);
  run<2, 1>(combined_face_orientation);
  run<2, 2>(combined_face_orientation);
  run<2, 3>(combined_face_orientation);
  run<2, 4>(combined_face_orientation);
}
