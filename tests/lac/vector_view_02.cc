// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// check VectorView::checkReinit1

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_view.h>
#include <cmath>
#include <fstream>
#include <iomanip>

template<typename number, typename size_type>
  void
  checkReinit1(const size_type N, const bool fast = false)
  {
    deallog << "Reinit with const size and fast" << std::endl;

    deallog
        << "Creating Vector<number> of size N+10 and filling with values 1 to N+10"
        << std::endl;

    Vector < number > V(N + 10);
    for (unsigned int i = 0; i < V.size(); i++)
      V(i) = i + 1;

    deallog
        << "Creating VectorView<number> of size N+10 pointing to Vector<number>"
        << std::endl;
    VectorView<number> VV(V.size(), V.begin());

    deallog << "Printing Vector<number>" << std::endl;
    for (unsigned int i = 0; i < V.size(); ++i)
      deallog << V(i) << '\t';
    deallog << std::endl;

    deallog << "Printing VectorView<number> pointing to Vector<number>"
        << std::endl;
    for (unsigned int i = 0; i < VV.size(); ++i)
      deallog << VV(i) << '\t';
    deallog << std::endl;

    deallog << "Reinit VectorView<number> to size N from N+10 with fast="
        << fast << std::endl;
    VV.reinit(N, fast);

    deallog << "Printing Vector<number>" << std::endl;
    for (unsigned int i = 0; i < V.size(); ++i)
      deallog << V(i) << '\t';
    deallog << std::endl;

    deallog << "Printing VectorView<number> pointing to Vector<number>"
        << std::endl;
    for (unsigned int i = 0; i < VV.size(); ++i)
      deallog << VV(i) << '\t';
    deallog << std::endl;
  }

int
main()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  checkReinit1<double, int>(10, false);
  checkReinit1<double, int>(10, true);
}

