// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II authors
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

// Test linear_operator constructor for full and reduced interface:

#include <deal.II/lac/linear_operator.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_memory.templates.h>

#include "../tests.h"

class MyMatrix1
{
public:
  size_t
  m() const
  {
    return 1;
  };
  size_t
  n() const
  {
    return 1;
  };

  void
  vmult(Vector<double> &, const Vector<double> &) const
  {
    deallog << "MyMatrix1::vmult" << std::endl;
  }

  void
  Tvmult(Vector<double> &, const Vector<double> &) const
  {
    deallog << "MyMatrix1::Tvmult" << std::endl;
  }
};

class MyMatrix2
{
public:
  size_t
  m() const
  {
    return 1;
  };
  size_t
  n() const
  {
    return 1;
  };

  void
  vmult(Vector<double> &, const Vector<double> &) const
  {
    deallog << "MyMatrix2::vmult" << std::endl;
  }

  void
  Tvmult(Vector<double> &, const Vector<double> &) const
  {
    deallog << "MyMatrix2::Tvmult" << std::endl;
  }

  void
  vmult_add(Vector<double> &, const Vector<double> &) const
  {
    deallog << "MyMatrix2::vmult_add" << std::endl;
  }

  void
  Tvmult_add(Vector<double> &, const Vector<double> &) const
  {
    deallog << "MyMatrix2::Tvmult_add" << std::endl;
  }
};

class MyMatrix3
{
public:
  size_t
  m() const
  {
    return 1;
  };
  size_t
  n() const
  {
    return 1;
  };

  template <typename OutVector, typename InVector>
  void
  vmult(OutVector &, const InVector &) const
  {
    deallog << "MyMatrix3::vmult" << std::endl;
  }

  template <typename OutVector, typename InVector>
  void
  Tvmult(OutVector &, const InVector &) const
  {
    deallog << "MyMatrix3::Tvmult" << std::endl;
  }

  template <typename OutVector, typename InVector>
  void
  vmult_add(OutVector &, const InVector &) const
  {
    deallog << "MyMatrix3::vmult_add" << std::endl;
  }

  template <typename OutVector, typename InVector>
  void
  Tvmult_add(OutVector &, const InVector &) const
  {
    deallog << "MyMatrix3::Tvmult_add" << std::endl;
  }
};

class MyMatrix4
{
public:
  size_t
  m() const
  {
    return 1;
  };
  size_t
  n() const
  {
    return 1;
  };

  template <typename OutVector, typename InVector>
  void
  vmult(OutVector &, const InVector &, bool = true) const
  {
    deallog << "MyMatrix4::vmult" << std::endl;
  }

  template <typename OutVector, typename InVector>
  void
  Tvmult(OutVector &, const InVector &, bool = true) const
  {
    deallog << "MyMatrix4::Tvmult" << std::endl;
  }

  template <typename OutVector, typename InVector>
  void
  vmult_add(OutVector &, const InVector &, bool = true) const
  {
    deallog << "MyMatrix4::vmult_add" << std::endl;
  }

  template <typename OutVector, typename InVector>
  void
  Tvmult_add(OutVector &, const InVector &, bool = true) const
  {
    deallog << "MyMatrix4::Tvmult_add" << std::endl;
  }
};

class MyMatrix5
{
public:
  size_t
  m() const
  {
    return 1;
  };
  size_t
  n() const
  {
    return 1;
  };

  template <typename number>
  void
  vmult(Vector<number> &, const Vector<number> &, bool = true) const
  {
    deallog << "MyMatrix5::vmult" << std::endl;
  }

  template <typename number>
  void
  Tvmult(Vector<number> &, const Vector<number> &, bool = true) const
  {
    deallog << "MyMatrix5::Tvmult" << std::endl;
  }

  template <typename number>
  void
  vmult_add(Vector<number> &, const Vector<number> &, bool = true) const
  {
    deallog << "MyMatrix5::vmult_add" << std::endl;
  }

  template <typename number>
  void
  Tvmult_add(Vector<number> &, const Vector<number> &, bool = true) const
  {
    deallog << "MyMatrix5::Tvmult_add" << std::endl;
  }
};

class MyMatrix6
{
public:
  size_t
  m() const
  {
    return 1;
  };
  size_t
  n() const
  {
    return 1;
  };

  template <typename number, typename number2>
  void
  vmult(Vector<number> &, const Vector<number2> &, bool = true) const
  {
    deallog << "MyMatrix6::vmult" << std::endl;
  }

  template <typename number, typename number2>
  void
  Tvmult(Vector<number> &, const Vector<number2> &, bool = true) const
  {
    deallog << "MyMatrix6::Tvmult" << std::endl;
  }

  template <typename number, typename number2>
  void
  vmult_add(Vector<number> &, const Vector<number2> &, bool = true) const
  {
    deallog << "MyMatrix6::vmult_add" << std::endl;
  }
};


int
main()
{
  initlog();
  deallog << std::setprecision(10);

  Vector<double> u(1);
  Vector<double> v(1);

  MyMatrix1 m1;
  MyMatrix2 m2;
  MyMatrix3 m3;
  MyMatrix4 m4;
  MyMatrix5 m5;
  MyMatrix6 m6;

  linear_operator(m1).vmult(v, u);
  linear_operator(m1).Tvmult(v, u);
  linear_operator(m1).vmult_add(v, u);
  linear_operator(m1).Tvmult_add(v, u);

  linear_operator(m2).vmult(v, u);
  linear_operator(m2).Tvmult(v, u);
  linear_operator(m2).vmult_add(v, u);
  linear_operator(m2).Tvmult_add(v, u);

  linear_operator(m3).vmult(v, u);
  linear_operator(m3).Tvmult(v, u);
  linear_operator(m3).vmult_add(v, u);
  linear_operator(m3).Tvmult_add(v, u);

  linear_operator(m4).vmult(v, u);
  linear_operator(m4).Tvmult(v, u);
  linear_operator(m4).vmult_add(v, u);
  linear_operator(m4).Tvmult_add(v, u);

  linear_operator(m5).vmult(v, u);
  linear_operator(m5).Tvmult(v, u);
  linear_operator(m5).vmult_add(v, u);
  linear_operator(m5).Tvmult_add(v, u);

  linear_operator(m6).vmult(v, u);
  linear_operator(m6).Tvmult(v, u);
  linear_operator(m6).vmult_add(v, u);
  linear_operator(m6).Tvmult_add(v, u);
}
