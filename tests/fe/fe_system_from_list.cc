// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// verify that the idea discussed in the multiple-element constructor FESystem
// really works

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <string>
#include <vector>

#include "../tests.h"

#define PRECISION 2

template <int dim>
struct MyFESystem
{
  MyFESystem(const std::vector<const FiniteElement<dim> *> &fes,
             const std::vector<unsigned int>               &multiplicities)
  {
    deallog << "Constructing FESystem from list." << std::endl;
  }
};



template <int dim>
class MySimulator
{
public:
  MySimulator(const unsigned int polynomial_degree);

private:
  MyFESystem<dim> fe;

  struct VectorElementDestroyer
  {
    const std::vector<const FiniteElement<dim> *> data;
    VectorElementDestroyer(
      const std::vector<const FiniteElement<dim> *> &pointers);
    ~VectorElementDestroyer(); // destructor to delete the pointers
    const std::vector<const FiniteElement<dim> *> &
    get_data() const;
  };

  static std::vector<const FiniteElement<dim> *>
  create_fe_list(const unsigned int polynomial_degree);

  static std::vector<unsigned int>
  create_fe_multiplicities();
};


template <int dim>
std::vector<const FiniteElement<dim> *>
MySimulator<dim>::create_fe_list(const unsigned int polynomial_degree)
{
  std::vector<const FiniteElement<dim> *> fe_list;
  fe_list.push_back(new FE_Q<dim>(1));
  fe_list.push_back(new FE_Q<dim>(2));
  fe_list.push_back(new FE_Q<dim>(3));
  fe_list.push_back(new FE_Q<dim>(4));
  deallog << "Created list of elements" << std::endl;
  return fe_list;
}

template <int dim>
std::vector<unsigned int>
MySimulator<dim>::create_fe_multiplicities()
{
  std::vector<unsigned int> multiplicities;
  multiplicities.push_back(1);
  multiplicities.push_back(2);
  multiplicities.push_back(3);
  multiplicities.push_back(4);
  return multiplicities;
}

template <int dim>
MySimulator<dim>::VectorElementDestroyer::VectorElementDestroyer(
  const std::vector<const FiniteElement<dim> *> &pointers)
  : data(pointers)
{}

template <int dim>
MySimulator<dim>::VectorElementDestroyer::~VectorElementDestroyer()
{
  for (unsigned int i = 0; i < data.size(); ++i)
    {
      deallog << "Destroying element " << data[i]->get_name() << std::endl;
      delete data[i];
    }
}

template <int dim>
const std::vector<const FiniteElement<dim> *> &
MySimulator<dim>::VectorElementDestroyer::get_data() const
{
  return data;
}


template <int dim>
MySimulator<dim>::MySimulator(const unsigned int polynomial_degree)
  : fe(VectorElementDestroyer(create_fe_list(polynomial_degree)).get_data(),
       create_fe_multiplicities())
{}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);

  MySimulator<1>(1);
  MySimulator<2>(1);
  MySimulator<3>(1);

  return 0;
}
