// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Checks that project_boundary_values_curl_conforming
// works for FESystems with FE_Nedelec and FE_Q mixed


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


const unsigned int dim = 3;

void
apply_boundary_values(DoFHandler<dim>           &dof_handler,
                      AffineConstraints<double> &constraints,
                      unsigned int               n_comps,
                      unsigned int               start_comp,
                      Mapping<dim>              &mapping,
                      Vector<double>            &dst)
{
  constraints.clear();
  for (unsigned int i = 0; i < 6; ++i)
    VectorTools::project_boundary_values_curl_conforming_l2(
      dof_handler,
      start_comp,
      Functions::ConstantFunction<dim>(1, n_comps),
      i,
      constraints,
      mapping);
  constraints.close();

  constraints.distribute(dst);
}



bool
test_boundary_values(DoFHandler<dim>    &dof_handler,
                     Mapping<dim>       &mapping,
                     FiniteElement<dim> &fe,
                     unsigned int        n_comps,
                     unsigned int        start_comp,
                     Vector<double>     &vec)
{
  // Initialize
  QGaussLobatto<dim - 1> quadrature(3);
  FEFaceValues<dim>      fe_values(mapping,
                              fe,
                              quadrature,
                              update_values | update_normal_vectors |
                                update_quadrature_points | update_JxW_values);

  unsigned                    num_cells = 0;
  std::vector<Vector<double>> local_averages(
    6, Vector<double>(dof_handler.n_dofs()));
  const FEValuesExtractors::Vector nedelec(start_comp);

  // For testing we take only one cell!
  Assert(dof_handler.n_dofs() ==
           dof_handler.begin_active()->get_fe().dofs_per_cell,
         ExcInternalError());

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // For testing we take only one cell! Otherwise, local to global is
      // missing
      assert(num_cells == 0);
      for (unsigned int face = 0;
           face < dealii::GeometryInfo<dim>::faces_per_cell;
           ++face)
        {
          if (cell->face(face)->at_boundary())
            {
              unsigned int   boundary_number = cell->face(face)->boundary_id();
              Tensor<1, dim> expected_value;
              expected_value[0] = 1;
              expected_value[1] = 1;
              expected_value[2] = 1;
              assert(boundary_number < 6);
              fe_values.reinit(cell, face);

              std::vector<Vector<double>> local_values(
                fe_values.n_quadrature_points, Vector<double>(n_comps));
              fe_values.get_function_values(vec, local_values);

              for (unsigned int q_point = 0;
                   q_point < fe_values.n_quadrature_points;
                   q_point++)
                {
                  // Compare n x v values
                  Tensor<1, dim> nedelec_values;
                  nedelec_values[0]     = local_values[q_point](start_comp);
                  nedelec_values[1]     = local_values[q_point](start_comp + 1);
                  nedelec_values[2]     = local_values[q_point](start_comp + 2);
                  Tensor<1, dim> normal = fe_values.normal_vector(q_point);

                  for (unsigned int i = 0; i < cell->get_fe().dofs_per_cell;
                       i++)
                    {
                      local_averages[boundary_number](i) +=
                        ((cross_product_3d(normal, nedelec_values) -
                          cross_product_3d(normal, expected_value)) *
                         cross_product_3d(
                           normal, fe_values[nedelec].value(i, q_point))) *
                        fe_values.JxW(q_point);
                    }
                }
            }
        }
      num_cells++;
    }
  // Test if ok
  bool ret = true;
  for (unsigned int bc = 0; bc < 6; ++bc)
    {
      for (unsigned int basis = 0; basis < dof_handler.n_dofs(); ++basis)
        {
          if (fabs(local_averages[bc](basis)) > 1.e-10)
            {
              deallog << "Wrong values on boundary " << bc
                      << " in basis function " << basis << " of size "
                      << fabs(local_averages[bc](basis)) << std::endl;
              ret = false;
            }
        }
    }
  return ret;
}


int
main(int /*argc*/, char ** /*argv*/)
{
  initlog();

  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation, 0, 1, true);
  MappingQ<dim>             mapping(1);
  AffineConstraints<double> constraints;
  Vector<double>            test_vec;
  Vector<double>            test_vec2;

  {
    // Test 1: Only FE_Nedelec
    deallog << "Testing FE_Nedelec(0): ";
    FE_Nedelec<dim> fe(0);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);
    test_vec2.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 3, 0, mapping, test_vec2);
    if (test_boundary_values(dof_handler, mapping, fe, 3, 0, test_vec2))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 2: FESystem(FE_Nedelec(0))
    deallog << "Testing FESystem(FE_Nedelec(0)): ";
    FESystem<dim>   fe_system(FE_Nedelec<dim>(0), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);
    test_vec.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 3, 0, mapping, test_vec);
    if (test_boundary_values(dof_handler, mapping, fe_system, 3, 0, test_vec))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
    // Now, in the case of only one element both initializations should give
    // identical vectors
    deallog << "Checking Consistency: ";
    bool equal = (test_vec.size() == test_vec2.size());
    if (equal)
      {
        for (unsigned int i = 0; i < test_vec.size(); ++i)
          {
            equal = equal && (test_vec(i) == test_vec2(i));
          }
      }
    if (equal)
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 1b: Only FE_Nedelec
    deallog << "Testing FE_Nedelec(1): ";
    FE_Nedelec<dim> fe(1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe);
    test_vec2.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 3, 0, mapping, test_vec2);
    if (test_boundary_values(dof_handler, mapping, fe, 3, 0, test_vec2))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 2b: FESystem(FE_Nedelec(1))
    deallog << "Testing FESystem(FE_Nedelec(1)): ";
    FESystem<dim>   fe_system(FE_Nedelec<dim>(1), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);
    test_vec.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 3, 0, mapping, test_vec);
    if (test_boundary_values(dof_handler, mapping, fe_system, 3, 0, test_vec))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
    // Now, in the case of only one element both initializations should give
    // identical vectors
    deallog << "Checking Consistency: ";
    bool equal = (test_vec.size() == test_vec2.size());
    if (equal)
      {
        for (unsigned int i = 0; i < test_vec.size(); ++i)
          {
            equal = equal && (test_vec(i) == test_vec2(i));
          }
      }
    if (equal)
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 3: FESystem(FE_Nedelec(0), FE_Nedelec(0))
    deallog << "Testing FESystem(FE_Nedelec(0), FE_Nedelec(0)): ";
    FESystem<dim>   fe_system(FE_Nedelec<dim>(0), 1, FE_Nedelec<dim>(0), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);
    test_vec.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 6, 3, mapping, test_vec);
    if (test_boundary_values(dof_handler, mapping, fe_system, 6, 3, test_vec))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 4: FESystem(FE_Nedelec(1), FE_Nedelec(0))
    deallog << "Testing FESystem(FE_Nedelec(1), FE_Nedelec(0)): ";
    FESystem<dim>   fe_system(FE_Nedelec<dim>(1), 1, FE_Nedelec<dim>(0), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);
    test_vec.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 6, 3, mapping, test_vec);
    if (test_boundary_values(dof_handler, mapping, fe_system, 6, 3, test_vec))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 5: FESystem(FE_Nedelec(0), FE_Nedelec(1))
    deallog << "Testing FESystem(FE_Nedelec(0), FE_Nedelec(1)): ";
    FESystem<dim>   fe_system(FE_Nedelec<dim>(0), 1, FE_Nedelec<dim>(1), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);
    test_vec.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 6, 3, mapping, test_vec);
    if (test_boundary_values(dof_handler, mapping, fe_system, 6, 3, test_vec))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 6: FESystem(FE_Nedelec(0), FE_Q(1))
    deallog << "Testing FESystem(FE_Nedelec(0), FE_Q(1)): ";
    FESystem<dim>   fe_system(FE_Nedelec<dim>(0), 1, FE_Q<dim>(1), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);
    test_vec.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 4, 0, mapping, test_vec);
    if (test_boundary_values(dof_handler, mapping, fe_system, 4, 0, test_vec))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 7: FESystem(FE_Q(1), FE_Nedelec(0))
    deallog << "Testing FESystem(FE_Q(1), FE_Nedelec(0)): ";
    FESystem<dim>   fe_system(FE_Q<dim>(1), 1, FE_Nedelec<dim>(0), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);
    test_vec.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 4, 1, mapping, test_vec);
    if (test_boundary_values(dof_handler, mapping, fe_system, 4, 1, test_vec))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 8: FESystem(FE_Nedelec(0), FE_Q(2))
    deallog << "Testing FESystem(FE_Nedelec(0), FE_Q(2)): ";
    FESystem<dim>   fe_system(FE_Nedelec<dim>(0), 1, FE_Q<dim>(2), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);
    test_vec.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 4, 0, mapping, test_vec);
    if (test_boundary_values(dof_handler, mapping, fe_system, 4, 0, test_vec))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 9: FESystem(FE_Q(2), FE_Nedelec(0))
    deallog << "Testing FESystem(FE_Q(2), FE_Nedelec(0)): ";
    FESystem<dim>   fe_system(FE_Q<dim>(2), 1, FE_Nedelec<dim>(0), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);
    test_vec.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 4, 1, mapping, test_vec);
    if (test_boundary_values(dof_handler, mapping, fe_system, 4, 1, test_vec))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 10: FESystem(FE_Nedelec(1), FE_Q(1))
    deallog << "Testing FESystem(FE_Nedelec(1), FE_Q(1)): ";
    FESystem<dim>   fe_system(FE_Nedelec<dim>(1), 1, FE_Q<dim>(1), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);
    test_vec.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 4, 0, mapping, test_vec);
    if (test_boundary_values(dof_handler, mapping, fe_system, 4, 0, test_vec))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }
  {
    // Test 11: FESystem(FE_Q(1), FE_Nedelec(1))
    deallog << "Testing FESystem(FE_Q(1), FE_Nedelec(1)): ";
    FESystem<dim>   fe_system(FE_Q<dim>(1), 1, FE_Nedelec<dim>(1), 1);
    DoFHandler<dim> dof_handler(triangulation);
    dof_handler.distribute_dofs(fe_system);
    test_vec.reinit(dof_handler.n_dofs());
    apply_boundary_values(dof_handler, constraints, 4, 1, mapping, test_vec);
    if (test_boundary_values(dof_handler, mapping, fe_system, 4, 1, test_vec))
      {
        deallog << "OK" << std::endl;
      }
    else
      {
        deallog << "Failed" << std::endl;
        abort();
      }
  }

  return 0;
}
