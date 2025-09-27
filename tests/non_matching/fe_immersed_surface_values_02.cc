#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_bernstein.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_fe_field.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/mapping_collection.h>

#include <deal.II/non_matching/fe_immersed_values.h>
#include <deal.II/non_matching/immersed_surface_quadrature.h>

#include <deal.II/numerics/vector_tools_interpolate.h>

#include "../tests.h"


using NonMatching::FEImmersedSurfaceValues;

// Set up a triangulation that we can construct an
// FEImmersedSurfaceValues object and compute known values when using a
// MappingFEField.


template <int dim>
class Test
{
public:
  Test();

  void
  run();

private:
  //  Construct a triangulation with a single  cell
  void
  setup_single_cell_triangulation();

  // Make the ImmersedSurfaceQuadraturewith contain a single point in the center
  // of the cell, with weight .5 and a normal different in all components.
  void
  setup_single_point_quadrature();

  // Test that the constant n_quadrature_points corresponds to what we sent in.
  void
  test_n_quadrature_points();

  // Test that we can use the get_quadrature function to get the quadrature
  // we passed to the FEImmersedSurfaceValues constructor.
  void
  test_get_quadrature();

  // Print the quadrature point in real space to make sure that it's mapped to
  // real space correctly.
  void
  test_point_mapped_correctly();

  // Print the normal in real space to make sure that it's mapped to real space
  // correctly.
  void
  test_normal();

  // Print JxW to check that the value is correct.
  void
  test_JxW();

  // Test that we can call the shape_surface_grad function.
  void
  test_shape_surface_grad();

  const FE_Bernstein<dim>                     element;
  const FESystem<dim>                         fe_system;
  Triangulation<dim>                          triangulation;
  DoFHandler<dim>                             dof_handler;
  DoFHandler<dim>                             euler_dof_handler;
  Vector<double>                              euler_vector;
  std::unique_ptr<MappingFEField<dim>>        euler_mapping;
  NonMatching::ImmersedSurfaceQuadrature<dim> quadrature;
};



template <int dim>
Test<dim>::Test()
  : element(1)
  , fe_system(element, dim)
  , dof_handler(triangulation)
  , euler_dof_handler(triangulation)
{}



template <int dim>
void
Test<dim>::run()
{
  setup_single_cell_triangulation();

  setup_single_point_quadrature();

  test_n_quadrature_points();
  test_get_quadrature();
  test_point_mapped_correctly();
  test_normal();
  test_JxW();
  test_shape_surface_grad();
}



template <int dim>
void
Test<dim>::setup_single_cell_triangulation()
{
  const Point<dim> lower_left;
  Point<dim>       upper_right;
  for (unsigned int d = 0; d < dim; ++d)
    {
      upper_right[d] = 2 + d;
    }

  GridGenerator::hyper_rectangle(triangulation, lower_left, upper_right);
  dof_handler.distribute_dofs(element);
  euler_dof_handler.distribute_dofs(fe_system);
  euler_vector.reinit(euler_dof_handler.n_dofs());

  const ComponentMask mask(dim, true);
  VectorTools::get_position_vector(euler_dof_handler, euler_vector, mask);
  euler_mapping = std::make_unique<MappingFEField<dim>>(euler_dof_handler,
                                                        euler_vector,
                                                        mask);
}



template <int dim>
void
Test<dim>::setup_single_point_quadrature()
{
  Point<dim>     point;
  const double   weight = .5;
  Tensor<1, dim> normal;
  for (unsigned int d = 0; d < dim; ++d)
    {
      point[d]  = .5;
      normal[d] = d + 1;
    }
  normal /= normal.norm();
  quadrature.push_back(point, weight, normal);
}



template <int dim>
void
Test<dim>::test_n_quadrature_points()
{
  FEImmersedSurfaceValues<dim> fe_values(*euler_mapping,
                                         element,
                                         quadrature,
                                         update_default);
  fe_values.reinit(triangulation.begin_active());

  deallog << "n_quadrature_points = " << fe_values.n_quadrature_points
          << std::endl;
}



template <int dim>
void
Test<dim>::test_get_quadrature()
{
  FEImmersedSurfaceValues<dim> fe_values(*euler_mapping,
                                         element,
                                         quadrature,
                                         update_default);
  fe_values.reinit(triangulation.begin_active());

  const NonMatching::ImmersedSurfaceQuadrature<dim> &stored_quadrature =
    fe_values.get_quadrature();

  for (unsigned int q = 0; q < stored_quadrature.size(); q++)
    {
      deallog << "(point, weight, normal) = ([" << stored_quadrature.point(q)
              << "], " << stored_quadrature.weight(q) << ", ["
              << stored_quadrature.normal_vector(q) << "])" << std::endl;
    }
}



template <int dim>
void
Test<dim>::test_point_mapped_correctly()
{
  FEImmersedSurfaceValues<dim> fe_values(*euler_mapping,
                                         element,
                                         quadrature,
                                         update_quadrature_points);
  fe_values.reinit(triangulation.begin_active());

  deallog << "point = " << fe_values.quadrature_point(0) << std::endl;
}



template <int dim>
void
Test<dim>::test_normal()
{
  FEImmersedSurfaceValues<dim> fe_values(*euler_mapping,
                                         element,
                                         quadrature,
                                         update_normal_vectors);
  fe_values.reinit(triangulation.begin_active());

  deallog << "normal = " << fe_values.normal_vector(0) << std::endl;
}



template <int dim>
void
Test<dim>::test_JxW()
{
  FEImmersedSurfaceValues<dim> fe_values(*euler_mapping,
                                         element,
                                         quadrature,
                                         update_JxW_values);
  fe_values.reinit(triangulation.begin_active());

  deallog << "JxW = " << fe_values.JxW(0) << std::endl;
}



template <int dim>
void
Test<dim>::test_shape_surface_grad()
{
  FEImmersedSurfaceValues<dim> fe_values(*euler_mapping,
                                         element,
                                         quadrature,
                                         update_gradients |
                                           update_normal_vectors);
  fe_values.reinit(dof_handler.begin_active());

  const unsigned int function_index = 0;
  const unsigned int q_index        = 0;
  const unsigned int component      = 0;

  deallog << "shape_surface_grad = "
          << fe_values.shape_surface_grad(function_index, q_index) << std::endl;

  deallog << "shape_surface_grad_component = "
          << fe_values.shape_surface_grad_component(function_index,
                                                    q_index,
                                                    component)
          << std::endl;
}



template <int dim>
void
run_test()
{
  deallog << "dim = " << dim << std::endl;
  Test<dim> test;
  test.run();
  deallog << std::endl;
}



int
main()
{
  initlog();
  run_test<1>();
  run_test<2>();
  run_test<3>();
}
