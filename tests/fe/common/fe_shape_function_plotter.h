// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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


#include <deal.II/base/function.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/fe/fe_nedelec_sz.h>
#include <deal.II/fe/fe_pyramid_p.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_wedge_p.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/reference_cell.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>
#include <string>

using namespace dealii;


/**
 * A helper struct that is used indicate whether we need to use
 * FEExtractors in order to retrieve shape function information
 * from the finite elements. This is typically required for finite
 * elements with non-primitive shape functions.
 */
template <class FE_Type>
struct RequiresFEExtractor : std::false_type
{};

/**
 * A specialisation of the above classes, for the FE_NedelecSZ class.
 */
template <int dim, int spacedim>
struct RequiresFEExtractor<FE_NedelecSZ<dim, spacedim>> : std::true_type
{};


/**
 * A specialisation of the above classes, for the FE_RaviartThomas class.
 */
template <int dim>
struct RequiresFEExtractor<FE_RaviartThomas<dim>> : std::true_type
{};


/**
 * A function that returns the shape function of a finite element.
 *
 * This class uses directly requests the shape function components
 * from the finite elements. Not all classes support this function
 * call, so we have to take this into consideration.
 */
template <typename FE_Type>
class FEComponentFunction : public Function<FE_Type::dimension>
{
public:
  static const unsigned int dim = FE_Type::dimension;

  FEComponentFunction(const FE_Type &      fe,
                      const unsigned int   basis_function,
                      const ReferenceCell &reference_cell)
    : Function<dim>(fe.n_components())
    , fe(fe)
    , bf(basis_function)
    , reference_cell(reference_cell)
  {
    Assert(basis_function < fe.n_dofs_per_cell(), ExcInternalError());
  }
  virtual ~FEComponentFunction() = default;

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    Assert(values.size() == this->n_components, ExcInternalError());

    // Some functionality does not yet exist for transition elements
    const bool is_pyramid = (reference_cell == ReferenceCells::Pyramid);
    const bool is_wedge   = (reference_cell == ReferenceCells::Wedge);
    if (!is_pyramid && !is_wedge)
      {
        Assert(reference_cell.contains_point(p), ExcInternalError());
      }

    // The main triangulation defines exactly the reference
    // element, so we can therefore directly ask the FE
    // for any shape functions etc.
    for (unsigned int i = 0; i < this->n_components; ++i)
      values[i] = fe.shape_value_component(bf, p, i);
  }

private:
  const FE_Type &      fe;
  const unsigned int   bf;
  const ReferenceCell &reference_cell;
};


/**
 * A function that returns the shape function of a finite element.
 *
 * This class uses FEExtractors to get the shape function information
 * from the finite element. Doing so is significantly more expensive
 * than the methodology used in the other variant of this class, so it
 * should be used only when necessary.
 */
template <typename FE_Type>
class FEExtractorFunction : public Function<FE_Type::dimension>
{
public:
  static const unsigned int dim = FE_Type::dimension;

  FEExtractorFunction(const FE_Type &      fe,
                      const unsigned int   basis_function,
                      const ReferenceCell &reference_cell)
    : Function<dim>(fe.n_components())
    , fe(fe)
    , bf(basis_function)
    , dof(tr)
    , fe_s(0)
    , fe_v(0)
  {
    Assert(basis_function < fe.n_dofs_per_cell(), ExcInternalError());

    GridGenerator::reference_cell(tr, reference_cell);
    dof.distribute_dofs(fe);
  }
  virtual ~FEExtractorFunction() = default;

  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const override
  {
    Assert(values.size() == this->n_components, ExcInternalError());
    Assert(GeometryInfo<dim>::is_inside_unit_cell(p), ExcInternalError());

    // For elements that don't support component-wise extraction, we use
    // the  FEValuesExtractors / FEValuesViews in conjunction with a
    // custom quadrature rule on the unit cell.
    const Quadrature<dim> quad(p);
    FEValues<dim>         fe_values(fe, quad, update_values);
    Assert(fe_values.get_quadrature().size() == 1, ExcInternalError());

    Assert(tr.n_active_cells() == 1, ExcInternalError());
    fe_values.reinit(dof.begin_active());
    constexpr unsigned int q_point = 0;

    if (this->n_components == 1)
      {
        const double Nx = fe_values[fe_s].value(bf, q_point);
        values[0]       = Nx;
      }
    else if (this->n_components == dim)
      {
        const Tensor<1, dim> Nx = fe_values[fe_v].value(bf, q_point);
        for (unsigned int i = 0; i < this->n_components; ++i)
          values[i] = Nx[i];
      }
    else
      {
        AssertThrow(false, ExcNotImplemented());
      }
  }

private:
  const FE_Type &    fe;
  const unsigned int bf;

  Triangulation<dim>         tr;
  DoFHandler<dim>            dof;
  FEValuesExtractors::Scalar fe_s;
  FEValuesExtractors::Vector fe_v;
};


/**
 * The main function that is used to plot finite element shape
 * functions. They are output in VTU format, so that they can
 * be visualised in Paraview and Visit.
 */
template <typename FE_Type>
void
plot_FE_shape_functions(const std::string    filename,
                        const unsigned int   order,
                        const ReferenceCell &reference_cell,
                        const unsigned int   n_div)
{
  const unsigned int dim = FE_Type::dimension;

  // Vector valued finite element to be visualised
  const FE_Type      fe_viz(order);
  const unsigned int n_components = fe_viz.n_components();
  Assert((n_components == 1) || (n_components == dim), ExcInternalError());

  // FESystem for the vector-valued FE to be visualised
  const bool is_hyper_cube = reference_cell.is_hyper_cube();
  const bool is_simplex    = reference_cell.is_simplex();
  const bool is_pyramid    = (reference_cell == ReferenceCells::Pyramid);
  const bool is_wedge      = (reference_cell == ReferenceCells::Wedge);
  if (!is_hyper_cube && !is_simplex && !is_pyramid)
    {
      (void)is_wedge;
      Assert(is_wedge, ExcInternalError());
    }

  const FESystem<dim> fe(
    is_hyper_cube ?
      FESystem<dim>(FE_Q<dim>(1), n_components) :
      (is_simplex ?
         FESystem<dim>(FE_SimplexP<dim>(1), n_components) :
         (is_pyramid ? FESystem<dim>(FE_PyramidP<dim>(1), n_components) :
                       FESystem<dim>(FE_WedgeP<dim>(1), n_components))));

  Triangulation<dim> tr;
  GridGenerator::reference_cell(tr, reference_cell);
  tr.refine_global(n_div);

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  // A collection of solution to store solution fields
  // characterising all vector valued basis functions
  std::vector<Vector<double>> bf_solns(fe_viz.n_dofs_per_cell(),
                                       Vector<double>(dof.n_dofs()));

  // Interpolate the vector valued shape functions onto a
  // continuous vector-valued FE space
  for (unsigned int bf = 0; bf < fe_viz.n_dofs_per_cell(); ++bf)
    if (RequiresFEExtractor<FE_Type>::value == false)
      {
        VectorTools::interpolate(dof,
                                 FEComponentFunction<FE_Type>(fe_viz,
                                                              bf,
                                                              reference_cell),
                                 bf_solns[bf]);
      }
    else
      {
        VectorTools::interpolate(dof,
                                 FEExtractorFunction<FE_Type>(fe_viz,
                                                              bf,
                                                              reference_cell),
                                 bf_solns[bf]);
      }

  // Plot the result, taking care of the fact that we might have either a
  // scalar- or vector-valued finite element
  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof);
  std::vector<DataComponentInterpretation::DataComponentInterpretation> dci(
    n_components,
    (n_components == 1 ?
       DataComponentInterpretation::component_is_scalar :
       DataComponentInterpretation::component_is_part_of_vector));

  for (unsigned int bf = 0; bf < fe_viz.n_dofs_per_cell(); ++bf)
    {
      const std::vector<std::string> bf_name(n_components,
                                             "basis_function_" +
                                               Utilities::to_string(bf, 2));
      data_out.add_data_vector(bf_solns[bf],
                               bf_name,
                               DataOut<dim>::type_dof_data,
                               dci);
    }

  data_out.build_patches();
  const std::string prefix = filename + "-order_" + Utilities::to_string(order);
  std::ofstream     output(dim == 2 ? prefix + "-2d.vtu" : prefix + "-3d.vtu");
  data_out.write_vtu(output);
}


template <typename FE_Type>
void
plot_one(const std::string    name,
         const unsigned int   lowest_order,
         const unsigned int   highest_order,
         const ReferenceCell &reference_cell,
         const unsigned int   n_divs = 4)
{
  for (unsigned int order = lowest_order; order <= highest_order; ++order)
    {
      std::cout << "  Plotting element type " << name << " of order " << order
                << "..." << std::endl;
      plot_FE_shape_functions<FE_Type>(name, order, reference_cell, n_divs);
    }
}
