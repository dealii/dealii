// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <memory>
#include <sstream>


DEAL_II_NAMESPACE_OPEN

namespace
{
  unsigned int
  count_nonzeros(const std::vector<unsigned int> &vec)
  {
    return std::count_if(vec.begin(), vec.end(), [](const unsigned int i) {
      return i > 0;
    });
  }
} // namespace
/* ----------------------- FESystem::InternalData ------------------- */


template <int dim, int spacedim>
FESystem<dim, spacedim>::InternalData::InternalData(
  const unsigned int n_base_elements)
  : base_fe_datas(n_base_elements)
  , base_fe_output_objects(n_base_elements)
{}



template <int dim, int spacedim>
FESystem<dim, spacedim>::InternalData::~InternalData()
{
  // delete pointers and set them to zero to avoid inadvertent use
  for (unsigned int i = 0; i < base_fe_datas.size(); ++i)
    base_fe_datas[i].reset();
}


template <int dim, int spacedim>
typename FiniteElement<dim, spacedim>::InternalDataBase &
FESystem<dim, spacedim>::InternalData::get_fe_data(
  const unsigned int base_no) const
{
  AssertIndexRange(base_no, base_fe_datas.size());
  return *base_fe_datas[base_no];
}



template <int dim, int spacedim>
void
FESystem<dim, spacedim>::InternalData::set_fe_data(
  const unsigned int base_no,
  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase> ptr)
{
  AssertIndexRange(base_no, base_fe_datas.size());
  base_fe_datas[base_no] = std::move(ptr);
}



template <int dim, int spacedim>
internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim> &
FESystem<dim, spacedim>::InternalData::get_fe_output_object(
  const unsigned int base_no) const
{
  AssertIndexRange(base_no, base_fe_output_objects.size());
  return base_fe_output_objects[base_no];
}



/* ---------------------------------- FESystem ------------------- */


template <int dim, int spacedim>
const unsigned int FESystem<dim, spacedim>::invalid_face_number;


template <int dim, int spacedim>
FESystem<dim, spacedim>::FESystem(const FiniteElement<dim, spacedim> &fe,
                                  const unsigned int n_elements)
  : FiniteElement<dim, spacedim>(
      FETools::Compositing::multiply_dof_numbers(&fe, n_elements),
      FETools::Compositing::compute_restriction_is_additive_flags(&fe,
                                                                  n_elements),
      FETools::Compositing::compute_nonzero_components(&fe, n_elements))
  , base_elements((n_elements > 0))
{
  std::vector<const FiniteElement<dim, spacedim> *> fes;
  fes.push_back(&fe);
  std::vector<unsigned int> multiplicities;
  multiplicities.push_back(n_elements);
  initialize(fes, multiplicities);
}



template <int dim, int spacedim>
FESystem<dim, spacedim>::FESystem(const FiniteElement<dim, spacedim> &fe1,
                                  const unsigned int                  n1,
                                  const FiniteElement<dim, spacedim> &fe2,
                                  const unsigned int                  n2)
  : FiniteElement<dim, spacedim>(
      FETools::Compositing::multiply_dof_numbers(&fe1, n1, &fe2, n2),
      FETools::Compositing::compute_restriction_is_additive_flags(&fe1,
                                                                  n1,
                                                                  &fe2,
                                                                  n2),
      FETools::Compositing::compute_nonzero_components(&fe1, n1, &fe2, n2))
  , base_elements((n1 > 0) + (n2 > 0))
{
  std::vector<const FiniteElement<dim, spacedim> *> fes;
  fes.push_back(&fe1);
  fes.push_back(&fe2);
  std::vector<unsigned int> multiplicities;
  multiplicities.push_back(n1);
  multiplicities.push_back(n2);
  initialize(fes, multiplicities);
}



template <int dim, int spacedim>
FESystem<dim, spacedim>::FESystem(const FiniteElement<dim, spacedim> &fe1,
                                  const unsigned int                  n1,
                                  const FiniteElement<dim, spacedim> &fe2,
                                  const unsigned int                  n2,
                                  const FiniteElement<dim, spacedim> &fe3,
                                  const unsigned int                  n3)
  : FiniteElement<dim, spacedim>(
      FETools::Compositing::multiply_dof_numbers(&fe1, n1, &fe2, n2, &fe3, n3),
      FETools::Compositing::compute_restriction_is_additive_flags(&fe1,
                                                                  n1,
                                                                  &fe2,
                                                                  n2,
                                                                  &fe3,
                                                                  n3),
      FETools::Compositing::compute_nonzero_components(&fe1,
                                                       n1,
                                                       &fe2,
                                                       n2,
                                                       &fe3,
                                                       n3))
  , base_elements((n1 > 0) + (n2 > 0) + (n3 > 0))
{
  std::vector<const FiniteElement<dim, spacedim> *> fes;
  fes.push_back(&fe1);
  fes.push_back(&fe2);
  fes.push_back(&fe3);
  std::vector<unsigned int> multiplicities;
  multiplicities.push_back(n1);
  multiplicities.push_back(n2);
  multiplicities.push_back(n3);
  initialize(fes, multiplicities);
}



template <int dim, int spacedim>
FESystem<dim, spacedim>::FESystem(const FiniteElement<dim, spacedim> &fe1,
                                  const unsigned int                  n1,
                                  const FiniteElement<dim, spacedim> &fe2,
                                  const unsigned int                  n2,
                                  const FiniteElement<dim, spacedim> &fe3,
                                  const unsigned int                  n3,
                                  const FiniteElement<dim, spacedim> &fe4,
                                  const unsigned int                  n4)
  : FiniteElement<dim, spacedim>(
      FETools::Compositing::multiply_dof_numbers(&fe1,
                                                 n1,
                                                 &fe2,
                                                 n2,
                                                 &fe3,
                                                 n3,
                                                 &fe4,
                                                 n4),
      FETools::Compositing::compute_restriction_is_additive_flags(&fe1,
                                                                  n1,
                                                                  &fe2,
                                                                  n2,
                                                                  &fe3,
                                                                  n3,
                                                                  &fe4,
                                                                  n4),
      FETools::Compositing::compute_nonzero_components(&fe1,
                                                       n1,
                                                       &fe2,
                                                       n2,
                                                       &fe3,
                                                       n3,
                                                       &fe4,
                                                       n4))
  , base_elements((n1 > 0) + (n2 > 0) + (n3 > 0) + (n4 > 0))
{
  std::vector<const FiniteElement<dim, spacedim> *> fes;
  fes.push_back(&fe1);
  fes.push_back(&fe2);
  fes.push_back(&fe3);
  fes.push_back(&fe4);
  std::vector<unsigned int> multiplicities;
  multiplicities.push_back(n1);
  multiplicities.push_back(n2);
  multiplicities.push_back(n3);
  multiplicities.push_back(n4);
  initialize(fes, multiplicities);
}



template <int dim, int spacedim>
FESystem<dim, spacedim>::FESystem(const FiniteElement<dim, spacedim> &fe1,
                                  const unsigned int                  n1,
                                  const FiniteElement<dim, spacedim> &fe2,
                                  const unsigned int                  n2,
                                  const FiniteElement<dim, spacedim> &fe3,
                                  const unsigned int                  n3,
                                  const FiniteElement<dim, spacedim> &fe4,
                                  const unsigned int                  n4,
                                  const FiniteElement<dim, spacedim> &fe5,
                                  const unsigned int                  n5)
  : FiniteElement<dim, spacedim>(
      FETools::Compositing::
        multiply_dof_numbers(&fe1, n1, &fe2, n2, &fe3, n3, &fe4, n4, &fe5, n5),
      FETools::Compositing::compute_restriction_is_additive_flags(&fe1,
                                                                  n1,
                                                                  &fe2,
                                                                  n2,
                                                                  &fe3,
                                                                  n3,
                                                                  &fe4,
                                                                  n4,
                                                                  &fe5,
                                                                  n5),
      FETools::Compositing::compute_nonzero_components(&fe1,
                                                       n1,
                                                       &fe2,
                                                       n2,
                                                       &fe3,
                                                       n3,
                                                       &fe4,
                                                       n4,
                                                       &fe5,
                                                       n5))
  , base_elements((n1 > 0) + (n2 > 0) + (n3 > 0) + (n4 > 0) + (n5 > 0))
{
  std::vector<const FiniteElement<dim, spacedim> *> fes;
  fes.push_back(&fe1);
  fes.push_back(&fe2);
  fes.push_back(&fe3);
  fes.push_back(&fe4);
  fes.push_back(&fe5);
  std::vector<unsigned int> multiplicities;
  multiplicities.push_back(n1);
  multiplicities.push_back(n2);
  multiplicities.push_back(n3);
  multiplicities.push_back(n4);
  multiplicities.push_back(n5);
  initialize(fes, multiplicities);
}



template <int dim, int spacedim>
FESystem<dim, spacedim>::FESystem(
  const std::vector<const FiniteElement<dim, spacedim> *> &fes,
  const std::vector<unsigned int> &                        multiplicities)
  : FiniteElement<dim, spacedim>(
      FETools::Compositing::multiply_dof_numbers(fes, multiplicities),
      FETools::Compositing::compute_restriction_is_additive_flags(
        fes,
        multiplicities),
      FETools::Compositing::compute_nonzero_components(fes, multiplicities))
  , base_elements(count_nonzeros(multiplicities))
{
  initialize(fes, multiplicities);
}



template <int dim, int spacedim>
std::string
FESystem<dim, spacedim>::get_name() const
{
  // note that the
  // FETools::get_fe_by_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;

  namebuf << "FESystem<" << Utilities::dim_string(dim, spacedim) << ">[";
  for (unsigned int i = 0; i < this->n_base_elements(); ++i)
    {
      namebuf << base_element(i).get_name();
      if (this->element_multiplicity(i) != 1)
        namebuf << '^' << this->element_multiplicity(i);
      if (i != this->n_base_elements() - 1)
        namebuf << '-';
    }
  namebuf << ']';

  return namebuf.str();
}



template <int dim, int spacedim>
std::unique_ptr<FiniteElement<dim, spacedim>>
FESystem<dim, spacedim>::clone() const
{
  std::vector<const FiniteElement<dim, spacedim> *> fes;
  std::vector<unsigned int>                         multiplicities;

  for (unsigned int i = 0; i < this->n_base_elements(); i++)
    {
      fes.push_back(&base_element(i));
      multiplicities.push_back(this->element_multiplicity(i));
    }
  return std::make_unique<FESystem<dim, spacedim>>(fes, multiplicities);
}



template <int dim, int spacedim>
const FiniteElement<dim, spacedim> &
FESystem<dim, spacedim>::get_sub_fe(
  const unsigned int first_component,
  const unsigned int n_selected_components) const
{
  Assert(first_component + n_selected_components <= this->n_components(),
         ExcMessage("Invalid arguments (not a part of this FiniteElement)."));

  const unsigned int base_index =
    this->component_to_base_table[first_component].first.first;
  const unsigned int component_in_base =
    this->component_to_base_table[first_component].first.second;
  const unsigned int base_components =
    this->base_element(base_index).n_components();

  // Only select our child base_index if that is all the user wanted. Error
  // handling will be done inside the recursion.
  if (n_selected_components <= base_components)
    return this->base_element(base_index)
      .get_sub_fe(component_in_base, n_selected_components);

  Assert(n_selected_components == this->n_components(),
         ExcMessage("You can not select a part of a FiniteElement."));
  return *this;
}



template <int dim, int spacedim>
double
FESystem<dim, spacedim>::shape_value(const unsigned int i,
                                     const Point<dim> & p) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  Assert(this->is_primitive(i),
         (typename FiniteElement<dim, spacedim>::ExcShapeFunctionNotPrimitive(
           i)));

  return (base_element(this->system_to_base_table[i].first.first)
            .shape_value(this->system_to_base_table[i].second, p));
}



template <int dim, int spacedim>
double
FESystem<dim, spacedim>::shape_value_component(
  const unsigned int i,
  const Point<dim> & p,
  const unsigned int component) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, this->n_components());

  // if this value is supposed to be
  // zero, then return right away...
  if (this->nonzero_components[i][component] == false)
    return 0;

  // ...otherwise: first find out to
  // which of the base elements this
  // desired component belongs, and
  // which component within this base
  // element it is
  const unsigned int base = this->component_to_base_index(component).first;
  const unsigned int component_in_base =
    this->component_to_base_index(component).second;

  // then get value from base
  // element. note that that will
  // throw an error should the
  // respective shape function not be
  // primitive; thus, there is no
  // need to check this here
  return (base_element(base).shape_value_component(
    this->system_to_base_table[i].second, p, component_in_base));
}



template <int dim, int spacedim>
Tensor<1, dim>
FESystem<dim, spacedim>::shape_grad(const unsigned int i,
                                    const Point<dim> & p) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  Assert(this->is_primitive(i),
         (typename FiniteElement<dim, spacedim>::ExcShapeFunctionNotPrimitive(
           i)));

  return (base_element(this->system_to_base_table[i].first.first)
            .shape_grad(this->system_to_base_table[i].second, p));
}



template <int dim, int spacedim>
Tensor<1, dim>
FESystem<dim, spacedim>::shape_grad_component(
  const unsigned int i,
  const Point<dim> & p,
  const unsigned int component) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, this->n_components());

  // if this value is supposed to be zero, then return right away...
  if (this->nonzero_components[i][component] == false)
    return Tensor<1, dim>();

  // ...otherwise: first find out to which of the base elements this desired
  // component belongs, and which component within this base element it is
  const unsigned int base = this->component_to_base_index(component).first;
  const unsigned int component_in_base =
    this->component_to_base_index(component).second;

  // then get value from base element. note that that will throw an error
  // should the respective shape function not be primitive; thus, there is no
  // need to check this here
  return (base_element(base).shape_grad_component(
    this->system_to_base_table[i].second, p, component_in_base));
}



template <int dim, int spacedim>
Tensor<2, dim>
FESystem<dim, spacedim>::shape_grad_grad(const unsigned int i,
                                         const Point<dim> & p) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  Assert(this->is_primitive(i),
         (typename FiniteElement<dim, spacedim>::ExcShapeFunctionNotPrimitive(
           i)));

  return (base_element(this->system_to_base_table[i].first.first)
            .shape_grad_grad(this->system_to_base_table[i].second, p));
}



template <int dim, int spacedim>
Tensor<2, dim>
FESystem<dim, spacedim>::shape_grad_grad_component(
  const unsigned int i,
  const Point<dim> & p,
  const unsigned int component) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, this->n_components());

  // if this value is supposed to be zero, then return right away...
  if (this->nonzero_components[i][component] == false)
    return Tensor<2, dim>();

  // ...otherwise: first find out to which of the base elements this desired
  // component belongs, and which component within this base element it is
  const unsigned int base = this->component_to_base_index(component).first;
  const unsigned int component_in_base =
    this->component_to_base_index(component).second;

  // then get value from base element. note that that will throw an error
  // should the respective shape function not be primitive; thus, there is no
  // need to check this here
  return (base_element(base).shape_grad_grad_component(
    this->system_to_base_table[i].second, p, component_in_base));
}



template <int dim, int spacedim>
Tensor<3, dim>
FESystem<dim, spacedim>::shape_3rd_derivative(const unsigned int i,
                                              const Point<dim> & p) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  Assert(this->is_primitive(i),
         (typename FiniteElement<dim, spacedim>::ExcShapeFunctionNotPrimitive(
           i)));

  return (base_element(this->system_to_base_table[i].first.first)
            .shape_3rd_derivative(this->system_to_base_table[i].second, p));
}



template <int dim, int spacedim>
Tensor<3, dim>
FESystem<dim, spacedim>::shape_3rd_derivative_component(
  const unsigned int i,
  const Point<dim> & p,
  const unsigned int component) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, this->n_components());

  // if this value is supposed to be zero, then return right away...
  if (this->nonzero_components[i][component] == false)
    return Tensor<3, dim>();

  // ...otherwise: first find out to which of the base elements this desired
  // component belongs, and which component within this base element it is
  const unsigned int base = this->component_to_base_index(component).first;
  const unsigned int component_in_base =
    this->component_to_base_index(component).second;

  // then get value from base element. note that that will throw an error
  // should the respective shape function not be primitive; thus, there is no
  // need to check this here
  return (base_element(base).shape_3rd_derivative_component(
    this->system_to_base_table[i].second, p, component_in_base));
}



template <int dim, int spacedim>
Tensor<4, dim>
FESystem<dim, spacedim>::shape_4th_derivative(const unsigned int i,
                                              const Point<dim> & p) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  Assert(this->is_primitive(i),
         (typename FiniteElement<dim, spacedim>::ExcShapeFunctionNotPrimitive(
           i)));

  return (base_element(this->system_to_base_table[i].first.first)
            .shape_4th_derivative(this->system_to_base_table[i].second, p));
}



template <int dim, int spacedim>
Tensor<4, dim>
FESystem<dim, spacedim>::shape_4th_derivative_component(
  const unsigned int i,
  const Point<dim> & p,
  const unsigned int component) const
{
  AssertIndexRange(i, this->n_dofs_per_cell());
  AssertIndexRange(component, this->n_components());

  // if this value is supposed to be zero, then return right away...
  if (this->nonzero_components[i][component] == false)
    return Tensor<4, dim>();

  // ...otherwise: first find out to which of the base elements this desired
  // component belongs, and which component within this base element it is
  const unsigned int base = this->component_to_base_index(component).first;
  const unsigned int component_in_base =
    this->component_to_base_index(component).second;

  // then get value from base element. note that that will throw an error
  // should the respective shape function not be primitive; thus, there is no
  // need to check this here
  return (base_element(base).shape_4th_derivative_component(
    this->system_to_base_table[i].second, p, component_in_base));
}



template <int dim, int spacedim>
void
FESystem<dim, spacedim>::get_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  FullMatrix<double> &                interpolation_matrix) const
{
  // check that the size of the matrices is correct. for historical
  // reasons, if you call matrix.reinit(8,0), it sets the sizes
  // to m==n==0 internally. this may happen when we use a FE_Nothing,
  // so write the test in a more lenient way
  Assert((interpolation_matrix.m() == this->n_dofs_per_cell()) ||
           (x_source_fe.n_dofs_per_cell() == 0),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              this->n_dofs_per_cell()));
  Assert((interpolation_matrix.n() == x_source_fe.n_dofs_per_cell()) ||
           (this->n_dofs_per_cell() == 0),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_cell()));

  // there are certain conditions that the two elements have to satisfy so
  // that this can work.
  //
  // condition 1: the other element must also be a system element

  AssertThrow(
    (x_source_fe.get_name().find("FESystem<") == 0) ||
      (dynamic_cast<const FESystem<dim, spacedim> *>(&x_source_fe) != nullptr),
    (typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented()));

  // ok, source is a system element, so we may be able to do the work
  const FESystem<dim, spacedim> &source_fe =
    dynamic_cast<const FESystem<dim, spacedim> &>(x_source_fe);

  // condition 2: same number of basis elements
  AssertThrow(
    this->n_base_elements() == source_fe.n_base_elements(),
    (typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented()));

  // condition 3: same number of basis elements
  for (unsigned int i = 0; i < this->n_base_elements(); ++i)
    AssertThrow(
      this->element_multiplicity(i) == source_fe.element_multiplicity(i),
      (typename FiniteElement<dim,
                              spacedim>::ExcInterpolationNotImplemented()));

  // ok, so let's try whether it works:

  // first let's see whether all the basis elements actually generate their
  // interpolation matrices. if we get past the following loop, then
  // apparently none of the called base elements threw an exception, so we're
  // fine continuing and assembling the one big matrix from the small ones of
  // the base elements
  std::vector<FullMatrix<double>> base_matrices(this->n_base_elements());
  for (unsigned int i = 0; i < this->n_base_elements(); ++i)
    {
      base_matrices[i].reinit(base_element(i).n_dofs_per_cell(),
                              source_fe.base_element(i).n_dofs_per_cell());
      base_element(i).get_interpolation_matrix(source_fe.base_element(i),
                                               base_matrices[i]);
    }

  // first clear big matrix, to make sure that entries that would couple
  // different bases (or multiplicity indices) are really zero. then assign
  // entries
  interpolation_matrix = 0;
  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    for (unsigned int j = 0; j < source_fe.n_dofs_per_cell(); ++j)
      if (this->system_to_base_table[i].first ==
          source_fe.system_to_base_table[j].first)
        interpolation_matrix(i, j) =
          (base_matrices[this->system_to_base_table[i].first.first](
            this->system_to_base_table[i].second,
            source_fe.system_to_base_table[j].second));
}



template <int dim, int spacedim>
const FullMatrix<double> &
FESystem<dim, spacedim>::get_restriction_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  AssertIndexRange(refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  Assert(refinement_case != RefinementCase<dim>::no_refinement,
         ExcMessage(
           "Restriction matrices are only available for refined cells!"));
  AssertIndexRange(child, GeometryInfo<dim>::n_children(refinement_case));

  // initialization upon first request
  if (this->restriction[refinement_case - 1][child].n() == 0)
    {
      std::lock_guard<std::mutex> lock(this->mutex);

      // check if updated while waiting for lock
      if (this->restriction[refinement_case - 1][child].n() ==
          this->n_dofs_per_cell())
        return this->restriction[refinement_case - 1][child];

      // Check if some of the matrices of the base elements are void.
      bool do_restriction = true;

      // shortcut for accessing local restrictions further down
      std::vector<const FullMatrix<double> *> base_matrices(
        this->n_base_elements());

      for (unsigned int i = 0; i < this->n_base_elements(); ++i)
        {
          base_matrices[i] =
            &base_element(i).get_restriction_matrix(child, refinement_case);
          if (base_matrices[i]->n() != base_element(i).n_dofs_per_cell())
            do_restriction = false;
        }
      Assert(do_restriction,
             (typename FiniteElement<dim, spacedim>::ExcProjectionVoid()));

      // if we did not encounter void matrices, initialize the matrix sizes
      if (do_restriction)
        {
          FullMatrix<double> restriction(this->n_dofs_per_cell(),
                                         this->n_dofs_per_cell());

          // distribute the matrices of the base finite elements to the
          // matrices of this object. for this, loop over all degrees of
          // freedom and take the respective entry of the underlying base
          // element.
          //
          // note that we by definition of a base element, they are
          // independent, i.e. do not couple. only DoFs that belong to the
          // same instance of a base element may couple
          for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
            for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
              {
                // first find out to which base element indices i and j
                // belong, and which instance thereof in case the base element
                // has a multiplicity greater than one. if they should not
                // happen to belong to the same instance of a base element,
                // then they cannot couple, so go on with the next index
                if (this->system_to_base_table[i].first !=
                    this->system_to_base_table[j].first)
                  continue;

                // so get the common base element and the indices therein:
                const unsigned int base =
                  this->system_to_base_table[i].first.first;

                const unsigned int base_index_i =
                                     this->system_to_base_table[i].second,
                                   base_index_j =
                                     this->system_to_base_table[j].second;

                // if we are sure that DoFs i and j may couple, then copy
                // entries of the matrices:
                restriction(i, j) =
                  (*base_matrices[base])(base_index_i, base_index_j);
              }

          restriction.swap(const_cast<FullMatrix<double> &>(
            this->restriction[refinement_case - 1][child]));
        }
    }

  return this->restriction[refinement_case - 1][child];
}



template <int dim, int spacedim>
const FullMatrix<double> &
FESystem<dim, spacedim>::get_prolongation_matrix(
  const unsigned int         child,
  const RefinementCase<dim> &refinement_case) const
{
  AssertIndexRange(refinement_case,
                   RefinementCase<dim>::isotropic_refinement + 1);
  Assert(refinement_case != RefinementCase<dim>::no_refinement,
         ExcMessage(
           "Restriction matrices are only available for refined cells!"));
  AssertIndexRange(child, GeometryInfo<dim>::n_children(refinement_case));

  // initialization upon first request, construction completely analogous to
  // restriction matrix
  if (this->prolongation[refinement_case - 1][child].n() == 0)
    {
      std::lock_guard<std::mutex> lock(this->mutex);

      if (this->prolongation[refinement_case - 1][child].n() ==
          this->n_dofs_per_cell())
        return this->prolongation[refinement_case - 1][child];

      bool                                    do_prolongation = true;
      std::vector<const FullMatrix<double> *> base_matrices(
        this->n_base_elements());
      for (unsigned int i = 0; i < this->n_base_elements(); ++i)
        {
          base_matrices[i] =
            &base_element(i).get_prolongation_matrix(child, refinement_case);
          if (base_matrices[i]->n() != base_element(i).n_dofs_per_cell())
            do_prolongation = false;
        }
      Assert(do_prolongation,
             (typename FiniteElement<dim, spacedim>::ExcEmbeddingVoid()));

      if (do_prolongation)
        {
          FullMatrix<double> prolongate(this->n_dofs_per_cell(),
                                        this->n_dofs_per_cell());

          for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
            for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
              {
                if (this->system_to_base_table[i].first !=
                    this->system_to_base_table[j].first)
                  continue;
                const unsigned int base =
                  this->system_to_base_table[i].first.first;

                const unsigned int base_index_i =
                                     this->system_to_base_table[i].second,
                                   base_index_j =
                                     this->system_to_base_table[j].second;
                prolongate(i, j) =
                  (*base_matrices[base])(base_index_i, base_index_j);
              }
          prolongate.swap(const_cast<FullMatrix<double> &>(
            this->prolongation[refinement_case - 1][child]));
        }
    }

  return this->prolongation[refinement_case - 1][child];
}


template <int dim, int spacedim>
unsigned int
FESystem<dim, spacedim>::face_to_cell_index(const unsigned int face_dof_index,
                                            const unsigned int face,
                                            const bool         face_orientation,
                                            const bool         face_flip,
                                            const bool face_rotation) const
{
  // we need to ask the base elements how they want to translate
  // the DoFs within their own numbering. thus, translate to
  // the base element numbering and then back
  const std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
    face_base_index = this->face_system_to_base_index(face_dof_index, face);

  const unsigned int base_face_to_cell_index =
    this->base_element(face_base_index.first.first)
      .face_to_cell_index(face_base_index.second,
                          face,
                          face_orientation,
                          face_flip,
                          face_rotation);

  // it would be nice if we had a base_to_system_index function, but
  // all that exists is a component_to_system_index function. we can't do
  // this here because it won't work for non-primitive elements. consequently,
  // simply do a loop over all dofs till we find whether it corresponds
  // to the one we're interested in -- crude, maybe, but works for now
  const std::pair<std::pair<unsigned int, unsigned int>, unsigned int> target =
    std::make_pair(face_base_index.first, base_face_to_cell_index);
  for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
    if (this->system_to_base_index(i) == target)
      return i;

  Assert(false, ExcInternalError());
  return numbers::invalid_unsigned_int;
}



//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------



template <int dim, int spacedim>
UpdateFlags
FESystem<dim, spacedim>::requires_update_flags(const UpdateFlags flags) const
{
  UpdateFlags out = update_default;
  // generate maximal set of flags
  // that are necessary
  for (unsigned int base_no = 0; base_no < this->n_base_elements(); ++base_no)
    out |= base_element(base_no).requires_update_flags(flags);
  return out;
}



template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FESystem<dim, spacedim>::get_data(
  const UpdateFlags             flags,
  const Mapping<dim, spacedim> &mapping,
  const Quadrature<dim> &       quadrature,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    & /*output_data*/) const
{
  // create an internal data object and set the update flags we will need
  // to deal with. the current object does not make use of these flags,
  // but we need to nevertheless set them correctly since we look
  // into the update_each flag of base elements in fill_fe_values,
  // and so the current object's update_each flag needs to be
  // correct in case the current FESystem is a base element for another,
  // higher-level FESystem itself.
  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
        data_ptr   = std::make_unique<InternalData>(this->n_base_elements());
  auto &data       = dynamic_cast<InternalData &>(*data_ptr);
  data.update_each = requires_update_flags(flags);

  // get data objects from each of the base elements and store
  // them. one might think that doing this in parallel (over the
  // base elements) would be a good idea, but this turns out to
  // be wrong because we would then run these jobs on different
  // threads/processors and this allocates memory in different
  // NUMA domains; this has large detrimental effects when later
  // writing into these objects in fill_fe_*_values. all of this
  // is particularly true when using FEValues objects in
  // WorkStream contexts where we explicitly make sure that
  // every function only uses objects previously allocated
  // in the same NUMA context and on the same thread as the
  // function is called
  for (unsigned int base_no = 0; base_no < this->n_base_elements(); ++base_no)
    {
      internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
        &base_fe_output_object = data.get_fe_output_object(base_no);
      base_fe_output_object.initialize(
        quadrature.size(),
        base_element(base_no),
        flags | base_element(base_no).requires_update_flags(flags));

      // let base objects produce their scratch objects. they may
      // also at this time write into the output objects we provide
      // for them; it would be nice if we could already copy something
      // out of the base output object into the system output object,
      // but we can't because we can't know what the elements already
      // copied and/or will want to update on every cell
      auto base_fe_data = base_element(base_no).get_data(flags,
                                                         mapping,
                                                         quadrature,
                                                         base_fe_output_object);

      data.set_fe_data(base_no, std::move(base_fe_data));
    }

  return data_ptr;
}

// The following function is a clone of get_data, with the exception
// that get_face_data of the base elements is called.

template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FESystem<dim, spacedim>::get_face_data(
  const UpdateFlags               flags,
  const Mapping<dim, spacedim> &  mapping,
  const hp::QCollection<dim - 1> &quadrature,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    & /*output_data*/) const
{
  // create an internal data object and set the update flags we will need
  // to deal with. the current object does not make use of these flags,
  // but we need to nevertheless set them correctly since we look
  // into the update_each flag of base elements in fill_fe_values,
  // and so the current object's update_each flag needs to be
  // correct in case the current FESystem is a base element for another,
  // higher-level FESystem itself.
  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
        data_ptr   = std::make_unique<InternalData>(this->n_base_elements());
  auto &data       = dynamic_cast<InternalData &>(*data_ptr);
  data.update_each = requires_update_flags(flags);

  // get data objects from each of the base elements and store
  // them. one might think that doing this in parallel (over the
  // base elements) would be a good idea, but this turns out to
  // be wrong because we would then run these jobs on different
  // threads/processors and this allocates memory in different
  // NUMA domains; this has large detrimental effects when later
  // writing into these objects in fill_fe_*_values. all of this
  // is particularly true when using FEValues objects in
  // WorkStream contexts where we explicitly make sure that
  // every function only uses objects previously allocated
  // in the same NUMA context and on the same thread as the
  // function is called
  for (unsigned int base_no = 0; base_no < this->n_base_elements(); ++base_no)
    {
      internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
        &base_fe_output_object = data.get_fe_output_object(base_no);
      base_fe_output_object.initialize(
        quadrature.max_n_quadrature_points(),
        base_element(base_no),
        flags | base_element(base_no).requires_update_flags(flags));

      // let base objects produce their scratch objects. they may
      // also at this time write into the output objects we provide
      // for them; it would be nice if we could already copy something
      // out of the base output object into the system output object,
      // but we can't because we can't know what the elements already
      // copied and/or will want to update on every cell
      auto base_fe_data = base_element(base_no).get_face_data(
        flags, mapping, quadrature, base_fe_output_object);

      data.set_fe_data(base_no, std::move(base_fe_data));
    }

  return data_ptr;
}



// The following function is a clone of get_data, with the exception
// that get_subface_data of the base elements is called.

template <int dim, int spacedim>
std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
FESystem<dim, spacedim>::get_subface_data(
  const UpdateFlags             flags,
  const Mapping<dim, spacedim> &mapping,
  const Quadrature<dim - 1> &   quadrature,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    & /*output_data*/) const
{
  // create an internal data object and set the update flags we will need
  // to deal with. the current object does not make use of these flags,
  // but we need to nevertheless set them correctly since we look
  // into the update_each flag of base elements in fill_fe_values,
  // and so the current object's update_each flag needs to be
  // correct in case the current FESystem is a base element for another,
  // higher-level FESystem itself.
  std::unique_ptr<typename FiniteElement<dim, spacedim>::InternalDataBase>
        data_ptr = std::make_unique<InternalData>(this->n_base_elements());
  auto &data     = dynamic_cast<InternalData &>(*data_ptr);

  data.update_each = requires_update_flags(flags);

  // get data objects from each of the base elements and store
  // them. one might think that doing this in parallel (over the
  // base elements) would be a good idea, but this turns out to
  // be wrong because we would then run these jobs on different
  // threads/processors and this allocates memory in different
  // NUMA domains; this has large detrimental effects when later
  // writing into these objects in fill_fe_*_values. all of this
  // is particularly true when using FEValues objects in
  // WorkStream contexts where we explicitly make sure that
  // every function only uses objects previously allocated
  // in the same NUMA context and on the same thread as the
  // function is called
  for (unsigned int base_no = 0; base_no < this->n_base_elements(); ++base_no)
    {
      internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
        &base_fe_output_object = data.get_fe_output_object(base_no);
      base_fe_output_object.initialize(
        quadrature.size(),
        base_element(base_no),
        flags | base_element(base_no).requires_update_flags(flags));

      // let base objects produce their scratch objects. they may
      // also at this time write into the output objects we provide
      // for them; it would be nice if we could already copy something
      // out of the base output object into the system output object,
      // but we can't because we can't know what the elements already
      // copied and/or will want to update on every cell
      auto base_fe_data = base_element(base_no).get_subface_data(
        flags, mapping, quadrature, base_fe_output_object);

      data.set_fe_data(base_no, std::move(base_fe_data));
    }

  return data_ptr;
}



template <int dim, int spacedim>
void
FESystem<dim, spacedim>::fill_fe_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const CellSimilarity::Similarity                            cell_similarity,
  const Quadrature<dim> &                                     quadrature,
  const Mapping<dim, spacedim> &                              mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase &   mapping_internal,
  const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                     spacedim>
    &                                                            mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  compute_fill(mapping,
               cell,
               invalid_face_number,
               invalid_face_number,
               hp::QCollection<dim>(quadrature),
               cell_similarity,
               mapping_internal,
               fe_internal,
               mapping_data,
               output_data);
}



template <int dim, int spacedim>
void
FESystem<dim, spacedim>::fill_fe_face_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const hp::QCollection<dim - 1> &                            quadrature,
  const Mapping<dim, spacedim> &                              mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase &   mapping_internal,
  const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                     spacedim>
    &                                                            mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  compute_fill(mapping,
               cell,
               face_no,
               invalid_face_number,
               quadrature,
               CellSimilarity::none,
               mapping_internal,
               fe_internal,
               mapping_data,
               output_data);
}



template <int dim, int spacedim>
void
FESystem<dim, spacedim>::fill_fe_subface_values(
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          sub_no,
  const Quadrature<dim - 1> &                                 quadrature,
  const Mapping<dim, spacedim> &                              mapping,
  const typename Mapping<dim, spacedim>::InternalDataBase &   mapping_internal,
  const dealii::internal::FEValuesImplementation::MappingRelatedData<dim,
                                                                     spacedim>
    &                                                            mapping_data,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  dealii::internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                     spacedim>
    &output_data) const
{
  compute_fill(mapping,
               cell,
               face_no,
               sub_no,
               hp::QCollection<dim - 1>(quadrature),
               CellSimilarity::none,
               mapping_internal,
               fe_internal,
               mapping_data,
               output_data);
}



template <int dim, int spacedim>
template <int dim_1>
void
FESystem<dim, spacedim>::compute_fill(
  const Mapping<dim, spacedim> &                              mapping,
  const typename Triangulation<dim, spacedim>::cell_iterator &cell,
  const unsigned int                                          face_no,
  const unsigned int                                          sub_no,
  const hp::QCollection<dim_1> &                              quadrature,
  const CellSimilarity::Similarity                            cell_similarity,
  const typename Mapping<dim, spacedim>::InternalDataBase &   mapping_internal,
  const typename FiniteElement<dim, spacedim>::InternalDataBase &fe_internal,
  const internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
    &mapping_data,
  internal::FEValuesImplementation::FiniteElementRelatedData<dim, spacedim>
    &output_data) const
{
  // convert data object to internal
  // data for this class. fails with
  // an exception if that is not
  // possible
  Assert(dynamic_cast<const InternalData *>(&fe_internal) != nullptr,
         ExcInternalError());
  const InternalData &fe_data = static_cast<const InternalData &>(fe_internal);

  // Either dim_1==dim
  // (fill_fe_values) or dim_1==dim-1
  // (fill_fe_(sub)face_values)
  Assert(dim_1 == dim || dim_1 == dim - 1, ExcInternalError());
  const UpdateFlags flags = fe_data.update_each;


  // loop over the base elements, let them compute what they need to compute,
  // and then copy what is necessary.
  //
  // one may think that it would be a good idea to parallelize this over
  // base elements, but it turns out to be not worthwhile: doing so lets
  // multiple threads access data objects that were created by the current
  // thread, leading to many NUMA memory access inefficiencies. we specifically
  // want to avoid this if this class is called in a WorkStream context where
  // we very carefully allocate objects only on the thread where they
  // will actually be used; spawning new tasks here would be counterproductive
  if (flags & (update_values | update_gradients | update_hessians |
               update_3rd_derivatives))
    for (unsigned int base_no = 0; base_no < this->n_base_elements(); ++base_no)
      {
        const FiniteElement<dim, spacedim> &base_fe = base_element(base_no);
        typename FiniteElement<dim, spacedim>::InternalDataBase &base_fe_data =
          fe_data.get_fe_data(base_no);
        internal::FEValuesImplementation::FiniteElementRelatedData<dim,
                                                                   spacedim>
          &base_data = fe_data.get_fe_output_object(base_no);

        // fill_fe_face_values needs argument Quadrature<dim-1> for both cases
        // dim_1==dim-1 and dim_1=dim. Hence the following workaround
        const Quadrature<dim> *         cell_quadrature     = nullptr;
        const hp::QCollection<dim - 1> *face_quadrature     = nullptr;
        const Quadrature<dim - 1> *     sub_face_quadrature = nullptr;
        const unsigned int              n_q_points = quadrature.size() == 1 ?
                                          quadrature[0].size() :
                                          quadrature[face_no].size();

        // static cast to the common base class of quadrature being either
        // Quadrature<dim> or Quadrature<dim-1>:

        if (face_no == invalid_face_number)
          {
            const Subscriptor *quadrature_base_pointer = &quadrature[0];
            Assert(dim_1 == dim, ExcDimensionMismatch(dim_1, dim));
            Assert(dynamic_cast<const Quadrature<dim> *>(
                     quadrature_base_pointer) != nullptr,
                   ExcInternalError());

            cell_quadrature =
              static_cast<const Quadrature<dim> *>(quadrature_base_pointer);
          }
        else if (sub_no == invalid_face_number)
          {
            const Subscriptor *quadrature_base_pointer = &quadrature;
            Assert(dim_1 == dim - 1, ExcDimensionMismatch(dim_1, dim - 1));
            Assert(dynamic_cast<const hp::QCollection<dim - 1> *>(
                     quadrature_base_pointer) != nullptr,
                   ExcInternalError());

            face_quadrature = static_cast<const hp::QCollection<dim - 1> *>(
              quadrature_base_pointer);
          }
        else
          {
            const Subscriptor *quadrature_base_pointer = &quadrature[0];
            Assert(dim_1 == dim - 1, ExcDimensionMismatch(dim_1, dim - 1));
            Assert(dynamic_cast<const Quadrature<dim - 1> *>(
                     quadrature_base_pointer) != nullptr,
                   ExcInternalError());

            sub_face_quadrature =
              static_cast<const Quadrature<dim - 1> *>(quadrature_base_pointer);
          }


        // Make sure that in the case of fill_fe_values the data is only
        // copied from base_data to data if base_data is changed. therefore
        // use fe_fe_data.current_update_flags()
        //
        // for the case of fill_fe_(sub)face_values the data needs to be
        // copied from base_data to data on each face, therefore use
        // base_fe_data.update_flags.
        if (face_no == invalid_face_number)
          base_fe.fill_fe_values(cell,
                                 cell_similarity,
                                 *cell_quadrature,
                                 mapping,
                                 mapping_internal,
                                 mapping_data,
                                 base_fe_data,
                                 base_data);
        else if (sub_no == invalid_face_number)
          base_fe.fill_fe_face_values(cell,
                                      face_no,
                                      *face_quadrature,
                                      mapping,
                                      mapping_internal,
                                      mapping_data,
                                      base_fe_data,
                                      base_data);
        else
          base_fe.fill_fe_subface_values(cell,
                                         face_no,
                                         sub_no,
                                         *sub_face_quadrature,
                                         mapping,
                                         mapping_internal,
                                         mapping_data,
                                         base_fe_data,
                                         base_data);

        // now data has been generated, so copy it. we used to work by
        // looping over all base elements (i.e. this outer loop), then over
        // multiplicity, then over the shape functions from that base
        // element, but that requires that we can infer the global number of
        // a shape function from its number in the base element. for that we
        // had the component_to_system_table.
        //
        // however, this does of course no longer work since we have
        // non-primitive elements. so we go the other way round: loop over
        // all shape functions of the composed element, and here only treat
        // those shape functions that belong to a given base element
        // TODO: Introduce the needed table and loop only over base element
        // shape functions. This here is not efficient at all AND very bad style
        const UpdateFlags base_flags = base_fe_data.update_each;

        // some base element might involve values that depend on the shape
        // of the geometry, so we always need to copy the shape values around
        // also in case we detected a cell similarity (but no heavy work will
        // be done inside the individual elements in case we have a
        // translation and simple elements).
        for (unsigned int system_index = 0;
             system_index < this->n_dofs_per_cell();
             ++system_index)
          if (this->system_to_base_table[system_index].first.first == base_no)
            {
              const unsigned int base_index =
                this->system_to_base_table[system_index].second;
              Assert(base_index < base_fe.n_dofs_per_cell(),
                     ExcInternalError());

              // now copy. if the shape function is primitive, then there
              // is only one value to be copied, but for non-primitive
              // elements, there might be more values to be copied
              //
              // so, find out from which index to take this one value, and
              // to which index to put
              unsigned int out_index = 0;
              for (unsigned int i = 0; i < system_index; ++i)
                out_index += this->n_nonzero_components(i);
              unsigned int in_index = 0;
              for (unsigned int i = 0; i < base_index; ++i)
                in_index += base_fe.n_nonzero_components(i);

              // then loop over the number of components to be copied
              Assert(this->n_nonzero_components(system_index) ==
                       base_fe.n_nonzero_components(base_index),
                     ExcInternalError());

              if (base_flags & update_values)
                for (unsigned int s = 0;
                     s < this->n_nonzero_components(system_index);
                     ++s)
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    output_data.shape_values[out_index + s][q] =
                      base_data.shape_values(in_index + s, q);

              if (base_flags & update_gradients)
                for (unsigned int s = 0;
                     s < this->n_nonzero_components(system_index);
                     ++s)
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    output_data.shape_gradients[out_index + s][q] =
                      base_data.shape_gradients[in_index + s][q];

              if (base_flags & update_hessians)
                for (unsigned int s = 0;
                     s < this->n_nonzero_components(system_index);
                     ++s)
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    output_data.shape_hessians[out_index + s][q] =
                      base_data.shape_hessians[in_index + s][q];

              if (base_flags & update_3rd_derivatives)
                for (unsigned int s = 0;
                     s < this->n_nonzero_components(system_index);
                     ++s)
                  for (unsigned int q = 0; q < n_q_points; ++q)
                    output_data.shape_3rd_derivatives[out_index + s][q] =
                      base_data.shape_3rd_derivatives[in_index + s][q];
            }
      }
}


template <int dim, int spacedim>
void
FESystem<dim, spacedim>::build_interface_constraints()
{
  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  // check whether all base elements implement their interface constraint
  // matrices. if this is not the case, then leave the interface costraints of
  // this composed element empty as well; however, the rest of the element is
  // usable
  for (unsigned int base = 0; base < this->n_base_elements(); ++base)
    if (base_element(base).constraints_are_implemented() == false)
      return;

  this->interface_constraints.TableBase<2, double>::reinit(
    this->interface_constraints_size());

  // the layout of the constraints matrix is described in the FiniteElement
  // class. you may want to look there first before trying to understand the
  // following, especially the mapping of the @p{m} index.
  //
  // in order to map it to the fe-system class, we have to know which base
  // element a degree of freedom within a vertex, line, etc belongs to. this
  // can be accomplished by the system_to_component_index function in
  // conjunction with the numbers first_{line,quad,...}_index
  for (unsigned int n = 0; n < this->interface_constraints.n(); ++n)
    for (unsigned int m = 0; m < this->interface_constraints.m(); ++m)
      {
        // for the pair (n,m) find out which base element they belong to and
        // the number therein
        //
        // first for the n index. this is simple since the n indices are in
        // the same order as they are usually on a face. note that for the
        // data type, first value in pair is (base element,instance of base
        // element), second is index within this instance
        const std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
          n_index = this->face_system_to_base_table[face_no][n];

        // likewise for the m index. this is more complicated due to the
        // strange ordering we have for the dofs on the refined faces.
        std::pair<std::pair<unsigned int, unsigned int>, unsigned int> m_index;
        switch (dim)
          {
            case 1:
              {
                // we should never get here!  (in 1d, the constraints matrix
                // should be of size zero)
                Assert(false, ExcInternalError());
                break;
              }

            case 2:
              {
                // the indices m=0..d_v-1 are from the center vertex.  their
                // order is the same as for the first vertex of the whole cell,
                // so we can use the system_to_base_table variable (using the
                // face_s_t_base_t function would yield the same)
                if (m < this->n_dofs_per_vertex())
                  m_index = this->system_to_base_table[m];
                else
                  // then come the two sets of line indices
                  {
                    const unsigned int index_in_line =
                      (m - this->n_dofs_per_vertex()) % this->n_dofs_per_line();
                    const unsigned int sub_line =
                      (m - this->n_dofs_per_vertex()) / this->n_dofs_per_line();
                    Assert(sub_line < 2, ExcInternalError());

                    // from this information, try to get base element and
                    // instance of base element. we do so by constructing the
                    // corresponding face index of m in the present element,
                    // then use face_system_to_base_table
                    const unsigned int tmp1 =
                      2 * this->n_dofs_per_vertex() + index_in_line;
                    m_index.first =
                      this->face_system_to_base_table[face_no][tmp1].first;

                    // what we are still missing is the index of m within the
                    // base elements interface_constraints table
                    //
                    // here, the second value of face_system_to_base_table can
                    // help: it denotes the face index of that shape function
                    // within the base element. since we know that it is a line
                    // dof, we can construct the rest: tmp2 will denote the
                    // index of this shape function among the line shape
                    // functions:
                    Assert(
                      this->face_system_to_base_table[face_no][tmp1].second >=
                        2 *
                          base_element(m_index.first.first).n_dofs_per_vertex(),
                      ExcInternalError());
                    const unsigned int tmp2 =
                      this->face_system_to_base_table[face_no][tmp1].second -
                      2 * base_element(m_index.first.first).n_dofs_per_vertex();
                    Assert(tmp2 < base_element(m_index.first.first)
                                    .n_dofs_per_line(),
                           ExcInternalError());
                    m_index.second =
                      base_element(m_index.first.first).n_dofs_per_vertex() +
                      base_element(m_index.first.first).n_dofs_per_line() *
                        sub_line +
                      tmp2;
                  }
                break;
              }

            case 3:
              {
                // same way as above, although a little more complicated...

                // the indices m=0..5*d_v-1 are from the center and the four
                // subline vertices.  their order is the same as for the first
                // vertex of the whole cell, so we can use the simple arithmetic
                if (m < 5 * this->n_dofs_per_vertex())
                  m_index = this->system_to_base_table[m];
                else
                  // then come the 12 sets of line indices
                  if (m < 5 * this->n_dofs_per_vertex() +
                            12 * this->n_dofs_per_line())
                  {
                    // for the meaning of all this, see the 2d part
                    const unsigned int index_in_line =
                      (m - 5 * this->n_dofs_per_vertex()) %
                      this->n_dofs_per_line();
                    const unsigned int sub_line =
                      (m - 5 * this->n_dofs_per_vertex()) /
                      this->n_dofs_per_line();
                    Assert(sub_line < 12, ExcInternalError());

                    const unsigned int tmp1 =
                      4 * this->n_dofs_per_vertex() + index_in_line;
                    m_index.first =
                      this->face_system_to_base_table[face_no][tmp1].first;

                    Assert(
                      this->face_system_to_base_table[face_no][tmp1].second >=
                        4 *
                          base_element(m_index.first.first).n_dofs_per_vertex(),
                      ExcInternalError());
                    const unsigned int tmp2 =
                      this->face_system_to_base_table[face_no][tmp1].second -
                      4 * base_element(m_index.first.first).n_dofs_per_vertex();
                    Assert(tmp2 < base_element(m_index.first.first)
                                    .n_dofs_per_line(),
                           ExcInternalError());
                    m_index.second =
                      5 *
                        base_element(m_index.first.first).n_dofs_per_vertex() +
                      base_element(m_index.first.first).n_dofs_per_line() *
                        sub_line +
                      tmp2;
                  }
                else
                  // on one of the four sub-quads
                  {
                    // for the meaning of all this, see the 2d part
                    const unsigned int index_in_quad =
                      (m - 5 * this->n_dofs_per_vertex() -
                       12 * this->n_dofs_per_line()) %
                      this->n_dofs_per_quad(face_no);
                    Assert(index_in_quad < this->n_dofs_per_quad(face_no),
                           ExcInternalError());
                    const unsigned int sub_quad =
                      ((m - 5 * this->n_dofs_per_vertex() -
                        12 * this->n_dofs_per_line()) /
                       this->n_dofs_per_quad(face_no));
                    Assert(sub_quad < 4, ExcInternalError());

                    const unsigned int tmp1 = 4 * this->n_dofs_per_vertex() +
                                              4 * this->n_dofs_per_line() +
                                              index_in_quad;
                    Assert(tmp1 <
                             this->face_system_to_base_table[face_no].size(),
                           ExcInternalError());
                    m_index.first =
                      this->face_system_to_base_table[face_no][tmp1].first;

                    Assert(
                      this->face_system_to_base_table[face_no][tmp1].second >=
                        4 * base_element(m_index.first.first)
                              .n_dofs_per_vertex() +
                          4 *
                            base_element(m_index.first.first).n_dofs_per_line(),
                      ExcInternalError());
                    const unsigned int tmp2 =
                      this->face_system_to_base_table[face_no][tmp1].second -
                      4 *
                        base_element(m_index.first.first).n_dofs_per_vertex() -
                      4 * base_element(m_index.first.first).n_dofs_per_line();
                    Assert(tmp2 < base_element(m_index.first.first)
                                    .n_dofs_per_quad(face_no),
                           ExcInternalError());
                    m_index.second =
                      5 *
                        base_element(m_index.first.first).n_dofs_per_vertex() +
                      12 * base_element(m_index.first.first).n_dofs_per_line() +
                      base_element(m_index.first.first)
                          .n_dofs_per_quad(face_no) *
                        sub_quad +
                      tmp2;
                  }

                break;
              }

            default:
              Assert(false, ExcNotImplemented());
          }

        // now that we gathered all information: use it to build the
        // matrix. note that if n and m belong to different base elements or
        // instances, then there definitely will be no coupling
        if (n_index.first == m_index.first)
          this->interface_constraints(m, n) =
            (base_element(n_index.first.first)
               .constraints()(m_index.second, n_index.second));
      }
}



template <int dim, int spacedim>
void
FESystem<dim, spacedim>::initialize(
  const std::vector<const FiniteElement<dim, spacedim> *> &fes,
  const std::vector<unsigned int> &                        multiplicities)
{
  Assert(fes.size() == multiplicities.size(),
         ExcDimensionMismatch(fes.size(), multiplicities.size()));
  Assert(fes.size() > 0,
         ExcMessage("Need to pass at least one finite element."));
  Assert(count_nonzeros(multiplicities) > 0,
         ExcMessage("You only passed FiniteElements with multiplicity 0."));

  // Note that we need to skip every FE with multiplicity 0 in the following
  // block of code

  this->base_to_block_indices.reinit(0, 0);

  for (unsigned int i = 0; i < fes.size(); i++)
    if (multiplicities[i] > 0)
      this->base_to_block_indices.push_back(multiplicities[i]);

  {
    Threads::TaskGroup<> clone_base_elements;

    unsigned int ind = 0;
    for (unsigned int i = 0; i < fes.size(); i++)
      if (multiplicities[i] > 0)
        {
          clone_base_elements += Threads::new_task([&, i, ind]() {
            base_elements[ind] = {fes[i]->clone(), multiplicities[i]};
          });
          ++ind;
        }
    Assert(ind > 0, ExcInternalError());

    // wait for all of these clone operations to finish
    clone_base_elements.join_all();
  }


  {
    // If the system is not primitive, these have not been initialized by
    // FiniteElement
    this->system_to_component_table.resize(this->n_dofs_per_cell());

    FETools::Compositing::build_cell_tables(this->system_to_base_table,
                                            this->system_to_component_table,
                                            this->component_to_base_table,
                                            *this);

    this->face_system_to_component_table.resize(this->n_unique_faces());

    for (unsigned int face_no = 0; face_no < this->n_unique_faces(); ++face_no)
      {
        this->face_system_to_component_table[0].resize(
          this->n_dofs_per_face(face_no));

        FETools::Compositing::build_face_tables(
          this->face_system_to_base_table[face_no],
          this->face_system_to_component_table[face_no],
          *this,
          true,
          face_no);
      }
  }

  // now initialize interface constraints, support points, and other tables.
  // (restriction and prolongation matrices are only built on demand.) do
  // this in parallel

  Threads::TaskGroup<> init_tasks;

  init_tasks +=
    Threads::new_task([&]() { this->build_interface_constraints(); });

  init_tasks += Threads::new_task([&]() {
    // if one of the base elements has no support points, then it makes no sense
    // to define support points for the composed element, so return an empty
    // array to demonstrate that fact. Note that we ignore FE_Nothing in this
    // logic.
    for (unsigned int base_el = 0; base_el < this->n_base_elements(); ++base_el)
      if (!base_element(base_el).has_support_points() &&
          base_element(base_el).n_dofs_per_cell() != 0)
        {
          this->unit_support_points.resize(0);
          return;
        }

    // generate unit support points from unit support points of sub elements
    this->unit_support_points.resize(this->n_dofs_per_cell());

    for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
      {
        const unsigned int base = this->system_to_base_table[i].first.first,
                           base_index = this->system_to_base_table[i].second;
        Assert(base < this->n_base_elements(), ExcInternalError());
        Assert(base_index < base_element(base).unit_support_points.size(),
               ExcInternalError());
        this->unit_support_points[i] =
          base_element(base).unit_support_points[base_index];
      }
  });

  // initialize face support points (for dim==2,3). same procedure as above
  if (dim > 1)
    init_tasks += Threads::new_task([&]() {
      for (unsigned int face_no = 0; face_no < this->n_unique_faces();
           ++face_no)
        {
          // if one of the base elements has no support points, then it makes
          // no sense to define support points for the composed element. In
          // that case, return an empty array to demonstrate that fact (note
          // that we ask whether the base element has no support points at
          // all, not only none on the face!)
          //
          // on the other hand, if there is an element that simply has no
          // degrees of freedom on the face at all, then we don't care whether
          // it has support points or not. this is, for example, the case for
          // the stable Stokes element Q(p)^dim \times DGP(p-1).
          bool flag_has_no_support_points = false;

          for (unsigned int base_el = 0; base_el < this->n_base_elements();
               ++base_el)
            if (!base_element(base_el).has_support_points() &&
                (base_element(base_el).n_dofs_per_face(face_no) > 0))
              {
                this->unit_face_support_points[face_no].resize(0);
                flag_has_no_support_points = true;
                break;
              }


          if (flag_has_no_support_points)
            continue;

          // generate unit face support points from unit support points of sub
          // elements
          this->unit_face_support_points[face_no].resize(
            this->n_dofs_per_face(face_no));

          for (unsigned int i = 0; i < this->n_dofs_per_face(face_no); ++i)
            {
              const unsigned int base_i =
                this->face_system_to_base_table[face_no][i].first.first;
              const unsigned int index_in_base =
                this->face_system_to_base_table[face_no][i].second;

              Assert(
                index_in_base <
                  base_element(base_i).unit_face_support_points[face_no].size(),
                ExcInternalError());

              this->unit_face_support_points[face_no][i] =
                base_element(base_i)
                  .unit_face_support_points[face_no][index_in_base];
            }
        }
    });

  // Initialize generalized support points and an (internal) index table
  init_tasks += Threads::new_task([&]() {
    // Iterate over all base elements, extract a representative set of
    // _unique_ generalized support points and store the information how
    // generalized support points of base elements are mapped to this list
    // of representatives. Complexity O(n^2), where n is the number of
    // generalized support points.

    generalized_support_points_index_table.resize(this->n_base_elements());

    for (unsigned int base = 0; base < this->n_base_elements(); ++base)
      {
        // If the current base element does not have generalized support
        // points, ignore it. Note that
        // * FESystem::convert_generalized_support_point_values_to_dof_values
        //   will simply skip such non-interpolatory base elements by
        //   assigning NaN to all dofs.
        // * If this routine does not pick up any generalized support
        //   points the corresponding vector will be empty and
        //   FiniteElement::has_generalized_support_points will return
        //   false.
        if (!base_element(base).has_generalized_support_points())
          continue;

        for (const auto &point :
             base_element(base).get_generalized_support_points())
          {
            // Is point already an element of generalized_support_points?
            const auto p =
              std::find(std::begin(this->generalized_support_points),
                        std::end(this->generalized_support_points),
                        point);

            if (p == std::end(this->generalized_support_points))
              {
                // If no, update the table and add the point to the vector
                const auto n = this->generalized_support_points.size();
                generalized_support_points_index_table[base].push_back(n);
                this->generalized_support_points.push_back(point);
              }
            else
              {
                // If yes, just add the correct index to the table.
                const auto n = p - std::begin(this->generalized_support_points);
                generalized_support_points_index_table[base].push_back(n);
              }
          }
      }

#ifdef DEBUG
    // check generalized_support_points_index_table for consistency
    for (unsigned int i = 0; i < base_elements.size(); ++i)
      {
        if (!base_element(i).has_generalized_support_points())
          continue;

        const auto &points =
          base_elements[i].first->get_generalized_support_points();
        for (unsigned int j = 0; j < points.size(); ++j)
          {
            const auto n = generalized_support_points_index_table[i][j];
            Assert(this->generalized_support_points[n] == points[j],
                   ExcInternalError());
          }
      }
#endif /* DEBUG */
  });

  // initialize quad dof index permutation in 3d and higher
  if (dim >= 3)
    init_tasks += Threads::new_task([&]() {
      for (unsigned int face_no = 0; face_no < this->n_unique_faces();
           ++face_no)
        {
          // the array into which we want to write should have the correct size
          // already.
          Assert(this->adjust_quad_dof_index_for_face_orientation_table[face_no]
                     .n_elements() == 8 * this->n_dofs_per_quad(face_no),
                 ExcInternalError());

          // to obtain the shifts for this composed element, copy the shift
          // information of the base elements
          unsigned int index = 0;
          for (unsigned int b = 0; b < this->n_base_elements(); ++b)
            {
              const Table<2, int> &temp =
                this->base_element(b)
                  .adjust_quad_dof_index_for_face_orientation_table[face_no];
              for (unsigned int c = 0; c < this->element_multiplicity(b); ++c)
                {
                  for (unsigned int i = 0; i < temp.size(0); ++i)
                    for (unsigned int j = 0; j < 8; ++j)
                      this->adjust_quad_dof_index_for_face_orientation_table
                        [face_no](index + i, j) = temp(i, j);
                  index += temp.size(0);
                }
            }
          Assert(index == this->n_dofs_per_quad(face_no), ExcInternalError());
        }

      // additionally compose the permutation information for lines
      Assert(this->adjust_line_dof_index_for_line_orientation_table.size() ==
               this->n_dofs_per_line(),
             ExcInternalError());
      unsigned int index = 0;
      for (unsigned int b = 0; b < this->n_base_elements(); ++b)
        {
          const std::vector<int> &temp2 =
            this->base_element(b)
              .adjust_line_dof_index_for_line_orientation_table;
          for (unsigned int c = 0; c < this->element_multiplicity(b); ++c)
            {
              std::copy(
                temp2.begin(),
                temp2.end(),
                this->adjust_line_dof_index_for_line_orientation_table.begin() +
                  index);
              index += temp2.size();
            }
        }
      Assert(index == this->n_dofs_per_line(), ExcInternalError());
    });

  // wait for all of this to finish
  init_tasks.join_all();
}



template <int dim, int spacedim>
bool
FESystem<dim, spacedim>::hp_constraints_are_implemented() const
{
  for (unsigned int b = 0; b < this->n_base_elements(); ++b)
    if (base_element(b).hp_constraints_are_implemented() == false)
      return false;

  return true;
}



template <int dim, int spacedim>
void
FESystem<dim, spacedim>::get_face_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  FullMatrix<double> &                interpolation_matrix,
  const unsigned int                  face_no) const
{
  Assert(interpolation_matrix.n() == this->n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.n(),
                              this->n_dofs_per_face(face_no)));
  Assert(interpolation_matrix.m() == x_source_fe.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_face(face_no)));

  // since dofs for each base are independent, we only have to stack things up
  // from base element to base element
  //
  // the problem is that we have to work with two FEs (this and
  // fe_other). only deal with the case that both are FESystems and that they
  // both have the same number of bases (counting multiplicity) each of which
  // match in their number of components. this covers
  // FESystem(FE_Q(p),1,FE_Q(q),2) vs FESystem(FE_Q(r),2,FE_Q(s),1), but not
  // FESystem(FE_Q(p),1,FE_Q(q),2) vs
  // FESystem(FESystem(FE_Q(r),2),1,FE_Q(s),1)
  if (const auto *fe_other_system =
        dynamic_cast<const FESystem<dim, spacedim> *>(&x_source_fe))
    {
      // clear matrix, since we will not get to set all elements
      interpolation_matrix = 0;

      // loop over all the base elements of this and the other element, counting
      // their multiplicities
      unsigned int base_index = 0, base_index_other = 0;
      unsigned int multiplicity = 0, multiplicity_other = 0;

      FullMatrix<double> base_to_base_interpolation;

      while (true)
        {
          const FiniteElement<dim, spacedim> &base = base_element(base_index),
                                             &base_other =
                                               fe_other_system->base_element(
                                                 base_index_other);

          Assert(base.n_components() == base_other.n_components(),
                 ExcNotImplemented());

          // get the interpolation from the bases
          base_to_base_interpolation.reinit(base_other.n_dofs_per_face(face_no),
                                            base.n_dofs_per_face(face_no));
          base.get_face_interpolation_matrix(base_other,
                                             base_to_base_interpolation,
                                             face_no);

          // now translate entries. we'd like to have something like
          // face_base_to_system_index, but that doesn't exist. rather, all we
          // have is the reverse. well, use that then
          for (unsigned int i = 0; i < this->n_dofs_per_face(face_no); ++i)
            if (this->face_system_to_base_index(i, face_no).first ==
                std::make_pair(base_index, multiplicity))
              for (unsigned int j = 0;
                   j < fe_other_system->n_dofs_per_face(face_no);
                   ++j)
                if (fe_other_system->face_system_to_base_index(j, face_no)
                      .first ==
                    std::make_pair(base_index_other, multiplicity_other))
                  interpolation_matrix(j, i) = base_to_base_interpolation(
                    fe_other_system->face_system_to_base_index(j, face_no)
                      .second,
                    this->face_system_to_base_index(i, face_no).second);

          // advance to the next base element for this and the other fe_system;
          // see if we can simply advance the multiplicity by one, or if have to
          // move on to the next base element
          ++multiplicity;
          if (multiplicity == this->element_multiplicity(base_index))
            {
              multiplicity = 0;
              ++base_index;
            }
          ++multiplicity_other;
          if (multiplicity_other ==
              fe_other_system->element_multiplicity(base_index_other))
            {
              multiplicity_other = 0;
              ++base_index_other;
            }

          // see if we have reached the end of the present element. if so, we
          // should have reached the end of the other one as well
          if (base_index == this->n_base_elements())
            {
              Assert(base_index_other == fe_other_system->n_base_elements(),
                     ExcInternalError());
              break;
            }

          // if we haven't reached the end of this element, we shouldn't have
          // reached the end of the other one either
          Assert(base_index_other != fe_other_system->n_base_elements(),
                 ExcInternalError());
        }
    }
  else
    {
      // repeat the cast to make the exception message more useful
      AssertThrow(
        (dynamic_cast<const FESystem<dim, spacedim> *>(&x_source_fe) !=
         nullptr),
        (typename FiniteElement<dim,
                                spacedim>::ExcInterpolationNotImplemented()));
    }
}



template <int dim, int spacedim>
void
FESystem<dim, spacedim>::get_subface_interpolation_matrix(
  const FiniteElement<dim, spacedim> &x_source_fe,
  const unsigned int                  subface,
  FullMatrix<double> &                interpolation_matrix,
  const unsigned int                  face_no) const
{
  AssertThrow(
    (x_source_fe.get_name().find("FESystem<") == 0) ||
      (dynamic_cast<const FESystem<dim, spacedim> *>(&x_source_fe) != nullptr),
    (typename FiniteElement<dim, spacedim>::ExcInterpolationNotImplemented()));

  Assert(interpolation_matrix.n() == this->n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.n(),
                              this->n_dofs_per_face(face_no)));
  Assert(interpolation_matrix.m() == x_source_fe.n_dofs_per_face(face_no),
         ExcDimensionMismatch(interpolation_matrix.m(),
                              x_source_fe.n_dofs_per_face(face_no)));

  // since dofs for each base are independent, we only have to stack things up
  // from base element to base element
  //
  // the problem is that we have to work with two FEs (this and
  // fe_other). only deal with the case that both are FESystems and that they
  // both have the same number of bases (counting multiplicity) each of which
  // match in their number of components. this covers
  // FESystem(FE_Q(p),1,FE_Q(q),2) vs FESystem(FE_Q(r),2,FE_Q(s),1), but not
  // FESystem(FE_Q(p),1,FE_Q(q),2) vs
  // FESystem(FESystem(FE_Q(r),2),1,FE_Q(s),1)
  const FESystem<dim, spacedim> *fe_other_system =
    dynamic_cast<const FESystem<dim, spacedim> *>(&x_source_fe);
  if (fe_other_system != nullptr)
    {
      // clear matrix, since we will not get to set all elements
      interpolation_matrix = 0;

      // loop over all the base elements of this and the other element, counting
      // their multiplicities
      unsigned int base_index = 0, base_index_other = 0;
      unsigned int multiplicity = 0, multiplicity_other = 0;

      FullMatrix<double> base_to_base_interpolation;

      while (true)
        {
          const FiniteElement<dim, spacedim> &base = base_element(base_index),
                                             &base_other =
                                               fe_other_system->base_element(
                                                 base_index_other);

          Assert(base.n_components() == base_other.n_components(),
                 ExcNotImplemented());

          // get the interpolation from the bases
          base_to_base_interpolation.reinit(base_other.n_dofs_per_face(face_no),
                                            base.n_dofs_per_face(face_no));
          base.get_subface_interpolation_matrix(base_other,
                                                subface,
                                                base_to_base_interpolation,
                                                face_no);

          // now translate entries. we'd like to have something like
          // face_base_to_system_index, but that doesn't exist. rather, all we
          // have is the reverse. well, use that then
          for (unsigned int i = 0; i < this->n_dofs_per_face(face_no); ++i)
            if (this->face_system_to_base_index(i, face_no).first ==
                std::make_pair(base_index, multiplicity))
              for (unsigned int j = 0;
                   j < fe_other_system->n_dofs_per_face(face_no);
                   ++j)
                if (fe_other_system->face_system_to_base_index(j, face_no)
                      .first ==
                    std::make_pair(base_index_other, multiplicity_other))
                  interpolation_matrix(j, i) = base_to_base_interpolation(
                    fe_other_system->face_system_to_base_index(j, face_no)
                      .second,
                    this->face_system_to_base_index(i, face_no).second);

          // advance to the next base element for this and the other fe_system;
          // see if we can simply advance the multiplicity by one, or if have to
          // move on to the next base element
          ++multiplicity;
          if (multiplicity == this->element_multiplicity(base_index))
            {
              multiplicity = 0;
              ++base_index;
            }
          ++multiplicity_other;
          if (multiplicity_other ==
              fe_other_system->element_multiplicity(base_index_other))
            {
              multiplicity_other = 0;
              ++base_index_other;
            }

          // see if we have reached the end of the present element. if so, we
          // should have reached the end of the other one as well
          if (base_index == this->n_base_elements())
            {
              Assert(base_index_other == fe_other_system->n_base_elements(),
                     ExcInternalError());
              break;
            }

          // if we haven't reached the end of this element, we shouldn't have
          // reached the end of the other one either
          Assert(base_index_other != fe_other_system->n_base_elements(),
                 ExcInternalError());
        }
    }
  else
    {
      // we should have caught this at the start, but check again anyway
      Assert(
        fe_other_system != nullptr,
        (typename FiniteElement<dim,
                                spacedim>::ExcInterpolationNotImplemented()));
    }
}



template <int dim, int spacedim>
template <int structdim>
std::vector<std::pair<unsigned int, unsigned int>>
FESystem<dim, spacedim>::hp_object_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  face_no) const
{
  // since dofs on each subobject (vertex, line, ...) are ordered such that
  // first come all from the first base element all multiplicities, then
  // second base element all multiplicities, etc., we simply have to stack all
  // the identities after each other
  //
  // the problem is that we have to work with two FEs (this and
  // fe_other). only deal with the case that both are FESystems and that they
  // both have the same number of bases (counting multiplicity) each of which
  // match in their number of components. this covers
  // FESystem(FE_Q(p),1,FE_Q(q),2) vs FESystem(FE_Q(r),2,FE_Q(s),1), but not
  // FESystem(FE_Q(p),1,FE_Q(q),2) vs
  // FESystem(FESystem(FE_Q(r),2),1,FE_Q(s),1)
  if (const FESystem<dim, spacedim> *fe_other_system =
        dynamic_cast<const FESystem<dim, spacedim> *>(&fe_other))
    {
      // loop over all the base elements of this and the other element,
      // counting their multiplicities
      unsigned int base_index = 0, base_index_other = 0;
      unsigned int multiplicity = 0, multiplicity_other = 0;

      // we also need to keep track of the number of dofs already treated for
      // each of the elements
      unsigned int dof_offset = 0, dof_offset_other = 0;

      std::vector<std::pair<unsigned int, unsigned int>> identities;

      while (true)
        {
          const FiniteElement<dim, spacedim> &base = base_element(base_index),
                                             &base_other =
                                               fe_other_system->base_element(
                                                 base_index_other);

          Assert(base.n_components() == base_other.n_components(),
                 ExcNotImplemented());

          // now translate the identities returned by the base elements to the
          // indices of this system element
          std::vector<std::pair<unsigned int, unsigned int>> base_identities;
          switch (structdim)
            {
              case 0:
                base_identities = base.hp_vertex_dof_identities(base_other);
                break;
              case 1:
                base_identities = base.hp_line_dof_identities(base_other);
                break;
              case 2:
                base_identities =
                  base.hp_quad_dof_identities(base_other, face_no);
                break;
              default:
                Assert(false, ExcNotImplemented());
            }

          for (const auto &base_identity : base_identities)
            identities.emplace_back(base_identity.first + dof_offset,
                                    base_identity.second + dof_offset_other);

          // record the dofs treated above as already taken care of
          dof_offset += base.template n_dofs_per_object<structdim>();
          dof_offset_other +=
            base_other.template n_dofs_per_object<structdim>();

          // advance to the next base element for this and the other
          // fe_system; see if we can simply advance the multiplicity by one,
          // or if have to move on to the next base element
          ++multiplicity;
          if (multiplicity == this->element_multiplicity(base_index))
            {
              multiplicity = 0;
              ++base_index;
            }
          ++multiplicity_other;
          if (multiplicity_other ==
              fe_other_system->element_multiplicity(base_index_other))
            {
              multiplicity_other = 0;
              ++base_index_other;
            }

          // see if we have reached the end of the present element. if so, we
          // should have reached the end of the other one as well
          if (base_index == this->n_base_elements())
            {
              Assert(base_index_other == fe_other_system->n_base_elements(),
                     ExcInternalError());
              break;
            }

          // if we haven't reached the end of this element, we shouldn't have
          // reached the end of the other one either
          Assert(base_index_other != fe_other_system->n_base_elements(),
                 ExcInternalError());
        }

      return identities;
    }
  else
    {
      Assert(false, ExcNotImplemented());
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FESystem<dim, spacedim>::hp_vertex_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  return hp_object_dof_identities<0>(fe_other);
}

template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FESystem<dim, spacedim>::hp_line_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other) const
{
  return hp_object_dof_identities<1>(fe_other);
}



template <int dim, int spacedim>
std::vector<std::pair<unsigned int, unsigned int>>
FESystem<dim, spacedim>::hp_quad_dof_identities(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  face_no) const
{
  return hp_object_dof_identities<2>(fe_other, face_no);
}



template <int dim, int spacedim>
FiniteElementDomination::Domination
FESystem<dim, spacedim>::compare_for_domination(
  const FiniteElement<dim, spacedim> &fe_other,
  const unsigned int                  codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face/cell domination
  // --------------------------------
  // at present all we can do is to compare with other FESystems that have the
  // same number of components and bases
  if (const FESystem<dim, spacedim> *fe_sys_other =
        dynamic_cast<const FESystem<dim, spacedim> *>(&fe_other))
    {
      Assert(this->n_components() == fe_sys_other->n_components(),
             ExcNotImplemented());
      Assert(this->n_base_elements() == fe_sys_other->n_base_elements(),
             ExcNotImplemented());

      FiniteElementDomination::Domination domination =
        FiniteElementDomination::no_requirements;

      // loop over all base elements and do some sanity checks
      for (unsigned int b = 0; b < this->n_base_elements(); ++b)
        {
          Assert(this->base_element(b).n_components() ==
                   fe_sys_other->base_element(b).n_components(),
                 ExcNotImplemented());
          Assert(this->element_multiplicity(b) ==
                   fe_sys_other->element_multiplicity(b),
                 ExcNotImplemented());

          // for this pair of base elements, check who dominates and combine
          // with previous result
          const FiniteElementDomination::Domination base_domination =
            (this->base_element(b).compare_for_domination(
              fe_sys_other->base_element(b), codim));
          domination = domination & base_domination;
        }

      return domination;
    }

  Assert(false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}



template <int dim, int spacedim>
const FiniteElement<dim, spacedim> &
FESystem<dim, spacedim>::base_element(const unsigned int index) const
{
  AssertIndexRange(index, base_elements.size());
  return *base_elements[index].first;
}



template <int dim, int spacedim>
bool
FESystem<dim, spacedim>::has_support_on_face(
  const unsigned int shape_index,
  const unsigned int face_index) const
{
  return (base_element(this->system_to_base_index(shape_index).first.first)
            .has_support_on_face(this->system_to_base_index(shape_index).second,
                                 face_index));
}



template <int dim, int spacedim>
Point<dim>
FESystem<dim, spacedim>::unit_support_point(const unsigned int index) const
{
  AssertIndexRange(index, this->n_dofs_per_cell());
  Assert((this->unit_support_points.size() == this->n_dofs_per_cell()) ||
           (this->unit_support_points.size() == 0),
         (typename FiniteElement<dim, spacedim>::ExcFEHasNoSupportPoints()));

  // let's see whether we have the information pre-computed
  if (this->unit_support_points.size() != 0)
    return this->unit_support_points[index];
  else
    // no. ask the base element whether it would like to provide this
    // information
    return (base_element(this->system_to_base_index(index).first.first)
              .unit_support_point(this->system_to_base_index(index).second));
}



template <int dim, int spacedim>
Point<dim - 1>
FESystem<dim, spacedim>::unit_face_support_point(
  const unsigned int index,
  const unsigned int face_no) const
{
  AssertIndexRange(index, this->n_dofs_per_face(face_no));
  Assert(
    (this->unit_face_support_points[this->n_unique_faces() == 1 ? 0 : face_no]
       .size() == this->n_dofs_per_face(face_no)) ||
      (this->unit_face_support_points[this->n_unique_faces() == 1 ? 0 : face_no]
         .size() == 0),
    (typename FiniteElement<dim, spacedim>::ExcFEHasNoSupportPoints()));

  // let's see whether we have the information pre-computed
  if (this->unit_face_support_points[this->n_unique_faces() == 1 ? 0 : face_no]
        .size() != 0)
    return this
      ->unit_face_support_points[this->n_unique_faces() == 1 ? 0 : face_no]
                                [index];
  else
    // no. ask the base element whether it would like to provide this
    // information
    return (
      base_element(this->face_system_to_base_index(index, face_no).first.first)
        .unit_face_support_point(
          this->face_system_to_base_index(index, face_no).second, face_no));
}



template <int dim, int spacedim>
std::pair<Table<2, bool>, std::vector<unsigned int>>
FESystem<dim, spacedim>::get_constant_modes() const
{
  // Note that this->n_components() is actually only an estimate of how many
  // constant modes we will need. There might be more than one such mode
  // (e.g. FE_Q_DG0).
  Table<2, bool> constant_modes(this->n_components(), this->n_dofs_per_cell());
  std::vector<unsigned int> components;
  for (unsigned int i = 0; i < base_elements.size(); ++i)
    {
      const std::pair<Table<2, bool>, std::vector<unsigned int>> base_table =
        base_elements[i].first->get_constant_modes();
      AssertDimension(base_table.first.n_rows(), base_table.second.size());
      const unsigned int element_multiplicity = this->element_multiplicity(i);

      // there might be more than one constant mode for some scalar elements,
      // so make sure the table actually fits: Create a new table with more
      // rows
      const unsigned int comp = components.size();
      if (constant_modes.n_rows() <
          comp + base_table.first.n_rows() * element_multiplicity)
        {
          Table<2, bool> new_constant_modes(comp + base_table.first.n_rows() *
                                                     element_multiplicity,
                                            constant_modes.n_cols());
          for (unsigned int r = 0; r < comp; ++r)
            for (unsigned int c = 0; c < this->n_dofs_per_cell(); ++c)
              new_constant_modes(r, c) = constant_modes(r, c);
          constant_modes.swap(new_constant_modes);
        }

      // next, fill the constant modes from the individual components as well
      // as the component numbers corresponding to the constant mode rows
      for (unsigned int k = 0; k < this->n_dofs_per_cell(); ++k)
        {
          std::pair<std::pair<unsigned int, unsigned int>, unsigned int> ind =
            this->system_to_base_index(k);
          if (ind.first.first == i)
            for (unsigned int c = 0; c < base_table.first.n_rows(); ++c)
              constant_modes(comp +
                               ind.first.second * base_table.first.n_rows() + c,
                             k) = base_table.first(c, ind.second);
        }
      for (unsigned int r = 0; r < element_multiplicity; ++r)
        for (const unsigned int c : base_table.second)
          components.push_back(
            comp + r * this->base_elements[i].first->n_components() + c);
    }
  AssertDimension(components.size(), constant_modes.n_rows());
  return std::pair<Table<2, bool>, std::vector<unsigned int>>(constant_modes,
                                                              components);
}



template <int dim, int spacedim>
void
FESystem<dim, spacedim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &point_values,
  std::vector<double> &              dof_values) const
{
  Assert(this->has_generalized_support_points(),
         ExcMessage("The FESystem does not have generalized support points"));

  AssertDimension(point_values.size(),
                  this->get_generalized_support_points().size());
  AssertDimension(dof_values.size(), this->n_dofs_per_cell());

  std::vector<double>         base_dof_values;
  std::vector<Vector<double>> base_point_values;

  // loop over all base elements (respecting multiplicity) and let them do
  // the work on their share of the input argument

  unsigned int current_vector_component = 0;
  for (unsigned int base = 0; base < base_elements.size(); ++base)
    {
      // We need access to the base_element, its multiplicity, the
      // number of generalized support points (n_base_points) and the
      // number of components we're dealing with.
      const auto &       base_element      = this->base_element(base);
      const unsigned int multiplicity      = this->element_multiplicity(base);
      const unsigned int n_base_dofs       = base_element.n_dofs_per_cell();
      const unsigned int n_base_components = base_element.n_components();

      // If the number of base degrees of freedom is zero, there is nothing
      // to do, skip the rest of the body in this case and continue with
      // the next element
      if (n_base_dofs == 0)
        {
          current_vector_component += multiplicity * n_base_components;
          continue;
        }

      if (base_element.has_generalized_support_points())
        {
          const unsigned int n_base_points =
            base_element.get_generalized_support_points().size();

          base_dof_values.resize(n_base_dofs);
          base_point_values.resize(n_base_points);

          for (unsigned int m = 0; m < multiplicity;
               ++m, current_vector_component += n_base_components)
            {
              // populate base_point_values for a recursive call to
              // convert_generalized_support_point_values_to_dof_values
              for (unsigned int j = 0; j < base_point_values.size(); ++j)
                {
                  base_point_values[j].reinit(n_base_components, false);

                  const auto n =
                    generalized_support_points_index_table[base][j];

                  // we have to extract the correct slice out of the global
                  // vector of values:
                  const auto begin =
                    std::begin(point_values[n]) + current_vector_component;
                  const auto end = begin + n_base_components;
                  std::copy(begin, end, std::begin(base_point_values[j]));
                }

              base_element
                .convert_generalized_support_point_values_to_dof_values(
                  base_point_values, base_dof_values);

              // Finally put these dof values back into global dof values
              // vector.

              // To do this, we could really use a base_to_system_index()
              // function, but that doesn't exist -- just do it by using the
              // reverse table -- the amount of work done here is not worth
              // trying to optimizing this.
              for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
                if (this->system_to_base_index(i).first ==
                    std::make_pair(base, m))
                  dof_values[i] =
                    base_dof_values[this->system_to_base_index(i).second];
            } /*for*/
        }
      else
        {
          // If the base element is non-interpolatory, assign NaN to all
          // DoFs associated to it.

          // To do this, we could really use a base_to_system_index()
          // function, but that doesn't exist -- just do it by using the
          // reverse table -- the amount of work done here is not worth
          // trying to optimizing this.
          for (unsigned int m = 0; m < multiplicity; ++m)
            for (unsigned int i = 0; i < this->n_dofs_per_cell(); ++i)
              if (this->system_to_base_index(i).first ==
                  std::make_pair(base, m))
                dof_values[i] = std::numeric_limits<double>::signaling_NaN();

          current_vector_component += multiplicity * n_base_components;
        }
    } /*for*/
}



template <int dim, int spacedim>
std::size_t
FESystem<dim, spacedim>::memory_consumption() const
{
  // neglect size of data stored in @p{base_elements} due to some problems
  // with the compiler. should be neglectable after all, considering the size
  // of the data of the subelements
  std::size_t mem = (FiniteElement<dim, spacedim>::memory_consumption() +
                     sizeof(base_elements));
  for (unsigned int i = 0; i < base_elements.size(); ++i)
    mem += MemoryConsumption::memory_consumption(*base_elements[i].first);
  return mem;
}



// explicit instantiations
#include "fe_system.inst"

DEAL_II_NAMESPACE_CLOSE
