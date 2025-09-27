// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_tools.h>

#include <memory>
#include <sstream>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace FE_DGPMonomial
  {
    namespace
    {
      // storage of hand-chosen support
      // points
      //
      // For dim=2, dofs_per_cell of
      // FE_DGPMonomial(k) is given by
      // 0.5(k+1)(k+2), i.e.
      //
      // k    0  1  2  3  4  5  6  7
      // dofs 1  3  6 10 15 21 28 36
      //
      // indirect access of unit points:
      // the points for degree k are
      // located at
      //
      // points[start_index[k]..start_index[k+1]-1]
      const unsigned int start_index2d[6] = {0, 1, 4, 10, 20, 35};
      const double       points2d[35][2]  = {
        {0, 0},       {0, 0},       {1, 0},       {0, 1},       {0, 0},
        {1, 0},       {0, 1},       {1, 1},       {0.5, 0},     {0, 0.5},
        {0, 0},       {1, 0},       {0, 1},       {1, 1},       {1. / 3., 0},
        {2. / 3., 0}, {0, 1. / 3.}, {0, 2. / 3.}, {0.5, 1},     {1, 0.5},
        {0, 0},       {1, 0},       {0, 1},       {1, 1},       {0.25, 0},
        {0.5, 0},     {0.75, 0},    {0, 0.25},    {0, 0.5},     {0, 0.75},
        {1. / 3., 1}, {2. / 3., 1}, {1, 1. / 3.}, {1, 2. / 3.}, {0.5, 0.5}};

      // For dim=3, dofs_per_cell of
      // FE_DGPMonomial(k) is given by
      // 1./6.(k+1)(k+2)(k+3), i.e.
      //
      // k    0  1  2  3  4  5  6   7
      // dofs 1  4 10 20 35 56 84 120
      const unsigned int start_index3d[6] = {0, 1, 5, 15 /*,35*/};
      const double       points3d[35][3]  = {{0, 0, 0},
                                             {0, 0, 0},
                                             {1, 0, 0},
                                             {0, 1, 0},
                                             {0, 0, 1},
                                             {0, 0, 0},
                                             {1, 0, 0},
                                             {0, 1, 0},
                                             {0, 0, 1},
                                             {0.5, 0, 0},
                                             {0, 0.5, 0},
                                             {0, 0, 0.5},
                                             {1, 1, 0},
                                             {1, 0, 1},
                                             {0, 1, 1}};


      template <int dim>
      void
      generate_unit_points(const unsigned int, std::vector<Point<dim>> &);

      template <>
      void
      generate_unit_points(const unsigned int k, std::vector<Point<1>> &p)
      {
        Assert(p.size() == k + 1, ExcDimensionMismatch(p.size(), k + 1));
        const double h = 1. / k;
        for (unsigned int i = 0; i < p.size(); ++i)
          p[i][0] = i * h;
      }

      template <>
      void
      generate_unit_points(const unsigned int k, std::vector<Point<2>> &p)
      {
        Assert(k <= 4, ExcNotImplemented());
        Assert(p.size() == start_index2d[k + 1] - start_index2d[k],
               ExcInternalError());
        for (unsigned int i = 0; i < p.size(); ++i)
          {
            p[i][0] = points2d[start_index2d[k] + i][0];
            p[i][1] = points2d[start_index2d[k] + i][1];
          }
      }

      template <>
      void
      generate_unit_points(const unsigned int k, std::vector<Point<3>> &p)
      {
        Assert(k <= 2, ExcNotImplemented());
        Assert(p.size() == start_index3d[k + 1] - start_index3d[k],
               ExcInternalError());
        for (unsigned int i = 0; i < p.size(); ++i)
          {
            p[i][0] = points3d[start_index3d[k] + i][0];
            p[i][1] = points3d[start_index3d[k] + i][1];
            p[i][2] = points3d[start_index3d[k] + i][2];
          }
      }
    } // namespace
  }   // namespace FE_DGPMonomial
} // namespace internal



template <int dim>
FE_DGPMonomial<dim>::FE_DGPMonomial(const unsigned int degree)
  : FE_Poly<dim>(PolynomialsP<dim>(degree),
                 FiniteElementData<dim>(get_dpo_vector(degree),
                                        1,
                                        degree,
                                        FiniteElementData<dim>::L2),
                 std::vector<bool>(
                   FiniteElementData<dim>(get_dpo_vector(degree), 1, degree)
                     .n_dofs_per_cell(),
                   true),
                 std::vector<ComponentMask>(
                   FiniteElementData<dim>(get_dpo_vector(degree), 1, degree)
                     .n_dofs_per_cell(),
                   ComponentMask(std::vector<bool>(1, true))))
{
  Assert(this->poly_space->n() == this->n_dofs_per_cell(), ExcInternalError());
  Assert(this->poly_space->degree() == this->degree, ExcInternalError());

  // DG doesn't have constraints, so
  // leave them empty

  // Reinit the vectors of
  // restriction and prolongation
  // matrices to the right sizes
  this->reinit_restriction_and_prolongation_matrices();
  // Fill prolongation matrices with embedding operators
  FETools::compute_embedding_matrices(*this, this->prolongation);
  // Fill restriction matrices with L2-projection
  FETools::compute_projection_matrices(*this, this->restriction);
}



template <int dim>
std::string
FE_DGPMonomial<dim>::get_name() const
{
  // note that the
  // FETools::get_fe_by_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_DGPMonomial<" << dim << ">(" << this->degree << ")";

  return namebuf.str();
}



template <int dim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_DGPMonomial<dim>::clone() const
{
  return std::make_unique<FE_DGPMonomial<dim>>(*this);
}



// TODO: Remove this function and use the one in FETools, if needed
template <int dim>
void
FE_DGPMonomial<dim>::get_interpolation_matrix(
  const FiniteElement<dim> &source_fe,
  FullMatrix<double>       &interpolation_matrix) const
{
  const FE_DGPMonomial<dim> *source_dgp_monomial =
    dynamic_cast<const FE_DGPMonomial<dim> *>(&source_fe);

  if (source_dgp_monomial)
    {
      // ok, source_fe is a DGP_Monomial
      // element. Then, the interpolation
      // matrix is simple
      const unsigned int m = interpolation_matrix.m();
      const unsigned int n = interpolation_matrix.n();
      (void)m;
      (void)n;
      Assert(m == this->n_dofs_per_cell(),
             ExcDimensionMismatch(m, this->n_dofs_per_cell()));
      Assert(n == source_dgp_monomial->n_dofs_per_cell(),
             ExcDimensionMismatch(n, source_dgp_monomial->n_dofs_per_cell()));

      const unsigned int min_mn =
        interpolation_matrix.m() < interpolation_matrix.n() ?
          interpolation_matrix.m() :
          interpolation_matrix.n();

      for (unsigned int i = 0; i < min_mn; ++i)
        interpolation_matrix(i, i) = 1.;
    }
  else
    {
      std::vector<Point<dim>> unit_points(this->n_dofs_per_cell());
      internal::FE_DGPMonomial::generate_unit_points(this->degree, unit_points);

      FullMatrix<double> source_fe_matrix(unit_points.size(),
                                          source_fe.n_dofs_per_cell());
      for (unsigned int j = 0; j < source_fe.n_dofs_per_cell(); ++j)
        for (unsigned int k = 0; k < unit_points.size(); ++k)
          source_fe_matrix(k, j) = source_fe.shape_value(j, unit_points[k]);

      FullMatrix<double> this_matrix(this->n_dofs_per_cell(),
                                     this->n_dofs_per_cell());
      for (unsigned int j = 0; j < this->n_dofs_per_cell(); ++j)
        for (unsigned int k = 0; k < unit_points.size(); ++k)
          this_matrix(k, j) =
            this->poly_space->compute_value(j, unit_points[k]);

      this_matrix.gauss_jordan();

      this_matrix.mmult(interpolation_matrix, source_fe_matrix);
    }
}



template <int dim>
void
FE_DGPMonomial<dim>::initialize_restriction()
{
  DEAL_II_NOT_IMPLEMENTED();
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim>
std::vector<unsigned int>
FE_DGPMonomial<dim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim + 1, 0U);
  dpo[dim] = deg + 1;
  for (unsigned int i = 1; i < dim; ++i)
    {
      dpo[dim] *= deg + 1 + i;
      dpo[dim] /= i + 1;
    }
  return dpo;
}


template <int dim>
void
FE_DGPMonomial<dim>::get_face_interpolation_matrix(
  const FiniteElement<dim> &x_source_fe,
  FullMatrix<double>       &interpolation_matrix,
  const unsigned int) const
{
  // this is only implemented, if the source
  // FE is also a DGPMonomial element. in that case,
  // both elements have no dofs on their
  // faces and the face interpolation matrix
  // is necessarily empty -- i.e. there isn't
  // much we need to do here.
  (void)interpolation_matrix;
  AssertThrow((x_source_fe.get_name().find("FE_DGPMonomial<") == 0) ||
                (dynamic_cast<const FE_DGPMonomial<dim> *>(&x_source_fe) !=
                 nullptr),
              typename FiniteElement<dim>::ExcInterpolationNotImplemented());

  Assert(interpolation_matrix.m() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0,
         ExcDimensionMismatch(interpolation_matrix.n(), 0));
}



template <int dim>
void
FE_DGPMonomial<dim>::get_subface_interpolation_matrix(
  const FiniteElement<dim> &x_source_fe,
  const unsigned int,
  FullMatrix<double> &interpolation_matrix,
  const unsigned int) const
{
  // this is only implemented, if the source
  // FE is also a DGPMonomial element. in that case,
  // both elements have no dofs on their
  // faces and the face interpolation matrix
  // is necessarily empty -- i.e. there isn't
  // much we need to do here.
  (void)interpolation_matrix;
  AssertThrow((x_source_fe.get_name().find("FE_DGPMonomial<") == 0) ||
                (dynamic_cast<const FE_DGPMonomial<dim> *>(&x_source_fe) !=
                 nullptr),
              typename FiniteElement<dim>::ExcInterpolationNotImplemented());

  Assert(interpolation_matrix.m() == 0,
         ExcDimensionMismatch(interpolation_matrix.m(), 0));
  Assert(interpolation_matrix.n() == 0,
         ExcDimensionMismatch(interpolation_matrix.n(), 0));
}



template <int dim>
bool
FE_DGPMonomial<dim>::hp_constraints_are_implemented() const
{
  return true;
}



template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGPMonomial<dim>::hp_vertex_dof_identities(
  const FiniteElement<dim> &fe_other) const
{
  // there are no such constraints for DGPMonomial
  // elements at all
  if (dynamic_cast<const FE_DGPMonomial<dim> *>(&fe_other) != nullptr)
    return std::vector<std::pair<unsigned int, unsigned int>>();
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}



template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGPMonomial<dim>::hp_line_dof_identities(
  const FiniteElement<dim> &fe_other) const
{
  // there are no such constraints for DGPMonomial
  // elements at all
  if (dynamic_cast<const FE_DGPMonomial<dim> *>(&fe_other) != nullptr)
    return std::vector<std::pair<unsigned int, unsigned int>>();
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}



template <int dim>
std::vector<std::pair<unsigned int, unsigned int>>
FE_DGPMonomial<dim>::hp_quad_dof_identities(const FiniteElement<dim> &fe_other,
                                            const unsigned int) const
{
  // there are no such constraints for DGPMonomial
  // elements at all
  if (dynamic_cast<const FE_DGPMonomial<dim> *>(&fe_other) != nullptr)
    return std::vector<std::pair<unsigned int, unsigned int>>();
  else
    {
      DEAL_II_NOT_IMPLEMENTED();
      return std::vector<std::pair<unsigned int, unsigned int>>();
    }
}



template <int dim>
FiniteElementDomination::Domination
FE_DGPMonomial<dim>::compare_for_domination(const FiniteElement<dim> &fe_other,
                                            const unsigned int codim) const
{
  Assert(codim <= dim, ExcImpossibleInDim(dim));

  // vertex/line/face domination
  // ---------------------------
  if (codim > 0)
    // this is a discontinuous element, so by definition there will
    // be no constraints wherever this element comes together with
    // any other kind of element
    return FiniteElementDomination::no_requirements;

  // cell domination
  // ---------------
  if (const FE_DGPMonomial<dim> *fe_monomial_other =
        dynamic_cast<const FE_DGPMonomial<dim> *>(&fe_other))
    {
      if (this->degree < fe_monomial_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_monomial_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (const FE_Nothing<dim> *fe_nothing =
             dynamic_cast<const FE_Nothing<dim> *>(&fe_other))
    {
      if (fe_nothing->is_dominating())
        return FiniteElementDomination::other_element_dominates;
      else
        // the FE_Nothing has no degrees of freedom and it is typically used
        // in a context where we don't require any continuity along the
        // interface
        return FiniteElementDomination::no_requirements;
    }

  DEAL_II_NOT_IMPLEMENTED();
  return FiniteElementDomination::neither_element_dominates;
}



template <>
bool
FE_DGPMonomial<1>::has_support_on_face(const unsigned int,
                                       const unsigned int face_index) const
{
  return face_index == 1 || (face_index == 0 && this->degree == 0);
}



template <>
bool
FE_DGPMonomial<2>::has_support_on_face(const unsigned int shape_index,
                                       const unsigned int face_index) const
{
  bool support_on_face = false;
  if (face_index == 1 || face_index == 2)
    support_on_face = true;
  else
    {
      auto *const polynomial_space_p =
        dynamic_cast<PolynomialsP<2> *>(this->poly_space.get());
      Assert(polynomial_space_p != nullptr, ExcInternalError());
      const std::array<unsigned int, 2> degrees =
        polynomial_space_p->directional_degrees(shape_index);

      if ((face_index == 0 && degrees[1] == 0) ||
          (face_index == 3 && degrees[0] == 0))
        support_on_face = true;
    }
  return support_on_face;
}



template <>
bool
FE_DGPMonomial<3>::has_support_on_face(const unsigned int shape_index,
                                       const unsigned int face_index) const
{
  bool support_on_face = false;
  if (face_index == 1 || face_index == 3 || face_index == 4)
    support_on_face = true;
  else
    {
      auto *const polynomial_space_p =
        dynamic_cast<PolynomialsP<3> *>(this->poly_space.get());
      Assert(polynomial_space_p != nullptr, ExcInternalError());
      const std::array<unsigned int, 3> degrees =
        polynomial_space_p->directional_degrees(shape_index);

      if ((face_index == 0 && degrees[1] == 0) ||
          (face_index == 2 && degrees[2] == 0) ||
          (face_index == 5 && degrees[0] == 0))
        support_on_face = true;
    }
  return support_on_face;
}



template <int dim>
std::size_t
FE_DGPMonomial<dim>::memory_consumption() const
{
  DEAL_II_NOT_IMPLEMENTED();
  return 0;
}



// explicit instantiations
#include "fe/fe_dgp_monomial.inst"


DEAL_II_NAMESPACE_CLOSE
