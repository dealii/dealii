// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/polynomials_p.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_bdm.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <iostream>
#include <memory>
#include <sstream>


DEAL_II_NAMESPACE_OPEN

// TODO: implement the adjust_quad_dof_index_for_face_orientation_table and
// adjust_line_dof_index_for_line_orientation_table fields, and write tests
// similar to bits/face_orientation_and_fe_q_*

template <int dim>
FE_BDM<dim>::FE_BDM(const unsigned int deg)
  : FE_PolyTensor<dim>(
      PolynomialsBDM<dim>(deg),
      FiniteElementData<dim>(get_dpo_vector(deg),
                             dim,
                             deg + 1,
                             FiniteElementData<dim>::Hdiv),
      get_ria_vector(deg),
      std::vector<ComponentMask>(PolynomialsBDM<dim>::n_polynomials(deg),
                                 ComponentMask(std::vector<bool>(dim, true))))
{
  Assert(dim >= 2, ExcImpossibleInDim(dim));
  Assert(
    deg > 0,
    ExcMessage(
      "Lowest order BDM element are degree 1, but you asked for degree 0"));

  const unsigned int n_dofs = this->n_dofs_per_cell();

  this->mapping_kind = {mapping_bdm};
  // These must be done first, since
  // they change the evaluation of
  // basis functions

  // Set up the generalized support
  // points
  initialize_support_points(deg);

  // Now compute the inverse node matrix, generating the correct
  // basis functions from the raw ones. For a discussion of what
  // exactly happens here, see FETools::compute_node_matrix.
  const FullMatrix<double> M = FETools::compute_node_matrix(*this);
  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
  this->inverse_node_matrix.invert(M);
  // From now on, the shape functions provided by FiniteElement::shape_value
  // and similar functions will be the correct ones, not
  // the raw shape functions from the polynomial space anymore.

  // Embedding errors become pretty large, so we just replace the
  // regular threshold in both "computing_..." functions by 1.
  this->reinit_restriction_and_prolongation_matrices(true, true);
  FETools::compute_embedding_matrices(*this, this->prolongation, true, 1.);

  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  FullMatrix<double> face_embeddings[GeometryInfo<dim>::max_children_per_face];
  for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_face; ++i)
    face_embeddings[i].reinit(this->n_dofs_per_face(face_no),
                              this->n_dofs_per_face(face_no));
  FETools::compute_face_embedding_matrices(
    *this, make_array_view(face_embeddings), 0, 0, 1.);
  this->interface_constraints.reinit((1 << (dim - 1)) *
                                       this->n_dofs_per_face(face_no),
                                     this->n_dofs_per_face(face_no));
  unsigned int target_row = 0;
  for (unsigned int d = 0; d < GeometryInfo<dim>::max_children_per_face; ++d)
    for (unsigned int i = 0; i < face_embeddings[d].m(); ++i)
      {
        for (unsigned int j = 0; j < face_embeddings[d].n(); ++j)
          this->interface_constraints(target_row, j) = face_embeddings[d](i, j);
        ++target_row;
      }

  // We need to initialize the dof permutation table and the one for the sign
  // change.
  initialize_quad_dof_index_permutation_and_sign_change();
}


template <int dim>
void
FE_BDM<dim>::initialize_quad_dof_index_permutation_and_sign_change()
{
  // for 1d and 2d, do nothing
  if (dim < 3)
    return;

  // TODO: Implement this for this class
  return;
}


template <int dim>
std::string
FE_BDM<dim>::get_name() const
{
  // note that the
  // FETools::get_fe_by_name
  // function depends on the
  // particular format of the string
  // this function returns, so they
  // have to be kept in synch

  // note that this->degree is the maximal
  // polynomial degree and is thus one higher
  // than the argument given to the
  // constructor
  std::ostringstream namebuf;
  namebuf << "FE_BDM<" << dim << ">(" << this->degree - 1 << ")";

  return namebuf.str();
}


template <int dim>
std::unique_ptr<FiniteElement<dim, dim>>
FE_BDM<dim>::clone() const
{
  return std::make_unique<FE_BDM<dim>>(*this);
}



template <int dim>
void
FE_BDM<dim>::convert_generalized_support_point_values_to_dof_values(
  const std::vector<Vector<double>> &support_point_values,
  std::vector<double>               &nodal_values) const
{
  Assert(support_point_values.size() == this->generalized_support_points.size(),
         ExcDimensionMismatch(support_point_values.size(),
                              this->generalized_support_points.size()));
  AssertDimension(support_point_values[0].size(), dim);
  Assert(nodal_values.size() == this->n_dofs_per_cell(),
         ExcDimensionMismatch(nodal_values.size(), this->n_dofs_per_cell()));

  // First do interpolation on faces. There, the component evaluated
  // depends on the face direction and orientation.

  // The index of the first dof on this face or the cell
  unsigned int dbase = 0;
  // The index of the first generalized support point on this face or the cell
  unsigned int pbase = 0;
  for (auto f : GeometryInfo<dim>::face_indices())
    {
      // Old version with no moments in 2d. See comment below in
      // initialize_support_points()
      if (test_values_face.empty())
        {
          for (unsigned int i = 0; i < this->n_dofs_per_face(f); ++i)
            nodal_values[dbase + i] =
              support_point_values[pbase + i]
                                  [GeometryInfo<dim>::unit_normal_direction[f]];
          pbase += this->n_dofs_per_face(f);
        }
      else
        {
          for (unsigned int i = 0; i < this->n_dofs_per_face(f); ++i)
            {
              double s = 0.;
              for (unsigned int k = 0; k < test_values_face.size(); ++k)
                s +=
                  support_point_values
                    [pbase + k][GeometryInfo<dim>::unit_normal_direction[f]] *
                  test_values_face[k][i];
              nodal_values[dbase + i] = s;
            }
          pbase += test_values_face.size();
        }
      dbase += this->n_dofs_per_face(f);
    }

  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;
  (void)face_no;

  AssertDimension(dbase,
                  this->n_dofs_per_face(face_no) *
                    GeometryInfo<dim>::faces_per_cell);
  AssertDimension(pbase,
                  this->generalized_support_points.size() -
                    test_values_cell.size());

  // Done for BDM1
  if (dbase == this->n_dofs_per_cell())
    return;

  // What's missing are the interior
  // degrees of freedom. In each
  // point, we take all components of
  // the solution.
  Assert((this->n_dofs_per_cell() - dbase) % dim == 0, ExcInternalError());

  for (unsigned int d = 0; d < dim; ++d, dbase += test_values_cell[0].size())
    {
      for (unsigned int i = 0; i < test_values_cell[0].size(); ++i)
        {
          double s = 0.;
          for (unsigned int k = 0; k < test_values_cell.size(); ++k)
            s += support_point_values[pbase + k][d] * test_values_cell[k][i];
          nodal_values[dbase + i] = s;
        }
    }

  Assert(dbase == this->n_dofs_per_cell(), ExcInternalError());
}



template <int dim>
std::vector<unsigned int>
FE_BDM<dim>::get_dpo_vector(const unsigned int deg)
{
  // compute the number of unknowns per cell interior/face/edge
  //
  // for the number of interior dofs, this is the number of
  // polynomials up to degree deg-2 in dim dimensions.
  //
  // the element is face-based and we have as many degrees of freedom
  // on the faces as there are polynomials of degree up to
  // deg. Observe the odd convention of
  // PolynomialSpace::n_polynomials()!

  std::vector<unsigned int> dpo(dim + 1, 0u);
  dpo[dim] =
    (deg > 1 ? dim * PolynomialSpace<dim>::n_polynomials(deg - 1) : 0u);
  dpo[dim - 1] = PolynomialSpace<dim - 1>::n_polynomials(deg + 1);

  return dpo;
}



template <int dim>
std::vector<bool>
FE_BDM<dim>::get_ria_vector(const unsigned int deg)
{
  if (dim == 1)
    {
      Assert(false, ExcImpossibleInDim(1));
      return std::vector<bool>();
    }

  const unsigned int dofs_per_cell = PolynomialsBDM<dim>::n_polynomials(deg);
  const unsigned int dofs_per_face =
    PolynomialSpace<dim - 1>::n_polynomials(deg + 1);

  Assert(GeometryInfo<dim>::faces_per_cell * dofs_per_face <= dofs_per_cell,
         ExcInternalError());

  // all dofs need to be
  // non-additive, since they have
  // continuity requirements.
  // however, the interior dofs are
  // made additive
  std::vector<bool> ret_val(dofs_per_cell, false);
  for (unsigned int i = GeometryInfo<dim>::faces_per_cell * dofs_per_face;
       i < dofs_per_cell;
       ++i)
    ret_val[i] = true;

  return ret_val;
}


namespace internal
{
  namespace FE_BDM
  {
    namespace
    {
      // This function sets up the values of the polynomials we want to
      // take moments with in the quadrature points. In fact, we multiply
      // those by the weights, such that the sum of function values and
      // test_values over quadrature points yields the interpolated degree
      // of freedom.
      template <int dim>
      void
      initialize_test_values(std::vector<std::vector<double>> &test_values,
                             const Quadrature<dim>            &quadrature,
                             const unsigned int                deg)
      {
        PolynomialsP<dim>           poly(deg);
        std::vector<Tensor<1, dim>> dummy1;
        std::vector<Tensor<2, dim>> dummy2;
        std::vector<Tensor<3, dim>> dummy3;
        std::vector<Tensor<4, dim>> dummy4;

        test_values.resize(quadrature.size());

        for (unsigned int k = 0; k < quadrature.size(); ++k)
          {
            test_values[k].resize(poly.n());
            poly.evaluate(quadrature.point(k),
                          test_values[k],
                          dummy1,
                          dummy2,
                          dummy3,
                          dummy4);
            for (unsigned int i = 0; i < poly.n(); ++i)
              {
                test_values[k][i] *= quadrature.weight(k);
              }
          }
      }

      // This specialization only serves to avoid error messages. Nothing
      // useful can be computed in dimension zero and thus the vector
      // length stays zero.
      template <>
      void
      initialize_test_values(std::vector<std::vector<double>> &,
                             const Quadrature<0> &,
                             const unsigned int)
      {}
    } // namespace
  }   // namespace FE_BDM
} // namespace internal


template <int dim>
void
FE_BDM<dim>::initialize_support_points(const unsigned int deg)
{
  // Our support points are quadrature points on faces and inside the
  // cell. First on the faces, we have to test polynomials of degree
  // up to deg, which means we need dg+1 points in each direction. The
  // fact that we do not have tensor product polynomials will be
  // considered later. In 2d, we can use point values.
  const QGauss<dim - 1> face_points(deg + 1);

  // TODO: the implementation makes the assumption that all faces have the
  // same number of dofs
  AssertDimension(this->n_unique_faces(), 1);
  const unsigned int face_no = 0;

  // Copy the quadrature formula to the face points.
  this->generalized_face_support_points[face_no].resize(face_points.size());
  for (unsigned int k = 0; k < face_points.size(); ++k)
    this->generalized_face_support_points[face_no][k] = face_points.point(k);

  // In the interior, we only test with polynomials of degree up to
  // deg-2, thus we use deg points. Note that deg>=1 and the lowest
  // order element has no points in the cell, such that we have to
  // distinguish this case.
  const QGauss<dim> cell_points(deg == 1 ? 0 : deg);

  // Compute the size of the whole support point set
  const unsigned int npoints =
    GeometryInfo<dim>::faces_per_cell * face_points.size() + cell_points.size();

  this->generalized_support_points.resize(npoints);

  const Quadrature<dim> faces =
    QProjector<dim>::project_to_all_faces(this->reference_cell(), face_points);

  for (unsigned int face_no = 0; face_no < GeometryInfo<dim>::faces_per_cell;
       ++face_no)
    {
      const auto offset = QProjector<dim>::DataSetDescriptor::face(
        this->reference_cell(),
        face_no,
        numbers::default_geometric_orientation,
        face_points.size());
      for (unsigned int face_point = 0; face_point < face_points.size();
           ++face_point)
        this->generalized_support_points[face_points.size() * face_no +
                                         face_point] =
          faces.point(offset + face_point);
    }

  // Currently, for backward compatibility, we do not use moments, but
  // point values on faces in 2d. In 3d, this is impossible, since the
  // moments are only taken with respect to PolynomialsP.
  if (dim > 2)
    internal::FE_BDM::initialize_test_values(test_values_face,
                                             face_points,
                                             deg);

  if (deg <= 1)
    return;

  // Remember where interior points start
  const unsigned int ibase =
    face_points.size() * GeometryInfo<dim>::faces_per_cell;
  for (unsigned int k = 0; k < cell_points.size(); ++k)
    {
      this->generalized_support_points[ibase + k] = cell_points.point(k);
    }
  // Finally, compute the values of
  // the test functions in the
  // interior quadrature points

  internal::FE_BDM::initialize_test_values(test_values_cell,
                                           cell_points,
                                           deg - 2);
}



/*-------------- Explicit Instantiations -------------------------------*/
#include "fe/fe_bdm.inst"

DEAL_II_NAMESPACE_CLOSE
