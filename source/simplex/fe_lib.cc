// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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

#include <deal.II/base/config.h>

#include <deal.II/simplex/fe_lib.h>

DEAL_II_NAMESPACE_OPEN

namespace Simplex
{
  namespace
  {
    /**
     * Helper function to set up the dpo vector of FE_P for a given @p dim and
     * @p degree.
     */
    std::vector<unsigned int>
    get_dpo_vector_fe_p(const unsigned int dim, const unsigned int degree)
    {
      std::vector<unsigned int> dpo(dim + 1, 0U);

      if (degree == 1)
        {
          // one dof at each vertex
          dpo[0] = 1;
        }
      else if (degree == 2)
        {
          // one dof at each vertex and in the middle of each line
          dpo[0] = 1;
          dpo[1] = 1;
          dpo[2] = 0;
        }
      else
        {
          Assert(false, ExcNotImplemented());
        }

      return dpo;
    }

    /**
     * Helper function to set up the dpo vector of FE_DGP for a given @p dim and
     * @p degree.
     */
    std::vector<unsigned int>
    get_dpo_vector_fe_dgp(const unsigned int dim, const unsigned int degree)
    {
      std::vector<unsigned int> dpo(dim + 1, 0U);

      // all dofs are internal
      if (dim == 2 && degree == 1)
        dpo[dim] = 3;
      else if (dim == 2 && degree == 2)
        dpo[dim] = 6;
      else if (dim == 3 && degree == 1)
        dpo[dim] = 4;
      else if (dim == 3 && degree == 2)
        dpo[dim] = 10;
      else
        {
          Assert(false, ExcNotImplemented());
        }

      return dpo;
    }
  } // namespace



  template <int dim, int spacedim>
  FE_Poly<dim, spacedim>::FE_Poly(const unsigned int               degree,
                                  const std::vector<unsigned int> &dpo_vector)
    : dealii::FE_Poly<dim, spacedim>(
        Simplex::ScalarPolynomial<dim>(degree),
        FiniteElementData<dim>(dpo_vector,
                               dim == 2 ? ReferenceCell::Type::Tri :
                                          ReferenceCell::Type::Tet,
                               1,
                               degree,
                               FiniteElementData<dim>::L2),
        std::vector<bool>(FiniteElementData<dim>(dpo_vector,
                                                 dim == 2 ?
                                                   ReferenceCell::Type::Tri :
                                                   ReferenceCell::Type::Tet,
                                                 1,
                                                 degree)
                            .dofs_per_cell,
                          true),
        std::vector<ComponentMask>(
          FiniteElementData<dim>(dpo_vector,
                                 dim == 2 ? ReferenceCell::Type::Tri :
                                            ReferenceCell::Type::Tet,
                                 1,
                                 degree)
            .dofs_per_cell,
          std::vector<bool>(1, true)))
  {
    this->unit_support_points.clear();

    if (dim == 2)
      {
        if (degree == 1)
          {
            this->unit_support_points.emplace_back(0.0, 0.0);
            this->unit_support_points.emplace_back(1.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 1.0);

            // TODO
            this->unit_face_support_points[0].emplace_back(0.0);
            this->unit_face_support_points[0].emplace_back(1.0);
          }
        else if (degree == 2)
          {
            this->unit_support_points.emplace_back(0.0, 0.0);
            this->unit_support_points.emplace_back(1.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 1.0);
            this->unit_support_points.emplace_back(0.5, 0.0);
            this->unit_support_points.emplace_back(0.5, 0.5);
            this->unit_support_points.emplace_back(0.0, 0.5);

            // TODO
            this->unit_face_support_points[0].emplace_back(0.0);
            this->unit_face_support_points[0].emplace_back(1.0);
            this->unit_face_support_points[0].emplace_back(0.5);
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
    else if (dim == 3)
      {
        if (degree == 1)
          {
            this->unit_support_points.emplace_back(0.0, 0.0, 0.0);
            this->unit_support_points.emplace_back(1.0, 0.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 1.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 0.0, 1.0);

            // TODO
            this->unit_face_support_points[0].emplace_back(1.0, 0.0);
            this->unit_face_support_points[0].emplace_back(0.0, 1.0);
            this->unit_face_support_points[0].emplace_back(0.0, 0.0);
          }
        else if (degree == 2)
          {
            this->unit_support_points.emplace_back(0.0, 0.0, 0.0);
            this->unit_support_points.emplace_back(1.0, 0.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 1.0, 0.0);
            this->unit_support_points.emplace_back(0.0, 0.0, 1.0);
            this->unit_support_points.emplace_back(0.5, 0.0, 0.0);
            this->unit_support_points.emplace_back(0.5, 0.5, 0.0);
            this->unit_support_points.emplace_back(0.0, 0.5, 0.0);
            this->unit_support_points.emplace_back(0.0, 0.0, 0.5);
            this->unit_support_points.emplace_back(0.5, 0.0, 0.5);
            this->unit_support_points.emplace_back(0.0, 0.5, 0.5);

            // TODO
            this->unit_face_support_points[0].emplace_back(1.0, 0.0);
            this->unit_face_support_points[0].emplace_back(0.0, 1.0);
            this->unit_face_support_points[0].emplace_back(0.0, 0.0);
            this->unit_face_support_points[0].emplace_back(0.5, 0.5);
            this->unit_face_support_points[0].emplace_back(0.0, 0.5);
            this->unit_face_support_points[0].emplace_back(0.5, 0.0);
          }
        else
          {
            Assert(false, ExcNotImplemented());
          }
      }
    else
      {
        Assert(false, ExcNotImplemented());
      }
  }



  template <int dim, int spacedim>
  void
  FE_Poly<dim, spacedim>::
    convert_generalized_support_point_values_to_dof_values(
      const std::vector<Vector<double>> &support_point_values,
      std::vector<double> &              nodal_values) const
  {
    AssertDimension(support_point_values.size(),
                    this->get_unit_support_points().size());
    AssertDimension(support_point_values.size(), nodal_values.size());
    AssertDimension(this->dofs_per_cell, nodal_values.size());

    for (unsigned int i = 0; i < this->dofs_per_cell; ++i)
      {
        AssertDimension(support_point_values[i].size(), 1);

        nodal_values[i] = support_point_values[i](0);
      }
  }



  template <int dim, int spacedim>
  FE_P<dim, spacedim>::FE_P(const unsigned int degree)
    : FE_Poly<dim, spacedim>(degree, get_dpo_vector_fe_p(dim, degree))
  {}



  template <int dim, int spacedim>
  std::unique_ptr<FiniteElement<dim, spacedim>>
  FE_P<dim, spacedim>::clone() const
  {
    return std::make_unique<FE_P<dim, spacedim>>(*this);
  }



  template <int dim, int spacedim>
  std::string
  FE_P<dim, spacedim>::get_name() const
  {
    std::ostringstream namebuf;
    namebuf << "FE_P<" << dim << ">(" << this->degree << ")";

    return namebuf.str();
  }


  template <int dim, int spacedim>
  FE_DGP<dim, spacedim>::FE_DGP(const unsigned int degree)
    : FE_Poly<dim, spacedim>(degree, get_dpo_vector_fe_dgp(dim, degree))
  {}



  template <int dim, int spacedim>
  std::unique_ptr<FiniteElement<dim, spacedim>>
  FE_DGP<dim, spacedim>::clone() const
  {
    return std::make_unique<FE_DGP<dim, spacedim>>(*this);
  }



  template <int dim, int spacedim>
  std::string
  FE_DGP<dim, spacedim>::get_name() const
  {
    std::ostringstream namebuf;
    namebuf << "FE_DGP<" << dim << ">(" << this->degree << ")";

    return namebuf.str();
  }



} // namespace Simplex

// explicit instantiations
#include "fe_lib.inst"

DEAL_II_NAMESPACE_CLOSE
