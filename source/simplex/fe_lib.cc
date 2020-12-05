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

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

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

    /**
     * Helper function to set up the dpo vector of FE_WedgeP for a given @p degree.
     */
    internal::GenericDoFsPerObject
    get_dpo_vector_fe_wedge(const unsigned int degree)
    {
      internal::GenericDoFsPerObject dpo;

      // all dofs are internal
      if (degree == 1)
        {
          dpo.dofs_per_object_exclusive  = {{1}, {0}, {0, 0, 0, 0, 0}, {0}};
          dpo.dofs_per_object_inclusive  = {{1}, {2}, {3, 3, 4, 4, 4}, {6}};
          dpo.object_index               = {{}, {6}, {6}, {6}};
          dpo.first_object_index_on_face = {{},
                                            {3, 3, 4, 4, 4},
                                            {3, 3, 4, 4, 4}};
        }
      else if (degree == 2)
        {
          dpo.dofs_per_object_exclusive = {{1}, {1}, {0, 0, 1, 1, 1}, {0}};
          dpo.dofs_per_object_inclusive = {{1}, {3}, {6, 6, 9, 9, 9}, {18}};
          dpo.object_index              = {{}, {6}, {15, 15, 15, 16, 17}, {18}};
          dpo.first_object_index_on_face = {{},
                                            {3, 3, 4, 4, 4},
                                            {6, 6, 8, 8, 8}};
        }
      else
        {
          Assert(false, ExcNotImplemented());
        }

      return dpo;
    }

    /**
     * Helper function to set up the dpo vector of FE_WedgeP for a given @p degree.
     */
    internal::GenericDoFsPerObject
    get_dpo_vector_fe_pyramid(const unsigned int degree)
    {
      internal::GenericDoFsPerObject dpo;

      // all dofs are internal
      if (degree == 1)
        {
          dpo.dofs_per_object_exclusive  = {{1}, {0}, {0, 0, 0, 0, 0}, {0}};
          dpo.dofs_per_object_inclusive  = {{1}, {2}, {4, 3, 3, 3, 3}, {5}};
          dpo.object_index               = {{}, {5}, {5}, {5}};
          dpo.first_object_index_on_face = {{},
                                            {4, 3, 3, 3, 3},
                                            {4, 3, 3, 3, 3}};
        }
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
  FiniteElementDomination::Domination
  FE_P<dim, spacedim>::compare_for_domination(
    const FiniteElement<dim, spacedim> &fe_other,
    const unsigned int                  codim) const
  {
    (void)fe_other;
    (void)codim;

    Assert((dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
           ExcNotImplemented());
    AssertDimension(dim, 2);
    AssertDimension(this->degree, fe_other.tensor_degree());

    return FiniteElementDomination::either_element_can_dominate;
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_P<dim, spacedim>::hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const
  {
    (void)fe_other;

    Assert((dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
           ExcNotImplemented());
    AssertDimension(dim, 2);
    AssertDimension(this->degree, fe_other.tensor_degree());

    return {{0, 0}};
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_P<dim, spacedim>::hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const
  {
    (void)fe_other;

    Assert((dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
           ExcNotImplemented());
    AssertDimension(dim, 2);
    AssertDimension(this->degree, fe_other.tensor_degree());

    std::vector<std::pair<unsigned int, unsigned int>> result;

    for (unsigned int i = 0; i < this->degree - 1; ++i)
      result.emplace_back(i, i);

    return result;
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


  template <int dim, int spacedim>
  FiniteElementDomination::Domination
  FE_DGP<dim, spacedim>::compare_for_domination(
    const FiniteElement<dim, spacedim> &fe_other,
    const unsigned int                  codim) const
  {
    (void)fe_other;
    (void)codim;

    Assert((dynamic_cast<const FE_DGQ<dim, spacedim> *>(&fe_other)),
           ExcNotImplemented());
    AssertDimension(dim, 2);
    AssertDimension(this->degree, fe_other.tensor_degree());

    return FiniteElementDomination::either_element_can_dominate;
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_DGP<dim, spacedim>::hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const
  {
    (void)fe_other;

    return {};
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_DGP<dim, spacedim>::hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const
  {
    (void)fe_other;

    return {};
  }



  template <int dim, int spacedim>
  FE_WedgeP<dim, spacedim>::FE_WedgeP(const unsigned int degree)
    : dealii::FE_Poly<dim, spacedim>(
        Simplex::ScalarWedgePolynomial<dim>(degree),
        FiniteElementData<dim>(get_dpo_vector_fe_wedge(degree),
                               ReferenceCell::Type::Wedge,
                               1,
                               degree,
                               FiniteElementData<dim>::L2),
        std::vector<bool>(FiniteElementData<dim>(get_dpo_vector_fe_wedge(
                                                   degree),
                                                 ReferenceCell::Type::Wedge,
                                                 1,
                                                 degree)
                            .dofs_per_cell,
                          true),
        std::vector<ComponentMask>(
          FiniteElementData<dim>(get_dpo_vector_fe_wedge(degree),
                                 ReferenceCell::Type::Wedge,
                                 1,
                                 degree)
            .dofs_per_cell,
          std::vector<bool>(1, true)))
  {
    AssertDimension(dim, 3);


    if (degree == 1)
      {
        this->unit_support_points.emplace_back(1.0, 0.0, 0.0);
        this->unit_support_points.emplace_back(0.0, 1.0, 0.0);
        this->unit_support_points.emplace_back(0.0, 0.0, 0.0);
        this->unit_support_points.emplace_back(1.0, 0.0, 1.0);
        this->unit_support_points.emplace_back(0.0, 1.0, 1.0);
        this->unit_support_points.emplace_back(0.0, 0.0, 1.0);

        // TODO
        this->unit_face_support_points[0].emplace_back(1.0, 0.0);
        this->unit_face_support_points[0].emplace_back(0.0, 1.0);
        this->unit_face_support_points[0].emplace_back(0.0, 0.0);
      }
  }



  template <int dim, int spacedim>
  std::unique_ptr<FiniteElement<dim, spacedim>>
  FE_WedgeP<dim, spacedim>::clone() const
  {
    return std::make_unique<FE_WedgeP<dim, spacedim>>(*this);
  }



  template <int dim, int spacedim>
  std::string
  FE_WedgeP<dim, spacedim>::get_name() const
  {
    std::ostringstream namebuf;
    namebuf << "FE_WedgeP<" << dim << ">(" << this->degree << ")";

    return namebuf.str();
  }



  template <int dim, int spacedim>
  FiniteElementDomination::Domination
  FE_WedgeP<dim, spacedim>::compare_for_domination(
    const FiniteElement<dim, spacedim> &fe_other,
    const unsigned int                  codim) const
  {
    (void)fe_other;
    (void)codim;

    Assert((dynamic_cast<const Simplex::FE_P<dim, spacedim> *>(&fe_other)) ||
             (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
           ExcNotImplemented());

    return FiniteElementDomination::this_element_dominates;
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_WedgeP<dim, spacedim>::hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const
  {
    (void)fe_other;

    Assert((dynamic_cast<const Simplex::FE_P<dim, spacedim> *>(&fe_other)) ||
             (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
           ExcNotImplemented());

    return {{0, 0}};
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_WedgeP<dim, spacedim>::hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const
  {
    (void)fe_other;

    Assert((dynamic_cast<const Simplex::FE_P<dim, spacedim> *>(&fe_other)) ||
             (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
           ExcNotImplemented());

    std::vector<std::pair<unsigned int, unsigned int>> result;

    for (unsigned int i = 0; i < this->degree - 1; ++i)
      result.emplace_back(i, i);

    return result;
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_WedgeP<dim, spacedim>::hp_quad_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other,
    const unsigned int                  face_no) const
  {
    (void)fe_other;


    AssertIndexRange(face_no, 5);

    if (face_no < 2)
      {
        Assert((dynamic_cast<const Simplex::FE_P<dim, spacedim> *>(&fe_other)),
               ExcNotImplemented());
      }
    else
      {
        Assert((dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
               ExcNotImplemented());
      }

    std::vector<std::pair<unsigned int, unsigned int>> result;

    for (unsigned int i = 0; i < this->n_dofs_per_quad(face_no); ++i)
      result.emplace_back(i, i);

    return result;
  }



  template <int dim, int spacedim>
  FE_PyramidP<dim, spacedim>::FE_PyramidP(const unsigned int degree)
    : dealii::FE_Poly<dim, spacedim>(
        Simplex::ScalarPyramidPolynomial<dim>(degree),
        FiniteElementData<dim>(get_dpo_vector_fe_pyramid(degree),
                               ReferenceCell::Type::Pyramid,
                               1,
                               degree,
                               FiniteElementData<dim>::H1),
        std::vector<bool>(FiniteElementData<dim>(get_dpo_vector_fe_pyramid(
                                                   degree),
                                                 ReferenceCell::Type::Pyramid,
                                                 1,
                                                 degree)
                            .dofs_per_cell,
                          true),
        std::vector<ComponentMask>(
          FiniteElementData<dim>(get_dpo_vector_fe_pyramid(degree),
                                 ReferenceCell::Type::Pyramid,
                                 1,
                                 degree)
            .dofs_per_cell,
          std::vector<bool>(1, true)))
  {
    AssertDimension(dim, 3);


    if (degree == 1)
      {
        this->unit_support_points.emplace_back(-1.0, -1.0, 0.0);
        this->unit_support_points.emplace_back(+1.0, -1.0, 0.0);
        this->unit_support_points.emplace_back(-1.0, +1.0, 0.0);
        this->unit_support_points.emplace_back(+1.0, +1.0, 0.0);
        this->unit_support_points.emplace_back(+0.0, +0.0, 1.0);
      }
  }



  template <int dim, int spacedim>
  std::unique_ptr<FiniteElement<dim, spacedim>>
  FE_PyramidP<dim, spacedim>::clone() const
  {
    return std::make_unique<FE_PyramidP<dim, spacedim>>(*this);
  }



  template <int dim, int spacedim>
  std::string
  FE_PyramidP<dim, spacedim>::get_name() const
  {
    std::ostringstream namebuf;
    namebuf << "FE_PyramidP<" << dim << ">(" << this->degree << ")";

    return namebuf.str();
  }



  template <int dim, int spacedim>
  FiniteElementDomination::Domination
  FE_PyramidP<dim, spacedim>::compare_for_domination(
    const FiniteElement<dim, spacedim> &fe_other,
    const unsigned int                  codim) const
  {
    (void)fe_other;
    (void)codim;

    Assert((dynamic_cast<const Simplex::FE_P<dim, spacedim> *>(&fe_other)) ||
             (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
           ExcNotImplemented());

    return FiniteElementDomination::this_element_dominates;
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_PyramidP<dim, spacedim>::hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const
  {
    (void)fe_other;

    Assert((dynamic_cast<const Simplex::FE_P<dim, spacedim> *>(&fe_other)) ||
             (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
           ExcNotImplemented());

    return {{0, 0}};
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_PyramidP<dim, spacedim>::hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const
  {
    (void)fe_other;

    Assert((dynamic_cast<const Simplex::FE_P<dim, spacedim> *>(&fe_other)) ||
             (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
           ExcNotImplemented());

    std::vector<std::pair<unsigned int, unsigned int>> result;

    for (unsigned int i = 0; i < this->degree - 1; ++i)
      result.emplace_back(i, i);

    return result;
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_PyramidP<dim, spacedim>::hp_quad_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other,
    const unsigned int                  face_no) const
  {
    (void)fe_other;


    AssertIndexRange(face_no, 5);

    if (face_no == 0)
      {
        Assert((dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other)),
               ExcNotImplemented());
      }
    else
      {
        Assert((dynamic_cast<const Simplex::FE_P<dim, spacedim> *>(&fe_other)),
               ExcNotImplemented());
      }

    std::vector<std::pair<unsigned int, unsigned int>> result;

    for (unsigned int i = 0; i < this->n_dofs_per_quad(face_no); ++i)
      result.emplace_back(i, i);

    return result;
  }



} // namespace Simplex

// explicit instantiations
#include "fe_lib.inst"

DEAL_II_NAMESPACE_CLOSE
