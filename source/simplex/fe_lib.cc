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
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>

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
     * Set up a vector that contains the unit (reference) cell support points
     * for FE_Poly and sufficiently similar elements.
     */
    template <int dim>
    std::vector<Point<dim>>
    unit_support_points_fe_poly(const unsigned int degree)
    {
      std::vector<Point<dim>> unit_points;
      if (dim == 1)
        {
          // We don't really have dim = 1 support for simplex elements yet, but
          // its convenient for populating the face array
          Assert(degree <= 2, ExcNotImplemented());
          if (degree >= 1)
            {
              unit_points.emplace_back(0.0);
              unit_points.emplace_back(1.0);

              if (degree == 2)
                unit_points.emplace_back(0.5);
            }
        }
      else if (dim == 2)
        {
          Assert(degree <= 2, ExcNotImplemented());
          if (degree >= 1)
            {
              unit_points.emplace_back(0.0, 0.0);
              unit_points.emplace_back(1.0, 0.0);
              unit_points.emplace_back(0.0, 1.0);

              if (degree == 2)
                {
                  unit_points.emplace_back(0.5, 0.0);
                  unit_points.emplace_back(0.5, 0.5);
                  unit_points.emplace_back(0.0, 0.5);
                }
            }
        }
      else if (dim == 3)
        {
          Assert(degree <= 2, ExcNotImplemented());
          if (degree >= 1)
            {
              unit_points.emplace_back(0.0, 0.0, 0.0);
              unit_points.emplace_back(1.0, 0.0, 0.0);
              unit_points.emplace_back(0.0, 1.0, 0.0);
              unit_points.emplace_back(0.0, 0.0, 1.0);

              if (degree == 2)
                {
                  unit_points.emplace_back(0.5, 0.0, 0.0);
                  unit_points.emplace_back(0.5, 0.5, 0.0);
                  unit_points.emplace_back(0.0, 0.5, 0.0);
                  unit_points.emplace_back(0.0, 0.0, 0.5);
                  unit_points.emplace_back(0.5, 0.0, 0.5);
                  unit_points.emplace_back(0.0, 0.5, 0.5);
                }
            }
        }
      else
        {
          Assert(false, ExcNotImplemented());
        }

      return unit_points;
    }

    /**
     * Set up a vector that contains the unit (reference) cell's faces support
     * points for FE_Poly and sufficiently similar elements.
     */
    template <int dim>
    std::vector<std::vector<Point<dim - 1>>>
    unit_face_support_points_fe_poly(const unsigned int degree)
    {
      Assert(dim == 2 || dim == 3, ExcNotImplemented());
      const auto &info = ReferenceCell::internal::Info::get_cell(
        dim == 2 ? ReferenceCell::Type::Tri : ReferenceCell::Type::Tet);
      std::vector<std::vector<Point<dim - 1>>> unit_face_points;

      // all faces have the same support points
      for (auto face_n : info.face_indices())
        {
          (void)face_n;
          unit_face_points.emplace_back(
            unit_support_points_fe_poly<dim - 1>(degree));
        }

      return unit_face_points;
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
    get_dpo_vector_fe_wedge_p(const unsigned int degree)
    {
      internal::GenericDoFsPerObject dpo;

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
     * Helper function to set up the dpo vector of FE_WedgeDGP for a given @p degree.
     */
    internal::GenericDoFsPerObject
    get_dpo_vector_fe_wedge_dgp(const unsigned int degree)
    {
      unsigned int n_dofs = 0;

      if (degree == 1)
        n_dofs = 6;
      else if (degree == 2)
        n_dofs = 18;
      else
        Assert(false, ExcNotImplemented());

      return internal::expand(3,
                              {{0, 0, 0, n_dofs}},
                              ReferenceCell::Type::Wedge);
    }

    /**
     * Helper function to set up the dpo vector of FE_PyramidP for a given @p degree.
     */
    internal::GenericDoFsPerObject
    get_dpo_vector_fe_pyramid_p(const unsigned int degree)
    {
      internal::GenericDoFsPerObject dpo;

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

    /**
     * Helper function to set up the dpo vector of FE_PyramidDGP for a given @p degree.
     */
    internal::GenericDoFsPerObject
    get_dpo_vector_fe_pyramid_dgp(const unsigned int degree)
    {
      unsigned int n_dofs = 0;

      if (degree == 1)
        n_dofs = 5;
      else
        Assert(false, ExcNotImplemented());

      return internal::expand(3,
                              {{0, 0, 0, n_dofs}},
                              ReferenceCell::Type::Pyramid);
    }
  } // namespace



  template <int dim, int spacedim>
  FE_Poly<dim, spacedim>::FE_Poly(
    const unsigned int                                degree,
    const std::vector<unsigned int> &                 dpo_vector,
    const typename FiniteElementData<dim>::Conformity conformity)
    : dealii::FE_Poly<dim, spacedim>(
        Simplex::ScalarPolynomial<dim>(degree),
        FiniteElementData<dim>(dpo_vector,
                               dim == 2 ? ReferenceCell::Type::Tri :
                                          ReferenceCell::Type::Tet,
                               1,
                               degree,
                               conformity),
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
    this->unit_support_points = unit_support_points_fe_poly<dim>(degree);
    // Discontinuous elements don't have face support points
    if (conformity == FiniteElementData<dim>::Conformity::H1)
      this->unit_face_support_points =
        unit_face_support_points_fe_poly<dim>(degree);
  }



  template <int dim, int spacedim>
  std::pair<Table<2, bool>, std::vector<unsigned int>>
  FE_Poly<dim, spacedim>::get_constant_modes() const
  {
    Table<2, bool> constant_modes(1, this->n_dofs_per_cell());
    constant_modes.fill(true);
    return std::pair<Table<2, bool>, std::vector<unsigned int>>(
      constant_modes, std::vector<unsigned int>(1, 0));
  }



  template <int dim, int spacedim>
  const FullMatrix<double> &
  FE_Poly<dim, spacedim>::get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case) const
  {
    Assert(refinement_case == RefinementCase<dim>::isotropic_refinement,
           ExcNotImplemented());
    AssertDimension(dim, spacedim);

    // initialization upon first request
    if (this->prolongation[refinement_case - 1][child].n() == 0)
      {
        std::lock_guard<std::mutex> lock(this->mutex);

        // if matrix got updated while waiting for the lock
        if (this->prolongation[refinement_case - 1][child].n() ==
            this->n_dofs_per_cell())
          return this->prolongation[refinement_case - 1][child];

        // now do the work. need to get a non-const version of data in order to
        // be able to modify them inside a const function
        auto &this_nonconst = const_cast<FE_Poly<dim, spacedim> &>(*this);

        std::vector<std::vector<FullMatrix<double>>> isotropic_matrices(
          RefinementCase<dim>::isotropic_refinement);
        isotropic_matrices.back().resize(
          GeometryInfo<dim>::n_children(RefinementCase<dim>(refinement_case)),
          FullMatrix<double>(this->n_dofs_per_cell(), this->n_dofs_per_cell()));

        FETools::compute_embedding_matrices(*this, isotropic_matrices, true);

        this_nonconst.prolongation[refinement_case - 1].swap(
          isotropic_matrices.back());
      }

    // finally return the matrix
    return this->prolongation[refinement_case - 1][child];
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
    : FE_Poly<dim, spacedim>(degree,
                             get_dpo_vector_fe_p(dim, degree),
                             FiniteElementData<dim>::H1)
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
    Assert(codim <= dim, ExcImpossibleInDim(dim));

    // vertex/line/face domination
    // (if fe_other is derived from FE_DGP)
    // ------------------------------------
    if (codim > 0)
      if (dynamic_cast<const FE_DGP<dim, spacedim> *>(&fe_other) != nullptr)
        // there are no requirements between continuous and discontinuous
        // elements
        return FiniteElementDomination::no_requirements;

    // vertex/line/face domination
    // (if fe_other is not derived from FE_DGP)
    // & cell domination
    // ----------------------------------------
    if (const FE_P<dim, spacedim> *fe_p_other =
          dynamic_cast<const FE_P<dim, spacedim> *>(&fe_other))
      {
        if (this->degree < fe_p_other->degree)
          return FiniteElementDomination::this_element_dominates;
        else if (this->degree == fe_p_other->degree)
          return FiniteElementDomination::either_element_can_dominate;
        else
          return FiniteElementDomination::other_element_dominates;
      }
    else if (const FE_Q<dim, spacedim> *fe_q_other =
               dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other))
      {
        if (this->degree < fe_q_other->degree)
          return FiniteElementDomination::this_element_dominates;
        else if (this->degree == fe_q_other->degree)
          return FiniteElementDomination::either_element_can_dominate;
        else
          return FiniteElementDomination::other_element_dominates;
      }
    else if (const FE_Nothing<dim, spacedim> *fe_nothing =
               dynamic_cast<const FE_Nothing<dim, spacedim> *>(&fe_other))
      {
        if (fe_nothing->is_dominating())
          return FiniteElementDomination::other_element_dominates;
        else
          // the FE_Nothing has no degrees of freedom and it is typically used
          // in a context where we don't require any continuity along the
          // interface
          return FiniteElementDomination::no_requirements;
      }

    Assert(false, ExcNotImplemented());
    return FiniteElementDomination::neither_element_dominates;
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_P<dim, spacedim>::hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const
  {
    AssertDimension(dim, 2);

    if (dynamic_cast<const FE_P<dim, spacedim> *>(&fe_other) != nullptr)
      {
        // there should be exactly one single DoF of each FE at a vertex, and
        // they should have identical value
        return {{0U, 0U}};
      }
    else if (dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other) != nullptr)
      {
        // there should be exactly one single DoF of each FE at a vertex, and
        // they should have identical value
        return {{0U, 0U}};
      }
    else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
      {
        // the FE_Nothing has no degrees of freedom, so there are no
        // equivalencies to be recorded
        return {};
      }
    else if (fe_other.n_unique_faces() == 1 && fe_other.n_dofs_per_face(0) == 0)
      {
        // if the other element has no elements on faces at all,
        // then it would be impossible to enforce any kind of
        // continuity even if we knew exactly what kind of element
        // we have -- simply because the other element declares
        // that it is discontinuous because it has no DoFs on
        // its faces. in that case, just state that we have no
        // constraints to declare
        return {};
      }
    else
      {
        Assert(false, ExcNotImplemented());
        return {};
      }
  }



  template <int dim, int spacedim>
  std::vector<std::pair<unsigned int, unsigned int>>
  FE_P<dim, spacedim>::hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const
  {
    AssertDimension(dim, 2);
    Assert(this->degree <= 2, ExcNotImplemented());

    if (const FE_P<dim, spacedim> *fe_p_other =
          dynamic_cast<const FE_P<dim, spacedim> *>(&fe_other))
      {
        // dofs are located along lines, so two dofs are identical if they are
        // located at identical positions.
        // Therefore, read the points in unit_support_points for the
        // first coordinate direction. For FE_P, they are currently hard-coded
        // and we iterate over points on the first line which begin after the 3
        // vertex points in the complete list of unit support points

        Assert(fe_p_other->degree <= 2, ExcNotImplemented());

        std::vector<std::pair<unsigned int, unsigned int>> identities;

        for (unsigned int i = 0; i < this->degree - 1; ++i)
          for (unsigned int j = 0; j < fe_p_other->degree - 1; ++j)
            if (std::fabs(this->unit_support_points[i + 3][0] -
                          fe_p_other->unit_support_points[i + 3][0]) < 1e-14)
              identities.emplace_back(i, j);

        return identities;
      }
    else if (const FE_Q<dim, spacedim> *fe_q_other =
               dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other))
      {
        // dofs are located along lines, so two dofs are identical if they are
        // located at identical positions. if we had only equidistant points, we
        // could simply check for similarity like (i+1)*q == (j+1)*p, but we
        // might have other support points (e.g. Gauss-Lobatto
        // points). Therefore, read the points in unit_support_points for the
        // first coordinate direction. For FE_Q, we take the lexicographic
        // ordering of the line support points in the first direction (i.e.,
        // x-direction), which we access between index 1 and p-1 (index 0 and p
        // are vertex dofs). For FE_P, they are currently hard-coded and we
        // iterate over points on the first line which begin after the 3 vertex
        // points in the complete list of unit support points

        const std::vector<unsigned int> &index_map_inverse_q_other =
          fe_q_other->get_poly_space_numbering_inverse();

        std::vector<std::pair<unsigned int, unsigned int>> identities;

        for (unsigned int i = 0; i < this->degree - 1; ++i)
          for (unsigned int j = 0; j < fe_q_other->degree - 1; ++j)
            if (std::fabs(this->unit_support_points[i + 3][0] -
                          fe_q_other->get_unit_support_points()
                            [index_map_inverse_q_other[j + 1]][0]) < 1e-14)
              identities.emplace_back(i, j);

        return identities;
      }
    else if (dynamic_cast<const FE_Nothing<dim> *>(&fe_other) != nullptr)
      {
        // the FE_Nothing has no degrees of freedom, so there are no
        // equivalencies to be recorded
        return {};
      }
    else if (fe_other.n_unique_faces() == 1 && fe_other.n_dofs_per_face(0) == 0)
      {
        // if the other element has no elements on faces at all,
        // then it would be impossible to enforce any kind of
        // continuity even if we knew exactly what kind of element
        // we have -- simply because the other element declares
        // that it is discontinuous because it has no DoFs on
        // its faces. in that case, just state that we have no
        // constraints to declare
        return {};
      }
    else
      {
        Assert(false, ExcNotImplemented());
        return {};
      }
  }



  template <int dim, int spacedim>
  FE_DGP<dim, spacedim>::FE_DGP(const unsigned int degree)
    : FE_Poly<dim, spacedim>(degree,
                             get_dpo_vector_fe_dgp(dim, degree),
                             FiniteElementData<dim>::L2)
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
    if (const FE_DGP<dim, spacedim> *fe_dgp_other =
          dynamic_cast<const FE_DGP<dim, spacedim> *>(&fe_other))
      {
        if (this->degree < fe_dgp_other->degree)
          return FiniteElementDomination::this_element_dominates;
        else if (this->degree == fe_dgp_other->degree)
          return FiniteElementDomination::either_element_can_dominate;
        else
          return FiniteElementDomination::other_element_dominates;
      }
    else if (const FE_DGQ<dim, spacedim> *fe_dgq_other =
               dynamic_cast<const FE_DGQ<dim, spacedim> *>(&fe_other))
      {
        if (this->degree < fe_dgq_other->degree)
          return FiniteElementDomination::this_element_dominates;
        else if (this->degree == fe_dgq_other->degree)
          return FiniteElementDomination::either_element_can_dominate;
        else
          return FiniteElementDomination::other_element_dominates;
      }
    else if (const FE_Nothing<dim, spacedim> *fe_nothing =
               dynamic_cast<const FE_Nothing<dim, spacedim> *>(&fe_other))
      {
        if (fe_nothing->is_dominating())
          return FiniteElementDomination::other_element_dominates;
        else
          // the FE_Nothing has no degrees of freedom and it is typically used
          // in a context where we don't require any continuity along the
          // interface
          return FiniteElementDomination::no_requirements;
      }

    Assert(false, ExcNotImplemented());
    return FiniteElementDomination::neither_element_dominates;
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
  FE_Wedge<dim, spacedim>::FE_Wedge(
    const unsigned int                                degree,
    const internal::GenericDoFsPerObject &            dpos,
    const typename FiniteElementData<dim>::Conformity conformity)
    : dealii::FE_Poly<dim, spacedim>(
        Simplex::ScalarWedgePolynomial<dim>(degree),
        FiniteElementData<dim>(dpos,
                               ReferenceCell::Type::Wedge,
                               1,
                               degree,
                               conformity),
        std::vector<bool>(
          FiniteElementData<dim>(dpos, ReferenceCell::Type::Wedge, 1, degree)
            .dofs_per_cell,
          true),
        std::vector<ComponentMask>(
          FiniteElementData<dim>(dpos, ReferenceCell::Type::Wedge, 1, degree)
            .dofs_per_cell,
          std::vector<bool>(1, true)))
  {
    AssertDimension(dim, 3);

    if (degree == 1)
      {
        this->unit_support_points.emplace_back(0.0, 0.0, 0.0);
        this->unit_support_points.emplace_back(1.0, 0.0, 0.0);
        this->unit_support_points.emplace_back(0.0, 1.0, 0.0);
        this->unit_support_points.emplace_back(0.0, 0.0, 1.0);
        this->unit_support_points.emplace_back(1.0, 0.0, 1.0);
        this->unit_support_points.emplace_back(0.0, 1.0, 1.0);
      }
  }



  template <int dim, int spacedim>
  FE_WedgeP<dim, spacedim>::FE_WedgeP(const unsigned int degree)
    : FE_Wedge<dim, spacedim>(degree,
                              get_dpo_vector_fe_wedge_p(degree),
                              FiniteElementData<dim>::H1)
  {}



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
    Assert(codim <= dim, ExcImpossibleInDim(dim));

    // vertex/line/face domination
    // (if fe_other is derived from FE_DGP)
    // ------------------------------------
    if (codim > 0)
      if (dynamic_cast<const FE_DGP<dim, spacedim> *>(&fe_other) != nullptr)
        // there are no requirements between continuous and discontinuous
        // elements
        return FiniteElementDomination::no_requirements;


    // vertex/line/face domination
    // (if fe_other is not derived from FE_DGP)
    // & cell domination
    // ----------------------------------------
    if (const FE_WedgeP<dim, spacedim> *fe_wp_other =
          dynamic_cast<const FE_WedgeP<dim, spacedim> *>(&fe_other))
      {
        if (this->degree < fe_wp_other->degree)
          return FiniteElementDomination::this_element_dominates;
        else if (this->degree == fe_wp_other->degree)
          return FiniteElementDomination::either_element_can_dominate;
        else
          return FiniteElementDomination::other_element_dominates;
      }
    else if (const FE_P<dim, spacedim> *fe_p_other =
               dynamic_cast<const FE_P<dim, spacedim> *>(&fe_other))
      {
        if (this->degree < fe_p_other->degree)
          return FiniteElementDomination::this_element_dominates;
        else if (this->degree == fe_p_other->degree)
          return FiniteElementDomination::either_element_can_dominate;
        else
          return FiniteElementDomination::other_element_dominates;
      }
    else if (const FE_Q<dim, spacedim> *fe_q_other =
               dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other))
      {
        if (this->degree < fe_q_other->degree)
          return FiniteElementDomination::this_element_dominates;
        else if (this->degree == fe_q_other->degree)
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

    Assert(false, ExcNotImplemented());
    return FiniteElementDomination::neither_element_dominates;
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
  FE_WedgeDGP<dim, spacedim>::FE_WedgeDGP(const unsigned int degree)
    : FE_Wedge<dim, spacedim>(degree,
                              get_dpo_vector_fe_wedge_dgp(degree),
                              FiniteElementData<dim>::L2)
  {}



  template <int dim, int spacedim>
  std::unique_ptr<FiniteElement<dim, spacedim>>
  FE_WedgeDGP<dim, spacedim>::clone() const
  {
    return std::make_unique<FE_WedgeDGP<dim, spacedim>>(*this);
  }



  template <int dim, int spacedim>
  std::string
  FE_WedgeDGP<dim, spacedim>::get_name() const
  {
    std::ostringstream namebuf;
    namebuf << "FE_WedgeDGP<" << dim << ">(" << this->degree << ")";

    return namebuf.str();
  }



  template <int dim, int spacedim>
  FE_Pyramid<dim, spacedim>::FE_Pyramid(
    const unsigned int                                degree,
    const internal::GenericDoFsPerObject &            dpos,
    const typename FiniteElementData<dim>::Conformity conformity)
    : dealii::FE_Poly<dim, spacedim>(
        Simplex::ScalarPyramidPolynomial<dim>(degree),
        FiniteElementData<dim>(dpos,
                               ReferenceCell::Type::Pyramid,
                               1,
                               degree,
                               conformity),
        std::vector<bool>(
          FiniteElementData<dim>(dpos, ReferenceCell::Type::Pyramid, 1, degree)
            .dofs_per_cell,
          true),
        std::vector<ComponentMask>(
          FiniteElementData<dim>(dpos, ReferenceCell::Type::Pyramid, 1, degree)
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
  FE_PyramidP<dim, spacedim>::FE_PyramidP(const unsigned int degree)
    : FE_Pyramid<dim, spacedim>(degree,
                                get_dpo_vector_fe_pyramid_p(degree),
                                FiniteElementData<dim>::H1)
  {}



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
    Assert(codim <= dim, ExcImpossibleInDim(dim));

    // vertex/line/face domination
    // (if fe_other is derived from FE_DGP)
    // ------------------------------------
    if (codim > 0)
      if (dynamic_cast<const FE_DGP<dim, spacedim> *>(&fe_other) != nullptr)
        // there are no requirements between continuous and discontinuous
        // elements
        return FiniteElementDomination::no_requirements;

    // vertex/line/face domination
    // (if fe_other is not derived from FE_DGP)
    // & cell domination
    // ----------------------------------------
    if (const FE_PyramidP<dim, spacedim> *fe_pp_other =
          dynamic_cast<const FE_PyramidP<dim, spacedim> *>(&fe_other))
      {
        if (this->degree < fe_pp_other->degree)
          return FiniteElementDomination::this_element_dominates;
        else if (this->degree == fe_pp_other->degree)
          return FiniteElementDomination::either_element_can_dominate;
        else
          return FiniteElementDomination::other_element_dominates;
      }
    else if (const FE_P<dim, spacedim> *fe_p_other =
               dynamic_cast<const FE_P<dim, spacedim> *>(&fe_other))
      {
        if (this->degree < fe_p_other->degree)
          return FiniteElementDomination::this_element_dominates;
        else if (this->degree == fe_p_other->degree)
          return FiniteElementDomination::either_element_can_dominate;
        else
          return FiniteElementDomination::other_element_dominates;
      }
    else if (const FE_Q<dim, spacedim> *fe_q_other =
               dynamic_cast<const FE_Q<dim, spacedim> *>(&fe_other))
      {
        if (this->degree < fe_q_other->degree)
          return FiniteElementDomination::this_element_dominates;
        else if (this->degree == fe_q_other->degree)
          return FiniteElementDomination::either_element_can_dominate;
        else
          return FiniteElementDomination::other_element_dominates;
      }
    else if (const FE_Nothing<dim, spacedim> *fe_nothing =
               dynamic_cast<const FE_Nothing<dim, spacedim> *>(&fe_other))
      {
        if (fe_nothing->is_dominating())
          return FiniteElementDomination::other_element_dominates;
        else
          // the FE_Nothing has no degrees of freedom and it is typically used
          // in a context where we don't require any continuity along the
          // interface
          return FiniteElementDomination::no_requirements;
      }

    Assert(false, ExcNotImplemented());
    return FiniteElementDomination::neither_element_dominates;
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



  template <int dim, int spacedim>
  FE_PyramidDGP<dim, spacedim>::FE_PyramidDGP(const unsigned int degree)
    : FE_Pyramid<dim, spacedim>(degree,
                                get_dpo_vector_fe_pyramid_dgp(degree),
                                FiniteElementData<dim>::L2)
  {}



  template <int dim, int spacedim>
  std::unique_ptr<FiniteElement<dim, spacedim>>
  FE_PyramidDGP<dim, spacedim>::clone() const
  {
    return std::make_unique<FE_PyramidDGP<dim, spacedim>>(*this);
  }



  template <int dim, int spacedim>
  std::string
  FE_PyramidDGP<dim, spacedim>::get_name() const
  {
    std::ostringstream namebuf;
    namebuf << "FE_PyramidDGP<" << dim << ">(" << this->degree << ")";

    return namebuf.str();
  }



} // namespace Simplex

// explicit instantiations
#include "fe_lib.inst"

DEAL_II_NAMESPACE_CLOSE
