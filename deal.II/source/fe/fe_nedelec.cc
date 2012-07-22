#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <sstream>
#include <iostream>

//TODO: implement the adjust_quad_dof_index_for_face_orientation_table and
//adjust_line_dof_index_for_line_orientation_table fields, and write tests
//similar to bits/face_orientation_and_fe_q_*


DEAL_II_NAMESPACE_OPEN

//#define DEBUG_NEDELEC


template <int dim>
FE_Nedelec<dim>::FE_Nedelec (const unsigned int p) :
FE_PolyTensor<PolynomialsNedelec<dim>, dim>
(p,
 FiniteElementData<dim> (get_dpo_vector (p), dim, p + 1,
                         FiniteElementData<dim>::Hcurl, 1),
 std::vector<bool> (PolynomialsNedelec<dim>::compute_n_pols (p), true),
 std::vector<std::vector<bool> >
 (PolynomialsNedelec<dim>::compute_n_pols (p),
  std::vector<bool> (dim, true)))
{
#ifdef DEBUG_NEDELEC
  deallog << get_name() << std::endl;
#endif

  Assert (dim >= 2, ExcImpossibleInDim(dim));

  const unsigned int n_dofs = this->dofs_per_cell;

  this->mapping_type = mapping_nedelec;
                                   // First, initialize the
                                   // generalized support points and
                                   // quadrature weights, since they
                                   // are required for interpolation.
  initialize_support_points (p);
  this->inverse_node_matrix.reinit (n_dofs, n_dofs);
  this->inverse_node_matrix.fill
    (FullMatrix<double> (IdentityMatrix (n_dofs)));
                                   // From now on, the shape functions
                                   // will be the correct ones, not
                                   // the raw shape functions anymore.

                                   // Reinit the vectors of
                                   // restriction and prolongation
                                   // matrices to the right sizes.
                                   // Restriction only for isotropic
                                   // refinement
#ifdef DEBUG_NEDELEC
  deallog << "Embedding" << std::endl;
#endif
  this->reinit_restriction_and_prolongation_matrices ();
                                   // Fill prolongation matrices with embedding operators
  FETools::compute_embedding_matrices (*this, this->prolongation, true);
#ifdef DEBUG_NEDELEC
  deallog << "Restriction" << std::endl;
#endif
  initialize_restriction ();

#ifdef DEBUG_NEDELEC
  deallog << "Face Embedding" << std::endl;
#endif
  FullMatrix<double> face_embeddings[GeometryInfo<dim>::max_children_per_face];

  for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_face; ++i)
     face_embeddings[i].reinit (this->dofs_per_face, this->dofs_per_face);

  FETools::compute_face_embedding_matrices<dim,double>
    (*this, face_embeddings, 0, 0);

  switch (dim)
    {
      case 1:
        {
          this->interface_constraints.reinit (0, 0);
          break;
        }

      case 2:
        {
          this->interface_constraints.reinit (2 * this->dofs_per_face,
                                              this->dofs_per_face);

          for (unsigned int i = 0; i < GeometryInfo<2>::max_children_per_face;
               ++i)
            for (unsigned int j = 0; j < this->dofs_per_face; ++j)
              for (unsigned int k = 0; k < this->dofs_per_face; ++k)
                this->interface_constraints (i * this->dofs_per_face + j, k)
                  = face_embeddings[i] (j, k);

          break;
        }

      case 3:
        {
          this->interface_constraints.reinit
          (4 * (this->dofs_per_face - this->degree), this->dofs_per_face);

          unsigned int target_row = 0;

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = this->degree; j < 2 * this->degree;
                 ++j, ++target_row)
              for (unsigned int k = 0; k < this->dofs_per_face; ++k)
                this->interface_constraints (target_row, k)
                  = face_embeddings[2 * i] (j, k);

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 3 * this->degree;
                 j < GeometryInfo<3>::lines_per_face * this->degree;
                 ++j, ++target_row)
              for (unsigned int k = 0; k < this->dofs_per_face; ++k)
                this->interface_constraints (target_row, k)
                  = face_embeddings[i] (j, k);

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              for (unsigned int k = i * this->degree;
                   k < (i + 1) * this->degree; ++k, ++target_row)
                for (unsigned int l = 0; l < this->dofs_per_face; ++l)
                  this->interface_constraints (target_row, l)
                    = face_embeddings[i + 2 * j] (k, l);

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              for (unsigned int k = (i + 2) * this->degree;
                   k < (i + 3) * this->degree; ++k, ++target_row)
                for (unsigned int l = 0; l < this->dofs_per_face; ++l)
                  this->interface_constraints (target_row, l)
                    = face_embeddings[2 * i + j] (k, l);

          for (unsigned int i = 0; i < GeometryInfo<3>::max_children_per_face;
               ++i)
            for (unsigned int
                   j = GeometryInfo<3>::lines_per_face * this->degree;
                 j < this->dofs_per_face; ++j, ++target_row)
              for (unsigned int k = 0; k < this->dofs_per_face; ++k)
                this->interface_constraints (target_row, k)
                  = face_embeddings[i] (j, k);

          break;
        }

      default:
        Assert (false, ExcNotImplemented ());
    }
}



template <int dim>
std::string
FE_Nedelec<dim>::get_name () const
{
                                   // note that the
                                   // FETools::get_fe_from_name
                                   // function depends on the
                                   // particular format of the string
                                   // this function returns, so they
                                   // have to be kept in synch

  std::ostringstream namebuf;
  namebuf << "FE_Nedelec<" << dim << ">(" << this->degree-1 << ")";

  return namebuf.str();
}


template <int dim>
FiniteElement<dim>
*FE_Nedelec<dim>::clone () const
{
  return new FE_Nedelec<dim> (*this);
}

//---------------------------------------------------------------------------
// Auxiliary and internal functions
//---------------------------------------------------------------------------



                   // Set the generalized support
                   // points and precompute the
                   // parts of the projection-based
                   // interpolation, which does
                   // not depend on the interpolated
                   // function.
template <>
void
FE_Nedelec<1>::initialize_support_points (const unsigned int)
{
  Assert (false, ExcNotImplemented ());
}



template <>
void
FE_Nedelec<2>::initialize_support_points (const unsigned int degree)
{
  const int dim = 2;

                       // Create polynomial basis.
  const std::vector<Polynomials::Polynomial<double> >& lobatto_polynomials
    = Polynomials::Lobatto::generate_complete_basis (degree + 1);
  std::vector<Polynomials::Polynomial<double> >
    lobatto_polynomials_grad (degree + 1);

  for (unsigned int i = 0; i < lobatto_polynomials_grad.size (); ++i)
    lobatto_polynomials_grad[i] = lobatto_polynomials[i + 1].derivative ();

                                   // Initialize quadratures to obtain
                                   // quadrature points later on.
  const QGauss<dim - 1> reference_edge_quadrature (degree + 1);
  const unsigned int&
    n_edge_points = reference_edge_quadrature.size ();
  const unsigned int n_boundary_points
    = GeometryInfo<dim>::lines_per_cell * n_edge_points;
  const Quadrature<dim>& edge_quadrature
    = QProjector<dim>::project_to_all_faces (reference_edge_quadrature);

  this->generalized_face_support_points.resize (n_edge_points);

                                   // Create face support points.
  for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
    this->generalized_face_support_points[q_point]
      = reference_edge_quadrature.point (q_point);

  if (degree > 0)
    {
                                       // If the polynomial degree is positive
                                       // we have support points on the faces
                                       // and in the interior of a cell.
      const QGauss<dim> quadrature (degree + 1);
      const unsigned int& n_interior_points = quadrature.size ();

      this->generalized_support_points.resize
        (n_boundary_points + n_interior_points);
      boundary_weights.reinit (n_edge_points, degree);

      for (unsigned int q_point = 0; q_point < n_edge_points;
           ++q_point)
        {
          for (unsigned int line = 0;
               line < GeometryInfo<dim>::lines_per_cell; ++line)
            this->generalized_support_points[line * n_edge_points
                                             + q_point]
              = edge_quadrature.point
              (QProjector<dim>::DataSetDescriptor::face
               (line, true, false, false, n_edge_points) + q_point);

          for (unsigned int i = 0; i < degree; ++i)
            boundary_weights (q_point, i)
              = reference_edge_quadrature.weight (q_point)
              * lobatto_polynomials_grad[i + 1].value
              (this->generalized_face_support_points[q_point] (0));
        }

      for (unsigned int q_point = 0; q_point < n_interior_points;
           ++q_point)
        this->generalized_support_points[q_point + n_boundary_points]
          = quadrature.point (q_point);
    }

  else
    {
                                       // In this case we only need support points
                                       // on the faces of a cell.
      const Quadrature<dim>& edge_quadrature
        = QProjector<dim>::project_to_all_faces
        (reference_edge_quadrature);

      this->generalized_support_points.resize (n_boundary_points);

      for (unsigned int line = 0;
           line < GeometryInfo<dim>::lines_per_cell; ++line)
        for (unsigned int q_point = 0; q_point < n_edge_points;
             ++q_point)
          this->generalized_support_points[line * n_edge_points
                                           + q_point]
            = edge_quadrature.point
            (QProjector<dim>::DataSetDescriptor::face
             (line, true, false, false, n_edge_points) + q_point);
    }
}



template <>
void
FE_Nedelec<3>::initialize_support_points (const unsigned int degree)
{
  const int dim = 3;

                       // Create polynomial basis.
  const std::vector<Polynomials::Polynomial<double> >& lobatto_polynomials
    = Polynomials::Lobatto::generate_complete_basis (degree + 1);
  std::vector<Polynomials::Polynomial<double> >
    lobatto_polynomials_grad (degree + 1);

  for (unsigned int i = 0; i < lobatto_polynomials_grad.size (); ++i)
    lobatto_polynomials_grad[i] = lobatto_polynomials[i + 1].derivative ();

                                   // Initialize quadratures to obtain
                                   // quadrature points later on.
  const QGauss<1> reference_edge_quadrature (degree + 1);
  const unsigned int& n_edge_points = reference_edge_quadrature.size ();
  const Quadrature<dim - 1>& edge_quadrature
    = QProjector<dim - 1>::project_to_all_faces
    (reference_edge_quadrature);

  if (degree > 0)
    {
                                       // If the polynomial degree is positive
                                       // we have support points on the edges,
                                       // faces and in the interior of a cell.
      const QGauss<dim - 1> reference_face_quadrature (degree + 1);
      const unsigned int& n_face_points
        = reference_face_quadrature.size ();
      const unsigned int n_boundary_points
        = GeometryInfo<dim>::lines_per_cell * n_edge_points
        + GeometryInfo<dim>::faces_per_cell * n_face_points;
      const QGauss<dim> quadrature (degree + 1);
      const unsigned int& n_interior_points = quadrature.size ();

      boundary_weights.reinit (n_edge_points + n_face_points,
                               2 * (degree + 1) * degree);
      this->generalized_face_support_points.resize
        (4 * n_edge_points + n_face_points);
      this->generalized_support_points.resize
        (n_boundary_points + n_interior_points);

                                       // Create support points on edges.
      for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
        {
          for (unsigned int line = 0;
               line < GeometryInfo<dim - 1>::lines_per_cell; ++line)
            this->generalized_face_support_points[line * n_edge_points
                                                  + q_point]
              = edge_quadrature.point
              (QProjector<dim - 1>::DataSetDescriptor::face
               (line, true, false, false, n_edge_points) + q_point);

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              {
                this->generalized_support_points
                  [q_point + (i + 4 * j) * n_edge_points]
                  = Point<dim>
                  (i, reference_edge_quadrature.point (q_point) (0),
                   j);
                this->generalized_support_points
                  [q_point + (i + 4 * j + 2) * n_edge_points]
                  = Point<dim>
                  (reference_edge_quadrature.point (q_point) (0),
                   i, j);
                this->generalized_support_points
                  [q_point + (i + 2 * (j + 4)) * n_edge_points]
                  = Point<dim>
                  (i, j,
                   reference_edge_quadrature.point (q_point) (0));
              }

          for (unsigned int i = 0; i < degree; ++i)
            boundary_weights (q_point, i)
              = reference_edge_quadrature.weight (q_point)
              * lobatto_polynomials_grad[i + 1].value
              (this->generalized_face_support_points[q_point] (1));
        }

                                       // Create support points on faces.
      for (unsigned int q_point = 0; q_point < n_face_points;
           ++q_point)
        {
          this->generalized_face_support_points[q_point
                                                + 4 * n_edge_points]
            = reference_face_quadrature.point (q_point);

          for (unsigned int i = 0; i <= degree; ++i)
            for (unsigned int j = 0; j < degree; ++j)
              {
                boundary_weights (q_point + n_edge_points,
                                  2 * (i * degree + j))
                  = reference_face_quadrature.weight (q_point)
                  * lobatto_polynomials_grad[i].value
                  (this->generalized_face_support_points
                   [q_point + 4 * n_edge_points] (0))
                  * lobatto_polynomials[j + 2].value
                  (this->generalized_face_support_points
                   [q_point + 4 * n_edge_points] (1));
                boundary_weights (q_point + n_edge_points,
                                  2 * (i * degree + j) + 1)
                  = reference_face_quadrature.weight (q_point)
                  * lobatto_polynomials_grad[i].value
                  (this->generalized_face_support_points
                   [q_point + 4 * n_edge_points] (1))
                  * lobatto_polynomials[j + 2].value
                  (this->generalized_face_support_points
                   [q_point + 4 * n_edge_points] (0));
              }
        }

      const Quadrature<dim>& face_quadrature
        = QProjector<dim>::project_to_all_faces
        (reference_face_quadrature);

      for (unsigned int face = 0;
           face < GeometryInfo<dim>::faces_per_cell; ++face)
        for (unsigned int q_point = 0; q_point < n_face_points;
             ++q_point)
          {
            this->generalized_support_points
              [face * n_face_points + q_point
               + GeometryInfo<dim>::lines_per_cell * n_edge_points]
              = face_quadrature.point
              (QProjector<dim>::DataSetDescriptor::face
               (face, true, false, false, n_face_points) + q_point);
          }

                                       // Create support points in the interior.
      for (unsigned int q_point = 0; q_point < n_interior_points;
           ++q_point)
        this->generalized_support_points[q_point + n_boundary_points]
          = quadrature.point (q_point);
    }

  else
    {
      this->generalized_face_support_points.resize (4 * n_edge_points);
      this->generalized_support_points.resize
        (GeometryInfo<dim>::lines_per_cell * n_edge_points);

      for (unsigned int q_point = 0; q_point < n_edge_points;
           ++q_point)
        {
          for (unsigned int line = 0;
               line < GeometryInfo<dim - 1>::lines_per_cell; ++line)
            this->generalized_face_support_points[line * n_edge_points
                                                  + q_point]
              = edge_quadrature.point
              (QProjector<dim - 1>::DataSetDescriptor::face
               (line, true, false, false, n_edge_points) + q_point);

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              {
                this->generalized_support_points
                  [q_point + (i + 4 * j) * n_edge_points]
                  = Point<dim>
                  (i, reference_edge_quadrature.point (q_point) (0),
                   j);
                this->generalized_support_points
                  [q_point + (i + 4 * j + 2) * n_edge_points]
                  = Point<dim>
                  (reference_edge_quadrature.point (q_point) (0),
                   i, j);
                this->generalized_support_points
                  [q_point + (i + 2 * (j + 4)) * n_edge_points]
                  = Point<dim>
                  (i, j,
                   reference_edge_quadrature.point (q_point) (0));
              }
        }
    }
}



                   // Set the restriction matrices.
template <>
void
FE_Nedelec<1>::initialize_restriction ()
{
                                   // there is only one refinement case in 1d,
                                   // which is the isotropic one
  for (unsigned int i = 0; i < GeometryInfo<1>::max_children_per_cell; ++i)
    this->restriction[0][i].reinit(0, 0);
}





                   // Restriction operator
template <int dim>
void
FE_Nedelec<dim>::initialize_restriction ()
{
                                   // This function does the same as the
                                   // function interpolate further below.
                                   // But since the functions, which we
                                   // interpolate here, are discontinuous
                                   // we have to use more quadrature
                                   // points as in interpolate.
  const QGauss<1> edge_quadrature (2 * this->degree);
  const std::vector<Point<1> >& edge_quadrature_points
    = edge_quadrature.get_points ();
  const unsigned int&
    n_edge_quadrature_points = edge_quadrature.size ();
  const unsigned int
    index = RefinementCase<dim>::isotropic_refinement - 1;
  const unsigned int deg = this->degree-1;
  const std::vector<Polynomials::Polynomial<double> >&
    legendre_polynomials = Polynomials::Legendre::generate_complete_basis (deg);

  switch (dim)
    {
      case 2:
        {
                        // First interpolate the shape
                        // functions of the child cells
                        // to the face shape functions
                        // of the parent cell.
          for (unsigned int dof = 0; dof < this->dofs_per_cell; ++dof)
            for (unsigned int q_point = 0; q_point < n_edge_quadrature_points;
                 ++q_point)
              {
                if (edge_quadrature_points[q_point] (0) < 0.5)
                  for (unsigned int i = 0; i <= deg; ++i)
                    {
                      const double weight
                        = 2.0 * edge_quadrature.weight (q_point)
                              * legendre_polynomials[i].value (edge_quadrature_points[q_point] (0));
                      Point<dim> quadrature_point (0.0,
                                                   2.0 * edge_quadrature_points[q_point] (0));

                      this->restriction[index][0] (i, dof) += weight
                                                           * this->shape_value_component
                                                             (dof,
                                                              quadrature_point,
                                                              1);
                      quadrature_point (0) = 1.0;
                      this->restriction[index][1] (i + this->degree, dof)
                        += weight * this->shape_value_component (dof,
                                                                 quadrature_point,
                                                                 1);
                      quadrature_point (0) = quadrature_point (1);
                      quadrature_point (1) = 0.0;
                      this->restriction[index][0] (i + 2 * this->degree, dof)
                        += weight * this->shape_value_component (dof,
                                                                 quadrature_point,
                                                                 0);
                      quadrature_point (1) = 1.0;
                      this->restriction[index][2] (i + 3 * this->degree, dof)
                        += weight * this->shape_value_component (dof,
                                                                 quadrature_point,
                                                                 0);
                    }

                else
                  for (unsigned int i = 0; i <= deg; ++i)
                    {
                      const double weight
                        = 2.0 * edge_quadrature.weight (q_point)
                              * legendre_polynomials[i].value (edge_quadrature_points[q_point] (0));
                      Point<dim> quadrature_point (0.0,
                                                   2.0 * edge_quadrature_points[q_point] (0)
                                                       - 1.0);

                      this->restriction[index][2] (i, dof) += weight
                                                           * this->shape_value_component
                                                             (dof,
                                                              quadrature_point,
                                                              1);
                      quadrature_point (0) = 1.0;
                      this->restriction[index][3] (i + this->degree, dof)
                        += weight * this->shape_value_component (dof,
                                                                 quadrature_point,
                                                                 1);
                      quadrature_point (0) = quadrature_point (1);
                      quadrature_point (1) = 0.0;
                      this->restriction[index][1] (i + 2 * this->degree, dof)
                        += weight * this->shape_value_component (dof,
                                                                 quadrature_point,
                                                                 0);
                      quadrature_point (1) = 1.0;
                      this->restriction[index][3] (i + 3 * this->degree, dof)
                        += weight * this->shape_value_component (dof,
                                                                 quadrature_point,
                                                                 0);
                    }
              }

                        // Then interpolate the shape functions
                        // of the child cells to the interior
                        // shape functions of the parent cell.
          if (this->degree > 1)
            {
              const QGauss<dim> quadrature (2 * this->degree);
              const std::vector<Point<dim> >&
                quadrature_points = quadrature.get_points ();
              const std::vector<Polynomials::Polynomial<double> >&
                lobatto_polynomials
                  = Polynomials::Lobatto::generate_complete_basis
                    (this->degree);
              const unsigned int n_boundary_dofs
                = GeometryInfo<dim>::faces_per_cell * this->degree;

              for (unsigned int dof = 0; dof < this->dofs_per_cell; ++dof)
                for (unsigned int q_point = 0; q_point < quadrature.size ();
                     ++q_point)
                  {
                    const double weight = 2.0 * quadrature.weight (q_point);
                    Point<dim> quadrature_point (2.0 * quadrature_points[q_point] (0),
                                                 2.0 * quadrature_points[q_point] (1));
                    unsigned int k;

                    if (quadrature_points[q_point] (0) < 0.5)
                      {
                        if (quadrature_points[q_point] (1) < 0.5)
                          k = 0;
                        
                        else
                          {
                            quadrature_point (1) -= 1.0;
                            k = 1;
                          }
                      }
                    
                    else
                      if (quadrature_points[q_point] (1) < 0.5)
                        {
                          quadrature_point (0) -= 1.0;
                          k = 2;
                        }
                      
                      else
                        {
                          quadrature_point (0) -= 1.0;
                          quadrature_point (1) -= 1.0;
                          k = 3;
                        }
                                  
                    for (unsigned int i = 0; i < this->degree; ++i)
                      {
                        const double L_i_0 = weight
                                             * legendre_polynomials[i].value
                                               (quadrature_points[q_point] (0));
                        const double L_i_1 = weight
                                             * legendre_polynomials[i].value
                                               (quadrature_points[q_point] (1));

                        for (unsigned int j = 0; j < deg; ++j)
                          {
                            this->restriction[index][k] (i * deg + j + n_boundary_dofs, dof)
                              += L_i_0 * (this->shape_grad_component (dof, quadrature_point, 0)[1]
                                          * legendre_polynomials[j + 1].value
                                            (quadrature_points[q_point] (1))
                                          + this->shape_value_component (dof, quadrature_point, 0)
                                          * lobatto_polynomials[j + 2].value
                                            (quadrature_points[q_point] (1)));
                            this->restriction[index][k] (i + (deg + j) * this->degree
                                                           + n_boundary_dofs, dof)
                              += L_i_1 * (this->shape_grad_component (dof, quadrature_point, 1)[0]
                                          * legendre_polynomials[j + 1].value
                                            (quadrature_points[q_point] (0))
                                          + this->shape_value_component (dof, quadrature_point, 1)
                                          * lobatto_polynomials[j + 2].value
                                            (quadrature_points[q_point] (0)));
                          }
                      }
                  }
            }

          break;
        }

      case 3:
        {
                        // First interpolate the shape
                        // functions of the child cells
                        // to the edge shape functions
                        // of the parent cell.
          for (unsigned int dof = 0; dof < this->dofs_per_cell; ++dof)
            for (unsigned int q_point = 0; q_point < n_edge_quadrature_points;
                 ++q_point)
              {
                const double weight = 2.0 * edge_quadrature.weight (q_point);
                
                if (edge_quadrature_points[q_point] (0) < 0.5)
                  for (unsigned int i = 0; i < this->degree; ++i)
                    {
                      const double L_i
                        = weight * legendre_polynomials[i].value
                                   (edge_quadrature_points[q_point] (0));

                      for (unsigned int j = 0; j < 2; ++j)
                        for (unsigned int k = 0; k < 2; ++k)
                          {
                            Point<dim> quadrature_point (j,
                                                         2.0 * edge_quadrature_points[q_point] (0),
                                                         k);

                            this->restriction[index][j + 4 * k]
                            (i + (j + 4 * k) * this->degree, dof)
                              += L_i * this->shape_value_component (dof,
                                                                    quadrature_point,
                                                                    1);
                            quadrature_point
                              = Point<dim> (2.0 * edge_quadrature_points[q_point] (0),
                                            j, k);
                            this->restriction[index][2 * (j + 2 * k)]
                            (i + (j + 4 * k + 2) * this->degree, dof)
                              += L_i * this->shape_value_component (dof,
                                                                    quadrature_point,
                                                                    0);
                            quadrature_point = Point<dim> (j, k,
                                                           2.0 * edge_quadrature_points[q_point] (0));
                            this->restriction[index][j + 2 * k]
                            (i + (j + 2 * (k + 4)) * this->degree, dof)
                              += L_i * this->shape_value_component (dof,
                                                                    quadrature_point,
                                                                    2);
                          }
                    }
                  
                else
                  for (unsigned int i = 0; i < this->degree; ++i)
                    {
                      const double L_i
                        = weight * legendre_polynomials[i].value
                                   (edge_quadrature_points[q_point] (0));

                      for (unsigned int j = 0; j < 2; ++j)
                        for (unsigned int k = 0; k < 2; ++k)
                          {
                            Point<dim> quadrature_point (j,
                                                         2.0 * edge_quadrature_points[q_point] (0)
                                                             - 1.0, k);

                            this->restriction[index][j + 4 * k + 2]
                            (i + (j + 4 * k) * this->degree, dof)
                              += L_i * this->shape_value_component (dof,
                                                                    quadrature_point,
                                                                    1);
                            quadrature_point
                              = Point<dim> (2.0 * edge_quadrature_points[q_point] (0) - 1.0,
                                            j, k);
                            this->restriction[index][2 * (j + 2 * k) + 1]
                            (i + (j + 4 * k + 2) * this->degree, dof)
                              += L_i * this->shape_value_component (dof,
                                                                    quadrature_point,
                                                                    0);
                            quadrature_point = Point<dim> (j, k,
                                                           2.0 * edge_quadrature_points[q_point] (0)
                                                               - 1.0);
                            this->restriction[index][j + 2 * (k + 2)]
                            (i + (j + 2 * (k + 4)) * this->degree, dof)
                              += L_i * this->shape_value_component (dof,
                                                                    quadrature_point,
                                                                    2);
                          }
                    }
              }

                        // Then interpolate the shape functions
                        // of the child cells to the face
                        // and interior shape functions of
                        // the parent cell.
          if (this->degree > 1)
            {
              const QGauss<2> face_quadrature (2 * this->degree);
              const std::vector<Point<2> >& face_quadrature_points
                = face_quadrature.get_points ();
              const std::vector<Polynomials::Polynomial<double> >&
                lobatto_polynomials
                  = Polynomials::Lobatto::generate_complete_basis
                    (this->degree);
              const QGauss<dim> quadrature (2 * this->degree);
              const std::vector<Point<dim> >&
                quadrature_points = quadrature.get_points ();
              const unsigned int n_edge_dofs
                = GeometryInfo<dim>::lines_per_cell * this->degree;
              const unsigned int n_boundary_dofs
                = 2 * GeometryInfo<dim>::faces_per_cell * deg * this->degree
                    + n_edge_dofs;
              const unsigned int& n_face_quadrature_points
                = face_quadrature.size ();
              const unsigned int& n_quadrature_points = quadrature.size ();
              
                            // First, the interpolation to
                            // the face shape functions.
              for (unsigned int dof = 0; dof < this->dofs_per_cell; ++dof)
                {
                  for (unsigned int q_point = 0; q_point < n_face_quadrature_points;
                       ++q_point)
                    {
                      const double weight = 2.0 * face_quadrature.weight (q_point);
                      
                      for (unsigned int i = 0; i < 2; ++i)
                        {
                          Point<dim> quadrature_point (i,
                                                       2.0 * face_quadrature_points[q_point] (0),
                                                       2.0 * face_quadrature_points[q_point] (1));
                          unsigned int l;
                          unsigned int m;
                          
                          if (face_quadrature_points[q_point] (0) < 0.5)
                            {
                              m = 0;
                              
                              if (face_quadrature_points[q_point] (1) < 0.5)
                                l = 0;
                              
                              else
                                {
                                  quadrature_point (2) -= 1.0;
                                  l = 1;
                                }
                            }
                        
                          else
                            {
                              quadrature_point (1) -= 1.0;
                              m = 1;
                              
                              if (face_quadrature_points[q_point] (1) < 0.5)
                                l = 0;
                              
                              else
                                {
                                  quadrature_point (2) -= 1.0;
                                  l = 1;
                                }
                            }
                          
                          for (unsigned int j = 0; j < this->degree; ++j)
                            {
                              const double L_j_0
                                = weight * legendre_polynomials[j].value
                                           (face_quadrature_points[q_point] (0));
                              const double L_j_1
                                = weight * legendre_polynomials[j].value
                                           (face_quadrature_points[q_point] (1));
                            
                              for (unsigned int k = 0; k < deg; ++k)
                                {
                                  const double Le_k_0
                                    = legendre_polynomials[k + 1].value
                                      (face_quadrature_points[q_point] (0));
                                  const double Le_k_1
                                    = legendre_polynomials[k + 1].value
                                      (face_quadrature_points[q_point] (1));
                                  const double lo_k_0
                                    = lobatto_polynomials[k + 2].value
                                      (face_quadrature_points[q_point] (0));
                                  const double lo_k_1
                                    = lobatto_polynomials[k + 2].value
                                      (face_quadrature_points[q_point] (1));
                                
                                  this->restriction[index][i + 2 * (2 * l + m)]
                                  ((2 * i * this->degree + j) * deg + k + n_edge_dofs, dof)
                                    += L_j_0 * (this->shape_grad_component (dof, quadrature_point, 1)[2]
                                                * Le_k_1
                                                + this->shape_value_component (dof, quadrature_point, 1)
                                                * lo_k_1);
                                  this->restriction[index][i + 2 * (2 * l + m)]
                                  (((2 * i + 1) * deg + k) * this->degree + j + n_edge_dofs, dof)
                                    += L_j_1 * (this->shape_grad_component (dof, quadrature_point, 2)[1]
                                                * Le_k_0
                                                + this->shape_value_component (dof, quadrature_point, 2)
                                                * lo_k_0);
                                  this->restriction[index][2 * (i + 2 * l) + m]
                                  ((2 * (i + 2) * this->degree + j) * deg + k + n_edge_dofs, dof)
                                    += L_j_0 * (this->shape_grad_component (dof, quadrature_point, 2)[0]
                                                * Le_k_1
                                                + this->shape_value_component (dof, quadrature_point, 2)
                                                * lo_k_1);
                                  this->restriction[index][2 * (i + 2 * l) + m]
                                  (((2 * i + 5) * deg + k) * this->degree + j + n_edge_dofs, dof)
                                    += L_j_1 * (this->shape_grad_component (dof, quadrature_point, 0)[2]
                                                * Le_k_0
                                                + this->shape_value_component (dof, quadrature_point, 0)
                                                * lo_k_0);
                                  this->restriction[index][2 * (2 * i + l) + m]
                                  ((2 * (i + 4) * this->degree + j) * deg + k + n_edge_dofs, dof)
                                    += L_j_0 * (this->shape_grad_component (dof, quadrature_point, 0)[1]
                                                * Le_k_1
                                                + this->shape_value_component (dof, quadrature_point, 0)
                                                * lo_k_1);
                                  this->restriction[index][2 * (2 * i + l) + m]
                                  (((2 * i + 9) * deg + k) * this->degree + j + n_edge_dofs, dof)
                                    += L_j_1 * (this->shape_grad_component (dof, quadrature_point, 1)[0]
                                                * Le_k_0
                                                + this->shape_value_component (dof, quadrature_point, 1)
                                                * lo_k_0);
                                }
                            }
                        }
                    }
                  
                            // Then, the interpolation to
                            // the interior shape functions.
                  for (unsigned int q_point = 0; q_point < n_quadrature_points;
                       ++q_point)
                    {
                      Point<dim> quadrature_point = 2.0 * quadrature_points[q_point];
                      unsigned int l;
                      unsigned int m;
                      unsigned int n;
                      
                      if (quadrature_points[q_point] (0) < 0.5)
                        {
                          n = 0;
                          
                          if (quadrature_points[q_point] (1) < 0.5)
                            {
                              m = 0;
                              
                              if (quadrature_points[q_point] (2) < 0.5)
                                l = 0;
                              
                              else
                                {
                                  quadrature_point (2) -= 1.0;
                                  l = 1;
                                }
                            }
                          
                          else
                            {
                              quadrature_point (1) -= 1.0;
                              m = 1;
                              
                              if (quadrature_points[q_point] (2) < 0.5)
                                l = 0;
                              
                              else
                                {
                                  quadrature_point (2) -= 1.0;
                                  l = 1;
                                }
                            }
                        }
                      
                      else
                        {
                          quadrature_point (0) -= 1.0;
                          n = 1;
                          
                          if (quadrature_points[q_point] (1) < 0.5)
                            {
                              m = 0;
                              
                              if (quadrature_points[q_point] (2) < 0.5)
                                l = 0;
                              
                              else
                                {
                                  quadrature_point (2) -= 1.0;
                                  l = 1;
                                }
                            }
                          
                          else
                            {
                              quadrature_point (1) -= 1.0;
                              m = 1;
                              
                              if (quadrature_points[q_point] (2) < 0.5)
                                l = 0;
                              
                              else
                                {
                                  quadrature_point (2) -= 1.0;
                                  l = 1;
                                }
                            }
                        }
                      
                      const double weight = 2.0 * quadrature.weight (q_point);
                      
                      for (unsigned int i = 0; i < this->degree; ++i)
                        {
                          const double L_i_0
                            = weight * legendre_polynomials[i].value (quadrature_points[q_point] (0));
                          const double L_i_1
                            = weight * legendre_polynomials[i].value (quadrature_points[q_point] (1));
                          const double L_i_2
                            = weight * legendre_polynomials[i].value (quadrature_points[q_point] (2));
                          
                          for (unsigned int j = 0; j < deg; ++j)
                            {
                              const double Le_j_0
                                = legendre_polynomials[j + 1].value (quadrature_points[q_point] (0));
                              const double Le_j_1
                                = legendre_polynomials[j + 1].value (quadrature_points[q_point] (1));
                              const double lo_j_0
                                = lobatto_polynomials[j + 2].value (quadrature_points[q_point] (0));
                              const double lo_j_1
                                = lobatto_polynomials[j + 2].value (quadrature_points[q_point] (1));
                              
                              for (unsigned int k = 0; k < deg; ++k)
                                {
                                  const double L_k
                                    = legendre_polynomials[k + 1].value (quadrature_points[q_point] (2));
                                  const double l_k_1
                                    = lobatto_polynomials[k + 2].value (quadrature_points[q_point] (1));
                                  const double l_k_2
                                    = lobatto_polynomials[k + 2].value (quadrature_points[q_point] (2));
                                  
                                  this->restriction[index][2 * (2 * l + m) + n]
                                  ((i * deg + j) * deg + k + n_boundary_dofs, dof)
                                    += L_i_0 * (this->shape_grad_component (dof, quadrature_point, 0)[2]
                                                * lo_j_1 * L_k
                                                + this->shape_grad_component (dof, quadrature_point, 0)[1]
                                                * Le_j_1 * l_k_2
                                                + this->shape_value_component (dof, quadrature_point, 0)
                                                * lo_j_1 * l_k_2);
                                  this->restriction[index][2 * (2 * l + m) + n]
                                  ((i + (j + deg) * this->degree) * deg + k + n_boundary_dofs, dof)
                                    += L_i_1 * (this->shape_grad_component (dof, quadrature_point, 1)[2]
                                                * lo_j_0 * L_k
                                                + this->shape_grad_component (dof, quadrature_point, 1)[0]
                                                * Le_j_0 * l_k_2
                                                + this->shape_value_component (dof, quadrature_point, 1)
                                                * lo_j_0 * l_k_2);
                                  this->restriction[index][2 * (2 * l + m) + n]
                                  (i + (k + (j + 2 * deg) * deg) * this->degree + n_boundary_dofs, dof)
                                    += L_i_2 * (this->shape_grad_component (dof, quadrature_point, 2)[1]
                                                * lo_j_0
                                                * legendre_polynomials[k + 1].value (quadrature_points[q_point] (1))
                                                + this->shape_grad_component (dof, quadrature_point, 2)[0]
                                                * Le_j_0 * l_k_1
                                                + this->shape_value_component (dof, quadrature_point, 2)
                                                * lo_j_0 * l_k_1);
                                }
                            }
                        }
                    }
                }
            }

          break;
        }

      default:
        Assert (false, ExcNotImplemented ());
    }
}


template <int dim>
std::vector<unsigned int>
FE_Nedelec<dim>::get_dpo_vector (const unsigned int degree, bool dg)
{
  std::vector<unsigned int> dpo (dim + 1);

  if (dg)
    {
      dpo[dim] = PolynomialsNedelec<dim>::compute_n_pols(degree);
    }
  else
    {
      dpo[0] = 0;
      dpo[1] = degree + 1;
      dpo[2] = 2 * degree * (degree + 1);

      if (dim == 3)
        dpo[3] = 3 * degree * degree * (degree + 1);
    }

  return dpo;
}


//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------

                                 // Chech wheter a given shape
                                 // function has support on a
                                 // given face.

                                 // We just switch through the
                                 // faces of the cell and return
                                 // true, if the shape function
                                 // has support on the face
                                 // and false otherwise.
template <int dim>
bool
FE_Nedelec<dim>::has_support_on_face (const unsigned int shape_index,
                                      const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
          ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
          ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));

  const unsigned int deg = this->degree-1;
  switch (dim)
    {
      case 2:
        switch (face_index)
          {
            case 0:
              if (!((shape_index > deg) && (shape_index < 2 * this->degree)))
                return true;

              else
                return false;

             case 1:
               if ((shape_index > deg) &&
                   (shape_index
                      < GeometryInfo<2>::lines_per_cell * this->degree))
                 return true;

               else
                 return false;

             case 2:
               if (shape_index < 3 * this->degree)
                 return true;

               else
                 return false;

             case 3:
               if (!((shape_index >= 2 * this->degree) &&
                     (shape_index < 3 * this->degree)))
                 return true;

               else
                 return false;

             default:
               {
                 Assert (false, ExcNotImplemented ());
                 return false;
               }
          }

      case 3:
        switch (face_index)
          {
            case 0:
              if (((shape_index > deg) && (shape_index < 2 * this->degree)) ||
                  ((shape_index >= 5 * this->degree) &&
                   (shape_index < 6 * this->degree)) ||
                  ((shape_index >= 9 * this->degree) &&
                   (shape_index < 10 * this->degree)) ||
                  ((shape_index >= 11 * this->degree) &&
                   (shape_index
                      < GeometryInfo<3>::lines_per_cell * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 2 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 4 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 5 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 6 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 7 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 9 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 10 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 11 * deg)
                        * this->degree)))
                return false;

              else
                return true;

            case 1:
              if (((shape_index > deg) && (shape_index < 4 * this->degree)) ||
                  ((shape_index >= 5 * this->degree) &&
                   (shape_index < 8 * this->degree)) ||
                  ((shape_index >= 9 * this->degree) &&
                   (shape_index < 10 * this->degree)) ||
                  ((shape_index >= 11 * this->degree) &&
                   (shape_index
                      < GeometryInfo<3>::lines_per_cell * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 2 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 5 * deg)
                        * this->degree)) ||
                   ((shape_index
                       >= (GeometryInfo<3>::lines_per_cell + 6 * deg)
                          * this->degree) &&
                    (shape_index
                       < (GeometryInfo<3>::lines_per_cell + 7 * deg)
                         * this->degree)) ||
                   ((shape_index
                       >= (GeometryInfo<3>::lines_per_cell + 9 * deg)
                          * this->degree) &&
                    (shape_index
                       < (GeometryInfo<3>::lines_per_cell + 10 * deg)
                         * this->degree)) ||
                   ((shape_index
                       >= (GeometryInfo<3>::lines_per_cell + 11 * deg)
                          * this->degree) &&
                    (shape_index
                       < (GeometryInfo<3>::lines_per_cell + 12 * deg)
                         * this->degree)))
                return true;

              else
                return false;

            case 2:
              if ((shape_index < 3 * this->degree) ||
                  ((shape_index >= 4 * this->degree) &&
                   (shape_index < 7 * this->degree)) ||
                  ((shape_index >= 8 * this->degree) &&
                   (shape_index < 10 * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 2 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 3 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 6 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 8 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 9 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 10 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 11 * deg)
                        * this->degree)))
                return true;

              else
                return false;

            case 3:
              if ((shape_index < 2 * this->degree) ||
                  ((shape_index >= 3 * this->degree) &&
                   (shape_index < 6 * this->degree)) ||
                  ((shape_index >= 7 * this->degree) &&
                   (shape_index < 8 * this->degree)) ||
                  ((shape_index >= 10 * this->degree) &&
                   (shape_index
                      < GeometryInfo<3>::lines_per_cell * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 2 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 3 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 4 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 6 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 9 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 10 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 11 * deg)
                        * this->degree)))
                return true;

              else
                return false;

            case 4:
              if ((shape_index < 4 * this->degree) ||
                  ((shape_index >= 8 * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 2 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 3 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 5 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 6 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 7 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 10 * deg)
                        * this->degree)))
                return true;

              else
                return false;

            case 5:
              if (((shape_index >= 4 * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 2 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 3 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 5 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 6 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 7 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 8 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 10 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 12 * deg)
                        * this->degree)))
                return true;

              else
                return false;

            default:
              {
                Assert (false, ExcNotImplemented ());
                return false;
              }
          }

      default:
        {
          Assert (false, ExcNotImplemented ());
          return false;
        }
    }
}

template <int dim>
FiniteElementDomination::Domination
FE_Nedelec<dim>::compare_for_face_domination (const FiniteElement<dim> &fe_other) const
{
  if (const FE_Nedelec<dim> *fe_nedelec_other
      = dynamic_cast<const FE_Nedelec<dim>*>(&fe_other))
    {
      if (this->degree < fe_nedelec_other->degree)
        return FiniteElementDomination::this_element_dominates;
      else if (this->degree == fe_nedelec_other->degree)
        return FiniteElementDomination::either_element_can_dominate;
      else
        return FiniteElementDomination::other_element_dominates;
    }
  else if (dynamic_cast<const FE_Nothing<dim>*>(&fe_other) != 0)
    {
                                       // the FE_Nothing has no
                                       // degrees of
                                       // freedom. nevertheless, we
                                       // say that the FE_Q element
                                       // dominates so that we don't
                                       // have to force the FE_Q side
                                       // to become a zero function
                                       // and rather allow the
                                       // function to be discontinuous
                                       // along the interface
      return FiniteElementDomination::other_element_dominates;
    }

  Assert (false, ExcNotImplemented());
  return FiniteElementDomination::neither_element_dominates;
}

template <int dim>
bool
FE_Nedelec<dim>::hp_constraints_are_implemented () const
{
  return true;
}

template <int dim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Nedelec<dim>::hp_vertex_dof_identities (const FiniteElement<dim>&)
const
{
                   // Nedelec elements do not have any dofs
                   // on vertices, hence return an empty vector.
  return std::vector<std::pair<unsigned int, unsigned int> > ();
}

template <int dim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Nedelec<dim>::hp_line_dof_identities (const FiniteElement<dim>& fe_other)
const
{
                                   // we can presently only compute these
                                   // identities if both FEs are
                                   // FE_Nedelec or if the other one is an
                                   // FE_Nothing
  if (const FE_Nedelec<dim> *fe_nedelec_other
        = dynamic_cast<const FE_Nedelec<dim>*> (&fe_other))
    {
                                       // dofs are located on lines, so
                                       // two dofs are identical, if their
                                       // edge shape functions have the
                                       // same polynomial degree.
      std::vector<std::pair<unsigned int, unsigned int> > identities;

      for (unsigned int i = 0;
           i < std::min (fe_nedelec_other->degree, this->degree); ++i)
        identities.push_back (std::make_pair (i, i));

      return identities;
    }

  else
    if (dynamic_cast<const FE_Nothing<dim>*> (&fe_other) != 0)
      {
                                       // the FE_Nothing has no
                                       // degrees of freedom, so there
                                       // are no equivalencies to be
                                       // recorded
        return std::vector<std::pair<unsigned int, unsigned int> > ();
      }

    else
      {
        Assert (false, ExcNotImplemented ());
        return std::vector<std::pair<unsigned int, unsigned int> > ();
      }
}

template <int dim>
std::vector<std::pair<unsigned int, unsigned int> >
FE_Nedelec<dim>::hp_quad_dof_identities (const FiniteElement<dim>& fe_other)
const
{
                                   // we can presently only compute
                                   // these identities if both FEs are
                                   // FE_Nedelec or if the other one is an
                                   // FE_Nothing
  if (const FE_Nedelec<dim> *fe_nedelec_other
        = dynamic_cast<const FE_Nedelec<dim>*> (&fe_other))
    {
                                       // dofs are located on the interior
                                       // of faces, so two dofs are identical,
                                       // if their face shape functions have
                                       // the same polynomial degree.
      const unsigned int p = fe_nedelec_other->degree;
      const unsigned int q = this->degree;
      const unsigned int p_min = std::min (p, q);
      std::vector<std::pair<unsigned int, unsigned int> > identities;

      for (unsigned int i = 0; i < p_min; ++i)
        for (unsigned int j = 0; j < p_min - 1; ++j)
          {
            identities.push_back (std::make_pair (i * (q - 1) + j,
                                                  i * (p - 1) + j));
            identities.push_back (std::make_pair (i + (j + q - 1) * q,
                                                  i + (j + p - 1) * p));
          }

      return identities;
    }

  else
    if (dynamic_cast<const FE_Nothing<dim>*> (&fe_other) != 0)
      {
                                       // the FE_Nothing has no
                                       // degrees of freedom, so there
                                       // are no equivalencies to be
                                       // recorded
        return std::vector<std::pair<unsigned int, unsigned int> > ();
      }

    else
      {
        Assert (false, ExcNotImplemented ());
        return std::vector<std::pair<unsigned int, unsigned int> > ();
      }
}

                   // In this function we compute the face
                   // interpolation matrix. This is usually
                   // done by projection-based interpolation,
                   // but, since one can compute the entries
                   // easy per hand, we save some computation
                   // time at this point and just fill in the
                   // correct values.
template <int dim>
void
FE_Nedelec<dim>::get_face_interpolation_matrix
  (const FiniteElement<dim>& source, FullMatrix<double>& interpolation_matrix)
const
{
                                   // this is only implemented, if the
                                   // source FE is also a
                                   // Nedelec element
  typedef FE_Nedelec<dim> FEN;
  typedef FiniteElement<dim> FEL;

  AssertThrow ((source.get_name ().find ("FE_Nedelec<") == 0) ||
               (dynamic_cast<const FEN*> (&source) != 0),
                typename FEL::ExcInterpolationNotImplemented());
  Assert (interpolation_matrix.m () == source.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m (),
                                source.dofs_per_face));
  Assert (interpolation_matrix.n () == this->dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.n (),
                                this->dofs_per_face));

                                   // ok, source is a Nedelec element, so
                                   // we will be able to do the work
  const FE_Nedelec<dim> &source_fe
    = dynamic_cast<const FE_Nedelec<dim>&> (source);

                                   // Make sure, that the element,
                   // for which the DoFs should be
                   // constrained is the one with
                   // the higher polynomial degree.
                   // Actually the procedure will work
                   // also if this assertion is not
                   // satisfied. But the matrices
                   // produced in that case might
                                   // lead to problems in the
                   // hp procedures, which use this
                                   // method.
  Assert (this->dofs_per_face <= source_fe.dofs_per_face,
          typename FEL::ExcInterpolationNotImplemented ());
  interpolation_matrix = 0;

                   // On lines we can just identify
                   // all degrees of freedom.
  for (unsigned int i = 0; i <this->degree; ++i)
    interpolation_matrix (i, i) = 1.0;

                   // In 3d we have some lines more
                   // and a face. The procedure stays
                   // the same as above, but we have
                   // to take a bit more care of the
                   // indices of the degrees of
                   // freedom.
  if (dim == 3)
    for (unsigned int i = 0; i <this->degree; ++i)
      {
        for (int j = 1; j < (int) GeometryInfo<dim>::lines_per_face; ++j)
          interpolation_matrix (j * source_fe.degree + i,
                                j * this->degree + i) = 1.0;

        for (unsigned int j = 0; j < this->degree-1; ++j)
          {
            interpolation_matrix
            ((i + GeometryInfo<dim>::lines_per_face) * (source_fe.degree - 1)
             + j + GeometryInfo<dim>::lines_per_face,
             (i + GeometryInfo<dim>::lines_per_face) * (this->degree-1) + j
             + GeometryInfo<dim>::lines_per_face)
              = 1.0;
            interpolation_matrix
            (i + (j + GeometryInfo<dim>::lines_per_face + source_fe.degree - 1)
               * source_fe.degree,
             i + (j + GeometryInfo<dim>::lines_per_face + this->degree-1) * this->degree)
              = 1.0;
          }
      }
}



template <>
void
FE_Nedelec<1>::get_subface_interpolation_matrix(
  const FiniteElement<1,1>&,
  const unsigned int,
  FullMatrix<double>&) const
{
  Assert (false, ExcNotImplemented ());
}



                   // In this function we compute the
                   // subface interpolation matrix.
template <int dim>
void
FE_Nedelec<dim>::get_subface_interpolation_matrix(
  const FiniteElement<dim>& source,
  const unsigned int subface,
  FullMatrix<double>& interpolation_matrix) const
{
                   // this is only implemented, if the
                   // source FE is also a
                   // Nedelec element
  typedef FE_Nedelec<dim> FEN;
  typedef FiniteElement<dim> FEL;

  AssertThrow ((source.get_name ().find ("FE_Nedelec<") == 0) ||
               (dynamic_cast<const FEN*> (&source) != 0),
                typename FEL::ExcInterpolationNotImplemented ());
  Assert (interpolation_matrix.m () == source.dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.m (),
                                source.dofs_per_face));
  Assert (interpolation_matrix.n () == this->dofs_per_face,
          ExcDimensionMismatch (interpolation_matrix.n (),
                                this->dofs_per_face));

                   // ok, source is a Nedelec element, so
                   // we will be able to do the work
  const FE_Nedelec<dim> &source_fe
    = dynamic_cast<const FE_Nedelec<dim>&> (source);

                   // Make sure, that the element,
                   // for which the DoFs should be
                   // constrained is the one with
                   // the higher polynomial degree.
                   // Actually the procedure will work
                   // also if this assertion is not
                   // satisfied. But the matrices
                   // produced in that case might
                   // lead to problems in the
                   // hp procedures, which use this
                   // method.
  Assert (this->dofs_per_face <= source_fe.dofs_per_face,
          typename FEL::ExcInterpolationNotImplemented ());
  interpolation_matrix = 0.0;
                   // Set the degrees of freedom
                   // by interpolation as usual.
  const QGauss<1> edge_quadrature (source_fe.degree);
  const std::vector<Point<1> >&
    edge_quadrature_points = edge_quadrature.get_points ();
  const std::vector<Polynomials::Polynomial<double> >&
    legendre_polynomials
      = Polynomials::Legendre::generate_complete_basis (source_fe.degree - 1);
  const unsigned int& n_edge_quadrature_points = edge_quadrature.size ();

  switch (dim)
    {
      case 2:
        {
          for (unsigned int q_point = 0; q_point < n_edge_quadrature_points;
               ++q_point)
            {
              const Point<dim> quadrature_point (0.0,
                                                 0.5 * (edge_quadrature_points[q_point] (0)
                                                        + subface));
              const double weight = 0.5 * edge_quadrature.weight (q_point);
                
              for (unsigned int dof = 0; dof < this->dofs_per_face; ++dof)
                {
                  const double shape_value = weight * this->shape_value_component (dof, quadrature_point, 1);
                  
                  for (unsigned int i = 0; i < source_fe.dofs_per_face; ++i)
                    interpolation_matrix (i, dof) += shape_value
                                                  * legendre_polynomials[i].value (edge_quadrature_points[q_point] (0));
                }
            }

          break;
        }

      case 3:
        {
          const double shifts[4][2] = { { 0.0, 0.0 }, { 1.0, 0.0 },
                                        { 0.0, 1.0 }, { 1.0, 1.0 } };

          for (unsigned int q_point = 0; q_point < n_edge_quadrature_points;
               ++q_point)
            {
              const double weight = 0.5 * edge_quadrature.weight (q_point);

              for (unsigned int i = 0; i < 2; ++i)
                {
                  Point<dim>
                    quadrature_point (0.5 * (i + shifts[subface][0]),
                                      0.5 * (edge_quadrature_points[q_point] (0)
                                             + shifts[subface][1]),
                                      0.0);
                    
                  for (unsigned int dof = 0; dof < this->dofs_per_face; ++dof)
                    {
                      const double shape_value
                        = weight * this->shape_value_component (this->face_to_cell_index (dof, 4),
                                                                quadrature_point, 1);
                      
                      for (unsigned int j = 0; j < source_fe.degree; ++j)
                        interpolation_matrix (i * source_fe.degree + j, dof)
                          += shape_value
                          * legendre_polynomials[j].value (edge_quadrature_points[q_point] (0));
                    }
                  
                  quadrature_point
                    = Point<dim> (0.5 * (edge_quadrature_points[q_point] (0)
                                         + shifts[subface][0]),
                                  0.5 * (i + shifts[subface][1]), 0.0);
                  
                  for (unsigned int dof = 0; dof < this->dofs_per_face; ++dof)
                    {
                      const double shape_value
                        = weight * this->shape_value_component (this->face_to_cell_index (dof, 4),
                                                                quadrature_point, 0);
                      
                      for (unsigned int j = 0; j < source_fe.degree; ++j)
                        interpolation_matrix ((i + 2) * source_fe.degree + j, dof)
                          += shape_value
                          * legendre_polynomials[j].value (edge_quadrature_points[q_point] (0));
                    }
                }
            }

          if (source_fe.degree > 1)
            {
              const QGauss<2> quadrature (source_fe.degree);
              const std::vector<Point<2> >&
                quadrature_points = quadrature.get_points ();
              const std::vector<Polynomials::Polynomial<double> >&
                lobatto_polynomials
                  = Polynomials::Lobatto::generate_complete_basis
                    (source_fe.degree);
              const unsigned int n_boundary_dofs
                = GeometryInfo<dim>::lines_per_face * source_fe.degree;
              const unsigned int& n_quadrature_points = quadrature.size ();
              FullMatrix<double>
                system_matrix_inv (source_fe.degree * (source_fe.degree - 1),
                                   source_fe.degree * (source_fe.degree - 1));

              {
                FullMatrix<double>
                  assembling_matrix (source_fe.degree * (source_fe.degree - 1),
                                     n_quadrature_points);

                for (unsigned int q_point = 0; q_point < n_quadrature_points;
                     ++q_point)
                  {
                    const double weight = std::sqrt (quadrature.weight (q_point));

                    for (unsigned int i = 0; i < source_fe.degree; ++i)
                      {
                        const double L_i = weight
                                           * legendre_polynomials[i].value
                                             (quadrature_points[q_point] (0));

                        for (unsigned int j = 0; j < source_fe.degree - 1; ++j)
                          assembling_matrix (i * (source_fe.degree - 1) + j,
                                             q_point)
                            = L_i * lobatto_polynomials[j + 2].value
                                    (quadrature_points[q_point] (1));
                      }
                  }

                FullMatrix<double> system_matrix (assembling_matrix.m (),
                                                  assembling_matrix.m ());

                assembling_matrix.mTmult (system_matrix, assembling_matrix);
                system_matrix_inv.invert (system_matrix);
              }

              FullMatrix<double> solution (system_matrix_inv.m (), 2);
              FullMatrix<double> system_rhs (system_matrix_inv.m (), 2);
              Vector<double> tmp (2);

              for (unsigned int dof = 0; dof < this->dofs_per_face; ++dof)
                {
                  system_rhs = 0.0;

                  for (unsigned int q_point = 0;
                       q_point < n_quadrature_points; ++q_point)
                    {
                      Point<dim>
                        quadrature_point
                        (0.5 * (quadrature_points[q_point] (0)
                                + shifts[subface][0]),
                         0.5 * (quadrature_points[q_point] (1)
                                + shifts[subface][1]), 0.0);
                      tmp (0) = 0.5 * this->shape_value_component
                                      (this->face_to_cell_index (dof, 4),
                                       quadrature_point, 0);
                      tmp (1) = 0.5 * this->shape_value_component
                                      (this->face_to_cell_index (dof, 4),
                                       quadrature_point, 1);
                      quadrature_point
                        = Point<dim> (quadrature_points[q_point] (0),
                                      quadrature_points[q_point] (1), 0.0);

                      for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int j = 0; j < source_fe.degree; ++j)
                          {
                            tmp (0) -= interpolation_matrix
                                       ((i + 2) * source_fe.degree + j, dof)
                                    * source_fe.shape_value_component
                                      ((i + 2) * source_fe.degree + j,
                                       quadrature_point, 0);
                            tmp (1) -= interpolation_matrix
                                       (i * source_fe.degree + j, dof)
                                    * source_fe.shape_value_component
                                      (i * source_fe.degree + j,
                                       quadrature_point, 1);
                          }

                      tmp *= quadrature.weight (q_point);

                      for (unsigned int i = 0; i < source_fe.degree; ++i)
                        {
                          const double L_i_0 = legendre_polynomials[i].value
                                               (quadrature_points[q_point] (0));
                          const double L_i_1 = legendre_polynomials[i].value
                                               (quadrature_points[q_point] (1));

                          for (unsigned int j = 0; j < source_fe.degree - 1; ++j)
                            {
                              system_rhs (i * (source_fe.degree - 1) + j, 0)
                                += tmp (0) * L_i_0
                                           * lobatto_polynomials[j + 2].value
                                             (quadrature_points[q_point] (1));
                              system_rhs (i * (source_fe.degree - 1) + j, 1)
                                += tmp (1) * L_i_1
                                           * lobatto_polynomials[j + 2].value
                                             (quadrature_points[q_point] (0));
                            }
                        }
                    }

                  system_matrix_inv.mmult (solution, system_rhs);

                  for (unsigned int i = 0; i < source_fe.degree; ++i)
                    for (unsigned int j = 0; j < source_fe.degree - 1; ++j)
                      {
                        if (std::abs (solution (i * (source_fe.degree - 1) + j, 0))
                              > 1e-14)
                          interpolation_matrix (i * (source_fe.degree - 1)
                                                  + j + n_boundary_dofs, dof)
                            = solution (i * (source_fe.degree - 1) + j, 0);

                        if (std::abs (solution (i * (source_fe.degree - 1) + j, 1))
                              > 1e-14)
                          interpolation_matrix (i + (j + source_fe.degree - 1)
                                                  * source_fe.degree
                                                  + n_boundary_dofs, dof)
                            = solution (i * (source_fe.degree - 1) + j, 1);
                      }
                }
            }

          break;
        }

      default:
        Assert (false, ExcNotImplemented ());
    }
}



                   // Since this is a vector valued element,
                   // we cannot interpolate a scalar function.
template <int dim>
void FE_Nedelec<dim>::interpolate (std::vector<double>&, const std::vector<double>&) const {
   Assert(false, ExcNotImplemented ());
}


                   // Interpolate a function, which is given by
                   // its values at the generalized support
                   // points in the finite element space on the
                   // reference cell.
                   // The interpolation on the edge degrees of
                   // freedom is done by direct calculation.
                   // For the interpolation on the remaining
                   // degrees of freedom we use a projection-
                   // based interpolation scheme. Therefore
                   // the remaining part of the interpolated
                   // function is projected on the face shape
                   // functions and the interior shape
                   // functions (if they exist).
template <int dim>
void
FE_Nedelec<dim>::interpolate (std::vector<double>& local_dofs,
                              const std::vector<Vector<double> >& values,
                              unsigned int offset) const
{
  Assert (values.size () == this->generalized_support_points.size (),
          ExcDimensionMismatch (values.size (),
                                this->generalized_support_points.size ()));
  Assert (local_dofs.size () == this->dofs_per_cell,
          ExcDimensionMismatch (local_dofs.size (),this->dofs_per_cell));
  Assert (values[0].size () >= offset + this->n_components (),
          ExcDimensionMismatch (values[0].size (),
                                offset + this->n_components ()));
  std::fill (local_dofs.begin (), local_dofs.end (), 0.);

  if (offset < dim)
    {
      const unsigned int deg = this->degree-1;
      const std::vector<Polynomials::Polynomial<double> >&
        legendre_polynomials
        = Polynomials::Legendre::generate_complete_basis (deg);
      const QGauss<1> edge_quadrature (this->degree);
      const std::vector<Point<1> >&
        edge_quadrature_points = edge_quadrature.get_points ();
      const unsigned int& n_edge_points = edge_quadrature.size ();
      
      switch (dim)
        {
          case 2:
            {
                                  // Let us begin with the
                                  // edge degrees of freedom.
              for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
                {
                  const double weight = edge_quadrature.weight (q_point);
                  
                  for (unsigned int i = 0; i < this->degree; ++i)
                    {
                      const double L_i
                        = weight * legendre_polynomials[i].value (edge_quadrature_points[q_point] (0));
                      
                      for (unsigned int j = 0; j < 2; ++j)
                        {
                          local_dofs[i + j * this->degree]
                            += L_i * values[q_point + j * n_edge_points] (1);
                          
                          if (offset == 0)
                            local_dofs[i + (j + 2) * this->degree]
                              += L_i * values[q_point + (j + 2) * n_edge_points] (0);
                        }
                    }
                }
            
                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
              for (unsigned int i = 0; i < GeometryInfo<dim>::lines_per_cell * this->degree; ++i)
                if (std::abs (local_dofs[i]) < 1e-14)
                  local_dofs[i] = 0.0;

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some interpolations onto
                                  // the interior shape
                                  // functions left.
              if (this->degree > 1)
                {
                  const std::vector<Polynomials::Polynomial<double> >&
                    lobatto_polynomials
                      = Polynomials::Lobatto::generate_complete_basis
                        (this->degree);
                  std::vector<Polynomials::Polynomial<double> >
                    lobatto_polynomials_grad (this->degree);

                  for (unsigned int i = 0; i < lobatto_polynomials_grad.size ();
                       ++i)
                    lobatto_polynomials_grad[i]
                      = lobatto_polynomials[i + 1].derivative ();

                                  // We set up the system
                                  // matrix and use it for
                                  // both, the horizontal
                                  // and the vertical
                                  // interior shape
                                  // functions.
                  const QGauss<dim> reference_quadrature (this->degree);
                  const unsigned int& n_interior_points
                    = reference_quadrature.size ();

                  FullMatrix<double> system_matrix (deg * this->degree,
                                                    deg * this->degree);

                  for (unsigned int i = 0; i < this->degree; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      for (unsigned int k = 0; k < this->degree; ++k)
                        for (unsigned int l = 0; l < deg; ++l)
                          for (unsigned int q_point = 0;
                               q_point < n_interior_points; ++q_point)
                            system_matrix (i * deg + j, k * deg + l)
                              += reference_quadrature.weight (q_point)
                                 * legendre_polynomials[i].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points]
                                    (0))
                                 * lobatto_polynomials[j + 2].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points]
                                    (1))
                                 * lobatto_polynomials_grad[k].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points]
                                    (0))
                                 * lobatto_polynomials[l + 2].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points]
                                    (1));

                  FullMatrix<double> system_matrix_inv (system_matrix.m (),
                                                        system_matrix.m ());
                
                  system_matrix_inv.invert (system_matrix);

                  Vector<double> solution (system_matrix_inv.m ());
                  Vector<double> system_rhs (system_matrix.m ());

                  if (offset == 0)
                    {
                                  // Set up the right hand side
                                  // for the horizontal shape
                                  // functions.
                      for (unsigned int q_point = 0;
                           q_point < n_interior_points; ++q_point)
                        {
                          double tmp
                            = values[q_point + GeometryInfo<dim>::lines_per_cell
                                     * n_edge_points] (0);

                          for (unsigned int i = 0; i < 2; ++i)
                            for (unsigned int j = 0; j < this->degree; ++j)
                              tmp -= local_dofs[(i + 2) * this->degree + j]
                                     * this->shape_value_component
                                       ((i + 2) * this->degree + j,
                                        this->generalized_support_points[q_point
                                                                         + GeometryInfo<dim>::lines_per_cell
                                                                         * n_edge_points],
                                        0);

                          for (unsigned int i = 0; i < this->degree; ++i)
                            for (unsigned int j = 0; j < this->degree-1; ++j)
                              system_rhs (i * (this->degree-1) + j)
                                += reference_quadrature.weight (q_point) * tmp
                                   * lobatto_polynomials_grad[i].value
                                     (this->generalized_support_points[q_point
                                                                       + GeometryInfo<dim>::lines_per_cell
                                                                       * n_edge_points]
                                      (0))
                                   * lobatto_polynomials[j + 2].value
                                     (this->generalized_support_points[q_point
                                                                       + GeometryInfo<dim>::lines_per_cell
                                                                       * n_edge_points]
                                      (1));
                        }

                      system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                      for (unsigned int i = 0; i < this->degree; ++i)
                        for (unsigned int j = 0; j < deg; ++j)
                          if (std::abs (solution (i * deg + j)) > 1e-14)
                            local_dofs[(i + GeometryInfo<dim>::lines_per_cell)
                                       * deg + j
                                       + GeometryInfo<dim>::lines_per_cell]
                              = solution (i * deg + j);
                    }

                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.
                  system_rhs = 0;

                  for (unsigned int q_point = 0; q_point < n_interior_points;
                       ++q_point)
                    {
                      double tmp
                        = values[q_point + GeometryInfo<dim>::lines_per_cell
                                 * n_edge_points] (1);

                      for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int j = 0; j < this->degree; ++j)
                          tmp -= local_dofs[i * this->degree + j]
                                 * this->shape_value_component
                                   (i * this->degree + j,
                                    this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points],
                                    1);

                      for (unsigned i = 0; i < this->degree; ++i)
                        for (unsigned int j = 0; j < deg; ++j)
                          system_rhs (i * deg + j)
                            += reference_quadrature.weight (q_point) * tmp
                               * lobatto_polynomials_grad[i].value
                                 (this->generalized_support_points[q_point
                                                                   + GeometryInfo<dim>::lines_per_cell
                                                                   * n_edge_points]
                                  (1))
                               * lobatto_polynomials[j + 2].value
                                 (this->generalized_support_points[q_point
                                                                   + GeometryInfo<dim>::lines_per_cell
                                                                   * n_edge_points]
                                  (0));
                    }

                  system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                  for (unsigned int i = 0; i < this->degree; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      if (std::abs (solution (i * deg + j)) > 1e-14)
                        local_dofs[i + (j + GeometryInfo<dim>::lines_per_cell
                                          + deg) * this->degree]
                          = solution (i * deg + j);
                }

              break;
            }

          case 3:
            {
                                  // Let us begin with the
                                  // edge degrees of freedom.
              for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
                {
                  const double weight = edge_quadrature.weight (q_point);
                  
                  for (unsigned int i = 0; i < this->degree; ++i)
                    {
                      const double L_i
                        = weight * legendre_polynomials[i].value (edge_quadrature_points[q_point] (0));
                    
                      for (unsigned int j = 0; j < 4; ++j)
                        local_dofs[i + (j + 8) * this->degree]
                          += L_i * values[q_point + (j + 8) * n_edge_points] (2);
                    
                      if (offset + 1 < dim)
                        {
                          for (unsigned int j = 0; j < 2; ++j)
                            for (unsigned int k = 0; k < 2; ++k)
                              local_dofs[i + (j + 4 * k) * this->degree]
                                += L_i * values[q_point + (j + 4 * k) * n_edge_points] (1);
                          
                          if (offset == 0)
                            for (unsigned int j = 0; j < 2; ++j)
                              for (unsigned int k = 0; k < 2; ++k)
                                local_dofs[i + (j + 4 * k + 2) * this->degree]
                                  += L_i * values[q_point + (j + 4 * k + 2) * n_edge_points] (0);
                        }
                    }
                }
            
                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
              for (unsigned int i = 0; i < GeometryInfo<dim>::lines_per_cell * this->degree; ++i)
                if (std::abs (local_dofs[i]) < 1e-14)
                  local_dofs[i] = 0.0;

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some interpolation to the
                                  // face and interior shape
                                  // functions left.
              if (this->degree > 1)
                {
                  const std::vector<Polynomials::Polynomial<double> >&
                    lobatto_polynomials
                      = Polynomials::Lobatto::generate_complete_basis
                        (this->degree);
                  FullMatrix<double> system_matrix (deg * this->degree,
                                                    deg * this->degree);
                  std::vector<Polynomials::Polynomial<double> >
                    lobatto_polynomials_grad (this->degree);

                  for (unsigned int i = 0; i < lobatto_polynomials_grad.size ();
                       ++i)
                    lobatto_polynomials_grad[i]
                      = lobatto_polynomials[i + 1].derivative ();

                                  // We set up the system
                                  // matrix and use it for
                                  // both, the horizontal
                                  // and the vertical, shape
                                  // functions.
                  const unsigned int
                    n_face_points = n_edge_points * n_edge_points;

                  for (unsigned int i = 0; i < this->degree; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      for (unsigned int k = 0; k < this->degree; ++k)
                        for (unsigned int l = 0; l < deg; ++l)
                          for (unsigned int q_point = 0; q_point < n_face_points;
                               ++q_point)
                            system_matrix (i * deg + j, k * deg + l)
                              += boundary_weights (q_point + n_edge_points,
                                                   2 * (k * deg + l))
                                 * legendre_polynomials[i].value
                                   (this->generalized_face_support_points[q_point
                                                                          + 4
                                                                          * n_edge_points]
                                    (0))
                                 * lobatto_polynomials[j + 2].value
                                   (this->generalized_face_support_points[q_point
                                                                          + 4
                                                                          * n_edge_points]
                                    (1));

                  FullMatrix<double> system_matrix_inv (system_matrix.m (),
                                                        system_matrix.m ());
                  
                  system_matrix_inv.invert (system_matrix);
                  
                  Vector<double> solution (system_matrix.m ());
                  Vector<double> system_rhs (system_matrix.m ());

                  for (unsigned int face = 0;
                       face < GeometryInfo<dim>::faces_per_cell; ++face)
                    {
                      switch (face)
                        {
                          case 0:
                            {
                              if (offset + 1 < dim)
                                {
                                      // Set up the right hand side
                                      // for the horizontal shape
                                      // functions.
                                  system_rhs = 0;

                                  for (unsigned int q_point = 0;
                                       q_point < n_face_points; ++q_point)
                                    {
                                      double tmp
                                        = values[q_point
                                                 + GeometryInfo<dim>::lines_per_cell
                                                 * n_edge_points] (1);

                                      for (unsigned int i = 0; i < 2; ++i)
                                        for (unsigned int j = 0; j < this->degree; ++j)
                                          tmp
                                            -= local_dofs[4 * i * this->degree
                                                            + j]
                                               * this->shape_value_component
                                                 (4 * i * this->degree + j,
                                                  this->generalized_support_points[q_point
                                                                                   + GeometryInfo<dim>::lines_per_cell
                                                                                   * n_edge_points],
                                                  1);

                                      for (unsigned int i = 0; i < this->degree; ++i)
                                        for (unsigned int j = 0; j < deg; ++j)
                                          system_rhs (i * deg + j)
                                            += boundary_weights
                                               (q_point + n_edge_points,
                                                2 * (i * deg + j)) * tmp;
                                    }

                                  system_matrix_inv.vmult (solution, system_rhs);

                                      // Add the computed values
                                      // to the resulting vector
                                      // only, if they are not
                                      // too small.
                                  for (unsigned int i = 0; i < this->degree; ++i)
                                    for (unsigned int j = 0; j < deg; ++j)
                                      if (std::abs (solution (i * deg + j))
                                            > 1e-14)
                                        local_dofs[(i
                                                    + GeometryInfo<dim>::lines_per_cell)
                                                   * deg + j
                                                   + GeometryInfo<dim>::lines_per_cell]
                                          = solution (i * deg + j);
                                }

                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.
                              system_rhs = 0;

                              for (unsigned int q_point = 0;
                                   q_point < n_face_points; ++q_point)
                                {
                                  double tmp
                                    = values[q_point
                                             + GeometryInfo<dim>::lines_per_cell
                                             * n_edge_points] (2);

                                  for (unsigned int i = 0; i < 2; ++i)
                                    for (unsigned int j = 0; j < this->degree; ++j)
                                      tmp -= local_dofs[2 * (i + 4)
                                                        * this->degree + j]
                                             * this->shape_value_component
                                               (2 * (i + 4) * this->degree + j,
                                                this->generalized_support_points[q_point
                                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                                 * n_edge_points],
                                                2);

                                  for (unsigned i = 0; i < this->degree; ++i)
                                    for (unsigned int j = 0; j < deg; ++j)
                                      system_rhs (i * deg + j)
                                        += boundary_weights
                                           (q_point + n_edge_points,
                                            2 * (i * deg + j) + 1)
                                           * tmp;
                                }

                              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                              for (unsigned int i = 0; i < this->degree; ++i)
                                for (unsigned int j = 0; j < deg; ++j)
                                  if (std::abs (solution (i * deg + j)) > 1e-14)
                                    local_dofs[i + (j + GeometryInfo<dim>::lines_per_cell
                                                      + deg)
                                                 * this->degree]
                                      = solution (i * deg + j);

                              break;
                            }

                          case 1:
                            {
                              if (offset + 1 < dim)
                                {
                                      // Set up the right hand side
                                      // for the horizontal shape
                                      // functions.
                                  system_rhs = 0;

                                  for (unsigned int q_point = 0;
                                       q_point < n_face_points; ++q_point)
                                    {
                                      double tmp
                                        = values[q_point
                                                 + GeometryInfo<dim>::lines_per_cell
                                                 * n_edge_points
                                                 + n_face_points] (1);

                                      for (unsigned int i = 0; i < 2; ++i)
                                        for (unsigned int j = 0; j < this->degree; ++j)
                                          tmp -= local_dofs[(4 * i + 1)
                                                            * this->degree + j]
                                                 * this->shape_value_component
                                                   ((4 * i + 1) * this->degree
                                                    + j,
                                                    this->generalized_support_points[q_point
                                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                                     * n_edge_points
                                                                                     + n_face_points],
                                                    1);

                                      for (unsigned int i = 0; i < this->degree; ++i)
                                        for (unsigned int j = 0; j < deg; ++j)
                                          system_rhs (i * deg + j)
                                            += boundary_weights
                                               (q_point + n_edge_points,
                                                2 * (i * deg + j)) * tmp;
                                    }

                                  system_matrix_inv.vmult (solution, system_rhs);

                                      // Add the computed values
                                      // to the resulting vector
                                      // only, if they are not
                                      // too small.
                                  for (unsigned int i = 0; i < this->degree; ++i)
                                    for (unsigned int j = 0; j < deg; ++j)
                                      if (std::abs (solution (i * deg + j))
                                            > 1e-14)
                                        local_dofs[(i + GeometryInfo<dim>::lines_per_cell
                                                      + 2 * this->degree) * deg + j
                                                      + GeometryInfo<dim>::lines_per_cell]
                                          = solution (i * deg + j);
                                }

                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.
                              system_rhs = 0;

                              for (unsigned int q_point = 0;
                                   q_point < n_face_points; ++q_point)
                                {
                                  double tmp
                                    = values[q_point
                                             + GeometryInfo<dim>::lines_per_cell
                                             * n_edge_points + n_face_points]
                                      (2);

                                  for (unsigned int i = 0; i < 2; ++i)
                                    for (unsigned int j = 0; j < this->degree; ++j)
                                      tmp -= local_dofs[(2 * (i + 4) + 1)
                                                        * this->degree + j]
                                             * this->shape_value_component
                                               ((2 * (i + 4) + 1) * this->degree
                                                + j,
                                                this->generalized_support_points[q_point
                                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                                 * n_edge_points
                                                                                 + n_face_points],
                                                2);

                                  for (unsigned i = 0; i < this->degree; ++i)
                                    for (unsigned int j = 0; j < deg; ++j)
                                      system_rhs (i * deg + j)
                                        += boundary_weights
                                           (q_point + n_edge_points,
                                            2 * (i * deg + j) + 1) * tmp;
                                }

                              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                              for (unsigned int i = 0; i < this->degree; ++i)
                                for (unsigned int j = 0; j < deg; ++j)
                                  if (std::abs (solution (i * deg + j)) > 1e-14)
                                    local_dofs[i + (j + GeometryInfo<dim>::lines_per_cell
                                                      + 3 * deg)
                                                 * this->degree]
                                      = solution (i * deg + j);

                              break;
                            }

                          case 2:
                            {
                              if (offset == 0)
                                {
                                      // Set up the right hand side
                                      // for the horizontal shape
                                      // functions.
                                  system_rhs = 0;

                                  for (unsigned int q_point = 0;
                                       q_point < n_face_points; ++q_point)
                                    {
                                      double tmp
                                        = values[q_point
                                                 + GeometryInfo<dim>::lines_per_cell
                                                 * n_edge_points + 2 * n_face_points]
                                          (2);

                                      for (unsigned int i = 0; i < 2; ++i)
                                        for (unsigned int j = 0; j < this->degree; ++j)
                                          tmp -= local_dofs[(i + 8) * this->degree
                                                                    + j]
                                              * this->shape_value_component
                                                ((i + 8) * this->degree + j,
                                                 this->generalized_support_points[q_point
                                                                                  + GeometryInfo<dim>::lines_per_cell
                                                                                  * n_edge_points
                                                                                  + 2
                                                                                  * n_face_points],
                                                 2);

                                      for (unsigned i = 0; i < this->degree; ++i)
                                        for (unsigned int j = 0; j < deg; ++j)
                                          system_rhs (i * deg + j)
                                            += boundary_weights
                                               (q_point + n_edge_points,
                                                2 * (i * deg + j)) * tmp;
                                    }

                                  system_matrix_inv.vmult (solution, system_rhs);

                                      // Add the computed values
                                      // to the resulting vector
                                      // only, if they are not
                                      // too small.
                                  for (unsigned int i = 0; i < this->degree; ++i)
                                    for (unsigned int j = 0; j < deg; ++j)
                                      if (std::abs (solution (i * deg + j))
                                            > 1e-14)
                                        local_dofs[(i + GeometryInfo<dim>::lines_per_cell
                                                      + 4 * this->degree) * deg
                                                      + j
                                                      + GeometryInfo<dim>::lines_per_cell]
                                          = solution (i * deg + j);
                                }

                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.
                              system_rhs = 0;

                              for (unsigned int q_point = 0;
                                   q_point < n_face_points; ++q_point)
                                {
                                  double tmp
                                    = values[q_point
                                             + GeometryInfo<dim>::lines_per_cell
                                             * n_edge_points
                                             + 2 * n_face_points] (0);

                                  for (unsigned int i = 0; i < 2; ++i)
                                    for (unsigned int j = 0; j < this->degree; ++j)
                                      tmp -= local_dofs[(4 * i + 2)
                                                        * this->degree + j]
                                          * this->shape_value_component
                                            ((4 * i + 2) * this->degree
                                                         + j,
                                             this->generalized_support_points[q_point
                                                                              + GeometryInfo<dim>::lines_per_cell
                                                                              * n_edge_points
                                                                              + 2
                                                                              * n_face_points],
                                             0);

                                  for (unsigned int i = 0; i < this->degree; ++i)
                                    for (unsigned int j = 0; j < deg; ++j)
                                      system_rhs (i * deg + j)
                                        += boundary_weights
                                           (q_point + n_edge_points,
                                            2 * (i * deg + j) + 1) * tmp;
                                }

                              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                              for (unsigned int i = 0; i < this->degree; ++i)
                                for (unsigned int j = 0; j < deg; ++j)
                                  if (std::abs (solution (i * deg + j)) > 1e-14)
                                    local_dofs[i + (j + GeometryInfo<dim>::lines_per_cell
                                                       + 5 * deg) * this->degree]
                                      = solution (i * deg + j);

                              break;
                            }

                          case 3:
                            {
                              if (offset == 0)
                                {
                                      // Set up the right hand side
                                      // for the horizontal shape
                                      // functions.
                                  system_rhs = 0;

                                  for (unsigned int q_point = 0;
                                       q_point < n_face_points; ++q_point)
                                    {
                                      double tmp
                                        = values[q_point
                                                 + GeometryInfo<dim>::lines_per_cell
                                                 * n_edge_points + 3 * n_face_points]
                                          (2);

                                      for (unsigned int i = 0; i < 2; ++i)
                                        for (unsigned int j = 0; j < this->degree; ++j)
                                          tmp -= local_dofs[(i + 10) * this->degree
                                                                     + j]
                                              * this->shape_value_component
                                                ((i + 10) * this->degree + j,
                                                 this->generalized_support_points[q_point
                                                                                  + GeometryInfo<dim>::lines_per_cell
                                                                                  * n_edge_points
                                                                                  + 3
                                                                                  * n_face_points],
                                                 2);

                                      for (unsigned i = 0; i < this->degree; ++i)
                                        for (unsigned int j = 0; j < deg; ++j)
                                          system_rhs (i * deg + j)
                                            += boundary_weights
                                               (q_point + n_edge_points,
                                                2 * (i * deg + j)) * tmp;
                                    }

                                  system_matrix_inv.vmult (solution, system_rhs);

                                      // Add the computed values
                                      // to the resulting vector
                                      // only, if they are not
                                      // too small.
                                  for (unsigned int i = 0; i < this->degree; ++i)
                                    for (unsigned int j = 0; j < deg; ++j)
                                      if (std::abs (solution (i * deg + j))
                                            > 1e-14)
                                        local_dofs[(i + GeometryInfo<dim>::lines_per_cell
                                                      + 6 * this->degree) * deg + j
                                                      + GeometryInfo<dim>::lines_per_cell]
                                          = solution (i * deg + j);
                                }

                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.
                              system_rhs = 0;

                              for (unsigned int q_point = 0;
                                   q_point < n_face_points; ++q_point)
                                {
                                  double tmp
                                    = values[q_point
                                             + GeometryInfo<dim>::lines_per_cell
                                             * n_edge_points + 3
                                             * n_face_points] (0);

                                  for (unsigned int i = 0; i < 2; ++i)
                                    for (unsigned int j = 0; j < this->degree; ++j)
                                      tmp -= local_dofs[(4 * i + 3)
                                                        * this->degree + j]
                                          * this->shape_value_component
                                            ((4 * i + 3) * this->degree
                                                         + j,
                                             this->generalized_support_points[q_point
                                                                              + GeometryInfo<dim>::lines_per_cell
                                                                              * n_edge_points
                                                                              + 3
                                                                              * n_face_points],
                                             0);

                                   for (unsigned int i = 0; i < this->degree; ++i)
                                     for (unsigned int j = 0; j < deg; ++j)
                                       system_rhs (i * deg + j)
                                         += boundary_weights
                                            (q_point + n_edge_points,
                                             2 * (i * deg + j) + 1) * tmp;
                                }

                              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                              for (unsigned int i = 0; i < this->degree; ++i)
                                for (unsigned int j = 0; j < deg; ++j)
                                  if (std::abs (solution (i * deg + j)) > 1e-14)
                                    local_dofs[i + (j + GeometryInfo<dim>::lines_per_cell
                                                      + 7 * deg) * this->degree]
                                      = solution (i * deg + j);

                              break;
                            }

                          case 4:
                            {
                              if (offset + 1 < dim)
                                {
                                      // Set up the right hand side
                                      // for the horizontal shape
                                      // functions.
                                  if (offset == 0)
                                    {
                                      system_rhs = 0;

                                      for (unsigned int q_point = 0;
                                           q_point < n_face_points; ++q_point)
                                        {
                                          double tmp
                                            = values[q_point
                                                     + GeometryInfo<dim>::lines_per_cell
                                                     * n_edge_points + 4
                                                     * n_face_points] (0);

                                          for (unsigned int i = 0; i < 2; ++i)
                                            for (unsigned int j = 0; j < this->degree; ++j)
                                              tmp -= local_dofs[(i + 2)
                                                                * this->degree
                                                                + j]
                                                     * this->shape_value_component
                                                       ((i + 2) * this->degree
                                                        + j,
                                                        this->generalized_support_points[q_point
                                                                                         + GeometryInfo<dim>::lines_per_cell
                                                                                         * n_edge_points
                                                                                         + 4
                                                                                         * n_face_points],
                                                        0);

                                          for (unsigned int i = 0; i < this->degree; ++i)
                                            for (unsigned int j = 0; j < deg; ++j)
                                              system_rhs (i * deg + j)
                                                += boundary_weights
                                                   (q_point + n_edge_points,
                                                    2 * (i * deg + j)) * tmp;
                                        }

                                      system_matrix_inv.vmult
                                        (solution, system_rhs);

                                          // Add the computed values
                                          // to the resulting vector
                                          // only, if they are not
                                          // too small.
                                      for (unsigned int i = 0; i < this->degree; ++i)
                                        for (unsigned int j = 0; j < deg; ++j)
                                          if (std::abs (solution (i * deg + j))
                                                > 1e-14)
                                            local_dofs[(i + GeometryInfo<dim>::lines_per_cell
                                                          + 8 * this->degree) * deg
                                                          + j
                                                          + GeometryInfo<dim>::lines_per_cell]
                                              = solution (i * deg + j);
                                    }

                                      // Set up the right hand side
                                      // for the vertical shape
                                      // functions.
                                  system_rhs = 0;

                                  for (unsigned int q_point = 0;
                                       q_point < n_face_points; ++q_point)
                                    {
                                      double tmp
                                        = values[q_point
                                                 + GeometryInfo<dim>::lines_per_cell
                                                 * n_edge_points + 4
                                                 * n_face_points] (1);

                                      for (unsigned int i = 0; i < 2; ++i)
                                        for (unsigned int j = 0; j < this->degree; ++j)
                                          tmp -= local_dofs[i * this->degree + j]
                                                 * this->shape_value_component
                                                   (i * this->degree + j,
                                                    this->generalized_support_points[q_point
                                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                                     * n_edge_points
                                                                                     + 4
                                                                                     * n_face_points],
                                                    1);

                                      for (unsigned i = 0; i < this->degree; ++i)
                                        for (unsigned int j = 0; j < deg; ++j)
                                          system_rhs (i * deg + j)
                                            += boundary_weights
                                               (q_point + n_edge_points,
                                                2 * (i * deg + j) + 1) * tmp;
                                    }

                                  system_matrix_inv.vmult (solution, system_rhs);

                                      // Add the computed values
                                      // to the resulting vector
                                      // only, if they are not
                                      // too small.
                                  for (unsigned int i = 0; i < this->degree; ++i)
                                    for (unsigned int j = 0; j < deg; ++j)
                                      if (std::abs (solution (i * deg + j))
                                            > 1e-14)
                                        local_dofs[i + (j + GeometryInfo<dim>::lines_per_cell
                                                          + 9 * deg)
                                                     * this->degree]
                                          = solution (i * deg + j);
                                }

                              break;
                            }

                          default:
                            if (offset + 1 < dim)
                              {
                                  // Set up the right hand side
                                  // for the horizontal shape
                                  // functions.
                                if (offset == 0)
                                  {
                                    system_rhs = 0;

                                    for (unsigned int q_point = 0;
                                         q_point < n_face_points; ++q_point)
                                      {
                                        double tmp
                                          = values[q_point
                                                   + GeometryInfo<dim>::lines_per_cell
                                                   * n_edge_points
                                                   + 5 * n_face_points] (0);

                                        for (unsigned int i = 0; i < 2; ++i)
                                          for (unsigned int j = 0; j < this->degree; ++j)
                                            tmp -= local_dofs[(i + 6)
                                                              * this->degree + j]
                                                   * this->shape_value_component
                                                     ((i + 6) * this->degree + j,
                                                      this->generalized_support_points[q_point
                                                                                       + GeometryInfo<dim>::lines_per_cell
                                                                                       * n_edge_points
                                                                                       + 5
                                                                                       * n_face_points],
                                                      0);

                                        for (unsigned int i = 0; i < this->degree; ++i)
                                          for (unsigned int j = 0; j < deg; ++j)
                                            system_rhs (i * deg + j)
                                              += boundary_weights
                                                 (q_point + n_edge_points,
                                                  2 * (i * deg + j)) * tmp;
                                      }

                                    system_matrix_inv.vmult
                                      (solution, system_rhs);

                                      // Add the computed values
                                      // to the resulting vector
                                      // only, if they are not
                                      // too small.
                                    for (unsigned int i = 0; i < this->degree; ++i)
                                      for (unsigned int j = 0; j < deg; ++j)
                                        if (std::abs (solution (i * deg + j))
                                              > 1e-14)
                                          local_dofs[(i + GeometryInfo<dim>::lines_per_cell
                                                        + 10 * this->degree)
                                                     * deg + j
                                                     + GeometryInfo<dim>::lines_per_cell]
                                            = solution (i * deg + j);
                                  }

                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.
                                system_rhs = 0;

                                for (unsigned int q_point = 0;
                                     q_point < n_face_points; ++q_point)
                                  {
                                    double tmp
                                      = values[q_point
                                               + GeometryInfo<dim>::lines_per_cell
                                               * n_edge_points + 5
                                               * n_face_points] (1);

                                    for (unsigned int i = 0; i < 2; ++i)
                                      for (unsigned int j = 0; j < this->degree; ++j)
                                        tmp -= local_dofs[(i + 4)
                                                          * this->degree + j]
                                               * this->shape_value_component
                                                 ((i + 4) * this->degree + j,
                                                  this->generalized_support_points[q_point
                                                                                   + GeometryInfo<dim>::lines_per_cell
                                                                                   * n_edge_points
                                                                                   + 5
                                                                                   * n_face_points],
                                                  1);

                                    for (unsigned i = 0; i < this->degree; ++i)
                                      for (unsigned int j = 0; j < deg; ++j)
                                        system_rhs (i * deg + j)
                                          += boundary_weights
                                             (q_point + n_edge_points,
                                              2 * (i * deg + j) + 1) * tmp;
                                  }

                                system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                                for (unsigned int i = 0; i < this->degree; ++i)
                                  for (unsigned int j = 0; j < deg; ++j)
                                    if (std::abs (solution (i * deg + j))
                                          > 1e-14)
                                      local_dofs[i + (j + GeometryInfo<dim>::lines_per_cell
                                                        + 11 * deg) * this->degree]
                                        = solution (i * deg + j);
                              }
                        }
                    }

                                  // Finally we project
                                  // the remaining parts
                                  // of the function on
                                  // the interior shape
                                  // functions.
                  const QGauss<dim> reference_quadrature (this->degree);
                  const unsigned int&
                    n_interior_points = reference_quadrature.size ();

                                  // We create the
                                  // system matrix.
                  system_matrix.reinit (this->degree * deg * deg,
                                        this->degree * deg * deg);
                  system_matrix = 0;

                  for (unsigned int i = 0; i < this->degree; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      for (unsigned int k = 0; k < deg; ++k)
                        for (unsigned int l = 0; l < this->degree; ++l)
                          for (unsigned int m = 0; m < deg; ++m)
                            for (unsigned int n = 0; n < deg; ++n)
                              for (unsigned int q_point = 0;
                                   q_point < n_interior_points; ++q_point)
                                system_matrix ((i * deg + j) * deg + k,
                                               (l * deg + m) * deg + n)
                                  += reference_quadrature.weight (q_point)
                                     * legendre_polynomials[i].value
                                       (this->generalized_support_points[q_point
                                                                         + GeometryInfo<dim>::lines_per_cell
                                                                         * n_edge_points
                                                                         + GeometryInfo<dim>::faces_per_cell
                                                                         * n_face_points]
                                        (0)) * lobatto_polynomials[j + 2].value
                                               (this->generalized_support_points[q_point
                                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                                 * n_edge_points
                                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                                 * n_face_points]
                                                (1))
                                             * lobatto_polynomials[k + 2].value
                                               (this->generalized_support_points[q_point
                                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                                 * n_edge_points
                                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                                 * n_face_points]
                                                (2))
                                             * lobatto_polynomials_grad[l].value
                                               (this->generalized_support_points[q_point
                                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                                 * n_edge_points
                                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                                 * n_face_points]
                                                (0))
                                             * lobatto_polynomials[m + 2].value
                                               (this->generalized_support_points[q_point
                                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                                 * n_edge_points
                                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                                 * n_face_points]
                                                (1))
                                             * lobatto_polynomials[n + 2].value
                                               (this->generalized_support_points[q_point
                                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                                 * n_edge_points
                                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                                 * n_face_points]
                                                (2));

                  system_matrix_inv.reinit (system_matrix.m (),
                                            system_matrix.m ());
                  system_matrix_inv.invert (system_matrix);
                  system_rhs.reinit (system_matrix_inv.m ());
                  solution.reinit (system_matrix.m ());

                  if (offset + 1 < dim)
                    {
                      if (offset == 0)
                        {
                                  // Set up the right hand side.
                          system_rhs = 0;

                          for (unsigned int q_point = 0;
                               q_point < n_interior_points; ++q_point)
                            {
                              double tmp
                                = values[q_point
                                         + GeometryInfo<dim>::lines_per_cell
                                         * n_edge_points
                                         + GeometryInfo<dim>::faces_per_cell
                                         * n_face_points] (0);

                              for (unsigned int i = 0; i < this->degree; ++i)
                                {
                                  for (unsigned int j = 0; j < 2; ++j)
                                    for (unsigned int k = 0; k < 2; ++k)
                                      tmp -= local_dofs[i + (j + 4 * k + 2)
                                                          * this->degree]
                                             * this->shape_value_component
                                               (i + (j + 4 * k + 2)
                                                  * this->degree,
                                                this->generalized_support_points[q_point
                                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                                 * n_edge_points
                                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                                 * n_face_points],
                                                0);

                                  for (unsigned int j = 0; j < deg; ++j)
                                    for (unsigned int k = 0; k < 4; ++k)
                                      tmp -= local_dofs[(i + 2 * (k + 2)
                                                           * this->degree
                                                           + GeometryInfo<dim>::lines_per_cell)
                                                        * deg + j
                                                        + GeometryInfo<dim>::lines_per_cell]
                                             * this->shape_value_component
                                               ((i + 2 * (k + 2) * this->degree
                                                   + GeometryInfo<dim>::lines_per_cell)
                                                * deg + j
                                                + GeometryInfo<dim>::lines_per_cell,
                                                this->generalized_support_points[q_point
                                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                                 * n_edge_points
                                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                                 * n_face_points],
                                                0);
                                }

                              for (unsigned int i = 0; i < this->degree; ++i)
                                for (unsigned int j = 0; j < deg; ++j)
                                  for (unsigned int k = 0; k < deg; ++k)
                                    system_rhs ((i * deg + j) * deg + k)
                                      += reference_quadrature.weight (q_point)
                                         * tmp
                                         * lobatto_polynomials_grad[i].value
                                           (this->generalized_support_points[q_point
                                                                             + GeometryInfo<dim>::lines_per_cell
                                                                             * n_edge_points
                                                                             + GeometryInfo<dim>::faces_per_cell
                                                                             * n_face_points]
                                            (0))
                                           * lobatto_polynomials[j + 2].value
                                             (this->generalized_support_points[q_point
                                                                               + GeometryInfo<dim>::lines_per_cell
                                                                               * n_edge_points
                                                                               + GeometryInfo<dim>::faces_per_cell
                                                                               * n_face_points]
                                            (1))
                                           * lobatto_polynomials[k + 2].value
                                             (this->generalized_support_points[q_point
                                                                               + GeometryInfo<dim>::lines_per_cell
                                                                               * n_edge_points
                                                                               + GeometryInfo<dim>::faces_per_cell
                                                                               * n_face_points]
                                            (2));
                            }

                          system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                          for (unsigned int i = 0; i < this->degree; ++i)
                            for (unsigned int j = 0; j < deg; ++j)
                              for (unsigned int k = 0; k < deg; ++k)
                                if (std::abs (solution ((i * deg + j) * deg + k))
                                      > 1e-14)
                                  local_dofs[((i + 2
                                                 * GeometryInfo<dim>::faces_per_cell)
                                              * deg + j
                                              + GeometryInfo<dim>::lines_per_cell
                                              + 2
                                              * GeometryInfo<dim>::faces_per_cell)
                                             * deg + k
                                             + GeometryInfo<dim>::lines_per_cell]
                                    = solution ((i * deg + j) * deg + k);
                        }

                                  // Set up the right hand side.
                      system_rhs = 0;

                      for (unsigned int q_point = 0; q_point < n_interior_points;
                           ++q_point)
                        {
                          double tmp
                            = values[q_point + GeometryInfo<dim>::lines_per_cell
                                             * n_edge_points
                                             + GeometryInfo<dim>::faces_per_cell
                                             * n_face_points] (1);

                          for (unsigned int i = 0; i < this->degree; ++i)
                            for (unsigned int j = 0; j < 2; ++j)
                              {
                                for (unsigned int k = 0; k < 2; ++k)
                                  tmp -= local_dofs[i + (4 * j + k)
                                                      * this->degree]
                                         * this->shape_value_component
                                           (i + (4 * j + k) * this->degree,
                                            this->generalized_support_points[q_point
                                                                             + GeometryInfo<dim>::lines_per_cell
                                                                             * n_edge_points
                                                                             + GeometryInfo<dim>::faces_per_cell
                                                                             * n_face_points],
                                            1);

                                for (unsigned int k = 0; k < deg; ++k)
                                  tmp -= local_dofs[(i + 2 * j * this->degree
                                                       + GeometryInfo<dim>::lines_per_cell)
                                                    * deg + k
                                                    + GeometryInfo<dim>::lines_per_cell]
                                         * this->shape_value_component
                                           ((i + 2 * j * this->degree
                                               + GeometryInfo<dim>::lines_per_cell)
                                            * deg + k
                                            + GeometryInfo<dim>::lines_per_cell,
                                            this->generalized_support_points[q_point
                                                                             + GeometryInfo<dim>::lines_per_cell
                                                                             * n_edge_points
                                                                             + GeometryInfo<dim>::faces_per_cell
                                                                             * n_face_points],
                                            1)
                                         + local_dofs[i + ((2 * j + 9) * deg + k
                                                           + GeometryInfo<dim>::lines_per_cell)
                                                        * this->degree]
                                         * this->shape_value_component
                                           (i + ((2 * j + 9) * deg + k
                                                 + GeometryInfo<dim>::lines_per_cell)
                                              * this->degree,
                                            this->generalized_support_points[q_point
                                                                             + GeometryInfo<dim>::lines_per_cell
                                                                             * n_edge_points
                                                                             + GeometryInfo<dim>::faces_per_cell
                                                                             * n_face_points],
                                            1);
                              }

                          for (unsigned int i = 0; i < this->degree; ++i)
                            for (unsigned int j = 0; j < deg; ++j)
                              for (unsigned int k = 0; k < deg; ++k)
                                system_rhs ((i * deg + j) * deg + k)
                                  += reference_quadrature.weight (q_point) * tmp
                                     * lobatto_polynomials_grad[i].value
                                       (this->generalized_support_points[q_point
                                                                         + GeometryInfo<dim>::lines_per_cell
                                                                         * n_edge_points
                                                                         + GeometryInfo<dim>::faces_per_cell
                                                                         * n_face_points]
                                        (1))
                                     * lobatto_polynomials[j + 2].value
                                       (this->generalized_support_points[q_point
                                                                         + GeometryInfo<dim>::lines_per_cell
                                                                         * n_edge_points
                                                                         + GeometryInfo<dim>::faces_per_cell
                                                                         * n_face_points]
                                        (0))
                                     * lobatto_polynomials[k + 2].value
                                       (this->generalized_support_points[q_point
                                                                         + GeometryInfo<dim>::lines_per_cell
                                                                         * n_edge_points
                                                                         + GeometryInfo<dim>::faces_per_cell
                                                                         * n_face_points]
                                        (2));
                        }

                      system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                      for (unsigned int i = 0; i < this->degree; ++i)
                        for (unsigned int j = 0; j < deg; ++j)
                          for (unsigned int k = 0; k < deg; ++k)
                            if (std::abs (solution ((i * deg + j) * deg + k))
                                  > 1e-14)
                              local_dofs[((i + this->degree + 2
                                             * GeometryInfo<dim>::faces_per_cell)
                                          * deg + j
                                          + GeometryInfo<dim>::lines_per_cell + 2
                                          * GeometryInfo<dim>::faces_per_cell)
                                         * deg + k
                                         + GeometryInfo<dim>::lines_per_cell]
                                = solution ((i * deg + j) * deg + k);
                    }

                                  // Set up the right hand side.
                  system_rhs = 0;

                  for (unsigned int q_point = 0; q_point < n_interior_points;
                       ++q_point)
                    {
                      double tmp
                        = values[q_point + GeometryInfo<dim>::lines_per_cell
                                         * n_edge_points
                                         + GeometryInfo<dim>::faces_per_cell
                                         * n_face_points] (2);

                      for (unsigned int i = 0; i < this->degree; ++i)
                        for (unsigned int j = 0; j < 4; ++j)
                          {
                            tmp -= local_dofs[i + (j + 8) * this->degree]
                                   * this->shape_value_component
                                     (i + (j + 8) * this->degree,
                                      this->generalized_support_points[q_point
                                                                       + GeometryInfo<dim>::lines_per_cell
                                                                       * n_edge_points
                                                                       + GeometryInfo<dim>::faces_per_cell
                                                                       * n_face_points],
                                      2);

                            for (unsigned int k = 0; k < deg; ++k)
                              tmp -= local_dofs[i + ((2 * j + 1) * deg + k
                                                     + GeometryInfo<dim>::lines_per_cell)
                                                  * this->degree]
                                     * this->shape_value_component
                                       (i + ((2 * j + 1) * deg + k
                                             + GeometryInfo<dim>::lines_per_cell)
                                          * this->degree,
                                        this->generalized_support_points[q_point
                                                                         + GeometryInfo<dim>::lines_per_cell
                                                                         * n_edge_points
                                                                         + GeometryInfo<dim>::faces_per_cell
                                                                         * n_face_points],
                                        2);
                          }

                      for (unsigned int i = 0; i < this->degree; ++i)
                        for (unsigned int j = 0; j < deg; ++j)
                          for (unsigned int k = 0; k < deg; ++k)
                            system_rhs ((i * deg + j) * deg + k)
                              += reference_quadrature.weight (q_point) * tmp
                                 * lobatto_polynomials_grad[i].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points]
                                    (2))
                                 * lobatto_polynomials[j + 2].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points]
                                    (0))
                                 * lobatto_polynomials[k + 2].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points]
                                    (1));
                    }

                  system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                  for (unsigned int i = 0; i < this->degree; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      for (unsigned int k = 0; k < deg; ++k)
                        if (std::abs (solution ((i * deg + j) * deg + k))
                              > 1e-14)
                          local_dofs[i + ((j + 2
                                             * (deg + GeometryInfo<dim>::faces_per_cell))
                                          * deg + k
                                          + GeometryInfo<dim>::lines_per_cell)
                                       * this->degree]
                            = solution ((i * deg + j) * deg + k);
                }

              break;
            }

          default:
            Assert (false, ExcNotImplemented ());
        }
    }
}


                   // Interpolate a function, which is given by
                   // its values at the generalized support
                   // points in the finite element space on the
                   // reference cell.
                   // This is done as usual by projection-based
                   // interpolation.
template <int dim>
void
FE_Nedelec<dim>::interpolate (std::vector<double>& local_dofs,
                              const VectorSlice<const std::vector<std::vector<double> > >& values)
const
{
  Assert (values.size () == this->n_components (),
          ExcDimensionMismatch (values.size (), this->n_components ()));
  Assert (values[0].size () == this->generalized_support_points.size (),
          ExcDimensionMismatch (values[0].size (),
                                this->generalized_support_points.size ()));
  Assert (local_dofs.size () == this->dofs_per_cell,
          ExcDimensionMismatch (local_dofs.size (), this->dofs_per_cell));
  std::fill (local_dofs.begin (), local_dofs.end (), 0.0);

  const unsigned int deg = this->degree-1;
  const std::vector<Polynomials::Polynomial<double> >&
    legendre_polynomials
      = Polynomials::Legendre::generate_complete_basis (deg);
  const QGauss<1> edge_quadrature (this->degree);
  const std::vector<Point<1> >&
    edge_quadrature_points = edge_quadrature.get_points ();
  const unsigned int& n_edge_points = edge_quadrature.size ();
      
  switch (dim)
    {
      case 2:
        {
                                  // Let us begin with the
                                  // edge degrees of freedom.
          for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
            {
              const double weight = edge_quadrature.weight (q_point);
              
              for (unsigned int i = 0; i < this->degree; ++i)
                {
                  const double L_i
                    = weight * legendre_polynomials[i].value (edge_quadrature_points[q_point] (0));
                  
                  for (unsigned int j = 0; j < 2; ++j)
                    for (unsigned int k = 0; k < 2; ++k)
                      local_dofs[i + (j + 2 * k) * this->degree]
                        += L_i * values[1 - k][q_point + (j + 2 * k) * n_edge_points];
                }
            }
            
                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
          for (unsigned int i = 0; i < GeometryInfo<dim>::lines_per_cell * this->degree; ++i)
            if (std::abs (local_dofs[i]) < 1e-14)
              local_dofs[i] = 0.0;

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some interpolations onto
                                  // the interior shape
                                  // functions left.
          if (this->degree > 1)
            {
              const std::vector<Polynomials::Polynomial<double> >&
                lobatto_polynomials
                  = Polynomials::Lobatto::generate_complete_basis
                    (this->degree);
              std::vector<Polynomials::Polynomial<double> >
                lobatto_polynomials_grad (this->degree);

              for (unsigned int i = 0; i < lobatto_polynomials_grad.size ();
                   ++i)
                lobatto_polynomials_grad[i]
                  = lobatto_polynomials[i + 1].derivative ();

                                  // We set up the system
                                  // matrix and use it for
                                  // both, the horizontal
                                  // and the vertical
                                  // interior shape
                                  // functions.
              const QGauss<dim> reference_quadrature (this->degree);
              const unsigned int& n_interior_points
                = reference_quadrature.size ();

              FullMatrix<double> system_matrix (deg * this->degree,
                                                deg * this->degree);

              for (unsigned int i = 0; i < this->degree; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k < this->degree; ++k)
                    for (unsigned int l = 0; l < deg; ++l)
                      for (unsigned int q_point = 0;
                           q_point < n_interior_points; ++q_point)
                        system_matrix (i * deg + j, k * deg + l)
                          += reference_quadrature.weight (q_point)
                             * legendre_polynomials[i].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points]
                                (0))
                             * lobatto_polynomials[j + 2].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points]
                                (1))
                             * lobatto_polynomials_grad[k].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points]
                                (0))
                             * lobatto_polynomials[l + 2].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points]
                                (1));

              FullMatrix<double> system_matrix_inv (system_matrix.m (),
                                                    system_matrix.m ());
                
              system_matrix_inv.invert (system_matrix);

              Vector<double> solution (system_matrix_inv.m ());
              Vector<double> system_rhs (system_matrix.m ());

              for (unsigned int q_point = 0; q_point < n_interior_points;
                   ++q_point)
                {
                  double tmp
                    = values[0][q_point + GeometryInfo<dim>::lines_per_cell
                                * n_edge_points];

                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int j = 0; j <= deg; ++j)
                      tmp -= local_dofs[(i + 2) * this->degree + j]
                             * this->shape_value_component
                               ((i + 2) * this->degree + j,
                                this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points],
                                0);

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      system_rhs (i * deg + j)
                        += reference_quadrature.weight (q_point) * tmp
                           * lobatto_polynomials_grad[i].value
                             (this->generalized_support_points[q_point
                                                               + GeometryInfo<dim>::lines_per_cell
                                                               * n_edge_points]
                              (0))
                           * lobatto_polynomials[j + 2].value
                             (this->generalized_support_points[q_point
                                                               + GeometryInfo<dim>::lines_per_cell
                                                               * n_edge_points]
                              (1));
                }

              solution.reinit (system_matrix.m ());
              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  if (std::abs (solution (i * deg + j)) > 1e-14)
                     local_dofs[(i + GeometryInfo<dim>::lines_per_cell) * deg
                                + j + GeometryInfo<dim>::lines_per_cell]
                       = solution (i * deg + j);

              system_rhs = 0;
                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.

              for (unsigned int q_point = 0; q_point < n_interior_points;
                   ++q_point)
                {
                  double tmp
                    = values[1][q_point + GeometryInfo<dim>::lines_per_cell
                                * n_edge_points];

                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int j = 0; j <= deg; ++j)
                      tmp -= local_dofs[i * this->degree + j]
                             * this->shape_value_component
                               (i * this->degree + j,
                                this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points],
                                1);

                  for (unsigned i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      system_rhs (i * deg + j)
                        += reference_quadrature.weight (q_point) * tmp
                           * lobatto_polynomials_grad[i].value
                             (this->generalized_support_points[q_point
                                                               + GeometryInfo<dim>::lines_per_cell
                                                               * n_edge_points]
                              (1))
                           * lobatto_polynomials[j + 2].value
                             (this->generalized_support_points[q_point
                                                               + GeometryInfo<dim>::lines_per_cell
                                                               * n_edge_points]
                              (0));
                }

              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  if (std::abs (solution (i * deg + j)) > 1e-14)
                     local_dofs[i + (j + GeometryInfo<dim>::lines_per_cell
                                     + deg) * this->degree]
                       = solution (i * deg + j);
            }

          break;
        }

      case 3:
        {
                                  // Let us begin with the
                                  // edge degrees of freedom.
          for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
            {
              const double weight = edge_quadrature.weight (q_point);
              
              for (unsigned int i = 0; i < this->degree; ++i)
                {
                  const double L_i
                    = weight * legendre_polynomials[i].value (edge_quadrature_points[q_point] (0));
                    
                  for (unsigned int j = 0; j < 4; ++j)
                    local_dofs[i + (j + 8) * this->degree]
                      += L_i * values[2][q_point + (j + 8) * n_edge_points];
                  
                  for (unsigned int j = 0; j < 2; ++j)
                    for (unsigned int k = 0; k < 2; ++k)
                      for (unsigned int l = 0; l < 2; ++l)
                        local_dofs[i + (j + 4 * k + 2 * l) * this->degree]
                          += L_i * values[1 - l][q_point + (j + 4 * k + 2 * l) * n_edge_points];
                }
            }
            
                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
          for (unsigned int i = 0; i < GeometryInfo<dim>::lines_per_cell * this->degree; ++i)
            if (std::abs (local_dofs[i]) < 1e-14)
              local_dofs[i] = 0.0;

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some interpolation to the
                                  // face and interior shape
                                  // functions left.
          if (this->degree > 1)
            {
              const std::vector<Polynomials::Polynomial<double> >&
                lobatto_polynomials
                  = Polynomials::Lobatto::generate_complete_basis
                    (this->degree);
              FullMatrix<double> system_matrix (deg * this->degree,
                                                deg * this->degree);
              std::vector<Polynomials::Polynomial<double> >
                lobatto_polynomials_grad (this->degree);

              for (unsigned int i = 0; i < lobatto_polynomials_grad.size ();
                   ++i)
                lobatto_polynomials_grad[i]
                  = lobatto_polynomials[i + 1].derivative ();

                                  // We set up the system
                                  // matrix and use it for
                                  // both, the horizontal
                                  // and the vertical, shape
                                  // functions.
              const unsigned int
                n_face_points = n_edge_points * n_edge_points;

              for (unsigned int i = 0; i < this->degree; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k < this->degree; ++k)
                    for (unsigned int l = 0; l < deg; ++l)
                      for (unsigned int q_point = 0; q_point < n_face_points;
                           ++q_point)
                        system_matrix (i * deg + j, k * deg + l)
                          += boundary_weights (q_point + n_edge_points,
                                               2 * (k * deg + l))
                             * legendre_polynomials[i].value
                               (this->generalized_face_support_points[q_point
                                                                      + 4
                                                                      * n_edge_points]
                                (0))
                             * lobatto_polynomials[j + 2].value
                               (this->generalized_face_support_points[q_point
                                                                      + 4
                                                                      * n_edge_points]
                                (1));

              FullMatrix<double> system_matrix_inv (system_matrix.m (),
                                                    system_matrix.m ());
              
              system_matrix_inv.invert (system_matrix);
              
              const unsigned int
                face_coordinates[GeometryInfo<3>::faces_per_cell][2]
                  = {{1, 2}, {1, 2}, {2, 0}, {2, 0}, {0, 1}, {0, 1}};
              const unsigned int
                edge_indices[GeometryInfo<3>::faces_per_cell][GeometryInfo<3>::lines_per_face]
                  = {{0, 4, 8, 10}, {1, 5, 9, 11}, {8, 9, 2, 6},
                     {10, 11, 3, 7}, {2, 3, 0, 1}, {6, 7, 4, 5}};
              Vector<double> solution (system_matrix.m ());
              Vector<double> system_rhs (system_matrix.m ());

              for (unsigned int face = 0;
                   face < GeometryInfo<dim>::faces_per_cell; ++face)
                {
                                  // Set up the right hand side
                                  // for the horizontal shape
                                  // functions.
                  system_rhs = 0;

                  for (unsigned int q_point = 0; q_point < n_face_points;
                       ++q_point)
                    {
                      double tmp
                        = values[face_coordinates[face][0]][q_point
                                                            + GeometryInfo<dim>::lines_per_cell
                                                            * n_edge_points];

                      for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int j = 0; j <= deg; ++j)
                          tmp -= local_dofs[edge_indices[face][i]
                                            * this->degree + j]
                                 * this->shape_value_component
                                   (edge_indices[face][i] * this->degree + j,
                                    this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points],
                                    face_coordinates[face][0]);

                      for (unsigned int i = 0; i <= deg; ++i)
                        for (unsigned int j = 0; j < deg; ++j)
                          system_rhs (i * deg + j)
                            += boundary_weights (q_point + n_edge_points,
                                                 2 * (i * deg + j)) * tmp;
                    }

                  system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      if (std::abs (solution (i * deg + j)) > 1e-14)
                        local_dofs[(2 * face * this->degree + i
                                    + GeometryInfo<dim>::lines_per_cell) * deg
                                   + j + GeometryInfo<dim>::lines_per_cell]
                          = solution (i * deg + j);

                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.
                  system_rhs = 0;

                  for (unsigned int q_point = 0; q_point < n_face_points;
                       ++q_point)
                    {
                      double tmp
                        = values[face_coordinates[face][1]][q_point
                                                            + GeometryInfo<dim>::lines_per_cell
                                                            * n_edge_points];

                      for (int i = 2; i < (int) GeometryInfo<dim>::lines_per_face; ++i)
                        for (unsigned int j = 0; j <= deg; ++j)
                          tmp -= local_dofs[edge_indices[face][i]
                                            * this->degree + j]
                                 * this->shape_value_component
                                   (edge_indices[face][i] * this->degree + j,
                                    this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points],
                                    face_coordinates[face][1]);

                      for (unsigned i = 0; i <= deg; ++i)
                        for (unsigned int j = 0; j < deg; ++j)
                          system_rhs (i * deg + j)
                            += boundary_weights (q_point + n_edge_points,
                                                 2 * (i * deg + j) + 1)
                               * tmp;
                    }

                  system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      if (std::abs (solution (i * deg + j)) > 1e-14)
                        local_dofs[((2 * face + 1) * deg + j + GeometryInfo<dim>::lines_per_cell)
                                   * this->degree + i]
                          = solution (i * deg + j);
                }

                                  // Finally we project
                                  // the remaining parts
                                  // of the function on
                                  // the interior shape
                                  // functions.
              const QGauss<dim> reference_quadrature (this->degree);
              const unsigned int
                n_interior_points = reference_quadrature.size ();

                                  // We create the
                                  // system matrix.
              system_matrix.reinit (this->degree * deg * deg,
                                    this->degree * deg * deg);
              system_matrix = 0;

              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k < deg; ++k)
                    for (unsigned int l = 0; l <= deg; ++l)
                      for (unsigned int m = 0; m < deg; ++m)
                        for (unsigned int n = 0; n < deg; ++n)
                          for (unsigned int q_point = 0;
                               q_point < n_interior_points; ++q_point)
                            system_matrix ((i * deg + j) * deg + k,
                                           (l * deg + m) * deg + n)
                              += reference_quadrature.weight (q_point)
                                 * legendre_polynomials[i].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points]
                                    (0))
                                 * lobatto_polynomials[j + 2].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points]
                                    (1))
                                 * lobatto_polynomials[k + 2].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points]
                                    (2))
                                 * lobatto_polynomials_grad[l].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points]
                                    (0))
                                 * lobatto_polynomials[m + 2].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points]
                                    (1))
                                 * lobatto_polynomials[n + 2].value
                                   (this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points]
                                    (2));

              system_matrix_inv.reinit (system_matrix.m (),
                                        system_matrix.m ());
              system_matrix_inv.invert (system_matrix);
                                  // Set up the right hand side.
              system_rhs.reinit (system_matrix.m ());
              system_rhs = 0;

              for (unsigned int q_point = 0; q_point < n_interior_points;
                   ++q_point)
                {
                  double tmp
                    = values[0][q_point + GeometryInfo<dim>::lines_per_cell
                                * n_edge_points
                                + GeometryInfo<dim>::faces_per_cell
                                * n_face_points];

                  for (unsigned int i = 0; i <= deg; ++i)
                    {
                      for (unsigned int j = 0; j < 2; ++j)
                        for (unsigned int k = 0; k < 2; ++k)
                          tmp -= local_dofs[i + (j + 4 * k + 2) * this->degree]
                                 * this->shape_value_component
                                   (i + (j + 4 * k + 2) * this->degree,
                                    this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points],
                                    0);

                      for (unsigned int j = 0; j < deg; ++j)
                        for (unsigned int k = 0; k < 4; ++k)
                          tmp -= local_dofs[(i + 2 * (k + 2) * this->degree
                                             + GeometryInfo<dim>::lines_per_cell)
                                            * deg + j
                                            + GeometryInfo<dim>::lines_per_cell]
                                 * this->shape_value_component
                                   ((i + 2 * (k + 2) * this->degree
                                     + GeometryInfo<dim>::lines_per_cell)
                                    * deg + j
                                    + GeometryInfo<dim>::lines_per_cell,
                                   this->generalized_support_points[q_point
                                                                    + GeometryInfo<dim>::lines_per_cell
                                                                    * n_edge_points
                                                                    + GeometryInfo<dim>::faces_per_cell
                                                                    * n_face_points],
                                   0);
                    }

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      for (unsigned int k = 0; k < deg; ++k)
                        system_rhs ((i * deg + j) * deg + k)
                          += reference_quadrature.weight (q_point) * tmp
                             * lobatto_polynomials_grad[i].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points
                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                 * n_face_points]
                                (0))
                             * lobatto_polynomials[j + 2].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points
                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                 * n_face_points]
                                (1))
                             * lobatto_polynomials[k + 2].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points
                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                 * n_face_points]
                                (2));
                }

              solution.reinit (system_rhs.size ());
              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k < deg; ++k)
                    if (std::abs (solution ((i * deg + j) * deg + k)) > 1e-14)
                      local_dofs[((i + 2 * GeometryInfo<dim>::faces_per_cell)
                                  * deg + j + GeometryInfo<dim>::lines_per_cell
                                  + 2 * GeometryInfo<dim>::faces_per_cell)
                                 * deg + k + GeometryInfo<dim>::lines_per_cell]
                        = solution ((i * deg + j) * deg + k);

                                  // Set up the right hand side.
              system_rhs = 0;

              for (unsigned int q_point = 0; q_point < n_interior_points;
                   ++q_point)
                {
                  double tmp
                    = values[1][q_point + GeometryInfo<dim>::lines_per_cell
                                * n_edge_points
                                + GeometryInfo<dim>::faces_per_cell
                                * n_face_points];

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < 2; ++j)
                      {
                        for (unsigned int k = 0; k < 2; ++k)
                          tmp -= local_dofs[i + (4 * j + k) * this->degree]
                                 * this->shape_value_component
                                   (i + (4 * j + k) * this->degree,
                                    this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points],
                                    1);

                        for (unsigned int k = 0; k < deg; ++k)
                          tmp -= local_dofs[(i + 2 * j * this->degree
                                             + GeometryInfo<dim>::lines_per_cell)
                                            * deg + k
                                            + GeometryInfo<dim>::lines_per_cell]
                                 * this->shape_value_component
                                   ((i + 2 * j * this->degree
                                     + GeometryInfo<dim>::lines_per_cell)
                                    * deg + k
                                    + GeometryInfo<dim>::lines_per_cell,
                                    this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points],
                                    1)
                                 + local_dofs[i + ((2 * j + 9) * deg + k
                                                   + GeometryInfo<dim>::lines_per_cell)
                                              * this->degree]
                                 * this->shape_value_component
                                   (i + ((2 * j + 9) * deg + k
                                          + GeometryInfo<dim>::lines_per_cell)
                                    * this->degree,
                                    this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points],
                                    1);
                      }

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      for (unsigned int k = 0; k < deg; ++k)
                        system_rhs ((i * deg + j) * deg + k)
                          += reference_quadrature.weight (q_point) * tmp
                             * lobatto_polynomials_grad[i].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points
                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                 * n_face_points]
                                (1))
                             * lobatto_polynomials[j + 2].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points
                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                 * n_face_points]
                                (0))
                             * lobatto_polynomials[k + 2].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points
                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                 * n_face_points]
                                (2));
                }

              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k < deg; ++k)
                    if (std::abs (solution ((i * deg + j) * deg + k)) > 1e-14)
                      local_dofs[((i + this->degree + 2
                                   * GeometryInfo<dim>::faces_per_cell) * deg
                                  + j + GeometryInfo<dim>::lines_per_cell + 2
                                  * GeometryInfo<dim>::faces_per_cell) * deg
                                  + k + GeometryInfo<dim>::lines_per_cell]
                        = solution ((i * deg + j) * deg + k);

                                  // Set up the right hand side.
              system_rhs = 0;

              for (unsigned int q_point = 0; q_point < n_interior_points;
                   ++q_point)
                {
                  double tmp
                    = values[2][q_point + GeometryInfo<dim>::lines_per_cell
                                * n_edge_points
                                + GeometryInfo<dim>::faces_per_cell
                                * n_face_points];

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < 4; ++j)
                      {
                        tmp -= local_dofs[i + (j + 8) * this->degree]
                               * this->shape_value_component
                                 (i + (j + 8) * this->degree,
                                  this->generalized_support_points[q_point
                                                                   + GeometryInfo<dim>::lines_per_cell
                                                                   * n_edge_points
                                                                   + GeometryInfo<dim>::faces_per_cell
                                                                   * n_face_points],
                                  2);

                        for (unsigned int k = 0; k < deg; ++k)
                          tmp -= local_dofs[i + ((2 * j + 1) * deg + k
                                                 + GeometryInfo<dim>::lines_per_cell)
                                            * this->degree]
                                 * this->shape_value_component
                                   (i + ((2 * j + 1) * deg + k
                                         + GeometryInfo<dim>::lines_per_cell)
                                    * this->degree,
                                    this->generalized_support_points[q_point
                                                                     + GeometryInfo<dim>::lines_per_cell
                                                                     * n_edge_points
                                                                     + GeometryInfo<dim>::faces_per_cell
                                                                     * n_face_points],
                                    2);
                      }

                  for (unsigned int i = 0; i <= deg; ++i)
                    for (unsigned int j = 0; j < deg; ++j)
                      for (unsigned int k = 0; k < deg; ++k)
                        system_rhs ((i * deg + j) * deg + k)
                          += reference_quadrature.weight (q_point) * tmp
                             * lobatto_polynomials_grad[i].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points
                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                 * n_face_points]
                                (2))
                             * lobatto_polynomials[j + 2].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points
                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                 * n_face_points]
                                (0))
                             * lobatto_polynomials[k + 2].value
                               (this->generalized_support_points[q_point
                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                 * n_edge_points
                                                                 + GeometryInfo<dim>::faces_per_cell
                                                                 * n_face_points]
                                (1));
                }

              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k < deg; ++k)
                    if (std::abs (solution ((i * deg + j) * deg + k)) > 1e-14)
                      local_dofs[i + ((j + 2 * (deg
                                                + GeometryInfo<dim>::faces_per_cell))
                                      * deg + k
                                      + GeometryInfo<dim>::lines_per_cell)
                                 * this->degree]
                        = solution ((i * deg + j) * deg + k);
            }

          break;
        }

      default:
        Assert (false, ExcNotImplemented ());
    }
}


template <int dim>
std::size_t
FE_Nedelec<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}


//----------------------------------------------------------------------//


// explicit instantiations
#include "fe_nedelec.inst"


DEAL_II_NAMESPACE_CLOSE
