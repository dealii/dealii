#include <base/logstream.h>
#include <base/utilities.h>
#include <base/quadrature.h>
#include <base/quadrature_lib.h>
#include <base/qprojector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/mapping.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_nothing.h>
#include <fe/fe_values.h>
#include <fe/fe_tools.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <sstream>
#include <iostream>

//TODO: implement the adjust_quad_dof_index_for_face_orientation_table and
//adjust_line_dof_index_for_line_orientation_table fields, and write tests
//similar to bits/face_orientation_and_fe_q_*


DEAL_II_NAMESPACE_OPEN


template <int dim>
FE_Nedelec<dim>::FE_Nedelec (const unsigned int p) :
FE_PolyTensor<PolynomialsNedelec<dim>, dim>
(p,
 FiniteElementData<dim> (get_dpo_vector (p), dim, p + 1,
                         FiniteElementData<dim>::Hcurl, 1),
 std::vector<bool> (PolynomialsNedelec<dim>::compute_n_pols (p), true),
 std::vector<std::vector<bool> >
 (PolynomialsNedelec<dim>::compute_n_pols (p),
  std::vector<bool> (dim, true))),
deg (p)
{
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
  this->reinit_restriction_and_prolongation_matrices ();
				   // Fill prolongation matrices with embedding operators
  FETools::compute_embedding_matrices (*this, this->prolongation);
  initialize_restriction ();
  
  switch (dim)
    {
      case 2:
	    this->interface_constraints.reinit (
	      GeometryInfo<dim>::max_children_per_face * this->dofs_per_face,
	      this->dofs_per_face);
	    break;
      case 3:
	    this->interface_constraints.reinit (
	      GeometryInfo<dim>::max_children_per_face * this->dofs_per_face
	      - 4 * this->dofs_per_line,
	      this->dofs_per_face);
	    break;
      default:
	    Assert(false, ExcNotImplemented());
	    break;
    }

  if (dim==2)
    {
      FullMatrix<double> face_embeddings[GeometryInfo<dim>::max_children_per_face];
      
      for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_face; ++i)
	face_embeddings[i].reinit (this->dofs_per_face, this->dofs_per_face);
      
      FETools::compute_face_embedding_matrices<dim,double>
	(*this, face_embeddings, 0, 0);
      unsigned int target_row = 0;
      
      for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_face; ++i)
	for (unsigned int j = 0; j < face_embeddings[i].m (); ++j)
	  {
	    for (unsigned int k = 0; k < face_embeddings[i].n (); ++k)
	      this->interface_constraints (target_row, k)
		= face_embeddings[i] (j, k);
	    
	    ++target_row;
	  }
    }
  else
    this->interface_constraints.reinit(0,0);
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
  namebuf << "FE_Nedelec<" << dim << ">(" << this->tensor_degree()-1 << ")";

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


#if deal_II_dimension == 1

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

#elif deal_II_dimension == 2

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

#else

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

#endif


#if deal_II_dimension == 1

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

#endif

                   // Restriction operator
template <int dim>
void
FE_Nedelec<dim>::initialize_restriction ()
{
				   // To save some computation time we just
				   // put in the correct values, which can
				   // be calculated by projection-based
				   // interpolation.
  switch (dim)
    {
      case 2:
        {
          const unsigned int n_boundary_dofs
            = GeometryInfo<dim>::lines_per_cell * this->degree;

          for (unsigned int ref = RefinementCase<dim>::cut_x;
               ref <= RefinementCase<dim>::isotropic_refinement; ++ref)
            {
              const unsigned int index = ref - 1;

              switch (ref)
                {
                  case RefinementCase<dim>::cut_x:
                    {
                      for (unsigned int i = 0; i <= deg; ++i)
                        {
                          for (unsigned int j = 0; j < 2; ++j)
                            this->restriction[index][j] (i + j * this->degree,
                                                         i + j * this->degree)
                              = 2.0;

                          for (unsigned int j = 2;
                               j < GeometryInfo<dim>::lines_per_cell; ++j)
                            for (unsigned int k = 0; k < 2; ++k)
                              this->restriction[index][k]
                              (i + j * this->degree, i + j * this->degree)
                                = 1.0;

                          for (unsigned int j = 0; j < deg; ++j)
                            for (unsigned int k = 0; k < 2; ++k)
                              for (unsigned int child = 0;
                                   child < GeometryInfo<dim>::n_children
                                           (RefinementCase<dim> (ref));
                                   ++ child)
                                this->restriction[index][child]
                                ((i + k * this->degree) * deg + j
                                 + n_boundary_dofs,
                                 (i + k * this->degree) * deg + j
                                 + n_boundary_dofs) = 1.0;
                        }

                      break;
                    }

                  case RefinementCase<dim>::cut_y:
                    {
                      for (unsigned int i = 0; i < this->degree; ++i)
                        {
                          for (unsigned int j = 0; j < 2; ++j)
                            {
                              for (unsigned int k = 0; k < 2; ++k)
                                this->restriction[index][k]
                                (i + j * this->degree, i + j * this->degree)
                                  = 1.0;

                              this->restriction[index][j]
                              (i + (j + 2) * this->degree,
                               i + (j + 2) * this->degree) = 2.0;
                            }

                          for (unsigned int j = 0; j < deg; ++j)
                            for (unsigned int k = 0; k < 2; ++k)
                              for (unsigned int child = 0;
                                   child < GeometryInfo<dim>::n_children
                                           (RefinementCase<dim> (ref));
                                   ++ child)
                                this->restriction[index][child]
                                ((i + k * this->degree) * deg + j
                                 + n_boundary_dofs,
                                 (i + k * this->degree) * deg + j
                                 + n_boundary_dofs) = 1.0;
                        }

                      break;
                    }

                  case RefinementCase<dim>::isotropic_refinement:
                    {
                      for (unsigned int i = 0; i < this->degree; ++i)
                        {
                          for (unsigned int j = 0; j < 2; ++j)
                            {
                              this->restriction[index][j]
                              (i + j * this->degree, i + j * this->degree)
                                = 1.0;
                              this->restriction[index][j]
                              (i + 2 * this->degree, i + 2 * this->degree)
                                = 1.0;
                              this->restriction[index][j + 2]
                              (i + j * this->degree, i + j * this->degree)
                                 = 1.0;
                              this->restriction[index][j + 2]
                              (i + 3 * this->degree, i + 3 * this->degree)
                                 = 1.0;
                            }

                          for (unsigned int j = 0; j < deg; ++j)
                            for (unsigned int k = 0; k < 2; ++k)
                              for (unsigned int child = 0;
                                   child < GeometryInfo<dim>::n_children
                                           (RefinementCase<dim> (ref));
                                   ++ child)
                                this->restriction[index][child]
                                ((i + k * this->degree) * deg + j
                                 + n_boundary_dofs,
                                 (i + k * this->degree) * deg + j
                                 + n_boundary_dofs) = 0.5;
                        }

                      break;
                    }

                  default:
                    Assert (false, ExcNotImplemented ());
                }
            }

          break;
        }

      case 3:
        {
          const unsigned int n_edge_dofs
            = GeometryInfo<dim>::lines_per_cell * deg;
          const unsigned int n_boundary_dofs
            = n_edge_dofs
              + 2 * GeometryInfo<dim>::faces_per_cell * deg * this->degree;

          for (unsigned int ref = RefinementCase<dim>::cut_x;
               ref <= RefinementCase<dim>::isotropic_refinement; ++ref)
            {
              const unsigned int index = ref - 1;

              switch (ref)
                {
                  case RefinementCase<3>::cut_x:
                    {
                      for (unsigned int i = 0; i <= deg; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          {
                            this->restriction[index][j] (i + j * this->degree,
                                                         i + j * this->degree)
                              = 2.0;
                            this->restriction[index][j] (i + 2 * this->degree,
                                                         i + 2 * this->degree)
                              = 1.0;
                            this->restriction[index][j] (i + 3 * this->degree,
                                                         i + 3 * this->degree)
                              = 1.0;
                            this->restriction[index][j]
                            (i + (j + 4) * this->degree,
                             i + (j + 4) * this->degree) = 2.0;
                            this->restriction[index][j] (i + 6 * this->degree,
                                                         i + 6 * this->degree)
                              = 1.0;
                            this->restriction[index][j] (i + 7 * this->degree,
                                                         i + 7 * this->degree)
                              = 1.0;
                            this->restriction[index][j]
                            (i + (j + 8) * this->degree,
                             i + (j + 8) * this->degree) = 2.0;
                            this->restriction[index][j]
                            (i + (j + 10) * this->degree,
                             i + (j + 10) * this->degree) = 2.0;
                          }

                      for (unsigned int i = 0; i < 2 * this->degree * deg; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          {
                            this->restriction[index][j]
                            (i + j * this->degree * deg + n_edge_dofs,
                             i + j * this->degree + deg + n_edge_dofs) = 2.0;

                            for (unsigned int k = 0; k < 4; ++k)
                              this->restriction[index][j]
                              (i + (2 * k + 4) * this->degree * deg
                                 + n_edge_dofs,
                               i + (2 * k + 4) * this->degree * deg
                                 + n_edge_dofs) = 1.0;
                          }

                      break;
                    }

                  case RefinementCase<3>::cut_y:
                    {
                      for (unsigned int i = 0; i <= deg; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          {
                            this->restriction[index][j] (i, i) = 1.0;
                            this->restriction[index][j] (i + this->degree,
                                                         i + this->degree)
                              = 1.0;
                            this->restriction[index][j]
                            (i + (j + 2) * this->degree,
                             i + (j + 2) * this->degree) = 2.0;
                            this->restriction[index][j] (i + 4 * this->degree,
                                                         i + 4 * this->degree)
                              = 1.0;
                            this->restriction[index][j] (i + 5 * this->degree,
                                                         i + 5 * this->degree)
                              = 1.0;

                            for (unsigned int k = 3; k < 6; ++k)
                              this->restriction[index][j]
                              (i + (j + 2 * k) * this->degree,
                               i + (j + 2 * k) * this->degree) = 2.0;
                          }

                      for (unsigned int i = 0; i < 2 * this->degree * deg; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          {
                            this->restriction[index][j] (i + n_edge_dofs,
                                                         i + n_edge_dofs)
                              = 1.0;
                            this->restriction[index][j]
                            (i + 2 * this->degree * deg + n_edge_dofs,
                             i + 2 * this->degree * deg + n_edge_dofs) = 1.0;
                            this->restriction[index][j]
                            (i + (2 * j + 4) * this->degree * deg
                               + n_edge_dofs,
                             i + (2 * j + 4) * this->degree * deg
                               + n_edge_dofs) = 2.0;
                            this->restriction[index][j]
                            (i + 8 * this->degree * deg + n_edge_dofs,
                             i + 8 * this->degree * deg + n_edge_dofs) = 1.0;
                            this->restriction[index][j]
                            (i + 10 * this->degree * deg + n_edge_dofs,
                             i + 10 * this->degree * deg + n_edge_dofs) = 1.0;
                          }

                      break;
                    }

                  case RefinementCase<3>::cut_xy:
                    {
                      for (unsigned int i = 0; i <= deg; ++i)
                        {
                          for (unsigned int j = 0; j < 2; ++j)
                            {
                              this->restriction[index][2 * j] (i, i) = 1.0;
                              this->restriction[index][2 * j + 1]
                              (i + this->degree, i + this->degree) = 1.0;
                              this->restriction[index][j]
                              (i + 2 * this->degree, i + 2 * this->degree)
                                = 1.0;
                              this->restriction[index][j + 2]
                              (i + 3 * this->degree, i + 3 * this->degree)
                                = 1.0;
                              this->restriction[index][2 * j]
                              (i + 4 * this->degree, i + 4 * this->degree)
                                = 1.0;
                              this->restriction[index][2 * j + 1]
                              (i + 5 * this->degree, i + 5 * this->degree)
                                = 1.0;
                              this->restriction[index][j]
                              (i + 6 * this->degree, i + 6 * this->degree)
                                = 1.0;
                              this->restriction[index][j + 2]
                              (i + 7 * this->degree, i + 7 * this->degree)
                                = 1.0;
                            }

                          for (unsigned int j = 0; j < 4; ++j)
                            this->restriction[index][j]
                            (i + (j + 8) * this->degree,
                             i + (j + 8) * this->degree) = 2.0;
                        }

                      for (unsigned int i = 0; i < 2 * this->degree * deg; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          {
                            this->restriction[index][2 * j] (i + n_edge_dofs,
                                                             i + n_edge_dofs)
                              = 1.0;
                            this->restriction[index][2 * j + 1]
                            (i + 2 * this->degree * deg + n_edge_dofs,
                             i + 2 * this->degree * deg + n_edge_dofs) = 1.0;
                            this->restriction[index][j]
                            (i + 4 * this->degree * deg + n_edge_dofs,
                             i + 4 * this->degree * deg + n_edge_dofs) = 1.0;
                            this->restriction[index][j + 2]
                            (i + 6 * this->degree * deg + n_edge_dofs,
                             i + 6 * this->degree * deg + n_edge_dofs) = 1.0;

                            for (unsigned int k = 0; k < 4; ++k)
                              this->restriction[index][k]
                              (i + 2 * (j + 4) * this->degree * deg
                                 + n_edge_dofs,
                               i + 2 * (j + 4) * this->degree * deg
                                 + n_edge_dofs) = 0.5;
                          }

                      break;
                    }

                  case RefinementCase<3>::cut_z:
                    {
                      for (unsigned int i = 0; i <= deg; ++i)
                        for (unsigned int j = 0; j < 4; ++j)
                          for (unsigned int k = 0; k < 2; ++k)
                            {
                              this->restriction[index][k]
                              (i + (j + 4 * k) * this->degree,
                               i + (j + 4 * k) * this->degree) = 2.0;
                              this->restriction[index][k]
                              (i + (j + 8) * this->degree,
                               i + (j + 8) * this->degree) = 1.0;
                            }

                      for (unsigned int i = 0; i < 2 * this->degree * deg; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          {
                            for (unsigned int k = 0; k < 3; ++k)
                              this->restriction[index][j]
                              (i + 2 * k * this->degree * deg + n_edge_dofs,
                               i + 2 * k * this->degree * deg + n_edge_dofs)
                                = 1.0;

                            this->restriction[index][j]
                            (i + 2 * (j + 4) * this->degree * deg
                               + n_edge_dofs,
                             i + 2 * (j + 4) * this->degree * deg
                               + n_edge_dofs) = 2.0;
                          }

                      break;
                    }

                  case RefinementCase<3>::cut_xz:
                    {
                      for (unsigned int i = 0; i <= deg; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          {
                            this->restriction[index][j] (i + j * this->degree,
                                                         i + j * this->degree)
                              = 2.0;
                            this->restriction[index][j + 2]
                            (i + (j + 4) * this->degree,
                             i + (j + 4) * this->degree) = 2.0;

                            for (unsigned int k = 0; k < 2; ++k)
                              {
                                this->restriction[index][j]
                                (i + (k + 2) * this->degree,
                                 i + (k + 2) * this->degree) = 1.0;
                                this->restriction[index][j + 2]
                                (i + (k + 6) * this->degree,
                                 i + (k + 6) * this->degree) = 1.0;
                                this->restriction[index][2 * j]
                                (i + 2 * (k + 4) * this->degree,
                                 i + 2 * (k + 4) * this->degree) = 1.0;
                                this->restriction[index][2 * j + 1]
                                (i + (2 * k + 9) * this->degree,
                                 i + (2 * k + 9) * this->degree) = 1.0;
                              }
                          }

                      for (unsigned int i = 0; i < 2 * this->degree * deg; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          {
                            this->restriction[index][2 * j] (i + n_edge_dofs,
                                                             i + n_edge_dofs)
                              = 1.0;
                            this->restriction[index][2 * j + 1]
                            (i + 2 * this->degree * deg + n_edge_dofs,
                             i + 2 * this->degree * deg + n_edge_dofs) = 1.0;

                            for (unsigned int k = 0; k < 2; ++k)
                              {
                                this->restriction[index][j + 2 * k]
                                (i + 4 * this->degree * deg + n_edge_dofs,
                                 i + 4 * this->degree * deg + n_edge_dofs)
                                   = 0.5;
                                this->restriction[index][j + 2 * k]
                                (i + 6 * this->degree * deg + n_edge_dofs,
                                 i + 6 * this->degree * deg + n_edge_dofs)
                                   = 0.5;
                                this->restriction[index][j + 2 * k]
                                (i + 2 * (k + 4) * this->degree * deg
                                   + n_edge_dofs,
                                 i + 2 * (k + 4) * this->degree * deg
                                   + n_edge_dofs) = 1.0;
                              }
                          }

                      break;
                    }

                  case RefinementCase<3>::cut_yz:
                    {
                      for (unsigned int i = 0; i <= deg; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          for (unsigned int k = 0; k < 2; ++k)
                            {
                          	  for (unsigned int l = 0; l < 2; ++l)
                          	    {
                                  this->restriction[index][j + 2 * l]
                                  (i + (k + 4 * l) * this->degree,
                                   i + (k + 4 * l) * this->degree) = 1.0;
                                  this->restriction[index][2 * j + l]
                                  (i + (k + 2 * (l + 4)) * this->degree,
                                   i + (k + 2 * (l + 4)) * this->degree) = 1.0;
                          	    }

                              this->restriction[index][j + 2 * k]
                              (i + (j + 4 * k + 2) * this->degree,
                               i + (j + 4 * k + 2) * this->degree) = 2.0;
                          }

                      for (unsigned int i = 0; i < 2 * this->degree * deg; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          {
                            for (unsigned int child = 0;
                                 child < GeometryInfo<dim>::n_children
                                         (RefinementCase<dim> (ref)); ++child)
                              this->restriction[index][child]
                              (i + 2 * j * this->degree * deg + n_edge_dofs,
                               i + 2 * j * this->degree * deg + n_edge_dofs)
                                = 0.5;

                            for (unsigned int k = 0; k < 2; ++k)
                              {
                                this->restriction[index][j + 2 * k]
                                (i + 2 * (j + 2) * this->degree * deg
                                   + n_edge_dofs,
                                 i + 2 * (j + 2) * this->degree * deg
                                   + n_edge_dofs) = 1.0;
                                this->restriction[index][2 * j + k]
                                (i + 2 * (j + 4) * this->degree * deg
                                   + n_edge_dofs,
                                 i + 2 * (j + 4) * this->degree * deg
                                   + n_edge_dofs) = 1.0;
                              }
                          }

                      break;
                    }

                  case RefinementCase<3>::isotropic_refinement:
                    {
                      for (unsigned int i = 0; i <= deg; ++i)
                        for (unsigned int j = 0; j < 2; ++j)
                          {
                            for (unsigned int k = 0; k < 2; ++k)
                              {
                                this->restriction[index][2 * j + k]
                                (i + k * this->degree, i + k * this->degree)
                                  = 1.0;
                                this->restriction[index][j + 2 * k]
                                (i + (k + 2) * this->degree,
                                 i + (k + 2) * this->degree) = 1.0;
                                this->restriction[index][2 * (j + 2) + k]
                                (i + (k + 4) * this->degree,
                                 i + (k + 4) * this->degree) = 1.0;
                                this->restriction[index][j + 2 * (k + 2)]
                                (i + (k + 6) * this->degree,
                                 i + (k + 6) * this->degree) = 1.0;
                                this->restriction[index][4 * j + k]
                                (i + (k + 8) * this->degree,
                                 i + (k + 8) * this->degree) = 1.0;
                              }

                            this->restriction[index][2 * (2 * j + 1)]
                            (i + 10 * this->degree, i + 10 * this->degree)
                              = 1.0;
                            this->restriction[index][4 * j + 3]
                            (i + 11 * this->degree, i + 11 * this->degree)
                              = 1.0;
                          }

                      for (unsigned int i = 0; i < 2 * this->degree * deg; ++i)
                        {
                          for (unsigned int j = 0; j < 4; ++j)
                            {
                              this->restriction[index][2 * j] (i + n_edge_dofs,
                                                               i + n_edge_dofs)
                                = 0.5;
                              this->restriction[index][2 * j + 1]
                              (i + 2 * this->degree * deg + n_edge_dofs,
                               i + 2 * this->degree * deg + n_edge_dofs) = 0.5;
                              this->restriction[index][j]
                              (i + 8 * this->degree * deg + n_edge_dofs,
                               i + 8 * this->degree * deg + n_edge_dofs) = 0.5;
                              this->restriction[index][j + 4]
                              (i + 10 * this->degree * deg + n_edge_dofs,
                               i + 10 * this->degree * deg + n_edge_dofs)
                                 = 0.5;
                            }

                          for (unsigned int j = 0; j < 2; ++j)
                            for (unsigned int k = 0; k < 2; ++k)
                              for (unsigned int l = 0; l < 2; ++l)
                                this->restriction[index][j + 2 * (2 * k + l)]
                                (i + 2 * (l + 2) * this->degree * deg
                                   + n_edge_dofs,
                                 i + 2 * (l + 2) * this->degree * deg
                                   + n_edge_dofs) = 0.5;
                        }

                      break;
                    }

                  default:
                    Assert (false, ExcNotImplemented ());
                }

              for (unsigned int i = 0; i < 3 * this->degree * deg * deg; ++i)
                for (unsigned int child = 0;
                     child < GeometryInfo<dim>::n_children
                             (RefinementCase<dim> (ref)); ++child)
                  this->restriction[index][child] (i + n_boundary_dofs,
                                                   i + n_boundary_dofs)
                    = 2.0 / GeometryInfo<dim>::n_children
                            (RefinementCase<dim> (ref));
            }

          break;
        }

      default:
        Assert (false, ExcNotImplemented ());
    }
}


#if deal_II_dimension == 1

template <>
std::vector<unsigned int>
FE_Nedelec<1>::get_dpo_vector (const unsigned int degree)
{
  std::vector<unsigned int> dpo (2);

  dpo[0] = 1;
  dpo[1] = degree;
  return dpo;
}

#endif


template <int dim>
std::vector<unsigned int>
FE_Nedelec<dim>::get_dpo_vector (const unsigned int degree)
{
  std::vector<unsigned int> dpo (dim + 1);

  dpo[0] = 0;
  dpo[1] = degree + 1;
  dpo[2] = 2 * degree * (degree + 1);

  if (dim == 3)
     dpo[3] = 3 * degree * degree * (degree + 1);

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
                      < (GeometryInfo<3>::lines_per_cell + 5 * deg)
                        * this->degree)) ||
                  ((shape_index
                      >= (GeometryInfo<3>::lines_per_cell + 6 * deg)
                         * this->degree) &&
                   (shape_index
                      < (GeometryInfo<3>::lines_per_cell + 7 * deg)
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
                       < (GeometryInfo<3>::lines_per_cell + 8 * deg)
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
                      >= (GeometryInfo<3>::lines_per_cell + 4 * deg)
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
                      >= (GeometryInfo<3>::lines_per_cell + 8 * deg)
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
                      >= (GeometryInfo<3>::lines_per_cell + 4 * deg)
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
bool
FE_Nedelec<dim>::hp_constraints_are_implemented () const
{
  return dim != 2;
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
            identities.push_back (std::make_pair ((i + 1) * (q + 1) + j,
                                                  (i + 1) * (p + 1) + j));
            identities.push_back (std::make_pair (i + (j + q + 2) * q,
                                                  i + (j + p + 2) * p));
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
  for (unsigned int i = 0; i <= deg; ++i)
    interpolation_matrix (i, i) = 1.0;

                   // In 3d we have some lines more
                   // and a face. The procedure stays
                   // the same as above, but we have
                   // to take a bit more care of the
                   // indices of the degrees of
                   // freedom.
  if (dim == 3)
    for (unsigned int i = 0; i <= deg; ++i)
      {
        for (int j = 1; j < (int) GeometryInfo<dim>::lines_per_face; ++j)
          interpolation_matrix (j * source_fe.degree + i,
                                j * this->degree + i) = 1.0;

        for (unsigned int j = 0; j < deg; ++j)
          {
            interpolation_matrix
              (i + (j + GeometryInfo<2>::lines_per_cell) * source_fe.degree,
               i + (j + GeometryInfo<2>::lines_per_cell) * this->degree)
              = 1.0;
            interpolation_matrix
              ((i * (source_fe.degree - 1)
               + GeometryInfo<2>::lines_per_cell) * source_fe.degree + j,
               (i * deg + GeometryInfo<2>::lines_per_cell) * this->degree)
              = 1.0;
          }
      }
}

#if deal_II_dimension == 1

template <int dim>
void
FE_Nedelec<dim>::get_subface_interpolation_matrix(
  const FiniteElement<dim>&,
  const unsigned int,
  FullMatrix<double>&) const
{
  Assert (false, ExcNotImplemented ());
}

#else

                   // In this function we compute the
                   // subface interpolation matrix.
                   // This is done by a projection-
                   // based interpolation. Therefore
                   // we first interpolate the
                   // shape functions of the higher
                   // order element on the lowest
                   // order edge shape functions.
                   // Then the remaining part of
                   // the interpolated shape
                   // functions is projected on the
                   // higher order edge shape
                   // functions, the face shape
                   // functions and the interior
                   // shape functions (if they all
                   // exist).
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
  interpolation_matrix = 0;
		           // Perform projection-based interpolation
				   // as usual.
  switch (dim)
    {
      case 2:
        {
          const QGauss<dim - 1> reference_edge_quadrature (this->degree);
          const Quadrature<dim - 1>& edge_quadrature
            = QProjector<dim - 1>::project_to_child
              (reference_edge_quadrature, subface);
          const unsigned int& n_edge_points = edge_quadrature.size ();
          const std::vector<Point<dim - 1> >&
            quadrature_points = edge_quadrature.get_points ();

                                  // Let us begin with the
                                  // interpolation part.
          for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
            {
              const double weight = edge_quadrature.weight (q_point);

              for (unsigned int dof = 0; dof < this->dofs_per_face; ++dof)
                interpolation_matrix (0, dof)
                  += weight
                     * this->shape_value_component
                       (dof, Point<dim> (0.0, quadrature_points[q_point] (0)),
                        1);
            }

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
          for (unsigned int dof = 0; dof < this->dofs_per_face; ++dof)
            if (std::abs (interpolation_matrix (0, dof)) < 1e-14)
              interpolation_matrix (0, dof) = 0.0;

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some higher order edge
                                  // shape functions to
                                  // consider.
                                  // Here the projection part
                                  // starts. The dof values
                                  // are obtained by solving
                                  // a linear system of
                                  // equations.
          if (deg > 0)
            {
            	                  // Shift value for scaling
            	                  // of quadrature points.
              const double shift[2] = {0.0, -1.0};
              const std::vector<Polynomials::Polynomial<double> >&
                lobatto_polynomials
                  = Polynomials::Lobatto::generate_complete_basis
                    (this->degree);
              FullMatrix<double> assembling_matrix (deg, n_edge_points);
              std::vector<Polynomials::Polynomial<double> >
                lobatto_polynomials_grad (this->degree);

              for (unsigned int i = 0; i < lobatto_polynomials_grad.size ();
                   ++i)
                lobatto_polynomials_grad[i]
                  = lobatto_polynomials[i + 1].derivative ();

                                  // Set up the system matrix
                                  // and right hand side
                                  // vector.
              for (unsigned int q_point = 0; q_point < n_edge_points;
                   ++q_point)
                {
               	  const double tmp = 2.0 * quadrature_points[q_point] (0)
               	                     + shift[subface];
               	  const double weight
               	    = std::sqrt (edge_quadrature.weight (q_point));

                  for (unsigned int i = 0; i < deg; ++i)
                    assembling_matrix (i, q_point)
                      = weight * lobatto_polynomials_grad[i + 1].value (tmp);
               	}

              FullMatrix<double> system_matrix (deg, deg);

              assembling_matrix.mTmult (system_matrix, assembling_matrix);

              FullMatrix<double> system_matrix_inv (deg, deg);

              system_matrix_inv.invert (system_matrix);

              Vector<double> solution (deg);
              Vector<double> system_rhs (deg);

              for (unsigned int dof = 0; dof < this->dofs_per_face; ++dof)
                {
                  system_rhs = 0;

                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                    {
                      const double tmp
                        = 2.0 * quadrature_points[q_point] (0)
                          + shift[subface];
                      const double weight
                        = edge_quadrature.weight (q_point)
                          * (this->shape_value_component
                             (dof, Point<dim> (0.0,
                                               quadrature_points[q_point] (0)),
                              1) - interpolation_matrix (0, dof));

                      for (unsigned int i = 0;  i < deg; ++i)
                        system_rhs (i)
                          += weight
                             * lobatto_polynomials_grad[i + 1].value (tmp);
                    }

                  system_matrix_inv.vmult (solution, system_rhs);

                  for (unsigned int i = 0; i < deg; ++i)
                    if (std::abs (solution (i)) > 1e-14)
                      interpolation_matrix (i + 1, dof) = solution (i);
                }
            }

          break;
        }

      case 3:
        {
          const QGauss<1> reference_edge_quadrature (this->degree);

          switch (subface)
            {
              case 0:
                {
                  const Quadrature<1>& edge_quadrature
                    = QProjector<1>::project_to_child
                      (reference_edge_quadrature, 0);
                  const unsigned int n_edge_points = edge_quadrature.size ();
                  const std::vector<Point<1> >&
                    edge_quadrature_points = edge_quadrature.get_points ();

                                  // Let us begin with the
                                  // interpolation part.
                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                    {
                      const double
                        weight = edge_quadrature.weight (q_point);

                      for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int dof = 0; dof < this->dofs_per_face;
                             ++dof)
                          {
                            interpolation_matrix (i * source_fe.degree, dof)
                              += weight
                                 * this->shape_value_component
                                   (this->face_to_cell_index (dof, 4),
                                    Point<dim>
                                    (0.5 * i,
                                     edge_quadrature_points[q_point] (0), 0.0),
                                     1);
                            interpolation_matrix ((i + 2) * source_fe.degree,
                                                  dof)
                              += weight
                                 * this->shape_value_component
                                   (this->face_to_cell_index (dof, 4),
                                    Point<dim>
                                    (edge_quadrature_points[q_point] (0),
                                     0.5 * i, 0.0), 0);
                          }
                    }

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int dof = 0; dof < this->dofs_per_face;
                         ++dof)
                      {
                        if (std::abs (interpolation_matrix
                                      (i * source_fe.degree, dof)) < 1e-14)
                          interpolation_matrix (i * source_fe.degree, dof)
                            = 0.0;

                        if (std::abs (interpolation_matrix
                                      ((i + 2) * source_fe.degree, dof))
                              < 1e-14)
                          interpolation_matrix ((i + 2) * source_fe.degree,
                                                dof) = 0.0;
                      }

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some higher order edge
                                  // shape functions to
                                  // consider.
                                  // Here the projection part
                                  // starts. The dof values
                                  // are obtained by solving
                                  // a linear system of
                                  // equations.
                  if (deg > 0)
                    {
                    			  // We start with projection
                    			  // on the higher order edge
                    			  // shape function.
                      const QGauss<dim - 1> reference_face_quadrature
                        (this->degree);
                      const Quadrature<dim - 1>& face_quadrature
                        = QProjector<dim - 1>::project_to_child
                          (reference_face_quadrature, 0);
                      const std::vector<Polynomials::Polynomial<double> >&
                        legendre_polynomials
                          = Polynomials::Legendre::generate_complete_basis
                            (deg);
                      const std::vector<Polynomials::Polynomial<double> >&
                        lobatto_polynomials
                          = Polynomials::Lobatto::generate_complete_basis
                            (this->degree);
                      const std::vector<Point<dim - 1> >&
                        face_quadrature_points = face_quadrature.get_points ();
                      const unsigned int& n_face_points
                        = face_quadrature.size ();
                      FullMatrix<double> assembling_matrix
                        (deg, n_edge_points);
                      FullMatrix<double> system_matrix (deg, deg);
                      FullMatrix<double> system_matrix_inv (deg, deg);
                      std::vector<Polynomials::Polynomial<double> >
                        lobatto_polynomials_grad (this->degree);

                      for (unsigned int i = 0; i <= deg; ++i)
                        lobatto_polynomials_grad[i]
                          = lobatto_polynomials[i + 1].derivative ();

//TODO:[Markus Buerg] We should not need those, since the projections
//on each face should just be copies of each other.

                                  // Shifted and scaled
                                  // quadrature points on
                                  // the four edges of a
                                  // face.
                      std::vector<std::vector<Point<dim> > >
                        edge_quadrature_points_full_dim
                        (GeometryInfo<dim>::lines_per_face,
                         std::vector<Point<dim> > (n_edge_points));

                      for (unsigned int q_point = 0; q_point < n_edge_points;
                           ++q_point)
                        {
                          edge_quadrature_points_full_dim[0][q_point]
                            = Point<dim> (0.0,
                                          edge_quadrature_points[q_point] (0),
                                          0.0);
                          edge_quadrature_points_full_dim[1][q_point]
                            = Point<dim> (0.5,
                                          edge_quadrature_points[q_point] (0),
                                          0.0);
                          edge_quadrature_points_full_dim[2][q_point]
                            = Point<dim> (edge_quadrature_points[q_point] (0),
                                          0.0, 0.0);
                          edge_quadrature_points_full_dim[3][q_point]
                            = Point<dim> (edge_quadrature_points[q_point] (0),
                                          0.5, 0.0);
                        }

                                  // Set up the system matrix.
                                  // This can be used for all
                                  // edges.
                  	  for (unsigned int q_point = 0;
                  	       q_point < n_edge_points; ++q_point)
                  	    {
                  	 	  const double tmp
                  	 	    = 2.0 * edge_quadrature_points[q_point] (0);
                  	 	  const double weight
                  	 	    = std::sqrt (edge_quadrature.weight (q_point));

                          for (unsigned int i = 0; i < deg; ++i)
                            assembling_matrix (i, q_point)
                               = weight * lobatto_polynomials_grad[i + 1].value
                                          (tmp);
                  	    }

                      assembling_matrix.mTmult (system_matrix,
                                                assembling_matrix);
                      system_matrix_inv.invert (system_matrix);

                      Vector<double> solution (deg);
                      Vector<double> system_rhs (deg);

                      for (unsigned int dof = 0; dof < this->dofs_per_face;
                           ++dof)
                          for (unsigned int line = 0;
                               line < GeometryInfo<dim>::lines_per_face;
                               ++line)
                            {
                                  // Set up the right hand side.
                     	      system_rhs = 0;

                              for (unsigned int q_point = 0;
                                   q_point < n_edge_points; ++q_point)
                                {
                                  const double right_hand_side_value
                                    = std::sqrt (edge_quadrature.weight
                                                       (q_point))
                                      * (this->shape_value_component
                                         (this->face_to_cell_index (dof, 4),
                                          edge_quadrature_points_full_dim[line][q_point],
                                          1)
                                         - interpolation_matrix
                                           (line * source_fe.degree, dof));
                                  const double tmp
                                    = 2.0 * edge_quadrature_points[q_point] (0);

                                  for (unsigned int i = 0; i < deg; ++i)
                                    system_rhs (i)
                                      += right_hand_side_value
                                         * lobatto_polynomials_grad[i + 1].value
                                           (tmp);
                                }

                              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                              for (unsigned int i = 0; i < deg; ++i)
                                if (std::abs (solution (i)) > 1e-14)
                                  interpolation_matrix
                                    (line * source_fe.degree + i + 1, dof)
                                    = solution (i);
                            }

                      assembling_matrix.reinit (deg * this->degree,
                                                n_face_points);

                      for (unsigned int q_point = 0; q_point < n_face_points;
                           ++q_point)
                        {
                          const Point<dim - 1> quadrature_point
                                (2.0 * face_quadrature_points[q_point] (0),
                                 2.0 * face_quadrature_points[q_point] (1));
                          const double weight
                            = std::sqrt (face_quadrature.weight (q_point));

                          for (unsigned int i = 0; i <= deg; ++i)
                            {
                              const double tmp
                                = weight * legendre_polynomials[i].value
                                           (quadrature_point (0));

                              for (unsigned int j = 0; j < deg; ++j)
                           	    assembling_matrix (i * deg + j, q_point)
                           	      = tmp * lobatto_polynomials[j + 2].value
                           	              (quadrature_point (1));
                            }
                        }

                      system_matrix.reinit (assembling_matrix.m (),
                                            assembling_matrix.m ());
                      assembling_matrix.mTmult (system_matrix,
                                                assembling_matrix);
                      system_matrix_inv.reinit (system_matrix.m (),
                                                system_matrix.m ());
                      system_matrix_inv.invert (system_matrix);
                      solution.reinit (system_matrix_inv.m ());
                      system_rhs.reinit (system_matrix_inv.m ());

                      for (unsigned int dof = 0; dof < this->dofs_per_face;
                           ++dof)
                        {
                          system_rhs = 0;

                                  // Now we project the remaining
                                  // part on the face shape
                                  // functions. First on the
                                  // horizontal ones, then on
                                  // the vertical ones.
                          for (unsigned int q_point = 0;
                               q_point < n_face_points; ++q_point)
                            {
                              const Point<dim> quadrature_point
                                (2.0 * face_quadrature_points[q_point] (0),
                                 2.0 * face_quadrature_points[q_point] (1),
                                 0.0);
                              double right_hand_side_value
                                = this->shape_value_component
                                  (this->face_to_cell_index (dof, 4),
                                   Point<dim>
                                   (face_quadrature_points[q_point] (0),
                                    face_quadrature_points[q_point] (1), 0.0),
                                   1);

                              for (unsigned int i = 0; i < 2; ++i)
                                for (unsigned int j = 0; j < source_fe.degree;
                                     ++j)
                                  right_hand_side_value
                                    -= interpolation_matrix
                                       (i * source_fe.degree + j, dof)
                                       * source_fe.shape_value_component
                                         (i * source_fe.degree + j,
                                          quadrature_point, 1);

                              right_hand_side_value
                                *= face_quadrature.weight (q_point);

                              for (unsigned int i = 0; i <= deg; ++i)
                                {
                                  const double tmp
                                    = right_hand_side_value
                                      * legendre_polynomials[i].value
                                        (quadrature_point (0));

                                  for (unsigned int j = 0; j < deg; ++j)
                                    system_rhs (i * deg + j)
                                      += tmp
                                         * lobatto_polynomials[j + 2].value
                           	               (quadrature_point (1));
                                }
                            }

                          system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                          for (unsigned int i = 0; i <= deg; ++i)
                            for (unsigned int j = 0; j < deg; ++j)
                              if (std::abs (solution (i * deg + j)) > 1e-14)
                                interpolation_matrix
                                  ((i + 4) * source_fe.degree + j - i, dof)
                                  = solution (i * deg + j);

                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.
                          system_rhs = 0;

                          for (unsigned int q_point = 0;
                               q_point < n_face_points; ++q_point)
                            {
                              const Point<dim> quadrature_point
                                (2.0 * face_quadrature_points[q_point] (0),
                                 2.0 * face_quadrature_points[q_point] (1),
                                 0.0);
                              double right_hand_side_value
                                = this->shape_value_component
                                  (this->face_to_cell_index (dof, 4),
                                   Point<dim>
                                   (face_quadrature_points[q_point] (0),
                                    face_quadrature_points[q_point] (1), 0.0),
                                   0);

                              for (unsigned int i = 0; i < 2; ++i)
                                for (unsigned int j = 0; j < source_fe.degree;
                                     ++j)
                                  right_hand_side_value
                                    -= interpolation_matrix
                                       ((i + 2) * source_fe.degree + j, dof)
                                       * source_fe.shape_value_component
                                         (i * source_fe.degree + j,
                                          quadrature_point, 0);

                              right_hand_side_value
                                *= face_quadrature.weight (q_point);

                              for (unsigned int i = 0; i <= deg; ++i)
                                {
                                  const double tmp
                                    = right_hand_side_value
                                      * legendre_polynomials[i].value
                                        (quadrature_point (0));

                                  for (unsigned int j = 0; j < deg; ++j)
                                    system_rhs (i * deg + j)
                                      += tmp
                                         * lobatto_polynomials[j + 2].value
                                           (quadrature_point (1));
                                }
                            }

                          system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                          for (unsigned int i = 0; i <= deg; ++i)
                            for (unsigned int j = 0; j < deg; ++j)
                              if (std::abs (solution (i * deg + j)) > 1e-14)
                                interpolation_matrix
                                (i + (j + source_fe.degree + 3)
                                 * source_fe.degree, dof) = solution (i * deg
                                                                      + j);
                        }
                    }

                  break;
                }

              case 1:
                {
                  const Quadrature<1>& edge_quadrature_x
                    = QProjector<1>::project_to_child
                      (reference_edge_quadrature, 1);
                  const Quadrature<1>& edge_quadrature_y
                    = QProjector<1>::project_to_child
                      (reference_edge_quadrature, 0);
                  const std::vector<Point<1> >&
                    edge_quadrature_x_points = edge_quadrature_x.get_points ();
                  const std::vector<Point<1> >&
                    edge_quadrature_y_points = edge_quadrature_y.get_points ();
                  const unsigned int& n_edge_points
                    = edge_quadrature_x.size ();

                                  // Let us begin with the
                                  // interpolation part.
                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                    {
                      const double weight
                        = edge_quadrature_x.weight (q_point);

                      for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int dof = 0; dof < this->dofs_per_face;
                             ++dof)
                          {
                            interpolation_matrix (i * source_fe.degree, dof)
                              += weight
                                 * this->shape_value_component
                                   (this->face_to_cell_index (dof, 4),
                                    Point<dim>
                                    (0.5 * (i + 1),
                                     edge_quadrature_y_points[q_point] (0), 0.0),
                                    1);
                            interpolation_matrix
                            ((i + 2) * source_fe.degree, dof)
                              += weight
                                 * this->shape_value_component
                                   (this->face_to_cell_index (dof, 4),
                                    Point<dim>
                                    (edge_quadrature_x_points[q_point] (0),
                                     0.5 * i, 0.0), 0);
                          }
                    }

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int dof = 0; dof < this->dofs_per_face;
                         ++dof)
                      {
                        if (std::abs (interpolation_matrix
                                      (i * source_fe.degree, dof)) < 1e-14)
                          interpolation_matrix (i * source_fe.degree, dof)
                            = 0.0;

                        if (std::abs (interpolation_matrix
                                      ((i + 2) * source_fe.degree, dof))
                              < 1e-14)
                          interpolation_matrix ((i + 2) * source_fe.degree,
                                                dof) = 0.0;
                      }

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some higher order edge
                                  // shape functions to
                                  // consider.
                                  // Here the projection part
                                  // starts. The dof values
                                  // are obtained by solving
                                  // a linear system of
                                  // equations.
                  if (deg > 0)
                    {
                    			  // We start with projection
                    			  // on the higher order edge
                    			  // shape function.
                      const QGauss<dim - 1> reference_face_quadrature
                        (this->degree);
                      const Quadrature<dim - 1>& face_quadrature
                        = QProjector<dim - 1>::project_to_child
                          (reference_face_quadrature, 1);
                      const std::vector<Point<dim - 1> >&
                        face_quadrature_points = face_quadrature.get_points ();
                      const std::vector<Polynomials::Polynomial<double> >&
                        legendre_polynomials
                          = Polynomials::Legendre::generate_complete_basis
                            (deg);
                      const std::vector<Polynomials::Polynomial<double> >&
                        lobatto_polynomials
                          = Polynomials::Lobatto::generate_complete_basis
                            (this->degree);
                      const unsigned int&
                        n_face_points = face_quadrature.size ();
                      FullMatrix<double> assembling_matrix (deg,
                                                            n_edge_points);
                      FullMatrix<double> system_matrix (deg, deg);
                      FullMatrix<double> system_matrix_inv (deg, deg);
                      std::vector<Polynomials::Polynomial<double> >
                        lobatto_polynomials_grad (this->degree);

                      for (unsigned int i = 0;
                           i < lobatto_polynomials_grad.size (); ++i)
                        lobatto_polynomials_grad[i]
                          = lobatto_polynomials[i + 1].derivative ();

                                  // Shifted and scaled
                                  // quadrature points and
                                  // weights on the four
                                  // edges of a face.
                      std::vector<std::vector<double> > edge_quadrature_points
                        (GeometryInfo<dim>::lines_per_face,
                         std::vector<double> (n_edge_points));
                      std::vector<std::vector<double> >
                        edge_quadrature_weights
                          (GeometryInfo<dim>::lines_per_face,
                           std::vector<double> (n_edge_points));
                      std::vector<std::vector<Point<dim> > >
                        edge_quadrature_points_full_dim
                          (GeometryInfo<dim>::lines_per_face,
                           std::vector<Point<dim> > (n_edge_points));

                      for (unsigned int q_point = 0; q_point < n_edge_points;
                           ++q_point)
                        {
                  	      edge_quadrature_points[0][q_point]
                  	        = 2.0 * edge_quadrature_y_points[q_point] (0);
                  	      edge_quadrature_points[1][q_point]
                  	        = edge_quadrature_points[0][q_point];
                  	      edge_quadrature_points[2][q_point]
                  	        = 2.0 * edge_quadrature_x_points[q_point] (0)
                  	          - 1.0;
                  	      edge_quadrature_points[3][q_point]
                  	        = edge_quadrature_points[2][q_point];
                          edge_quadrature_points_full_dim[0][q_point]
                            = Point<dim>
                              (0.5, edge_quadrature_y_points[q_point] (0),
                               0.0);
                          edge_quadrature_points_full_dim[1][q_point]
                            = Point<dim>
                              (1.0, edge_quadrature_y_points[q_point] (0),
                               0.0);
                          edge_quadrature_points_full_dim[2][q_point]
                            = Point<dim>
                              (edge_quadrature_x_points[q_point] (0), 0.0,
                               0.0);
                          edge_quadrature_points_full_dim[3][q_point]
                            = Point<dim>
                              (edge_quadrature_x_points[q_point] (0), 0.5,
                               0.0);
                          edge_quadrature_weights[0][q_point]
                            = std::sqrt (edge_quadrature_y.weight (q_point));
                          edge_quadrature_weights[1][q_point]
                            = edge_quadrature_weights[0][q_point];
                          edge_quadrature_weights[2][q_point]
                            = std::sqrt (edge_quadrature_x.weight (q_point));
                          edge_quadrature_weights[3][q_point]
                            = edge_quadrature_weights[2][q_point];
                        }

                                  // Set up the system matrix.
                                  // This can be used for all
                                  // edges.
                  	  for (unsigned int q_point = 0;
                  	       q_point < n_edge_points; ++q_point)
                  	    {
                  	 	  const double tmp
                  	 	    = 2.0 * edge_quadrature_y_points[q_point] (0);
                  	 	  const double weight
                  	 	    = std::sqrt (edge_quadrature_y.weight (q_point));

                          for (unsigned int i = 0; i < deg; ++i)
                            assembling_matrix (i, q_point)
                              = weight
                                * lobatto_polynomials_grad[i + 1].value (tmp);
                  	    }

                      assembling_matrix.mTmult (system_matrix,
                                                assembling_matrix);
                      system_matrix_inv.invert (system_matrix);

                      Vector<double> system_rhs (system_matrix.m ());
                      Vector<double> solution (system_rhs.size ());

                      for (unsigned int dof = 0; dof < this->dofs_per_face;
                          ++dof)
                          for (unsigned int line = 0;
                               line < GeometryInfo<dim - 1>::lines_per_cell;
                               ++line)
                            {
                                  // Set up the right hand side.
                              system_rhs = 0;

                              for (unsigned int q_point = 0;
                                   q_point < n_edge_points; ++q_point)
                                {
                                   const double right_hand_side_value
                                     = edge_quadrature_weights[line][q_point]
                                       * (this->shape_value_component
                                          (this->face_to_cell_index (dof, 4),
                                           edge_quadrature_points_full_dim[line][q_point],
                                           1) - interpolation_matrix
                                                (line * source_fe.degree,
                                                 dof));

                                   for (unsigned int i = 0; i < deg; ++i)
                                     system_rhs (i)
                                       += right_hand_side_value
                                          * lobatto_polynomials_grad[i + 1].value
                                            (edge_quadrature_points[line][q_point]);
                                }

                              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                              for (unsigned int i = 0; i < solution.size ();
                                   ++i)
                                if (std::abs (solution (i)) > 1e-14)
                                  interpolation_matrix
                                  (line * source_fe.degree + i + 1, dof)
                                    = solution (i);
                            }

                                  // Now we project the remaining
                                  // part on the face shape
                                  // functions. First on the
                                  // horizontal ones, then on
                                  // the vertical ones.
                      assembling_matrix.reinit (deg * this->degree,
                                                n_face_points);

                      for (unsigned int q_point = 0;
                           q_point < n_face_points; ++q_point)
                        {
                          const Point<dim - 1> quadrature_point
                                (2.0 * face_quadrature_points[q_point] (0)
                                 - 1.0,
                                 2.0 * face_quadrature_points[q_point] (1));
                          const double weight
                            = std::sqrt (face_quadrature.weight (q_point));

                          for (unsigned int i = 0; i <= deg; ++i)
                            {
                              const double tmp
                                = weight * legendre_polynomials[i].value
                                           (quadrature_point (0));

                              for (unsigned int j = 0; j < deg; ++j)
                           	    assembling_matrix (i * deg + j, q_point)
                           	      = tmp * lobatto_polynomials[j + 2].value
                           	              (quadrature_point (1));
                            }
                        }

                      system_matrix.reinit (assembling_matrix.m (),
                                            assembling_matrix.m ());
                      assembling_matrix.mTmult (system_matrix,
                                                assembling_matrix);
                      system_matrix_inv.reinit (system_matrix.m (),
                                                system_matrix.m ());
                      system_matrix_inv.invert (system_matrix);
                      solution.reinit (system_matrix_inv.m ());
                      system_rhs.reinit (assembling_matrix.m ());

                      for (unsigned int dof = 0; dof < this->dofs_per_face;
                          ++dof)
                        {
                          system_rhs = 0;

                          for (unsigned int q_point = 0;
                               q_point < n_face_points; ++q_point)
                            {
                              const Point<dim> quadrature_point
                                (2.0 * face_quadrature_points[q_point] (0)
                                 - 1.0,
                                 2.0 * face_quadrature_points[q_point] (1),
                                 0.0);
                              double right_hand_side_value
                                = this->shape_value_component
                                  (this->face_to_cell_index (dof, 4),
                                   Point<dim>
                                   (face_quadrature_points[q_point] (0),
                                    face_quadrature_points[q_point] (1), 0),
                                   1);

                              for (unsigned int i = 0; i < 2; ++i)
                                for (unsigned int j = 0; j < source_fe.degree;
                                     ++j)
                                  right_hand_side_value
                                    -= interpolation_matrix
                                       (i * source_fe.degree + j, dof)
                                       * source_fe.shape_value_component
                                         (i * source_fe.degree + j,
                                          quadrature_point, 1);

                              right_hand_side_value
                                *= face_quadrature.weight (q_point);

                              for (unsigned int i = 0; i <= deg; ++i)
                                {
                                  const double tmp
                                    = right_hand_side_value
                                      * legendre_polynomials[i].value
                                        (quadrature_point (0));

                                  for (unsigned int j = 0; j < deg; ++j)
                                    system_rhs (i * deg + j)
                                      = tmp
                                        * lobatto_polynomials[j + 2].value
                           	              (quadrature_point (1));
                                }
                            }

                          system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                          for (unsigned int i = 0; i <= deg; ++i)
                            for (unsigned int j = 0; j < deg; ++j)
                              if (std::abs (solution (i * deg + j)) > 1e-14)
                                interpolation_matrix
                                ((i + 4) * source_fe.degree + j - i, dof)
                                  = solution (i * deg + j);

                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.
                          system_rhs = 0;

                          for (unsigned int q_point = 0;
                               q_point < n_face_points; ++q_point)
                            {
                              const Point<dim> quadrature_point
                                (2.0 * face_quadrature_points[q_point] (0)
                                 - 1.0,
                                 2.0 * face_quadrature_points[q_point] (1),
                                 0.0);
                              double right_hand_side_value
                                = this->shape_value_component
                                  (this->face_to_cell_index (dof, 4),
                                   Point<dim>
                                   (face_quadrature_points[q_point] (0),
                                    face_quadrature_points[q_point] (1), 0),
                                   0);

                              for (unsigned int i = 0; i < 2; ++i)
                                for (unsigned int j = 0; j < source_fe.degree;
                                     ++j)
                                  right_hand_side_value
                                    -= interpolation_matrix
                                       ((i + 2) * source_fe.degree + j, dof)
                                       * source_fe.shape_value_component
                                         (i * source_fe.degree + j,
                                          quadrature_point, 0);

                              right_hand_side_value
                                *= face_quadrature.weight (q_point);

                              for (unsigned int i = 0; i <= deg; ++i)
                                {
                                  const double tmp
                                    = right_hand_side_value
                                      * legendre_polynomials[i].value
                                        (quadrature_point (0));

                                  for (unsigned int j = 0; j < deg; ++j)
                                    system_rhs (i * deg + j)
                                      += tmp
                                         * lobatto_polynomials[j + 2].value
                                           (quadrature_point (1));
                                }
                            }

                          system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                          for (unsigned int i = 0; i <= deg; ++i)
                            for (unsigned int j = 0; j < deg; ++j)
                              if (std::abs (solution (i * deg + j)) > 1e-14)
                                interpolation_matrix
                                (i + (j + source_fe.degree + 3)
                                     * source_fe.degree, dof)
                                  = solution (i * deg + j);
                        }
                    }

                  break;
                }

              case 2:
                {
                  const Quadrature<1>& edge_quadrature_x
                    = QProjector<1>::project_to_child
                      (reference_edge_quadrature, 0);
                  const Quadrature<1>& edge_quadrature_y
                    = QProjector<1>::project_to_child
                      (reference_edge_quadrature, 1);
                  const unsigned int& n_edge_points
                    = edge_quadrature_x.size ();
                  const std::vector<Point<1> >&
                    edge_quadrature_x_points = edge_quadrature_x.get_points ();
                  const std::vector<Point<1> >&
                    edge_quadrature_y_points = edge_quadrature_y.get_points ();

                                  // Let us begin with the
                                  // interpolation part.
                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                    {
                      const double weight
                        = edge_quadrature_x.weight (q_point);

                      for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int dof = 0; dof < this->dofs_per_face;
                             ++dof)
                          {
                            interpolation_matrix (i * source_fe.degree, dof)
                              += weight
                                 * this->shape_value_component
                                   (this->face_to_cell_index (dof, 4),
                                    Point<dim>
                                    (0.5 * i,
                                     edge_quadrature_y_points[q_point] (0), 0.0),
                                    1);
                            interpolation_matrix ((i + 2) * source_fe.degree,
                                                  dof)
                              += weight
                                 * this->shape_value_component
                                   (this->face_to_cell_index (dof, 4),
                                    Point<dim>
                                    (edge_quadrature_x_points[q_point] (0),
                                     0.5 * (i + 1), 0.0), 0);
                          }
                    }

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int dof = 0; dof < this->dofs_per_face;
                         ++dof)
                      {
                        if (std::abs (interpolation_matrix
                                      (i * source_fe.degree, dof)) < 1e-14)
                          interpolation_matrix (i * source_fe.degree, dof)
                            = 0.0;

                        if (std::abs (interpolation_matrix
                                      ((i + 2) * source_fe.degree, dof))
                              < 1e-14)
                          interpolation_matrix ((i + 2) * source_fe.degree,
                                                dof) = 0.0;
                      }

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some higher order edge
                                  // shape functions to
                                  // consider.
                                  // Here the projection part
                                  // starts. The dof values
                                  // are obtained by solving
                                  // a linear system of
                                  // equations.
                  if (deg > 0)
                    {
                    			  // We start with projection
                    			  // on the higher order edge
                    			  // shape function.
                      const QGauss<dim - 1> reference_face_quadrature (this->degree);
                      const Quadrature<dim - 1>& face_quadrature
                        = QProjector<dim - 1>::project_to_child
                          (reference_face_quadrature, 2);
                      const std::vector<Point<dim - 1> >&
                        face_quadrature_points = face_quadrature.get_points ();
                      const std::vector<Polynomials::Polynomial<double> >& legendre_polynomials
                        = Polynomials::Legendre::generate_complete_basis (deg);
                      const std::vector<Polynomials::Polynomial<double> >& lobatto_polynomials
                        = Polynomials::Lobatto::generate_complete_basis (this->degree);
                      const unsigned int& n_face_points
                        = face_quadrature.size ();
                      FullMatrix<double> assembling_matrix (deg,
                                                            n_edge_points);
                      FullMatrix<double> system_matrix (deg, deg);
                      FullMatrix<double> system_matrix_inv (deg, deg);
                      std::vector<Polynomials::Polynomial<double> >
                        lobatto_polynomials_grad (this->degree);

                      for (unsigned int i = 0;
                           i < lobatto_polynomials_grad.size (); ++i)
                        lobatto_polynomials_grad[i]
                          = lobatto_polynomials[i + 1].derivative ();

                                  // Shifted and scaled
                                  // quadrature points and
                                  // weights on the four
                                  // edges of a face.
                      std::vector<std::vector<double> >
                        edge_quadrature_points
                        (GeometryInfo<dim>::lines_per_face,
                         std::vector<double> (n_edge_points));
                      std::vector<std::vector<double> >
                        edge_quadrature_weights
                        (GeometryInfo<dim>::lines_per_face,
                         std::vector<double> (n_edge_points));
                      std::vector<std::vector<Point<dim> > >
                        edge_quadrature_points_full_dim
                        (GeometryInfo<dim>::lines_per_face,
                         std::vector<Point<dim> > (n_edge_points));

                      for (unsigned int q_point = 0; q_point < n_edge_points;
                           ++q_point)
                        {
                  	      edge_quadrature_points[0][q_point]
                  	        = 2.0 * edge_quadrature_y_points[q_point] (0)
                  	          - 1.0;
                  	      edge_quadrature_points[1][q_point]
                  	        = edge_quadrature_points[0][q_point];
                  	      edge_quadrature_points[2][q_point]
                  	        = 2.0 * edge_quadrature_x_points[q_point] (0);
                  	      edge_quadrature_points[3][q_point]
                  	        = edge_quadrature_points[2][q_point];
                          edge_quadrature_points_full_dim[0][q_point]
                            = Point<dim>
                              (0.0, edge_quadrature_y_points[q_point] (0),
                               0.0);
                          edge_quadrature_points_full_dim[1][q_point]
                            = Point<dim>
                              (0.5, edge_quadrature_y_points[q_point] (0),
                               0.0);
                          edge_quadrature_points_full_dim[2][q_point]
                            = Point<dim>
                              (edge_quadrature_x_points[q_point] (0), 0.5,
                               0.0);
                          edge_quadrature_points_full_dim[3][q_point]
                            = Point<dim>
                              (edge_quadrature_x_points[q_point] (0), 1.0,
                               0.0);
                          edge_quadrature_weights[0][q_point]
                            = std::sqrt (edge_quadrature_y.weight (q_point));
                          edge_quadrature_weights[1][q_point]
                            = edge_quadrature_weights[0][q_point];
                          edge_quadrature_weights[2][q_point]
                            = std::sqrt (edge_quadrature_x.weight (q_point));
                          edge_quadrature_weights[3][q_point]
                            = edge_quadrature_weights[2][q_point];
                        }

                                  // Set up the system matrix.
                                  // This can be used for all
                                  // edges.
                  	  for (unsigned int q_point = 0;
                  	       q_point < n_edge_points; ++q_point)
                  	    {
                  	 	  const double weight
                  	        = std::sqrt (edge_quadrature_y.weight
                  	                     (q_point));
               	   	      const double tmp
               	   	        = 2.0 * edge_quadrature_y_points[q_point] (0)
               	   	          - 1.0;

                          for (unsigned int i = 0; i < deg; ++i)
                            assembling_matrix (i, q_point)
                              = weight
                                * lobatto_polynomials_grad[i + 1].value
                                  (tmp);
                        }

                      assembling_matrix.mTmult (system_matrix,
                                               assembling_matrix);
                      system_matrix_inv.invert (system_matrix);

                      Vector<double> system_rhs (system_matrix.m ());
                      Vector<double> solution (system_rhs.size ());

                      for (unsigned int dof = 0; dof < this->dofs_per_face;
                           ++dof)
                          for (unsigned int line = 0;
                               line < GeometryInfo<dim - 1>::lines_per_cell;
                               ++line)
                            {
                                  // Set up the right hand side.
                              system_rhs = 0;

                              for (unsigned int q_point = 0;
                                   q_point < n_edge_points; ++q_point)
                                {
                                  const double right_hand_side_value
                                    = edge_quadrature_weights[line][q_point]
                                      * (this->shape_value_component
                                         (this->face_to_cell_index (dof, 4),
                                          edge_quadrature_points_full_dim[line][q_point],
                                          1) - interpolation_matrix
                                               (line * source_fe.degree, dof));

                                  for (unsigned int i = 0; i < deg; ++i)
                                    system_rhs (i)
                                      += right_hand_side_value
                                         * lobatto_polynomials_grad[i + 1].value
                                           (edge_quadrature_points[line][q_point]);
                                }

                              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                              for (unsigned int i = 0; i < solution.size ();
                                   ++i)
                                if (std::abs (solution (i)) > 1e-14)
                                  interpolation_matrix
                                  (line * source_fe.degree + i + 1, dof)
                                    = solution (i);
                            }

                      assembling_matrix.reinit (deg * this->degree,
                                                n_face_points);

                      for (unsigned int q_point = 0;
                           q_point < n_face_points; ++q_point)
                        {
                           const double weight
                             = std::sqrt (face_quadrature.weight (q_point));
                           const Point<dim - 1> quadrature_point
                             (2.0 * face_quadrature_points[q_point] (0),
                              2.0 * face_quadrature_points[q_point] (1) - 1.0);

                           for (unsigned int i = 0; i <= deg; ++i)
                             {
                               const double tmp
                                 = weight * legendre_polynomials[i].value
                                            (quadrature_point (0));

                               for (unsigned int j = 0; j < deg; ++j)
                           	     assembling_matrix (i * deg + j, q_point)
                           	       = tmp * lobatto_polynomials[j + 2].value
                           	               (quadrature_point (1));
                             }
                        }

                      system_matrix.reinit (assembling_matrix.m (),
                                            assembling_matrix.m ());
                      assembling_matrix.mTmult (system_matrix,
                                                assembling_matrix);
                      system_matrix_inv.reinit (system_matrix.m (),
                                                system_matrix.m ());
                      system_matrix_inv.invert (system_matrix);
                      solution.reinit (system_matrix_inv.m ());
                      system_rhs.reinit (assembling_matrix.m ());
                      system_rhs = 0;

                                  // Now we project the remaining
                                  // part on the face shape
                                  // functions. First on the
                                  // horizontal ones, then on
                                  // the vertical ones.
                      for (unsigned int dof = 0; dof < this->dofs_per_face;
                           ++dof)
                        {
                          for (unsigned int q_point = 0;
                               q_point < n_face_points; ++q_point)
                            {
                              const Point<dim> quadrature_point
                                (2.0 * face_quadrature_points[q_point] (0),
                                 2.0 * face_quadrature_points[q_point] (1)
                                 - 1.0, 0.0);
                              double right_hand_side_value
                                = this->shape_value_component
                                  (this->face_to_cell_index (dof, 4),
                                   Point<dim>
                                   (face_quadrature_points[q_point] (0),
                                    face_quadrature_points[q_point] (1), 0.0),
                                   1);

                              for (unsigned int i = 0; i < 2; ++i)
                                for (unsigned int j = 0; j < source_fe.degree;
                                     ++j)
                                  right_hand_side_value
                                    -= interpolation_matrix
                                       (i * source_fe.degree + j, dof)
                                       * source_fe.shape_value_component
                                         (i * source_fe.degree + j,
                                          quadrature_point, 1);

                              right_hand_side_value
                                *= face_quadrature.weight (q_point);

                              for (unsigned int i = 0; i <= deg; ++i)
                                {
                                  const double tmp
                                    = right_hand_side_value
                                      * legendre_polynomials[i].value
                                        (quadrature_point (0));

                                  for (unsigned int j = 0; j < deg; ++j)
                                      system_rhs (i * deg + j)
                                        +=  tmp
                                            * lobatto_polynomials[j + 2].value
                           	                  (quadrature_point (1));
                                }
                            }

                          system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                          for (unsigned int i = 0; i <= deg; ++i)
                            for (unsigned int j = 0; j < deg; ++j)
                              if (std::abs (solution (i * deg + j)) > 1e-14)
                                interpolation_matrix
                                ((i + 4) * source_fe.degree + j - i, dof)
                                  = solution (i * deg + j);

                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.
                          system_rhs = 0;

                          for (unsigned int q_point = 0;
                               q_point < n_face_points; ++q_point)
                            {
                              const Point<dim> quadrature_point
                                (2.0 * face_quadrature_points[q_point] (0),
                                 2.0 * face_quadrature_points[q_point] (1)
                                 - 1.0, 0.0);
                              double right_hand_side_value
                                = this->shape_value_component
                                  (this->face_to_cell_index (dof, 4),
                                   Point<dim>
                                   (face_quadrature_points[q_point] (0),
                                    face_quadrature_points[q_point] (1), 0.0),
                                   0);

                              for (unsigned int i = 0; i < 2; ++i)
                                for (unsigned int j = 0; j < source_fe.degree;
                                     ++j)
                                  right_hand_side_value
                                    -= interpolation_matrix
                                       ((i + 2) * source_fe.degree + j, dof)
                                       * source_fe.shape_value_component
                                         (i * source_fe.degree + j,
                                          quadrature_point, 0);

                              right_hand_side_value *= face_quadrature.weight
                                                       (q_point);

                              for (unsigned int i = 0; i <= deg; ++i)
                                {
                                  const double tmp
                                    = right_hand_side_value
                                       * legendre_polynomials[i].value
                                         (quadrature_point (0));

                                  for (unsigned int j = 0; j < deg; ++j)
                                    system_rhs (i * deg + j)
                                      += tmp
                                         * lobatto_polynomials[j + 2].value
                                           (quadrature_point (1));
                                }
                            }

                          system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                          for (unsigned int i = 0; i <= deg; ++i)
                            for (unsigned int j = 0; j < deg; ++j)
                              if (std::abs (solution (i * deg + j)) > 1e-14)
                                interpolation_matrix
                                (i + (j + source_fe.degree + 3)
                                 * source_fe.degree, dof)
                                  = solution (i * deg + j);
                        }
                    }

                  break;
                }

              case 3:
                {
                  const Quadrature<1>& edge_quadrature
                    = QProjector<1>::project_to_child
                      (reference_edge_quadrature, 1);
                  const unsigned int& n_edge_points = edge_quadrature.size ();
                  const std::vector<Point<1> >&
                    edge_quadrature_points = edge_quadrature.get_points ();

                                  // Let us begin with the
                                  // interpolation part.
                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                    {
                      const double weight = edge_quadrature.weight (q_point);

                      for (unsigned int i = 0; i < 2; ++i)
                        for (unsigned int dof = 0; dof < this->dofs_per_face;
                             ++dof)
                          {
                            interpolation_matrix (i * source_fe.degree, dof)
                              += weight
                                 * this->shape_value_component
                                   (this->face_to_cell_index (dof, 4),
                                    Point<dim>
                                    (0.5 * (i + 1),
                                     edge_quadrature_points[q_point] (0), 0.0),
                                    1);
                            interpolation_matrix ((i + 2) * source_fe.degree,
                                                  dof)
                              += weight
                                 * this->shape_value_component
                                   (this->face_to_cell_index (dof, 4),
                                    Point<dim>
                                    (edge_quadrature_points[q_point] (0),
                                     0.5 * (i + 1), 0.0), 0);
                          }
                    }

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int dof = 0; dof < this->dofs_per_face;
                         ++dof)
                      {
                        if (std::abs (interpolation_matrix
                                      (i * source_fe.degree, dof)) < 1e-14)
                          interpolation_matrix (i * source_fe.degree, dof)
                            = 0.0;

                        if (std::abs (interpolation_matrix
                                      ((i + 2) * source_fe.degree, dof))
                              < 1e-14)
                          interpolation_matrix ((i + 2) * source_fe.degree,
                                                dof) = 0.0;
                      }

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some higher order edge
                                  // shape functions to
                                  // consider.
                                  // Here the projection part
                                  // starts. The dof values
                                  // are obtained by solving
                                  // a linear system of
                                  // equations.
                  if (deg > 1)
                    {
                    			  // We start with projection
                    			  // on the higher order edge
                    			  // shape function.
                      const QGauss<dim - 1>
                        reference_face_quadrature (this->degree);
                      const Quadrature<dim - 1>& face_quadrature
                        = QProjector<dim - 1>::project_to_child
                          (reference_face_quadrature, 3);
                      const std::vector<Point<dim - 1> >&
                        face_quadrature_points = face_quadrature.get_points ();
                      const std::vector<Polynomials::Polynomial<double> >&
                        legendre_polynomials
                          = Polynomials::Legendre::generate_complete_basis
                            (deg);
                      const std::vector<Polynomials::Polynomial<double> >&
                        lobatto_polynomials
                          = Polynomials::Lobatto::generate_complete_basis
                            (this->degree);
                      const unsigned int& n_face_points
                        = face_quadrature.size ();
                      FullMatrix<double> assembling_matrix (deg,
                                                            n_edge_points);
                      FullMatrix<double> system_matrix (deg, deg);
                      FullMatrix<double> system_matrix_inv (deg, deg);
                      std::vector<Polynomials::Polynomial<double> >
                        lobatto_polynomials_grad (this->degree);

                      for (unsigned int i = 0;
                           i < lobatto_polynomials_grad.size (); ++i)
                        lobatto_polynomials_grad[i]
                          = lobatto_polynomials[i + 1].derivative ();

                                  // Shifted and scaled
                                  // quadrature points on
                                  // the four edges of a
                                  // face.
                      std::vector<std::vector<Point<dim> > >
                        edge_quadrature_points_full_dim
                        (GeometryInfo<dim>::lines_per_face,
                         std::vector<Point<dim> > (n_edge_points));

                      for (unsigned int q_point = 0; q_point < n_edge_points;
                           ++q_point)
                        {
                          edge_quadrature_points_full_dim[0][q_point]
                            = Point<dim>
                              (0.5, edge_quadrature_points[q_point] (0), 0.0);
                          edge_quadrature_points_full_dim[1][q_point]
                            = Point<dim>
                              (1.0, edge_quadrature_points[q_point] (0), 0.0);
                          edge_quadrature_points_full_dim[2][q_point]
                            = Point<dim> (edge_quadrature_points[q_point] (0),
                                          0.5, 0.0);
                          edge_quadrature_points_full_dim[3][q_point]
                            = Point<dim> (edge_quadrature_points[q_point] (0),
                                          1.0, 0.0);
                        }

                                  // Set up the system matrix.
                                  // This can be used for all
                                  // edges.
                  	  for (unsigned int q_point = 0;
                           q_point < n_edge_points; ++q_point)
              	        {
               	 	      const double tmp
               	 	        = 2.0 * edge_quadrature_points[q_point] (0)
               	 	          - 1.0;
               	 	      const double weight
                 	 	    = std::sqrt (edge_quadrature.weight (q_point));

                          for (unsigned int i = 0; i < deg; ++i)
                            assembling_matrix (i, q_point)
                              = weight
                                * lobatto_polynomials_grad[i + 1].value
                                  (tmp);
               	        }

                      assembling_matrix.mTmult (system_matrix,
                                                assembling_matrix);
                      system_matrix_inv.invert (system_matrix);

                      Vector<double> system_rhs (system_matrix.m ());
                      Vector<double> solution (system_rhs.size ());

                      for (unsigned int dof = 0; dof < this->dofs_per_face;
                           ++dof)
                          for (unsigned int line = 0;
                               line < GeometryInfo<dim - 1>::lines_per_cell;
                               ++line)
                            {
                                  // Set up the right hand side.
                              system_rhs = 0;

                              for (unsigned int q_point = 0;
                                   q_point < n_edge_points; ++q_point)
                                {
                                  const double right_hand_side_value
                                    = std::sqrt (edge_quadrature.weight
                                                 (q_point))
                                      * (this->shape_value_component
                                         (this->face_to_cell_index (dof, 4),
                                          edge_quadrature_points_full_dim[line][q_point],
                                          1) - interpolation_matrix
                                               (line * source_fe.degree, dof));
                                  const double tmp
                                    = 2.0 * edge_quadrature_points[q_point] (0)
                                      - 1.0;

                                  for (unsigned int i = 0; i < deg; ++i)
                                    system_rhs (i)
                                      += right_hand_side_value
                                         * lobatto_polynomials_grad[i + 1].value
                                           (tmp);
                                }

                              system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                              for (unsigned int i = 0; i < solution.size ();
                                   ++i)
                                if (std::abs (solution (i)) > 1e-14)
                                  interpolation_matrix
                                  (line * source_fe.degree + i + 1, dof)
                                    = solution (i);
                            }

                      assembling_matrix.reinit (deg * this->degree,
                                                n_face_points);

                      for (unsigned int q_point = 0;
                           q_point < n_face_points; ++q_point)
                        {
                          const double weight
                            = std::sqrt (face_quadrature.weight
                                         (q_point));
                          const Point<dim - 1> quadrature_point
                             (2.0 * face_quadrature_points[q_point] (0) - 1.0,
                              2.0 * face_quadrature_points[q_point] (1) - 1.0);

                          for (unsigned int i = 0; i <= deg; ++i)
                            {
                              const double tmp
                                = weight * legendre_polynomials[i].value
                                           (quadrature_point (0));

                              for (unsigned int j = 0; j < deg; ++j)
                           	    assembling_matrix (i * deg + j, q_point)
                           	      = tmp * lobatto_polynomials[j + 2].value
                           	              (quadrature_point (1));
                            }
                        }

                      system_matrix.reinit (assembling_matrix.m (),
                                            assembling_matrix.m ());
                      assembling_matrix.mTmult (system_matrix,
                                                assembling_matrix);
                      system_matrix_inv.reinit (system_matrix.m (),
                                                system_matrix.m ());
                      system_matrix_inv.invert (system_matrix);
                      solution.reinit (system_matrix.m ());
                      system_rhs.reinit (assembling_matrix.m ());
                      system_rhs = 0;

                      for (unsigned int dof = 0; dof < this->dofs_per_face;
                           ++dof)
                        {
                                  // Now we project the remaining
                                  // part on the face shape
                                  // functions. First on the
                                  // horizontal ones, then on
                                  // the vertical ones.
                          for (unsigned int q_point = 0;
                               q_point < n_face_points; ++q_point)
                            {
                              const Point<dim> quadrature_point
                                (2.0 * face_quadrature_points[q_point] (0)
                                 - 1.0,
                                 2.0 * face_quadrature_points[q_point] (1)
                                 - 1.0, 0.0);
                              double right_hand_side_value
                                = this->shape_value_component
                                  (this->face_to_cell_index (dof, 4),
                                   Point<dim>
                                   (face_quadrature_points[q_point] (0),
                                    face_quadrature_points[q_point] (1), 0.0),
                                   1);

                              for (unsigned int i = 0; i < 2; ++i)
                                for (unsigned int j = 0; j < source_fe.degree;
                                     ++j)
                                  right_hand_side_value
                                    -= interpolation_matrix
                                       (i * source_fe.degree + j, dof)
                                       * source_fe.shape_value_component
                                         (i * source_fe.degree + j,
                                          quadrature_point, 1);

                              right_hand_side_value
                                *= face_quadrature.weight (q_point);

                              for (unsigned int i = 0; i <= deg; ++i)
                                {
                                  const double tmp
                                    = right_hand_side_value
                                       * legendre_polynomials[i].value
                                         (quadrature_point (0));

                                  for (unsigned int j = 0; j < deg; ++j)
                                    system_rhs (i * deg + j)
                                      +=  tmp
                                          * lobatto_polynomials[j + 2].value
                           	                (quadrature_point (1));
                                }
                            }

                          system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                          for (unsigned int i = 0; i <= deg; ++i)
                            for (unsigned int j = 0; j < deg; ++j)
                              if (std::abs (solution (i * deg + j)) > 1e-14)
                                interpolation_matrix
                                ((i + 4) * source_fe.degree + j - i, dof)
                                  = solution (i * deg + j);

                                  // Set up the right hand side
                                  // for the vertical shape
                                  // functions.
                          system_rhs = 0;

                          for (unsigned int q_point = 0;
                               q_point < n_face_points; ++q_point)
                            {
                              const Point<dim> quadrature_point
                                (2.0 * face_quadrature_points[q_point] (0)
                                 - 1.0,
                                 2.0 * face_quadrature_points[q_point] (1)
                                 - 1.0, 0.0);
                              double right_hand_side_value
                                = this->shape_value_component
                                  (this->face_to_cell_index (dof, 4),
                                   Point<dim>
                                   (face_quadrature_points[q_point] (0),
                                    face_quadrature_points[q_point] (1), 0.0),
                                   0);

                              for (unsigned int i = 0; i < 2; ++i)
                                for (unsigned int j = 0; j < source_fe.degree;
                                     ++j)
                                  right_hand_side_value
                                    -= interpolation_matrix
                                       ((i + 2) * source_fe.degree + j, dof)
                                       * source_fe.shape_value_component
                                         (i * source_fe.degree + j,
                                          quadrature_point, 0);

                              right_hand_side_value
                                *= face_quadrature.weight (q_point);

                              for (unsigned int i = 0; i <= deg; ++i)
                                {
                                  const double L_i
                                    = legendre_polynomials[i].value
                                      (quadrature_point (0));
                                  const double tmp
                                    = right_hand_side_value * L_i;

                                  for (unsigned int j = 0; j < deg; ++j)
                                    system_rhs (i * deg + j)
                                      += tmp
                                         * lobatto_polynomials[j + 2].value
                                           (quadrature_point (1));
                                }
                            }

                          system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the interpolation
                                  // matrix only, if they are
                                  // not too small.
                          for (unsigned int i = 0; i <= deg; ++i)
                            for (unsigned int j = 0; j < deg; ++j)
                              if (std::abs (solution (i * deg + j)) > 1e-14)
                                interpolation_matrix
                                (i + (j + source_fe.degree + 3)
                                 * source_fe.degree, dof)
                                  = solution (i * deg + j);
                        }
                    }

                  break;
                }

              default:
                Assert (false, ExcNotImplemented ());
            }

          break;
        }

      default:
        Assert (false, ExcNotImplemented ());
    }
}

#endif

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
                   // This is done as usual by projection-based
                   // interpolation.
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
    switch (dim)
      {
        case 2:
          {
            const QGauss<1> reference_edge_quadrature (this->degree);
            const unsigned int& n_edge_points
              = reference_edge_quadrature.size ();

                                  // Let us begin with the
                                  // interpolation part.
            for (unsigned int i = 0; i < 2; ++i)
              {
                for (unsigned int q_point = 0; q_point < n_edge_points;
                     ++q_point)
                  local_dofs[i * this->degree]
                    += reference_edge_quadrature.weight (q_point)
                       * values[q_point + i * n_edge_points] (1);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                if (std::abs (local_dofs[i * this->degree]) < 1e-14)
                  local_dofs[i * this->degree] = 0.0;
              }

            if (offset == 0)
              for (unsigned int i = 0; i < 2; ++i)
                {
                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                    local_dofs[(i + 2) * this->degree]
                      += reference_edge_quadrature.weight (q_point)
                         * values[q_point + (i + 2) * n_edge_points] (0);

                  if (std::abs (local_dofs[(i + 2) * this->degree]) < 1e-14)
                    local_dofs[(i + 2) * this->degree] = 0.0;
                }

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some higher order edge
                                  // shape functions to
                                  // consider.
                                  // Here the projection part
                                  // starts. The dof values
                                  // are obtained by solving
                                  // a linear system of
                                  // equations.
            if (deg > 0)
              {
                    			  // We start with projection
                    			  // on the higher order edge
                    			  // shape function.
                const std::vector<Polynomials::Polynomial<double> >&
                  lobatto_polynomials
                    = Polynomials::Lobatto::generate_complete_basis
                      (this->degree);
                const unsigned int
                  line_coordinate[GeometryInfo<2>::lines_per_cell]
                    = {1, 1, 0, 0};
                std::vector<Polynomials::Polynomial<double> >
                  lobatto_polynomials_grad (this->degree);

                for (unsigned int i = 0; i < lobatto_polynomials_grad.size ();
                     ++i)
                  lobatto_polynomials_grad[i]
                    = lobatto_polynomials[i + 1].derivative ();

                                  // Set up the system matrix.
                                  // This can be used for all
                                  // edges.
                FullMatrix<double> system_matrix (deg, deg);

                for (unsigned int i = 0; i < system_matrix.m (); ++i)
                  for (unsigned int j = 0; j < system_matrix.n (); ++j)
                    for (unsigned int q_point = 0; q_point < n_edge_points;
                         ++q_point)
                      system_matrix (i, j)
                        += boundary_weights (q_point, j)
                           * lobatto_polynomials_grad[i + 1].value
                             (this->generalized_face_support_points[q_point]
                              (1));

                FullMatrix<double> system_matrix_inv (deg, deg);

                system_matrix_inv.invert (system_matrix);

                Vector<double> system_rhs (system_matrix.m ());
                Vector<double> solution (system_rhs.size ());

                for (unsigned int line = 0;
                     line < GeometryInfo<dim>::lines_per_cell; ++line)
                  if ((line < 2) || (offset == 0))
                    {
                                  // Set up the right hand side.
                      system_rhs = 0;

                      for (unsigned int q_point = 0; q_point < n_edge_points;
                           ++q_point)
                        {
                          const double tmp
                            = values[line * n_edge_points + q_point]
                              (line_coordinate[line])
                              - local_dofs[line * this->degree]
                              * this->shape_value_component
                                (line * this->degree,
                                 this->generalized_support_points[line
                                                                  * n_edge_points
                                                                  + q_point],
                                 line_coordinate[line]);

                          for (unsigned int i = 0; i < system_rhs.size ();
                               ++i)
                            system_rhs (i) += boundary_weights (q_point, i)
                                              * tmp;
                        }

                      system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                      for (unsigned int i = 0; i < solution.size (); ++i)
                        if (std::abs (solution (i)) > 1e-14)
                          local_dofs[line * this->degree + i + 1]
                            = solution (i);
                    }

                                  // Then we go on to the
                                  // interior shape
                                  // functions. Again we
                                  // set up the system
                                  // matrix and use it
                                  // for both, the
                                  // horizontal and the
                                  // vertical, interior
                                  // shape functions.
                const QGauss<dim> reference_quadrature (this->degree);
                const std::vector<Polynomials::Polynomial<double> >&
                  legendre_polynomials
                    = Polynomials::Legendre::generate_complete_basis (deg);
                const unsigned int& n_interior_points
                  = reference_quadrature.size ();

                system_matrix.reinit (deg * this->degree, deg * this->degree);
                system_matrix = 0;

                for (unsigned int i = 0; i <= deg; ++i)
                  for (unsigned int j = 0; j < deg; ++j)
                    for (unsigned int k = 0; k <= deg; ++k)
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

                system_matrix_inv.reinit (system_matrix.m (),
                                          system_matrix.m ());
                system_matrix_inv.invert (system_matrix);
                solution.reinit (system_matrix_inv.m ());
                system_rhs.reinit (system_matrix.m ());

                if (offset == 0)
                  {
                                  // Set up the right hand side
                                  // for the horizontal shape
                                  // functions.
                    system_rhs = 0;

                    for (unsigned int q_point = 0;
                         q_point < n_interior_points; ++q_point)
                      {
                        double tmp
                          = values[q_point + GeometryInfo<dim>::lines_per_cell
                                   * n_edge_points] (0);

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

                    system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                    for (unsigned int i = 0; i <= deg; ++i)
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
            const QGauss<1>
              reference_edge_quadrature (this->degree);
            const unsigned int&
              n_edge_points = reference_edge_quadrature.size ();

                                  // Let us begin with the
                                  // interpolation part.
            for (unsigned int i = 0; i < 4; ++i)
              {
                for (unsigned int q_point = 0; q_point < n_edge_points;
                     ++q_point)
                  local_dofs[(i + 8) * this->degree]
                    += reference_edge_quadrature.weight (q_point)
                       * values[q_point + (i + 8) * n_edge_points] (2);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                if (std::abs (local_dofs[(i + 8) * this->degree]) < 1e-14)
                  local_dofs[(i + 8) * this->degree] = 0.0;
              }

            if (offset + 1 < dim)
              {
                for (unsigned int i = 0; i < 2; ++i)
                  for (unsigned int j = 0; j < 2; ++j)
                    {
                      for (unsigned int q_point = 0; q_point < n_edge_points;
                           ++q_point)
                        local_dofs[(i + 4 * j) * this->degree]
                          += reference_edge_quadrature.weight (q_point)
                             * values[q_point + (i + 4 * j) * n_edge_points]
                               (1);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                      if (std::abs (local_dofs[(i + 4 * j) * this->degree])
                            < 1e-14)
                        local_dofs[(i + 4 * j) * this->degree] = 0.0;
                    }

                if (offset == 0)
                  for (unsigned int i = 0; i < 2; ++i)
                    for (unsigned int j = 0; j < 2; ++j)
                      {
                        for (unsigned int q_point = 0;
                             q_point < n_edge_points; ++q_point)
                          local_dofs[(i + 4 * j + 2) * this->degree]
                            += reference_edge_quadrature.weight (q_point)
                               * values[q_point + (i + 4 * j + 2)
                                        * n_edge_points] (0);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                        if (std::abs (local_dofs[(i + 4 * j + 2)
                                                 * this->degree]) < 1e-14)
                          local_dofs[(i + 4 * j + 2) * this->degree] = 0.0;
                      }
              }

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some higher order shape
                                  // functions to consider.
                                  // Here the projection part
                                  // starts. The dof values
                                  // are obtained by solving
                                  // a linear system of
                                  // equations.
            if (deg > 0)
              {
                    			  // We start with projection
                    			  // on the higher order edge
                    			  // shape function.
                const std::vector<Polynomials::Polynomial<double> >&
                  lobatto_polynomials
                    = Polynomials::Lobatto::generate_complete_basis
                      (this->degree);
                const unsigned int
                  line_coordinate[GeometryInfo<3>::lines_per_cell]
                    = {1, 1, 0, 0, 1, 1, 0, 0, 2, 2, 2, 2};
                FullMatrix<double> system_matrix (deg, deg);
                FullMatrix<double> system_matrix_inv (deg, deg);
                std::vector<Polynomials::Polynomial<double> >
                  lobatto_polynomials_grad (this->degree);

                for (unsigned int i = 0; i < lobatto_polynomials_grad.size ();
                     ++i)
                  lobatto_polynomials_grad[i]
                    = lobatto_polynomials[i + 1].derivative ();

                Vector<double> system_rhs (system_matrix.m ());
                Vector<double> solution (system_rhs.size ());

                                  // Set up the system matrix.
                                  // This can be used for all
                                  // edges.
                for (unsigned int i = 0; i < system_matrix.m (); ++i)
                  for (unsigned int j = 0; j < system_matrix.n (); ++j)
                    for (unsigned int q_point = 0; q_point < n_edge_points;
                         ++q_point)
                      system_matrix (i, j)
                        += boundary_weights (q_point, j)
                           * lobatto_polynomials_grad[i + 1].value
                             (this->generalized_face_support_points[q_point]
                              (1));

                system_matrix_inv.invert (system_matrix);

                for (unsigned int line = 0;
                     line < GeometryInfo<dim>::lines_per_cell; ++line)
                  {
                                  // Set up the right hand side.
                    system_rhs = 0;

                    if ((((line == 0) || (line == 1) || (line == 4) ||
                          (line == 5)) && (offset + 1 < dim)) ||
                        (((line == 2) || (line == 3) || (line == 6) ||
                          (line == 7)) && (offset == 0)) || (line > 7))
                      {
                        for (unsigned int q_point = 0; q_point < n_edge_points;
                             ++q_point)
                          {
                            double tmp
                              = values[line * n_edge_points + q_point]
                                (line_coordinate[line])
                                - local_dofs[line * this->degree]
                                * this->shape_value_component
                                  (line * this->degree,
                                   this->generalized_support_points[line
                                                                    * this->degree
                                                                    + q_point],
                                   line_coordinate[line]);

                            for (unsigned int i = 0; i < system_rhs.size ();
                                 ++i)
                              system_rhs (i)
                                += boundary_weights (q_point, i) * tmp;
                          }

                        system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                        for (unsigned int i = 0; i < solution.size (); ++i)
                          if (std::abs (solution (i)) > 1e-14)
                            local_dofs[line * this->degree + i + 1]
                              = solution (i);
                      }
                  }

                                  // Then we go on to the
                                  // face shape functions.
                                  // Again we set up the
                                  // system matrix and
                                  // use it for both, the
                                  // horizontal and the
                                  // vertical, shape
                                  // functions.
                const std::vector<Polynomials::Polynomial<double> >&
                  legendre_polynomials
                    = Polynomials::Legendre::generate_complete_basis (deg);
                const unsigned int
                  n_face_points = n_edge_points * n_edge_points;

                system_matrix.reinit (deg * this->degree, deg * this->degree);
                system_matrix = 0;

                for (unsigned int i = 0; i <= deg; ++i)
                  for (unsigned int j = 0; j < deg; ++j)
                    for (unsigned int k = 0; k <= deg; ++k)
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

                system_matrix_inv.reinit (system_matrix.m (),
                                          system_matrix.n ());
                system_matrix_inv.invert (system_matrix);
                solution.reinit (system_matrix.m ());
                system_rhs.reinit (system_matrix.m ());

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
                                      for (unsigned int j = 0; j <= deg; ++j)
                                        tmp
                                          -= local_dofs[4 * i * this->degree
                                                        + j]
                                             * this->shape_value_component
                                               (4 * i * this->degree + j,
                                                this->generalized_support_points[q_point
                                                                                 + GeometryInfo<dim>::lines_per_cell
                                                                                 * n_edge_points],
                                                1);

                                    for (unsigned int i = 0; i <= deg; ++i)
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
                                for (unsigned int i = 0; i <= deg; ++i)
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
                                  for (unsigned int j = 0; j <= deg; ++j)
                                    tmp -= local_dofs[2 * (i + 4)
                                                      * this->degree + j]
                                           * this->shape_value_component
                                             (2 * (i + 4) * this->degree + j,
                                              this->generalized_support_points[q_point
                                                                               + GeometryInfo<dim>::lines_per_cell
                                                                               * n_edge_points],
                                              2);

                                for (unsigned i = 0; i <= deg; ++i)
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
                            for (unsigned int i = 0; i <= deg; ++i)
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
                                      for (unsigned int j = 0; j <= deg; ++j)
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

                                    for (unsigned int i = 0; i <= deg; ++i)
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
                                for (unsigned int i = 0; i <= deg; ++i)
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
                                  for (unsigned int j = 0; j <= deg; ++j)
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

                                for (unsigned i = 0; i <= deg; ++i)
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
                            for (unsigned int i = 0; i <= deg; ++i)
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
                                     q_point < n_face_points;
                                     ++q_point)
                                  {
                                    double tmp
                                      = values[q_point
                                               + GeometryInfo<dim>::lines_per_cell
                                               * n_edge_points
                                               + 2 * n_face_points] (0);

                                    for (unsigned int i = 0; i < 2; ++i)
                                      for (unsigned int j = 0; j <= deg; ++j)
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

                                      for (unsigned int i = 0; i <= deg; ++i)
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
                                for (unsigned int i = 0; i <= deg; ++i)
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
                                           * n_edge_points + 2 * n_face_points]
                                    (2);

                                for (unsigned int i = 0; i < 2; ++i)
                                  for (unsigned int j = 0; j <= deg; ++j)
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

                                for (unsigned i = 0; i <= deg; ++i)
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
                            for (unsigned int i = 0; i <= deg; ++i)
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
                                               * n_edge_points + 3
                                               * n_face_points] (0);

                                    for (unsigned int i = 0; i < 2; ++i)
                                      for (unsigned int j = 0; j <= deg; ++j)
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

                                    for (unsigned int i = 0; i <= deg; ++i)
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
                                for (unsigned int i = 0; i <= deg; ++i)
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
                                           * n_edge_points + 3 * n_face_points]
                                    (2);

                                for (unsigned int i = 0; i < 2; ++i)
                                  for (unsigned int j = 0; j <= deg; ++j)
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

                                for (unsigned i = 0; i <= deg; ++i)
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
                            for (unsigned int i = 0; i <= deg; ++i)
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
                                          for (unsigned int j = 0; j <= deg; ++j)
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

                                        for (unsigned int i = 0; i <= deg; ++i)
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
                                    for (unsigned int i = 0; i <= deg; ++i)
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
                                      for (unsigned int j = 0; j <= deg; ++j)
                                        tmp -= local_dofs[i * this->degree + j]
                                               * this->shape_value_component
                                                 (i * this->degree + j,
                                                  this->generalized_support_points[q_point
                                                                                   + GeometryInfo<dim>::lines_per_cell
                                                                                   * n_edge_points
                                                                                   + 4
                                                                                   * n_face_points],
                                                  1);

                                    for (unsigned i = 0; i <= deg; ++i)
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
                                for (unsigned int i = 0; i <= deg; ++i)
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
                                        for (unsigned int j = 0; j <= deg; ++j)
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

                                      for (unsigned int i = 0; i <= deg; ++i)
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
                                  for (unsigned int i = 0; i <= deg; ++i)
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
                                    for (unsigned int j = 0; j <= deg; ++j)
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

                                  for (unsigned i = 0; i <= deg; ++i)
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
                              for (unsigned int i = 0; i <= deg; ++i)
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

                            for (unsigned int i = 0; i <= deg; ++i)
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

                            for (unsigned int i = 0; i <= deg; ++i)
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
                        for (unsigned int i = 0; i <= deg; ++i)
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

                        for (unsigned int i = 0; i <= deg; ++i)
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

  switch (dim)
    {
      case 2:
        {
                                  // Let us begin with the
                                  // interpolation part.
          const QGauss<dim - 1> reference_edge_quadrature (this->degree);
          const unsigned int&
            n_edge_points = reference_edge_quadrature.size ();

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              {
                for (unsigned int q_point = 0; q_point < n_edge_points;
                     ++q_point)
                  local_dofs[(i + 2 * j) * this->degree]
                    += reference_edge_quadrature.weight (q_point)
                       * values[1 - j][q_point + (i + 2 * j) * n_edge_points];

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                if (std::abs (local_dofs[(i + 2 * j) * this->degree]) < 1e-14)
                  local_dofs[(i + 2 * j) * this->degree] = 0.0;
              }

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some higher order edge
                                  // shape functions to
                                  // consider.
                                  // Here the projection part
                                  // starts. The dof values
                                  // are obtained by solving
                                  // a linear system of
                                  // equations.
          if (deg > 0)
            {
                    			  // We start with projection
                    			  // on the higher order edge
                    			  // shape function.
              const std::vector<Polynomials::Polynomial<double> >&
                lobatto_polynomials
                  = Polynomials::Lobatto::generate_complete_basis
                    (this->degree);
              FullMatrix<double> system_matrix (deg, deg);
              std::vector<Polynomials::Polynomial<double> >
                lobatto_polynomials_grad (this->degree);

              for (unsigned int i = 0; i < lobatto_polynomials_grad.size ();
                   ++i)
                lobatto_polynomials_grad[i]
                  = lobatto_polynomials[i + 1].derivative ();

                                  // Set up the system matrix.
                                  // This can be used for all
                                  // edges.
              for (unsigned int i = 0; i < system_matrix.m (); ++i)
                for (unsigned int j = 0; j < system_matrix.n (); ++j)
                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                     system_matrix (i, j)
                       += boundary_weights (q_point, j)
                          * lobatto_polynomials_grad[i + 1].value
                            (this->generalized_face_support_points[q_point]
                             (1));

              FullMatrix<double> system_matrix_inv (deg, deg);

              system_matrix_inv.invert (system_matrix);

              const unsigned int
                line_coordinate[GeometryInfo<2>::lines_per_cell]
                  = {1, 1, 0, 0};
              Vector<double> system_rhs (system_matrix.m ());
              Vector<double> solution (system_rhs.size ());

              for (unsigned int line = 0;
                   line < GeometryInfo<dim>::lines_per_cell; ++line)
                {
                                  // Set up the right hand side.
                  system_rhs = 0;

                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                    {
                      const double tmp
                        = values[line_coordinate[line]][line * n_edge_points
                                                        + q_point]
                          - local_dofs[line * this->degree]
                          * this->shape_value_component
                            (line * this->degree,
                             this->generalized_support_points[line
                                                              * n_edge_points
                                                              + q_point],
                             line_coordinate[line]);

                      for (unsigned int i = 0; i < system_rhs.size (); ++i)
                        system_rhs (i) += boundary_weights (q_point, i) * tmp;
                    }

                  system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                  for (unsigned int i = 0; i < solution.size (); ++i)
                    if (std::abs (solution (i)) > 1e-14)
                      local_dofs[line * this->degree + i + 1] = solution (i);
                }

                                  // Then we go on to the
                                  // interior shape
                                  // functions. Again we
                                  // set up the system
                                  // matrix and use it
                                  // for both, the
                                  // horizontal and the
                                  // vertical, interior
                                  // shape functions.
              const QGauss<dim> reference_quadrature (this->degree);
              const unsigned int&
                n_interior_points = reference_quadrature.size ();
              const std::vector<Polynomials::Polynomial<double> >&
                legendre_polynomials
                  = Polynomials::Legendre::generate_complete_basis (deg);

              system_matrix.reinit (deg * this->degree, deg * this->degree);
              system_matrix = 0;

              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k <= deg; ++k)
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

              system_matrix_inv.reinit (system_matrix.m (),
                                        system_matrix.m ());
              system_matrix_inv.invert (system_matrix);
                                  // Set up the right hand side
                                  // for the horizontal shape
                                  // functions.
              system_rhs.reinit (system_matrix_inv.m ());
              system_rhs = 0;

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
                                  // interpolation part.
          const QGauss<1> reference_edge_quadrature (this->degree);
          const unsigned int&
            n_edge_points = reference_edge_quadrature.size ();

          for (unsigned int q_point = 0; q_point < n_edge_points; ++q_point)
            {
              for (unsigned int i = 0; i < 4; ++i)
                local_dofs[(i + 8) * this->degree]
                  += reference_edge_quadrature.weight (q_point)
                     * values[2][q_point + (i + 8) * n_edge_points];

              for (unsigned int i = 0; i < 2; ++i)
                for (unsigned int j = 0; j < 2; ++j)
                  for (unsigned int k = 0; k < 2; ++k)
                    local_dofs[(i + 2 * (2 * j + k)) * this->degree]
                      += reference_edge_quadrature.weight (q_point)
                         * values[1 - k][q_point + (i + 2 * (2 * j + k))
                                         * n_edge_points];
            }

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
          for (unsigned int i = 0; i < 4; ++i)
            if (std::abs (local_dofs[(i + 8) * this->degree]) < 1e-14)
              local_dofs[(i + 8) * this->degree] = 0.0;

          for (unsigned int i = 0; i < 2; ++i)
            for (unsigned int j = 0; j < 2; ++j)
              for (unsigned int k = 0; k < 2; ++k)
                if (std::abs (local_dofs[(i + 2 * (2 * j + k)) * this->degree])
                      < 1e-14)
                  local_dofs[(i + 2 * (2 * j + k)) * this->degree] = 0.0;

                                  // If the degree is greater
                                  // than 0, then we have still
                                  // some higher order shape
                                  // functions to consider.
                                  // Here the projection part
                                  // starts. The dof values
                                  // are obtained by solving
                                  // a linear system of
                                  // equations.
          if (deg > 0)
            {
                    			  // We start with projection
                    			  // on the higher order edge
                    			  // shape function.
              const std::vector<Polynomials::Polynomial<double> >&
                lobatto_polynomials
                  = Polynomials::Lobatto::generate_complete_basis
                    (this->degree);
              FullMatrix<double> system_matrix (deg, deg);
              std::vector<Polynomials::Polynomial<double> >
                lobatto_polynomials_grad (this->degree);

              for (unsigned int i = 0; i < lobatto_polynomials_grad.size ();
                   ++i)
                lobatto_polynomials_grad[i]
                  = lobatto_polynomials[i + 1].derivative ();

                                  // Set up the system matrix.
                                  // This can be used for all
                                  // edges.
              for (unsigned int i = 0; i < system_matrix.m (); ++i)
                for (unsigned int j = 0; j < system_matrix.n (); ++j)
                  for (unsigned int q_point = 0; q_point < n_edge_points;
                       ++q_point)
                    system_matrix (i, j)
                      += boundary_weights (q_point, j)
                         * lobatto_polynomials_grad[i + 1].value
                           (this->generalized_face_support_points[q_point]
                            (1));

              FullMatrix<double> system_matrix_inv (deg, deg);

              system_matrix_inv.invert (system_matrix);

              const unsigned int
                line_coordinate[GeometryInfo<3>::lines_per_cell]
                  = {1, 1, 0, 0, 1, 1, 0, 0, 2, 2, 2, 2};
              Vector<double> system_rhs (system_matrix.m ());
              Vector<double> solution (system_rhs.size ());

              for (unsigned int line = 0;
                   line < GeometryInfo<dim>::lines_per_cell; ++line)
                {
                                  // Set up the right hand side.
                  system_rhs = 0;

                  for (unsigned int q_point = 0; q_point <= deg; ++q_point)
                    {
                      const double tmp
                        = values[line_coordinate[line]][line * this->degree
                                                        + q_point]
                          - local_dofs[line * this->degree]
                          * this->shape_value_component
                            (line * this->degree,
                             this->generalized_support_points[line
                                                              * this->degree
                                                              + q_point],
                             line_coordinate[line]);

                      for (unsigned int i = 0; i < system_rhs.size (); ++i)
                           system_rhs (i) += boundary_weights (q_point, i)
                                             * tmp;
                    }

                  system_matrix_inv.vmult (solution, system_rhs);

                                  // Add the computed values
                                  // to the resulting vector
                                  // only, if they are not
                                  // too small.
                  for (unsigned int i = 0; i < solution.size (); ++i)
                    if (std::abs (solution (i)) > 1e-14)
                      local_dofs[line * this->degree + i + 1] = solution (i);
                }

                                  // Then we go on to the
                                  // face shape functions.
                                  // Again we set up the
                                  // system matrix and
                                  // use it for both, the
                                  // horizontal and the
                                  // vertical, shape
                                  // functions.
              const std::vector<Polynomials::Polynomial<double> >&
                legendre_polynomials
                  = Polynomials::Legendre::generate_complete_basis (deg);
              const unsigned int n_face_points = n_edge_points * n_edge_points;

              system_matrix.reinit (deg * this->degree, deg * this->degree);
              system_matrix = 0;

              for (unsigned int i = 0; i <= deg; ++i)
                for (unsigned int j = 0; j < deg; ++j)
                  for (unsigned int k = 0; k <= deg; ++k)
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

              system_matrix_inv.reinit (system_matrix.m (),
                                        system_matrix.m ());
              system_matrix_inv.invert (system_matrix);
              solution.reinit (system_matrix.m ());
              system_rhs.reinit (system_matrix.m ());

              const unsigned int
                face_coordinates[GeometryInfo<3>::faces_per_cell][2]
                  = {{1, 2}, {1, 2}, {0, 2}, {0, 2}, {0, 1}, {0, 1}};
              const unsigned int
                edge_indices[GeometryInfo<3>::faces_per_cell][GeometryInfo<3>::lines_per_face]
                  = {{0, 4, 8, 10}, {1, 5, 9, 11}, {2, 6, 8, 9},
                     {3, 7, 10, 11}, {2, 3, 0, 1}, {6, 7, 4, 5}};

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
unsigned int
FE_Nedelec<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}


template class FE_Nedelec<deal_II_dimension>;

DEAL_II_NAMESPACE_CLOSE
