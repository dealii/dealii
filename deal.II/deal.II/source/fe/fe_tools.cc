//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <lac/full_matrix.h>
#include <lac/householder.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <fe/fe_tools.h>
#include <fe/fe.h>
#include <fe/fe_q.h>
#include <fe/fe_q_hierarchical.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_dgp_nonparametric.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_cartesian.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>

#include <boost/shared_ptr.hpp>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif



namespace 
{
  template <int dim>
  inline
  unsigned int
  max_dofs_per_cell (const DoFHandler<dim> &dh) 
  {
    return dh.get_fe().dofs_per_cell;
  }


  template <int dim>
  inline
  unsigned int
  max_dofs_per_face (const DoFHandler<dim> &dh) 
  {
    return dh.get_fe().dofs_per_face;
  }


  template <int dim>
  inline
  unsigned int
  n_components (const DoFHandler<dim> &dh) 
  {
    return dh.get_fe().n_components();
  }


                                   // forwarder function for
                                   // FE::get_interpolation_matrix. we
                                   // will want to call that function
                                   // for arbitrary FullMatrix<T>
                                   // types, but it only accepts
                                   // double arguments. since it is a
                                   // virtual function, this can also
                                   // not be changed. so have a
                                   // forwarder function that calls
                                   // that function directly if
                                   // T==double, and otherwise uses a
                                   // temporary
  template <int dim>
  inline
  void gim_forwarder (const FiniteElement<dim> &fe1,
                      const FiniteElement<dim> &fe2,
                      FullMatrix<double> &interpolation_matrix)
  {
    fe2.get_interpolation_matrix (fe1, interpolation_matrix);
  }

  
  template <int dim, typename number>
  inline
  void gim_forwarder (const FiniteElement<dim> &fe1,
                      const FiniteElement<dim> &fe2,
                      FullMatrix<number> &interpolation_matrix)
  {
    FullMatrix<double> tmp (interpolation_matrix.m(),
                            interpolation_matrix.n());
    fe2.get_interpolation_matrix (fe1, tmp);
    interpolation_matrix = tmp;
  }

				   // return true if the given pattern
				   // string matches the given name at
				   // the first position of the string
  inline
  bool
  match_at_string_start (const std::string &name,
			 const std::string &pattern)
  {
    if (pattern.size() > name.size())
      return false;

    for (unsigned int i=0; i<pattern.size(); ++i)
      if (pattern[i] != name[i])
	return false;

    return true;
  }



				   // read an integer at the position
				   // in "name" indicated by the
				   // second argument, and retun this
				   // integer together with how many
				   // characters it takes in the
				   // string
				   //
				   // if no integer can be read at the
				   // indicated position, return
				   // (-1,-1)
  inline
  std::pair<int, unsigned int>
  get_integer (const std::string &name,
	       const unsigned int position)
  {
    Assert (position < name.size(), ExcInternalError());
    
    const std::string test_string (name.begin()+position,
				   name.end());
    
#ifdef HAVE_STD_STRINGSTREAM
    std::istringstream str(test_string);
#else
    std::istrstream str(test_string.c_str());
#endif

    int i;
    if (str >> i)
      {
					 // compute the number of
					 // digits of i. assuming it
					 // is less than 6 is likely
					 // ok
	if (i<10)
	  return std::make_pair (i, 1U);
	else if (i<100)
	  return std::make_pair (i, 2U);
	else if (i<1000)
	  return std::make_pair (i, 3U);
	else if (i<10000)
	  return std::make_pair (i, 4U);
	else if (i<100000)
	  return std::make_pair (i, 5U);
	else
	  {
	    Assert (false, ExcNotImplemented());
	    return std::make_pair (-1, deal_II_numbers::invalid_unsigned_int);
	  }
      }
    else
      return std::make_pair (-1, deal_II_numbers::invalid_unsigned_int);
  }



				   // return how many characters
				   // starting at the given position
				   // of the string match either the
				   // generic string "<dim>" or the
				   // specialized string with "dim"
				   // replaced with the numeric value
				   // of the template argument
  template <int dim>
  inline
  unsigned int match_dimension (const std::string &name,
				const unsigned int position)
  {
    if (position >= name.size())
      return 0;

    if ((position+5 < name.size())
	&&
	(name[position] == '<')
	&&
	(name[position+1] == 'd')
	&&
	(name[position+2] == 'i')
	&&
	(name[position+3] == 'm')
	&&
	(name[position+4] == '>'))
      return 5;

    Assert (dim<10, ExcNotImplemented());
    const char dim_char = '0'+dim;
    
    if ((position+3 < name.size())
	&&
	(name[position] == '<')
	&&
	(name[position+1] == dim_char)
	&&
	(name[position+2] == '>'))
      return 3;

				     // some other string that doesn't
				     // match
    return 0;
  }
}



template<int dim>
void FETools::compute_component_wise(
  const FiniteElement<dim>& element,
  std::vector<unsigned int>& renumbering,
  std::vector<std::vector<unsigned int> >& comp_start)
{
  Assert(renumbering.size() == element.dofs_per_cell,
	 ExcDimensionMismatch(renumbering.size(),
			      element.dofs_per_cell));
  
  comp_start.resize(element.n_base_elements());
  
  unsigned int k=0;
  for (unsigned int i=0;i<comp_start.size();++i)
    {
      comp_start[i].resize(element.element_multiplicity(i));
      const unsigned int increment
	= element.base_element(i).dofs_per_cell;
      
      for (unsigned int j=0;j<comp_start[i].size();++j)
	{
	  comp_start[i][j] = k;
	  k += increment;
	}
    }
  
				   // For each index i of the
				   // unstructured cellwise
				   // numbering, renumbering
				   // contains the index of the
				   // cell-block numbering
  for (unsigned int i=0;i<element.dofs_per_cell;++i)
    {
      std::pair<std::pair<unsigned int, unsigned int>, unsigned int>
	indices = element.system_to_base_index(i);
      renumbering[i] = comp_start[indices.first.first][indices.first.second]
			     +indices.second;
    }
}



template <int dim, typename number>
void FETools::get_interpolation_matrix (const FiniteElement<dim> &fe1,
                                        const FiniteElement<dim> &fe2,
                                        FullMatrix<number> &interpolation_matrix)
{
  Assert (fe1.n_components() == fe2.n_components(),
	  ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
  Assert(interpolation_matrix.m()==fe2.dofs_per_cell &&
	 interpolation_matrix.n()==fe1.dofs_per_cell,
	 ExcMatrixDimensionMismatch(interpolation_matrix.m(),
				    interpolation_matrix.n(),
				    fe2.dofs_per_cell,
				    fe1.dofs_per_cell));

				   // first try the easy way: maybe
				   // the FE wants to implement things
				   // itself:
  bool fe_implements_interpolation = true;
  try 
    {
      gim_forwarder (fe1, fe2, interpolation_matrix);
    }
  catch (typename FiniteElementBase<dim>::ExcInterpolationNotImplemented &)
    {
                                       // too bad....
      fe_implements_interpolation = false;
    }
  if (fe_implements_interpolation == true)
    return;

				   // uh, so this was not the
				   // case. hm. then do it the hard
				   // way. note that this will only
				   // work if the element is
				   // primitive, so check this first
  Assert (fe1.is_primitive() == true, ExcFEMustBePrimitive());
  Assert (fe2.is_primitive() == true, ExcFEMustBePrimitive());

				   // Initialize FEValues for fe1 at
				   // the unit support points of the
				   // fe2 element.
  const std::vector<Point<dim> > &
    fe2_support_points = fe2.get_unit_support_points ();

  Assert(fe2_support_points.size()==fe2.dofs_per_cell,
	 typename FiniteElementBase<dim>::ExcFEHasNoSupportPoints());

  for (unsigned int i=0; i<fe2.dofs_per_cell; ++i)	
    {
      const unsigned int i1 = fe2.system_to_component_index(i).first;
      for (unsigned int j=0; j<fe1.dofs_per_cell; ++j)
	{
	  const unsigned int j1 = fe1.system_to_component_index(j).first;
	  if (i1==j1)
	    interpolation_matrix(i,j) = fe1.shape_value (j,fe2_support_points[i]);
	  else
	    interpolation_matrix(i,j)=0.;
	}  
    }
}



template <int dim, typename number>
void FETools::get_back_interpolation_matrix(const FiniteElement<dim> &fe1,
					    const FiniteElement<dim> &fe2,
					    FullMatrix<number> &interpolation_matrix)
{
  Assert (fe1.n_components() == fe2.n_components(),
	  ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
  Assert(interpolation_matrix.m()==fe1.dofs_per_cell &&
	 interpolation_matrix.n()==fe1.dofs_per_cell, 
	 ExcMatrixDimensionMismatch(interpolation_matrix.m(),
				    interpolation_matrix.n(),
				    fe1.dofs_per_cell,
				    fe1.dofs_per_cell));
    
  FullMatrix<number> first_matrix (fe2.dofs_per_cell, fe1.dofs_per_cell);
  FullMatrix<number> second_matrix(fe1.dofs_per_cell, fe2.dofs_per_cell);
  
  get_interpolation_matrix(fe1, fe2, first_matrix);
  get_interpolation_matrix(fe2, fe1, second_matrix);

				   // int_matrix=second_matrix*first_matrix
  second_matrix.mmult(interpolation_matrix, first_matrix);
}



template <int dim, typename number>
void FETools::get_interpolation_difference_matrix (const FiniteElement<dim> &fe1,
						   const FiniteElement<dim> &fe2,
						   FullMatrix<number> &difference_matrix)
{
  Assert (fe1.n_components() == fe2.n_components(),
	  ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
  Assert(difference_matrix.m()==fe1.dofs_per_cell &&
	 difference_matrix.n()==fe1.dofs_per_cell, 
	 ExcMatrixDimensionMismatch(difference_matrix.m(),
				    difference_matrix.n(),
				    fe1.dofs_per_cell,
				    fe1.dofs_per_cell));
   
  FullMatrix<number> interpolation_matrix(fe1.dofs_per_cell);
  get_back_interpolation_matrix(fe1, fe2, interpolation_matrix);
  
  for (unsigned int i=0; i<fe1.dofs_per_cell; ++i)
    difference_matrix(i,i) = 1.;
  
				   // compute difference
  difference_matrix.add (-1, interpolation_matrix);
}



template <int dim, typename number>
void FETools::get_projection_matrix (const FiniteElement<dim> &fe1,
				     const FiniteElement<dim> &fe2,
				     FullMatrix<number> &matrix)
{
  Assert (fe1.n_components() == 1, ExcNotImplemented());
  Assert (fe1.n_components() == fe2.n_components(),
	  ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
  Assert(matrix.m()==fe2.dofs_per_cell && matrix.n()==fe1.dofs_per_cell,
	 ExcMatrixDimensionMismatch(matrix.m(), matrix.n(),
				    fe2.dofs_per_cell,
				    fe1.dofs_per_cell));
  matrix = 0;
  
  unsigned int n1 = fe1.dofs_per_cell;
  unsigned int n2 = fe2.dofs_per_cell;

  				   // First, create a local mass matrix for
  				   // the unit cell
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  
				   // Choose a quadrature rule
				   // Gauss is exact up to degree 2n-1
  const unsigned int degree = std::max(fe1.tensor_degree(), fe2.tensor_degree());
  Assert (degree != deal_II_numbers::invalid_unsigned_int,
	  ExcNotImplemented());
  
  QGauss<dim> quadrature(degree+1);
				   // Set up FEValues.
  const UpdateFlags flags = update_values | update_q_points | update_JxW_values;
  FEValues<dim> val1 (fe1, quadrature, update_values);
  val1.reinit (tr.begin_active());
  FEValues<dim> val2 (fe2, quadrature, flags);
  val2.reinit (tr.begin_active());

				   // Integrate and invert mass matrix
				   // This happens in the target space
  FullMatrix<double> mass (n2, n2);

  for (unsigned int k=0;k<quadrature.n_quadrature_points;++k)
    {
      const double w = val2.JxW(k);
      for (unsigned int i=0;i<n2;++i)
	{
	  const double v = val2.shape_value(i,k);
	  for (unsigned int j=0;j<n2;++j)
	    mass(i,j) += w*v * val2.shape_value(j,k);
	}
    }
				   // Gauss-Jordan should be
				   // sufficient since we expect the
				   // mass matrix to be
				   // well-conditioned
  mass.gauss_jordan();
  
				   // Now, test every function of fe1
				   // with test functions of fe2 and
				   // compute the projection of each
				   // unit vector.
  Vector<double> b(n2);
  Vector<double> x(n2);
  
  for (unsigned int j=0;j<n1;++j)
    {
      b = 0.;
      for (unsigned int i=0;i<n2;++i)
        for (unsigned int k=0;k<quadrature.n_quadrature_points;++k)
          {
            const double w = val2.JxW(k);
            const double u = val1.shape_value(j,k);
            const double v = val2.shape_value(i,k);
            b(i) += u*v*w;
          }
      
				       // Multiply by the inverse
      mass.vmult(x,b);
      for (unsigned int i=0;i<n2;++i)
	matrix(i,j) = x(i);
    }
}


template<int dim, typename number>
void
FETools::compute_embedding_matrices(const FiniteElement<dim>& fe,
				    FullMatrix<number>* matrices)
{
  const unsigned int nc = GeometryInfo<dim>::children_per_cell;
  const unsigned int n  = fe.dofs_per_cell;
  const unsigned int nd = fe.n_components();
  const unsigned int degree = fe.degree;
  
  for (unsigned int i=0;i<nc;++i)
    {
      Assert(matrices[i].n() == n, ExcDimensionMismatch(matrices[i].n(),n));
      Assert(matrices[i].m() == n, ExcDimensionMismatch(matrices[i].m(),n));
    }
  
                                   // Set up a meshes, one with a single
                                   // reference cell and refine it once
  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria, 0, 1);
  tria.refine_global(1);

  MappingCartesian<dim> mapping;
  QGauss<dim> q_fine(degree+1);
  const unsigned int nq = q_fine.n_quadrature_points;
  
  FEValues<dim> fine (mapping, fe, q_fine,
		      update_q_points | update_JxW_values | update_values);
  
				   // We search for the polynomial on
				   // the small cell, being equal to
				   // the coarse polynomial in all
				   // quadrature points.
				    
				   // First build the matrix for this
				   // least squares problem. This
				   // contains the values of the fine
				   // cell polynomials in the fine
				   // cell grid points.
				    
				   // This matrix is the same for all
				   // children.
  fine.reinit(tria.begin_active());
  FullMatrix<number> A(nq*nd, n);
  for (unsigned int d=0;d<nd;++d)
    for (unsigned int k=0;k<nq;++k)
      for (unsigned int j=0;j<n;++j)
	A(k*nd+d,j) = fine.shape_value_component(j,k,d);

  Householder<double> H(A);
  
  Vector<number> v_coarse(nq*nd);
  Vector<number> v_fine(n);
  
  unsigned int cell_number = 0;
  for (typename Triangulation<dim>::active_cell_iterator fine_cell
         = tria.begin_active();
       fine_cell != tria.end(); ++fine_cell, ++cell_number)
    {
      fine.reinit(fine_cell);

                                       // evaluate on the coarse cell (which
                                       // is the first -- inactive -- cell on
                                       // the lowest level of the
                                       // triangulation we have created)
      const Quadrature<dim> q_coarse (fine.get_quadrature_points(),
                                      fine.get_JxW_values());
      FEValues<dim> coarse (mapping, fe, q_coarse, update_values);
      coarse.reinit(tria.begin(0));

      FullMatrix<double> &this_matrix = matrices[cell_number];
      
				       // Compute this once for each
				       // coarse grid basis function
      for (unsigned int i=0;i<n;++i)
	{
					   // The right hand side of
					   // the least squares
					   // problem consists of the
					   // function values of the
					   // coarse grid function in
					   // each quadrature point.
	  for (unsigned int d=0;d<nd;++d)
	    for (unsigned int k=0;k<nq;++k)
	      v_coarse(k*nd+d) = coarse.shape_value_component(i,k,d);

					   // solve the least squares
					   // problem.
	  const double result = H.least_squares(v_fine, v_coarse);
	  Assert (result < 1.e-12, ExcLeastSquaresError(result));
	    
					   // Copy into the result
					   // matrix. Since the matrix
					   // maps a coarse grid
					   // function to a fine grid
					   // function, the columns
					   // are fine grid.
	  for (unsigned int j=0;j<n;++j)
	    this_matrix(j,i) = v_fine(j);
	}
				       // Remove small entries from
				       // the matrix
      for (unsigned int i=0; i<this_matrix.m(); ++i)
	for (unsigned int j=0; j<this_matrix.n(); ++j)
	  if (std::fabs(this_matrix(i,j)) < 1e-12)
	    this_matrix(i,j) = 0.;
    }
  Assert (cell_number == GeometryInfo<dim>::children_per_cell,
          ExcInternalError());
}



//TODO[GK]: this function does not work yet.
template<int dim, typename number>
void
FETools::compute_projection_matrices(const FiniteElement<dim>& fe,
				     FullMatrix<number>* matrices)
{
  Assert(false, ExcNotImplemented());
  const unsigned int nc = GeometryInfo<dim>::children_per_cell;
  const unsigned int n  = fe.dofs_per_cell;
  const unsigned int nd = fe.n_components();
  const unsigned int degree = fe.degree;
  
  for (unsigned int i=0;i<nc;++i)
    {
      Assert(matrices[i].n() == n, ExcDimensionMismatch(matrices[i].n(),n));
      Assert(matrices[i].m() == n, ExcDimensionMismatch(matrices[i].m(),n));
    }
  
  Triangulation<dim> tr;
  GridGenerator::hyper_cube (tr, 0, 1);
  tr.refine_global(1);
  
  MappingCartesian<dim> mapping;
  QGauss<dim> q_fine(degree+1);
  const unsigned int nq = q_fine.n_quadrature_points;
  
  FEValues<dim> coarse (mapping, fe, q_fine,
			update_q_points | update_JxW_values | update_values);
  FEValues<dim> fine (mapping, fe, q_fine,
		      update_q_points | update_JxW_values | update_values);
  
  typename Triangulation<dim>::cell_iterator coarse_cell
    = tr.begin(0);
  typename Triangulation<dim>::cell_iterator fine_cell;

				   // Compute the coarse level mass
				   // matrix
  coarse.reinit(coarse_cell);
  FullMatrix<number> A(n, n);
  for (unsigned int k=0;k<nq;++k)
    for (unsigned int i=0;i<n;++i)
      for (unsigned int j=0;j<n;++j)
	if (fe.is_primitive())
	  A(i,j) = coarse.JxW(k)
		   * coarse.shape_value(i,k)
		   * coarse.shape_value(j,k);
	else
	  for (unsigned int d=0;d<nd;++d)
	    A(i,j) = coarse.JxW(k)
		     * coarse.shape_value_component(i,k,d)
		     * coarse.shape_value_component(j,k,d);
  
  Householder<double> H(A);
  
  Vector<number> v_coarse(n);
  Vector<number> v_fine(n);
  
  for (unsigned int cell_number=0;cell_number<GeometryInfo<dim>::children_per_cell;++cell_number)
    {
      FullMatrix<double> &this_matrix = matrices[cell_number];
      
				       // Compute right hand side,
				       // which is a fine level basis
				       // function tested with the
				       // coarse level functions.
      fine.reinit(coarse_cell->child(cell_number));
      Quadrature<dim> q_coarse (fine.get_quadrature_points(),
				fine.get_JxW_values());
      FEValues<dim> coarse (mapping, fe, q_coarse, update_values);
      coarse.reinit(coarse_cell);
      
				       // Build RHS

				       // Outer loop over all fine
				       // grid shape functions phi_j
      for (unsigned int j=0;j<fe.dofs_per_cell;++j)
	{
					   // Loop over all quadrature points
	  for (unsigned int k=0;k<fine.n_quadrature_points;++k)
	    {
					       // integrate the scalar
					       // product
					       // (phi_i,phi_j) for
					       // all coarse shape
					       // functions to get the
					       // right hand side
	      for (unsigned int i=0;i<fe.dofs_per_cell;++i)
		{
		  if (fe.is_primitive())
		    v_fine(i) += fine.JxW(k)
				 * coarse.shape_value(i,k)
				 * fine.shape_value(j,k);
		  else
		    for (unsigned int d=0;d<nd;++d)
		      v_fine(i) += fine.JxW(k)
				   * coarse.shape_value_component(i,k,d)
				   * fine.shape_value_component(j,k,d);
		}
	    }
					   // RHS ready. Solve system
					   // and enter row into
					   // matrix
	  H.least_squares(v_coarse, v_fine);
	  for (unsigned int i=0;i<fe.dofs_per_cell;++i)
	    this_matrix(j,i) = v_coarse(i);
	}
      
				       // Remove small entries from
				       // the matrix
      for (unsigned int i=0; i<this_matrix.m(); ++i)
 	for (unsigned int j=0; j<this_matrix.n(); ++j)
 	  if (std::fabs(this_matrix(i,j)) < 1e-12)
 	    this_matrix(i,j) = 0.;
    }
}


template <int dim,
          class InVector, class OutVector>
void
FETools::interpolate(const DoFHandler<dim> &dof1,
                     const InVector        &u1,
                     const DoFHandler<dim> &dof2,
                     OutVector             &u2)
{
  ConstraintMatrix dummy;
  dummy.close();
  interpolate(dof1, u1, dof2, dummy, u2);
}

  
template <int dim, class InVector, class OutVector>
void
FETools::interpolate(const DoFHandler<dim>  &dof1,
                     const InVector         &u1,
                     const DoFHandler<dim>  &dof2,
                     const ConstraintMatrix &constraints,
                     OutVector              &u2)
{
  Assert(&dof1.get_tria() == &dof2.get_tria(), ExcTriangulationMismatch());
  Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(),
                              dof2.get_fe().n_components()));
  Assert(u1.size()==dof1.n_dofs(),
         ExcDimensionMismatch(u1.size(), dof1.n_dofs()));  
  Assert(u2.size()==dof2.n_dofs(),
         ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

				   // for continuous elements on grids
				   // with hanging nodes we need
				   // hanging node
				   // constraints. Consequentely, if
				   // there are no constraints then
				   // hanging nodes are not allowed.
  const bool hanging_nodes_not_allowed
    = ((dof2.get_fe().dofs_per_vertex != 0) &&
       (constraints.n_constraints() == 0));
  
  const unsigned int dofs_per_cell1=dof1.get_fe().dofs_per_cell;
  const unsigned int dofs_per_cell2=dof2.get_fe().dofs_per_cell;
  
  Vector<typename OutVector::value_type> u1_local(dofs_per_cell1);
  Vector<typename OutVector::value_type> u2_local(dofs_per_cell2);

  FullMatrix<double> interpolation_matrix(dofs_per_cell2,
					  dofs_per_cell1);
  FETools::get_interpolation_matrix(dof1.get_fe(), dof2.get_fe(),
				    interpolation_matrix);
  
  typename DoFHandler<dim>::active_cell_iterator cell1 = dof1.begin_active(),
						 endc1 = dof1.end(),
						 cell2 = dof2.begin_active(),
						 endc2 = dof2.end();

  std::vector<unsigned int> touch_count(dof2.n_dofs(),0);
  std::vector<unsigned int> dofs (dofs_per_cell2);
  u2 = 0;

  for (; cell1!=endc1; ++cell1, ++cell2) 
    {
      if (hanging_nodes_not_allowed)
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  Assert (cell1->at_boundary(face) ||
		  cell1->neighbor(face)->level() == cell1->level(),
		  ExcHangingNodesNotAllowed(0));
      
      cell1->get_dof_values(u1, u1_local);
      interpolation_matrix.vmult(u2_local, u1_local);
      cell2->get_dof_indices(dofs);
      for (unsigned int i=0; i<dofs_per_cell2; ++i)
	{
	  u2(dofs[i])+=u2_local(i);
	  ++touch_count[dofs[i]];
	}
    }
                                   // cell1 is at the end, so should
                                   // be cell2
  Assert (cell2 == endc2, ExcInternalError());
  
				   // when a discontinuous element is
				   // interpolated to a continuous
				   // one, we take the mean values.
  for (unsigned int i=0; i<dof2.n_dofs(); ++i)
    {
      Assert(touch_count[i]!=0, ExcInternalError());
      u2(i) /= touch_count[i];
    }

				   // Apply hanging node constraints.
  constraints.distribute(u2);
}



template <int dim, class InVector, class OutVector>
void
FETools::back_interpolate(const DoFHandler<dim>    &dof1,
                          const InVector           &u1,
                          const FiniteElement<dim> &fe2,
                          OutVector                &u1_interpolated)
{
  Assert(dof1.get_fe().n_components() == fe2.n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), fe2.n_components()));
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));  
  Assert(u1_interpolated.size()==dof1.n_dofs(),
	 ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));
  
				   // For continuous elements on grids
				   // with hanging nodes we need
				   // hanging node
				   // constraints. Consequently, when
				   // the elements are continuous no
				   // hanging node constraints are
				   // allowed.
  const bool hanging_nodes_not_allowed=
    (dof1.get_fe().dofs_per_vertex != 0) || (fe2.dofs_per_vertex != 0);

  const unsigned int dofs_per_cell1=dof1.get_fe().dofs_per_cell;

  Vector<typename OutVector::value_type> u1_local(dofs_per_cell1);
  Vector<typename OutVector::value_type> u1_int_local(dofs_per_cell1);
  
  typename DoFHandler<dim>::active_cell_iterator cell = dof1.begin_active(),
						 endc = dof1.end();

  FullMatrix<double> interpolation_matrix(dofs_per_cell1, dofs_per_cell1);
  FETools::get_back_interpolation_matrix(dof1.get_fe(), fe2,
					 interpolation_matrix);
  for (; cell!=endc; ++cell) 
    {
      if (hanging_nodes_not_allowed)
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  Assert (cell->at_boundary(face) ||
		  cell->neighbor(face)->level() == cell->level(),
		  ExcHangingNodesNotAllowed(0));

      cell->get_dof_values(u1, u1_local);
      interpolation_matrix.vmult(u1_int_local, u1_local);
      cell->set_dof_values(u1_int_local, u1_interpolated);
    }
}


  
template <int dim, class InVector, class OutVector>
void FETools::back_interpolate(const DoFHandler<dim> &dof1,
			       const ConstraintMatrix &constraints1,
			       const InVector &u1,
			       const DoFHandler<dim> &dof2,
			       const ConstraintMatrix &constraints2,
			       OutVector &u1_interpolated)
{
				   // For discontinuous elements
				   // without constraints take the
				   // simpler version of the
				   // back_interpolate function.
  if (dof1.get_fe().dofs_per_vertex==0 && dof2.get_fe().dofs_per_vertex==0
      && constraints1.n_constraints()==0 && constraints2.n_constraints()==0)
    back_interpolate(dof1, u1, dof2.get_fe(), u1_interpolated);
  else
    {
      Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
	     ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
      Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));  
      Assert(u1_interpolated.size()==dof1.n_dofs(),
	     ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));
      
				       // For continuous elements
				       // first interpolate to dof2,
				       // taking into account
				       // constraints2, and then
				       // interpolate back to dof1
				       // taking into account
				       // constraints1
      Vector<typename OutVector::value_type> u2(dof2.n_dofs());
      interpolate(dof1, u1, dof2, constraints2, u2);
      interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
    }
}


  
template <int dim, class InVector, class OutVector>
void FETools::interpolation_difference (const DoFHandler<dim> &dof1,
					const InVector &u1,
					const FiniteElement<dim> &fe2,
					OutVector &u1_difference)
{
  Assert(dof1.get_fe().n_components() == fe2.n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), fe2.n_components()));
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));  
  Assert(u1_difference.size()==dof1.n_dofs(),
	 ExcDimensionMismatch(u1_difference.size(), dof1.n_dofs()));
  
				   // For continuous elements on grids
				   // with hanging nodes we need
				   // hnaging node
				   // constraints. Consequently, when
				   // the elements are continuous no
				   // hanging node constraints are
				   // allowed.
  const bool hanging_nodes_not_allowed=
    (dof1.get_fe().dofs_per_vertex != 0) || (fe2.dofs_per_vertex != 0);

  const unsigned int dofs_per_cell=dof1.get_fe().dofs_per_cell;

  Vector<typename OutVector::value_type> u1_local(dofs_per_cell);
  Vector<typename OutVector::value_type> u1_diff_local(dofs_per_cell);
  
  FullMatrix<double> difference_matrix(dofs_per_cell, dofs_per_cell);
  FETools::get_interpolation_difference_matrix(dof1.get_fe(), fe2,
					       difference_matrix);
  
  typename DoFHandler<dim>::active_cell_iterator cell = dof1.begin_active(),
						 endc = dof1.end();
  
  for (; cell!=endc; ++cell)
    {
      if (hanging_nodes_not_allowed)
	for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	  Assert (cell->at_boundary(face) ||
		  cell->neighbor(face)->level() == cell->level(),
		  ExcHangingNodesNotAllowed(0));

      cell->get_dof_values(u1, u1_local);
      difference_matrix.vmult(u1_diff_local, u1_local);
      cell->set_dof_values(u1_diff_local, u1_difference);
    }
}



template <int dim, class InVector, class OutVector>
void FETools::interpolation_difference(const DoFHandler<dim> &dof1,
				       const ConstraintMatrix &constraints1,
				       const InVector &u1,
				       const DoFHandler<dim> &dof2,
				       const ConstraintMatrix &constraints2,
				       OutVector &u1_difference)
{
 				   // For discontinuous elements
				   // without constraints take the
				   // cheaper version of the
				   // interpolation_difference function.
  if (dof1.get_fe().dofs_per_vertex==0 && dof2.get_fe().dofs_per_vertex==0
      && constraints1.n_constraints()==0 && constraints2.n_constraints()==0)
    interpolation_difference(dof1, u1, dof2.get_fe(), u1_difference);
  else
    {
      back_interpolate(dof1, constraints1, u1, dof2, constraints2, u1_difference);
      u1_difference.sadd(-1, u1);
    }
}

  
template <int dim, class InVector, class OutVector>
void FETools::project_dg(const DoFHandler<dim> &dof1,
			 const InVector &u1,
			 const DoFHandler<dim> &dof2,
			 OutVector &u2)
{
  Assert(&dof1.get_tria()==&dof2.get_tria(), ExcTriangulationMismatch());
  Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));  
  Assert(u2.size()==dof2.n_dofs(), ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

  typename DoFHandler<dim>::active_cell_iterator cell1 = dof1.begin_active();
  typename DoFHandler<dim>::active_cell_iterator cell2 = dof2.begin_active();
  typename DoFHandler<dim>::active_cell_iterator end = dof2.end();

  const unsigned int n1 = dof1.get_fe().dofs_per_cell;
  const unsigned int n2 = dof2.get_fe().dofs_per_cell;
  
  Vector<double> u1_local(n1);
  Vector<double> u2_local(n2);
  std::vector<unsigned int> dofs(n2);
  
  FullMatrix<double> matrix(n2,n1);
  get_projection_matrix(dof1.get_fe(), dof2.get_fe(), matrix);
  
  while (cell2 != end)
    {
      cell1->get_dof_values(u1, u1_local);
      matrix.vmult(u2_local, u1_local);
      cell2->get_dof_indices(dofs);
      for (unsigned int i=0; i<n2; ++i)
	{
	  u2(dofs[i])+=u2_local(i);
	}
     
      ++cell1;
      ++cell2;
    }
}

  
template <int dim, class InVector, class OutVector>
void FETools::extrapolate(const DoFHandler<dim> &dof1,
			  const InVector &u1,
			  const DoFHandler<dim> &dof2,
			  OutVector &u2)
{
  ConstraintMatrix dummy;
  dummy.close();
  extrapolate(dof1, u1, dof2, dummy, u2);
}



template <int dim, class InVector, class OutVector>
void FETools::extrapolate(const DoFHandler<dim> &dof1,
			  const InVector &u1,
			  const DoFHandler<dim> &dof2,
			  const ConstraintMatrix &constraints,
			  OutVector &u2)
{
  Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
  Assert(&dof1.get_tria()==&dof2.get_tria(), ExcTriangulationMismatch());
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
  Assert(u2.size()==dof2.n_dofs(), ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

  OutVector u3;
  u3.reinit(u2);
  interpolate(dof1, u1, dof2, constraints, u3);

  const unsigned int dofs_per_cell  = dof2.get_fe().dofs_per_cell;
  Vector<typename OutVector::value_type> dof_values(dofs_per_cell);

				   // make sure that each cell on the
				   // coarsest level is at least once
				   // refined. otherwise, we can't
				   // treat these cells and would
				   // generate a bogus result
  {
    typename DoFHandler<dim>::cell_iterator cell = dof2.begin(0),
					    endc = dof2.end(0);
    for (; cell!=endc; ++cell)
      Assert (cell->has_children(), ExcGridNotRefinedAtLeastOnce());
  } 

				   // then traverse grid bottom up
  for (unsigned int level=0; level<dof1.get_tria().n_levels()-1; ++level)
    {
      typename DoFHandler<dim>::cell_iterator cell=dof2.begin(level),
					      endc=dof2.end(level);

      for (; cell!=endc; ++cell)
	if (!cell->active())
	  {
					     // check whether this
					     // cell has active
					     // children
	    bool active_children=false;
	    for (unsigned int child_n=0;
		 child_n<GeometryInfo<dim>::children_per_cell; ++child_n)
	      if (cell->child(child_n)->active())
		{
		  active_children=true;
		  break;
		}

					     // if there are active
					     // children, the we have
					     // to work on this
					     // cell. get the data
					     // from the one vector
					     // and set it on the
					     // other
	    if (active_children)
	      {
		cell->get_interpolated_dof_values(u3, dof_values);
		cell->set_dof_values_by_interpolation(dof_values, u2);
	      }
	  }
    }

				   // Apply hanging node constraints.
  constraints.distribute(u2);
}



template <int dim>
FiniteElement<dim> *
FETools::get_fe_from_name (const std::string &name)
{
				   // get the finite element that
				   // would be created from the given
				   // string at the first position, as
				   // well as how many characters of
				   // the string were eaten
  const std::pair<FiniteElement<dim>*, unsigned int>
    tmp = get_fe_from_name_aux<dim> (name);
				   // make sure that we took all
				   // characters in the name,
				   // i.e. that there is not some junk
				   // left over at the end of which we
				   // didn't know what to do
				   // with. make sure we don't create
				   // a memory leak here
  if (tmp.second != name.size())
    {
      delete tmp.first;
      AssertThrow (false, ExcInvalidFEName (name));
    }

				   // otherwise, return the just
				   // created pointer
  return tmp.first;
}



template <int dim>
std::pair<FiniteElement<dim> *, unsigned int>
FETools::get_fe_from_name_aux (const std::string &name)
{  
				   // so, let's see what's at position
				   // 0 of this string, and create a
				   // respective finite element
				   //
				   // start with the longest names, to
				   // make sure we don't match FE_Q
				   // when it's actually a
				   // FE_Q_Hierarchic
  if (match_at_string_start (name, std::string("FE_Q_Hierarchical")))
    {
      unsigned int position = std::string("FE_Q_Hierarchical").size();
				       // as described in the
				       // documentation, at this point
				       // we may have either a) no
				       // dimension specification, b)
				       // the correct dimension (like
				       // <2>), or c) a generic
				       // <dim>. check how many of
				       // these characters match, and
				       // advance the cursor by the
				       // respective amount
      position += match_dimension<dim> (name, position);
      
				       // make sure the next character
				       // is an opening parenthesis
      AssertThrow (name[position] == '(', ExcInvalidFEName(name));
      ++position;
				       // next thing is to parse the
				       // degree of the finite element
      const std::pair<int,unsigned int> tmp = get_integer (name, position);
      AssertThrow (tmp.first>=0, ExcInvalidFEName(name));
      position += tmp.second;

				       // make sure the next character
				       // is an closing parenthesis
      AssertThrow (name[position] == ')', ExcInvalidFEName(name));
      ++position;

				       // ok, everything seems
				       // good. so create finite
				       // element and return position
				       // count
      return std::make_pair (static_cast<FiniteElement<dim>*>
			     (new FE_Q_Hierarchical<dim>(tmp.first)),
			     position);
    }
				   // check other possibilities in
				   // exactly the same way
  else if (match_at_string_start (name, std::string("FE_RaviartThomas")))
    {
      unsigned int position = std::string("FE_RaviartThomas").size();
      position += match_dimension<dim> (name, position);
      AssertThrow (name[position] == '(', ExcInvalidFEName(name));
      ++position;
      const std::pair<int,unsigned int> tmp = get_integer (name, position);
      AssertThrow (tmp.first>=0, ExcInvalidFEName(name));
      position += tmp.second;
      AssertThrow (name[position] == ')', ExcInvalidFEName(name));
      ++position;
      return std::make_pair (static_cast<FiniteElement<dim>*>
			     (new FE_RaviartThomas<dim>(tmp.first)),
			     position);
    }
  else if (match_at_string_start (name, std::string("FE_Nedelec")))
    {
      unsigned int position = std::string("FE_Nedelec").size();
      position += match_dimension<dim> (name, position);
      AssertThrow (name[position] == '(', ExcInvalidFEName(name));
      ++position;
      const std::pair<int,unsigned int> tmp = get_integer (name, position);
      AssertThrow (tmp.first>=0, ExcInvalidFEName(name));
      position += tmp.second;
      AssertThrow (name[position] == ')', ExcInvalidFEName(name));
      ++position;
      return std::make_pair (static_cast<FiniteElement<dim>*>
			     (new FE_Nedelec<dim>(tmp.first)),
			     position);
    }
  else if (match_at_string_start (name, std::string("FE_DGPNonparametric")))
    {
      unsigned int position = std::string("FE_DGPNonparametric").size();
      position += match_dimension<dim> (name, position);
      AssertThrow (name[position] == '(', ExcInvalidFEName(name));
      ++position;
      const std::pair<int,unsigned int> tmp = get_integer (name, position);
      AssertThrow (tmp.first>=0, ExcInvalidFEName(name));
      position += tmp.second;
      AssertThrow (name[position] == ')', ExcInvalidFEName(name));
      ++position;
      return std::make_pair (static_cast<FiniteElement<dim>*>
			     (new FE_DGPNonparametric<dim>(tmp.first)),
			     position);
    }
  else if (match_at_string_start (name, std::string("FE_DGP")))
    {
      unsigned int position = std::string("FE_DGP").size();
      position += match_dimension<dim> (name, position);
      AssertThrow (name[position] == '(', ExcInvalidFEName(name));
      ++position;
      const std::pair<int,unsigned int> tmp = get_integer (name, position);
      AssertThrow (tmp.first>=0, ExcInvalidFEName(name));
      position += tmp.second;
      AssertThrow (name[position] == ')', ExcInvalidFEName(name));
      ++position;
      return std::make_pair (static_cast<FiniteElement<dim>*>(new FE_DGP<dim>(tmp.first)),
			     position);
    }
  else if (match_at_string_start (name, std::string("FE_DGQ")))
    {
      unsigned int position = std::string("FE_DGQ").size();
      position += match_dimension<dim> (name, position);
      AssertThrow (name[position] == '(', ExcInvalidFEName(name));
      ++position;
      const std::pair<int,unsigned int> tmp = get_integer (name, position);
      AssertThrow (tmp.first>=0, ExcInvalidFEName(name));
      position += tmp.second;
      AssertThrow (name[position] == ')', ExcInvalidFEName(name));
      ++position;
      return std::make_pair (static_cast<FiniteElement<dim>*> (new FE_DGQ<dim>(tmp.first)),
			     position);
    }
  else if (match_at_string_start (name, std::string("FE_Q")))
    {
      unsigned int position = std::string("FE_Q").size();
      position += match_dimension<dim> (name, position);
      AssertThrow (name[position] == '(', ExcInvalidFEName(name));
      ++position;
      const std::pair<int,unsigned int> tmp = get_integer (name, position);
      AssertThrow (tmp.first>=0, ExcInvalidFEName(name));
      position += tmp.second;
      AssertThrow (name[position] == ')', ExcInvalidFEName(name));
      ++position;
      return std::make_pair (static_cast<FiniteElement<dim>*>(new FE_Q<dim>(tmp.first)),
			     position);
    }
  
  else
				     // now things get a little more
				     // complicated: FESystem. it's
				     // more complicated, since we
				     // have to figure out what the
				     // base elements are. this can
				     // only be done recursively
    if (match_at_string_start (name, std::string("FESystem")))
      {
	unsigned int position = std::string("FESystem").size();
	position += match_dimension<dim> (name, position);

					 // FESystem puts the names of
					 // the basis elements into
					 // square brackets
	AssertThrow (name[position] == '[', ExcInvalidFEName(name));
	++position;

					 // next we have to get at the
					 // base elements. start with
					 // the first. wrap the whole
					 // block into try-catch to
					 // make sure we destroy the
					 // pointers we got from
					 // recursive calls if one of
					 // these calls should throw
					 // an exception
	std::vector<FiniteElement<dim>*> base_fes;
	std::vector<unsigned int>        base_multiplicities;
	try
	  {
	    do
	      {
						 // first get the
						 // element at this
						 // position of the
						 // string, i.e. of
						 // the substring
						 // ranging from the
						 // present position
						 // to the end
		const std::pair<FiniteElement<dim> *, unsigned int> tmp_x
		  = get_fe_from_name_aux<dim> (std::string(name.begin()+position,
							   name.end()));
		base_fes.push_back (tmp_x.first);
		position += tmp_x.second;
		
						 // next check whether
						 // FESystem placed a
						 // multiplicity after
						 // the element name
		if (name[position] == '^')
		  {
						     // yes. this is
						     // the case. move
						     // the cursor
						     // beyond the '^'
						     // and read this
						     // multiplicity
		    ++position;
		    const std::pair<int,unsigned int> tmp = get_integer (name,
									 position);
		    AssertThrow (tmp.first>=0, ExcInvalidFEName(name));
		    position += tmp.second;
		    base_multiplicities.push_back (tmp.first);
		  }
		else
						   // no, so
						   // multiplicity is
						   // 1
		  base_multiplicities.push_back (1);

						 // so that's it for
						 // this base
						 // element. base
						 // elements are
						 // separated by '-',
						 // and the list is
						 // terminated by ']',
						 // so loop while the
						 // next character is
						 // '-'
	      }
	    while (name[position++] == '-');

					     // so we got to the end
					     // of the '-' separated
					     // list. make sure that
					     // we actually had a ']'
					     // there
	    AssertThrow (name[position-1] == ']', ExcInvalidFEName(name));

					     // just one more sanity check
	    Assert ((base_fes.size() == base_multiplicities.size())
		    &&
		    (base_fes.size() > 0),
		    ExcInternalError());
	    
					     // ok, apparently
					     // everything went ok. so
					     // generate the composed
					     // element
	    FiniteElement<dim> *system_element = 0;
	    switch (base_fes.size())
	      {
		case 1:
		{
		  system_element = new FESystem<dim>(*base_fes[0],
						     base_multiplicities[0]);
		  break;
		}

		case 2:
		{
		  system_element = new FESystem<dim>(*base_fes[0],
						     base_multiplicities[0],
						     *base_fes[1],
						     base_multiplicities[1]);
		  break;
		}
		
		case 3:
		{
		  system_element = new FESystem<dim>(*base_fes[0],
						     base_multiplicities[0],
						     *base_fes[1],
						     base_multiplicities[1],
						     *base_fes[2],
						     base_multiplicities[2]);
		  break;
		}

		default:
		      Assert (false, ExcNotImplemented());
	      }

					     // now we don't need the
					     // list of base elements
					     // any more
	    for (unsigned int i=0; i<base_fes.size(); ++i)
	      delete base_fes[i];
	    
					     // finally return our findings
	    return std::make_pair (system_element, position);
	  }
	catch (...)
	  {
					     // ups, some exception
					     // was thrown. prevent a
					     // memory leak, and then
					     // pass on the exception
					     // to the caller
	    for (unsigned int i=0; i<base_fes.size(); ++i)
	      delete base_fes[i];
	    throw;
	  }

					 // this is a place where we
					 // should really never get,
					 // since above we have either
					 // returned from the
					 // try-clause, or have
					 // re-thrown in the catch
					 // clause. check that we
					 // never get here
	Assert (false, ExcInternalError());
      }
  
    
				   // hm, if we have come thus far, we
				   // didn't know what to do with the
				   // string we got. so do as the docs
				   // say: raise an exception
  AssertThrow (false, ExcInvalidFEName(name));

				   // make some compilers happy that
				   // do not realize that we can't get
				   // here after throwing
  return std::pair<FiniteElement<dim> *, unsigned int> (0,0);
}



template <int dim>
void
FETools::
compute_projection_from_quadrature_points_matrix (const FiniteElement<dim> &fe,
                                                  const Quadrature<dim>    &lhs_quadrature,
                                                  const Quadrature<dim>    &rhs_quadrature,
                                                  FullMatrix<double>       &X)
{
  Assert (fe.n_components() == 1, ExcNotImplemented());

                                   // first build the matrices M and Q
                                   // described in the documentation
  FullMatrix<double> M (fe.dofs_per_cell, fe.dofs_per_cell);
  FullMatrix<double> Q (fe.dofs_per_cell, rhs_quadrature.n_quadrature_points);

  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
      for (unsigned int q=0; q<lhs_quadrature.n_quadrature_points; ++q)
        M(i,j) += fe.shape_value (i, lhs_quadrature.point(q)) *
                  fe.shape_value (j, lhs_quadrature.point(q)) *
                  lhs_quadrature.weight(q);

  for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
    for (unsigned int q=0; q<rhs_quadrature.n_quadrature_points; ++q)
      Q(i,q) += fe.shape_value (i, rhs_quadrature.point(q)) *
                rhs_quadrature.weight(q);

                                   // then invert M
  FullMatrix<double> M_inverse (fe.dofs_per_cell, fe.dofs_per_cell);
  M_inverse.invert (M);

                                   // finally compute the result
  X.reinit (fe.dofs_per_cell, rhs_quadrature.n_quadrature_points);
  M_inverse.mmult (X, Q);

  Assert (X.m() == fe.dofs_per_cell, ExcInternalError());
  Assert (X.n() == rhs_quadrature.n_quadrature_points, ExcInternalError());
}



template <int dim>
void
FETools::
compute_interpolation_to_quadrature_points_matrix (const FiniteElement<dim> &fe,
                                                   const Quadrature<dim>    &quadrature,
                                                   FullMatrix<double>       &I_q)
{
  Assert (fe.n_components() == 1, ExcNotImplemented());
  Assert (I_q.m() == quadrature.n_quadrature_points,
          ExcMessage ("Wrong matrix size"));
  Assert (I_q.n() == fe.dofs_per_cell, ExcMessage ("Wrong matrix size"));

  for (unsigned int q=0; q<quadrature.n_quadrature_points; ++q)
    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      I_q(q,i) = fe.shape_value (i, quadrature.point(q));
}




/*-------------- Explicit Instantiations -------------------------------*/


template
void FETools::compute_component_wise(
  const FiniteElement<deal_II_dimension>& element,
  std::vector<unsigned int>&, std::vector<std::vector<unsigned int> >&);
template
void FETools::get_interpolation_matrix<deal_II_dimension>
(const FiniteElement<deal_II_dimension> &,
 const FiniteElement<deal_II_dimension> &,
 FullMatrix<double> &);
template
void FETools::get_back_interpolation_matrix<deal_II_dimension>
(const FiniteElement<deal_II_dimension> &,
 const FiniteElement<deal_II_dimension> &,
 FullMatrix<double> &);
template
void FETools::get_interpolation_difference_matrix<deal_II_dimension>
(const FiniteElement<deal_II_dimension> &,
 const FiniteElement<deal_II_dimension> &,
 FullMatrix<double> &);
template
void FETools::get_interpolation_matrix<deal_II_dimension>
(const FiniteElement<deal_II_dimension> &,
 const FiniteElement<deal_II_dimension> &,
 FullMatrix<float> &);
template
void FETools::get_back_interpolation_matrix<deal_II_dimension>
(const FiniteElement<deal_II_dimension> &,
 const FiniteElement<deal_II_dimension> &,
 FullMatrix<float> &);
template
void FETools::get_interpolation_difference_matrix<deal_II_dimension>
(const FiniteElement<deal_II_dimension> &,
 const FiniteElement<deal_II_dimension> &,
 FullMatrix<float> &);

template
void FETools::get_projection_matrix<deal_II_dimension>
(const FiniteElement<deal_II_dimension> &,
 const FiniteElement<deal_II_dimension> &,
 FullMatrix<double> &);

template
void FETools::compute_embedding_matrices<deal_II_dimension>
(const FiniteElement<deal_II_dimension> &, FullMatrix<double>*);

template
void FETools::compute_projection_matrices<deal_II_dimension>
(const FiniteElement<deal_II_dimension> &, FullMatrix<double>*);

template
void FETools::interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<double> &,
 const DoFHandler<deal_II_dimension> &, Vector<double> &);
template
void FETools::interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<double> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 Vector<double> &);
template
void FETools::back_interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<double> &,
 const FiniteElement<deal_II_dimension> &, Vector<double> &);
template
void FETools::back_interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 const Vector<double> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 Vector<double> &);
template
void FETools::interpolation_difference<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<double> &,
 const FiniteElement<deal_II_dimension> &, Vector<double> &);
template
void FETools::interpolation_difference<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 const Vector<double> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 Vector<double> &);
template
void FETools::project_dg<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<double> &,
 const DoFHandler<deal_II_dimension> &, Vector<double> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<double> &,
 const DoFHandler<deal_II_dimension> &, Vector<double> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<double> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 Vector<double> &);


template
void FETools::interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<float> &,
 const DoFHandler<deal_II_dimension> &, Vector<float> &);
template
void FETools::interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<float> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 Vector<float> &);
template
void FETools::back_interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<float> &,
 const FiniteElement<deal_II_dimension> &, Vector<float> &);
template
void FETools::back_interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 const Vector<float> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 Vector<float> &);
template
void FETools::interpolation_difference<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<float> &,
 const FiniteElement<deal_II_dimension> &, Vector<float> &);
template
void FETools::interpolation_difference<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 const Vector<float> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 Vector<float> &);
template
void FETools::project_dg<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<float> &,
 const DoFHandler<deal_II_dimension> &, Vector<float> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<float> &,
 const DoFHandler<deal_II_dimension> &, Vector<float> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const Vector<float> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 Vector<float> &);


template
void FETools::interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<double> &,
 const DoFHandler<deal_II_dimension> &, BlockVector<double> &);
template
void FETools::interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<double> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 BlockVector<double> &);
template
void FETools::back_interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<double> &,
 const FiniteElement<deal_II_dimension> &, BlockVector<double> &);
template
void FETools::back_interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 const BlockVector<double> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 BlockVector<double> &);
template
void FETools::interpolation_difference<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<double> &,
 const FiniteElement<deal_II_dimension> &, BlockVector<double> &);
template
void FETools::interpolation_difference<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 const BlockVector<double> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 BlockVector<double> &);
template
void FETools::project_dg<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<double> &,
 const DoFHandler<deal_II_dimension> &, BlockVector<double> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<double> &,
 const DoFHandler<deal_II_dimension> &, BlockVector<double> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<double> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 BlockVector<double> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<double> &,
 const DoFHandler<deal_II_dimension> &, Vector<double> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<double> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 Vector<double> &);


template
void FETools::interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<float> &,
 const DoFHandler<deal_II_dimension> &, BlockVector<float> &);
template
void FETools::interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<float> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 BlockVector<float> &);
template
void FETools::back_interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<float> &,
 const FiniteElement<deal_II_dimension> &, BlockVector<float> &);
template
void FETools::back_interpolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 const BlockVector<float> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 BlockVector<float> &);
template
void FETools::interpolation_difference<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<float> &,
 const FiniteElement<deal_II_dimension> &, BlockVector<float> &);
template
void FETools::interpolation_difference<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 const BlockVector<float> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 BlockVector<float> &);
template
void FETools::project_dg<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<float> &,
 const DoFHandler<deal_II_dimension> &, BlockVector<float> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<float> &,
 const DoFHandler<deal_II_dimension> &, BlockVector<float> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<float> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 BlockVector<float> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<float> &,
 const DoFHandler<deal_II_dimension> &, Vector<float> &);
template
void FETools::extrapolate<deal_II_dimension>
(const DoFHandler<deal_II_dimension> &, const BlockVector<float> &,
 const DoFHandler<deal_II_dimension> &, const ConstraintMatrix &,
 Vector<float> &);

template
FiniteElement<deal_II_dimension> *
FETools::get_fe_from_name<deal_II_dimension> (const std::string &);

template
void
FETools::
compute_projection_from_quadrature_points_matrix (const FiniteElement<deal_II_dimension> &fe,
                                                  const Quadrature<deal_II_dimension>    &lhs_quadrature,
                                                  const Quadrature<deal_II_dimension>    &rhs_quadrature,
                                                  FullMatrix<double>       &X);

template
void
FETools::
compute_interpolation_to_quadrature_points_matrix (const FiniteElement<deal_II_dimension> &fe,
                                                   const Quadrature<deal_II_dimension>    &quadrature,
                                                   FullMatrix<double>       &I_q);


/*----------------------------   fe_tools.cc     ---------------------------*/
