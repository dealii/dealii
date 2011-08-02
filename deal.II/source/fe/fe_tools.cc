//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/householder.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_tools.h>
#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/fe_dgp_nonparametric.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/hp/dof_handler.h>

#include <deal.II/base/std_cxx1x/shared_ptr.h>

#include <deal.II/base/index_set.h>

#include <iostream>


DEAL_II_NAMESPACE_OPEN

namespace FETools
{
				// Not implemented in the general case.
  template <class FE>
  FiniteElement<FE::dimension, FE::dimension>*
  FEFactory<FE>::get (const Quadrature<1> &) const
  {
    Assert(false, ExcNotImplemented());
    return 0;
  }

				// Specializations for FE_Q.
  template <>
  FiniteElement<1, 1>*
  FEFactory<FE_Q<1, 1> >::get (const Quadrature<1> &quad) const
  {
    return new FE_Q<1>(quad);
  }
  template <>
  FiniteElement<2, 2>*
  FEFactory<FE_Q<2, 2> >::get (const Quadrature<1> &quad) const
  {
    return new FE_Q<2>(quad);
  }
  template <>
  FiniteElement<3, 3>*
  FEFactory<FE_Q<3, 3> >::get (const Quadrature<1> &quad) const
  {
    return new FE_Q<3>(quad);
  }


				// Specializations for FE_DGQArbitraryNodes.
  template <>
  FiniteElement<1, 1>*
  FEFactory<FE_DGQ<1> >::get (const Quadrature<1> &quad) const
  {
    return new FE_DGQArbitraryNodes<1>(quad);
  }
  template <>
  FiniteElement<2, 2>*
  FEFactory<FE_DGQ<2> >::get (const Quadrature<1> &quad) const
  {
    return new FE_DGQArbitraryNodes<2>(quad);
  }
  template <>
  FiniteElement<3, 3>*
  FEFactory<FE_DGQ<3> >::get (const Quadrature<1> &quad) const
  {
    return new FE_DGQArbitraryNodes<3>(quad);
  }
}

namespace
{
                                   // use a shared pointer to factory
                                   // objects, to ensure that they get
                                   // deleted at the end of the
                                   // program run and don't end up as
                                   // apparent memory leaks to
                                   // programs like valgrind.
                                   // a function that returns the
                                   // default set of finite element
                                   // names and factory objects for
                                   // them. used to initialize
                                   // fe_name_map below
  template <int dim>
#ifdef DEAL_II_ANON_NAMESPACE_BUG
    static
#endif
  std::map<std::string,std_cxx1x::shared_ptr<const FETools::FEFactoryBase<dim> > >
  get_default_fe_names ()
  {
    std::map<std::string,
             std_cxx1x::shared_ptr<const FETools::FEFactoryBase<dim> > > default_map;

    typedef std_cxx1x::shared_ptr<const FETools::FEFactoryBase<dim> >
      FEFactoryPointer;

    default_map["FE_Q_Hierarchical"]
      = FEFactoryPointer(new FETools::FEFactory<FE_Q_Hierarchical<dim> >);
    default_map["FE_ABF"]
      = FEFactoryPointer(new FETools::FEFactory<FE_RaviartThomas<dim> >);
    default_map["FE_RaviartThomas"]
      = FEFactoryPointer(new FETools::FEFactory<FE_RaviartThomas<dim> >);
    default_map["FE_RaviartThomasNodal"]
      = FEFactoryPointer(new FETools::FEFactory<FE_RaviartThomasNodal<dim> >);
    default_map["FE_Nedelec"]
      = FEFactoryPointer(new FETools::FEFactory<FE_Nedelec<dim> >);
    default_map["FE_DGPNonparametric"]
      = FEFactoryPointer(new FETools::FEFactory<FE_DGPNonparametric<dim> >);
    default_map["FE_DGP"]
      = FEFactoryPointer(new FETools::FEFactory<FE_DGP<dim> >);
    default_map["FE_DGPMonomial"]
      = FEFactoryPointer(new FETools::FEFactory<FE_DGPMonomial<dim> >);
    default_map["FE_DGQ"]
      = FEFactoryPointer(new FETools::FEFactory<FE_DGQ<dim> >);
    default_map["FE_DGQArbitraryNodes"]
      = FEFactoryPointer(new FETools::FEFactory<FE_DGQ<dim> >);
    default_map["FE_Q"]
      = FEFactoryPointer(new FETools::FEFactory<FE_Q<dim> >);

    return default_map;
  }



                                   // have a lock that guarantees that
                                   // at most one thread is changing
                                   // and accessing the fe_name_map
                                   // variable. make this lock local
                                   // to this file.
				   //
				   // this and the next variable are
				   // declared static (even though
				   // they're in an anonymous
				   // namespace) in order to make icc
				   // happy (which otherwise reports a
				   // multiply defined symbol when
				   // linking libraries for more than
				   // one space dimension together
  static
  Threads::ThreadMutex fe_name_map_lock;

                                   // This is the map used by
                                   // FETools::get_fe_from_name and
                                   // FETools::add_fe_name. Since
                                   // FEFactoryBase has a template
                                   // parameter dim, it could not be a
                                   // member variable of FETools. On
                                   // the other hand, it is only
                                   // accessed by functions in this
                                   // file, so it is safe to make it a
                                   // static variable here. It must be
                                   // static so that we can link
                                   // several dimensions together.
                                   //
                                   // it is initialized at program start time
                                   // using the function above. because at
                                   // this time there are no threads running,
                                   // there are no thread-safety issues
                                   // here. since this is compiled for all
                                   // dimensions at once, need to create
                                   // objects for each dimension and then
                                   // separate between them further down
  static
  std::map<std::string,
	   std_cxx1x::shared_ptr<const FETools::FEFactoryBase<1> > >
  fe_name_map_1d
  = get_default_fe_names<1> ();
  static
  std::map<std::string,
	   std_cxx1x::shared_ptr<const FETools::FEFactoryBase<2> > >
  fe_name_map_2d
  = get_default_fe_names<2> ();
  static
  std::map<std::string,
	   std_cxx1x::shared_ptr<const FETools::FEFactoryBase<3> > >
  fe_name_map_3d
  = get_default_fe_names<3> ();
}






namespace
{

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
  template <int dim, int spacedim>
  inline
  void gim_forwarder (const FiniteElement<dim,spacedim> &fe1,
                      const FiniteElement<dim,spacedim> &fe2,
                      FullMatrix<double> &interpolation_matrix)
  {
    fe2.get_interpolation_matrix (fe1, interpolation_matrix);
  }


  template <int dim, typename number, int spacedim>
  inline
  void gim_forwarder (const FiniteElement<dim,spacedim> &fe1,
                      const FiniteElement<dim,spacedim> &fe2,
                      FullMatrix<number> &interpolation_matrix)
  {
    FullMatrix<double> tmp (interpolation_matrix.m(),
                            interpolation_matrix.n());
    fe2.get_interpolation_matrix (fe1, tmp);
    interpolation_matrix = tmp;
  }



				   // return how many characters
				   // starting at the given position
				   // of the string match either the
				   // generic string "<dim>" or the
				   // specialized string with "dim"
				   // replaced with the numeric value
				   // of the template argument
  template <int dim, int spacedim>
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


namespace FETools
{
  template <int dim, int spacedim>
  FEFactoryBase<dim,spacedim>::~FEFactoryBase()
  {}


  template<int dim, int spacedim>
  void compute_component_wise(
    const FiniteElement<dim,spacedim>& element,
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



  template<int dim, int spacedim>
  void compute_block_renumbering (
    const FiniteElement<dim,spacedim>& element,
    std::vector<unsigned int>& renumbering,
    std::vector<unsigned int>& block_data,
    bool return_start_indices)
  {
    Assert(renumbering.size() == element.dofs_per_cell,
	   ExcDimensionMismatch(renumbering.size(),
				element.dofs_per_cell));
    Assert(block_data.size() == element.n_blocks(),
	   ExcDimensionMismatch(block_data.size(),
				element.n_blocks()));

    unsigned int k=0;
    unsigned int i=0;
    for (unsigned int b=0;b<element.n_base_elements();++b)
      for (unsigned int m=0;m<element.element_multiplicity(b);++m)
	{
	  block_data[i++] = (return_start_indices)
			    ? k
			    : (element.base_element(b).n_dofs_per_cell());
	  k += element.base_element(b).n_dofs_per_cell();
	}
    Assert (i == element.n_blocks(), ExcInternalError());

    std::vector<unsigned int> start_indices(block_data.size());
    k = 0;
    for (unsigned int i=0;i<block_data.size();++i)
      if (return_start_indices)
	start_indices[i] = block_data[i];
      else
	{
	  start_indices[i] = k;
	  k += block_data[i];
	}

//TODO:[GK] This does not work for a single RT
    for (unsigned int i=0;i<element.dofs_per_cell;++i)
      {
	std::pair<unsigned int, unsigned int>
	  indices = element.system_to_block_index(i);
	renumbering[i] = start_indices[indices.first]
			 +indices.second;
      }
  }



  template <int dim, typename number, int spacedim>
  void get_interpolation_matrix (const FiniteElement<dim,spacedim> &fe1,
				 const FiniteElement<dim,spacedim> &fe2,
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
    catch (typename FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented &)
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
    Assert (fe1.is_primitive() == true, ExcFENotPrimitive());
    Assert (fe2.is_primitive() == true, ExcFENotPrimitive());

				     // Initialize FEValues for fe1 at
				     // the unit support points of the
				     // fe2 element.
    const std::vector<Point<dim> > &
      fe2_support_points = fe2.get_unit_support_points ();

    typedef FiniteElement<dim,spacedim> FEL;
    Assert(fe2_support_points.size()==fe2.dofs_per_cell,
	   typename FEL::ExcFEHasNoSupportPoints());

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



  template <int dim, typename number, int spacedim>
  void get_back_interpolation_matrix(const FiniteElement<dim,spacedim> &fe1,
				     const FiniteElement<dim,spacedim> &fe2,
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



  template <int dim, typename number, int spacedim>
  void get_interpolation_difference_matrix (const FiniteElement<dim,spacedim> &fe1,
					    const FiniteElement<dim,spacedim> &fe2,
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



  template <int dim, typename number, int spacedim>
  void get_projection_matrix (const FiniteElement<dim,spacedim> &fe1,
			      const FiniteElement<dim,spacedim> &fe2,
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
    Triangulation<dim,spacedim> tr;
    GridGenerator::hyper_cube(tr);

				     // Choose a quadrature rule
				     // Gauss is exact up to degree 2n-1
    const unsigned int degree = std::max(fe1.tensor_degree(), fe2.tensor_degree());
    Assert (degree != numbers::invalid_unsigned_int,
	    ExcNotImplemented());

    QGauss<dim> quadrature(degree+1);
				     // Set up FEValues.
    const UpdateFlags flags = update_values | update_quadrature_points | update_JxW_values;
    FEValues<dim> val1 (fe1, quadrature, update_values);
    val1.reinit (tr.begin_active());
    FEValues<dim> val2 (fe2, quadrature, flags);
    val2.reinit (tr.begin_active());

				     // Integrate and invert mass matrix
				     // This happens in the target space
    FullMatrix<double> mass (n2, n2);

    for (unsigned int k=0;k<quadrature.size();++k)
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
	  for (unsigned int k=0;k<quadrature.size();++k)
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


  template<int dim, int spacedim>
  void
  compute_node_matrix(
    FullMatrix<double>& N,
    const FiniteElement<dim,spacedim>& fe)
  {
    const unsigned int n_dofs = fe.dofs_per_cell;
    Assert (fe.has_generalized_support_points(), ExcNotInitialized());
    Assert (N.n()==n_dofs, ExcDimensionMismatch(N.n(), n_dofs));
    Assert (N.m()==n_dofs, ExcDimensionMismatch(N.m(), n_dofs));

    const std::vector<Point<dim> >& points = fe.get_generalized_support_points();

				     // We need the values of the
				     // polynomials in all generalized
				     // support points.
    std::vector<std::vector<double> >
      values (dim, std::vector<double>(points.size()));

				     // In this vector, we store the
				     // result of the interpolation
    std::vector<double> local_dofs(n_dofs);

				     // One row per shape
				     // function. Remember that these
				     // are the 'raw' shape functions
				     // where the inverse node matrix is
				     // empty. Otherwise, this would
				     // yield identity.
    for (unsigned int i=0;i<n_dofs;++i)
      {
	for (unsigned int k=0;k<values[0].size();++k)
	  for (unsigned int d=0;d<dim;++d)
	    values[d][k] = fe.shape_value_component(i,points[k],d);
	fe.interpolate(local_dofs, values);
					 // Enter the interpolated dofs
					 // into the matrix
	for (unsigned int j=0;j<n_dofs;++j)
	  N(j,i) = local_dofs[j];
      }
  }



  template<>
  void
  compute_embedding_matrices(const FiniteElement<1,2> &,
			     std::vector<std::vector<FullMatrix<double> > > &,
			     const bool)
  {
    Assert(false, ExcNotImplemented());
  }



  template<>
  void
  compute_embedding_matrices(const FiniteElement<2,3>&,
			     std::vector<std::vector<FullMatrix<double> > >&,
			     const bool)
  {
    Assert(false, ExcNotImplemented());
  }



  namespace {
    template<int dim, typename number, int spacedim>
    void
    compute_embedding_for_shape_function (
      const unsigned int i,
      const FiniteElement<dim, spacedim>& fe,
      const FEValues<dim>& coarse,
      const Householder<double>& H,
      FullMatrix<number>& this_matrix)
    {
      const unsigned int n  = fe.dofs_per_cell;
      const unsigned int nd = fe.n_components ();
      const unsigned int nq = coarse.n_quadrature_points;

      Vector<number> v_coarse(nq*nd);
      Vector<number> v_fine(n);

				       // The right hand side of
				       // the least squares
				       // problem consists of the
				       // function values of the
				       // coarse grid function in
				       // each quadrature point.
      if (fe.is_primitive ())
	{
	  const unsigned int
	    d = fe.system_to_component_index (i).first;
	  const double* phi_i = &coarse.shape_value (i, 0);

	  for (unsigned int k = 0; k < nq; ++k)
	    v_coarse (k * nd + d) = phi_i[k];
	}

      else
	for (unsigned int d = 0; d < nd; ++d)
	  for (unsigned int k = 0; k < nq; ++k)
	    v_coarse (k * nd + d) = coarse.shape_value_component (i, k, d);

				       // solve the least squares
				       // problem.
      const double result = H.least_squares (v_fine, v_coarse);
      Assert (result < 1.e-12, ExcLeastSquaresError (result));

				       // Copy into the result
				       // matrix. Since the matrix
				       // maps a coarse grid
				       // function to a fine grid
				       // function, the columns
				       // are fine grid.
      for (unsigned int j = 0; j < n; ++j)
	this_matrix(j, i) = v_fine(j);
    }


    template<int dim, typename number, int spacedim>
    void
    compute_embedding_matrices_for_refinement_case (
      const FiniteElement<dim, spacedim>& fe,
      std::vector<FullMatrix<number> >& matrices,
      const unsigned int ref_case)
    {
      const unsigned int n  = fe.dofs_per_cell;
      const unsigned int nc = GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));
      for (unsigned int i = 0; i < nc; ++i)
	{
	  Assert(matrices[i].n() == n, ExcDimensionMismatch(matrices[i].n (), n));
	  Assert(matrices[i].m() == n, ExcDimensionMismatch(matrices[i].m (), n));
	}

				       // Set up meshes, one with a single
				       // reference cell and refine it once
      Triangulation<dim,spacedim> tria;
      GridGenerator::hyper_cube (tria, 0, 1);
      tria.begin_active()->set_refine_flag (RefinementCase<dim>(ref_case));
      tria.execute_coarsening_and_refinement ();

      MappingCartesian<dim> mapping;
      const unsigned int degree = fe.degree;
      QGauss<dim> q_fine (degree+1);
      const unsigned int nq = q_fine.size();

      FEValues<dim> fine (mapping, fe, q_fine,
			  update_quadrature_points |
			  update_JxW_values |
			  update_values);

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
      fine.reinit (tria.begin_active ());
      const unsigned int nd = fe.n_components ();
      FullMatrix<number> A (nq*nd, n);

      for (unsigned int j = 0; j < n; ++j)
	for (unsigned int d = 0; d < nd; ++d)
	  for (unsigned int k = 0; k < nq; ++k)
	    A (k * nd + d, j) = fine.shape_value_component (j, k, d);

      Householder<double> H (A);
      unsigned int cell_number = 0;

      Threads::TaskGroup<void> task_group;

      for (typename Triangulation<dim>::active_cell_iterator
	     fine_cell = tria.begin_active (); fine_cell != tria.end ();
	   ++fine_cell, ++cell_number)
	{
	  fine.reinit (fine_cell);

					   // evaluate on the coarse cell (which
					   // is the first -- inactive -- cell on
					   // the lowest level of the
					   // triangulation we have created)
	  const Quadrature<dim> q_coarse (fine.get_quadrature_points (),
					  fine.get_JxW_values ());
	  FEValues<dim> coarse (mapping, fe, q_coarse, update_values);

	  coarse.reinit (tria.begin (0));

	  FullMatrix<double> &this_matrix = matrices[cell_number];

					   // Compute this once for each
					   // coarse grid basis function. can
					   // spawn subtasks if n is
					   // sufficiently large so that there
					   // are more than about 5000
					   // operations in the inner loop
					   // (which is basically const * n^2
					   // operations).
	  if (n > 30)
	    {
	      for (unsigned int i = 0; i < n; ++i)
		{
		  task_group +=
		    Threads::new_task (&compute_embedding_for_shape_function<dim, number, spacedim>,
				       i, fe, coarse, H, this_matrix);
		}
	      task_group.join_all();
	    }
	  else
	    {
	      for (unsigned int i = 0; i < n; ++i)
		{
		  compute_embedding_for_shape_function<dim, number, spacedim>
		    (i, fe, coarse, H, this_matrix);
		}
	    }

					   // Remove small entries from
					   // the matrix
	  for (unsigned int i = 0; i < this_matrix.m (); ++i)
	    for (unsigned int j = 0; j < this_matrix.n (); ++j)
	      if (std::fabs (this_matrix (i, j)) < 1e-12)
		this_matrix (i, j) = 0.;
	}

      Assert (cell_number == GeometryInfo<dim>::n_children (RefinementCase<dim> (ref_case)),
	      ExcInternalError ());
    }
  }


  template <int dim, typename number, int spacedim>
  void
  compute_embedding_matrices(const FiniteElement<dim,spacedim>& fe,
			     std::vector<std::vector<FullMatrix<number> > >& matrices,
			     const bool isotropic_only)
  {
    Threads::TaskGroup<void> task_group;

				     // loop over all possible refinement cases
    unsigned int ref_case = (isotropic_only)
			    ? RefinementCase<dim>::isotropic_refinement
			    : RefinementCase<dim>::cut_x;

    for (;ref_case <= RefinementCase<dim>::isotropic_refinement; ++ref_case)
      task_group += Threads::new_task (&compute_embedding_matrices_for_refinement_case<dim, number, spacedim>,
				       fe, matrices[ref_case-1], ref_case);

    task_group.join_all ();
  }



  template <int dim, typename number, int spacedim>
  void
  compute_face_embedding_matrices(const FiniteElement<dim,spacedim>& fe,
				  FullMatrix<number> (&matrices)[GeometryInfo<dim>::max_children_per_face],
				  const unsigned int face_coarse,
				  const unsigned int face_fine)
  {
    Assert(face_coarse==0, ExcNotImplemented());
    Assert(face_fine==0, ExcNotImplemented());

    const unsigned int nc = GeometryInfo<dim>::max_children_per_face;
    const unsigned int n  = fe.dofs_per_face;
    const unsigned int nd = fe.n_components();
    const unsigned int degree = fe.degree;

    const bool normal = fe.conforms(FiniteElementData<dim>::Hdiv);
    const bool tangential = fe.conforms(FiniteElementData<dim>::Hcurl);

    for (unsigned int i=0;i<nc;++i)
      {
	Assert(matrices[i].n() == n, ExcDimensionMismatch(matrices[i].n(),n));
	Assert(matrices[i].m() == n, ExcDimensionMismatch(matrices[i].m(),n));
      }

				     // In order to make the loops below
				     // simpler, we introduce vectors
				     // containing for indices 0-n the
				     // number of the corresponding
				     // shape value on the cell.
    std::vector<unsigned int> face_c_dofs(n);
    std::vector<unsigned int> face_f_dofs(n);
    unsigned int k=0;
    for (unsigned int i=0;i<GeometryInfo<dim>::vertices_per_face;++i)
      {
	const unsigned int offset_c = GeometryInfo<dim>::face_to_cell_vertices(face_coarse, i)
				      *fe.dofs_per_vertex;
	const unsigned int offset_f = GeometryInfo<dim>::face_to_cell_vertices(face_fine, i)
				      *fe.dofs_per_vertex;
	for (unsigned int j=0;j<fe.dofs_per_vertex;++j)
	  {
	    face_c_dofs[k] = offset_c + j;
	    face_f_dofs[k] = offset_f + j;
	    ++k;
	  }
      }
    for (unsigned int i=1;i<=GeometryInfo<dim>::lines_per_face;++i)
      {
	const unsigned int offset_c = fe.first_line_index
				      + GeometryInfo<dim>::face_to_cell_lines(face_coarse, i-1)
				      *fe.dofs_per_line;
	const unsigned int offset_f = fe.first_line_index
				      + GeometryInfo<dim>::face_to_cell_lines(face_fine, i-1)
				      *fe.dofs_per_line;
	for (unsigned int j=0;j<fe.dofs_per_line;++j)
	  {
	    face_c_dofs[k] = offset_c + j;
	    face_f_dofs[k] = offset_f + j;
	    ++k;
	  }
      }
    for (unsigned int i=1;i<=GeometryInfo<dim>::quads_per_face;++i)
      {
	const unsigned int offset_c = fe.first_quad_index
				      + face_coarse
				      *fe.dofs_per_quad;
	const unsigned int offset_f = fe.first_quad_index
				      + face_fine
				      *fe.dofs_per_quad;
	for (unsigned int j=0;j<fe.dofs_per_quad;++j)
	  {
	    face_c_dofs[k] = offset_c + j;
	    face_f_dofs[k] = offset_f + j;
	    ++k;
	  }
      }
    Assert (k == fe.dofs_per_face, ExcInternalError());

				     // Set up meshes, one with a single
				     // reference cell and refine it once
    Triangulation<dim,spacedim> tria;
    GridGenerator::hyper_cube (tria, 0, 1);
    tria.refine_global(1);
    MappingCartesian<dim> mapping;

				     // Setup quadrature and FEValues
				     // for a face. We cannot use
				     // FEFaceValues and
				     // FESubfaceValues because of
				     // some nifty handling of
				     // refinement cases. Guido stops
				     // disliking and instead starts
				     // hating the anisotropic implementation
    QGauss<dim-1> q_gauss(degree+1);
    const Quadrature<dim> q_fine = QProjector<dim>::project_to_face(q_gauss, face_fine);
    const unsigned int nq = q_fine.size();

    FEValues<dim> fine (mapping, fe, q_fine,
			update_quadrature_points | update_JxW_values | update_values);

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
    for (unsigned int j=0;j<n;++j)
      for (unsigned int k=0;k<nq;++k)
	if (nd != dim)
	  for (unsigned int d=0;d<nd;++d)
	    A(k*nd+d,j) = fine.shape_value_component(face_f_dofs[j],k,d);
	else
	  {
	    if (normal)
	      A(k*nd,j) = fine.shape_value_component(face_f_dofs[j],k,0);
	    if (tangential)
	      for (unsigned int d=1;d<dim;++d)
		A(k*nd+d,j) = fine.shape_value_component(face_f_dofs[j],k,d);
	  }

    Householder<double> H(A);

    Vector<number> v_coarse(nq*nd);
    Vector<number> v_fine(n);



    for (unsigned int cell_number = 0; cell_number < GeometryInfo<dim>::max_children_per_face;
	 ++cell_number)
      {
	const Quadrature<dim> q_coarse
	  = QProjector<dim>::project_to_subface(q_gauss, face_coarse, cell_number);
	FEValues<dim> coarse (mapping, fe, q_coarse, update_values);

	typename Triangulation<dim,spacedim>::active_cell_iterator fine_cell
	  = tria.begin(0)->child(GeometryInfo<dim>::child_cell_on_face(
				   tria.begin(0)->refinement_case(), face_coarse, cell_number));
	fine.reinit(fine_cell);
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
	    for (unsigned int k=0;k<nq;++k)
	      if (nd != dim)
		for (unsigned int d=0;d<nd;++d)
		  v_coarse(k*nd+d) = coarse.shape_value_component (face_c_dofs[i],k,d);
	      else
	  {
	    if (normal)
	      v_coarse(k*nd) = coarse.shape_value_component(face_c_dofs[i],k,0);
	    if (tangential)
	      for (unsigned int d=1;d<dim;++d)
		v_coarse(k*nd+d) = coarse.shape_value_component(face_c_dofs[i],k,d);
	  }
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
  }



  template <>
  void
  compute_projection_matrices(const FiniteElement<1,2>&,
			      std::vector<std::vector<FullMatrix<double> > >&, bool)
  {
    Assert(false, ExcNotImplemented());
  }



  template <>
  void
  compute_projection_matrices(const FiniteElement<2,3>&,
			      std::vector<std::vector<FullMatrix<double> > >&, bool)
  {
    Assert(false, ExcNotImplemented());
  }



  template <int dim, typename number, int spacedim>
  void
  compute_projection_matrices(const FiniteElement<dim,spacedim>& fe,
			      std::vector<std::vector<FullMatrix<number> > >& matrices,
			      const bool isotropic_only)
  {
    const unsigned int n  = fe.dofs_per_cell;
    const unsigned int nd = fe.n_components();
    const unsigned int degree = fe.degree;

				     // prepare FEValues, quadrature etc on
				     // coarse cell
    MappingCartesian<dim> mapping;
    QGauss<dim> q_fine(degree+1);
    const unsigned int nq = q_fine.size();

				     // create mass matrix on coarse cell.
    FullMatrix<number> mass(n, n);
    {
				       // set up a triangulation for coarse cell
      Triangulation<dim,spacedim> tr;
      GridGenerator::hyper_cube (tr, 0, 1);

      FEValues<dim> coarse (mapping, fe, q_fine,
			    update_JxW_values | update_values);

      typename Triangulation<dim,spacedim>::cell_iterator coarse_cell
	= tr.begin(0);
      coarse.reinit (coarse_cell);

      const std::vector<double> & JxW = coarse.get_JxW_values();
      for (unsigned int i=0;i<n;++i)
	for (unsigned int j=0;j<n;++j)
	  if (fe.is_primitive())
	    {
	      const double * coarse_i = &coarse.shape_value(i,0);
	      const double * coarse_j = &coarse.shape_value(j,0);
	      double mass_ij = 0;
	      for (unsigned int k=0;k<nq;++k)
		mass_ij += JxW[k] * coarse_i[k] * coarse_j[k];
	      mass(i,j) = mass_ij;
	    }
	  else
	    {
	      double mass_ij = 0;
	      for (unsigned int d=0;d<nd;++d)
		for (unsigned int k=0;k<nq;++k)
		  mass_ij += JxW[k] * coarse.shape_value_component(i,k,d)
			     * coarse.shape_value_component(j,k,d);
	      mass(i,j) = mass_ij;
	    }

				       // invert mass matrix
      mass.gauss_jordan();
    }

				     // loop over all possible
				     // refinement cases
    unsigned int ref_case = (isotropic_only)
			    ? RefinementCase<dim>::isotropic_refinement
			    : RefinementCase<dim>::cut_x;
    for (;ref_case <= RefinementCase<dim>::isotropic_refinement; ++ref_case)
      {
	const unsigned int
	  nc = GeometryInfo<dim>::n_children(RefinementCase<dim>(ref_case));

	for (unsigned int i=0;i<nc;++i)
	  {
	    Assert(matrices[ref_case-1][i].n() == n,
		   ExcDimensionMismatch(matrices[ref_case-1][i].n(),n));
	    Assert(matrices[ref_case-1][i].m() == n,
		   ExcDimensionMismatch(matrices[ref_case-1][i].m(),n));
	  }

					 // create a respective refinement on the
					 // triangulation
	Triangulation<dim,spacedim> tr;
	GridGenerator::hyper_cube (tr, 0, 1);
	tr.begin_active()->set_refine_flag(RefinementCase<dim>(ref_case));
	tr.execute_coarsening_and_refinement();

	FEValues<dim> fine (mapping, fe, q_fine,
			    update_quadrature_points | update_JxW_values |
			    update_values);

	typename Triangulation<dim,spacedim>::cell_iterator coarse_cell
	  = tr.begin(0);

	Vector<number> v_coarse(n);
	Vector<number> v_fine(n);

	for (unsigned int cell_number=0;cell_number<nc;++cell_number)
	  {
	    FullMatrix<double> &this_matrix = matrices[ref_case-1][cell_number];

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

	    const std::vector<double> & JxW = fine.get_JxW_values();

					     // Outer loop over all fine
					     // grid shape functions phi_j
	    for (unsigned int j=0;j<fe.dofs_per_cell;++j)
	      {
		for (unsigned int i=0; i<fe.dofs_per_cell;++i)
		  {
		    if (fe.is_primitive())
		      {
			const double * coarse_i = &coarse.shape_value(i,0);
			const double * fine_j = &fine.shape_value(j,0);

			double update = 0;
			for (unsigned int k=0; k<nq; ++k)
			  update += JxW[k] * coarse_i[k] * fine_j[k];
			v_fine(i) = update;
		      }
		    else
		      {
			double update = 0;
			for (unsigned int d=0; d<nd; ++d)
			  for (unsigned int k=0; k<nq; ++k)
			    update += JxW[k] * coarse.shape_value_component(i,k,d)
				      * fine.shape_value_component(j,k,d);
			v_fine(i) = update;
		      }
		  }

						 // RHS ready. Solve system
						 // and enter row into
						 // matrix
		mass.vmult (v_coarse, v_fine);
		for (unsigned int i=0;i<fe.dofs_per_cell;++i)
		  this_matrix(i,j) = v_coarse(i);
	      }

					     // Remove small entries from
					     // the matrix
	    for (unsigned int i=0; i<this_matrix.m(); ++i)
	      for (unsigned int j=0; j<this_matrix.n(); ++j)
		if (std::fabs(this_matrix(i,j)) < 1e-12)
		  this_matrix(i,j) = 0.;
	  }
      }
  }


  template <int dim, int spacedim,
	    template <int, int> class DH1,
	    template <int, int> class DH2,
	    class InVector, class OutVector>
  void
  interpolate(const DH1<dim, spacedim> &dof1,
	      const InVector           &u1,
	      const DH2<dim, spacedim> &dof2,
	      OutVector                &u2)
	    {
	      ConstraintMatrix dummy;
	      dummy.close();
	      interpolate(dof1, u1, dof2, dummy, u2);
	    }



  template <int dim, int spacedim,
            template <int, int> class DH1,
            template <int, int> class DH2,
            class InVector, class OutVector>
  void
  interpolate (const DH1<dim, spacedim> &dof1,
	       const InVector           &u1,
	       const DH2<dim, spacedim> &dof2,
	       const ConstraintMatrix   &constraints,
	       OutVector                &u2)
  {
    Assert(&dof1.get_tria()==&dof2.get_tria(), ExcTriangulationMismatch());

    Assert(u1.size()==dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size()==dof2.n_dofs(),
           ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

#ifdef DEAL_II_USE_PETSC
    if (dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
                                   // if u1 is a parallel distributed
                                   // PETSc vector, we check the local
                                   // size of u1 for safety
          Assert(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size() == dof1.locally_owned_dofs().n_elements(),
                 ExcDimensionMismatch(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size(), dof1.locally_owned_dofs().n_elements()));
        };

    if (dynamic_cast<PETScWrappers::MPI::Vector *>(&u2) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof2) != 0)
        {
          Assert(dynamic_cast<PETScWrappers::MPI::Vector *>(&u2)->local_size() == dof2.locally_owned_dofs().n_elements(),
                 ExcDimensionMismatch(dynamic_cast<PETScWrappers::MPI::Vector *>(&u2)->local_size(), dof2.locally_owned_dofs().n_elements()));
        };
#endif
                                   // allocate vectors at maximal
                                   // size. will be reinited in inner
                                   // cell, but Vector makes sure that
                                   // this does not lead to
                                   // reallocation of memory
    Vector<typename OutVector::value_type> u1_local(DoFTools::max_dofs_per_cell(dof1));
    Vector<typename OutVector::value_type> u2_local(DoFTools::max_dofs_per_cell(dof2));

                                   // have a map for interpolation
                                   // matrices. shared_ptr make sure
                                   // that memory is released again
    std::map<const FiniteElement<dim,spacedim> *,
      std::map<const FiniteElement<dim,spacedim> *,
      std_cxx1x::shared_ptr<FullMatrix<double> > > >
      interpolation_matrices;

    typename DH1<dim,spacedim>::active_cell_iterator cell1 = dof1.begin_active(),
                                                     endc1 = dof1.end();
    typename DH2<dim,spacedim>::active_cell_iterator cell2 = dof2.begin_active(),
                                                     endc2 = dof2.end();

    std::vector<unsigned int> dofs;
    dofs.reserve (DoFTools::max_dofs_per_cell (dof2));

    u2 = 0;
    OutVector touch_count(u2);
    touch_count = 0;

                                   // for distributed triangulations,
                                   // we can only interpolate u1 on
                                   // a cell, which this processor owns,
				   // so we have to know the subdomain_id
    const types::subdomain_id_t subdomain_id =
      dof1.get_tria().locally_owned_subdomain();

    for (; cell1!=endc1; ++cell1, ++cell2)
      if ((cell1->subdomain_id() == subdomain_id)
          ||
          (subdomain_id == types::invalid_subdomain_id))
        {
          Assert(cell1->get_fe().n_components() == cell2->get_fe().n_components(),
                 ExcDimensionMismatch (cell1->get_fe().n_components(),
                                       cell2->get_fe().n_components()));

                                       // for continuous elements on
                                       // grids with hanging nodes we
                                       // need hanging node
                                       // constraints. Consequentely,
                                       // if there are no constraints
                                       // then hanging nodes are not
                                       // allowed.
          const bool hanging_nodes_not_allowed
            = ((cell2->get_fe().dofs_per_vertex != 0) &&
               (constraints.n_constraints() == 0));

          if (hanging_nodes_not_allowed)
	    for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
	      Assert (cell1->at_boundary(face) ||
		      cell1->neighbor(face)->level() == cell1->level(),
		      ExcHangingNodesNotAllowed(0));


          const unsigned int dofs_per_cell1 = cell1->get_fe().dofs_per_cell;
          const unsigned int dofs_per_cell2 = cell2->get_fe().dofs_per_cell;
          u1_local.reinit (dofs_per_cell1);
          u2_local.reinit (dofs_per_cell2);

                                       // check if interpolation
                                       // matrix for this particular
                                       // pair of elements is already
                                       // there
          if (interpolation_matrices[&cell1->get_fe()][&cell2->get_fe()].get() == 0)
            {
              std_cxx1x::shared_ptr<FullMatrix<double> >
                interpolation_matrix (new FullMatrix<double> (dofs_per_cell2,
                                                              dofs_per_cell1));
              interpolation_matrices[&cell1->get_fe()][&cell2->get_fe()]
                = interpolation_matrix;

              get_interpolation_matrix(cell1->get_fe(),
				       cell2->get_fe(),
				       *interpolation_matrix);
            }

          cell1->get_dof_values(u1, u1_local);
          interpolation_matrices[&cell1->get_fe()][&cell2->get_fe()]
            ->vmult(u2_local, u1_local);

          dofs.resize (dofs_per_cell2);
          cell2->get_dof_indices(dofs);

          for (unsigned int i=0; i<dofs_per_cell2; ++i)
            {
              u2(dofs[i])+=u2_local(i);
              touch_count(dofs[i]) += 1;
            }
        }
                                   // cell1 is at the end, so should
                                   // be cell2
    Assert (cell2 == endc2, ExcInternalError());

    u2.compress();
    touch_count.compress();

                                   // if we work on parallel distributed
                                   // vectors, we have to ensure, that we only
                                   // work on dofs this processor owns.
    IndexSet  locally_owned_dofs = dof2.locally_owned_dofs();

				   // when a discontinuous element is
				   // interpolated to a continuous
				   // one, we take the mean values.
				   // for parallel vectors check,
				   // if this component is owned by
				   // this processor.
    for (unsigned int i=0; i<dof2.n_dofs(); ++i)
      if (locally_owned_dofs.is_element(i))
        {
          Assert(touch_count(i)!=0, ExcInternalError());
          u2(i) /= touch_count(i);
        }

				   // finish the work on parallel vectors
    u2.compress();
				   // Apply hanging node constraints.
    constraints.distribute(u2);

				   // and finally update ghost values
    u2.update_ghost_values();
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void
  back_interpolate(const DoFHandler<dim,spacedim>    &dof1,
		   const InVector           &u1,
		   const FiniteElement<dim,spacedim> &fe2,
		   OutVector                &u1_interpolated)
  {
    Assert(dof1.get_fe().n_components() == fe2.n_components(),
	   ExcDimensionMismatch(dof1.get_fe().n_components(), fe2.n_components()));
    Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u1_interpolated.size()==dof1.n_dofs(),
	   ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));

#ifdef DEAL_II_USE_PETSC
    if (dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
                                   // if u1 is a parallel distributed
                                   // PETSc vector, we check the local
                                   // size of u1 for safety
          Assert(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size() == dof1.locally_owned_dofs().n_elements(),
                 ExcDimensionMismatch(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size(), dof1.locally_owned_dofs().n_elements()));
        };

    if (dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
          Assert(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated)->local_size() == dof1.locally_owned_dofs().n_elements(),
                 ExcDimensionMismatch(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated)->local_size(), dof1.locally_owned_dofs().n_elements()));
        };
#endif

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

    const types::subdomain_id_t subdomain_id =
      dof1.get_tria().locally_owned_subdomain();

    typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof1.begin_active(),
							    endc = dof1.end();

    FullMatrix<double> interpolation_matrix(dofs_per_cell1, dofs_per_cell1);
    get_back_interpolation_matrix(dof1.get_fe(), fe2,
				  interpolation_matrix);
    for (; cell!=endc; ++cell)
      if ((cell->subdomain_id() == subdomain_id)
          ||
          (subdomain_id == types::invalid_subdomain_id))
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

				   // if we work on a parallel PETSc vector
				   // we have to finish the work and
				   // update ghost values
    u1_interpolated.compress();
    u1_interpolated.update_ghost_values();
  }



  template <int dim,
	    template <int> class DH,
	    class InVector, class OutVector, int spacedim>
  void
  back_interpolate(const DH<dim>            &dof1,
		   const InVector           &u1,
		   const FiniteElement<dim,spacedim> &fe2,
		   OutVector                &u1_interpolated)
  {
    Assert(u1.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u1_interpolated.size() == dof1.n_dofs(),
           ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));

#ifdef DEAL_II_USE_PETSC
    if (dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
                                   // if u1 is a parallel distributed
                                   // PETSc vector, we check the local
                                   // size of u1 for safety
          Assert(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size() == dof1.locally_owned_dofs().n_elements(),
                 ExcDimensionMismatch(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size(), dof1.locally_owned_dofs().n_elements()));
        };

    if (dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
          Assert(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated)->local_size() == dof1.locally_owned_dofs().n_elements(),
                 ExcDimensionMismatch(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_interpolated)->local_size(), dof1.locally_owned_dofs().n_elements()));
        };
#endif

    Vector<typename OutVector::value_type> u1_local(DoFTools::max_dofs_per_cell(dof1));
    Vector<typename OutVector::value_type> u1_int_local(DoFTools::max_dofs_per_cell(dof1));

    const types::subdomain_id_t subdomain_id =
      dof1.get_tria().locally_owned_subdomain();

    typename DH<dim>::active_cell_iterator cell = dof1.begin_active(),
                                           endc = dof1.end();

					       // map from possible fe objects in
					       // dof1 to the back_interpolation
					       // matrices
    std::map<const FiniteElement<dim> *,
       std_cxx1x::shared_ptr<FullMatrix<double> > > interpolation_matrices;

    for (; cell!=endc; ++cell)
      if ((cell->subdomain_id() == subdomain_id)
          ||
          (subdomain_id == types::invalid_subdomain_id))
        {
          Assert(cell->get_fe().n_components() == fe2.n_components(),
                 ExcDimensionMismatch(cell->get_fe().n_components(),
                 fe2.n_components()));

				     // For continuous elements on
				     // grids with hanging nodes we
				     // need hanging node
				     // constraints. Consequently,
				     // when the elements are
				     // continuous no hanging node
				     // constraints are allowed.
          const bool hanging_nodes_not_allowed=
            (cell->get_fe().dofs_per_vertex != 0) || (fe2.dofs_per_vertex != 0);

          if (hanging_nodes_not_allowed)
            for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
              Assert (cell->at_boundary(face) ||
                      cell->neighbor(face)->level() == cell->level(),
                      ExcHangingNodesNotAllowed(0));

          const unsigned int dofs_per_cell1 = cell->get_fe().dofs_per_cell;

				     // make sure back_interpolation
				     // matrix is available
          if (interpolation_matrices[&cell->get_fe()] != 0)
            {
              interpolation_matrices[&cell->get_fe()] =
                std_cxx1x::shared_ptr<FullMatrix<double> >
                  (new FullMatrix<double>(dofs_per_cell1, dofs_per_cell1));
              get_back_interpolation_matrix(dof1.get_fe(), fe2,
                                            *interpolation_matrices[&cell->get_fe()]);
            }

          u1_local.reinit (dofs_per_cell1);
          u1_int_local.reinit (dofs_per_cell1);

          cell->get_dof_values(u1, u1_local);
          interpolation_matrices[&cell->get_fe()]->vmult(u1_int_local, u1_local);
          cell->set_dof_values(u1_int_local, u1_interpolated);
        };

				   // if we work on a parallel PETSc vector
				   // we have to finish the work and
				   // update ghost values
    u1_interpolated.compress();
    u1_interpolated.update_ghost_values();
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void back_interpolate(const DoFHandler<dim,spacedim> &dof1,
			const ConstraintMatrix &constraints1,
			const InVector &u1,
			const DoFHandler<dim,spacedim> &dof2,
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
#ifdef DEAL_II_USE_PETSC
        if (dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1) != 0)
          {
            AssertThrow (dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1_interpolated) != 0,
                         ExcMessage("Interpolation can only be done between vectors of same type!"));

					 // if u1 is a parallel distributed
					 // PETSc vector, we create a vector u2
					 // with based on the sets of locally owned
					 // and relevant dofs of dof2
            IndexSet  dof2_locally_owned_dofs = dof2.locally_owned_dofs();
            IndexSet  dof2_locally_relevant_dofs;
            DoFTools::extract_locally_relevant_dofs (dof2,
                                                     dof2_locally_relevant_dofs);

            PETScWrappers::MPI::Vector  u2 (dynamic_cast<const PETScWrappers::MPI::Vector *> (&u1)->get_mpi_communicator(),
                                            dof2_locally_owned_dofs,
                                            dof2_locally_relevant_dofs);
            interpolate(dof1, u1, dof2, constraints2, u2);
            interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
          }
        else
#endif
          {
            Vector<typename OutVector::value_type> u2(dof2.n_dofs());
            interpolate(dof1, u1, dof2, constraints2, u2);
            interpolate(dof2, u2, dof1, constraints1, u1_interpolated);
          }
      }
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void interpolation_difference (const DoFHandler<dim,spacedim> &dof1,
				 const InVector &u1,
				 const FiniteElement<dim,spacedim> &fe2,
				 OutVector &u1_difference)
  {
    Assert(dof1.get_fe().n_components() == fe2.n_components(),
	   ExcDimensionMismatch(dof1.get_fe().n_components(), fe2.n_components()));
    Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u1_difference.size()==dof1.n_dofs(),
	   ExcDimensionMismatch(u1_difference.size(), dof1.n_dofs()));

#ifdef DEAL_II_USE_PETSC
    if (dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
                                   // if u1 is a parallel distributed
                                   // PETSc vector, we check the local
                                   // size of u1 for safety
          Assert(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size() == dof1.locally_owned_dofs().n_elements(),
                 ExcDimensionMismatch(dynamic_cast<const PETScWrappers::MPI::Vector *>(&u1)->local_size(), dof1.locally_owned_dofs().n_elements()));
        };

    if (dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_difference) != 0)
      if (dynamic_cast<const DoFHandler<dim>*>(&dof1) != 0)
        {
          Assert(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_difference)->local_size() == dof1.locally_owned_dofs().n_elements(),
                 ExcDimensionMismatch(dynamic_cast<PETScWrappers::MPI::Vector *>(&u1_difference)->local_size(), dof1.locally_owned_dofs().n_elements()));
        };
#endif

				     // For continuous elements on grids
				     // with hanging nodes we need
				     // hanging node
				     // constraints. Consequently, when
				     // the elements are continuous no
				     // hanging node constraints are
				     // allowed.
    const bool hanging_nodes_not_allowed=
      (dof1.get_fe().dofs_per_vertex != 0) || (fe2.dofs_per_vertex != 0);

    const unsigned int dofs_per_cell=dof1.get_fe().dofs_per_cell;

    Vector<typename OutVector::value_type> u1_local(dofs_per_cell);
    Vector<typename OutVector::value_type> u1_diff_local(dofs_per_cell);

    const types::subdomain_id_t subdomain_id =
      dof1.get_tria().locally_owned_subdomain();

    FullMatrix<double> difference_matrix(dofs_per_cell, dofs_per_cell);
    get_interpolation_difference_matrix(dof1.get_fe(), fe2,
					difference_matrix);

    typename DoFHandler<dim,spacedim>::active_cell_iterator cell = dof1.begin_active(),
							    endc = dof1.end();

    for (; cell!=endc; ++cell)
      if ((cell->subdomain_id() == subdomain_id)
          ||
          (subdomain_id == types::invalid_subdomain_id))
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

				   // if we work on a parallel PETSc vector
				   // we have to finish the work and
				   // update ghost values
      u1_difference.compress();
      u1_difference.update_ghost_values();
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void interpolation_difference(const DoFHandler<dim,spacedim> &dof1,
				const ConstraintMatrix &constraints1,
				const InVector &u1,
				const DoFHandler<dim,spacedim> &dof2,
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

				   // if we work on a parallel PETSc vector
				   // we have to finish the work and
				   // update ghost values
        u1_difference.compress();
        u1_difference.update_ghost_values();
      }
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void project_dg(const DoFHandler<dim,spacedim> &dof1,
		  const InVector &u1,
		  const DoFHandler<dim,spacedim> &dof2,
		  OutVector &u2)
  {
    Assert(&dof1.get_tria()==&dof2.get_tria(), ExcTriangulationMismatch());
    Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
	   ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
    Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
    Assert(u2.size()==dof2.n_dofs(), ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

    typename DoFHandler<dim,spacedim>::active_cell_iterator cell1 = dof1.begin_active();
    typename DoFHandler<dim,spacedim>::active_cell_iterator cell2 = dof2.begin_active();
    typename DoFHandler<dim,spacedim>::active_cell_iterator end = dof2.end();

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


  template <int dim, class InVector, class OutVector, int spacedim>
  void extrapolate(const DoFHandler<dim,spacedim> &dof1,
		   const InVector &u1,
		   const DoFHandler<dim,spacedim> &dof2,
		   OutVector &u2)
  {
    ConstraintMatrix dummy;
    dummy.close();
    extrapolate(dof1, u1, dof2, dummy, u2);
  }



  template <int dim, class InVector, class OutVector, int spacedim>
  void extrapolate(const DoFHandler<dim,spacedim> &dof1,
		   const InVector &u1,
		   const DoFHandler<dim,spacedim> &dof2,
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
      typename DoFHandler<dim,spacedim>::cell_iterator cell = dof2.begin(0),
						       endc = dof2.end(0);
      for (; cell!=endc; ++cell)
	Assert (cell->has_children(), ExcGridNotRefinedAtLeastOnce());
    }

				     // then traverse grid bottom up
    for (unsigned int level=0; level<dof1.get_tria().n_levels()-1; ++level)
      {
	typename DoFHandler<dim,spacedim>::cell_iterator cell=dof2.begin(level),
							 endc=dof2.end(level);

	for (; cell!=endc; ++cell)
	  if (!cell->active())
	    {
					       // check whether this
					       // cell has active
					       // children
	      bool active_children=false;
	      for (unsigned int child_n=0; child_n<cell->n_children(); ++child_n)
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


  template <>
  void
  add_fe_name(const std::string& parameter_name,
	      const FEFactoryBase<1,1>* factory)
  {
				     // Erase everything after the
				     // actual class name
    std::string name = parameter_name;
    unsigned int name_end =
      name.find_first_not_of(std::string("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"));
    if (name_end < name.size())
      name.erase(name_end);
				     // first make sure that no other
				     // thread intercepts the
				     // operation of this function;
				     // for this, acquire the lock
				     // until we quit this function
    Threads::ThreadMutex::ScopedLock lock(fe_name_map_lock);

    Assert(fe_name_map_1d.find(name) == fe_name_map_1d.end(),
	   ExcMessage("Cannot change existing element in finite element name list"));

				     // Insert the normalized name into
				     // the map
    fe_name_map_1d[name] =
      std_cxx1x::shared_ptr<const FETools::FEFactoryBase<1> > (factory);
  }



  template <>
  void
  add_fe_name(const std::string& parameter_name,
	      const FEFactoryBase<2,2>* factory)
  {
				     // Erase everything after the
				     // actual class name
    std::string name = parameter_name;
    unsigned int name_end =
      name.find_first_not_of(std::string("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"));
    if (name_end < name.size())
      name.erase(name_end);
				     // first make sure that no other
				     // thread intercepts the
				     // operation of this function;
				     // for this, acquire the lock
				     // until we quit this function
    Threads::ThreadMutex::ScopedLock lock(fe_name_map_lock);

    Assert(fe_name_map_2d.find(name) == fe_name_map_2d.end(),
	   ExcMessage("Cannot change existing element in finite element name list"));

				     // Insert the normalized name into
				     // the map
    fe_name_map_2d[name] =
      std_cxx1x::shared_ptr<const FETools::FEFactoryBase<2> > (factory);
  }


  template <>
  void
  add_fe_name(const std::string& parameter_name,
	      const FEFactoryBase<3,3>* factory)
  {
				     // Erase everything after the
				     // actual class name
    std::string name = parameter_name;
    unsigned int name_end =
      name.find_first_not_of(std::string("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"));
    if (name_end < name.size())
      name.erase(name_end);
				     // first make sure that no other
				     // thread intercepts the
				     // operation of this function;
				     // for this, acquire the lock
				     // until we quit this function
    Threads::ThreadMutex::ScopedLock lock(fe_name_map_lock);

    Assert(fe_name_map_3d.find(name) == fe_name_map_3d.end(),
	   ExcMessage("Cannot change existing element in finite element name list"));

				     // Insert the normalized name into
				     // the map
    fe_name_map_3d[name] =
      std_cxx1x::shared_ptr<const FETools::FEFactoryBase<3> > (factory);
  }


  namespace internal
  {
    namespace
    {
				// TODO: this encapsulates the call to the
				// dimension-dependent fe_name_map so that we
				// have a unique interface. could be done
				// smarter?
      template <int dim, int spacedim>
      FiniteElement<dim,spacedim>*
      get_fe_from_name_ext (std::string &name,
			    const std::map<std::string,
			    std_cxx1x::shared_ptr<const FETools::FEFactoryBase<dim> > >
			    &fe_name_map)
      {
					 // Extract the name of the
					 // finite element class, which only
					 // contains characters, numbers and
					 // underscores.
	unsigned int name_end =
	  name.find_first_not_of(std::string("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_"));
	const std::string name_part(name, 0, name_end);
	name.erase(0, name_part.size());

					 // now things get a little more
					 // complicated: FESystem. it's
					 // more complicated, since we
					 // have to figure out what the
					 // base elements are. this can
					 // only be done recursively
	if (name_part == "FESystem")
	  {
					     // next we have to get at the
					     // base elements. start with
					     // the first. wrap the whole
					     // block into try-catch to
					     // make sure we destroy the
					     // pointers we got from
					     // recursive calls if one of
					     // these calls should throw
					     // an exception
	    std::vector<FiniteElement<dim,spacedim>*> base_fes;
	    std::vector<unsigned int>        base_multiplicities;
	    try
	      {
						 // Now, just the [...]
						 // part should be left.
		if (name.size() == 0 || name[0] != '[')
		  throw (std::string("Invalid first character in ") + name);
		do
		  {
						     // Erase the
						     // leading '[' or '-'
		    name.erase(0,1);
						     // Now, the name of the
						     // first base element is
						     // first... Let's get it
		    base_fes.push_back (get_fe_from_name_ext<dim,spacedim> (name,
									    fe_name_map));
						     // next check whether
						     // FESystem placed a
						     // multiplicity after
						     // the element name
		    if (name[0] == '^')
		      {
							 // yes. Delete the '^'
							 // and read this
							 // multiplicity
			name.erase(0,1);

			const std::pair<int,unsigned int> tmp
			  = Utilities::get_integer_at_position (name, 0);
			name.erase(0, tmp.second);
							 // add to length,
							 // including the '^'
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
		while (name[0] == '-');

						 // so we got to the end
						 // of the '-' separated
						 // list. make sure that
						 // we actually had a ']'
						 // there
		if (name.size() == 0 || name[0] != ']')
		  throw (std::string("Invalid first character in ") + name);
		name.erase(0,1);
						 // just one more sanity check
		Assert ((base_fes.size() == base_multiplicities.size())
			&&
			(base_fes.size() > 0),
			ExcInternalError());

						 // ok, apparently
						 // everything went ok. so
						 // generate the composed
						 // element
		FiniteElement<dim,spacedim> *system_element = 0;
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
			  AssertThrow (false, ExcNotImplemented());
		  }

						 // now we don't need the
						 // list of base elements
						 // any more
		for (unsigned int i=0; i<base_fes.size(); ++i)
		  delete base_fes[i];

						 // finally return our
						 // findings
						 // Add the closing ']' to
						 // the length
		return system_element;

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
	else
	  {
					     // Make sure no other thread
					     // is just adding an element
	    Threads::ThreadMutex::ScopedLock lock (fe_name_map_lock);

	    AssertThrow (fe_name_map.find(name_part) != fe_name_map.end(),
			 ExcInvalidFEName(name));
					     // Now, just the (degree)
					     // or (Quadrature<1>(degree+1))
					     // part should be left.
	    if (name.size() == 0 || name[0] != '(')
	      throw (std::string("Invalid first character in ") + name);
	    name.erase(0,1);
	    if (name[0] != 'Q')
	      {
		const std::pair<int,unsigned int> tmp
		  = Utilities::get_integer_at_position (name, 0);
		name.erase(0, tmp.second+1);
		return fe_name_map.find(name_part)->second->get(tmp.first);
	      }
	    else
	      {
		unsigned int position = name.find('(');
		const std::string quadrature_name(name, 0, position);
		name.erase(0,position+1);
		if (quadrature_name.compare("QGaussLobatto") == 0)
		  {
		    const std::pair<int,unsigned int> tmp
		      = Utilities::get_integer_at_position (name, 0);
				// delete "))"
		    name.erase(0, tmp.second+2);
		    return fe_name_map.find(name_part)->second->get(QGaussLobatto<1>(tmp.first));
		  }
		else
		  {
		    AssertThrow (false,ExcNotImplemented());
		  }
	      }
	  }


					 // hm, if we have come thus far, we
					 // didn't know what to do with the
					 // string we got. so do as the docs
					 // say: raise an exception
	AssertThrow (false, ExcInvalidFEName(name));

					 // make some compilers happy that
					 // do not realize that we can't get
					 // here after throwing
	return 0;
      }



				// need to work around problem with different
				// name map names for different dimensions
				// TODO: should be possible to do nicer!
      template <int dim,int spacedim>
      FiniteElement<dim,spacedim>* get_fe_from_name (std::string &name);

      template <>
      FiniteElement<1,1>*
      get_fe_from_name<1> (std::string &name)
      {
	return get_fe_from_name_ext<1,1> (name, fe_name_map_1d);
      }

      template <>
      FiniteElement<2,2>*
      get_fe_from_name<2> (std::string &name)
      {
	return get_fe_from_name_ext<2,2> (name, fe_name_map_2d);
      }

      template <>
      FiniteElement<3,3>*
      get_fe_from_name<3> (std::string &name)
      {
	return get_fe_from_name_ext<3,3> (name, fe_name_map_3d);
      }
    }
  }





  template <int dim>
  FiniteElement<dim, dim> *
  get_fe_from_name (const std::string &parameter_name)
  {
				     // Create a version of the name
				     // string where all template
				     // parameters are eliminated.
    std::string name = parameter_name;
    for (unsigned int pos1 = name.find('<');
	 pos1 < name.size();
	 pos1 = name.find('<'))
      {

	const unsigned int pos2 = name.find('>');
					 // If there is only a single
					 // character between those two,
					 // it should be 'd' or the number
					 // representing the dimension.
	if (pos2-pos1 == 2)
	  {
	    const char dimchar = '0' + dim;
	    if (name.at(pos1+1) != 'd')
	      Assert (name.at(pos1+1) == dimchar,
		      ExcInvalidFEDimension(name.at(pos1+1), dim));
	  }
	else
	  Assert(pos2-pos1 == 4, ExcInvalidFEName(name));

					 // If pos1==pos2, then we are
					 // probably at the end of the
					 // string
	if (pos2 != pos1)
	  name.erase(pos1, pos2-pos1+1);
      }
				     // Replace all occurences of "^dim"
				     // by "^d" to be handled by the
				     // next loop
    for (unsigned int pos = name.find("^dim");
	 pos < name.size();
	 pos = name.find("^dim"))
      name.erase(pos+2, 2);

				     // Replace all occurences of "^d"
				     // by using the actual dimension
    for (unsigned int pos = name.find("^d");
	 pos < name.size();
	 pos = name.find("^d"))
      name.at(pos+1) = '0' + dim;

    try
      {
	FiniteElement<dim,dim> *fe = internal::get_fe_from_name<dim,dim> (name);

					 // Make sure the auxiliary function
					 // ate up all characters of the name.
	AssertThrow (name.size() == 0,
		     ExcInvalidFEName(parameter_name
				      + std::string(" extra characters after "
						    "end of name")));
	return fe;
      }
    catch (const std::string &errline)
      {
	AssertThrow(false, ExcInvalidFEName(parameter_name
					    + std::string(" at ")
					    + errline));
	return 0;
      }
  }


// template <int dim>
// FiniteElement<dim> *
// get_fe_from_name (const std::string &parameter_name)
// {
//     return internal::get_fe_from_name<dim,dim>(parameter_name);
// }


  template <int dim, int spacedim>
  void

  compute_projection_from_quadrature_points_matrix (const FiniteElement<dim,spacedim> &fe,
						    const Quadrature<dim>    &lhs_quadrature,
						    const Quadrature<dim>    &rhs_quadrature,
						    FullMatrix<double>       &X)
  {
    Assert (fe.n_components() == 1, ExcNotImplemented());

				     // first build the matrices M and Q
				     // described in the documentation
    FullMatrix<double> M (fe.dofs_per_cell, fe.dofs_per_cell);
    FullMatrix<double> Q (fe.dofs_per_cell, rhs_quadrature.size());

    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
	for (unsigned int q=0; q<lhs_quadrature.size(); ++q)
	  M(i,j) += fe.shape_value (i, lhs_quadrature.point(q)) *
		    fe.shape_value (j, lhs_quadrature.point(q)) *
		    lhs_quadrature.weight(q);

    for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
      for (unsigned int q=0; q<rhs_quadrature.size(); ++q)
	Q(i,q) += fe.shape_value (i, rhs_quadrature.point(q)) *
		  rhs_quadrature.weight(q);

				     // then invert M
    FullMatrix<double> M_inverse (fe.dofs_per_cell, fe.dofs_per_cell);
    M_inverse.invert (M);

				     // finally compute the result
    X.reinit (fe.dofs_per_cell, rhs_quadrature.size());
    M_inverse.mmult (X, Q);

    Assert (X.m() == fe.dofs_per_cell, ExcInternalError());
    Assert (X.n() == rhs_quadrature.size(), ExcInternalError());
  }



  template <int dim, int spacedim>
  void
  compute_interpolation_to_quadrature_points_matrix (const FiniteElement<dim,spacedim> &fe,
						     const Quadrature<dim>    &quadrature,
						     FullMatrix<double>       &I_q)
  {
    Assert (fe.n_components() == 1, ExcNotImplemented());
    Assert (I_q.m() == quadrature.size(),
	    ExcMessage ("Wrong matrix size"));
    Assert (I_q.n() == fe.dofs_per_cell, ExcMessage ("Wrong matrix size"));

    for (unsigned int q=0; q<quadrature.size(); ++q)
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	I_q(q,i) = fe.shape_value (i, quadrature.point(q));
  }



  template <int dim>
  void
  compute_projection_from_quadrature_points(
    const FullMatrix<double>                &projection_matrix,
    const std::vector< Tensor<1, dim > >    &vector_of_tensors_at_qp,
    std::vector< Tensor<1, dim > >          &vector_of_tensors_at_nodes)
  {

				     // check that the number columns of the projection_matrix
				     // matches the size of the vector_of_tensors_at_qp
    Assert(projection_matrix.n_cols() == vector_of_tensors_at_qp.size(),
	   ExcDimensionMismatch(projection_matrix.n_cols(),
				vector_of_tensors_at_qp.size()));

				     // check that the number rows of the projection_matrix
				     // matches the size of the vector_of_tensors_at_nodes
    Assert(projection_matrix.n_rows() == vector_of_tensors_at_nodes.size(),
	   ExcDimensionMismatch(projection_matrix.n_rows(),
				vector_of_tensors_at_nodes.size()));

				     // number of support points (nodes) to project to
    const unsigned int n_support_points = projection_matrix.n_rows();
				     // number of quadrature points to project from
    const unsigned int n_quad_points = projection_matrix.n_cols();

				     // component projected to the nodes
    Vector<double> component_at_node(n_support_points);
				     // component at the quadrature point
    Vector<double> component_at_qp(n_quad_points);

    for (unsigned int ii = 0; ii < dim; ++ii) {

      component_at_qp = 0;

				       // populate the vector of components at the qps
				       // from vector_of_tensors_at_qp
				       // vector_of_tensors_at_qp data is in form:
				       //      columns:        0, 1, ...,  dim
				       //      rows:           0,1,....,  n_quad_points
				       // so extract the ii'th column of vector_of_tensors_at_qp
      for (unsigned int q = 0; q < n_quad_points; ++q) {
	component_at_qp(q) = vector_of_tensors_at_qp[q][ii];
      }

				       // project from the qps -> nodes
				       // component_at_node = projection_matrix_u * component_at_qp
      projection_matrix.vmult(component_at_node, component_at_qp);

				       // rewrite the projection of the components
				       // back into the vector of tensors
      for (unsigned int nn =0; nn <n_support_points; ++nn) {
	vector_of_tensors_at_nodes[nn][ii] = component_at_node(nn);
      }
    }
  }



  template <int dim>
  void
  compute_projection_from_quadrature_points(
    const FullMatrix<double>                        &projection_matrix,
    const std::vector< SymmetricTensor<2, dim > >   &vector_of_tensors_at_qp,
    std::vector< SymmetricTensor<2, dim > >         &vector_of_tensors_at_nodes)
  {

				     // check that the number columns of the projection_matrix
				     // matches the size of the vector_of_tensors_at_qp
    Assert(projection_matrix.n_cols() == vector_of_tensors_at_qp.size(),
	   ExcDimensionMismatch(projection_matrix.n_cols(),
				vector_of_tensors_at_qp.size()));

				     // check that the number rows of the projection_matrix
				     // matches the size of the vector_of_tensors_at_nodes
    Assert(projection_matrix.n_rows() == vector_of_tensors_at_nodes.size(),
	   ExcDimensionMismatch(projection_matrix.n_rows(),
				vector_of_tensors_at_nodes.size()));

				     // number of support points (nodes)
    const unsigned int n_support_points = projection_matrix.n_rows();
				     // number of quadrature points to project from
    const unsigned int n_quad_points = projection_matrix.n_cols();

				     // number of unique entries in a symmetric second-order tensor
    const unsigned int n_independent_components =
      SymmetricTensor<2, dim >::n_independent_components;

				     // component projected to the nodes
    Vector<double> component_at_node(n_support_points);
				     // component at the quadrature point
    Vector<double> component_at_qp(n_quad_points);

				     // loop over the number of unique dimensions of the tensor
    for (unsigned int ii = 0; ii < n_independent_components; ++ii) {

      component_at_qp = 0;

				       // row-column entry of tensor corresponding the unrolled index
      TableIndices<2>  row_column_index = SymmetricTensor< 2, dim >::unrolled_to_component_indices(ii);
      const unsigned int row = row_column_index[0];
      const unsigned int column = row_column_index[1];

				       //  populate the vector of components at the qps
				       //  from vector_of_tensors_at_qp
				       //  vector_of_tensors_at_qp is in form:
				       //      columns:       0, 1, ..., n_independent_components
				       //      rows:           0,1,....,  n_quad_points
				       //  so extract the ii'th column of vector_of_tensors_at_qp
      for (unsigned int q = 0; q < n_quad_points; ++q) {
	component_at_qp(q) = (vector_of_tensors_at_qp[q])[row][column];
      }

				       // project from the qps -> nodes
				       // component_at_node = projection_matrix_u * component_at_qp
      projection_matrix.vmult(component_at_node, component_at_qp);

				       // rewrite the projection of the components back into the vector of tensors
      for (unsigned int nn =0; nn <n_support_points; ++nn) {
	(vector_of_tensors_at_nodes[nn])[row][column] = component_at_node(nn);
      }
    }
  }



  template <int dim, int spacedim>
  void
  compute_projection_from_face_quadrature_points_matrix (const FiniteElement<dim, spacedim> &fe,
							 const Quadrature<dim-1>    &lhs_quadrature,
							 const Quadrature<dim-1>    &rhs_quadrature,
							 const typename DoFHandler<dim, spacedim>::active_cell_iterator & cell,
							 unsigned int face,
							 FullMatrix<double>       &X)
  {
    Assert (fe.n_components() == 1, ExcNotImplemented());
    Assert (lhs_quadrature.size () > fe.degree, ExcNotGreaterThan (lhs_quadrature.size (), fe.degree));



				     // build the matrices M and Q
				     // described in the documentation
    FullMatrix<double> M (fe.dofs_per_cell, fe.dofs_per_cell);
    FullMatrix<double> Q (fe.dofs_per_cell, rhs_quadrature.size());

    {
				       // need an FEFaceValues object to evaluate shape function
				       // values on the specified face.
      FEFaceValues <dim> fe_face_values (fe, lhs_quadrature, update_values);
      fe_face_values.reinit (cell, face); // setup shape_value on this face.

      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	for (unsigned int j=0; j<fe.dofs_per_cell; ++j)
	  for (unsigned int q=0; q<lhs_quadrature.size(); ++q)
	    M(i,j) += fe_face_values.shape_value (i, q) *
		      fe_face_values.shape_value (j, q) *
		      lhs_quadrature.weight(q);
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	{
	  M(i,i) = (M(i,i) == 0 ? 1 : M(i,i));
	}
    }

    {
      FEFaceValues <dim> fe_face_values (fe, rhs_quadrature, update_values);
      fe_face_values.reinit (cell, face); // setup shape_value on this face.

      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	for (unsigned int q=0; q<rhs_quadrature.size(); ++q)
	  Q(i,q) += fe_face_values.shape_value (i, q) *
		    rhs_quadrature.weight(q);
    }
				     // then invert M
    FullMatrix<double> M_inverse (fe.dofs_per_cell, fe.dofs_per_cell);
    M_inverse.invert (M);

				     // finally compute the result
    X.reinit (fe.dofs_per_cell, rhs_quadrature.size());
    M_inverse.mmult (X, Q);

    Assert (X.m() == fe.dofs_per_cell, ExcInternalError());
    Assert (X.n() == rhs_quadrature.size(), ExcInternalError());
  }



  template <int dim>
  void
  hierarchic_to_lexicographic_numbering (const FiniteElementData<dim> &fe,
					 std::vector<unsigned int> &h2l)
  {
    Assert (h2l.size() == fe.dofs_per_cell,
	    ExcDimensionMismatch (h2l.size(), fe.dofs_per_cell));
    h2l = hierarchic_to_lexicographic_numbering (fe);
  }



  template <int dim>
  std::vector<unsigned int>
  hierarchic_to_lexicographic_numbering (const FiniteElementData<dim> &fe)
  {
    Assert (fe.n_components() == 1, ExcInvalidFE());

    std::vector<unsigned int> h2l (fe.dofs_per_cell);

    const unsigned int dofs_per_cell = fe.dofs_per_cell;
				     // polynomial degree
    const unsigned int degree = fe.dofs_per_line+1;
				     // number of grid points in each
				     // direction
    const unsigned int n = degree+1;

				     // the following lines of code are
				     // somewhat odd, due to the way the
				     // hierarchic numbering is
				     // organized. if someone would
				     // really want to understand these
				     // lines, you better draw some
				     // pictures where you indicate the
				     // indices and orders of vertices,
				     // lines, etc, along with the
				     // numbers of the degrees of
				     // freedom in hierarchical and
				     // lexicographical order
    switch (dim)
      {
	case 1:
	{
	  h2l[0] = 0;
	  h2l[1] = dofs_per_cell-1;
	  for (unsigned int i=2; i<dofs_per_cell; ++i)
	    h2l[i] = i-1;

	  break;
	}

	case 2:
	{
	  unsigned int next_index = 0;
					   // first the four vertices
	  h2l[next_index++] = 0;
	  h2l[next_index++] = n-1;
	  h2l[next_index++] = n*(n-1);
	  h2l[next_index++] = n*n-1;

					   // left   line
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (1+i)*n;

					   // right  line
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (2+i)*n-1;

					   // bottom line
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = 1+i;

					   // top    line
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n*(n-1)+i+1;

					   // inside quad
	  Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
		  ExcInternalError());
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = n*(i+1)+j+1;

	  Assert (next_index == fe.dofs_per_cell, ExcInternalError());

	  break;
	}

	case 3:
	{
	  unsigned int next_index = 0;
					   // first the eight vertices
	  h2l[next_index++] = 0;                 // 0
	  h2l[next_index++] = (      1)*degree;  // 1
	  h2l[next_index++] = (    n  )*degree;  // 2
	  h2l[next_index++] = (    n+1)*degree;  // 3
	  h2l[next_index++] = (n*n    )*degree;  // 4
	  h2l[next_index++] = (n*n  +1)*degree;  // 5
	  h2l[next_index++] = (n*n+n  )*degree;  // 6
	  h2l[next_index++] = (n*n+n+1)*degree;  // 7

					   // line 0
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (i+1)*n;
					   // line 1
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n-1+(i+1)*n;
					   // line 2
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = 1+i;
					   // line 3
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = 1+i+n*(n-1);

					   // line 4
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (n-1)*n*n+(i+1)*n;
					   // line 5
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (n-1)*(n*n+1)+(i+1)*n;
					   // line 6
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n*n*(n-1)+i+1;
					   // line 7
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n*n*(n-1)+i+1+n*(n-1);

					   // line 8
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (i+1)*n*n;
					   // line 9
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n-1+(i+1)*n*n;
					   // line 10
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (i+1)*n*n+n*(n-1);
					   // line 11
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n-1+(i+1)*n*n+n*(n-1);


					   // inside quads
	  Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
		  ExcInternalError());
					   // face 0
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = (i+1)*n*n+n*(j+1);
					   // face 1
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = (i+1)*n*n+n-1+n*(j+1);
					   // face 2, note the orientation!
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = (j+1)*n*n+i+1;
					   // face 3, note the orientation!
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = (j+1)*n*n+n*(n-1)+i+1;
					   // face 4
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = n*(i+1)+j+1;
					   // face 5
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = (n-1)*n*n+n*(i+1)+j+1;

					   // inside hex
	  Assert (fe.dofs_per_hex == fe.dofs_per_quad*fe.dofs_per_line,
		  ExcInternalError());
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      for (unsigned int k=0; k<fe.dofs_per_line; ++k)
		h2l[next_index++]	= n*n*(i+1)+n*(j+1)+k+1;

	  Assert (next_index == fe.dofs_per_cell, ExcInternalError());

	  break;
	}

	default:
	      Assert (false, ExcNotImplemented());
      }

    return h2l;
  }



  template <int dim>
  void
  lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe,
					 std::vector<unsigned int>    &l2h)
  {
    l2h = lexicographic_to_hierarchic_numbering (fe);
  }



  template <int dim>
  std::vector<unsigned int>
  lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe)
  {
    return Utilities::invert_permutation(hierarchic_to_lexicographic_numbering (fe));
  }

}



/*-------------- Explicit Instantiations -------------------------------*/
#include "fe_tools.inst"


/*----------------------------   fe_tools.cc     ---------------------------*/

DEAL_II_NAMESPACE_CLOSE
