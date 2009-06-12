//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <lac/vector.h>
#include <lac/petsc_vector.h>
#include <lac/trilinos_vector.h>

#include <grid/grid_refinement.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria.h>

#include <numeric>
#include <algorithm>
#include <cmath>
#include <functional>
#include <fstream>

DEAL_II_NAMESPACE_OPEN


namespace 
{
  template <typename number>
  inline
  number
  max_element (const Vector<number> &criteria)
  {
    return *std::max_element(criteria.begin(), criteria.end());
  }

  

  template <typename number>
  inline
  number
  min_element (const Vector<number> &criteria)
  {
    return *std::min_element(criteria.begin(), criteria.end());
  }


#ifdef DEAL_II_USE_PETSC  
  PetscScalar
  max_element (const PETScWrappers::Vector &criteria)
  {
                                     // this is horribly slow (since we have
                                     // to get the array of values from PETSc
                                     // in every iteration), but works
    PetscScalar m = 0;
    for (unsigned int i=0; i<criteria.size(); ++i)
      m = std::max (m, criteria(i));
    return m;
  }


  PetscScalar
  min_element (const PETScWrappers::Vector &criteria)
  {
                                     // this is horribly slow (since we have
                                     // to get the array of values from PETSc
                                     // in every iteration), but works
    PetscScalar m = criteria(0);
    for (unsigned int i=1; i<criteria.size(); ++i)
      m = std::min (m, criteria(i));
    return m;
  }
#endif


#ifdef DEAL_II_USE_TRILINOS  
  TrilinosScalar
  max_element (const TrilinosWrappers::Vector &criteria)
  {
    TrilinosScalar m = 0;
    criteria.trilinos_vector().MaxValue(&m);
    return m;
  }


  TrilinosScalar
  min_element (const TrilinosWrappers::Vector &criteria)
  {
    TrilinosScalar m = 0;
    criteria.trilinos_vector().MinValue(&m);
    return m;
  }
#endif
}


namespace
{
				   /**
				    * Sorts the vector @p ind as an
				    * index vector of @p a in
				    * increasing order.  This
				    * implementation of quicksort
				    * seems to be faster than the
				    * STL version and is needed in
				    * @p refine_and_coarsen_optimize
				    */

  template <class Vector>
  void qsort_index (const Vector              &a,
		    std::vector<unsigned int> &ind,
		    int                        l,
		    int                        r)
  {
    int i,j;
    typename Vector::value_type v;
  
    if (r<=l)
      return;

    v = a(ind[r]);
    i = l-1;
    j = r;
    do
      {
	do 
	  {
	    ++i;
	  }
	while ((a(ind[i])>v) && (i<r));
	do
	  {
	    --j;
	  }
	while ((a(ind[j])<v) && (j>0));
      
	if (i<j)
	  std::swap (ind[i], ind[j]);
	else
	  std::swap (ind[i], ind[r]);
      }
    while (i<j);
    qsort_index(a,ind,l,i-1);
    qsort_index(a,ind,i+1,r);  
  }
}




template <int dim, class Vector, int spacedim>
void GridRefinement::refine (Triangulation<dim,spacedim> &tria,
			     const Vector       &criteria,
			     const double        threshold)
{
  Assert (criteria.size() == tria.n_active_cells(),
	  ExcDimensionMismatch(criteria.size(), tria.n_active_cells()));
  Assert (criteria.is_non_negative (), ExcNegativeCriteria());
  
				   // when all indicators are zero we
				   // do not need to refine but only
				   // to coarsen
  if (criteria.all_zero())
    return;
  
  typename Triangulation<dim,spacedim>::active_cell_iterator cell = tria.begin_active();
  const unsigned int n_cells = criteria.size();

//TODO: This is undocumented, looks fishy and seems unnecessary
  
  double new_threshold=threshold;
				   // when threshold==0 find the
				   // smallest value in criteria
				   // greater 0
  if (new_threshold==0)
    {
      new_threshold = criteria(0);
      for (unsigned int index=1; index<n_cells; ++index)
	if (criteria(index)>0
	    && (criteria(index)<new_threshold))
	new_threshold=criteria(index);
    }
  
  for (unsigned int index=0; index<n_cells; ++cell, ++index)
    if (std::fabs(criteria(index)) >= new_threshold)
      cell->set_refine_flag();
}



template <int dim, class Vector, int spacedim>
void GridRefinement::coarsen (Triangulation<dim,spacedim> &tria,
			      const Vector       &criteria,
			      const double        threshold)
{
  Assert (criteria.size() == tria.n_active_cells(),
	  ExcDimensionMismatch(criteria.size(), tria.n_active_cells()));
  Assert (criteria.is_non_negative (), ExcNegativeCriteria());

  typename Triangulation<dim,spacedim>::active_cell_iterator cell = tria.begin_active();
  const unsigned int n_cells = criteria.size();
  
  for (unsigned int index=0; index<n_cells; ++cell, ++index)
    if (std::fabs(criteria(index)) <= threshold)
      if (!cell->refine_flag_set())
        cell->set_coarsen_flag();
}



template <int dim, class Vector, int spacedim>
void
GridRefinement::refine_and_coarsen_fixed_number (Triangulation<dim,spacedim> &tria,
						 const Vector       &criteria,
						 const double        top_fraction,
						 const double        bottom_fraction,
						 const unsigned int  max_n_cells)
{
				   // correct number of cells is
				   // checked in @p{refine}
  Assert ((top_fraction>=0) && (top_fraction<=1), ExcInvalidParameterValue());
  Assert ((bottom_fraction>=0) && (bottom_fraction<=1), ExcInvalidParameterValue());
  Assert (top_fraction+bottom_fraction <= 1, ExcInvalidParameterValue());
  Assert (criteria.is_non_negative (), ExcNegativeCriteria());

  int refine_cells  = static_cast<int>(top_fraction*criteria.size());
  int coarsen_cells = static_cast<int>(bottom_fraction*criteria.size());

				   // first we have to see whether we
				   // currently already exceed the target
				   // number of cells
  if (tria.n_active_cells() >= max_n_cells)
    {
				       // if yes, then we need to stop
				       // refining cells and instead try to
				       // only coarsen as many as it would
				       // take to get to the target
      
				       // as we have no information on cells
				       // being refined isotropically or
				       // anisotropically, assume isotropic
				       // refinement here, though that may
				       // result in a worse approximation
      refine_cells  = 0;
      coarsen_cells = (tria.n_active_cells() - max_n_cells) *
		      GeometryInfo<dim>::max_children_per_cell /
		      (GeometryInfo<dim>::max_children_per_cell - 1);
    }
				   // otherwise, see if we would exceed the
				   // maximum desired number of cells with the
				   // number of cells that are likely going to
				   // result from refinement. here, each cell
				   // to be refined is replaced by
				   // C=GeometryInfo<dim>::max_children_per_cell
				   // new cells, i.e. there will be C-1 more
				   // cells than before. similarly, C cells
				   // will be replaced by 1

				   // again, this is true for isotropically
				   // refined cells. we take this as an
				   // approximation of a mixed refinement.
  else if (static_cast<unsigned int>
	   (tria.n_active_cells()
	    + refine_cells * (GeometryInfo<dim>::max_children_per_cell - 1)
	    - (coarsen_cells *
	       (GeometryInfo<dim>::max_children_per_cell - 1) /
	       GeometryInfo<dim>::max_children_per_cell))
	   >
	   max_n_cells)
    {
				       // we have to adjust the
				       // fractions. assume we want
				       // alpha*refine_fraction and
				       // alpha*coarsen_fraction as new
				       // fractions and the resulting number
				       // of cells to be equal to
				       // max_n_cells. this leads to the
				       // following equation for lambda
      const double alpha
	=
	1. *
	(max_n_cells - tria.n_active_cells())
	/
	(refine_cells * (GeometryInfo<dim>::max_children_per_cell - 1)
	 - (coarsen_cells *
	    (GeometryInfo<dim>::max_children_per_cell - 1) /
	    GeometryInfo<dim>::max_children_per_cell));
      refine_cells  = static_cast<int> (refine_cells * alpha);
      coarsen_cells = static_cast<int> (coarsen_cells * alpha);
    }
  
  if (refine_cells || coarsen_cells)
    {
      dealii::Vector<typename Vector::value_type> tmp(criteria);
      if (refine_cells)
	{
	  std::nth_element (tmp.begin(), tmp.begin()+refine_cells,
			    tmp.end(),
			    std::greater<double>());
	  refine (tria, criteria, *(tmp.begin() + refine_cells));
	}

      if (coarsen_cells)
	{
	  std::nth_element (tmp.begin(), tmp.begin()+tmp.size()-coarsen_cells,
			    tmp.end(),
			    std::greater<double>());
	  coarsen (tria, criteria,
		   *(tmp.begin() + tmp.size() - coarsen_cells));
	}
    }
}



template <int dim, class Vector, int spacedim>
void
GridRefinement::refine_and_coarsen_fixed_fraction (Triangulation<dim,spacedim> &tria,
						   const Vector       &criteria,
						   const double        top_fraction,
						   const double        bottom_fraction,
						   const unsigned int  max_n_cells)
{
				   // correct number of cells is
				   // checked in @p{refine}
  Assert ((top_fraction>=0) && (top_fraction<=1), ExcInvalidParameterValue());
  Assert ((bottom_fraction>=0) && (bottom_fraction<=1), ExcInvalidParameterValue());
  Assert (top_fraction+bottom_fraction <= 1, ExcInvalidParameterValue());
  Assert (criteria.is_non_negative (), ExcNegativeCriteria());
  
				   // let tmp be the cellwise square of the
				   // error, which is what we have to sum
				   // up and compare with
				   // @p{fraction_of_error*total_error}.
  dealii::Vector<typename Vector::value_type> tmp(criteria);
  const double total_error = tmp.l1_norm();

				   // sort the largest criteria to the
				   // beginning of the vector
  std::sort (tmp.begin(), tmp.end(), std::greater<double>());

				   // compute thresholds
  typename dealii::Vector<typename Vector::value_type>::const_iterator
    pp=tmp.begin();
  for (double sum=0;
       (sum<top_fraction*total_error) && (pp!=(tmp.end()-1));
       ++pp)
    sum += *pp;
  double top_threshold = ( pp != tmp.begin () ?
			   (*pp+*(pp-1))/2 :
			   *pp );
  typename dealii::Vector<typename Vector::value_type>::const_iterator
    qq=(tmp.end()-1);
  for (double sum=0;
       (sum<bottom_fraction*total_error) && (qq!=tmp.begin());
       --qq)
    sum += *qq;
  double bottom_threshold = ( qq != (tmp.end()-1) ?
			      (*qq + *(qq+1))/2 :
			      0);

				   // we now have an idea how many cells we
				   // are going to refine and coarsen. we use
				   // this information to see whether we are
				   // over the limit and if so use a function
				   // that knows how to deal with this
				   // situation

				   // note, that at this point, we have no
				   // information about anisotropically refined
				   // cells, thus use the situation of purely
				   // isotropic refinement as guess for a mixed
				   // refinemnt as well.
  {
    const unsigned int refine_cells  = pp - tmp.begin(),
		       coarsen_cells = tmp.end() - qq;

    if (static_cast<unsigned int>
	(tria.n_active_cells()
	 + refine_cells * (GeometryInfo<dim>::max_children_per_cell - 1)
	 - (coarsen_cells *
	    (GeometryInfo<dim>::max_children_per_cell - 1) /
	    GeometryInfo<dim>::max_children_per_cell))
	>
	max_n_cells)
      {
	refine_and_coarsen_fixed_number (tria,
					 criteria,
					 1.*refine_cells/criteria.size(),
					 1.*coarsen_cells/criteria.size(),
					 max_n_cells);
	return;
      }
  }
  

				   // in some rare cases it may happen that
				   // both thresholds are the same (e.g. if
				   // there are many cells with the same
				   // error indicator). That would mean that
				   // all cells will be flagged for
				   // refinement or coarsening, but some will
				   // be flagged for both, namely those for
				   // which the indicator equals the
				   // thresholds. This is forbidden, however.
				   //
				   // In some rare cases with very few cells
				   // we also could get integer round off
				   // errors and get problems with
				   // the top and bottom fractions.
				   //
				   // In these case we arbitrarily reduce the
				   // bottom threshold by one permille below
				   // the top threshold
				   //
				   // Finally, in some cases
				   // (especially involving symmetric
				   // solutions) there are many cells
				   // with the same error indicator
				   // values. if there are many with
				   // indicator equal to the top
				   // threshold, no refinement will
				   // take place below; to avoid this
				   // case, we also lower the top
				   // threshold if it equals the
				   // largest indicator and the
				   // top_fraction!=1
  if ((top_threshold == max_element(criteria)) &&
      (top_fraction != 1))
    top_threshold *= 0.999;
  
  if (bottom_threshold>=top_threshold)
    bottom_threshold = 0.999*top_threshold;
  
				   // actually flag cells
  if (top_threshold < max_element(criteria))
    refine (tria, criteria, top_threshold);
  
  if (bottom_threshold > min_element(criteria))
    coarsen (tria, criteria, bottom_threshold);
}



template <int dim, class Vector, int spacedim>
void
GridRefinement::refine_and_coarsen_optimize (Triangulation<dim,spacedim> &tria,
					     const Vector       &criteria)
{
  Assert (criteria.size() == tria.n_active_cells(),
	  ExcDimensionMismatch(criteria.size(), tria.n_active_cells()));
  Assert (criteria.is_non_negative (), ExcNegativeCriteria());
  
				   // get an increasing order on
				   // the error indicator
  std::vector<unsigned int> tmp(criteria.size());
  for (unsigned int i=0;i<criteria.size();++i)
    tmp[i] = i;
  
  qsort_index(criteria,tmp,0,criteria.size()-1);
  
  double s0 = 0.75 * criteria(tmp[0]);
  double E  = criteria.l1_norm();
  
  unsigned int N = criteria.size();
  unsigned int M = 0;
  
				   // The first M cells are refined
				   // to minimize the expected error
				   // multiplied with the expected
				   // number of cells.
				   // We assume that the error is
				   // decreased by 3/4 a_K if the cell
				   // K with error indicator a_K is
				   // refined.
				   // The expected number of cells is
				   // N+3*M (N is the current number
				   // of cells)
  double min = (3.*(1.+M)+N) * (E-s0);
  
  unsigned int minArg = N-1;
  
  for (M=1;M<criteria.size();++M)
    {
      s0+= 0.75 * criteria(tmp[M]);
      
      if ( (3.*(1+M)+N)*(E-s0) <= min)
	{
	  min = (3.*(1+M)+N)*(E-s0);
	  minArg = M;
	}
    }
  refine(tria,criteria,criteria(tmp[minArg]));
}


// explicit instantiations
template
void
GridRefinement::
refine<deal_II_dimension,Vector<float>,deal_II_dimension>
       (Triangulation<deal_II_dimension> &,
        const Vector<float> &,
        const double);

template
void
GridRefinement::
refine<deal_II_dimension,Vector<double>,deal_II_dimension>
       (Triangulation<deal_II_dimension> &,
        const Vector<double> &,
        const double);

template
void
GridRefinement::
coarsen<deal_II_dimension,Vector<float>,deal_II_dimension>
       (Triangulation<deal_II_dimension> &,
         const Vector<float> &,
         const double);

template
void
GridRefinement::
coarsen<deal_II_dimension,Vector<double>,deal_II_dimension>
       (Triangulation<deal_II_dimension> &,
         const Vector<double> &,
         const double);


template
void
GridRefinement::
refine_and_coarsen_fixed_number<deal_II_dimension,Vector<double>,deal_II_dimension>
(Triangulation<deal_II_dimension> &,
 const Vector<double> &,
 const double,
 const double,
 const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_fixed_number<deal_II_dimension,Vector<float>,deal_II_dimension>
 (Triangulation<deal_II_dimension> &,
  const Vector<float> &,
  const double,
  const double,
  const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_fixed_fraction<deal_II_dimension,Vector<double>,deal_II_dimension>
 (Triangulation<deal_II_dimension> &,
  const Vector<double> &,
  const double,
  const double,
  const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_fixed_fraction<deal_II_dimension,Vector<float>,deal_II_dimension>
 (Triangulation<deal_II_dimension> &,
  const Vector<float> &,
  const double,
  const double,
  const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_optimize<deal_II_dimension,Vector<float>,deal_II_dimension>
(Triangulation<deal_II_dimension> &,
 const Vector<float> &);

template
void
GridRefinement::
refine_and_coarsen_optimize<deal_II_dimension,Vector<double>,deal_II_dimension>
(Triangulation<deal_II_dimension> &,
 const Vector<double> &);

#ifdef DEAL_II_USE_PETSC
template
void
GridRefinement::
refine<deal_II_dimension,PETScWrappers::Vector,deal_II_dimension>
 (Triangulation<deal_II_dimension> &,
  const PETScWrappers::Vector &,
  const double);

template
void
GridRefinement::
coarsen<deal_II_dimension,PETScWrappers::Vector,deal_II_dimension>
(Triangulation<deal_II_dimension> &,
 const PETScWrappers::Vector &,
 const double);


template
void
GridRefinement::
refine_and_coarsen_fixed_number<deal_II_dimension,PETScWrappers::Vector,deal_II_dimension>
(Triangulation<deal_II_dimension> &,
 const PETScWrappers::Vector &,
 const double,
 const double,
 const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_fixed_fraction<deal_II_dimension,PETScWrappers::Vector,deal_II_dimension>
(Triangulation<deal_II_dimension> &,
 const PETScWrappers::Vector &,
 const double,
 const double,
 const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_optimize<deal_II_dimension,PETScWrappers::Vector,deal_II_dimension>
(Triangulation<deal_II_dimension> &,
 const PETScWrappers::Vector &);
#endif

#ifdef DEAL_II_USE_TRILINOS
template
void
GridRefinement::
refine<deal_II_dimension,TrilinosWrappers::Vector,deal_II_dimension>
 (Triangulation<deal_II_dimension> &,
  const TrilinosWrappers::Vector &,
  const double);

template
void
GridRefinement::
coarsen<deal_II_dimension,TrilinosWrappers::Vector,deal_II_dimension>
(Triangulation<deal_II_dimension> &,
 const TrilinosWrappers::Vector &,
 const double);


template
void
GridRefinement::
refine_and_coarsen_fixed_number<deal_II_dimension,TrilinosWrappers::Vector,deal_II_dimension>
(Triangulation<deal_II_dimension> &,
 const TrilinosWrappers::Vector &,
 const double,
 const double,
 const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_fixed_fraction<deal_II_dimension,TrilinosWrappers::Vector,deal_II_dimension>
(Triangulation<deal_II_dimension> &,
 const TrilinosWrappers::Vector &,
 const double,
 const double,
 const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_optimize<deal_II_dimension,TrilinosWrappers::Vector,deal_II_dimension>
(Triangulation<deal_II_dimension> &,
 const TrilinosWrappers::Vector &);
#endif



#if deal_II_dimension != 3
// explicit instantiations
template
void
GridRefinement::
refine<deal_II_dimension,Vector<float>,deal_II_dimension+1>
       (Triangulation<deal_II_dimension,deal_II_dimension+1> &,
        const Vector<float> &,
        const double);

template
void
GridRefinement::
refine<deal_II_dimension,Vector<double>,deal_II_dimension+1>
       (Triangulation<deal_II_dimension,deal_II_dimension+1> &,
        const Vector<double> &,
        const double);

template
void
GridRefinement::
coarsen<deal_II_dimension,Vector<float>,deal_II_dimension+1>
       (Triangulation<deal_II_dimension,deal_II_dimension+1> &,
         const Vector<float> &,
         const double);

template
void
GridRefinement::
coarsen<deal_II_dimension,Vector<double>,deal_II_dimension+1>
       (Triangulation<deal_II_dimension,deal_II_dimension+1> &,
         const Vector<double> &,
         const double);


template
void
GridRefinement::
refine_and_coarsen_fixed_number<deal_II_dimension,Vector<double>,deal_II_dimension+1>
(Triangulation<deal_II_dimension,deal_II_dimension+1> &,
 const Vector<double> &,
 const double,
 const double,
 const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_fixed_number<deal_II_dimension,Vector<float>,deal_II_dimension+1>
 (Triangulation<deal_II_dimension,deal_II_dimension+1> &,
  const Vector<float> &,
  const double,
  const double,
  const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_fixed_fraction<deal_II_dimension,Vector<double>,deal_II_dimension+1>
 (Triangulation<deal_II_dimension,deal_II_dimension+1> &,
  const Vector<double> &,
  const double,
  const double,
  const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_fixed_fraction<deal_II_dimension,Vector<float>,deal_II_dimension+1>
 (Triangulation<deal_II_dimension,deal_II_dimension+1> &,
  const Vector<float> &,
  const double,
  const double,
  const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_optimize<deal_II_dimension,Vector<float>,deal_II_dimension+1>
(Triangulation<deal_II_dimension,deal_II_dimension+1> &,
 const Vector<float> &);

template
void
GridRefinement::
refine_and_coarsen_optimize<deal_II_dimension,Vector<double>,deal_II_dimension+1>
(Triangulation<deal_II_dimension,deal_II_dimension+1> &,
 const Vector<double> &);

#ifdef DEAL_II_USE_PETSC
template
void
GridRefinement::
refine<deal_II_dimension,PETScWrappers::Vector,deal_II_dimension+1>
 (Triangulation<deal_II_dimension,deal_II_dimension+1> &,
  const PETScWrappers::Vector &,
  const double);

template
void
GridRefinement::
coarsen<deal_II_dimension,PETScWrappers::Vector,deal_II_dimension+1>
(Triangulation<deal_II_dimension,deal_II_dimension+1> &,
 const PETScWrappers::Vector &,
 const double);


template
void
GridRefinement::
refine_and_coarsen_fixed_number<deal_II_dimension,PETScWrappers::Vector,deal_II_dimension+1>
(Triangulation<deal_II_dimension,deal_II_dimension+1> &,
 const PETScWrappers::Vector &,
 const double,
 const double,
 const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_fixed_fraction<deal_II_dimension,PETScWrappers::Vector,deal_II_dimension+1>
(Triangulation<deal_II_dimension,deal_II_dimension+1> &,
 const PETScWrappers::Vector &,
 const double,
 const double,
 const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_optimize<deal_II_dimension,PETScWrappers::Vector,deal_II_dimension+1>
(Triangulation<deal_II_dimension,deal_II_dimension+1> &,
 const PETScWrappers::Vector &);
#endif

#ifdef DEAL_II_USE_TRILINOS
template
void
GridRefinement::
refine<deal_II_dimension,TrilinosWrappers::Vector,deal_II_dimension+1>
 (Triangulation<deal_II_dimension,deal_II_dimension+1> &,
  const TrilinosWrappers::Vector &,
  const double);

template
void
GridRefinement::
coarsen<deal_II_dimension,TrilinosWrappers::Vector,deal_II_dimension+1>
(Triangulation<deal_II_dimension,deal_II_dimension+1> &,
 const TrilinosWrappers::Vector &,
 const double);


template
void
GridRefinement::
refine_and_coarsen_fixed_number<deal_II_dimension,TrilinosWrappers::Vector,deal_II_dimension+1>
(Triangulation<deal_II_dimension,deal_II_dimension+1> &,
 const TrilinosWrappers::Vector &,
 const double,
 const double,
 const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_fixed_fraction<deal_II_dimension,TrilinosWrappers::Vector,deal_II_dimension+1>
(Triangulation<deal_II_dimension,deal_II_dimension+1> &,
 const TrilinosWrappers::Vector &,
 const double,
 const double,
 const unsigned int);

template
void
GridRefinement::
refine_and_coarsen_optimize<deal_II_dimension,TrilinosWrappers::Vector,deal_II_dimension+1>
(Triangulation<deal_II_dimension,deal_II_dimension+1> &,
 const TrilinosWrappers::Vector &);
#endif


#endif

DEAL_II_NAMESPACE_CLOSE

