/* $Id$ */
/* Author: Wolfgang Bangerth, Texas A&M University, 2007 */

/*    $Id$       */
/*    Version: $Name$                                          */
/*                                                                */
/*    Copyright (C) 2007, 2008 by the deal.II authors */
/*                                                                */
/*    This file is subject to QPL and may not be  distributed     */
/*    without copyright and license information. Please refer     */
/*    to the file deal.II/doc/license.html for the  text  and     */
/*    further information on this license.                        */



#include <base/quadrature_lib.h>
#include <base/logstream.h>
#include <base/function.h>
#include <base/utilities.h>

#include <lac/block_vector.h>
#include <lac/full_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <lac/solver_cg.h>
#include <lac/precondition.h>
#include <lac/sparse_direct.h>

#include <grid/tria.h>
#include <grid/grid_generator.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/tria_boundary_lib.h>
#include <grid/grid_tools.h>
#include <grid/grid_refinement.h>

#include <dofs/dof_handler.h>
#include <dofs/dof_renumbering.h>
#include <dofs/dof_accessor.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_constraints.h>

#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_system.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>
#include <fe/mapping_c1.h>

#include <numerics/vectors.h>
#include <numerics/matrices.h>
#include <numerics/data_out.h>
#include <numerics/derivative_approximation.h>
#include <numerics/solution_transfer.h>

#include <fstream>
#include <sstream>

using namespace dealii;


                                 
template <int dim>
class BoussinesqFlowProblem 
{
  public:
    BoussinesqFlowProblem (const unsigned int degree);
    void run ();
    
  private:
    void setup_dofs (const bool setup_matrices);
    void assemble_system ();
    void assemble_rhs_T ();
    double get_maximal_velocity () const;
    void solve ();
    void output_results () const;
    void refine_mesh ();
    
    const unsigned int   degree;
    
    Triangulation<dim>   triangulation;
    FESystem<dim>        fe;
    DoFHandler<dim>      dof_handler;

    ConstraintMatrix     hanging_node_constraints;
    
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;

    double time_step;
    unsigned int timestep_number;
 
    BlockVector<double> solution;
    BlockVector<double> old_solution;
    BlockVector<double> system_rhs;

    boost::shared_ptr<SparseDirectUMFPACK> A_preconditioner;

    bool rebuild_matrices;
    bool rebuild_preconditioner;
};





template <int dim>
class PressureBoundaryValues : public Function<dim> 
{
  public:
    PressureBoundaryValues () : Function<dim>(1) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};


template <int dim>
double
PressureBoundaryValues<dim>::value (const Point<dim>  &/*p*/,
                                    const unsigned int /*component*/) const 
{
  return 0;
}



template <int dim>
class TemperatureBoundaryValues : public Function<dim> 
{
  public:
    TemperatureBoundaryValues () : Function<dim>(1) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;
};



template <int dim>
double
TemperatureBoundaryValues<dim>::value (const Point<dim> &p,
                                      const unsigned int /*component*/) const 
{
//TODO: leftover from olden times. replace by something sensible once we have
//diffusion in the temperature field
  if (p[0] == 0)
    return 1;
  else
    return 0;
}




template <int dim>
class InitialValues : public Function<dim> 
{
  public:
    InitialValues () : Function<dim>(dim+2) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p, 
                               Vector<double>   &value) const;
};


template <int dim>
double
InitialValues<dim>::value (const Point<dim>  &,
                           const unsigned int) const 
{
  return 0;
}


template <int dim>
void
InitialValues<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &values) const 
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = InitialValues<dim>::value (p, c);
}



template <int dim>
class RightHandSide : public Function<dim> 
{
  public:
    RightHandSide () : Function<dim>(dim+2) {}
    
    virtual double value (const Point<dim>   &p,
                          const unsigned int  component = 0) const;

    virtual void vector_value (const Point<dim> &p, 
                               Vector<double>   &value) const;
};


template <int dim>
double
RightHandSide<dim>::value (const Point<dim>  &p,
                           const unsigned int component) const 
{
  if (component == dim+1)
    return ((p.distance (Point<dim>(.3,.1)) < 1./32)
	    ||
	    (p.distance (Point<dim>(.45,.1)) < 1./32)
	    ||
	    (p.distance (Point<dim>(.75,.1)) < 1./32)
	    ?
	    1
	    :
	    0);
  else
    return 0;
}


template <int dim>
void
RightHandSide<dim>::vector_value (const Point<dim> &p,
                                  Vector<double>   &values) const 
{
  for (unsigned int c=0; c<this->n_components; ++c)
    values(c) = RightHandSide<dim>::value (p, c);
}



namespace internals
{
  namespace VectorTools
  {
				     /**
				      * A structure that stores the dim DoF
				      * indices that correspond to a
				      * vector-valued quantity at a single
				      * support point.
				      */
    template <int dim>
    struct VectorDoFTuple
    {
	unsigned int dof_indices[dim];

	bool operator < (const VectorDoFTuple<dim> &other) const
	  {
	    for (unsigned int i=0; i<dim; ++i)
	      if (dof_indices[i] < other.dof_indices[i])
		return true;
	      else
		if (dof_indices[i] > other.dof_indices[i])
		  return false;
	    return false;
	  }

	bool operator == (const VectorDoFTuple<dim> &other) const
	  {
	    for (unsigned int i=0; i<dim; ++i)
	      if (dof_indices[i] != other.dof_indices[i])
		return false;

	    return true;
	  }

	bool operator != (const VectorDoFTuple<dim> &other) const
	  {
	    return ! (*this == other);
	  }
    };
  }
}



template <int dim, template <int> class DH>
void
compute_no_normal_flux_constraints (const DH<dim>         &dof_handler,
				    const unsigned int     first_vector_component,
				    const std::set<unsigned char> &boundary_ids,
				    ConstraintMatrix      &constraints,
				    const Mapping<dim>    &mapping = StaticMappingQ1<dim>::mapping)
{
  Assert (dim > 1,
	  ExcMessage ("This function is not useful in 1d because it amounts "
		      "to imposing Dirichlet values on the vector-valued "
		      "quantity."));
	  
  const FiniteElement<dim> &fe = dof_handler.get_fe();

  std::vector<unsigned int> face_dofs (fe.dofs_per_face);
  std::vector<Point<dim> >  dof_locations  (fe.dofs_per_face);

				   // have a map that stores normal vectors
				   // for each vector-dof tuple we want to
				   // constrain. since we can get at the same
				   // vector dof tuple more than once (for
				   // example if it is located at a vertex
				   // that we visit from all adjacent cells),
				   // we will want to average later on the
				   // normal vectors computed on different
				   // cells as described in the documentation
				   // of this function. however, we can only
				   // average if the contributions came from
				   // different cells, whereas we want to
				   // constrain twice or more in case the
				   // contributions came from different faces
				   // of the same cell. consequently, we also
				   // have to store which cell a normal vector
				   // was computed on
  typedef 
    std::multimap< ::internals::VectorTools::VectorDoFTuple<dim>,
    std::pair<Tensor<1,dim>, typename DH<dim>::active_cell_iterator> >
    DoFToNormalsMap;

  DoFToNormalsMap dof_to_normals_map;

				   // now loop over all cells and all faces
  typename DH<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    for (unsigned int face_no=0; face_no < GeometryInfo<dim>::faces_per_cell;
	 ++face_no)
      if (boundary_ids.find(cell->face(face_no)->boundary_indicator())
	  != boundary_ids.end())
	{
	  typename DH<dim>::face_iterator face = cell->face(face_no);

					   // get the indices of the
					   // dofs on this cell...
	  face->get_dof_indices (face_dofs, cell->active_fe_index());

					   // ...and the normal
					   // vectors at the locations
					   // where they are defined:
	  const std::vector<Point<dim-1> > &
	    unit_support_points = fe.get_unit_face_support_points();  
	  Quadrature<dim-1> aux_quad (unit_support_points);
	  FEFaceValues<dim> fe_values (mapping, fe, aux_quad,
				       update_normal_vectors);
	  fe_values.reinit(cell, face_no);

					   // then identify which of
					   // them correspond to the
					   // selected set of vector
					   // components
	  for (unsigned int i=0; i<face_dofs.size(); ++i)
	    if (fe.face_system_to_component_index(i).first ==
		first_vector_component)
	      {
						 // find corresponding other
						 // components of vector
		::internals::VectorTools::VectorDoFTuple<dim> vector_dofs;
		vector_dofs.dof_indices[0] = face_dofs[i];
		
		for (unsigned int k=0; k<fe.dofs_per_face; ++k)
		  if ((k != i)
		      &&
		      (unit_support_points[k] == unit_support_points[i])
		      &&
		      (fe.face_system_to_component_index(k).first >=
		      first_vector_component)
		      &&
		      (fe.face_system_to_component_index(k).first <
		       first_vector_component + dim))
		    vector_dofs.dof_indices[fe.face_system_to_component_index(k).first]
		      = face_dofs[k];

						 // and enter the
						 // (dofs,(normal_vector,cell))
						 // entry into the map
		dof_to_normals_map
		  .insert (std::make_pair (vector_dofs,
					   std::make_pair (fe_values.normal_vector(i),
							   cell)));
	      }
	}

				   // Now do something with the
				   // collected information. To this
				   // end, loop through all sets of
				   // pairs (dofs,normal_vector) and
				   // identify which entries belong to
				   // the same set of dofs and then do
				   // as described in the
				   // documentation, i.e. either
				   // average the normal vector or
				   // don't for this particular set of
				   // dofs
  typename DoFToNormalsMap::const_iterator
    p = dof_to_normals_map.begin();
  
  while (p != dof_to_normals_map.end())
    {
				       // first find the range of entries in
				       // the multimap that corresponds to the
				       // same vector-dof tuple. as usual, we
				       // define the range half-open. the
				       // first entry of course is 'p'
      typename DoFToNormalsMap::const_iterator same_dof_range[2]
	= { p };
      for (++p; p != dof_to_normals_map.end(); ++p)
	if (p->first != same_dof_range[0]->first)
	  {
	    same_dof_range[1] = p;
	    break;
	  }
      if (p == dof_to_normals_map.end())
	same_dof_range[1] = dof_to_normals_map.end();

				       // now compute the reverse mapping: for
				       // each of the cells that contributed
				       // to the current set of vector dofs,
				       // add up the normal vectors. the
				       // values of the map are pairs of
				       // normal vectors and number of cells
				       // that have contributed
      typedef
	std::map
	<typename DH<dim>::active_cell_iterator,
	std::pair<Tensor<1,dim>, unsigned int> >
	CellToNormalsMap;

      CellToNormalsMap cell_to_normals_map;
      for (typename DoFToNormalsMap::const_iterator
	     q = same_dof_range[0];
	   q != same_dof_range[1]; ++q)
	if (cell_to_normals_map.find (q->second.second)
	    == cell_to_normals_map.end())
	  cell_to_normals_map[q->second.second]
	    = std::make_pair (q->second.first, 1U);
	else
	  {
	    const Tensor<1,dim> old_normal
	      = cell_to_normals_map[q->second.second].first;
	    const unsigned int old_count
	      = cell_to_normals_map[q->second.second].second;

	    Assert (old_count > 0, ExcInternalError());

					     // in the same entry,
					     // store again the now
					     // averaged normal vector
					     // and the new count
	    cell_to_normals_map[q->second.second]
	      = std::make_pair ((old_normal * old_count + q->second.first) / (old_count + 1),
				old_count + 1);
	  }

      Assert (cell_to_normals_map.size() >= 1, ExcInternalError());

				       // count the maximum number of
				       // contributions from each cell
      unsigned int max_n_contributions_per_cell = 1;
      for (typename CellToNormalsMap::const_iterator
	     x = cell_to_normals_map.begin();
	   x != cell_to_normals_map.end(); ++x)
	max_n_contributions_per_cell
	  = std::max (max_n_contributions_per_cell,
		      x->second.second);

				       // verify that each cell can have only
				       // contributed at most dim times, since
				       // that is the maximum number of faces
				       // that come together at a single place
      Assert (max_n_contributions_per_cell <= dim, ExcInternalError());

      switch (max_n_contributions_per_cell)
	{
					   // first deal with the case that a
					   // number of cells all have
					   // registered that they have a
					   // normal vector defined at the
					   // location of a given vector dof,
					   // and that each of them have
					   // encountered this vector dof
					   // exactly once while looping over
					   // all their faces. as stated in
					   // the documentation, this is the
					   // case where we want to simply
					   // average over all normal vectors
	  case 1:
	  {
	    
					     // compute the average
					     // normal vector from all
					     // the ones that have the
					     // same set of dofs. we
					     // could add them up and
					     // divide them by the
					     // number of additions,
					     // or simply normalize
					     // them right away since
					     // we want them to have
					     // unit length anyway
	    Tensor<1,dim> normal;
	    for (typename CellToNormalsMap::const_iterator
		   x = cell_to_normals_map.begin();
		 x != cell_to_normals_map.end(); ++x)
	      normal += x->second.first;
	    normal /= normal.norm();

	    const ::internals::VectorTools::VectorDoFTuple<dim> &
		dof_indices = same_dof_range[0]->first;
	    
	    if (cell_to_normals_map.size() == 2)
	      {
		std::cout << "XX " << same_dof_range[0]->first.dof_indices[0]
			  << ' ' <<  same_dof_range[0]->first.dof_indices[1]
			  << std::endl;
		std::cout << "   " << cell_to_normals_map.begin()->first
			  << ' '   << cell_to_normals_map.begin()->second.first
			  << std::endl;
		std::cout << "   " << (++cell_to_normals_map.begin())->first
			  << ' '   << (++cell_to_normals_map.begin())->second.first
			  << std::endl;
		std::cout << "   " << normal << std::endl;
	      }
	    if (cell_to_normals_map.size() == 1)
	      {
		std::cout << "YY " << same_dof_range[0]->first.dof_indices[0]
			  << ' ' <<  same_dof_range[0]->first.dof_indices[1]
			  << std::endl;
		std::cout << "   " << cell_to_normals_map.begin()->first
			  << ' '   << cell_to_normals_map.begin()->second.first
			  << std::endl;
		std::cout << "   " << normal << std::endl;
	      }
	    
	    
					     // then construct constraints
					     // from this. choose the DoF that
					     // has the largest component in
					     // the normal vector as the one
					     // to be constrained as this
					     // makes the process stable in
					     // cases where the normal vector
					     // has the form n=(1,0) or
					     // n=(0,1)
	    switch (dim)
	      {
		case 2:
		{
		  if (std::fabs(normal[0]) > std::fabs(normal[1]))
		    {
		      constraints.add_line (dof_indices.dof_indices[0]);
		      constraints.add_entry (dof_indices.dof_indices[0],
					     dof_indices.dof_indices[1],
					     -normal[1]/normal[0]);
		    }
		  else
		    {
		      constraints.add_line (dof_indices.dof_indices[1]);
		      constraints.add_entry (dof_indices.dof_indices[1],
					     dof_indices.dof_indices[0],
					     -normal[0]/normal[1]);
		    }
		  break;
		}

		case 3:
		{
		  if ((std::fabs(normal[0]) >= std::fabs(normal[1]))
		      &&
		      (std::fabs(normal[0]) >= std::fabs(normal[2])))
		    {
		      constraints.add_line (dof_indices.dof_indices[0]);
		      constraints.add_entry (dof_indices.dof_indices[0],
					     dof_indices.dof_indices[1],
					     -normal[1]/normal[0]);
		      constraints.add_entry (dof_indices.dof_indices[0],
					     dof_indices.dof_indices[2],
					     -normal[2]/normal[0]);
		    }
		  else
		    if ((std::fabs(normal[1]) >= std::fabs(normal[0]))
			&&
			(std::fabs(normal[1]) >= std::fabs(normal[2])))
		      {
			constraints.add_line (dof_indices.dof_indices[1]);
			constraints.add_entry (dof_indices.dof_indices[1],
					       dof_indices.dof_indices[0],
					       -normal[0]/normal[1]);
			constraints.add_entry (dof_indices.dof_indices[1],
					       dof_indices.dof_indices[2],
					       -normal[2]/normal[1]);
		      }
		    else
		      {
			constraints.add_line (dof_indices.dof_indices[2]);
			constraints.add_entry (dof_indices.dof_indices[2],
					       dof_indices.dof_indices[0],
					       -normal[0]/normal[2]);
			constraints.add_entry (dof_indices.dof_indices[2],
					       dof_indices.dof_indices[1],
					       -normal[1]/normal[2]);
		      }
	    
		  break;
		}

		default:
		      Assert (false, ExcNotImplemented());
	      }

	    break;
	  }


					   // this is the slightly
					   // more complicated case
					   // that a single cell has
					   // contributed with exactly
					   // DIM normal vectors to
					   // the same set of vector
					   // dofs. this is what
					   // happens in a corner in
					   // 2d and 3d (but not on an
					   // edge in 3d, where we
					   // have only 2, i.e. <DIM,
					   // contributions. Here we
					   // do not want to average
					   // the normal
					   // vectors. Since we have
					   // DIM contributions, let's
					   // assume (and verify) that
					   // they are in fact all
					   // linearly independent; in
					   // that case, all vector
					   // components are
					   // constrained and we need
					   // to set them to zero
	  case dim:
	  {
					     // assert that indeed
					     // only a single cell has
					     // contributed
	    Assert (cell_to_normals_map.size() == 1,
		    ExcInternalError());

					     // check linear
					     // independence by
					     // computing the
					     // determinant of the
					     // matrix created from
					     // all the normal
					     // vectors. if they are
					     // linearly independent,
					     // then the determinant
					     // is nonzero. if they
					     // are orthogonal, then
					     // the matrix is in fact
					     // equal to 1 (since they
					     // are all unit vectors);
					     // make sure the
					     // determinant is larger
					     // than 1e-3 to avoid
					     // cases where cells are
					     // degenerate
	    {
	      Tensor<2,dim> t;

	      typename DoFToNormalsMap::const_iterator x = same_dof_range[0];
	      for (unsigned int i=0; i<dim; ++i, ++x)
		for (unsigned int j=0; j<dim; ++j)
		  t[i][j] = x->second.first[j];
	      
	      Assert (std::fabs(determinant (t)) > 1e-3,
		      ExcMessage("Found a set of normal vectors that are nearly collinear."));
	    }
	    
					     // so all components of
					     // this vector dof are
					     // constrained. enter
					     // this into the
					     // constraint matrix
	    for (unsigned int i=0; i<dim; ++i)
	      {
		constraints.add_line (same_dof_range[0]->first.dof_indices[i]);
						 // no add_entries here
	      }
	    
	    std::cout << "DIM contributions at "
		      << same_dof_range[0]->first.dof_indices[0]
		      << ' '
		      << same_dof_range[0]->first.dof_indices[1]
		      << std::endl;
	    break;
	  }
	  

					   // this is the case of an
					   // edge contribution in 3d
	  default:
		Assert (dim >= 3, ExcInternalError());
		Assert (false, ExcNotImplemented());
	}
    }
}






template <int dim>
Tensor<1,dim>
extract_u (const FEValuesBase<dim> &fe_values,
           const unsigned int i,
           const unsigned int q)
{
  Tensor<1,dim> tmp;

  const unsigned int component
    = fe_values.get_fe().system_to_component_index(i).first;

  if (component < dim)
    tmp[component] = fe_values.shape_value (i,q);

  return tmp;
}



template <int dim>
double
extract_div_u (const FEValuesBase<dim> &fe_values,
               const unsigned int i,
               const unsigned int q)
{
  const unsigned int component
    = fe_values.get_fe().system_to_component_index(i).first;

  if (component < dim)
    return fe_values.shape_grad (i,q)[component];
  else
    return 0;
}



template <int dim>
Tensor<2,dim>
extract_grad_s_u (const FEValuesBase<dim> &fe_values,
		  const unsigned int i,
		  const unsigned int q)
{
  Tensor<2,dim> tmp;

  const unsigned int component
    = fe_values.get_fe().system_to_component_index(i).first;
  
  if (component < dim)
    {
      const Tensor<1,dim> grad_phi_over_2 = fe_values.shape_grad (i,q) / 2;
      
      for (unsigned int e=0; e<dim; ++e)
	tmp[component][e] += grad_phi_over_2[e];
      for (unsigned int d=0; d<dim; ++d)
	tmp[d][component] += grad_phi_over_2[d];
    }
  
  return tmp;
}


  
template <int dim>
double extract_p (const FEValuesBase<dim> &fe_values,
                  const unsigned int i,
                  const unsigned int q)
{
  const unsigned int component
    = fe_values.get_fe().system_to_component_index(i).first;

  if (component == dim)
    return fe_values.shape_value (i,q);
  else
    return 0;
}



template <int dim>
double extract_T (const FEValuesBase<dim> &fe_values,
                  const unsigned int i,
                  const unsigned int q)
{
  const unsigned int component
    = fe_values.get_fe().system_to_component_index(i).first;

  if (component == dim+1)
    return fe_values.shape_value (i,q);
  else
    return 0;
}



template <int dim>
Tensor<1,dim>
extract_grad_T (const FEValuesBase<dim> &fe_values,
                const unsigned int i,
                const unsigned int q)
{
  Tensor<1,dim> tmp;

  const unsigned int component
    = fe_values.get_fe().system_to_component_index(i).first;

  if (component == dim+1)
    tmp = fe_values.shape_grad (i,q);

  return tmp;
}




template <class Matrix, class Preconditioner>
class InverseMatrix : public Subscriptor
{
  public:
    InverseMatrix (const Matrix         &m,
		   const Preconditioner &preconditioner);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const Matrix> matrix;
    const Preconditioner &preconditioner;

    mutable GrowingVectorMemory<> vector_memory;    
};


template <class Matrix, class Preconditioner>
InverseMatrix<Matrix,Preconditioner>::InverseMatrix (const Matrix &m,
						     const Preconditioner &preconditioner)
                :
                matrix (&m),
		preconditioner (preconditioner)
{}


                                 
template <class Matrix, class Preconditioner>
void InverseMatrix<Matrix,Preconditioner>::vmult (Vector<double>       &dst,
						  const Vector<double> &src) const
{
  SolverControl solver_control (src.size(), 1e-6*src.l2_norm());
  SolverCG<> cg (solver_control, vector_memory);

  dst = 0;

  try
    {
      cg.solve (*matrix, dst, src, preconditioner);
    }
  catch (...)
    {
      Assert (false, ExcInternalError());
    }
}



template <class Preconditioner>
class SchurComplement : public Subscriptor
{
  public:
    SchurComplement (const BlockSparseMatrix<double> &A,
                     const InverseMatrix<SparseMatrix<double>,Preconditioner> &Minv);

    void vmult (Vector<double>       &dst,
                const Vector<double> &src) const;

  private:
    const SmartPointer<const BlockSparseMatrix<double> > system_matrix;
    const SmartPointer<const InverseMatrix<SparseMatrix<double>,Preconditioner> > m_inverse;
    
    mutable Vector<double> tmp1, tmp2;
};



template <class Preconditioner>
SchurComplement<Preconditioner>::
SchurComplement (const BlockSparseMatrix<double> &A,
                 const InverseMatrix<SparseMatrix<double>,Preconditioner> &Minv)
                :
                system_matrix (&A),
                m_inverse (&Minv),
                tmp1 (A.block(0,0).m()),
                tmp2 (A.block(0,0).m())
{}


template <class Preconditioner>
void SchurComplement<Preconditioner>::vmult (Vector<double>       &dst,
					     const Vector<double> &src) const
{
  system_matrix->block(0,1).vmult (tmp1, src);
  m_inverse->vmult (tmp2, tmp1);
  system_matrix->block(1,0).vmult (dst, tmp2);
}



template <int dim>
BoussinesqFlowProblem<dim>::BoussinesqFlowProblem (const unsigned int degree)
                :
                degree (degree),
                fe (FE_Q<dim>(degree+1), dim,
                    FE_Q<dim>(degree), 1,
                    FE_DGQ<dim>(degree-1), 1),
                dof_handler (triangulation),
                time_step (0),
		rebuild_matrices (true),
		rebuild_preconditioner (true)
{}




template <int dim>
void BoussinesqFlowProblem<dim>::setup_dofs (const bool setup_matrices)
{  
  dof_handler.distribute_dofs (fe); 
  DoFRenumbering::component_wise (dof_handler);

  hanging_node_constraints.clear ();
  DoFTools::make_hanging_node_constraints (dof_handler,
					   hanging_node_constraints);
  std::set<unsigned char> no_normal_flux_boundaries;
  no_normal_flux_boundaries.insert (0);
  compute_no_normal_flux_constraints (dof_handler, 0, no_normal_flux_boundaries,
				      hanging_node_constraints);
  hanging_node_constraints.close ();

  std::vector<unsigned int> dofs_per_component (dim+2);
  DoFTools::count_dofs_per_component (dof_handler, dofs_per_component);  
  const unsigned int n_u = dofs_per_component[0] * dim,
                     n_p = dofs_per_component[dim],
                     n_T = dofs_per_component[dim+1];

  std::cout << "Number of active cells: "
            << triangulation.n_active_cells()
            << std::endl
            << "Number of degrees of freedom: "
            << dof_handler.n_dofs()
            << " (" << n_u << '+' << n_p << '+'<< n_T <<')'
            << std::endl
            << std::endl;
  
  const unsigned int
    n_couplings = dof_handler.max_couplings_between_dofs();

  if (setup_matrices == true)
    {
      system_matrix.clear ();
      
      sparsity_pattern.reinit (3,3);
      sparsity_pattern.block(0,0).reinit (n_u, n_u, n_couplings);
      sparsity_pattern.block(1,0).reinit (n_p, n_u, n_couplings);
      sparsity_pattern.block(2,0).reinit (n_T, n_u, n_couplings);
      sparsity_pattern.block(0,1).reinit (n_u, n_p, n_couplings);
      sparsity_pattern.block(1,1).reinit (n_p, n_p, n_couplings);
      sparsity_pattern.block(2,1).reinit (n_T, n_p, n_couplings);
      sparsity_pattern.block(0,2).reinit (n_u, n_T, n_couplings);
      sparsity_pattern.block(1,2).reinit (n_p, n_T, n_couplings);
      sparsity_pattern.block(2,2).reinit (n_T, n_T, n_couplings);
  
      sparsity_pattern.collect_sizes();    

  
      DoFTools::make_sparsity_pattern (dof_handler, sparsity_pattern);
      hanging_node_constraints.condense (sparsity_pattern);
      sparsity_pattern.compress();
  
      system_matrix.reinit (sparsity_pattern);
    }
                                   
  solution.reinit (3);
  solution.block(0).reinit (n_u);
  solution.block(1).reinit (n_p);
  solution.block(2).reinit (n_T);
  solution.collect_sizes ();
  
  old_solution.reinit (3);
  old_solution.block(0).reinit (n_u);
  old_solution.block(1).reinit (n_p);
  old_solution.block(2).reinit (n_T);
  old_solution.collect_sizes ();
  
  system_rhs.reinit (3);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_p);
  system_rhs.block(2).reinit (n_T);
  system_rhs.collect_sizes ();
}


template <int dim>
double
scalar_product (const Tensor<2,dim> &a,
		const Tensor<2,dim> &b)
{
  double tmp = 0;
  for (unsigned int i=0; i<dim; ++i)
    for (unsigned int j=0; j<dim; ++j)
      tmp += a[i][j] * b[i][j];
  return tmp;
}



template <int dim>
void BoussinesqFlowProblem<dim>::assemble_system () 
{
  if (rebuild_matrices == true)
    system_matrix=0;

  system_rhs=0;
  
  QGauss<dim>   quadrature_formula(degree+2); 
  QGauss<dim-1> face_quadrature_formula(degree+2);

  FEValues<dim> fe_values (fe, quadrature_formula,
			   update_values    |
			   update_quadrature_points  |
			   update_JxW_values |
			   (rebuild_matrices == true
			    ?
			    update_gradients
			    :
			    UpdateFlags(0)));
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
                                    update_values    | update_normal_vectors |
                                    update_quadrature_points  | update_JxW_values);

  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  
  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int   n_face_q_points = face_quadrature_formula.size();

  FullMatrix<double>   local_matrix (dofs_per_cell, dofs_per_cell);
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<unsigned int> local_dof_indices (dofs_per_cell);
  
  const PressureBoundaryValues<dim> pressure_boundary_values;
  
  std::vector<double>               boundary_values (n_face_q_points);
  
  std::vector<Vector<double> >      old_solution_values(n_q_points, Vector<double>(dim+2));
  std::vector<std::vector<Tensor<1,dim> > >  old_solution_grads(n_q_points,
                                                                std::vector<Tensor<1,dim> > (dim+2));

  const double Raleigh_number = 10;
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    { 
      fe_values.reinit (cell);
      local_matrix = 0;
      local_rhs = 0;

      fe_values.get_function_values (old_solution, old_solution_values);

      for (unsigned int q=0; q<n_q_points; ++q)
	{
	  const double old_temperature = old_solution_values[q](dim+1);
	  
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    {

	      const Tensor<1,dim> phi_i_u      = extract_u (fe_values, i, q);

	      if (rebuild_matrices)
		{
		  const Tensor<2,dim> phi_i_grads_u= extract_grad_s_u (fe_values, i, q);
		  const double        div_phi_i_u  = extract_div_u (fe_values, i, q);
		  const double        phi_i_p      = extract_p (fe_values, i, q);
		  const double        phi_i_T      = extract_T (fe_values, i, q); 
		  const Tensor<1,dim> grad_phi_i_T = extract_grad_T(fe_values, i, q);
            
		  for (unsigned int j=0; j<dofs_per_cell; ++j)
		    {
		      const Tensor<2,dim> phi_j_grads_u     = extract_grad_s_u (fe_values, j, q);
		      const double        div_phi_j_u = extract_div_u (fe_values, j, q);
		      const double        phi_j_p     = extract_p (fe_values, j, q);
		      const double        phi_j_T     = extract_T (fe_values, j, q);
                
		      local_matrix(i,j) += (scalar_product(phi_i_grads_u, phi_j_grads_u)
					    - div_phi_i_u * phi_j_p
					    - phi_i_p * div_phi_j_u
					    + phi_i_p * phi_j_p
					    + phi_i_T * phi_j_T)
					   * fe_values.JxW(q);     
		    }
		}
	      
	      const Point<dim> gravity (0,1);
	      
	      local_rhs(i) += (Raleigh_number *
			       gravity * phi_i_u * old_temperature)*
			      fe_values.JxW(q);
          }
	}
      

      for (unsigned int face_no=0;
           face_no<GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        if (cell->at_boundary(face_no))
          {
            fe_face_values.reinit (cell, face_no);
            
            pressure_boundary_values
              .value_list (fe_face_values.get_quadrature_points(),
                           boundary_values);

            for (unsigned int q=0; q<n_face_q_points; ++q) 
              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  const Tensor<1,dim>
                    phi_i_u = extract_u (fe_face_values, i, q);

                  local_rhs(i) += -(phi_i_u *
                                    fe_face_values.normal_vector(q) *
                                    boundary_values[q] *
                                    fe_face_values.JxW(q));
                }
          }

      cell->get_dof_indices (local_dof_indices);

      if (rebuild_matrices == true)
	{
	  for (unsigned int i=0; i<dofs_per_cell; ++i)
	    for (unsigned int j=0; j<dofs_per_cell; ++j)
	      system_matrix.add (local_dof_indices[i],
				 local_dof_indices[j],
				 local_matrix(i,j));
	}
      
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += local_rhs(i);
    }

  if (rebuild_matrices == true)
    hanging_node_constraints.condense (system_matrix);
  
  hanging_node_constraints.condense (system_rhs);  
  
  if (rebuild_matrices == true)
    {
//       std::map<unsigned int,double> boundary_values;

//       typename DoFHandler<dim>::active_cell_iterator
// 	cell = dof_handler.begin_active(),
// 	emdc = dof_handler.end();
//       for (; cell!=endc; ++cell)
// 	for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_cell; ++v)
// 	  if (cell->vertex(v).distance(dim == 2
// 				       ?
// 				       Point<dim>(0,-1)
// 				       :
// 				       Point<dim>(0,0,-1)) < 1e-6)
// 	    {
// 	      std::cout << "Found cell and vertex: " << cell << ' '
// 			<< v << std::endl;
	      
// 	      boundary_values[cell->vertex_dof_index(v,0)] = 0;
// 	      break;
// 	    } 
      
//      std::vector<bool> component_mask (dim+2, true);
//       component_mask[dim] = component_mask[dim+1] = false;
//       VectorTools::interpolate_boundary_values (dof_handler,
// 						0,
// 						ZeroFunction<dim>(dim+2),
// 						boundary_values,
// 						component_mask);

//       MatrixTools::apply_boundary_values (boundary_values,
// 					  system_matrix,
// 					  solution,
// 					  system_rhs);  
    }

  if (rebuild_preconditioner == true)
    {
      Assert (rebuild_matrices == true,
	      ExcMessage ("There is no point in rebuilding the preconditioner "
			  "without a rebuilt matrix!"));
      
      std::cout << "   Rebuilding preconditioner..." << std::flush;
      
      A_preconditioner
	= boost::shared_ptr<SparseDirectUMFPACK>(new SparseDirectUMFPACK());
      A_preconditioner->initialize (system_matrix.block(0,0));

      std::cout << std::endl;
      
      rebuild_preconditioner = false;
    }

  rebuild_matrices = false;
}






template <int dim>
void BoussinesqFlowProblem<dim>::assemble_rhs_T () 
{  
  QGauss<dim>   quadrature_formula(degree+2); 
  QGauss<dim-1> face_quadrature_formula(degree+2);
  FEValues<dim> fe_values (fe, quadrature_formula, 
                           update_values    | update_gradients |
                           update_quadrature_points  | update_JxW_values);
  FEFaceValues<dim> fe_face_values (fe, face_quadrature_formula, 
                                    update_values    | update_normal_vectors |
                                    update_quadrature_points  | update_JxW_values);
  FESubfaceValues<dim> fe_subface_values (fe, face_quadrature_formula, 
					  update_values    | update_normal_vectors |
					  update_JxW_values);
  FEFaceValues<dim> fe_face_values_neighbor (fe, face_quadrature_formula, 
                                             update_values);
  FESubfaceValues<dim> fe_subface_values_neighbor (fe, face_quadrature_formula, 
						   update_values);
 
  const unsigned int   dofs_per_cell   = fe.dofs_per_cell;
  const unsigned int   n_q_points      = quadrature_formula.size();
  const unsigned int   n_face_q_points = face_quadrature_formula.size();
  
  Vector<double>       local_rhs (dofs_per_cell);

  std::vector<Vector<double> > old_solution_values(n_q_points, Vector<double>(dim+2));

  std::vector<Vector<double> > old_solution_values_face(n_face_q_points, Vector<double>(dim+2));
  std::vector<Vector<double> > old_solution_values_face_neighbor(n_face_q_points, Vector<double>(dim+2));
  std::vector<Vector<double> > present_solution_values(n_q_points, Vector<double>(dim+2));
  std::vector<Vector<double> > present_solution_values_face(n_face_q_points, Vector<double>(dim+2));

  std::vector<std::vector<Tensor<1,dim> > >
    present_solution_grads(n_q_points,
			   std::vector<Tensor<1,dim> >(dim+2));
	  
  
  std::vector<double> neighbor_temperature (n_face_q_points);
  std::vector<unsigned int> local_dof_indices (dofs_per_cell);

  TemperatureBoundaryValues<dim> temperature_boundary_values;
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

  for (; cell!=endc; ++cell)
    {
      local_rhs = 0;
      fe_values.reinit (cell);

      fe_values.get_function_values (old_solution, old_solution_values);
      fe_values.get_function_values (solution, present_solution_values);
      fe_values.get_function_gradients (solution, present_solution_grads);

      for (unsigned int q=0; q<n_q_points; ++q) 
        for (unsigned int i=0; i<dofs_per_cell; ++i)
          {
            const double old_T = old_solution_values[q](dim+1);
            Tensor<1,dim> present_u;
            for (unsigned int d=0; d<dim; ++d)
              present_u[d] = present_solution_values[q](d);

	    double present_div_u = 0;
            for (unsigned int d=0; d<dim; ++d)
              present_div_u += present_solution_grads[q][d][d];	    
	    
            const double        phi_i_T      = extract_T(fe_values, i, q);
            const Tensor<1,dim> grad_phi_i_T = extract_grad_T(fe_values, i, q);

	    const Point<dim> p = fe_values.quadrature_point(q);
	    
            local_rhs(i) += (time_step *
                             old_T *
                             (present_u *
			      grad_phi_i_T
			      +
			      present_div_u *
			      phi_i_T)
                             +
                             old_T * phi_i_T
			     +
			     time_step *
			     RightHandSide<dim>().value (p, dim+1)
			     * phi_i_T)
                            *
                            fe_values.JxW(q);
          }

      
//TODO: unify the code that actually does the assembly down below      
      for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
           ++face_no)
        if (cell->at_boundary(face_no)
	    ||
	    ((cell->neighbor(face_no)->has_children() == false)
	     &&
	     (cell->neighbor(face_no)->level() == cell->level())))
	  {
					     // cell either at
					     // boundary or with a
					     // neighbor that has the
					     // same refinement level
					     // and is not further
					     // refined
	    fe_face_values.reinit (cell, face_no);

	    fe_face_values.get_function_values (old_solution,
						old_solution_values_face);
	    fe_face_values.get_function_values (solution,
						present_solution_values_face);

	    if (cell->at_boundary(face_no))
	      temperature_boundary_values
		.value_list (fe_face_values.get_quadrature_points(),
			     neighbor_temperature);
	    else
	      {
		const typename DoFHandler<dim>::active_cell_iterator
		  neighbor = cell->neighbor(face_no);
		
		fe_face_values_neighbor.reinit (neighbor,
						cell->neighbor_of_neighbor(face_no));
             
		fe_face_values_neighbor
		  .get_function_values (old_solution,
					old_solution_values_face_neighbor);
             
		for (unsigned int q=0; q<n_face_q_points; ++q)
		  neighbor_temperature[q] = old_solution_values_face_neighbor[q](dim+1);
	      }

	    for (unsigned int q=0; q<n_face_q_points; ++q)
	      {
		Tensor<1,dim> present_u_face;
		for (unsigned int d=0; d<dim; ++d)
		  present_u_face[d] = present_solution_values_face[q](d);

		const double normal_flux = present_u_face *
					   fe_face_values.normal_vector(q);

		const bool is_outflow_q_point = (normal_flux >= 0);
                                     
		for (unsigned int i=0; i<dofs_per_cell; ++i)
		  local_rhs(i) -= time_step *
				  normal_flux *
				  (is_outflow_q_point == true
				   ?
				   old_solution_values_face[q](dim+1)
				   :
				   neighbor_temperature[q]) *
				  extract_T(fe_face_values,i,q) *
				  fe_face_values.JxW(q);
	      }
	  }
	else
	  if (cell->neighbor(face_no)->has_children())
	    {
					       // neighbor is further
					       // refined. loop over
					       // all sub faces
	      for (unsigned int subface_no=0;
		   subface_no<GeometryInfo<dim>::subfaces_per_face;
		   ++subface_no)
		{
		  fe_subface_values.reinit (cell, face_no, subface_no);

		  fe_subface_values.get_function_values (old_solution,
							 old_solution_values_face);
		  fe_subface_values.get_function_values (solution,
							 present_solution_values_face);

		  const typename DoFHandler<dim>::active_cell_iterator
		    neighbor = cell->neighbor_child_on_subface (face_no, subface_no);

		  fe_face_values_neighbor.reinit (neighbor,
						  cell->neighbor_of_neighbor(face_no));
             
		  fe_face_values_neighbor
		    .get_function_values (old_solution,
					  old_solution_values_face_neighbor);
		  
		  for (unsigned int q=0; q<n_face_q_points; ++q)
		    neighbor_temperature[q] = old_solution_values_face_neighbor[q](dim+1);

		  for (unsigned int q=0; q<n_face_q_points; ++q)
		    {
		      Tensor<1,dim> present_u_face;
		      for (unsigned int d=0; d<dim; ++d)
			present_u_face[d] = present_solution_values_face[q](d);

		      const double normal_flux = present_u_face *
						 fe_subface_values.normal_vector(q);

		      const bool is_outflow_q_point = (normal_flux >= 0);
                                     
		      for (unsigned int i=0; i<dofs_per_cell; ++i)
			local_rhs(i) -= time_step *
					normal_flux *
					(is_outflow_q_point == true
					 ?
					 old_solution_values_face[q](dim+1)
					 :
					 neighbor_temperature[q]) *
					extract_T(fe_face_values,i,q) *
					fe_face_values.JxW(q);
		    }
		}
	    }
	  else
	    {
					       // neighbor is less
					       // refined. we need to
					       // use a subface values
					       // object for the
					       // neighbor's subface
	      fe_face_values.reinit (cell, face_no);

	      fe_face_values.get_function_values (old_solution, old_solution_values_face);
	      fe_face_values.get_function_values (solution, present_solution_values_face);

	      const typename DoFHandler<dim>::active_cell_iterator
		neighbor = cell->neighbor (face_no);

	      const std::pair<unsigned int, unsigned int> faceno_subfaceno=
		cell->neighbor_of_coarser_neighbor(face_no);
	      const unsigned int neighbor_face_no    = faceno_subfaceno.first,
				 neighbor_subface_no = faceno_subfaceno.second;

	      Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
							   neighbor_subface_no)
		      == cell,
		      ExcInternalError());

	      fe_subface_values_neighbor.reinit (neighbor,
						 neighbor_face_no,
						 neighbor_subface_no);
	      
	      fe_subface_values_neighbor
		.get_function_values (old_solution,
				      old_solution_values_face_neighbor);
             
	      for (unsigned int q=0; q<n_face_q_points; ++q)
		neighbor_temperature[q] = old_solution_values_face_neighbor[q](dim+1);

	      for (unsigned int q=0; q<n_face_q_points; ++q)
		{
		  Tensor<1,dim> present_u_face;
		  for (unsigned int d=0; d<dim; ++d)
		    present_u_face[d] = present_solution_values_face[q](d);

		  const double normal_flux = present_u_face *
					     fe_face_values.normal_vector(q);

		  const bool is_outflow_q_point = (normal_flux >= 0);
                                     
		  for (unsigned int i=0; i<dofs_per_cell; ++i)
		    local_rhs(i) -= time_step *
				    normal_flux *
				    (is_outflow_q_point == true
				     ?
				     old_solution_values_face[q](dim+1)
				     :
				     neighbor_temperature[q]) *
				    extract_T(fe_face_values,i,q) *
				    fe_face_values.JxW(q);
		}
	    }      
      
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        system_rhs(local_dof_indices[i]) += local_rhs(i);       
    }
} 




template <int dim>
void BoussinesqFlowProblem<dim>::solve () 
{
  solution = old_solution;
  
  const InverseMatrix<SparseMatrix<double>,SparseDirectUMFPACK>
    A_inverse (system_matrix.block(0,0), *A_preconditioner);
  Vector<double> tmp (solution.block(0).size());
  Vector<double> schur_rhs (solution.block(1).size());
  Vector<double> tmp2 (solution.block(2).size());
  
  {
    A_inverse.vmult (tmp, system_rhs.block(0));
    system_matrix.block(1,0).vmult (schur_rhs, tmp);
    schur_rhs -= system_rhs.block(1);

    
    SchurComplement<SparseDirectUMFPACK>
      schur_complement (system_matrix, A_inverse);
    
    SolverControl solver_control (system_matrix.block(0,0).m(),
                                  1e-6*schur_rhs.l2_norm());
    SolverCG<>    cg (solver_control);
    
    PreconditionSSOR<> preconditioner;
    preconditioner.initialize (system_matrix.block(1,1), 1.2);

    InverseMatrix<SparseMatrix<double>,PreconditionSSOR<> >
      m_inverse (system_matrix.block(1,1), preconditioner);
    
    try
      {
	cg.solve (schur_complement, solution.block(1), schur_rhs,
		  m_inverse);
      }
    catch (...)
      {
	abort ();
      }

				     // produce a consistent flow field
    hanging_node_constraints.distribute (solution);
  
    std::cout << "   "
              << solver_control.last_step()
              << " CG Schur complement iterations for pressure."
              << std::endl;    
  }

  {
    system_matrix.block(0,1).vmult (tmp, solution.block(1));
    tmp *= -1;
    tmp += system_rhs.block(0);
    
    A_inverse.vmult (solution.block(0), tmp);

				     // produce a consistent pressure field
    hanging_node_constraints.distribute (solution);
  }

				   // for DGQ1 needs to be /15
  time_step = GridTools::minimal_cell_diameter(triangulation) /
              std::max (get_maximal_velocity(), .05) / 2;

  assemble_rhs_T ();
  {
    
    SolverControl solver_control (system_matrix.block(2,2).m(),
                                  1e-8*system_rhs.block(2).l2_norm());
    SolverCG<>   cg (solver_control);
    PreconditionJacobi<> preconditioner;
    preconditioner.initialize (system_matrix.block(2,2));
    
    try
      {
	cg.solve (system_matrix.block(2,2), solution.block(2), system_rhs.block(2),
		  preconditioner);
      }
    catch (...)
      {
	abort ();
      }
	
				     // produce a consistent temperature field
    hanging_node_constraints.distribute (solution);
                
    std::cout << "   "
              << solver_control.last_step()
              << " CG iterations for temperature."
              << std::endl;
    std::cout << "   Max temperature: "
	      << *std::max_element (solution.block(2).begin(),
				    solution.block(2).end())
	      << std::endl;
  } 
}
                                 


template <int dim>
void BoussinesqFlowProblem<dim>::output_results ()  const
{
  if (timestep_number % 10 != 0)
    return;
  
  std::vector<std::string> solution_names (dim, "velocity");
  solution_names.push_back ("p");
  solution_names.push_back ("T");
  
  DataOut<dim> data_out;

  data_out.attach_dof_handler (dof_handler);

  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation
    (dim+2, DataComponentInterpretation::component_is_scalar);
  for (unsigned int i=0; i<dim; ++i)
    data_component_interpretation[i]
      = DataComponentInterpretation::component_is_part_of_vector;
  
  data_out.add_data_vector (solution, solution_names,
			    DataOut<dim>::type_dof_data,
			    data_component_interpretation);
  
  data_out.build_patches ();
  
  std::ostringstream filename;
  filename << "solution-" << Utilities::int_to_string(timestep_number, 4) << ".vtk";

  std::ofstream output (filename.str().c_str());
  data_out.write_vtk (output);
}



template <int dim>
void
BoussinesqFlowProblem<dim>::refine_mesh () 
{
  Vector<float> estimated_error_per_cell (triangulation.n_active_cells());

//TODO do this better      
  DerivativeApproximation::approximate_gradient (dof_handler,
						 old_solution,
						 estimated_error_per_cell,
						 dim+1);

  typename Triangulation<dim>::active_cell_iterator
    cell = triangulation.begin_active(),
    endc = triangulation.end();
  for (unsigned int cell_index=0; cell!=endc; ++cell, ++cell_index)
    estimated_error_per_cell(cell_index) *= cell->diameter();

  GridRefinement::refine_and_coarsen_fixed_fraction (triangulation,
						     estimated_error_per_cell,
						     0.3, 0.03,
						     static_cast<unsigned int>
						     (triangulation.n_active_cells()*1.1));

  SolutionTransfer<dim, double> soltrans(dof_handler);

  triangulation.prepare_coarsening_and_refinement();

  Vector<double> x_old_solution (dof_handler.n_dofs());
  x_old_solution = old_solution;
  
  soltrans.prepare_for_coarsening_and_refinement(x_old_solution);

  triangulation.execute_coarsening_and_refinement ();
  setup_dofs (true);

  Vector<double> tmp (dof_handler.n_dofs());
  soltrans.interpolate(x_old_solution, tmp);

  rebuild_matrices       = true;
  rebuild_preconditioner = true;

  old_solution = tmp;
}



template <int dim>
double
BoussinesqFlowProblem<dim>::get_maximal_velocity () const
{
  QGauss<dim>   quadrature_formula(degree+2); 
  const unsigned int   n_q_points
    = quadrature_formula.size();

  FEValues<dim> fe_values (fe, quadrature_formula, 
                           update_values);
  std::vector<Vector<double> > solution_values(n_q_points,
                                               Vector<double>(dim+2));
  double max_velocity = 0;
  
  typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      fe_values.get_function_values (solution, solution_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          Tensor<1,dim> velocity;
          for (unsigned int i=0; i<dim; ++i)
            velocity[i] = solution_values[q](i);          
          
          max_velocity = std::max (max_velocity,
                                   velocity.norm());
        }
    }
  
  return max_velocity;
}



template <int dim>
void BoussinesqFlowProblem<dim>::run () 
{
  switch (dim)
    {
      case 2:
      {
	GridGenerator::hyper_cube (triangulation);
	
	triangulation.refine_global (5);

	break;
      }

      case 3:
      {
	GridGenerator::hyper_shell (triangulation,
				    Point<dim>(), 0.5, 1.0);
	
	static HyperShellBoundary<dim> boundary;
	triangulation.set_boundary (0, boundary);
	
	triangulation.refine_global (2);

	break;
      }

      default:
	    Assert (false, ExcNotImplemented());
    }
  

  const bool do_adaptivity = false;

  if (do_adaptivity)
    {
      setup_dofs(false);

      VectorTools::project (dof_handler,
			    hanging_node_constraints,
			    QGauss<dim>(degree+2),
			    InitialValues<dim>(),
			    old_solution);
  
      for (unsigned int pre_refinement=0; pre_refinement<4-dim; ++pre_refinement)
	{
	  refine_mesh ();

	  VectorTools::project (dof_handler,
				hanging_node_constraints,
				QGauss<dim>(degree+2),
				InitialValues<dim>(),
				old_solution);
	}
    }
  else
    {
      setup_dofs(true);

      VectorTools::project (dof_handler,
			    hanging_node_constraints,
			    QGauss<dim>(degree+2),
			    InitialValues<dim>(),
			    old_solution);
    }
  
  timestep_number = 0;
  double time = 0;
  
  do
    { 
      std::cout << "Timestep " << timestep_number
		<< ":  t=" << time
		<< ", dt=" << time_step
                << std::endl; 

      std::cout << "   Assembling..." << std::endl;
      assemble_system ();      

      std::cout << "   Solving..." << std::endl;
      solve ();
      
      output_results ();

      time += time_step;
      ++timestep_number;
   
      old_solution = solution; 

      std::cout << std::endl;

      if (do_adaptivity)
	if (timestep_number % 10 == 0)
	  refine_mesh ();
    }
  while (time <= 500);
}

    

int main () 
{
  try
    {
      deallog.depth_console (0);

      BoussinesqFlowProblem<2> flow_problem(1);
      flow_problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      
      return 1;
    }
  catch (...) 
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
