//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/config.h>

#ifdef DEAL_II_USE_P4EST

#include <lac/vector.h>
#include <lac/petsc_vector.h>
#include <lac/trilinos_vector.h>
#include <lac/block_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_block_vector.h>

#include <distributed/solution_transfer.h>
#include <distributed/tria.h>
#include <dofs/dof_tools.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>

#include <base/std_cxx1x/bind.h>

DEAL_II_NAMESPACE_OPEN

namespace parallel
{
  namespace distributed
  {

    namespace
    {
      template <class VECTOR>
      void compress_vector_insert(VECTOR & vec)
      {
	vec.compress();
      }

#ifdef DEAL_II_USE_TRILINOS
      void compress_vector_insert(TrilinosWrappers::Vector & vec)
      {
	vec.compress(Insert);
      }

      void compress_vector_insert(TrilinosWrappers::BlockVector & vec)
      {
	for (unsigned int i=0;i<vec.n_blocks();++i)
	  vec.block(i).compress(Insert);
      }

      void compress_vector_insert(TrilinosWrappers::MPI::Vector & vec)
      {
	vec.compress(Insert);
      }

      void compress_vector_insert(TrilinosWrappers::MPI::BlockVector & vec)
      {
	for (unsigned int i=0;i<vec.n_blocks();++i)
	  vec.block(i).compress(Insert);
      }
#endif
    }



    template<int dim, typename VECTOR, class DH>
    SolutionTransfer<dim, VECTOR, DH>::SolutionTransfer(const DH &dof)
		    :
		    dof_handler(&dof, typeid(*this).name())
    {}



    template<int dim, typename VECTOR, class DH>
    SolutionTransfer<dim, VECTOR, DH>::~SolutionTransfer()
    {}



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    prepare_for_coarsening_and_refinement (const std::vector<const VECTOR*> &all_in)
    {
      Assert(all_in.size() > 0, ExcMessage("Please transfer atleast one vector!"));
      input_vectors = all_in;
      SolutionTransfer<dim, VECTOR, DH> *ptr = this;

//TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim> * tria
	= (dynamic_cast<parallel::distributed::Triangulation<dim>*>
	   (const_cast<dealii::Triangulation<dim>*>
	    (&dof_handler->get_tria())));
      Assert (tria != 0, ExcInternalError());

      offset
	= tria->register_data_attach(static_cast<size_t>
				     (get_data_size() * input_vectors.size()),
				     std_cxx1x::bind(&SolutionTransfer<dim, VECTOR, DH>::pack_callback,
						     ptr,
						     _1,
						     _2,
						     _3));
    }



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    prepare_for_coarsening_and_refinement (const VECTOR &in)
    {
      std::vector<const VECTOR*> all_in(1, &in);
      prepare_for_coarsening_and_refinement(all_in);
    }



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    interpolate (std::vector<VECTOR*> &all_out)
    {
      Assert(input_vectors.size()==all_out.size(),
	     ExcDimensionMismatch(input_vectors.size(), all_out.size()) );

//TODO: casting away constness is bad
      parallel::distributed::Triangulation<dim> * tria
	= (dynamic_cast<parallel::distributed::Triangulation<dim>*>
	   (const_cast<dealii::Triangulation<dim>*>
	    (&dof_handler->get_tria())));
      Assert (tria != 0, ExcInternalError());

      tria->notify_ready_to_unpack(offset,
				   std_cxx1x::bind(&SolutionTransfer<dim, VECTOR, DH>::unpack_callback,
						   this,
						   _1,
						   _2,
						   _3,
						   std_cxx1x::ref(all_out)));

      for (typename std::vector<VECTOR*>::iterator it=all_out.begin();
	   it !=all_out.end();
	   ++it)
	compress_vector_insert(*(*it));

      input_vectors.clear();
    }



    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    interpolate (VECTOR &out)
    {
      std::vector<VECTOR*> all_out(1, &out);
      interpolate(all_out);
    }



    template<int dim, typename VECTOR, class DH>
    unsigned int
    SolutionTransfer<dim, VECTOR, DH>::
    get_data_size() const
    {
      return sizeof(double)* DoFTools::max_dofs_per_cell(*dof_handler);
    }


    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    pack_callback(const typename Triangulation<dim,dim>::cell_iterator & cell_,
		  const typename Triangulation<dim,dim>::CellStatus /*status*/,
		  void* data)
    {
      double *data_store = reinterpret_cast<double *>(data);

      typename DH::cell_iterator cell(&dof_handler->get_tria(), cell_->level(), cell_->index(), dof_handler);

      const unsigned int dofs_per_cell=cell->get_fe().dofs_per_cell;
      Vector<double> dofvalues(dofs_per_cell);
      for (typename std::vector<const VECTOR*>::iterator it=input_vectors.begin();
	   it !=input_vectors.end();
	   ++it)
	{
	  cell->get_interpolated_dof_values(*(*it), dofvalues);
	  std::memcpy(data_store, &dofvalues(0), sizeof(double)*dofs_per_cell);
	  data_store += dofs_per_cell;
	}
    }


    template<int dim, typename VECTOR, class DH>
    void
    SolutionTransfer<dim, VECTOR, DH>::
    unpack_callback(const typename Triangulation<dim,dim>::cell_iterator & cell_,
		    const typename Triangulation<dim,dim>::CellStatus /*status*/,
		    const void* data,
		    std::vector<VECTOR*> &all_out)
    {
      typename DH::cell_iterator
	cell(&dof_handler->get_tria(), cell_->level(), cell_->index(), dof_handler);

      const unsigned int dofs_per_cell=cell->get_fe().dofs_per_cell;
      Vector<double> dofvalues(dofs_per_cell);
      const double *data_store = reinterpret_cast<const double *>(data);

      for (typename std::vector<VECTOR*>::iterator it = all_out.begin();
	   it != all_out.end();
	   ++it)
	{
	  std::memcpy(&dofvalues(0), data_store, sizeof(double)*dofs_per_cell);
	  cell->set_dof_values_by_interpolation(dofvalues, *(*it));
	  data_store += dofs_per_cell;
	}
    }


  }
}


// explicit instantiations

namespace parallel
{
  namespace distributed
  {
#if deal_II_dimension > 1
    template class SolutionTransfer<deal_II_dimension,Vector<double>, DoFHandler<deal_II_dimension> >;


#ifdef DEAL_II_USE_PETSC
    template class SolutionTransfer<deal_II_dimension, PETScWrappers::Vector, DoFHandler<deal_II_dimension> >;

    template class SolutionTransfer<deal_II_dimension, PETScWrappers::BlockVector, DoFHandler<deal_II_dimension> >;

    template class SolutionTransfer<deal_II_dimension, PETScWrappers::MPI::Vector, DoFHandler<deal_II_dimension> >;

    template class SolutionTransfer<deal_II_dimension, PETScWrappers::MPI::BlockVector, DoFHandler<deal_II_dimension> >;

#endif

#ifdef DEAL_II_USE_TRILINOS
    template class SolutionTransfer<deal_II_dimension, TrilinosWrappers::Vector, DoFHandler<deal_II_dimension> >;

    template class SolutionTransfer<deal_II_dimension, TrilinosWrappers::BlockVector, DoFHandler<deal_II_dimension> >;

    template class SolutionTransfer<deal_II_dimension, TrilinosWrappers::MPI::Vector, DoFHandler<deal_II_dimension> >;

    template class SolutionTransfer<deal_II_dimension, TrilinosWrappers::MPI::BlockVector, DoFHandler<deal_II_dimension> >;
#endif
#endif
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
