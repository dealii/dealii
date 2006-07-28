//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/exceptions.h>
#include <base/memory_consumption.h>
#include <dofs/dof_objects.h>
#include <dofs/dof_handler.h>
#include <fe/fe.h>

namespace internal
{
  namespace DoFHandler
  {
    template<int dim>
    unsigned int
    DoFObjects<dim>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (dofs));
    }



    template <int dim>
    template <int spacedim>
    unsigned int
    DoFObjects<dim>::
    get_dof_index (const ::DoFHandler<spacedim> &dof_handler,
		   const unsigned int       obj_index,
		   const unsigned int       fe_index,
		   const unsigned int       local_index) const
    {
      Assert (fe_index == ::DoFHandler<spacedim>::default_fe_index,
	      ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
      Assert (local_index<dof_handler.get_fe().template n_dofs_per_object<dim>(),
	      ExcIndexRange (local_index, 0, dof_handler.get_fe().template n_dofs_per_object<dim>()));
      Assert (obj_index * dof_handler.get_fe().template n_dofs_per_object<dim>()+local_index
	      <
	      dofs.size(),
	      ExcInternalError());
      
      return dofs[obj_index * dof_handler.get_fe().template n_dofs_per_object<dim>() + local_index];
    }

    

    template <int dim>
    template <int spacedim>
    void
    DoFObjects<dim>::
    set_dof_index (const ::DoFHandler<spacedim> &dof_handler,
		   const unsigned int       obj_index,
		   const unsigned int       fe_index,
		   const unsigned int       local_index,
		   const unsigned int       global_index)
    {
      Assert (fe_index == ::DoFHandler<spacedim>::default_fe_index,
	      ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
      Assert (local_index<dof_handler.get_fe().template n_dofs_per_object<dim>(),
	      ExcIndexRange (local_index, 0, dof_handler.get_fe().template n_dofs_per_object<dim>()));
      Assert (obj_index * dof_handler.get_fe().template n_dofs_per_object<dim>()+local_index
	      <
	      dofs.size(),
	      ExcInternalError());

      dofs[obj_index * dof_handler.get_fe().template n_dofs_per_object<dim>() + local_index] = global_index;
    }

// explicit instantiations

    template class DoFObjects<1>;

    template
    unsigned int
    DoFObjects<1>::
    get_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
		   const unsigned int       obj_index,
		   const unsigned int       fe_index,
		   const unsigned int       local_index) const;
    
    template
    void
    DoFObjects<1>::
    set_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
		   const unsigned int       obj_index,
		   const unsigned int       fe_index,
		   const unsigned int       local_index,
		   const unsigned int       global_index);

#if deal_II_dimension >= 2

    template class DoFObjects<2>;

    template
    unsigned int
    DoFObjects<2>::
    get_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
		   const unsigned int       obj_index,
		   const unsigned int       fe_index,
		   const unsigned int       local_index) const;
    
    template
    void
    DoFObjects<2>::
    set_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
		   const unsigned int       obj_index,
		   const unsigned int       fe_index,
		   const unsigned int       local_index,
		   const unsigned int       global_index);

#endif

#if deal_II_dimension >= 3
    
    template class DoFObjects<3>;

    template
    unsigned int
    DoFObjects<3>::
    get_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
		   const unsigned int       obj_index,
		   const unsigned int       fe_index,
		   const unsigned int       local_index) const;
    
    template
    void
    DoFObjects<3>::
    set_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
		   const unsigned int       obj_index,
		   const unsigned int       fe_index,
		   const unsigned int       local_index,
		   const unsigned int       global_index);

#endif

  }
}
