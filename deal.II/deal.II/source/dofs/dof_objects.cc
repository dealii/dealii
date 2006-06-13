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
    inline
    unsigned int
    DoFObjects<dim>::n_active_fe_indices (const ::DoFHandler<spacedim> &,
						    const unsigned) const
    {
      return 1;
    }
    

    template <int dim>
    template <int spacedim>
    inline
    bool
    DoFObjects<dim>::fe_index_is_active (const ::DoFHandler<spacedim> &,
						   const unsigned int,
						   const unsigned int fe_index) const
    {
      Assert (fe_index == 0,
              ExcMessage ("Only zero fe_index values are allowed for "
                          "non-hp DoFHandlers."));
      return true;
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
      unsigned int dofs_per_obj;
      switch (dim)
	{
	  case 1 :
		dofs_per_obj = dof_handler.get_fe().dofs_per_line;
		break;
	  case 2 :
		dofs_per_obj = dof_handler.get_fe().dofs_per_quad;
		break;
	  case 3 :
		dofs_per_obj = dof_handler.get_fe().dofs_per_hex;
	}

      Assert (fe_index == ::DoFHandler<spacedim>::default_fe_index,
	      ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
      Assert (local_index<dofs_per_obj,
	      ExcIndexRange (local_index, 0, dofs_per_obj));
      Assert (obj_index * dofs_per_obj+local_index
	      <
	      dofs.size(),
	      ExcInternalError());
      
      return dofs[obj_index * dofs_per_obj + local_index];
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
      unsigned int dofs_per_obj;
      switch (dim)
	{
	  case 1 :
		dofs_per_obj = dof_handler.get_fe().dofs_per_line;
		break;
	  case 2 :
		dofs_per_obj = dof_handler.get_fe().dofs_per_quad;
		break;
	  case 3 :
		dofs_per_obj = dof_handler.get_fe().dofs_per_hex;
	}
            
      Assert (fe_index == ::DoFHandler<spacedim>::default_fe_index,
	      ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
      Assert (local_index<dofs_per_obj,
	      ExcIndexRange (local_index, 0, dof_handler.get_fe().dofs_per_line));
      Assert (obj_index * dofs_per_obj+local_index
	      <
	      dofs.size(),
	      ExcInternalError());

      dofs[obj_index * dofs_per_obj + local_index] = global_index;
    }

// explicit instantiations
    template
    unsigned int
    DoFObjects<1>::
    memory_consumption () const;
    

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

    template 
    unsigned int
    DoFObjects<1>::
    n_active_fe_indices (const ::DoFHandler<deal_II_dimension> &,
			 const unsigned) const;

    template
    bool
    DoFObjects<1>::
    fe_index_is_active (const ::DoFHandler<deal_II_dimension> &,
			const unsigned int,
			const unsigned int fe_index) const;

#if deal_II_dimension >= 2

    template
    unsigned int
    DoFObjects<2>::
    memory_consumption () const;
    

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

    template 
    unsigned int
    DoFObjects<2>::
    n_active_fe_indices (const ::DoFHandler<deal_II_dimension> &,
			 const unsigned) const;

    template
    bool
    DoFObjects<2>::
    fe_index_is_active (const ::DoFHandler<deal_II_dimension> &,
			const unsigned int,
			const unsigned int fe_index) const;
    
#endif

#if deal_II_dimension >= 3
    
    template
    unsigned int
    DoFObjects<3>::
    memory_consumption () const;
    

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

    template 
    unsigned int
    DoFObjects<3>::
    n_active_fe_indices (const ::DoFHandler<deal_II_dimension> &,
			 const unsigned) const;

    template
    bool
    DoFObjects<3>::
    fe_index_is_active (const ::DoFHandler<deal_II_dimension> &,
			const unsigned int,
			const unsigned int fe_index) const;

#endif

  }
}
