//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/exceptions.h>
#include <base/memory_consumption.h>
#include <dofs/dof_levels.h>
#include <dofs/dof_handler.h>
#include <fe/fe.h>


namespace internal
{
  namespace DoFHandler
  {
    template <int dim>
    unsigned int
    DoFLevel<1>::
    get_line_dof_index (const ::DoFHandler<dim> &dof_handler,
                        const unsigned int       line_index,
                        const unsigned int       fe_index,
                        const unsigned int       local_index) const
    {
      Assert (fe_index == ::DoFHandler<dim>::default_fe_index,
	      ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
      Assert (local_index<dof_handler.get_fe().dofs_per_line,
	      ExcIndexRange (local_index, 0, dof_handler.get_fe().dofs_per_line));
      Assert (line_index * dof_handler.get_fe().dofs_per_line+local_index
	      <
	      line_dofs.size(),
	      ExcInternalError());
      
      return line_dofs[line_index * dof_handler.get_fe().dofs_per_line+local_index];
    }



    template <int dim>
    void
    DoFLevel<1>::
    set_line_dof_index (const ::DoFHandler<dim> &dof_handler,
                        const unsigned int       line_index,
                        const unsigned int       fe_index,
                        const unsigned int       local_index,
                        const unsigned int       global_index)
    {
      Assert (fe_index == ::DoFHandler<dim>::default_fe_index,
	      ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
      Assert (local_index<dof_handler.get_fe().dofs_per_line,
	      ExcIndexRange (local_index, 0, dof_handler.get_fe().dofs_per_line));
      Assert (line_index * dof_handler.get_fe().dofs_per_line+local_index
	      <
	      line_dofs.size(),
	      ExcInternalError());
      
      line_dofs[line_index * dof_handler.get_fe().dofs_per_line+local_index] = global_index;
    }



    template <int dim>
    unsigned int
    DoFLevel<2>::
    get_quad_dof_index (const ::DoFHandler<dim> &dof_handler,
                        const unsigned int       quad_index,
                        const unsigned int       fe_index,
                        const unsigned int       local_index) const
    {
      Assert (fe_index == ::DoFHandler<dim>::default_fe_index,
	      ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
      Assert (local_index<dof_handler.get_fe().dofs_per_quad,
	      ExcIndexRange (local_index, 0, dof_handler.get_fe().dofs_per_quad));
      Assert (quad_index * dof_handler.get_fe().dofs_per_quad+local_index
	      <
	      quad_dofs.size(),
	      ExcInternalError());
      
      return quad_dofs[quad_index * dof_handler.get_fe().dofs_per_quad+local_index];
    }



    template <int dim>
    void
    DoFLevel<2>::
    set_quad_dof_index (const ::DoFHandler<dim> &dof_handler,
                        const unsigned int       quad_index,
                        const unsigned int       fe_index,
                        const unsigned int       local_index,
                        const unsigned int       global_index)
    {
      Assert (fe_index == ::DoFHandler<dim>::default_fe_index,
	      ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
      Assert (local_index<dof_handler.get_fe().dofs_per_quad,
	      ExcIndexRange (local_index, 0, dof_handler.get_fe().dofs_per_quad));
      Assert (quad_index * dof_handler.get_fe().dofs_per_quad+local_index
	      <
	      quad_dofs.size(),
	      ExcInternalError());
      
      quad_dofs[quad_index * dof_handler.get_fe().dofs_per_quad+local_index] = global_index;
    }



    template <int dim>
    unsigned int
    DoFLevel<3>::
    get_hex_dof_index (const ::DoFHandler<dim> &dof_handler,
                       const unsigned int       hex_index,
                       const unsigned int       fe_index,
                       const unsigned int       local_index) const
    {
      Assert (fe_index == ::DoFHandler<dim>::default_fe_index,
	      ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
      Assert (local_index<dof_handler.get_fe().dofs_per_hex,
	      ExcIndexRange (local_index, 0, dof_handler.get_fe().dofs_per_hex));
      Assert (hex_index * dof_handler.get_fe().dofs_per_hex+local_index
	      <
	      hex_dofs.size(),
	      ExcInternalError());
      
      return hex_dofs[hex_index * dof_handler.get_fe().dofs_per_hex+local_index];
    }



    template <int dim>
    void
    DoFLevel<3>::
    set_hex_dof_index (const ::DoFHandler<dim> &dof_handler,
                       const unsigned int       hex_index,
                       const unsigned int       fe_index,
                       const unsigned int       local_index,
                       const unsigned int       global_index)
    {
      Assert (fe_index == ::DoFHandler<dim>::default_fe_index,
	      ExcMessage ("Only the default FE index is allowed for non-hp DoFHandler objects"));
      Assert (local_index<dof_handler.get_fe().dofs_per_hex,
	      ExcIndexRange (local_index, 0, dof_handler.get_fe().dofs_per_hex));
      Assert (hex_index * dof_handler.get_fe().dofs_per_hex+local_index
	      <
	      hex_dofs.size(),
	      ExcInternalError());
      
      hex_dofs[hex_index * dof_handler.get_fe().dofs_per_hex+local_index] = global_index;
    }
    
    

// explicit instantiations
    template
    unsigned int
    DoFLevel<1>::
    get_line_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
                        const unsigned int       line_index,
                        const unsigned int       fe_index,
                        const unsigned int       local_index) const;
    
    template
    void
    DoFLevel<1>::
    set_line_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
                        const unsigned int       line_index,
                        const unsigned int       fe_index,
                        const unsigned int       local_index,
                        const unsigned int       global_index);

#if deal_II_dimension >= 2

    template
    unsigned int
    DoFLevel<2>::
    get_quad_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
                        const unsigned int       quad_index,
                        const unsigned int       fe_index,
                        const unsigned int       local_index) const;
    
    template
    void
    DoFLevel<2>::
    set_quad_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
                        const unsigned int       quad_index,
                        const unsigned int       fe_index,
                        const unsigned int       local_index,
                        const unsigned int       global_index);

#endif

#if deal_II_dimension >= 3
    
    template
    unsigned int
    DoFLevel<3>::
    get_hex_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
                       const unsigned int       hex_index,
                       const unsigned int       fe_index,
                       const unsigned int       local_index) const;
    
    template
    void
    DoFLevel<3>::
    set_hex_dof_index (const ::DoFHandler<deal_II_dimension> &dof_handler,
                       const unsigned int       hex_index,
                       const unsigned int       fe_index,
                       const unsigned int       local_index,
                       const unsigned int       global_index);

#endif
    
  }
  
}
