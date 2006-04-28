//----------------------------  hp_dof_levels.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  hp_dof_levels.cc  ------------------------


#include <base/memory_consumption.h>
#include <dofs/hp_dof_levels.h>


namespace internal
{
  namespace hp
  {
    unsigned int
    DoFLevel<1>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (line_dofs) +
              MemoryConsumption::memory_consumption (dof_line_index_offset));
    }



    unsigned int
    DoFLevel<2>::memory_consumption () const
    {
      return (DoFLevel<1>::memory_consumption () +
              MemoryConsumption::memory_consumption (quad_dofs) +
              MemoryConsumption::memory_consumption (dof_quad_index_offset));
    }



    unsigned int
    DoFLevel<3>::memory_consumption () const
    {
      return (DoFLevel<2>::memory_consumption () +
              MemoryConsumption::memory_consumption (hex_dofs) +
              MemoryConsumption::memory_consumption (dof_hex_index_offset));
    }

    template <int dim>
    void
    DoFLevel<0>::
    set_hp_vertex_dof_index (const ::hp::FECollection<dim> &fe,
                             const unsigned int           fe_index,
                             const unsigned int          *start_of_list,
                             const unsigned int           local_index,
                             const unsigned int           global_index)
    {
      Assert (fe_index != hp::DoFHandler<dim>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working with hp DoFHandlers"));
      abort ();
    }


    template <int dim>
    unsigned int
    DoFLevel<0>::
    get_hp_vertex_dof_index (const ::hp::FECollection<dim> &fe,
                             const unsigned int           fe_index,
                             const unsigned int          *start_of_list,
                             const unsigned int           local_index)
    {
      Assert (fe_index != hp::DoFHandler<dim>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working with hp DoFHandlers"));
      Assert (local_index < fe[fe_index].dofs_per_vertex,
              ExcIndexRange(local_index, 0, fe[fe_index].dofs_per_vertex));
                                       // hop along the list of index
                                       // sets until we find the one
                                       // with the correct fe_index, and
                                       // then poke into that
                                       // part. trigger an exception if
                                       // we can't find a set for this
                                       // particular fe_index
      const unsigned int *pointer = start_of_list;
      while (true)
        {
          Assert (*pointer != deal_II_numbers::invalid_unsigned_int,
                  ExcInternalError());
          if (*pointer == fe_index)
            return *(pointer + 1 + local_index);
          else
            pointer += fe[*pointer].dofs_per_vertex;
        }
    }  


// explicit instantiations
    template
    void
    DoFLevel<0>::
    set_hp_vertex_dof_index (const ::hp::FECollection<deal_II_dimension> &fe,
                             const unsigned int           fe_index,
                             const unsigned int          *start_of_list,
                             const unsigned int           local_index,
                             const unsigned int           global_index);
    
    template
    unsigned int
    DoFLevel<0>::
    get_hp_vertex_dof_index (const ::hp::FECollection<deal_II_dimension> &fe,
                             const unsigned int           fe_index,
                             const unsigned int          *start_of_list,
                             const unsigned int           local_index);
  }
}
