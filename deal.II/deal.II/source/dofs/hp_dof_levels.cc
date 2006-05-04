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
#include <dofs/hp_dof_handler.h>
#include <fe/fe_collection.h>

//TODO[WB]: make a few of these functions inline


namespace internal
{
  namespace hp
  {
    template <int dim>
    unsigned int
    DoFLevel<1>::
    get_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
		   const unsigned int           line_index,
		   const unsigned int           fe_index,
		   const unsigned int           local_index,
		   internal::StructuralDimension<1> /*dummy*/) const
    {
      Assert (fe_index != ::hp::DoFHandler<dim>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));
      Assert (local_index < dof_handler.get_fe()[fe_index].dofs_per_line,
              ExcIndexRange(local_index, 0,
                            dof_handler.get_fe()[fe_index].dofs_per_line));
      Assert (line_index < line_dof_offsets.size(),
              ExcIndexRange (line_index, 0, line_dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (line_dof_offsets[line_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
      if (dim == 1)
        {
          Assert (fe_index == this->active_fe_indices[line_index],
                  ExcMessage ("FE index does not match that of the present cell"));
          return line_dofs[line_dof_offsets[line_index]+local_index];
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object. hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = line_dof_offsets[line_index];
          const unsigned int *pointer        = &line_dofs[starting_offset];
          while (true)
            {
              Assert (*pointer != deal_II_numbers::invalid_unsigned_int,
                      ExcInternalError());
              if (*pointer == fe_index)
                return *(pointer + 1 + local_index);
              else
                pointer += dof_handler.get_fe()[*pointer].dofs_per_line + 1;
            }
        }
    }



    template <int dim>
    void
    DoFLevel<1>::
    set_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
                   const unsigned int           line_index,
                   const unsigned int           fe_index,
                   const unsigned int           local_index,
                   const unsigned int           global_index,
                   internal::StructuralDimension<1> /*dummy*/)
    {
      Assert (fe_index != ::hp::DoFHandler<dim>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));
      Assert (local_index < dof_handler.get_fe()[fe_index].dofs_per_line,
              ExcIndexRange(local_index, 0,
                            dof_handler.get_fe()[fe_index].dofs_per_line));
      Assert (line_index < line_dof_offsets.size(),
              ExcIndexRange (line_index, 0, line_dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (line_dof_offsets[line_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
      if (dim == 1)
        {
          Assert (fe_index == this->active_fe_indices[line_index],
                  ExcMessage ("FE index does not match that of the present cell"));
          line_dofs[line_dof_offsets[line_index]+local_index] = global_index;
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object.  hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = line_dof_offsets[line_index];
          unsigned int      *pointer         = &line_dofs[starting_offset];
          while (true)
            {
              Assert (*pointer != deal_II_numbers::invalid_unsigned_int,
                      ExcInternalError());
              if (*pointer == fe_index)
                {
                  *(pointer + 1 + local_index) = global_index;
                  return;
                }
              else
                pointer += dof_handler.get_fe()[*pointer].dofs_per_line + 1;
            }
        }  
    }



    template <int dim>
    unsigned int
    DoFLevel<1>::
    n_active_fe_indices (const ::hp::DoFHandler<dim> &dof_handler,
                         const unsigned int           line_index,
                         internal::StructuralDimension<1> /*dummy*/) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (line_index < line_dof_offsets.size(),
              ExcIndexRange (line_index, 0, line_dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (line_dof_offsets[line_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
      if (dim == 1)
        {
          return 1;
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object. hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = line_dof_offsets[line_index];
          const unsigned int *pointer        = &line_dofs[starting_offset];
          unsigned int counter = 0;
          while (true)
            {
              if (*pointer == deal_II_numbers::invalid_unsigned_int)
                                                 // end of list reached
                return counter;
              else
                {
                  ++counter;
                  pointer += dof_handler.get_fe()[*pointer].dofs_per_line + 1;
                }
            }
        }
    }



    template <int dim>
    bool
    DoFLevel<1>::
    fe_index_is_active (const ::hp::DoFHandler<dim> &dof_handler,
                        const unsigned int           line_index,
                        const unsigned int           fe_index,
                        internal::StructuralDimension<1> /*dummy*/) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (line_index < line_dof_offsets.size(),
              ExcIndexRange (line_index, 0, line_dof_offsets.size()));
      Assert (fe_index != ::hp::DoFHandler<dim>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (line_dof_offsets[line_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
      if (dim == 1)
        {
          Assert (line_index < active_fe_indices.size(),
                  ExcInternalError());
          return (fe_index == active_fe_indices[line_index]);
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object. hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = line_dof_offsets[line_index];
          const unsigned int *pointer        = &line_dofs[starting_offset];
          while (true)
            {
              if (*pointer == deal_II_numbers::invalid_unsigned_int)
                                                 // end of list reached
                return false;
              else
                if (*pointer == fe_index)
                  return true;
                else
                  pointer += dof_handler.get_fe()[*pointer].dofs_per_line + 1;
            }
        }
    }
    


    template <int dim>
    unsigned int
    DoFLevel<2>::
    get_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
                   const unsigned int           quad_index,
                   const unsigned int           fe_index,
                   const unsigned int           local_index,
                   internal::StructuralDimension<2> /*dummy*/) const
    {
      Assert (dim >= 2, ExcMessage ("You can only access quads in 2d or higher"));
      Assert (fe_index != ::hp::DoFHandler<dim>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));
      Assert (local_index < dof_handler.get_fe()[fe_index].dofs_per_quad,
              ExcIndexRange(local_index, 0,
                            dof_handler.get_fe()[fe_index].dofs_per_quad));
      Assert (quad_index < quad_dof_offsets.size(),
              ExcIndexRange (quad_index, 0, quad_dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (quad_dof_offsets[quad_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
      if (dim == 2)
        {
          Assert (fe_index == this->active_fe_indices[quad_index],
                  ExcMessage ("FE index does not match that of the present cell"));
          return quad_dofs[quad_dof_offsets[quad_index]+local_index];
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object. hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = quad_dof_offsets[quad_index];
          const unsigned int *pointer        = &quad_dofs[starting_offset];
          while (true)
            {
              Assert (*pointer != deal_II_numbers::invalid_unsigned_int,
                      ExcInternalError());
              if (*pointer == fe_index)
                return *(pointer + 1 + local_index);
              else
                pointer += dof_handler.get_fe()[*pointer].dofs_per_quad + 1;
            }
        }
    }



    template <int dim>
    void
    DoFLevel<2>::
    set_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
                   const unsigned int           quad_index,
                   const unsigned int           fe_index,
                   const unsigned int           local_index,
                   const unsigned int           global_index,
                   internal::StructuralDimension<2> /*dummy*/)
    {
      Assert (dim >= 2, ExcMessage ("You can only access quads in 2d or higher"));
      Assert (fe_index != ::hp::DoFHandler<dim>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));
      Assert (local_index < dof_handler.get_fe()[fe_index].dofs_per_quad,
              ExcIndexRange(local_index, 0,
                            dof_handler.get_fe()[fe_index].dofs_per_quad));
      Assert (quad_index < quad_dof_offsets.size(),
              ExcIndexRange (quad_index, 0, quad_dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (quad_dof_offsets[quad_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
      if (dim == 2)
        {
          Assert (fe_index == this->active_fe_indices[quad_index],
                  ExcMessage ("FE index does not match that of the present cell"));
          quad_dofs[quad_dof_offsets[quad_index]+local_index] = global_index;
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object.  hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = quad_dof_offsets[quad_index];
          unsigned int      *pointer         = &quad_dofs[starting_offset];
          while (true)
            {
              Assert (*pointer != deal_II_numbers::invalid_unsigned_int,
                      ExcInternalError());
              if (*pointer == fe_index)
                {
                  *(pointer + 1 + local_index) = global_index;
                  return;
                }
              else
                pointer += dof_handler.get_fe()[*pointer].dofs_per_quad + 1;
            }
        }  
    }



    template <int dim>
    unsigned int
    DoFLevel<2>::
    n_active_fe_indices (const ::hp::DoFHandler<dim> &dof_handler,
                         const unsigned int           quad_index,
                         internal::StructuralDimension<2> /*dummy*/) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (quad_index < quad_dof_offsets.size(),
              ExcIndexRange (quad_index, 0, quad_dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (quad_dof_offsets[quad_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 2d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
      if (dim == 2)
        {
          return 1;
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object. hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = quad_dof_offsets[quad_index];
          const unsigned int *pointer        = &quad_dofs[starting_offset];
          unsigned int counter = 0;
          while (true)
            {
              if (*pointer == deal_II_numbers::invalid_unsigned_int)
                                                 // end of list reached
                return counter;
              else
                {
                  ++counter;
                  pointer += dof_handler.get_fe()[*pointer].dofs_per_quad + 1;
                }
            }
        }
    }



    template <int dim>
    bool
    DoFLevel<2>::
    fe_index_is_active (const ::hp::DoFHandler<dim> &dof_handler,
                        const unsigned int           quad_index,
                        const unsigned int           fe_index,
                        internal::StructuralDimension<2> /*dummy*/) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (quad_index < quad_dof_offsets.size(),
              ExcIndexRange (quad_index, 0, quad_dof_offsets.size()));
      Assert (fe_index != ::hp::DoFHandler<dim>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (quad_dof_offsets[quad_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 2d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
      if (dim == 2)
        {
          Assert (quad_index < active_fe_indices.size(),
                  ExcInternalError());
          return (fe_index == active_fe_indices[quad_index]);
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object. hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = quad_dof_offsets[quad_index];
          const unsigned int *pointer        = &quad_dofs[starting_offset];
          while (true)
            {
              if (*pointer == deal_II_numbers::invalid_unsigned_int)
                                                 // end of list reached
                return false;
              else
                if (*pointer == fe_index)
                  return true;
                else
                  pointer += dof_handler.get_fe()[*pointer].dofs_per_quad + 1;
            }
        }
    }



    template <int dim>
    unsigned int
    DoFLevel<3>::
    get_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
                   const unsigned int           hex_index,
                   const unsigned int           fe_index,
                   const unsigned int           local_index,
                   internal::StructuralDimension<3> /*dummy*/) const
    {
      Assert (dim >= 3, ExcMessage ("You can only access hexs in 3d or higher"));
      Assert (fe_index != ::hp::DoFHandler<dim>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));
      Assert (local_index < dof_handler.get_fe()[fe_index].dofs_per_hex,
              ExcIndexRange(local_index, 0,
                            dof_handler.get_fe()[fe_index].dofs_per_hex));
      Assert (hex_index < hex_dof_offsets.size(),
              ExcIndexRange (hex_index, 0, hex_dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (hex_dof_offsets[hex_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
      if (dim == 3)
        {
          Assert (fe_index == this->active_fe_indices[hex_index],
                  ExcMessage ("FE index does not match that of the present cell"));
          return hex_dofs[hex_dof_offsets[hex_index]+local_index];
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object. hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = hex_dof_offsets[hex_index];
          const unsigned int *pointer        = &hex_dofs[starting_offset];
          while (true)
            {
              Assert (*pointer != deal_II_numbers::invalid_unsigned_int,
                      ExcInternalError());
              if (*pointer == fe_index)
                return *(pointer + 1 + local_index);
              else
                pointer += dof_handler.get_fe()[*pointer].dofs_per_hex + 1;
            }
        }
    }



    template <int dim>
    void
    DoFLevel<3>::
    set_dof_index (const ::hp::DoFHandler<dim> &dof_handler,
                   const unsigned int           hex_index,
                   const unsigned int           fe_index,
                   const unsigned int           local_index,
                   const unsigned int           global_index,
                   internal::StructuralDimension<3> /*dummy*/)
    {
      Assert (dim >= 3, ExcMessage ("You can only access hexs in 3d or higher"));
      Assert (fe_index != ::hp::DoFHandler<dim>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));
      Assert (local_index < dof_handler.get_fe()[fe_index].dofs_per_hex,
              ExcIndexRange(local_index, 0,
                            dof_handler.get_fe()[fe_index].dofs_per_hex));
      Assert (hex_index < hex_dof_offsets.size(),
              ExcIndexRange (hex_index, 0, hex_dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (hex_dof_offsets[hex_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
      if (dim == 3)
        {
          Assert (fe_index == this->active_fe_indices[hex_index],
                  ExcMessage ("FE index does not match that of the present cell"));
          hex_dofs[hex_dof_offsets[hex_index]+local_index] = global_index;
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object.  hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = hex_dof_offsets[hex_index];
          unsigned int      *pointer         = &hex_dofs[starting_offset];
          while (true)
            {
              Assert (*pointer != deal_II_numbers::invalid_unsigned_int,
                      ExcInternalError());
              if (*pointer == fe_index)
                {
                  *(pointer + 1 + local_index) = global_index;
                  return;
                }
              else
                pointer += dof_handler.get_fe()[*pointer].dofs_per_hex + 1;
            }
        }  
    }



    template <int dim>
    unsigned int
    DoFLevel<3>::
    n_active_fe_indices (const ::hp::DoFHandler<dim> &dof_handler,
                         const unsigned int           hex_index,
                         internal::StructuralDimension<3> /*dummy*/) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (hex_index < hex_dof_offsets.size(),
              ExcIndexRange (hex_index, 0, hex_dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (hex_dof_offsets[hex_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 3d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
      if (dim == 3)
        {
          return 1;
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object. hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = hex_dof_offsets[hex_index];
          const unsigned int *pointer        = &hex_dofs[starting_offset];
          unsigned int counter = 0;
          while (true)
            {
              if (*pointer == deal_II_numbers::invalid_unsigned_int)
                                                 // end of list reached
                return counter;
              else
                {
                  ++counter;
                  pointer += dof_handler.get_fe()[*pointer].dofs_per_hex + 1;
                }
            }
        }
    }
    


    template <int dim>
    bool
    DoFLevel<3>::
    fe_index_is_active (const ::hp::DoFHandler<dim> &dof_handler,
                        const unsigned int           hex_index,
                        const unsigned int           fe_index,
                        internal::StructuralDimension<3> /*dummy*/) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (hex_index < hex_dof_offsets.size(),
              ExcIndexRange (hex_index, 0, hex_dof_offsets.size()));
      Assert (fe_index != ::hp::DoFHandler<dim>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (hex_dof_offsets[hex_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 2d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
      if (dim == 3)
        {
          Assert (hex_index < active_fe_indices.size(),
                  ExcInternalError());
          return (fe_index == active_fe_indices[hex_index]);
        }
      else
        {
                                           // we are in higher space
                                           // dimensions, so there may
                                           // be multiple finite
                                           // elements associated with
                                           // this object. hop along
                                           // the list of index sets
                                           // until we find the one
                                           // with the correct
                                           // fe_index, and then poke
                                           // into that part. trigger
                                           // an exception if we can't
                                           // find a set for this
                                           // particular fe_index
          const unsigned int starting_offset = hex_dof_offsets[hex_index];
          const unsigned int *pointer        = &hex_dofs[starting_offset];
          while (true)
            {
              if (*pointer == deal_II_numbers::invalid_unsigned_int)
                                                 // end of list reached
                return false;
              else
                if (*pointer == fe_index)
                  return true;
                else
                  pointer += dof_handler.get_fe()[*pointer].dofs_per_hex + 1;
            }
        }
    }

    

// explicit instantiations
    template
    unsigned int
    DoFLevel<1>::
    get_dof_index (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                   const unsigned int           line_index,
                   const unsigned int           fe_index,
                   const unsigned int           local_index,
                   internal::StructuralDimension<1> dummy) const;
    
    template
    void
    DoFLevel<1>::
    set_dof_index (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                   const unsigned int           line_index,
                   const unsigned int           fe_index,
                   const unsigned int           local_index,
                   const unsigned int           global_index,
                   internal::StructuralDimension<1> dummy);

    template
    unsigned int
    DoFLevel<1>::
    n_active_fe_indices (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                         const unsigned int           line_index,
                         internal::StructuralDimension<1> dummy) const;

    template
    bool
    DoFLevel<1>::
    fe_index_is_active (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                        const unsigned int           line_index,
                        const unsigned int           fe_index,
                        internal::StructuralDimension<1> dummy) const;
    
#if deal_II_dimension >= 2

    template
    unsigned int
    DoFLevel<2>::
    get_dof_index (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                   const unsigned int           quad_index,
                   const unsigned int           fe_index,
                   const unsigned int           local_index,
                   internal::StructuralDimension<2> dummy) const;
    
    template
    void
    DoFLevel<2>::
    set_dof_index (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                   const unsigned int           quad_index,
                   const unsigned int           fe_index,
                   const unsigned int           local_index,
                   const unsigned int           global_index,
                   internal::StructuralDimension<2> dummy);

    template
    unsigned int
    DoFLevel<2>::
    n_active_fe_indices (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                         const unsigned int           quad_index,
                         internal::StructuralDimension<2> dummy) const;
    

    template
    bool
    DoFLevel<2>::
    fe_index_is_active (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                        const unsigned int           line_index,
                        const unsigned int           fe_index,
                        internal::StructuralDimension<2> dummy) const;
#endif

#if deal_II_dimension >= 3
    
    template
    unsigned int
    DoFLevel<3>::
    get_dof_index (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                   const unsigned int           hex_index,
                   const unsigned int           fe_index,
                   const unsigned int           local_index,
                   internal::StructuralDimension<3> dummy) const;
    
    template
    void
    DoFLevel<3>::
    set_dof_index (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                   const unsigned int           hex_index,
                   const unsigned int           fe_index,
                   const unsigned int           local_index,
                   const unsigned int           global_index,
                   internal::StructuralDimension<3> dummy);

    template
    unsigned int
    DoFLevel<3>::
    n_active_fe_indices (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                         const unsigned int           hex_index,
                         internal::StructuralDimension<3> dummy) const;

    template
    bool
    DoFLevel<3>::
    fe_index_is_active (const ::hp::DoFHandler<deal_II_dimension> &dof_handler,
                        const unsigned int           line_index,
                        const unsigned int           fe_index,
                        internal::StructuralDimension<3> dummy) const;
#endif
  }
}
