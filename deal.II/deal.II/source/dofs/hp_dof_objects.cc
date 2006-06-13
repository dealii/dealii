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


#include <base/memory_consumption.h>
#include <dofs/hp_dof_objects.h>
#include <dofs/hp_dof_levels.h>
#include <dofs/hp_dof_handler.h>


namespace internal
{
  namespace hp
  {
    template <int structdim>
    unsigned int
    DoFObjects<structdim>::memory_consumption () const
    {
      return (MemoryConsumption::memory_consumption (dofs) +
              MemoryConsumption::memory_consumption (dof_offsets));
    }

#if deal_II_dimension == 1    
    template <>
    template <>
    unsigned int
    DoFObjects<1>::
    get_dof_index<1> (const ::hp::DoFHandler<1> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           obj_level) const
    {
      Assert (fe_index != ::hp::DoFHandler<1>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
      
      Assert (fe_index == dof_handler.levels[obj_level]->active_fe_indices[obj_index],
	      ExcMessage ("FE index does not match that of the present cell"));
      return dofs[dof_offsets[obj_index]+local_index];
    }


    template <>
    template <>
    void
    DoFObjects<1>::
    set_dof_index<1> (const ::hp::DoFHandler<1> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           global_index,
		      const unsigned int           obj_level)
    {
      Assert (fe_index != ::hp::DoFHandler<1>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
          Assert (fe_index == dof_handler.levels[obj_level]->active_fe_indices[obj_index],
                  ExcMessage ("FE index does not match that of the present cell"));
          dofs[dof_offsets[obj_index]+local_index] = global_index;
    }


    


    template <>
    template <>
    unsigned int
    DoFObjects<1>::
    n_active_fe_indices<1> (const ::hp::DoFHandler<1> &dof_handler,
			    const unsigned int           obj_index) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
         return 1;
    }



    template <>
    template <>
    bool
    DoFObjects<1>::
    fe_index_is_active<1> (const ::hp::DoFHandler<1> &dof_handler,
			   const unsigned int           obj_index,
			   const unsigned int           fe_index,
			   const unsigned int           obj_level) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));
      Assert (fe_index != ::hp::DoFHandler<1>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 1d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
          Assert (obj_index < dof_handler.levels[obj_level]->active_fe_indices.size(),
                  ExcInternalError());
          return (fe_index == dof_handler.levels[obj_level]->active_fe_indices[obj_index]);
    }
#endif
#if deal_II_dimension == 2    

    template <>
    template <>
    unsigned int
    DoFObjects<1>::
    get_dof_index<2> (const ::hp::DoFHandler<2> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           ) const
    {
      Assert (fe_index != ::hp::DoFHandler<2>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
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
          const unsigned int starting_offset = dof_offsets[obj_index];
          const unsigned int *pointer        = &dofs[starting_offset];
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


    template <>
    template <>
    void
    DoFObjects<1>::
    set_dof_index<2> (const ::hp::DoFHandler<2> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           global_index,
		      const unsigned int           )
    {
      Assert (fe_index != ::hp::DoFHandler<2>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

    
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
          const unsigned int starting_offset = dof_offsets[obj_index];
          unsigned int      *pointer         = &dofs[starting_offset];
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





    template <>
    template <>
    unsigned int
    DoFObjects<1>::
    n_active_fe_indices<2> (const ::hp::DoFHandler<2> &dof_handler,
			    const unsigned int           obj_index) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
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
          const unsigned int starting_offset = dof_offsets[obj_index];
          const unsigned int *pointer        = &dofs[starting_offset];
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



    template <>
    template <>
    bool
    DoFObjects<1>::
    fe_index_is_active<2> (const ::hp::DoFHandler<2> &dof_handler,
			   const unsigned int           obj_index,
			   const unsigned int           fe_index,
			   const unsigned int           ) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));
      Assert (fe_index != ::hp::DoFHandler<2>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
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
          const unsigned int starting_offset = dof_offsets[obj_index];
          const unsigned int *pointer        = &dofs[starting_offset];
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
#endif
#if deal_II_dimension == 3    

    template <>
    template <>
    unsigned int
    DoFObjects<1>::
    get_dof_index<3> (const ::hp::DoFHandler<3> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           ) const
    {
      Assert (fe_index != ::hp::DoFHandler<3>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
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
          const unsigned int starting_offset = dof_offsets[obj_index];
          const unsigned int *pointer        = &dofs[starting_offset];
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



    template <>
    template <>
    void
    DoFObjects<1>::
    set_dof_index<3> (const ::hp::DoFHandler<3> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           global_index,
		      const unsigned int           )
    {
      Assert (fe_index != ::hp::DoFHandler<3>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

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
          const unsigned int starting_offset = dof_offsets[obj_index];
          unsigned int      *pointer         = &dofs[starting_offset];
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



    template <>
    template <>
    unsigned int
    DoFObjects<1>::
    n_active_fe_indices<3> (const ::hp::DoFHandler<3> &dof_handler,
			    const unsigned int           obj_index) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
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
          const unsigned int starting_offset = dof_offsets[obj_index];
          const unsigned int *pointer        = &dofs[starting_offset];
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



    template <>
    template <>
    bool
    DoFObjects<1>::
    fe_index_is_active<3> (const ::hp::DoFHandler<3> &dof_handler,
			   const unsigned int           obj_index,
			   const unsigned int           fe_index,
			   const unsigned int           ) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));
      Assert (fe_index != ::hp::DoFHandler<3>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
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
          const unsigned int starting_offset = dof_offsets[obj_index];
          const unsigned int *pointer        = &dofs[starting_offset];
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
    
#endif
#if deal_II_dimension == 2    

    template <>
    template <>
    unsigned int
    DoFObjects<2>::
    get_dof_index<2> (const ::hp::DoFHandler<2> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           obj_level) const
    {
      Assert (fe_index != ::hp::DoFHandler<2>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

                                       // if we are in 2d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
          Assert (fe_index == dof_handler.levels[obj_level]->active_fe_indices[obj_index],
                  ExcMessage ("FE index does not match that of the present cell"));
          return dofs[dof_offsets[obj_index]+local_index];
    }



    template <>
    template <>
    void
    DoFObjects<2>::
    set_dof_index<2> (const ::hp::DoFHandler<2> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           global_index,
		      const unsigned int           obj_level)
    {
      Assert (fe_index != ::hp::DoFHandler<2>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

                                       // if we are in 2d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
      Assert (fe_index == dof_handler.levels[obj_level]->active_fe_indices[obj_index],
                  ExcMessage ("FE index does not match that of the present cell"));
          dofs[dof_offsets[obj_index]+local_index] = global_index;
    }



    template <>
    template <>
    unsigned int
    DoFObjects<2>::
    n_active_fe_indices<2> (const ::hp::DoFHandler<2> &dof_handler,
			    const unsigned int           obj_index) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 2d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
          return 1;
    }



    template <>
    template <>
    bool
    DoFObjects<2>::
    fe_index_is_active<2> (const ::hp::DoFHandler<2> &dof_handler,
			   const unsigned int           obj_index,
			   const unsigned int           fe_index,
			   const unsigned int           obj_level) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));
      Assert (fe_index != ::hp::DoFHandler<2>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 2d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
          Assert (obj_index < dof_handler.levels[obj_level]->active_fe_indices.size(),
                  ExcInternalError());
          return (fe_index == dof_handler.levels[obj_level]->active_fe_indices[obj_index]);
    }
#endif
    
#if deal_II_dimension == 3    
    template <>
    template <>
    unsigned int
    DoFObjects<2>::
    get_dof_index<3> (const ::hp::DoFHandler<3> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           ) const
    {
      Assert (fe_index != ::hp::DoFHandler<3>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

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
          const unsigned int starting_offset = dof_offsets[obj_index];
          const unsigned int *pointer        = &dofs[starting_offset];
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



    template <>
    template <>
    void
    DoFObjects<2>::
    set_dof_index<3> (const ::hp::DoFHandler<3> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           global_index,
		      const unsigned int           )
    {
      Assert (fe_index != ::hp::DoFHandler<3>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

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
          const unsigned int starting_offset = dof_offsets[obj_index];
          unsigned int      *pointer         = &dofs[starting_offset];
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



    template <>
    template <>
    unsigned int
    DoFObjects<2>::
    n_active_fe_indices<3> (const ::hp::DoFHandler<3> &dof_handler,
			    const unsigned int           obj_index) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
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
          const unsigned int starting_offset = dof_offsets[obj_index];
          const unsigned int *pointer        = &dofs[starting_offset];
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



    template <>
    template <>
    bool
    DoFObjects<2>::
    fe_index_is_active<3> (const ::hp::DoFHandler<3> &dof_handler,
			   const unsigned int           obj_index,
			   const unsigned int           fe_index,
			   const unsigned int           ) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));
      Assert (fe_index != ::hp::DoFHandler<3>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
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
          const unsigned int starting_offset = dof_offsets[obj_index];
          const unsigned int *pointer        = &dofs[starting_offset];
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



    template <>
    template <>
    unsigned int
    DoFObjects<3>::
    get_dof_index<3> (const ::hp::DoFHandler<3> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           obj_level) const
    {
      Assert (fe_index != ::hp::DoFHandler<3>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

                                       // if we are in 3d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
           Assert (fe_index == dof_handler.levels[obj_level]->active_fe_indices[obj_index],
                  ExcMessage ("FE index does not match that of the present cell"));
          return dofs[dof_offsets[obj_index]+local_index];
    }



    template <>
    template <>
    void
    DoFObjects<3>::
    set_dof_index<3> (const ::hp::DoFHandler<3> &dof_handler,
		      const unsigned int           obj_index,
		      const unsigned int           fe_index,
		      const unsigned int           local_index,
		      const unsigned int           global_index,
		      const unsigned int           obj_level)
    {
      Assert (fe_index != ::hp::DoFHandler<3>::default_fe_index,
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
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));

                                       // if we are in 3d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique. then
                                       // fe_index must be
                                       // active_fe_index
          Assert (fe_index == dof_handler.levels[obj_level]->active_fe_indices[obj_index],
                  ExcMessage ("FE index does not match that of the present cell"));
          dofs[dof_offsets[obj_index]+local_index] = global_index;
    }



    template <>
    template <>
    unsigned int
    DoFObjects<3>::
    n_active_fe_indices<3> (const ::hp::DoFHandler<3> &dof_handler,
			    const unsigned int           obj_index) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 3d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
           return 1;
    }
    


    template <>
    template <>
    bool
    DoFObjects<3>::
    fe_index_is_active<3> (const ::hp::DoFHandler<3> &dof_handler,
			   const unsigned int           obj_index,
			   const unsigned int           fe_index,
			   const unsigned int           obj_level) const
    {
      Assert (&dof_handler != 0,
              ExcMessage ("No DoFHandler is specified for this iterator"));
      Assert (&dof_handler.get_fe() != 0,
              ExcMessage ("No finite element collection is associated with "
                          "this DoFHandler"));
      Assert (obj_index < dof_offsets.size(),
              ExcIndexRange (obj_index, 0, dof_offsets.size()));
      Assert (fe_index != ::hp::DoFHandler<3>::default_fe_index,
              ExcMessage ("You need to specify a FE index when working "
                          "with hp DoFHandlers"));
      Assert (fe_index < dof_handler.get_fe().size(),
              ExcIndexRange (fe_index, 0, dof_handler.get_fe().size()));

                                       // make sure we are on an
                                       // object for which DoFs have
                                       // been allocated at all
      Assert (dof_offsets[obj_index] != deal_II_numbers::invalid_unsigned_int,
              ExcMessage ("You are trying to access degree of freedom "
                          "information for an object on which no such "
                          "information is available"));
      
                                       // if we are in 3d, then the
                                       // only set of indices we store
                                       // is the one for the cell,
                                       // which is unique
          Assert (obj_index < dof_handler.levels[obj_level]->active_fe_indices.size(),
                  ExcInternalError());
          return (fe_index == dof_handler.levels[obj_level]->active_fe_indices[obj_index]);
    }
#endif

				     // explicit instantiations
    template
    unsigned int
    DoFObjects<1>::
    memory_consumption () const;
    
#if deal_II_dimension >= 2

    template
    unsigned int
    DoFObjects<2>::
    memory_consumption () const;
    
#endif
    
    
#if deal_II_dimension >= 3

    template
    unsigned int
    DoFObjects<3>::
    memory_consumption () const;
    
#endif
  }
}
