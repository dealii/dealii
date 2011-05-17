//---------------------------------------------------------------------------
//    $Id: fe_update_flags.h,v 1.31 2005/10/24 04:33:03 guido Exp $
//    Version: $Name:  $
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__assemble_flags_h
#define __deal2__assemble_flags_h

#include <deal.II/base/config.h>
/**
 * The enum type given to the constructors of LocalAssembleBase objects, 
 * telling those objects which data to assemble on each mesh cell. 
 * When the GlobalAssembler calls the local one, it checks for each flag, 
 * and if it finds one, it assemble the corresponding object.
 *
 * By default, all flags are off, i.e. no procedure will be called.
 *
 * You can select more than one flag by concatenation
 * using the bitwise or operator|(AssembleFlags,AssembleFlags).
 */
enum AssembleFlags
  {
    //! No update
    assemble_default                      = 0,
    //! Cell term. 
    /**
     * Assemble the cell term. This is
     * usually needed, unless the matrix
     * is only a flux matrix.
     */
    assemble_cell                       = 0x0001,
    //! Assemble boundary term.
    /**
     * Calls the assemble_boundary_term for
     * each boundary face.	
     */
    assemble_boundary                  = 0x0002,
    //! Assemble face term.
    /** Call the assemble_face_term for each face of each cell in the
     * triangulation
     */
    assemble_face                      = 0x0004,
    /** Assemble rhs cell term. Used in assemble_rhs method.*/
    assemble_rhs_cell		       = 0x0008,
    /** Assemble rhs boundary terms. */
    assemble_rhs_boundary	       = 0x0010
  };





/**
 * Global operator which returns an object in which all bits are set
 * which are either set in the first or the second argument. This
 * operator exists since if it did not then the result of the bit-or
 * <tt>operator |</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * AssembleFlags.
 */
inline
AssembleFlags
operator | (AssembleFlags f1, AssembleFlags f2)
{
  return static_cast<AssembleFlags> (
				     static_cast<unsigned int> (f1) |
				     static_cast<unsigned int> (f2));
}




/**
 * Global operator which sets the bits from the second argument also
 * in the first one.
 */
inline
AssembleFlags &
operator |= (AssembleFlags &f1, AssembleFlags f2)
{
  f1 = f1 | f2;
  return f1;
}


/**
 * Global operator which returns an object in which all bits are set
 * which are set in the first as well as the second argument. This
 * operator exists since if it did not then the result of the bit-and
 * <tt>operator &</tt> would be an integer which would in turn trigger
 * a compiler warning when we tried to assign it to an object of type
 * AssembleFlags.
 */
inline
AssembleFlags
operator & (AssembleFlags f1, AssembleFlags f2)
{
  return static_cast<AssembleFlags> (
				     static_cast<unsigned int> (f1) &
				     static_cast<unsigned int> (f2));
}


/**
 * Global operator which clears all the bits in the first argument if
 * they are not also set in the second argument.
 */
inline
AssembleFlags &
operator &= (AssembleFlags &f1, AssembleFlags f2)
{
  f1 = f1 & f2;
  return f1;
}
#endif
