// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__fe_block_mask_h
#define __deal2__fe_block_mask_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>

#include <vector>
#include <iosfwd>

DEAL_II_NAMESPACE_OPEN



/**
 * This class represents a mask that can be used to select individual
 * vector blocks of a finite element (see also
 * @ref GlossBlockMask "this glossary entry"). It will typically have as many
 * elements as the finite element has blocks, and one can
 * use <code>operator[]</code> to query whether a particular block
 * has been selected.
 *
 * The semantics of this class are the same as the related ComponentMask
 * class, i.e., a default constructed mask represents all possible
 * blocks. See there for more information about these semantics.
 *
 * Objects of this kind are used in many places where one wants to restrict
 * operations to a certain subset of blocks, e.g. in DoFTools::extract_dofs.
 * These objects can
 * either be created by hand, or, simpler, by asking the finite element
 * to generate a block mask from certain selected blocks using
 * code such as this where we create a mask that only denotes the
 * velocity block of a Stokes element (see @ref vector_valued):
 * @code
 *   FESystem<dim> stokes_fe (FESystem<dim>(FE_Q<dim>(2), dim), 1,    // Q2 element for the velocities
 *                            FE_Q<dim>(1),                     1);     // Q1 element for the pressure
 *   FEValuesExtractors::Scalar pressure(dim);
 *   BlockMask pressure_mask = stokes_fe.block_mask (pressure);
 * @endcode
 * Note that by wrapping the velocity elements into a single FESystem
 * object we make sure that the overall element has only 2 blocks.
 * The result is a block mask that, in both 2d and 3d, would have values
 * <code>[false, true]</code>. (Compare this to the corresponding
 * component mask discussed in the ComponentMask documentation.)
 * Similarly, using
 * @code
 *   FEValuesExtractors::Vector velocities(0);
 *   BlockMask velocity_mask = stokes_fe.block_mask (velocities);
 * @endcode
 * would result in a mask <code>[true, false]</code> in both 2d and 3d.
 *
 * @ingroup fe
 * @author Wolfgang Bangerth
 * @date 2012
 * @ingroup vector_valued
 */
class BlockMask
{
public:
  /**
   * Initialize a block mask. The default is that a block
   * mask represents a set of blocks that are <i>all</i>
   * selected, i.e., calling this constructor results in
   * a block mask that always returns <code>true</code>
   * whenever asked whether a block is selected.
   */
  BlockMask ();

  /**
   * Initialize an object of this type with a set of selected
   * blocks specified by the argument.
   *
   * @param block_mask A vector of <code>true/false</code>
   * entries that determine which blocks of a finite element
   * are selected. If the length of the given vector is zero,
   * then this interpreted as the case where <i>every</i> block
   * is selected.
   */
  BlockMask (const std::vector<bool> &block_mask);

  /**
   * Initialize the block mask with a number of elements that
   * are either all true or false.
   *
   * @param n_blocks The number of elements of this mask
   * @param initializer The value each of these elements is supposed to have:
   *                    either true or false.
   */
  BlockMask (const unsigned int n_blocks,
             const bool         initializer);

  /**
   * If this block mask has been initialized with a mask of
   * size greater than zero, then return the size of the mask
   * represented by this object. On the other hand, if this
   * mask has been initialized as an empty object that represents
   * a mask that is true for every element (i.e., if this object
   * would return true when calling represents_the_all_selected_mask())
   * then return zero since no definite size is known.
   */
  unsigned int size () const;

  /**
   * Return whether a particular block is selected by this
   * mask. If this mask represents the case of an object that
   * selects <i>all blocks</i> (e.g. if it is created
   * using the default constructor or is converted from an
   * empty vector of type bool) then this function returns
   * true regardless of the given argument.
   *
   * @param block_index The index for which the function
   * should return whether the block is selected. If this
   * object represents a mask in which all blocks are always
   * selected then any index is allowed here. Otherwise, the
   * given index needs to be between zero and the number of
   * blocks that this mask represents.
   */
  bool operator[] (const unsigned int block_index) const;

  /**
   * Return whether this block mask represents a mask with
   * exactly <code>n</code> blocks. This is true if either
   * it was initilized with a vector with exactly <code>n</code>
   * entries of type <code>bool</code> (in this case, @p n must
   * equal the result of size()) or if it was initialized
   * with an empty vector (or using the default constructor) in
   * which case it can represent a mask with an arbitrary number
   * of blocks and will always say that a block is
   * selected.
   */
  bool
  represents_n_blocks (const unsigned int n) const;

  /**
   * Return the number of blocks that are selected by this mask.
   *
   * Since empty block masks represent a block mask that
   * would return <code>true</code> for every block, this
   * function may not know the true size of the block
   * mask and it therefore requires an argument that denotes the
   * overall number of blocks.
   *
   * If the object has been initialized with a non-empty mask (i.e.,
   * if the size() function returns something greater than zero,
   * or equivalently if represents_the_all_selected_mask() returns
   * false) then the argument can be omitted and the result of size()
   * is taken.
   */
  unsigned int
  n_selected_blocks (const unsigned int overall_number_of_blocks = numbers::invalid_unsigned_int) const;

  /**
   * Return the index of the first selected block. The argument
   * is there for the same reason it exists with the
   * n_selected_blocks() function.
   *
   * The function throws an exception if no block is selected at all.
   */
  unsigned int
  first_selected_block (const unsigned int overall_number_of_blocks = numbers::invalid_unsigned_int) const;

  /**
   * Return true if this mask represents a default
   * constructed mask that corresponds to one in which
   * all blocks are selected. If true, then the size()
   * function will return zero.
   */
  bool
  represents_the_all_selected_mask () const;

  /**
   * Return a block mask that contains the union of the
   * blocks selected by the current object and the one
   * passed as an argument.
   */
  BlockMask operator | (const BlockMask &mask) const;

  /**
   * Return a block mask that has only those elements set that
   * are set both in the current object as well as the one
   * passed as an argument.
   */
  BlockMask operator & (const BlockMask &mask) const;

  /**
   * Return whether this object and the argument are identical.
   */
  bool operator== (const BlockMask &mask) const;

  /**
   * Return whether this object and the argument are not identical.
   */
  bool operator!= (const BlockMask &mask) const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t
  memory_consumption () const;

private:
  /**
   * The actual block mask.
   */
  std::vector<bool> block_mask;

  // make the output operator a friend so it can access
  // the block_mask array
  friend
  std::ostream &operator << (std::ostream &out,
                             const BlockMask &mask);
};


/**
 * Write a block mask to an output stream. If the block
 * mask represents one where all blocks are selected without
 * specifying a particular size of the mask, then it
 * writes the string <code>[all blocks selected]</code> to the
 * stream. Otherwise, it prints the block mask in a form like
 * <code>[true,true,true,false]</code>.
 *
 * @param out The stream to write to.
 * @param mask The mask to write.
 * @return A reference to the first argument.
 */
std::ostream &operator << (std::ostream &out,
                           const BlockMask &mask);


// -------------------- inline functions ---------------------

inline
BlockMask::BlockMask()
{}


inline
BlockMask::BlockMask(const std::vector<bool> &block_mask)
  :
  block_mask (block_mask)
{}


inline
BlockMask::BlockMask(const unsigned int n_blocks,
                     const bool         initializer)
  :
  block_mask (n_blocks, initializer)
{}


inline
unsigned int
BlockMask::size () const
{
  return block_mask.size();
}


inline
bool
BlockMask::operator [](const unsigned int block_index) const
{
  // if the mask represents the all-block mask
  // then always return true
  if (block_mask.size() == 0)
    return true;
  else
    {
      // otherwise check the validity of the index and
      // return whatever is appropriate
      Assert (block_index < block_mask.size(),
              ExcIndexRange (block_index, 0, block_mask.size()));
      return block_mask[block_index];
    }
}


inline
bool
BlockMask::represents_n_blocks(const unsigned int n) const
{
  return ((block_mask.size() == 0)
          ||
          (block_mask.size() == n));
}


inline
unsigned int
BlockMask::n_selected_blocks(const unsigned int n) const
{
  if ((n != numbers::invalid_unsigned_int) && (size() > 0))
    AssertDimension (n, size());

  const unsigned int real_n = (n != numbers::invalid_unsigned_int
                               ?
                               n
                               :
                               size());
  if (block_mask.size() == 0)
    return real_n;
  else
    {
      AssertDimension (real_n, block_mask.size());
      unsigned int c = 0;
      for (unsigned int i=0; i<block_mask.size(); ++i)
        if (block_mask[i] == true)
          ++c;
      return c;
    }
}


inline
unsigned int
BlockMask::first_selected_block(const unsigned int n) const
{
  if ((n != numbers::invalid_unsigned_int) && (size() > 0))
    AssertDimension (n, size());

  if (block_mask.size() == 0)
    return 0;
  else
    {
      for (unsigned int c=0; c<block_mask.size(); ++c)
        if (block_mask[c] == true)
          return c;

      Assert (false, ExcMessage ("No block is selected at all!"));
      return numbers::invalid_unsigned_int;
    }
}



inline
bool
BlockMask::represents_the_all_selected_mask () const
{
  return (block_mask.size() == 0);
}



inline
BlockMask
BlockMask::operator | (const BlockMask &mask) const
{
  // if one of the two masks denotes the all-block mask,
  // then return the other one
  if (block_mask.size() == 0)
    return mask;
  else if (mask.block_mask.size() == 0)
    return *this;
  else
    {
      // if both masks have individual entries set, form
      // the combination of the two
      AssertDimension(block_mask.size(), mask.block_mask.size());
      std::vector<bool> new_mask (block_mask.size());
      for (unsigned int i=0; i<block_mask.size(); ++i)
        new_mask[i] = (block_mask[i] || mask.block_mask[i]);

      return new_mask;
    }
}


inline
BlockMask
BlockMask::operator & (const BlockMask &mask) const
{
  // if one of the two masks denotes the all-block mask,
  // then return the other one
  if (block_mask.size() == 0)
    return mask;
  else if (mask.block_mask.size() == 0)
    return *this;
  else
    {
      // if both masks have individual entries set, form
      // the combination of the two
      AssertDimension(block_mask.size(), mask.block_mask.size());
      std::vector<bool> new_mask (block_mask.size());
      for (unsigned int i=0; i<block_mask.size(); ++i)
        new_mask[i] = (block_mask[i] && mask.block_mask[i]);

      return new_mask;
    }
}


inline
bool
BlockMask::operator== (const BlockMask &mask) const
{
  return block_mask == mask.block_mask;
}


inline
bool
BlockMask::operator!= (const BlockMask &mask) const
{
  return block_mask != mask.block_mask;
}


DEAL_II_NAMESPACE_CLOSE

#endif
