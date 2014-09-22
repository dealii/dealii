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

#ifndef __deal2__fe_component_mask_h
#define __deal2__fe_component_mask_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_consumption.h>

#include <vector>
#include <iosfwd>

DEAL_II_NAMESPACE_OPEN



/**
 * This class represents a mask that can be used to select individual
 * vector components of a finite element (see also
 * @ref GlossComponentMask "this glossary entry"). It will typically have as many
 * elements as the finite element has vector components, and one can
 * use <code>operator[]</code> to query whether a particular component
 * has been selected.
 *
 * Objects of this kind are used in many places where one wants to restrict
 * operations to a certain subset of components, e.g. in DoFTools::make_zero_boundary_values
 * or VectorTools::interpolate_boundary_values. These objects can
 * either be created by hand, or, simpler, by asking the finite element
 * to generate a component mask from certain selected components using
 * code such as this where we create a mask that only denotes the
 * velocity components of a Stokes element (see @ref vector_valued):
 * @code
 *   FESystem<dim> stokes_fe (FE_Q<dim>(2), dim,    // Q2 element for the velocities
 *                            FE_Q<dim>(1), 1);     // Q1 element for the pressure
 *   FEValuesExtractors::Scalar pressure(dim);
 *   ComponentMask pressure_mask = stokes_fe.component_mask (pressure);
 * @endcode
 * The result is a component mask that, in 2d, would have values
 * <code>[false, false, true]</code>. Similarly, using
 * @code
 *   FEValuesExtractors::Vector velocities(0);
 *   ComponentMask velocity_mask = stokes_fe.component_mask (velocities);
 * @endcode
 * would result in a mask <code>[true, true, false]</code> in 2d. Of
 * course, in 3d, the result would be <code>[true, true, true, false]</code>.
 *
 * @ingroup fe
 * @author Wolfgang Bangerth
 * @date 2012
 * @ingroup vector_valued
 */
class ComponentMask
{
public:
  /**
   * Initialize a component mask. The default is that a component
   * mask represents a set of components that are <i>all</i>
   * selected, i.e., calling this constructor results in
   * a component mask that always returns <code>true</code>
   * whenever asked whether a component is selected.
   */
  ComponentMask ();

  /**
   * Initialize an object of this type with a set of selected
   * components specified by the argument.
   *
   * @param component_mask A vector of <code>true/false</code>
   * entries that determine which components of a finite element
   * are selected. If the length of the given vector is zero,
   * then this interpreted as the case where <i>every</i> component
   * is selected.
   */
  ComponentMask (const std::vector<bool> &component_mask);

  /**
   * Initialize the component mask with a number of elements that
   * are either all true or false.
   *
   * @param n_components The number of elements of this mask
   * @param initializer The value each of these elements is supposed to have:
   *                    either true or false.
   */
  ComponentMask (const unsigned int n_components,
                 const bool         initializer);

  /**
   * Set a particular entry in the mask to a value.
   */
  void set (const unsigned int index, const bool value);

  /**
   * If this component mask has been initialized with a mask of
   * size greater than zero, then return the size of the mask
   * represented by this object. On the other hand, if this
   * mask has been initialized as an empty object that represents
   * a mask that is true for every element (i.e., if this object
   * would return true when calling represents_the_all_selected_mask())
   * then return zero since no definite size is known.
   */
  unsigned int size () const;

  /**
   * Return whether a particular component is selected by this
   * mask. If this mask represents the case of an object that
   * selects <i>all components</i> (e.g. if it is created
   * using the default constructor or is converted from an
   * empty vector of type bool) then this function returns
   * true regardless of the given argument.
   *
   * @param component_index The index for which the function
   * should return whether the component is selected. If this
   * object represents a mask in which all components are always
   * selected then any index is allowed here. Otherwise, the
   * given index needs to be between zero and the number of
   * components that this mask represents.
   */
  bool operator[] (const unsigned int component_index) const;

  /**
   * Return whether this component mask represents a mask with
   * exactly <code>n</code> components. This is true if either
   * it was initilized with a vector with exactly <code>n</code>
   * entries of type <code>bool</code> (in this case, @p n must
   * equal the result of size()) or if it was initialized
   * with an empty vector (or using the default constructor) in
   * which case it can represent a mask with an arbitrary number
   * of components and will always say that a component is
   * selected.
   */
  bool
  represents_n_components (const unsigned int n) const;

  /**
   * Return the number of components that are selected by this mask.
   *
   * Since empty component masks represent a component mask that
   * would return <code>true</code> for every component, this
   * function may not know the true size of the component
   * mask and it therefore requires an argument that denotes the
   * overall number of components.
   *
   * If the object has been initialized with a non-empty mask (i.e.,
   * if the size() function returns something greater than zero,
   * or equivalently if represents_the_all_selected_mask() returns
   * false) then the argument can be omitted and the result of size()
   * is taken.
   */
  unsigned int
  n_selected_components (const unsigned int overall_number_of_components = numbers::invalid_unsigned_int) const;

  /**
   * Return the index of the first selected component. The argument
   * is there for the same reason it exists with the
   * n_selected_components() function.
   *
   * The function throws an exception if no component is selected at all.
   */
  unsigned int
  first_selected_component (const unsigned int overall_number_of_components = numbers::invalid_unsigned_int) const;

  /**
   * Return true if this mask represents a default
   * constructed mask that corresponds to one in which
   * all components are selected. If true, then the size()
   * function will return zero.
   */
  bool
  represents_the_all_selected_mask () const;

  /**
   * Return a component mask that contains the union of the
   * components selected by the current object and the one
   * passed as an argument.
   */
  ComponentMask operator | (const ComponentMask &mask) const;

  /**
   * Return a component mask that has only those elements set that
   * are set both in the current object as well as the one
   * passed as an argument.
   */
  ComponentMask operator & (const ComponentMask &mask) const;

  /**
   * Return whether this object and the argument are identical.
   */
  bool operator== (const ComponentMask &mask) const;

  /**
   * Return whether this object and the argument are not identical.
   */
  bool operator!= (const ComponentMask &mask) const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t
  memory_consumption () const;

  /**
   * Exception
   */
  DeclException0 (ExcNoComponentSelected);

private:
  /**
   * The actual component mask.
   */
  std::vector<bool> component_mask;

  // make the output operator a friend so it can access
  // the component_mask array
  friend
  std::ostream &operator << (std::ostream &out,
                             const ComponentMask &mask);
};


/**
 * Write a component mask to an output stream. If the component
 * mask represents one where all components are selected without
 * specifying a particular size of the mask, then it
 * writes the string <code>[all components selected]</code> to the
 * stream. Otherwise, it prints the component mask in a form like
 * <code>[true,true,true,false]</code>.
 *
 * @param out The stream to write to.
 * @param mask The mask to write.
 * @return A reference to the first argument.
 */
std::ostream &operator << (std::ostream &out,
                           const ComponentMask &mask);


// -------------------- inline functions ---------------------

inline
ComponentMask::ComponentMask()
{}


inline
ComponentMask::ComponentMask(const std::vector<bool> &component_mask)
  :
  component_mask (component_mask)
{}


inline
ComponentMask::ComponentMask(const unsigned int n_components,
                             const bool         initializer)
  :
  component_mask (n_components, initializer)
{}


inline
unsigned int
ComponentMask::size () const
{
  return component_mask.size();
}


inline
void
ComponentMask::set(const unsigned int index, const bool value)
{
  AssertIndexRange(index, component_mask.size());
  component_mask[index] = value;
}


inline
bool
ComponentMask::operator [](const unsigned int component_index) const
{
  // if the mask represents the all-component mask
  // then always return true
  if (component_mask.size() == 0)
    return true;
  else
    {
      // otherwise check the validity of the index and
      // return whatever is appropriate
      AssertIndexRange (component_index, component_mask.size());
      return component_mask[component_index];
    }
}


inline
bool
ComponentMask::represents_n_components(const unsigned int n) const
{
  return ((component_mask.size() == 0)
          ||
          (component_mask.size() == n));
}


inline
unsigned int
ComponentMask::n_selected_components(const unsigned int n) const
{
  if ((n != numbers::invalid_unsigned_int) && (size() > 0))
    AssertDimension (n, size());

  const unsigned int real_n = (n != numbers::invalid_unsigned_int
                               ?
                               n
                               :
                               size());
  if (component_mask.size() == 0)
    return real_n;
  else
    {
      AssertDimension (real_n, component_mask.size());
      unsigned int c = 0;
      for (unsigned int i=0; i<component_mask.size(); ++i)
        if (component_mask[i] == true)
          ++c;
      return c;
    }
}


inline
unsigned int
ComponentMask::first_selected_component(const unsigned int n) const
{
  if ((n != numbers::invalid_unsigned_int) && (size() > 0))
    AssertDimension (n, size());

  if (component_mask.size() == 0)
    return 0;
  else
    {
      for (unsigned int c=0; c<component_mask.size(); ++c)
        if (component_mask[c] == true)
          return c;

      Assert (false, ExcMessage ("No component is selected at all!"));
      return numbers::invalid_unsigned_int;
    }
}



inline
bool
ComponentMask::represents_the_all_selected_mask () const
{
  return (component_mask.size() == 0);
}



inline
ComponentMask
ComponentMask::operator | (const ComponentMask &mask) const
{
  // if one of the two masks denotes the all-component mask,
  // then return the other one
  if (component_mask.size() == 0)
    return mask;
  else if (mask.component_mask.size() == 0)
    return *this;
  else
    {
      // if both masks have individual entries set, form
      // the combination of the two
      AssertDimension(component_mask.size(), mask.component_mask.size());
      std::vector<bool> new_mask (component_mask.size());
      for (unsigned int i=0; i<component_mask.size(); ++i)
        new_mask[i] = (component_mask[i] || mask.component_mask[i]);

      return new_mask;
    }
}


inline
ComponentMask
ComponentMask::operator & (const ComponentMask &mask) const
{
  // if one of the two masks denotes the all-component mask,
  // then return the other one
  if (component_mask.size() == 0)
    return mask;
  else if (mask.component_mask.size() == 0)
    return *this;
  else
    {
      // if both masks have individual entries set, form
      // the combination of the two
      AssertDimension(component_mask.size(), mask.component_mask.size());
      std::vector<bool> new_mask (component_mask.size());
      for (unsigned int i=0; i<component_mask.size(); ++i)
        new_mask[i] = (component_mask[i] && mask.component_mask[i]);

      return new_mask;
    }
}


inline
bool
ComponentMask::operator== (const ComponentMask &mask) const
{
  return component_mask == mask.component_mask;
}


inline
bool
ComponentMask::operator!= (const ComponentMask &mask) const
{
  return component_mask != mask.component_mask;
}



DEAL_II_NAMESPACE_CLOSE

#endif
