// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_grid_refinement_choice_h
#define dealii_grid_refinement_choice_h


#include <deal.II/base/config.h>

#include <deal.II/base/exception_macros.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/numbers.h>

#include <cstdint>
#include <string>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
template <int dim>
class RefinementChoice;

namespace internal
{
  /**
   * A helper function to create a RefinementChoice object from an integer.
   *
   * RefinementChoice objects are "singletons" (actually, "multitons" -- there
   * are multiple, but they are only a handful and these are all that can be
   * used). What is then necessary is to have a way to create these with their
   * internal id to distinguish the few possible ones in existence. We could do
   * this via a public constructor of RefinementChoice, but that would allow
   * users to create ones outside the range we envision, and we don't want to do
   * that. Rather, the constructor that takes an integer is made `private` but
   * we have this one function in an internal namespace that is a friend of the
   * class and can be used to create the objects.
   */
  template <int dim>
  constexpr RefinementChoice<dim>
  make_refinement_choice_from_int(const std::uint8_t kind);
} // namespace internal



/**
 * @brief A class which enumerates refinement possibilities for all supported
 * ReferenceCell objects.
 *
 * Not every RefinementChoice can be used with every reference cell. For
 * example, RefinementChoices::CutX<dim> only makes sense for quadrilaterals and
 * hexahedra. See ReferenceCell::is_compatible() for more information.
 */
template <int dim>
class RefinementChoice
{
public:
  /**
   * Default constructor. Creates a RefinementChoices::None<dim>.
   */
  constexpr RefinementChoice();

  /**
   * Conversion operator to the underlying integral type.
   *
   * @note Since RefinementChoices::None<dim> is always 0, the subsequent
   * conversion from an integer to a boolean may be used to check if any
   * refinement should be done (e.g., in an if-statement).
   *
   * @note This function is intended for converting to a boolean, serialization,
   * and switch statements. In most applications (such as indices into arrays of
   * prolongation or restriction matrices) the more useful integer value is the
   * one provided by ReferenceCell::refinement_choice_index().
   */
  constexpr
  operator std::uint8_t() const;

  /**
   * Comparison operator.
   */
  constexpr bool
  operator==(const RefinementChoice<dim> &other_choice) const;

  /**
   * Comparison operator.
   */
  constexpr bool
  operator!=(const RefinementChoice<dim> &other_choice) const;

  /**
   * Return a string representation of the present object, e.g., "CutX".
   */
  std::string
  to_string() const;

  /**
   * Create a new RefinementChoice which is the union of the current object
   * and @p other_choice.
   *
   * @note This function is intended to be used with anisotropic refinement. As
   * all tetrahedral refinement possibilities are mutually exclusive, this
   * operation is only defined for tetrahedra when `other_choice == *this`.
   */
  constexpr RefinementChoice<dim>
  operator|(const RefinementChoice<dim> &other_choice) const;

  /**
   * Create a new RefinementChoice which is the intersection of the current
   * object and @p other_choice.
   *
   * @note This function is intended to be used with anisotropic refinement. As
   * all tetrahedral refinement possibilities are mutually exclusive, for those
   * choices if `other_choice == *this` then this returns `*this` and otherwise
   * it returns RefinementChoices::None<3>.
   */
  constexpr RefinementChoice<dim>
  operator&(const RefinementChoice<dim> &other_choice) const;

  /**
   * Invert the RefinementChoice.
   *
   * @note This function is intended to be used with anisotropic refinement.
   * Hence, even though RefinementChoices::None<3> is valid for both hexahedra
   * and tetrahedra, this function assumes the caller is refining a hexahedron
   * and will (for that value) return RefinementChoices::Isotropic<3> (which is
   * not valid for tetrahedra).
   */
  constexpr RefinementChoice<dim>
  operator~() const;

  /**
   * Return whether or not the refinement is anisotropic.
   *
   * @note Since tetrahedra have multiple isotropic refinements, generic code
   * should use this function instead of comparing with
   * RefinementChoices::Isotropic<dim> (which is not a valid refinement choice
   * for tetrahedra).
   */
  constexpr bool
  is_anisotropic() const;

  /**
   * Return the RefinementChoice corresponding to cutting a hypercube along the
   * axis @p i. For example, if `i == 0` and `dim == 3` then the returned value
   * is RefinementChoices::CutX<3>.
   */
  static constexpr RefinementChoice<dim>
  cut_axis(const unsigned int i);

  /**
   * Return whether or not the present object is one of the RefinementChoices
   * specific to tetrahedra.
   */
  constexpr bool
  is_tetrahedral() const;

  /**
   * Return the internal line for tetrahedral refinement. This function may only
   * be called when `*this` is one of RefinementChoices::Cut49<3>,
   * RefinementChoices::Cut57<3>, or RefinementChoices::Cut68<3>.
   */
  constexpr unsigned int
  new_tetrahedron_line() const;

private:
  /**
   * The value of RefinementChoice::kind corresponding to isotropic refinement.
   */
  static constexpr std::uint8_t isotropic_kind = dim == 0 ? 1 : (1 << dim) - 1;

  /**
   * The variable that stores what this object actually corresponds to.
   */
  std::uint8_t kind;

  /**
   * Constructor.
   */
  constexpr explicit RefinementChoice(const std::uint8_t kind);

  /**
   * A kind of constructor -- not quite private because it can be called by
   * anyone, but at least hidden in an internal namespace.
   */
  template <int dim_>
  friend constexpr RefinementChoice<dim_>
  internal::make_refinement_choice_from_int(const std::uint8_t);
};


/*----------------------------- inline functions -----------------------------*/


template <int dim>
inline constexpr RefinementChoice<dim>::RefinementChoice()
  : kind(0)
{}



template <int dim>
inline constexpr RefinementChoice<dim>::RefinementChoice(
  const std::uint8_t kind)
  : kind(kind)
{
  static_assert(0 <= dim && dim <= 3, "Only implemented for 0d - 3d.");
  // Prevent users from creating illegal combinations of dim and kind. These
  // checks duplicate the values of 'kind' encoded in RefinementChoices.
  if constexpr (dim == 0)
    Assert(kind <= isotropic_kind,
           ExcMessage("The only valid RefinementChoices in 0d are None and "
                      "Isotropic."));
  else if constexpr (dim == 1)
    Assert(kind <= isotropic_kind,
           ExcMessage("The only valid RefinementChoices in 1d are None and "
                      "Isotropic."));
  else if constexpr (dim == 2)
    Assert(kind <= isotropic_kind,
           ExcMessage("The only valid RefinementChoices in 2d are None, "
                      "CutX, CutY, and Isotropic."));
  else if constexpr (dim == 3)
    Assert(kind <= isotropic_kind || (8 <= kind && kind <= 10),
           ExcMessage("This RefinementChoice isn't implemented in 3d."));
}

namespace internal
{
  template <int dim>
  inline constexpr RefinementChoice<dim>
  make_refinement_choice_from_int(const std::uint8_t kind)
  {
    return RefinementChoice<dim>(kind);
  }
} // namespace internal

/**
 * A namespace in which we define objects that correspond to specific refinement
 * choices. The objects defined here are a complete enumeration of all possible
 * refinement choices that can be used in deal.II.
 */
namespace RefinementChoices
{
  /**
   * Do not refine the cell (the default).
   */
  template <int dim>
  constexpr RefinementChoice<dim> None =
    internal::make_refinement_choice_from_int<dim>(0);

  /**
   * Refine a quadrilateral or hexahedron by cutting in the $x$ direction.
   *
   * @note This may only be used in 2d or 3d. To refine a 1d cell use
   * RefinementChoices::Isotropic.
   */
  template <int dim>
  constexpr RefinementChoice<dim> CutX =
    internal::make_refinement_choice_from_int<dim>(1);

  /**
   * Refine a quadrilateral or hexahedron by cutting in the $y$ direction.
   *
   * @note This may only be used in 2d or 3d.
   */
  template <int dim>
  constexpr RefinementChoice<dim> CutY =
    internal::make_refinement_choice_from_int<dim>(2);

  /**
   * Refine a hexahedron by cutting in the $x$ and $y$ directions.
   *
   * @note This may only be used in 3d. To refine a 2d cell use
   * RefinementChoices::Isotropic.
   */
  template <int dim>
  constexpr RefinementChoice<dim> CutXY =
    internal::make_refinement_choice_from_int<dim>(3);

  /**
   * Refine a hexahedron by cutting in the $z$ direction.
   */
  template <int dim>
  constexpr RefinementChoice<dim> CutZ =
    internal::make_refinement_choice_from_int<dim>(4);

  /**
   * Refine a hexahedron by cutting in the $x$ and $z$ directions.
   */
  template <int dim>
  constexpr RefinementChoice<dim> CutXZ =
    internal::make_refinement_choice_from_int<dim>(5);

  /**
   * Refine a hexahedron by cutting in the $y$ and $z$ directions.
   */
  template <int dim>
  constexpr RefinementChoice<dim> CutYZ =
    internal::make_refinement_choice_from_int<dim>(6);

  /**
   * Refine a cell isotropically.
   *
   * @note There are three distinct ways to isotropically refine a tetrahedron
   * and, as such, this RefinementChoice cannot be used with tetrahedra. Use one
   * of RefinementChoices::Cut49<3>, RefinementChoices::Cut57<3>, or
   * RefinementChoices::Cut68<3> instead.
   */
  template <int dim>
  constexpr RefinementChoice<dim> Isotropic =
    internal::make_refinement_choice_from_int<dim>(dim == 0 ? 1 :
                                                              (1 << dim) - 1);

  /**
   * Refine a Tetrahedron isotropically, with the internal line placed between
   * the fourth and ninth vertices.
   */
  template <int dim>
  constexpr RefinementChoice<dim> Cut49 =
    internal::make_refinement_choice_from_int<dim>(8);

  /**
   * Refine a Tetrahedron isotropically, with the internal line placed between
   * the fifth and seventh vertices.
   */
  template <int dim>
  constexpr RefinementChoice<dim> Cut57 =
    internal::make_refinement_choice_from_int<dim>(9);

  /**
   * Refine a Tetrahedron isotropically, with the internal line placed between
   * the sixth and eighth vertices.
   */
  template <int dim>
  constexpr RefinementChoice<dim> Cut68 =
    internal::make_refinement_choice_from_int<dim>(10);
} // namespace RefinementChoices



template <int dim>
inline constexpr RefinementChoice<dim>::
operator std::uint8_t() const
{
  return this->kind;
}



template <int dim>
inline constexpr bool
RefinementChoice<dim>::operator==(const RefinementChoice<dim> &other) const
{
  return this->kind == other.kind;
}



template <int dim>
inline constexpr bool
RefinementChoice<dim>::operator!=(const RefinementChoice<dim> &other) const
{
  return !(*this == other);
}



template <int dim>
inline constexpr RefinementChoice<dim>
RefinementChoice<dim>::operator|(
  const RefinementChoice<dim> &other_choice) const
{
  // Canonicalize by ensuring this->kind <= other_choice.kind
  if (*this == other_choice)
    return *this;
  if (this->kind > other_choice.kind)
    return other_choice.operator|(*this);

  if (*this == RefinementChoices::None<dim>)
    return other_choice;

  // Hypercube refinements are stored as bit fields:
  if (dim < 3 || !other_choice.is_tetrahedral())
    return RefinementChoice<dim>(this->kind | other_choice.kind);

  // Tetrahedral refinements can only be combined with None (handled above) or
  // themselves (also handled above).
  Assert(this->is_tetrahedral(),
         ExcMessage("Tetrahedral refinements (e.g., Cut49) cannot be combined "
                    "with refinements of other cells (e.g., CutX or "
                    "Isotropic)."));
  Assert(*this == other_choice,
         ExcMessage("The union of distinct tetrahedral refinements, e.g., "
                    "Cut49 with Cut57, is not defined since these are mutually "
                    "exclusive operations."));
  return RefinementChoices::None<dim>;
}



template <int dim>
inline constexpr RefinementChoice<dim>
RefinementChoice<dim>::operator&(
  const RefinementChoice<dim> &other_choice) const
{
  // Like operator|(), canonicalize by ensuring this->kind <= other_choice.kind
  if (*this == other_choice)
    return *this;
  if (this->kind > other_choice.kind)
    return other_choice.operator&(*this);

  if (*this == RefinementChoices::None<dim>)
    return *this;

  // Hypercube refinements are stored as bit fields:
  if (dim < 3 || !other_choice.is_tetrahedral())
    return RefinementChoice<dim>(this->kind & other_choice.kind);

  // Tetrahedral refinements can only be combined with None (handled above) or
  // other tetrahedral refinements. Since these are mutually exclusive their
  // intersection is empty.
  Assert(this->is_tetrahedral(),
         ExcMessage("Tetrahedral refinements (e.g., Cut49) cannot be combined "
                    "with refinements of other cells (e.g., CutX or "
                    "Isotropic)."));
  Assert(*this != other_choice, ExcInternalError());
  return RefinementChoices::None<dim>;
}



template <int dim>
inline constexpr RefinementChoice<dim>
RefinementChoice<dim>::operator~() const
{
  if (this->is_tetrahedral())
    return RefinementChoices::None<dim>;
  else
    return RefinementChoice<dim>(~(this->kind) & isotropic_kind);
}



template <int dim>
inline constexpr bool
RefinementChoice<dim>::is_anisotropic() const
{
  return !(*this == RefinementChoices::None<dim> ||
           *this == RefinementChoices::Isotropic<dim> ||
           this->is_tetrahedral());
}



template <int dim>
inline constexpr RefinementChoice<dim>
RefinementChoice<dim>::cut_axis(const unsigned int i)
{
  // Since there are no axes in 0d this also prevents calling this function
  // in that case.
  AssertIndexRange(i, dim);
  // As with RefinementCase, cutting along an axis sets the bit corresponding to
  // that axis.
  return RefinementChoice<dim>(static_cast<std::uint8_t>(1u << i));
}



template <int dim>
inline constexpr bool
RefinementChoice<dim>::is_tetrahedral() const
{
  return dim == 3 && this->kind > isotropic_kind;
}



template <int dim>
inline constexpr unsigned int
RefinementChoice<dim>::new_tetrahedron_line() const
{
  Assert(is_tetrahedral(),
         ExcMessage("This function may only be called on tetrahedral "
                    "refinement choices."));
  switch (this->kind)
    {
      case RefinementChoices::Cut49<3>:
        return 2;
      case RefinementChoices::Cut57<3>:
        return 1;
      case RefinementChoices::Cut68<3>:
        return 0;
      default:
        DEAL_II_NOT_IMPLEMENTED();
        return numbers::invalid_unsigned_int;
    }
}


DEAL_II_NAMESPACE_CLOSE

#endif
