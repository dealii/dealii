// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2002 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_filtered_iterator_h
#define dealii_filtered_iterator_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/iterator_range.h>
#include <deal.II/base/types.h>

#include <deal.II/grid/tria_iterator_base.h>

#include <memory>
#include <set>
#include <tuple>

#ifdef DEAL_II_HAVE_CXX20
#  include <concepts>
#endif


DEAL_II_NAMESPACE_OPEN


/**
 * In this namespace a number of classes is declared that may be used as
 * filters in the FilteredIterator class. The filters either check for binary
 * information (for example, the IteratorFilters::Active filter class checks
 * whether the object pointed to is active), or for valued information by
 * comparison with prescribed values (for example, the LevelEqualTo filter
 * class checks whether the level of the object pointed to by the iterator
 * under consideration is equal to a value that was given to the filter upon
 * construction.
 *
 * For examples of use of these classes as well as requirements on filters see
 * the general description of the FilteredIterator class.
 *
 * @ingroup Iterators
 */
namespace IteratorFilters
{
  /**
   * Filter that evaluates to true if either the iterator points to an active
   * object or an iterator past the end.
   *
   * @ingroup Iterators
   */
  class Active
  {
  public:
    /**
     * Evaluate the iterator and return true if the object is active or past
     * the end.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };

  /**
   * Filter that evaluates to true if either the iterator points to an object
   * for which the user flag is set or an iterator past the end. See
   * @ref GlossUserFlags
   * for information about user flags.
   *
   * @ingroup Iterators
   */
  class UserFlagSet
  {
  public:
    /**
     * Evaluate the iterator and return true if the object has a set user flag
     * or past the end.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };


  /**
   * Filter that evaluates to true if either the iterator points to an object
   * for which the user flag is not set or an iterator past the end. Inverse
   * filter to the previous class.
   *
   * @ingroup Iterators
   */
  class UserFlagNotSet
  {
  public:
    /**
     * Evaluate the iterator and return true if the object has an unset user
     * flag or past the end.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };


  /**
   * Filter for iterators that evaluates to true if either the iterator is
   * past the end or the level of the object pointed to is equal to a value
   * given to the constructor.
   *
   * @ingroup Iterators
   */
  class LevelEqualTo
  {
  public:
    /**
     * Constructor. Store the level which iterators shall have to be evaluated
     * to true.
     */
    LevelEqualTo(const unsigned int level);

    /**
     * Evaluation operator. Returns true if either the level of the object
     * pointed to is equal to the stored value or the iterator is past the
     * end.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;

  protected:
    /**
     * Stored value to compare the level with.
     */
    const unsigned int level;
  };



  /**
   * Filter for iterators that evaluates to true if either the iterator is
   * past the end or the subdomain id of the object pointed to is equal to a
   * value given to the constructor, assuming that the iterator allows
   * querying for a subdomain id.
   *
   * @ingroup Iterators
   */
  class SubdomainEqualTo
  {
  public:
    /**
     * Constructor. Store the subdomain which iterators shall have to be
     * evaluated to true.
     */
    SubdomainEqualTo(const types::subdomain_id subdomain_id);

    /**
     * Evaluation operator. Returns true if either the subdomain of the object
     * pointed to is equal to the stored value or the iterator is past the
     * end.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;

  protected:
    /**
     * Stored value to compare the subdomain with.
     */
    const types::subdomain_id subdomain_id;
  };



  /**
   * Filter for iterators that evaluates to true if a cell is owned by the
   * current processor, i.e., if it is a
   * @ref GlossLocallyOwnedCell "locally owned cell".
   *
   * This class is used in step-32, in connection with the methods of the
   * @ref distributed
   * topic.
   *
   * @ingroup Iterators
   */
  class LocallyOwnedCell
  {
  public:
    /**
     * Evaluation operator. Returns true if the cell is locally owned.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };



  /**
   * Filter for iterators that evaluates to true if the level subdomain id of
   * a cell is equal to the current processor id.
   *
   * @ingroup Iterators
   */
  class LocallyOwnedLevelCell
  {
  public:
    /**
     * Evaluation operator. Returns true if the level subdomain id of the cell
     * is equal to the current processor id.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };


  /**
   * Filter for iterators that evaluates to true if the iterator of the object
   * pointed to is equal to a value or set of values given to the constructor,
   * assuming that the iterator allows querying for a material id.
   *
   *
   * @ingroup Iterators
   */
  class MaterialIdEqualTo
  {
  public:
    /**
     * Constructor. Store the material id which iterators shall have to be
     * evaluated to true and state if the iterator must be locally owned.
     */
    MaterialIdEqualTo(const types::material_id material_id,
                      const bool               only_locally_owned = false);

    /**
     * Constructor. Store a collection of material ids which iterators shall
     * have to be evaluated to true and state if the iterator must be locally
     * owned.
     */
    MaterialIdEqualTo(const std::set<types::material_id> &material_ids,
                      const bool only_locally_owned = false);

    /**
     * Evaluation operator. Returns true if the material id of the object
     * pointed to is equal within the stored set of value allowable values
     * and, if required, if the cell is locally owned.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;

  protected:
    /**
     * Stored value to compare the material id with.
     */
    const std::set<types::material_id> material_ids;
    /**
     * Flag stating whether only locally owned cells must return true.
     */
    const bool only_locally_owned;
  };

  /**
   * Filter for iterators that evaluates to true if the iterator of the object
   * pointed to is equal to a value or set of values given to the constructor,
   * assuming that the iterator allows querying for an active FE index.
   *
   *
   * @ingroup Iterators
   */
  class ActiveFEIndexEqualTo
  {
  public:
    /**
     * Constructor. Store the active FE index which iterators shall have to be
     * evaluated to true and state if the iterator must be locally owned.
     */
    ActiveFEIndexEqualTo(const unsigned int active_fe_index,
                         const bool         only_locally_owned = false);

    /**
     * Constructor. Store a collection of active FE indices which iterators
     * shall have to be evaluated to true and state if the iterator must be
     * locally owned.
     */
    ActiveFEIndexEqualTo(const std::set<unsigned int> &active_fe_indices,
                         const bool only_locally_owned = false);

    /**
     * Evaluation operator. Returns true if the active FE index of the object
     * pointed to is equal within the stored set of value allowable values
     * and, if required, if the cell is locally owned.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;

  protected:
    /**
     * Stored value to compare the material id with.
     */
    const std::set<unsigned int> active_fe_indices;
    /**
     * Flag stating whether only locally owned cells must return true.
     */
    const bool only_locally_owned;
  };

  /**
   * Filter for iterators that evaluates to true if the iterator of the object
   * pointed to is on the boundary.
   *
   *
   * @ingroup Iterators
   */
  class AtBoundary
  {
  public:
    /**
     * Evaluate the iterator and return true if the object at the boundary.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;
  };


  /**
   * Filter for iterators that evaluates to true if the iterator of the object
   * pointed to is equal to a value or set of values given to the constructor,
   * assuming that the iterator allows querying for a boundary id.
   *
   * @ingroup Iterators
   */
  class BoundaryIdEqualTo
  {
  public:
    /**
     * Constructor. Store the boundary id which iterators shall have to be
     * evaluated to true.
     */
    BoundaryIdEqualTo(const types::boundary_id boundary_id);

    /**
     * Constructor. Store a collection of boundary ids which iterators shall
     * have to be evaluated to true.
     */
    BoundaryIdEqualTo(const std::set<types::boundary_id> &boundary_ids);

    /**
     * Evaluation operator. Returns true if the boundary id of the object
     * pointed to is equal within the stored set of value allowable values.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;

  protected:
    /**
     * Stored value to compare the material id with.
     */
    const std::set<types::boundary_id> boundary_ids;
  };


  /**
   * Filter for iterators that evaluates to true if the iterator of the object
   * pointed to is equal to a value or set of values given to the constructor,
   * assuming that the iterator allows querying for a manifold id.
   *
   * @ingroup Iterators
   */
  class ManifoldIdEqualTo
  {
  public:
    /**
     * Constructor. Store the boundary id which iterators shall have to be
     * evaluated to true.
     */
    ManifoldIdEqualTo(const types::manifold_id manifold_id);

    /**
     * Constructor. Store a collection of boundary ids which iterators shall
     * have to be evaluated to true.
     */
    ManifoldIdEqualTo(const std::set<types::manifold_id> &manifold_ids);

    /**
     * Evaluation operator. Returns true if the boundary id of the object
     * pointed to is equal within the stored set of value allowable values.
     */
    template <class Iterator>
    bool
    operator()(const Iterator &i) const;

  protected:
    /**
     * Stored value to compare the material id with.
     */
    const std::set<types::manifold_id> manifold_ids;
  };
} // namespace IteratorFilters


/**
 * This class provides a certain view on a range of triangulation or
 * DoFHandler iterators by only iterating over elements that satisfy a given
 * filter (called a <em>predicate</em>, following the notation of the C++
 * standard library). Once initialized with a predicate and a value for the
 * iterator, a filtered iterator advances to the next or previous element that
 * satisfies the predicate if operators ++ or \-- are invoked. Intermediate
 * iterator values that lie in between but do not satisfy the predicate are
 * skipped. It is thus very simple to write loops over a certain class of
 * objects without the need to explicitly write down the condition they have
 * to satisfy in each loop iteration. This in particular is helpful if
 * functions are called with a pair of iterators denoting a range on which
 * they shall act, by choosing a filtered iterator instead of usual ones.
 * The use of this class is particularly useful when writing "ranged-based"
 * `for` loops of the form
 * @code
 *   for (const auto &cell : <some FilteredIterator object>)
 *     { ... }
 * @endcode
 * An example is to write
 * @code
 *   DoFHandler<dim> dof_handler;
 *   ...
 *   for (const auto &cell :
 *          dof_handler.active_cell_iterators()
 *          | IteratorFilters::LocallyOwnedCell()
 *          | IteratorFilters::AtBoundary())
 *     {
 *       fe_values.reinit (cell);
 *       ...do the local integration on 'cell'...;
 *     }
 * @endcode
 * where the resulting loop iterates over all active cells that are also
 * locally owned and that point to cells that are at the boundary of the
 * domain.
 *
 * This class implements functionality comparable to what was added to C++20
 * in the form of
 * [range adaptors](https://en.cppreference.com/w/cpp/ranges/filter_view)
 * and the filters and filter views that represent them.
 *
 * This class is used in step-32.
 *
 *
 * <h3>Predicates</h3>
 *
 * The object that represent the condition an iterator has to satisfy only
 * has to provide an interface that allows to call the evaluation operator,
 * i.e. <code>bool operator() (const BaseIterator&)</code>. This includes
 * function pointers as well as classes that implement an <code>bool operator
 * ()(const BaseIterator&)</code>. Then, the FilteredIterator will skip all
 * objects where the return value of this function is <code>false</code>.
 *
 *
 * An example of a simple valid predicate is the following: given the function
 * @code
 *   template <typename BIterator>
 *   bool level_equal_to_3 (const BIterator& c)
 *   {
 *     return (static_cast<unsigned int>(c->level()) == 3);
 *   };
 * @endcode
 * then
 * @code
 *   &level_equal_to_3<typename Triangulation<dim>::active_cell_iterator>
 * @endcode
 * is a valid predicate.
 *
 * Likewise, given the following binary function
 * @code
 *   template <typename BIterator>
 *   bool level_equal_to (const BIterator&     c,
 *                        const unsigned int level)
 *   {
 *     return (static_cast<unsigned int>(c->level()) == level);
 *   };
 * @endcode
 * then the lambda function
 * @code
 *   [](const BIterator& c) { return level_equal_to<active_cell_iterator>(c,
 * 3);}
 * @endcode
 * is another valid predicate (here: a function that returns true if either
 * the iterator is past the end or the level is equal to the second argument;
 * this second argument is taken considered fixed when creating the lambda
 * function).
 *
 * Finally, classes can be predicates. The following class is one:
 * @code
 *   class Active
 *   {
 *   public:
 *     template <class Iterator>
 *     bool operator () (const Iterator &i) const
 *     {
 *       return i->is_active();
 *     }
 *   };
 * @endcode
 * and objects of this type can be used as predicates. Likewise, this more
 * complicated one can also be used:
 * @code
 *   class SubdomainEqualTo
 *   {
 *   public:
 *     SubdomainEqualTo (const types::subdomain_id subdomain_id)
 *       : subdomain_id (subdomain_id)
 *     {};
 *
 *     template <class Iterator>
 *     bool operator () (const Iterator &i) const
 *     {
 *       return (i->subdomain_id() == subdomain_id);
 *     }
 *
 *   private:
 *     const types::subdomain_id subdomain_id;
 *   };
 * @endcode
 * Objects like <code>SubdomainEqualTo(3)</code> can then be used as
 * predicates.
 *
 * Since whenever a predicate is evaluated it is checked that the iterator
 * checked is actually valid (i.e. not past the end), no checks for this case
 * have to be performed inside predicates.
 *
 * A number of filter classes are already implemented in the IteratorFilters
 * namespace, but writing different ones is simple following the examples
 * above.
 *
 *
 * <h3>Initialization of filtered iterators</h3>
 *
 * Filtered iterators are given a predicate at construction time which cannot
 * be changed any more. This behavior would be expected if the predicate
 * would have been given as a template parameter to the class, but since that
 * would make the declaration of filtered iterators a nightmare, we rather
 * give the predicate as an unchangeable entity to the constructor. Note that
 * one can assign a filtered iterator with one predicate to another filtered
 * iterator with another type; yet, this does <em>not</em> change the
 * predicate of the assigned-to iterator, only the pointer indicating the
 * iterator is changed.
 *
 * If a filtered iterator is not assigned a value of the underlying
 * (unfiltered) iterator type, the default value is taken. If, however, a
 * value is given to the constructor, that value has either to be past the
 * end, or has to satisfy the predicate. For example, if the predicate only
 * evaluates to true if the level of an object is equal to three, then
 * <code>tria.begin_active(3)</code> would be a valid choice while
 * <code>tria.begin()</code> would not since the latter also returns iterators
 * to non-active cells which always start at level 0.
 *
 * Since one often only has some iterator and wants to set a filtered iterator
 * to the first one that satisfies a predicate (for example, the first one for
 * which the user flag is set, or the first one with a given subdomain id),
 * there are assignment functions #set_to_next_positive and
 * #set_to_previous_positive that assign the next or last previous iterator
 * that satisfies the predicate, i.e. they follow the list of iterators in
 * either direction until they find a matching one (or the past-the-end
 * iterator). Like the <code>operator=</code> they return the resulting value
 * of the filtered iterator.
 *
 *
 * <h3>Examples</h3>
 *
 * The following call counts the number of active cells that have a set user
 * flag:
 * @code
 *   FilteredIterator<typename Triangulation<dim>::active_cell_iterator>
 *     begin (IteratorFilters::UserFlagSet()),
 *     end (IteratorFilters::UserFlagSet());
 *   begin.set_to_next_positive(tria.begin_active());
 *   end = tria.end();
 *   n_flagged_cells = std::distance (begin, end);
 * @endcode
 * Note that by the @p set_to_next_positive call the first cell with a set
 * user flag was assigned to the @p begin iterator. For the end iterator, no
 * such call was necessary, since the past-the-end iterator always satisfies
 * all predicates.
 *
 * The same can be achieved by the following snippet, though harder to read:
 * @code
 *   using FI =
 *     FilteredIterator<typename Triangulation<dim>::active_cell_iterator>;
 *   n_flagged_cells =
 *     std::distance (
 *       FI(IteratorFilters::UserFlagSet()).set_to_next_positive(
 *         tria.begin_active()),
 *       FI(IteratorFilters::UserFlagSet(), tria.end()));
 * @endcode
 * It relies on the fact that if we create an unnamed filtered iterator with a
 * given predicate but no iterator value and assign it the next positive value
 * with respect to this predicate, it returns itself which is then used as the
 * first parameter to the @p std::distance function. This procedure is not
 * necessary for the end element to this function here, since the past-the-end
 * iterator always satisfies the predicate so that we can assign this value to
 * the filtered iterator directly in the constructor.
 *
 * Finally, the following loop only assembles the matrix on cells with
 * subdomain id equal to three:
 * @code
 * FilteredIterator<typename Triangulation<dim>::active_cell_iterator>
 *   cell (IteratorFilters::SubdomainEqualTo(3)),
 *   endc (IteratorFilters::SubdomainEqualTo(3), tria.end());
 * cell.set_to_next_positive (tria.begin_active());
 * for (; cell!=endc; ++cell)
 *   assemble_local_matrix (cell);
 * @endcode
 *
 * Since comparison between filtered and unfiltered iterators is defined, we
 * could as well have let the @p endc variable in the last example be of type
 * Triangulation::active_cell_iterator since it is unchanged and its value
 * does not depend on the filter.
 *
 * @ingroup grid
 * @ingroup Iterators
 */
template <typename BaseIterator>
class FilteredIterator : public BaseIterator
{
public:
  /**
   * Typedef to the accessor type of the underlying iterator.
   */
  using AccessorType = typename BaseIterator::AccessorType;

  /**
   * Constructor. Set the iterator to the default state and use the given
   * predicate for filtering subsequent assignment and iteration.
   *
   * @dealiiConceptRequires{(std::predicate<Predicate, BaseIterator>)}
   */
  template <typename Predicate>
  DEAL_II_CXX20_REQUIRES((std::predicate<Predicate, BaseIterator>))
  FilteredIterator(Predicate p);

  /**
   * Constructor. Use the given predicate for filtering and initialize the
   * iterator with the given value.
   *
   * If the initial value @p bi does not satisfy the predicate @p p then it is
   * advanced until we either hit the past-the-end iterator, or the
   * predicate is satisfied. This allows, for example, to write code like
   * @code
   * FilteredIterator<typename Triangulation<dim>::active_cell_iterator>
   *   cell (IteratorFilters::SubdomainEqualTo(13),
   *         triangulation.begin_active());
   * @endcode
   *
   * If the cell <code>triangulation.begin_active()</code> does not have a
   * subdomain_id equal to 13, then the iterator will automatically be
   * advanced to the first cell that has.
   *
   * @dealiiConceptRequires{(std::predicate<Predicate, BaseIterator>)}
   */
  template <typename Predicate>
  DEAL_II_CXX20_REQUIRES((std::predicate<Predicate, BaseIterator>))
  FilteredIterator(Predicate p, const BaseIterator &bi);

  /**
   * Copy constructor. Copy the predicate and iterator value of the given
   * argument.
   */
  FilteredIterator(const FilteredIterator &fi);

  /**
   * Assignment operator. Copy the iterator value of the argument, but as
   * discussed in the class documentation, the predicate of the argument is
   * not copied. The iterator value underlying the argument has to satisfy the
   * predicate of the object assigned to, as given at its construction time.
   */
  FilteredIterator &
  operator=(const FilteredIterator &fi);

  /**
   * Assignment operator. Copy the iterator value of the argument, and keep
   * the predicate of this object. The given iterator value has to satisfy the
   * predicate of the object assigned to, as given at its construction time.
   */
  FilteredIterator &
  operator=(const BaseIterator &fi);

  /**
   * Search for the next iterator from @p bi onwards that satisfies the
   * predicate of this object and assign it to this object.
   *
   * Since filtered iterators are automatically converted to the underlying
   * base iterator type, you can also give a filtered iterator as argument to
   * this function.
   */
  FilteredIterator &
  set_to_next_positive(const BaseIterator &bi);

  /**
   * As above, but search for the previous iterator from @p bi backwards that
   * satisfies the predicate of this object and assign it to this object.
   *
   * Since filtered iterators are automatically converted to the underlying
   * base iterator type, you can also give a filtered iterator as argument to
   * this function.
   */
  FilteredIterator &
  set_to_previous_positive(const BaseIterator &bi);

  /**
   * Compare for equality of the underlying iterator values of this and the
   * given object.
   *
   * We do not compare for equality of the predicates.
   */
  bool
  operator==(const FilteredIterator &fi) const;

  /**
   * Compare for equality of the underlying iterator value of this object with
   * the given object.
   *
   * The predicate of this object is irrelevant for this operation.
   */
  bool
  operator==(const BaseIterator &fi) const;

  /**
   * Compare for inequality of the underlying iterator values of this and the
   * given object.
   *
   * We do not compare for equality of the predicates.
   */
  bool
  operator!=(const FilteredIterator &fi) const;

  /**
   * Compare for inequality of the underlying iterator value of this object
   * with the given object.
   *
   * The predicate of this object is irrelevant for this operation.
   */
  bool
  operator!=(const BaseIterator &fi) const;

  /**
   * Compare for ordering of the underlying iterator values of this and the
   * given object.
   *
   * We do not compare the predicates.
   */
  bool
  operator<(const FilteredIterator &fi) const;

  /**
   * Compare for ordering of the underlying iterator value of this object with
   * the given object.
   *
   * The predicate of this object is irrelevant for this operation.
   */
  bool
  operator<(const BaseIterator &fi) const;

  /**
   * Prefix advancement operator: move to the next iterator value satisfying
   * the predicate and return the new iterator value.
   */
  FilteredIterator &
  operator++();

  /**
   * Postfix advancement operator: move to the next iterator value satisfying
   * the predicate and return the old iterator value.
   */
  FilteredIterator
  operator++(int);

  /**
   * Prefix decrement operator: move to the previous iterator value satisfying
   * the predicate and return the new iterator value.
   */
  FilteredIterator &
  operator--();

  /**
   * Postfix advancement operator: move to the previous iterator value
   * satisfying the predicate and return the old iterator value.
   */
  FilteredIterator
  operator--(int);

  /**
   * Exception.
   */
  DeclException1(
    ExcInvalidElement,
    BaseIterator,
    << "The element " << arg1
    << " with which you want to compare or which you want to"
    << " assign from is invalid since it does not satisfy the predicate.");

private:
  /**
   * Base class to encapsulate a predicate object. Since predicates can be of
   * different types and we do not want to code these types into the template
   * parameter list of the filtered iterator class, we use a base class with
   * an abstract function and templatized derived classes that implement the
   * use of actual predicate types through the virtual function.
   *
   * @ingroup Iterators
   */
  class PredicateBase
  {
  public:
    /**
     * Mark the destructor virtual to allow destruction through pointers to
     * the base class.
     */
    virtual ~PredicateBase() = default;

    /**
     * Abstract function which in derived classes denotes the evaluation of
     * the predicate on the given iterator.
     */
    virtual bool
    operator()(const BaseIterator &bi) const = 0;

    /**
     * Generate a copy of this object, i.e. of the actual type of this
     * pointer.
     */
    virtual std::unique_ptr<PredicateBase>
    clone() const = 0;
  };


  /**
   * Actual implementation of the above abstract base class. Use a template
   * parameter to denote the actual type of the predicate and store a copy of
   * it. When the virtual function is called evaluate the given iterator with
   * the stored copy of the predicate.
   *
   * @ingroup Iterators
   */
  template <typename Predicate>
  class PredicateTemplate : public PredicateBase
  {
  public:
    /**
     * Constructor. Take a predicate and store a copy of it.
     */
    PredicateTemplate(const Predicate &predicate);

    /**
     * Evaluate the iterator with the stored copy of the predicate.
     */
    virtual bool
    operator()(const BaseIterator &bi) const override;

    /**
     * Generate a copy of this object, i.e. of the actual type of this
     * pointer.
     */
    virtual std::unique_ptr<PredicateBase>
    clone() const override;

  private:
    /**
     * Copy of the predicate.
     */
    const Predicate predicate;
  };

  /**
   * Pointer to an object that encapsulated the actual data type of the
   * predicate given to the constructor.
   */
  std::unique_ptr<const PredicateBase> predicate;
};



/**
 * Create an object of type FilteredIterator given the base iterator and
 * predicate.  This function makes the creation of temporary objects (for
 * example as function arguments) a lot simpler because one does not have to
 * explicitly specify the type of the base iterator by hand -- it is deduced
 * automatically here.
 *
 * @relatesalso FilteredIterator
 *
 * @dealiiConceptRequires{(std::predicate<Predicate, BaseIterator>)}
 */
template <typename BaseIterator, typename Predicate>
DEAL_II_CXX20_REQUIRES((std::predicate<Predicate, BaseIterator>))
FilteredIterator<BaseIterator> make_filtered_iterator(const BaseIterator &i,
                                                      const Predicate    &p)
{
  FilteredIterator<BaseIterator> fi(p);
  fi.set_to_next_positive(i);
  return fi;
}



namespace internal
{
  namespace FilteredIteratorImplementation
  {
    // The following classes create a nested sequence of
    // FilteredIterator<FilteredIterator<...<BaseIterator>...>> with as many
    // levels of FilteredIterator classes as there are elements in the TypeList
    // if the latter is given as a std::tuple<Args...>.
    template <typename BaseIterator, typename TypeList>
    struct NestFilteredIterators;

    template <typename BaseIterator, typename Predicate>
    struct NestFilteredIterators<BaseIterator, std::tuple<Predicate>>
    {
      using type = ::dealii::FilteredIterator<BaseIterator>;
    };

    template <typename BaseIterator, typename Predicate, typename... Targs>
    struct NestFilteredIterators<BaseIterator, std::tuple<Predicate, Targs...>>
    {
      using type = ::dealii::FilteredIterator<
        typename NestFilteredIterators<BaseIterator,
                                       std::tuple<Targs...>>::type>;
    };
  } // namespace FilteredIteratorImplementation
} // namespace internal



/**
 * Filter the  given range of iterators using a Predicate. This allows to
 * replace:
 * @code
 *   DoFHandler<dim> dof_handler;
 *   ...
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     {
 *       if (cell->is_locally_owned())
 *         {
 *           fe_values.reinit (cell);
 *           ...do the local integration on 'cell'...;
 *         }
 *     }
 * @endcode
 * by:
 * @code
 *   DoFHandler<dim> dof_handler;
 *   ...
 *   const auto filtered_iterators_range =
 *     filter_iterators(dof_handler.active_cell_iterators(),
 *                      IteratorFilters::LocallyOwnedCell());
 *   for (const auto &cell : filtered_iterators_range)
 *     {
 *       fe_values.reinit (cell);
 *       ...do the local integration on 'cell'...;
 *     }
 * @endcode
 *
 * @note This function is obsolete and should be replaced by calling
 *   `operator|` instead. The latter is certainly more in the
 *   spirit of
 *   [C++20 range
 * adaptors](https://en.cppreference.com/w/cpp/ranges/filter_view) and results
 * in the following code instead of the one shown above:
 *   @code
 *   DoFHandler<dim> dof_handler;
 *   ...
 *   for (const auto &cell :
 *          dof_handler.active_cell_iterators()
 *          | IteratorFilters::LocallyOwnedCell())
 *     {
 *       fe_values.reinit (cell);
 *       ...do the local integration on 'cell'...;
 *     }
 *   @endcode
 *
 * @relatesalso FilteredIterator
 * @ingroup CPP11
 *
 * @dealiiConceptRequires{(std::predicate<Predicate, BaseIterator>)}
 */
template <typename BaseIterator, typename Predicate>
DEAL_II_CXX20_REQUIRES((std::predicate<Predicate, BaseIterator>))
inline IteratorRange<FilteredIterator<BaseIterator>> filter_iterators(
  IteratorRange<BaseIterator> i,
  const Predicate            &p)
{
  FilteredIterator<BaseIterator> fi(p, *(i.begin()));
  FilteredIterator<BaseIterator> fi_end(p, *(i.end()));

  return IteratorRange<FilteredIterator<BaseIterator>>(fi, fi_end);
}



/**
 * Filter the given range of iterators through an arbitrary number of
 * Predicates. This allows to replace:
 * @code
 *   DoFHandler<dim> dof_handler;
 *   ...
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     {
 *       if (cell->is_locally_owned())
 *         {
 *           if (cell->at_boundary())
 *             {
 *               fe_values.reinit (cell);
 *               ...do the local integration on 'cell'...;
 *             }
 *         }
 *     }
 * @endcode
 * by:
 * @code
 *   DoFHandler<dim> dof_handler;
 *   ...
 *   const auto filtered_iterators_range =
 *     filter_iterators(dof_handler.active_cell_iterators(),
 *                      IteratorFilters::LocallyOwnedCell(),
 *                      IteratorFilters::AtBoundary());
 *   for (const auto &cell : filter_iterators_range)
 *     {
 *       fe_values.reinit (cell);
 *       ...do the local integration on 'cell'...;
 *     }
 * @endcode
 *
 * @note This function is obsolete and should be replaced by calling
 *   `operator|` once or multiple times instead.
 *   The latter is certainly more in the spirit of
 *   [C++20 range
 * adaptors](https://en.cppreference.com/w/cpp/ranges/filter_view) and results
 * in the following code instead of the one shown above:
 *   @code
 *   DoFHandler<dim> dof_handler;
 *   ...
 *   for (const auto &cell :
 *          dof_handler.active_cell_iterators()
 *          | IteratorFilters::LocallyOwnedCell()
 *          | IteratorFilters::AtBoundary())
 *     {
 *       fe_values.reinit (cell);
 *       ...do the local integration on 'cell'...;
 *     }
 *   @endcode
 *
 * @relatesalso FilteredIterator
 * @ingroup CPP11
 *
 * @dealiiConceptRequires{(std::predicate<Predicate, BaseIterator>)}
 */
template <typename BaseIterator, typename Predicate, typename... Targs>
DEAL_II_CXX20_REQUIRES((std::predicate<Predicate, BaseIterator>))
IteratorRange<
  typename internal::FilteredIteratorImplementation::
    NestFilteredIterators<BaseIterator, std::tuple<Predicate, Targs...>>::
      type> filter_iterators(IteratorRange<BaseIterator> i,
                             const Predicate            &p,
                             const Targs... args)
{
  // Recursively create filtered iterators, one predicate at a time
  auto fi = filter_iterators(i, p);
  return filter_iterators(fi, args...);
}



/**
 * Filter the  given range of iterators using a predicate. This allows to
 * replace:
 * @code
 *   DoFHandler<dim> dof_handler;
 *   ...
 *   for (const auto &cell : dof_handler.active_cell_iterators())
 *     {
 *       if (cell->is_locally_owned())
 *         {
 *           fe_values.reinit (cell);
 *           ...do the local integration on 'cell'...;
 *         }
 *     }
 * @endcode
 * by:
 * @code
 *   DoFHandler<dim> dof_handler;
 *   ...
 *   const auto filtered_iterators_range =
 *     filter_iterators();
 *   for (const auto &cell :
 *           dof_handler.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
 *     {
 *       fe_values.reinit (cell);
 *       ...do the local integration on 'cell'...;
 *     }
 * @endcode
 * Here, the `operator|` is to be interpreted in the same way as is done in
 * the [range adaptors](https://en.cppreference.com/w/cpp/ranges) feature
 * that is part of [C++20](https://en.wikipedia.org/wiki/C%2B%2B20). It has
 * the same meaning as the `|` symbol on the command line: It takes what is
 * on its left as its inputs, and filters and transforms to produce some
 * output. In the example above, it "filters" all of the active cell iterators
 * and removes those that do not satisfy the predicate -- that is, it produces
 * a range of iterators that only contains those cells that are both active
 * and locally owned.
 *
 * @note The `|` operator can be applied more than once. This results in code
 *   such as the following:
 *   @code
 *   DoFHandler<dim> dof_handler;
 *   ...
 *   for (const auto &cell :
 *          dof_handler.active_cell_iterators()
 *          | IteratorFilters::LocallyOwnedCell()
 *          | IteratorFilters::AtBoundary())
 *     {
 *       fe_values.reinit (cell);
 *       ...do the local integration on 'cell'...;
 *     }
 *   @endcode
 *   In this code, the loop executes over all active cells that are both
 *   locally owned and located at the boundary.
 *
 * @relatesalso FilteredIterator
 * @ingroup CPP11
 *
 * @dealiiConceptRequires{(std::predicate<Predicate, BaseIterator>)}
 */
template <typename BaseIterator, typename Predicate>
DEAL_II_CXX20_REQUIRES((std::predicate<Predicate, BaseIterator>))
inline IteratorRange<FilteredIterator<BaseIterator>>
operator|(IteratorRange<BaseIterator> i, const Predicate &p)
{
  return filter_iterators(i, p);
}



/* ------------------ Inline functions and templates ------------ */


template <typename BaseIterator>
template <typename Predicate>
DEAL_II_CXX20_REQUIRES((std::predicate<Predicate, BaseIterator>))
inline FilteredIterator<BaseIterator>::FilteredIterator(Predicate p)
  : predicate(new PredicateTemplate<Predicate>(p))
{}



template <typename BaseIterator>
template <typename Predicate>
DEAL_II_CXX20_REQUIRES((std::predicate<Predicate, BaseIterator>))
inline FilteredIterator<BaseIterator>::FilteredIterator(Predicate           p,
                                                        const BaseIterator &bi)
  : BaseIterator(bi)
  , predicate(new PredicateTemplate<Predicate>(p))
{
  if ((this->state() == IteratorState::valid) && !(*predicate)(*this))
    set_to_next_positive(bi);
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator>::FilteredIterator(
  const FilteredIterator &fi)
  : // this construction looks strange, but without going through the
    // address of fi, GCC would not cast fi to the base class of type
    // BaseIterator but tries to go through constructing a new
    // BaseIterator with an Accessor.
  BaseIterator(*static_cast<const BaseIterator *>(&fi))
  , predicate(fi.predicate->clone())
{}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::operator=(const FilteredIterator &fi)
{
  // Using equivalent code to the one for 'operator=(const BaseIterator &bi)'
  // below, some compiler would not cast fi to the base class of type
  // BaseIterator but try to go through constructing a new Accessor from fi
  // which fails. Hence, we just use an explicit upcast and call the above-
  // mentioned method.
  const BaseIterator &bi = fi;
  return operator=(bi);
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::operator=(const BaseIterator &bi)
{
  Assert((bi.state() != IteratorState::valid) || (*predicate)(bi),
         ExcInvalidElement(bi));
  BaseIterator::operator=(bi);
  return *this;
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::set_to_next_positive(const BaseIterator &bi)
{
  BaseIterator::operator=(bi);
  while ((this->state() == IteratorState::valid) && (!(*predicate)(*this)))
    BaseIterator::operator++();

  return *this;
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::set_to_previous_positive(const BaseIterator &bi)
{
  BaseIterator::operator=(bi);
  while ((this->state() == IteratorState::valid) && (!(*predicate)(*this)))
    BaseIterator::operator--();

  return *this;
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator==(const FilteredIterator &fi) const
{
  return (static_cast<const BaseIterator &>(*this) ==
          static_cast<const BaseIterator &>(fi));
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator!=(const FilteredIterator &fi) const
{
  return (static_cast<const BaseIterator &>(*this) !=
          static_cast<const BaseIterator &>(fi));
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator<(const FilteredIterator &fi) const
{
  return (static_cast<const BaseIterator &>(*this) <
          static_cast<const BaseIterator &>(fi));
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator==(const BaseIterator &bi) const
{
  return (static_cast<const BaseIterator &>(*this) == bi);
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator!=(const BaseIterator &bi) const
{
  return (static_cast<const BaseIterator &>(*this) != bi);
}



template <typename BaseIterator>
inline bool
FilteredIterator<BaseIterator>::operator<(const BaseIterator &bi) const
{
  return (static_cast<const BaseIterator &>(*this) < bi);
}


template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::operator++()
{
  if (this->state() == IteratorState::valid)
    do
      BaseIterator::operator++();
    while ((this->state() == IteratorState::valid) && !(*predicate)(*this));
  return *this;
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator>
FilteredIterator<BaseIterator>::operator++(int)
{
  const FilteredIterator old_state = *this;

  if (this->state() == IteratorState::valid)
    do
      BaseIterator::operator++();
    while ((this->state() == IteratorState::valid) && !(*predicate)(*this));
  return old_state;
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::operator--()
{
  if (this->state() == IteratorState::valid)
    do
      BaseIterator::operator--();
    while ((this->state() == IteratorState::valid) && !(*predicate)(*this));
  return *this;
}



template <typename BaseIterator>
inline FilteredIterator<BaseIterator>
FilteredIterator<BaseIterator>::operator--(int)
{
  const FilteredIterator old_state = *this;

  if (this->state() == IteratorState::valid)
    do
      BaseIterator::operator--();
    while ((this->state() == IteratorState::valid) && !(*predicate)(*this));
  return old_state;
}



template <typename BaseIterator>
template <typename Predicate>
inline FilteredIterator<BaseIterator>::PredicateTemplate<
  Predicate>::PredicateTemplate(const Predicate &predicate)
  : predicate(predicate)
{}



template <typename BaseIterator>
template <typename Predicate>
bool
FilteredIterator<BaseIterator>::PredicateTemplate<Predicate>::operator()(
  const BaseIterator &bi) const
{
  return predicate(bi);
}



template <typename BaseIterator>
template <typename Predicate>
std::unique_ptr<typename FilteredIterator<BaseIterator>::PredicateBase>
FilteredIterator<BaseIterator>::PredicateTemplate<Predicate>::clone() const
{
  return std::make_unique<PredicateTemplate>(predicate);
}



namespace IteratorFilters
{
  // ---------------- IteratorFilters::Active ---------

  template <class Iterator>
  inline bool
  Active::operator()(const Iterator &i) const
  {
    return i->is_active();
  }


  // ---------------- IteratorFilters::UserFlagSet ---------

  template <class Iterator>
  inline bool
  UserFlagSet::operator()(const Iterator &i) const
  {
    return (i->user_flag_set());
  }


  // ---------------- IteratorFilters::UserFlagNotSet ---------

  template <class Iterator>
  inline bool
  UserFlagNotSet::operator()(const Iterator &i) const
  {
    return (!i->user_flag_set());
  }


  // ---------------- IteratorFilters::LevelEqualTo ---------
  inline LevelEqualTo::LevelEqualTo(const unsigned int level)
    : level(level)
  {}



  template <class Iterator>
  inline bool
  LevelEqualTo::operator()(const Iterator &i) const
  {
    return (static_cast<unsigned int>(i->level()) == level);
  }



  // ---------------- IteratorFilters::SubdomainEqualTo ---------
  inline SubdomainEqualTo::SubdomainEqualTo(
    const types::subdomain_id subdomain_id)
    : subdomain_id(subdomain_id)
  {}



  template <class Iterator>
  inline bool
  SubdomainEqualTo::operator()(const Iterator &i) const
  {
    return (i->subdomain_id() == subdomain_id);
  }



  // ---------------- IteratorFilters::LocallyOwnedCell ---------

  template <class Iterator>
  inline bool
  LocallyOwnedCell::operator()(const Iterator &i) const
  {
    return (i->is_locally_owned());
  }


  // ---------------- IteratorFilters::LocallyOwnedLevelCell ---------

  template <class Iterator>
  inline bool
  LocallyOwnedLevelCell::operator()(const Iterator &i) const
  {
    return (i->is_locally_owned_on_level());
  }



  // ---------------- IteratorFilters::MaterialIdEqualTo ---------
  inline MaterialIdEqualTo::MaterialIdEqualTo(
    const types::material_id material_id,
    const bool               only_locally_owned)
    : material_ids{material_id}
    , only_locally_owned(only_locally_owned)
  {}



  inline MaterialIdEqualTo::MaterialIdEqualTo(
    const std::set<types::material_id> &material_ids,
    const bool                          only_locally_owned)
    : material_ids(material_ids)
    , only_locally_owned(only_locally_owned)
  {}



  template <class Iterator>
  inline bool
  MaterialIdEqualTo::operator()(const Iterator &i) const
  {
    return only_locally_owned == true ?
             (material_ids.find(i->material_id()) != material_ids.end() &&
              i->is_locally_owned()) :
             material_ids.find(i->material_id()) != material_ids.end();
  }



  // ---------------- IteratorFilters::ActiveFEIndexEqualTo ---------
  inline ActiveFEIndexEqualTo::ActiveFEIndexEqualTo(
    const unsigned int active_fe_index,
    const bool         only_locally_owned)
    : active_fe_indices{active_fe_index}
    , only_locally_owned(only_locally_owned)
  {}



  inline ActiveFEIndexEqualTo::ActiveFEIndexEqualTo(
    const std::set<unsigned int> &active_fe_indices,
    const bool                    only_locally_owned)
    : active_fe_indices(active_fe_indices)
    , only_locally_owned(only_locally_owned)
  {}



  template <class Iterator>
  inline bool
  ActiveFEIndexEqualTo::operator()(const Iterator &i) const
  {
    return only_locally_owned == true ?
             (i->is_locally_owned() &&
              active_fe_indices.find(i->active_fe_index()) !=
                active_fe_indices.end()) :
             active_fe_indices.find(i->active_fe_index()) !=
               active_fe_indices.end();
  }



  // ---------------- IteratorFilters::AtBoundary ---------

  template <class Iterator>
  inline bool
  AtBoundary::operator()(const Iterator &i) const
  {
    return (i->at_boundary());
  }



  // ---------------- IteratorFilters::BoundaryIdEqualTo ---------
  inline BoundaryIdEqualTo::BoundaryIdEqualTo(
    const types::boundary_id boundary_id)
    : boundary_ids{boundary_id}
  {}



  inline BoundaryIdEqualTo::BoundaryIdEqualTo(
    const std::set<types::boundary_id> &boundary_ids)
    : boundary_ids(boundary_ids)
  {}



  template <class Iterator>
  inline bool
  BoundaryIdEqualTo::operator()(const Iterator &i) const
  {
    return boundary_ids.find(i->boundary_id()) != boundary_ids.end();
  }



  // ---------------- IteratorFilters::ManifoldIdEqualTo ---------
  inline ManifoldIdEqualTo::ManifoldIdEqualTo(
    const types::manifold_id manifold_id)
    : manifold_ids{manifold_id}
  {}



  inline ManifoldIdEqualTo::ManifoldIdEqualTo(
    const std::set<types::manifold_id> &manifold_ids)
    : manifold_ids(manifold_ids)
  {}



  template <class Iterator>
  inline bool
  ManifoldIdEqualTo::operator()(const Iterator &i) const
  {
    return manifold_ids.find(i->manifold_id()) != manifold_ids.end();
  }
} // namespace IteratorFilters


DEAL_II_NAMESPACE_CLOSE

#endif
