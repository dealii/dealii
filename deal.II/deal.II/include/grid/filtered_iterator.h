//------------------------ filtered_iterator.h -----------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//------------------------ filtered_iterator.h -----------------------
#ifndef __deal2__filtered_iterator_h
#define __deal2__filtered_iterator_h



/**
 * In this namespace a number of classes is declared that may be used
 * as filters in the @ref{FilteredIterator} class. The filters either
 * check for binary information (for example, the @p{Active} filter
 * class checks whether the object pointed to is active), or for
 * valued information by comparison with prescribed values (for
 * example, the @ref{LevelEqualTo} filter class checks whether the
 * level of the object pointed to by the iterator under consideration
 * is equal to a value that was given to the filter upon construction.
 *
 * For examples of use of these classes as well as requirements on
 * filters see the general description of the @ref{FilteredIterator}
 * class.
 *
 * @author Wolfgang Bangerth, 2002
 */
namespace IteratorFilters
{
				   /**
				    * Filter that evaluates to true if
				    * either the iterator points to an
				    * active object or an iterator
				    * past the end.
				    */
  class Active 
  {
    public:
				       /**
					* Evaluate the iterator and
					* return true if the object is
					* active or past the end.
					*/
      template <class Iterator>
      bool operator () (const Iterator &i) const;
  };
  
				   /**
				    * Filter that evaluates to true if
				    * either the iterator points to an
				    * object for which the user flag
				    * is set or an iterator past the
				    * end.
				    */
  class UserFlagSet
  {
    public:
				       /**
					* Evaluate the iterator and
					* return true if the object
					* has a set user flag or past
					* the end.
					*/
      template <class Iterator>
      bool operator () (const Iterator &i) const;
  };
  

				   /**
				    * Filter that evaluates to true if
				    * either the iterator points to an
				    * object for which the user flag
				    * is not set or an iterator past
				    * the end. Inverse filter to the
				    * previous class.
				    */
  class UserFlagNotSet
  {
    public:
				       /**
					* Evaluate the iterator and
					* return true if the object
					* has an unset user flag or
					* past the end.
					*/
      template <class Iterator>
      bool operator () (const Iterator &i) const;
  };

  
				   /**
				    * Filter for iterators that
				    * evaluates to true if either the
				    * iterator is past the end or the
				    * level of the object pointed to
				    * is equal to a value given to the
				    * constructor.
				    */
  class LevelEqualTo 
  {
    public:
				       /**
					* Constructor. Store the level
					* which iterators shall have
					* to be evaluated to true.
					*/
      LevelEqualTo (const unsigned int level);

				       /**
					* Evaluation operator. Returns
					* true if either the level of
					* the object pointed to is
					* equal to the stored value or
					* the iterator is past the
					* end.
					*/
      template <class Iterator>
      bool operator () (const Iterator &i) const;

    protected:
				       /**
					* Stored value to compare the
					* level with.
					*/
      const unsigned int level;
  };



				   /**
				    * Filter for iterators that
				    * evaluates to true if either the
				    * iterator is past the end or the
				    * subdomain id of the object
				    * pointed to is equal to a value
				    * given to the constructor,
				    * assuming that the iterator
				    * allows querying for a subdomain
				    * id).
				    */
  class SubdomainEqualTo 
  {
    public:
				       /**
					* Constructor. Store the
					* subdomain which iterators
					* shall have to be evaluated
					* to true.
					*/
      SubdomainEqualTo (const unsigned int subdomain_id);

				       /**
					* Evaluation operator. Returns
					* true if either the subdomain
					* of the object pointed to is
					* equal to the stored value or
					* the iterator is past the
					* end.
					*/
      template <class Iterator>
      bool operator () (const Iterator &i) const;

    protected:
				       /**
					* Stored value to compare the
					* subdomain with.
					*/
      const unsigned int subdomain_id;
  };
};


/**
 * This class provides a certain view on a range of triangulation or
 * DoFHandler iterators by only iterating over elements that satisfy a
 * given filter (called a @em{predicate}, following the notation of
 * the C++ standard library). Once initialized with a predicate and a
 * value for the iterator, a filtered iterator hops to the next or
 * previous element that satisfies the predicate if operators ++ or --
 * are invoked. Intermediate iterator values that lie in between but
 * do not satisfy the predicate are skipped. It is thus very simple to
 * write loops over a certain class of objects without the need to
 * explicitely write down the condition they have to satisfy in each
 * loop iteration. This in particular is helpful if functions are
 * called with a pair of iterators denoting a range on which they
 * shall act, by choosing a filtered iterator instead of usual ones.
 *
 *
 * @sect3{Predicates}
 *
 * The object that represent the condition an iterator has to satisfy
 * only have to provide an interface that allows to call the
 * evaluation operator, i.e. @p{operator()}. This includes function
 * pointers as well as classes that implement an @p{operator ()}.
 *
 * An example of a simple valid predicate is the following: given the function
 * @begin{verbatim}
 *   template <typename Iterator>
 *   bool level_equal_to_3 (const Iterator c)
 *   {
 *     return (static_cast<unsigned int>(c->level()) == 3);
 *   };
 * @end{verbatim} 
 * then
 * @begin{verbatim}
 *   &level_equal_to_3<typename Triangulation<dim>::active_cell_iterator>
 * @end{verbatim}
 * is a valid predicate.
 *
 * Likewise, given the following binary function
 * @begin{verbatim}
 *   template <typename Iterator>
 *   bool level_equal_to (const Iterator     c,
 *                        const unsigned int level)
 *   {
 *     return (static_cast<unsigned int>(c->level()) == level);
 *   };
 * @end{verbatim} 
 * then
 * @begin{verbatim}
 *   std::bind2nd (std::ptr_fun(&level_equal_to<active_cell_iterator>), 3)
 * @end{verbatim}
 * is another valid predicate (here: a function that returns true if
 * either the iterator is past the end or the level is equal to the
 * second argument; this second argument is bound to a fixed value
 * using the @p{std::bind2nd} function).
 *
 * Finally, classes can be predicates. The following class is one:
 * @begin{verbatim}
 *   class Active 
 *   {
 *     public:
 *       template <class Iterator>
 *       bool operator () (const Iterator &i) const {
 *         return (i->active());
 *       }
 *   };
 * @end{verbatim}
 * and objects of this type can be used as predicates. Likewise, this
 * more complicated one can also be used:
 * @begin{verbatim}
 *   class SubdomainEqualTo
 *   {
 *     public:
 *       SubdomainEqualTo (const unsigned int subdomain_id)
 *                   : subdomain_id (subdomain_id) {};
 *
 *       template <class Iterator>
 *       bool operator () (const Iterator &i) const {
 *         return (i->subdomain_id() == subdomain_id);
 *       }
 *
 *     private:
 *       const unsigned int subdomain_id;
 *   };
 * @end{verbatim}
 * Objects like @p{SubdomainEqualTo(3)} can then be used as predicates.
 *
 * Since whenever a predicate is evaluated it is checked that the
 * iterator checked is actually valid (i.e. not past the end), no
 * checks for this case have to be performed inside predicates.
 *
 * A number of filter classes are already implemented in the
 * @ref{IteratorFilters} namespace, but writing different ones is
 * simple following the examples above.
 *
 *
 * @sect3{Initialization of filtered iterators}
 *
 * Filtered iterators are given a predicate at construction time which
 * cannot be changed any more. This behaviour would be expected if the
 * predicate would have been given as a template parameter to the
 * class, but since that would make the declaration of filtered
 * iterators a nightmare, we rather give the predicate as an
 * unchangeable entity to the constructor. Note that one can assign a
 * filtered iterator with one predicate to another filtered iterator
 * with another type; yet, this does @em{not} change the predicate of
 * the assigned-to iterator, only the pointer indicating the iterator
 * is changed.
 *
 * If a filtered iterator is not assigned a value of the underlying
 * (unfiltered) iterator type, the default value is taken. If,
 * however, a value is given to the constructor, that value has either
 * to be past the end, or has to satisfy the predicate. For example,
 * if the predicate only evaluates to true if the level of an object
 * is equal to three, then @p{tria.begin_active(3)} would be a valid
 * choice while @p{tria.begin()} would not since the latter also
 * returns iterators to non-active cells which always start at level
 * 0.
 *
 * Since one often only has some iterator and wants to set a filtered
 * iterator to the first one that satisfies a predicate (for example,
 * the first one for which the user flag is set, or the first one with
 * a given subdomain id), there are assignement functions
 * @p{set_to_next_positive} and @p{set_to_previous_positive} that
 * assign the next or last previous iterator that satisfies the
 * predicate, i.e. they follow the list of iterators in either
 * direction until they find a matching one (or the past-the-end
 * iterator). Like the @p{operator=} they return the resulting value
 * of the filtered iterator.
 *
 *
 * @sect3{Examples}
 *
 * The following call counts the number of active cells that
 * have a set user flag:
 * @begin{verbatim}
 *   FilteredIterator<typename Triangulation<dim>::active_cell_iterator>
 *      begin (IteratorFilters::UserFlagSet()),
 *      end (IteratorFilters::UserFlagSet());
 *   begin.set_to_next_positive(tria.begin_active());
 *   end = tria.end();
 *   n_flagged_cells = std::distance (begin, end);
 * @begin{verbatim}
 * Note that by the @p{set_to_next_positive} call the first cell with
 * a set user flag was assigned to the @p{begin} iterator. For the
 * @{end} iterator, no such call was necessary, since the past-the-end
 * iterator always satisfies all predicates.
 *
 * The same can be achieved by the following snippet, though harder to read:
 * @begin{verbatim}
 *   typedef FilteredIterator<typename Triangulation<dim>::active_cell_iterator> FI;
 *   n_flagged_cells =
 *      std::distance (FI(IteratorFilters::UserFlagSet())
 *                            .set_to_next_positive(tria.begin_active()),
 *                     FI(IteratorFilters::UserFlagSet(), tria.end()));
 * @begin{verbatim}
 * It relies on the fact that if we create an unnamed filtered
 * iterator with a given predicate but no iterator value and assign it
 * the next positive value with respect to this predicate, it returns
 * itself which is then used as the first parameter to the
 * @p{std::distance} function. This procedure is not necessary for the
 * end element to this function here, since the past-the-end iterator
 * always satisfies the predicate so that we can assign this value to
 * the filtered iterator directly in the constructor.
 *
 * Finally, the following loop only assembles the matrix on cells with
 * subdomain id equal to three:
 * @begin{verbatim}
 * FilteredIterator<typename Triangulation<dim>::active_cell_iterator>
 *   cell (FilteredIterator::SubdomainEqualTo(3)),
 *   endc (FilteredIterator::SubdomainEqualTo(3), tria.end());
 * cell.set_to_next_positive (tria.begin_active());
 * for (; cell!=endc; ++cell)
 *   assemble_local_matrix (cell);
 * @end{verbatim}
 *
 * Since comparison between filtered and unfiltered iterators is
 * defined, we could as well have let the @p{endc} variable in the
 * last example be of type
 * @p{Triangulation<dim>::active_cell_iterator} since it is unchanged
 * and its value does not depend on the filter.
 *
 * @author Wolfgang Bangerth, 2002
 */
template <typename BaseIterator>
class FilteredIterator : public BaseIterator
{
  public:
				     /**
				      * Typedef to the accessor type
				      * of the underlying iterator.
				      */
    typedef typename BaseIterator::AccessorType AccessorType;

				     /**
				      * Constructor. Set the iterator
				      * to the default state and use
				      * the given predicate for
				      * filtering subsequent
				      * assignement and iteration.
				      */
    template <typename Predicate>
    FilteredIterator (Predicate p);

				     /**
				      * Constructor. Use the given
				      * predicate for filtering and
				      * initialize the iterator with
				      * the given value. This initial
				      * value has to satisfy the
				      * predicate.
				      */
    template <typename Predicate>    
    FilteredIterator (Predicate           p,
		      const BaseIterator &bi);

				     /**
				      * Copy constructor. Copy the
				      * predicate and iterator value
				      * of the given argument.
				      */
    FilteredIterator (const FilteredIterator &fi);

				     /**
				      * Destructor.
				      */
    ~FilteredIterator ();
    
				     /**
				      * Assignment operator. Copy the
				      * iterator value of the
				      * argument, but as discussed in
				      * the class documentation, the
				      * predicate of the argument is
				      * not copied. The iterator value
				      * underlying the argument has to
				      * satisfy the predicate of the
				      * object assigned to, as given
				      * at its construction time.
				      */
    FilteredIterator & operator = (const FilteredIterator &fi);

				     /**
				      * Assignment operator. Copy the
				      * iterator value of the
				      * argument, and keep the
				      * predicate of this object. The
				      * given iterator value has to
				      * satisfy the predicate of the
				      * object assigned to, as given
				      * at its construction time.
				      */
    FilteredIterator & operator = (const BaseIterator &fi);

				     /**
				      * Search for the next iterator
				      * from @p{bi} onwards that
				      * satisfies the predicate of
				      * this object and assign it to
				      * this object.
				      *
				      * Since filtered iterators are
				      * automatically converted to the
				      * underlying base iterator type,
				      * you can also give a filtered
				      * iterator as argument to this
				      * function.
				      */
    FilteredIterator &
    set_to_next_positive (const BaseIterator &bi);
    
				     /**
				      * As above, but search for the
				      * previous iterator from @p{bi}
				      * backwards that satisfies the
				      * predicate of this object and
				      * assign it to this object.
				      *
				      * Since filtered iterators are
				      * automatically converted to the
				      * underlying base iterator type,
				      * you can also give a filtered
				      * iterator as argument to this
				      * function.
				      */
    FilteredIterator &
    set_to_previous_positive (const BaseIterator &bi);
    
				     /**
				      * Compare for equality of the
				      * underlying iterator values of
				      * this and the given object.
				      *
				      * We do not compare for equality
				      * of the predicates.
				      */
    bool operator == (const FilteredIterator &fi) const;

				     /**
				      * Compare for equality of the
				      * underlying iterator value of
				      * this object with the given
				      * object.
				      *
				      * The predicate of this object
				      * is irrelevant for this
				      * operation.
				      */
    bool operator == (const BaseIterator &fi) const;

				     /**
				      * Compare for inequality of the
				      * underlying iterator values of
				      * this and the given object.
				      *
				      * We do not compare for equality
				      * of the predicates.
				      */
    bool operator != (const FilteredIterator &fi) const;

				     /**
				      * Compare for inequality of the
				      * underlying iterator value of
				      * this object with the given
				      * object.
				      *
				      * The predicate of this object
				      * is irrelevant for this
				      * operation.
				      */
    bool operator != (const BaseIterator &fi) const;

				     /**
				      * Compare for ordering of the
				      * underlying iterator values of
				      * this and the given object.
				      *
				      * We do not compare the
				      * predicates.
				      */    
    bool operator <  (const FilteredIterator &fi) const;

				     /**
				      * Compare for ordering of the
				      * underlying iterator value of
				      * this object with the given
				      * object.
				      *
				      * The predicate of this object
				      * is irrelevant for this
				      * operation.
				      */
    bool operator <  (const BaseIterator &fi) const;

				     /**
				      * Prefix advancement operator:
				      * move to the next iterator
				      * value satisfying the predicate
				      * and return the new iterator
				      * value.
				      */
    FilteredIterator & operator ++ ();

				     /**
				      * Postfix advancement operator:
				      * move to the next iterator
				      * value satisfying the predicate
				      * and return the old iterator
				      * value.
				      */
    FilteredIterator   operator ++ (int);

				     /**
				      * Prefix decrement operator:
				      * move to the previous iterator
				      * value satisfying the predicate
				      * and return the new iterator
				      * value.
				      */
    FilteredIterator & operator -- ();

				     /**
				      * Postfix advancement operator:
				      * move to the previous iterator
				      * value satisfying the predicate
				      * and return the old iterator
				      * value.
				      */
    FilteredIterator   operator -- (int);

				     /**
				      * Exception.
				      */
    DeclException1 (ExcInvalidElement,
		    BaseIterator,
		    << "The element " << arg1
		    << " with which you want to compare or which you want to"
		    << " assign is invalid since it does not satisfy the predicate.");
    
  private:

				     /**
				      * Base class to encapsulate a
				      * predicate object. Since
				      * predicates can be of different
				      * types and we do not want to
				      * code these types into the
				      * template parameter list of the
				      * filtered iterator class, we
				      * use a base class with an
				      * abstract function and
				      * templatized derived classes
				      * that implement the use of
				      * actual predicate types through
				      * the virtual function.
				      */
    class PredicateBase
    {
      public:
					 /**
					  * Mark the destructor
					  * virtual to allow
					  * destruction through
					  * pointers to the base
					  * class.
					  */
	virtual ~PredicateBase () {};

					 /**
					  * Abstract function which in
					  * derived classes denotes
					  * the evaluation of the
					  * predicate on the give
					  * iterator.
					  */
	virtual bool operator () (const BaseIterator &bi) const = 0;

					 /**
					  * Generate a copy of this
					  * object, i.e. of the actual
					  * type of this pointer.
					  */
	virtual PredicateBase * clone () const = 0;
    };


				     /**
				      * Actual implementation of the
				      * above abstract base class. Use
				      * a template parameter to denote
				      * the actual type of the
				      * predicate and store a copy of
				      * it. When the virtual function
				      * is called evaluate the given
				      * iterator with the stored copy
				      * of the predicate.
				      */
    template <typename Predicate>
    class PredicateTemplate : public PredicateBase
    {
      public:
					 /**
					  * Constructor. Take a
					  * predicate and store a copy
					  * of it.
					  */
	PredicateTemplate (const Predicate &predicate);

					 /**
					  * Evaluate the iterator with
					  * the stored copy of the
					  * predicate.
					  */
	virtual bool operator () (const BaseIterator &bi) const;

					 /**
					  * Generate a copy of this
					  * object, i.e. of the actual
					  * type of this pointer.
					  */
	virtual PredicateBase * clone () const;
	
      private:
					 /**
					  * Copy of the predicate.
					  */
	const Predicate predicate;
    };

				     /**
				      * Pointer to an object that
				      * encapsulated the actual data
				      * type of the predicate given to
				      * the constructor.
				      */
    const PredicateBase * predicate;
};



/* ------------------ Inline functions and templates ------------ */


template <typename BaseIterator>
template <typename Predicate>
inline
FilteredIterator<BaseIterator>::
FilteredIterator (Predicate p)
		:
		predicate (new PredicateTemplate<Predicate>(p))
{};



template <typename BaseIterator>
template <typename Predicate>
inline
FilteredIterator<BaseIterator>::
FilteredIterator (Predicate          p,
		  const BaseIterator &bi)
		:
		BaseIterator (bi),
		predicate (new PredicateTemplate<Predicate>(p))
{
  Assert ((state() != IteratorState::valid) || (*predicate) (*this),
	  ExcInvalidElement(bi));
};



template <typename BaseIterator>
inline
FilteredIterator<BaseIterator>::
FilteredIterator (const FilteredIterator &fi)
		:
		BaseIterator (static_cast<BaseIterator>(fi)),
		predicate (fi.predicate->clone ())
{};



template <typename BaseIterator>
inline
FilteredIterator<BaseIterator>::
~FilteredIterator ()
{
  delete predicate;
  predicate = 0;
};



template <typename BaseIterator>
inline
FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::
operator = (const FilteredIterator &fi)
{
  Assert ((fi.state() != IteratorState::valid) || (*predicate)(fi),
	  ExcInvalidElement(fi));
  BaseIterator::operator = (fi);
  return *this;
};



template <typename BaseIterator>
inline
FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::
operator = (const BaseIterator &bi)
{
  Assert ((bi.state() != IteratorState::valid) || (*predicate)(bi),
	  ExcInvalidElement(bi));  
  BaseIterator::operator = (bi);
  return *this;
};



template <typename BaseIterator>
inline
FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::
set_to_next_positive (const BaseIterator &bi)
{
  BaseIterator::operator = (bi);
  while ((state() == IteratorState::valid) &&
	 ( ! (*predicate)(*this)))
    BaseIterator::operator++ ();
  
  return *this;
};



template <typename BaseIterator>
inline
FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::
set_to_previous_positive (const BaseIterator &bi)
{
  BaseIterator::operator = (bi);
  while ((state() == IteratorState::valid) &&
	 ( ! predicate(*this)))
    BaseIterator::operator-- ();
  
  return *this;
};



template <typename BaseIterator>
inline
bool
FilteredIterator<BaseIterator>::
operator == (const FilteredIterator &fi) const
{
  return (static_cast<const BaseIterator &>(*this)
	  ==
	  static_cast<const BaseIterator &>(fi));
};



template <typename BaseIterator>
inline
bool
FilteredIterator<BaseIterator>::
operator != (const FilteredIterator &fi) const
{
  return (static_cast<const BaseIterator &>(*this)
	  !=
	  static_cast<const BaseIterator &>(fi));
};



template <typename BaseIterator>
inline
bool
FilteredIterator<BaseIterator>::
operator < (const FilteredIterator &fi) const
{
  return (static_cast<const BaseIterator &>(*this)
	  <
	  static_cast<const BaseIterator &>(fi));
};




template <typename BaseIterator>
inline
bool
FilteredIterator<BaseIterator>::
operator == (const BaseIterator &bi) const
{
  return (static_cast<const BaseIterator &>(*this) == bi);
};



template <typename BaseIterator>
inline
bool
FilteredIterator<BaseIterator>::
operator != (const BaseIterator &bi) const
{
  return (static_cast<const BaseIterator &>(*this) != bi);
};



template <typename BaseIterator>
inline
bool
FilteredIterator<BaseIterator>::
operator < (const BaseIterator &bi) const
{
  return (static_cast<const BaseIterator &>(*this) < bi);
};


template <typename BaseIterator>
inline
FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::
operator ++ ()
{
  if (state() == IteratorState::valid)
    do
      BaseIterator::operator++ ();
    while ((state() == IteratorState::valid) &&
	   !(*predicate) (*this));
  return *this;
};



template <typename BaseIterator>
inline
FilteredIterator<BaseIterator>
FilteredIterator<BaseIterator>::
operator ++ (int)
{
  const FilteredIterator old_state = *this;
  
  if (state() == IteratorState::valid)
    do
      BaseIterator::operator++ ();
    while ((state() == IteratorState::valid) &&
	   !(*predicate) (*this));
  return old_state;
};




template <typename BaseIterator>
inline
FilteredIterator<BaseIterator> &
FilteredIterator<BaseIterator>::
operator -- ()
{
  if (state() == IteratorState::valid)
    do
      BaseIterator::operator-- ();
    while ((state() == IteratorState::valid) &&
	   !(*predicate) (*this));
  return *this;
};



template <typename BaseIterator>
inline
FilteredIterator<BaseIterator>
FilteredIterator<BaseIterator>::
operator -- (int)
{
  const FilteredIterator old_state = *this;
  
  if (state() == IteratorState::valid)
    do
      BaseIterator::operator-- ();
    while ((state() == IteratorState::valid) &&
	   !(*predicate) (*this));
  return old_state;
};



template <typename BaseIterator>
template <typename Predicate>
inline
FilteredIterator<BaseIterator>::PredicateTemplate<Predicate>::
PredicateTemplate (const Predicate &predicate)
		:
		predicate (predicate)
{};



template <typename BaseIterator>
template <typename Predicate>
bool
FilteredIterator<BaseIterator>::PredicateTemplate<Predicate>::
operator () (const BaseIterator &bi) const
{
  return predicate(bi);
};



template <typename BaseIterator>
template <typename Predicate>
FilteredIterator<BaseIterator>::PredicateBase *
FilteredIterator<BaseIterator>::PredicateTemplate<Predicate>::
clone () const
{
  return new PredicateTemplate (predicate);
};



namespace IteratorFilters 
{

// ---------------- IteratorFilters::Active ---------

  template <class Iterator>
  inline
  bool
  Active::operator () (const Iterator &i) const
  {
    return (i->active());
  };


// ---------------- IteratorFilters::UserFlagSet ---------

  template <class Iterator>
  inline
  bool
  UserFlagSet::operator () (const Iterator &i) const
  {
    return (i->user_flag_set());
  };


// ---------------- IteratorFilters::UserFlagNotSet ---------

  template <class Iterator>
  inline
  bool
  UserFlagNotSet::operator () (const Iterator &i) const
  {
    return (! i->user_flag_set());
  };


// ---------------- IteratorFilters::LevelEqualTo ---------  
  inline
  LevelEqualTo::LevelEqualTo (const unsigned int level)
		  :
		  level (level)
  {};



  template <class Iterator>
  inline
  bool
  LevelEqualTo::operator () (const Iterator &i) const
  {
    return (static_cast<unsigned int>(i->level()) == level);
  };



// ---------------- IteratorFilters::SubdomainEqualTo ---------
  inline
  SubdomainEqualTo::SubdomainEqualTo (const unsigned int subdomain_id)
		  :
		  subdomain_id (subdomain_id)
  {};



  template <class Iterator>
  inline
  bool
  SubdomainEqualTo::operator () (const Iterator &i) const
  {
    return (static_cast<unsigned int>(i->subdomain_id()) == subdomain_id);
  };
};


/*------------------------- filtered_iterator.h ------------------------*/
#endif
/*------------------------- filtered_iterator.h ------------------------*/


