/*
    Copyright 2005-2013 Intel Corporation.  All Rights Reserved.

    This file is part of Threading Building Blocks.

    Threading Building Blocks is free software; you can redistribute it
    and/or modify it under the terms of the GNU General Public License
    version 2 as published by the Free Software Foundation.

    Threading Building Blocks is distributed in the hope that it will be
    useful, but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Threading Building Blocks; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

    As a special exception, you may use this file as part of a free software
    library without restriction.  Specifically, if other files instantiate
    templates or use macros or inline functions from this file, or you compile
    this file and link it with other files to produce an executable, this
    file does not by itself cause the resulting executable to be covered by
    the GNU General Public License.  This exception does not however
    invalidate any other reasons why the executable file might be covered by
    the GNU General Public License.
*/

/* Container implementations in this header are based on PPL implementations
   provided by Microsoft. */

#ifndef __TBB_concurrent_unordered_set_H
#define __TBB_concurrent_unordered_set_H

#include "internal/_concurrent_unordered_impl.h"

namespace tbb
{

namespace interface5 {

// Template class for hash set traits
template<typename Key, typename Hash_compare, typename Allocator, bool Allow_multimapping>
class concurrent_unordered_set_traits
{
protected:
    typedef Key value_type;
    typedef Key key_type;
    typedef Hash_compare hash_compare;
    typedef typename Allocator::template rebind<value_type>::other allocator_type;
    enum { allow_multimapping = Allow_multimapping };

    concurrent_unordered_set_traits() : my_hash_compare() {}
    concurrent_unordered_set_traits(const hash_compare& hc) : my_hash_compare(hc) {}

    typedef hash_compare value_compare;

    static const Key& get_key(const value_type& value) {
        return value;
    }

    hash_compare my_hash_compare; // the comparator predicate for keys
};

template <typename Key, typename Hasher = tbb::tbb_hash<Key>, typename Key_equality = std::equal_to<Key>, typename Allocator = tbb::tbb_allocator<Key> >
class concurrent_unordered_set : public internal::concurrent_unordered_base< concurrent_unordered_set_traits<Key, internal::hash_compare<Key, Hasher, Key_equality>, Allocator, false> >
{
    // Base type definitions
    typedef internal::hash_compare<Key, Hasher, Key_equality> hash_compare;
    typedef internal::concurrent_unordered_base< concurrent_unordered_set_traits<Key, hash_compare, Allocator, false> > base_type;
    typedef concurrent_unordered_set_traits<Key, internal::hash_compare<Key, Hasher, Key_equality>, Allocator, false> traits_type;
    using traits_type::my_hash_compare;
#if __TBB_EXTRA_DEBUG
public:
#endif
    using traits_type::allow_multimapping;
public:
    using base_type::end;
    using base_type::find;
    using base_type::insert;

    // Type definitions
    typedef Key key_type;
    typedef typename base_type::value_type value_type;
    typedef Key mapped_type;
    typedef Hasher hasher;
    typedef Key_equality key_equal;
    typedef hash_compare key_compare;

    typedef typename base_type::allocator_type allocator_type;
    typedef typename base_type::pointer pointer;
    typedef typename base_type::const_pointer const_pointer;
    typedef typename base_type::reference reference;
    typedef typename base_type::const_reference const_reference;

    typedef typename base_type::size_type size_type;
    typedef typename base_type::difference_type difference_type;

    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;
    typedef typename base_type::iterator local_iterator;
    typedef typename base_type::const_iterator const_local_iterator;

    // Construction/destruction/copying
    explicit concurrent_unordered_set(size_type n_of_buckets = 8, const hasher& a_hasher = hasher(),
        const key_equal& a_keyeq = key_equal(), const allocator_type& a = allocator_type())
        : base_type(n_of_buckets, key_compare(a_hasher, a_keyeq), a)
    {
    }

    concurrent_unordered_set(const Allocator& a) : base_type(8, key_compare(), a)
    {
    }

    template <typename Iterator>
    concurrent_unordered_set(Iterator first, Iterator last, size_type n_of_buckets = 8, const hasher& a_hasher = hasher(),
        const key_equal& a_keyeq = key_equal(), const allocator_type& a = allocator_type())
        : base_type(n_of_buckets, key_compare(a_hasher, a_keyeq), a)
    {
        for (; first != last; ++first)
            base_type::insert(*first);
    }

    concurrent_unordered_set(const concurrent_unordered_set& table) : base_type(table)
    {
    }

    concurrent_unordered_set(const concurrent_unordered_set& table, const Allocator& a)
        : base_type(table, a)
    {
    }

    concurrent_unordered_set& operator=(const concurrent_unordered_set& table)
    {
        base_type::operator=(table);
        return (*this);
    }

    iterator unsafe_erase(const_iterator where)
    {
        return base_type::unsafe_erase(where);
    }

    size_type unsafe_erase(const key_type& key)
    {
        return base_type::unsafe_erase(key);
    }

    iterator unsafe_erase(const_iterator first, const_iterator last)
    {
        return base_type::unsafe_erase(first, last);
    }

    void swap(concurrent_unordered_set& table)
    {
        base_type::swap(table);
    }

    // Observers
    hasher hash_function() const
    {
        return my_hash_compare.my_hash_object;
    }

    key_equal key_eq() const
    {
        return my_hash_compare.my_key_compare_object;
    }
};

template <typename Key, typename Hasher = tbb::tbb_hash<Key>, typename Key_equality = std::equal_to<Key>,
         typename Allocator = tbb::tbb_allocator<Key> >
class concurrent_unordered_multiset :
    public internal::concurrent_unordered_base< concurrent_unordered_set_traits<Key,
    internal::hash_compare<Key, Hasher, Key_equality>, Allocator, true> >
{
public:
    // Base type definitions
    typedef internal::hash_compare<Key, Hasher, Key_equality> hash_compare;
    typedef concurrent_unordered_set_traits<Key, hash_compare, Allocator, true> traits_type;
    typedef internal::concurrent_unordered_base< traits_type > base_type;
    using traits_type::allow_multimapping;
    using traits_type::my_hash_compare;

    // Type definitions
    typedef Key key_type;
    typedef typename base_type::value_type value_type;
    typedef Key mapped_type;
    typedef Hasher hasher;
    typedef Key_equality key_equal;
    typedef hash_compare key_compare;

    typedef typename base_type::allocator_type allocator_type;
    typedef typename base_type::pointer pointer;
    typedef typename base_type::const_pointer const_pointer;
    typedef typename base_type::reference reference;
    typedef typename base_type::const_reference const_reference;

    typedef typename base_type::size_type size_type;
    typedef typename base_type::difference_type difference_type;

    typedef typename base_type::iterator iterator;
    typedef typename base_type::const_iterator const_iterator;
    typedef typename base_type::iterator local_iterator;
    typedef typename base_type::const_iterator const_local_iterator;

    // Construction/destruction/copying
    explicit concurrent_unordered_multiset(size_type n_of_buckets = 8,
        const hasher& _Hasher = hasher(), const key_equal& _Key_equality = key_equal(),
        const allocator_type& a = allocator_type())
        : base_type(n_of_buckets, key_compare(_Hasher, _Key_equality), a)
    {
    }

    concurrent_unordered_multiset(const Allocator& a) : base_type(8, key_compare(), a)
    {
    }

    template <typename Iterator>
    concurrent_unordered_multiset(Iterator first, Iterator last, size_type n_of_buckets = 8,
        const hasher& _Hasher = hasher(), const key_equal& _Key_equality = key_equal(),
        const allocator_type& a = allocator_type())
        : base_type(n_of_buckets, key_compare(_Hasher, _Key_equality), a)
    {
        for (; first != last; ++first)
        {
            base_type::insert(*first);
        }
    }

    concurrent_unordered_multiset(const concurrent_unordered_multiset& table) : base_type(table)
    {
    }

    concurrent_unordered_multiset(const concurrent_unordered_multiset& table, const Allocator& a) : base_type(table, a)
    {
    }

    concurrent_unordered_multiset& operator=(const concurrent_unordered_multiset& table)
    {
        base_type::operator=(table);
        return (*this);
    }

    // Modifiers
    std::pair<iterator, bool> insert(const value_type& value)
    {
        return base_type::insert(value);
    }

    iterator insert(const_iterator where, const value_type& value)
    {
        return base_type::insert(where, value);
    }

    template<class Iterator>
    void insert(Iterator first, Iterator last)
    {
        base_type::insert(first, last);
    }

    iterator unsafe_erase(const_iterator where)
    {
        return base_type::unsafe_erase(where);
    }

    size_type unsafe_erase(const key_type& key)
    {
        return base_type::unsafe_erase(key);
    }

    iterator unsafe_erase(const_iterator first, const_iterator last)
    {
        return base_type::unsafe_erase(first, last);
    }

    void swap(concurrent_unordered_multiset& table)
    {
        base_type::swap(table);
    }

    // Observers
    hasher hash_function() const
    {
        return my_hash_compare.my_hash_object;
    }

    key_equal key_eq() const
    {
        return my_hash_compare.my_key_compare_object;
    }
};
} // namespace interface5

using interface5::concurrent_unordered_set;
using interface5::concurrent_unordered_multiset;

} // namespace tbb

#endif// __TBB_concurrent_unordered_set_H
