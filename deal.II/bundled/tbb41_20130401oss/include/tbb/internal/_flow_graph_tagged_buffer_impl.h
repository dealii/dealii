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

// tagged buffer that can expand, and can support as many deletions as additions
// list-based, with elements of list held in std::vector (for destruction management),
// multiplicative hashing (like ets).  No synchronization built-in.
//

#ifndef __TBB__flow_graph_tagged_buffer_impl_H
#define __TBB__flow_graph_tagged_buffer_impl_H

#ifndef __TBB_flow_graph_H
#error Do not #include this internal file directly; use public TBB headers instead.
#endif

template<typename TagType, typename ValueType, size_t NoTagMark>
struct buffer_element {
    TagType t;
    ValueType v;
    buffer_element *next;
    buffer_element() : t(NoTagMark), next(NULL) {}
};

template
    <
     typename TagType, 
     typename ValueType, 
     size_t   NoTagMark = 0,
     typename Allocator=tbb::cache_aligned_allocator< buffer_element<TagType,ValueType,NoTagMark> >
    >
class tagged_buffer {
public:
    static const size_t INITIAL_SIZE = 8;  // initial size of the hash pointer table
    static const TagType NO_TAG = TagType(NoTagMark);
    typedef ValueType value_type;
    typedef buffer_element<TagType,ValueType, NO_TAG> element_type;
    typedef value_type *pointer_type;
    typedef std::vector<element_type, Allocator> list_array_type;
    typedef typename Allocator::template rebind<element_type*>::other pointer_array_allocator_type;
    typedef typename Allocator::template rebind<list_array_type>::other list_array_allocator;
private:

    size_t my_size;
    size_t nelements;
    element_type** array;
    std::vector<element_type, Allocator> *lists;
    element_type* free_list;

    size_t mask() { return my_size - 1; }

    static size_t hash(TagType t) {
        return uintptr_t(t)*tbb::internal::select_size_t_constant<0x9E3779B9,0x9E3779B97F4A7C15ULL>::value;
    }

    void set_up_free_list( element_type **p_free_list, list_array_type *la, size_t sz) {
        for(size_t i=0; i < sz - 1; ++i ) {  // construct free list
            (*la)[i].next = &((*la)[i+1]);
            (*la)[i].t = NO_TAG;
        }
        (*la)[sz-1].next = NULL;
        *p_free_list = &((*la)[0]);
    }

    void grow_array() {
        // make the pointer array larger
        element_type **new_array;
        element_type **old_array = array;
        size_t old_size = my_size;
        my_size *=2;
        new_array = pointer_array_allocator_type().allocate(my_size);
        for(size_t i=0; i < my_size; ++i) new_array[i] = NULL;
        list_array_type *new_list_array = new list_array_type(old_size, element_type(), Allocator());
        set_up_free_list(&free_list, new_list_array, old_size );

        for(size_t i=0; i < old_size; ++i) {
            for( element_type* op = old_array[i]; op; op = op->next) {
                internal_tagged_insert(new_array, my_size, op->t, op->v);
            }
        }
        pointer_array_allocator_type().deallocate(old_array, old_size);

        delete lists;  // destroy and deallocate instead
        array = new_array;
        lists = new_list_array;
    }

    void internal_tagged_insert( element_type **ar, size_t sz, TagType t, value_type v) {
        size_t l_mask = sz-1;
        size_t h = hash(t) & l_mask;
        __TBB_ASSERT(free_list, "Error: free list not set up.");
        element_type* my_elem = free_list; free_list = free_list->next;
        my_elem->t = t;
        my_elem->v = v;
        my_elem->next = ar[h];
        ar[h] = my_elem;
    }

    void internal_initialize_buffer() {
        array = pointer_array_allocator_type().allocate(my_size);
        for(size_t i = 0; i < my_size; ++i) array[i] = NULL;
        lists = new list_array_type(INITIAL_SIZE/2, element_type(), Allocator());
        set_up_free_list(&free_list, lists, INITIAL_SIZE/2);
    }

    void internal_free_buffer() {
        if(array) {
            pointer_array_allocator_type().deallocate(array, my_size); 
            array = NULL;
        }
        if(lists) {
            delete lists;
            lists = NULL;
        }
        my_size = INITIAL_SIZE;
        nelements = 0;
    }

public:
    tagged_buffer() : my_size(INITIAL_SIZE), nelements(0) {
        internal_initialize_buffer();
    }

    ~tagged_buffer() {
        internal_free_buffer();
    }

    void reset() {
        internal_free_buffer();
        internal_initialize_buffer();
    }

    bool tagged_insert(TagType t, value_type v) {
        pointer_type p;
        if(tagged_find_ref(t, p)) {
            *p = v;  // replace the value
            return false;
        }
        ++nelements;
        if(nelements*2 > my_size) grow_array();
        internal_tagged_insert(array, my_size, t, v);
        return true;
    }

    // returns reference to array element.v
    bool tagged_find_ref(TagType t, pointer_type &v) {
        size_t i = hash(t) & mask();
        for(element_type* p = array[i]; p; p = p->next) {
            if(p->t == t) {
                v = &(p->v);
                return true;
            }
        }
        return false;
    }

    bool tagged_find( TagType t, value_type &v) {
        value_type *p;
        if(tagged_find_ref(t, p)) {
            v = *p;
            return true;
        }
        else
            return false;
    }

    void tagged_delete(TagType t) {
        size_t h = hash(t) & mask();
        element_type* prev = NULL;
        for(element_type* p = array[h]; p; prev = p, p = p->next) {
            if(p->t == t) {
                p->t = NO_TAG;
                if(prev) prev->next = p->next;
                else array[h] = p->next;
                p->next = free_list;
                free_list = p;
                --nelements;
                return;
            }
        }
        __TBB_ASSERT(false, "tag not found for delete");
    }

    // search for v in the array; if found {set t, return true} else return false
    // we use this in join_node_FE to find if a tag's items are all available.
    bool find_value_tag( TagType &t, value_type v) {
        for(size_t i= 0; i < my_size / 2; ++i) {  // remember the vector is half the size of the hash array
            if( (*lists)[i].t != NO_TAG && (*lists)[i].v == v) {
                t = (*lists)[i].t;
                return true;
            }
        }
        return false;
    }
};
#endif // __TBB__flow_graph_tagged_buffer_impl_H
