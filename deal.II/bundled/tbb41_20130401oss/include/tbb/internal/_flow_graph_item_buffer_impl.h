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

#ifndef __TBB__flow_graph_item_buffer_impl_H
#define __TBB__flow_graph_item_buffer_impl_H

#ifndef __TBB_flow_graph_H
#error Do not #include this internal file directly; use public TBB headers instead.
#endif

    //! Expandable buffer of items.  The possible operations are push, pop,
    //* tests for empty and so forth.  No mutual exclusion is built in.
    template <typename T, typename A=cache_aligned_allocator<T> >
    class item_buffer {
    public:
        typedef T input_type;
        typedef T output_type;
    protected:
        typedef size_t size_type;
        typedef std::pair< T, bool > item_type;
        typedef typename A::template rebind<item_type>::other allocator_type;

        item_type *my_array;
        size_type my_array_size;
        static const size_type initial_buffer_size = 4;
        size_type my_head;
        size_type my_tail;

        bool buffer_empty() { return my_head == my_tail; }

        item_type &item(size_type i) { return my_array[i & (my_array_size - 1) ]; } // may not be marked valid

        bool item_valid(size_type i) { return item(i).second; }

        void fetch_front(T &v) { __TBB_ASSERT(item_valid(my_head), "front not valid"); v = item(my_head).first; }
        void fetch_back(T &v) { __TBB_ASSERT(item_valid(my_tail-1), "back not valid"); v = item(my_tail-1).first; }

        void invalidate(size_type i) { __TBB_ASSERT(item_valid(i), "Item not valid"); item(i).second = false; }
        void validate(size_type i) { __TBB_ASSERT(!item_valid(i), "Item already valid"); item(i).second = true; }

        void invalidate_front() { invalidate(my_head); }
        void validate_front() { validate(my_head); }
        void invalidate_back() { invalidate(my_tail-1); }

        size_type size() { return my_tail - my_head; }
        size_type capacity() { return my_array_size; }
        bool buffer_full() { return size() == capacity(); }

        //! Grows the internal array.
        void grow_my_array( size_t minimum_size ) {
            size_type old_size = my_array_size;
            size_type new_size = old_size ? 2*old_size : initial_buffer_size;
            while( new_size<minimum_size )
                new_size*=2;

            item_type* new_array = allocator_type().allocate(new_size);
            item_type* old_array = my_array;

            for( size_type i=0; i<new_size; ++i ) {
                new (&(new_array[i].first)) input_type;
                new_array[i].second = false;
            }

            size_t t=my_head;
            for( size_type i=0; i<old_size; ++i, ++t )
                new_array[t&(new_size-1)] = old_array[t&(old_size-1)];
            my_array = new_array;
            my_array_size = new_size;
            if( old_array ) {
                for( size_type i=0; i<old_size; ++i, ++t )
                    old_array[i].first.~input_type();
                allocator_type().deallocate(old_array,old_size);
            }
        }

        bool push_back(T &v) {
            if(buffer_full()) {
                grow_my_array(size() + 1);
            }
            item(my_tail) = std::make_pair( v, true );
            ++my_tail;
            return true;
        }

        bool pop_back(T &v) {
            if (!item_valid(my_tail-1)) {
                return false;
            }
            fetch_back(v);
            invalidate_back();
            --my_tail;
            return true;
        }

        bool pop_front(T &v) {
            if(!item_valid(my_head)) {
                return false;
            }
            fetch_front(v);
            invalidate_front();
            ++my_head;
            return true;
        }

        void clean_up_buffer() {
            if (my_array) {
                for( size_type i=0; i<my_array_size; ++i ) {
                    my_array[i].first.~input_type();
                }
                allocator_type().deallocate(my_array,my_array_size); 
            }
            my_array = NULL;
            my_head = my_tail = my_array_size = 0;
        }

    public:
        //! Constructor
        item_buffer( ) : my_array(NULL), my_array_size(0),
            my_head(0), my_tail(0) {
            grow_my_array(initial_buffer_size);
        }

        ~item_buffer() {
            clean_up_buffer();
        }

        void reset() { clean_up_buffer(); grow_my_array(initial_buffer_size); }

    };

    //! item_buffer with reservable front-end.  NOTE: if reserving, do not
    //* complete operation with pop_front(); use consume_front().  
    //* No synchronization built-in.
    template<typename T, typename A=cache_aligned_allocator<T> >
    class reservable_item_buffer : public item_buffer<T, A> {
    protected:
        using item_buffer<T, A>::buffer_empty;
        using item_buffer<T, A>::fetch_front;
        using item_buffer<T, A>::invalidate_front;
        using item_buffer<T, A>::validate_front;
        using item_buffer<T, A>::item_valid;
        using item_buffer<T, A>::my_head;

    public:
        reservable_item_buffer() : item_buffer<T, A>(), my_reserved(false) {}
        void reset() {my_reserved = false; item_buffer<T,A>::reset(); }
    protected:

        bool reserve_front(T &v) {
            if(my_reserved || !item_valid(my_head)) return false;
            my_reserved = true;
            // reserving the head
            fetch_front(v);
            // invalidate the head, but don't commit until consume is called
            invalidate_front();
            return true;
        }

        void consume_front() {
            __TBB_ASSERT(my_reserved, "Attempt to consume a non-reserved item");
            ++my_head;
            my_reserved = false;
        }

        void release_front() {
            __TBB_ASSERT(my_reserved, "Attempt to release a non-reserved item");
            validate_front();
            my_reserved = false;
        }

        bool my_reserved;
    };

#endif // __TBB__flow_graph_item_buffer_impl_H
