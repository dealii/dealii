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

#ifndef __TBB_partitioner_H
#define __TBB_partitioner_H

#ifndef __TBB_INITIAL_CHUNKS
#define __TBB_INITIAL_CHUNKS 2
#endif
#ifndef __TBB_RANGE_POOL_CAPACITY
#define __TBB_RANGE_POOL_CAPACITY 8
#endif
#ifndef __TBB_INIT_DEPTH
#define __TBB_INIT_DEPTH 5
#endif

#include "task.h"
#include "aligned_space.h"
#include "atomic.h"

#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    // Workaround for overzealous compiler warnings
    #pragma warning (push)
    #pragma warning (disable: 4244)
#endif

namespace tbb {

class auto_partitioner;
class simple_partitioner;
class affinity_partitioner;
namespace interface6 {
    namespace internal {
        class affinity_partition_type;
    }
}

namespace internal {
size_t __TBB_EXPORTED_FUNC get_initial_auto_partitioner_divisor();

//! Defines entry point for affinity partitioner into tbb run-time library.
class affinity_partitioner_base_v3: no_copy {
    friend class tbb::affinity_partitioner;
    friend class tbb::interface6::internal::affinity_partition_type;
    //! Array that remembers affinities of tree positions to affinity_id.
    /** NULL if my_size==0. */
    affinity_id* my_array;
    //! Number of elements in my_array.
    size_t my_size;
    //! Zeros the fields.
    affinity_partitioner_base_v3() : my_array(NULL), my_size(0) {}
    //! Deallocates my_array.
    ~affinity_partitioner_base_v3() {resize(0);}
    //! Resize my_array.
    /** Retains values if resulting size is the same. */
    void __TBB_EXPORTED_METHOD resize( unsigned factor );
};

//! Provides backward-compatible methods for partition objects without affinity.
class partition_type_base {
public:
    void set_affinity( task & ) {}
    void note_affinity( task::affinity_id ) {}
    task* continue_after_execute_range() {return NULL;}
    bool decide_whether_to_delay() {return false;}
    void spawn_or_delay( bool, task& b ) {
        task::spawn(b);
    }
};

template<typename Range, typename Body, typename Partitioner> class start_scan;

} // namespace internal
//! @endcond

namespace serial {
namespace interface6 {
template<typename Range, typename Body, typename Partitioner> class start_for;
}
}

namespace interface6 {
//! @cond INTERNAL
namespace internal {
using namespace tbb::internal;
template<typename Range, typename Body, typename Partitioner> class start_for;
template<typename Range, typename Body, typename Partitioner> class start_reduce;

//! Join task node that contains shared flag for stealing feedback
class flag_task: public task {
public:
    tbb::atomic<bool> my_child_stolen;
    flag_task() { my_child_stolen = false; }
    task* execute() { return NULL; }
    static void mark_task_stolen(task &t) {
        tbb::atomic<bool> &flag = static_cast<flag_task*>(t.parent())->my_child_stolen;
#if TBB_USE_THREADING_TOOLS
        // Threading tools respect lock prefix but report false-positive data-race via plain store
        flag.fetch_and_store<release>(true);
#else
        flag = true;
#endif //TBB_USE_THREADING_TOOLS
    }
    static bool is_peer_stolen(task &t) {
        return static_cast<flag_task*>(t.parent())->my_child_stolen;
    }
};

//! Task to signal the demand without carrying the work
class signal_task: public task {
public:
    task* execute() {
        if( is_stolen_task() ) {
            flag_task::mark_task_stolen(*this);
        }
        return NULL;
    }
};

//! Depth is a relative depth of recursive division inside a range pool. Relative depth allows
//! infinite absolute depth of the recursion for heavily imbalanced workloads with range represented
//! by a number that cannot fit into machine word.
typedef unsigned char depth_t;

//! Range pool stores ranges of type T in a circular buffer with MaxCapacity
template <typename T, depth_t MaxCapacity>
class range_vector {
    depth_t my_head;
    depth_t my_tail;
    depth_t my_size;
    depth_t my_depth[MaxCapacity]; // relative depths of stored ranges
    tbb::aligned_space<T, MaxCapacity> my_pool;

public:
    //! initialize via first range in pool
    range_vector(const T& elem) : my_head(0), my_tail(0), my_size(1) {
        my_depth[0] = 0;
        new( my_pool.begin() ) T(elem);//TODO: std::move?
    }
    ~range_vector() {
        while( !empty() ) pop_back();
    }
    bool empty() const { return my_size == 0; }
    depth_t size() const { return my_size; }
    //! Populates range pool via ranges up to max depth or while divisible
    //! max_depth starts from 0, e.g. value 2 makes 3 ranges in the pool up to two 1/4 pieces
    void split_to_fill(depth_t max_depth) {
        while( my_size < MaxCapacity && my_depth[my_head] < max_depth
          && my_pool.begin()[my_head].is_divisible() ) {
            depth_t prev = my_head;
            my_head = (my_head + 1) % MaxCapacity;
            new(my_pool.begin()+my_head) T(my_pool.begin()[prev]); // copy TODO: std::move?
            my_pool.begin()[prev].~T(); // instead of assignment
            new(my_pool.begin()+prev) T(my_pool.begin()[my_head], split()); // do 'inverse' split
            my_depth[my_head] = ++my_depth[prev];
            my_size++;
        }
    }
    void pop_back() {
        __TBB_ASSERT(my_size > 0, "range_vector::pop_back() with empty size");
        my_pool.begin()[my_head].~T();
        my_size--;
        my_head = (my_head + MaxCapacity - 1) % MaxCapacity;
    }
    void pop_front() {
        __TBB_ASSERT(my_size > 0, "range_vector::pop_front() with empty size");
        my_pool.begin()[my_tail].~T();
        my_size--;
        my_tail = (my_tail + 1) % MaxCapacity;
    }
    T& back() {
        __TBB_ASSERT(my_size > 0, "range_vector::back() with empty size");
        return my_pool.begin()[my_head];
    }
    T& front() {
        __TBB_ASSERT(my_size > 0, "range_vector::front() with empty size");
        return my_pool.begin()[my_tail];
    }
    //! similarly to front(), returns depth of the first range in the pool
    depth_t front_depth() {
        __TBB_ASSERT(my_size > 0, "range_vector::front_depth() with empty size");
        return my_depth[my_tail];
    }
};

//! Provides default methods for partition objects and common algorithm blocks.
template <typename Partition>
struct partition_type_base {
    // decision makers
    void set_affinity( task & ) {}
    void note_affinity( task::affinity_id ) {}
    bool check_being_stolen(task &) { return false; } // part of old should_execute_range()
    bool check_for_demand(task &) { return false; }
    bool divisions_left() { return true; } // part of old should_execute_range()
    bool should_create_trap() { return false; }
    depth_t max_depth() { return 0; }
    void align_depth(depth_t) { }
    // common function blocks
    Partition& derived() { return *static_cast<Partition*>(this); }
    template<typename StartType>
    flag_task* split_work(StartType &start) {
        flag_task* parent_ptr = start.create_continuation(); // the type here is to express expectation
        start.set_parent(parent_ptr);
        parent_ptr->set_ref_count(2);
        StartType& right_work = *new( parent_ptr->allocate_child() ) StartType(start, split());
        start.spawn(right_work);
        return parent_ptr;
    }
    template<typename StartType, typename Range>
    void execute(StartType &start, Range &range) {
        // The algorithm in a few words ([]-denotes calls to decision methods of partitioner):
        // [If this task is stolen, adjust depth and divisions if necessary, set flag].
        // If range is divisible {
        //    Spread the work while [initial divisions left];
        //    Create trap task [if necessary];
        // }
        // If not divisible or [max depth is reached], execute, else do the range pool part
        task* parent_ptr = start.parent();
        if( range.is_divisible() ) {
            if( derived().divisions_left() )
                do parent_ptr = split_work(start); // split until divisions_left()
                while( range.is_divisible() && derived().divisions_left() );
            if( derived().should_create_trap() ) { // only for range pool
                if( parent_ptr->ref_count() > 1 ) { // create new parent if necessary
                    parent_ptr = start.create_continuation();
                    start.set_parent(parent_ptr);
                } else __TBB_ASSERT(parent_ptr->ref_count() == 1, NULL);
                parent_ptr->set_ref_count(2); // safe because parent has only one reference
                signal_task& right_signal = *new( parent_ptr->allocate_child() ) signal_task();
                start.spawn(right_signal); // pure signal is to avoid deep recursion in the end
            }
        }
        if( !range.is_divisible() || !derived().max_depth() )
            start.run_body( range ); // simple partitioner goes always here
        else { // do range pool
            internal::range_vector<Range, Partition::range_pool_size> range_pool(range);
            do {
                range_pool.split_to_fill(derived().max_depth()); // fill range pool
                if( derived().check_for_demand( start ) ) {
                    if( range_pool.size() > 1 ) {
                        parent_ptr = start.create_continuation();
                        start.set_parent(parent_ptr);
                        parent_ptr->set_ref_count(2);
                        StartType& right_work = *new( parent_ptr->allocate_child() ) StartType(start, range_pool.front(), range_pool.front_depth());
                        start.spawn(right_work);
                        range_pool.pop_front();
                        continue;
                    }
                    if( range_pool.back().is_divisible() ) // was not enough depth to fork a task
                        continue; // note: check_for_demand() should guarantee increasing max_depth() next time
                }
                start.run_body( range_pool.back() );
                range_pool.pop_back();
            } while( !range_pool.empty() && !start.is_cancelled() );
        }
    }
};

//! Provides default methods for auto (adaptive) partition objects.
template <typename Partition>
struct auto_partition_type_base : partition_type_base<Partition> {
    size_t my_divisor;
    depth_t my_max_depth;
    auto_partition_type_base() : my_max_depth(__TBB_INIT_DEPTH) {
        my_divisor = tbb::internal::get_initial_auto_partitioner_divisor()*__TBB_INITIAL_CHUNKS/4;
        __TBB_ASSERT(my_divisor, "initial value of get_initial_auto_partitioner_divisor() is not valid");
    }
    auto_partition_type_base(auto_partition_type_base &src, split) {
        my_max_depth = src.my_max_depth;
#if __TBB_INITIAL_TASK_IMBALANCE
        if( src.my_divisor <= 1 ) my_divisor = 0;
        else my_divisor = src.my_divisor = (src.my_divisor+1u) / 2u;
#else
        my_divisor = src.my_divisor / 2u;
        src.my_divisor = src.my_divisor - my_divisor; // TODO: check the effect separately
        if(my_divisor) src.my_max_depth += static_cast<depth_t>(__TBB_Log2(src.my_divisor/my_divisor));
#endif
    }
    bool check_being_stolen( task &t) { // part of old should_execute_range()
        if( !my_divisor ) { // if not from the top P tasks of binary tree
            my_divisor = 1; // TODO: replace by on-stack flag (partition_state's member)?
            if( t.is_stolen_task() ) {
#if TBB_USE_EXCEPTIONS
                // RTTI is available, check whether the cast is valid
                __TBB_ASSERT(dynamic_cast<flag_task*>(t.parent()), 0);
                // correctness of the cast relies on avoiding the root task for which:
                // - initial value of my_divisor != 0 (protected by separate assertion)
                // - is_stolen_task() always returns false for the root task.
#endif
                flag_task::mark_task_stolen(t);
                my_max_depth++;
                return true;
            }
        }
        return false;
    }
    bool divisions_left() { // part of old should_execute_range()
        if( my_divisor > 1 ) return true;
        if( my_divisor && my_max_depth > 1 ) { // can split the task and once more internally. TODO: on-stack flag instead
            // keep same fragmentation while splitting for the local task pool
            my_max_depth--;
            my_divisor = 0; // decrease max_depth once per task
            return true;
        } else return false;
    }
    bool should_create_trap() {
        return my_divisor > 0;
    }
    bool check_for_demand(task &t) {
        if( flag_task::is_peer_stolen(t) ) {
            my_max_depth++;
            return true;
        } else return false;
    }
    void align_depth(depth_t base) {
        __TBB_ASSERT(base <= my_max_depth, 0);
        my_max_depth -= base;
    }
    depth_t max_depth() { return my_max_depth; }
};

//! Provides default methods for affinity (adaptive) partition objects.
class affinity_partition_type : public auto_partition_type_base<affinity_partition_type> {
    static const unsigned factor_power = 4;
    static const unsigned factor = 1<<factor_power;
    bool my_delay;
    unsigned map_begin, map_end, map_mid;
    tbb::internal::affinity_id* my_array;
    void set_mid() {
        unsigned d = (map_end - map_begin)/2; // we could add 1 but it is rather for LIFO affinity
        if( d > factor )
            d &= 0u-factor;
        map_mid = map_end - d;
    }
public:
    affinity_partition_type( tbb::internal::affinity_partitioner_base_v3& ap ) {
        __TBB_ASSERT( (factor&(factor-1))==0, "factor must be power of two" );
        ap.resize(factor);
        my_array = ap.my_array;
        map_begin = 0;
        map_end = unsigned(ap.my_size);
        set_mid();
        my_delay = true;
        my_divisor /= __TBB_INITIAL_CHUNKS; // let exactly P tasks to be distributed across workers
        my_max_depth = factor_power+1; // the first factor_power ranges will be spawned, and >=1 ranges should be left
        __TBB_ASSERT( my_max_depth < __TBB_RANGE_POOL_CAPACITY, 0 );
    }
    affinity_partition_type(affinity_partition_type& p, split)
        : auto_partition_type_base<affinity_partition_type>(p, split()), my_array(p.my_array) {
        __TBB_ASSERT( p.map_end-p.map_begin<factor || (p.map_end-p.map_begin)%factor==0, NULL );
        map_end = p.map_end;
        map_begin = p.map_end = p.map_mid;
        set_mid(); p.set_mid();
        my_delay = p.my_delay;
    }
    void set_affinity( task &t ) {
        if( map_begin<map_end )
            t.set_affinity( my_array[map_begin] );
    }
    void note_affinity( task::affinity_id id ) {
        if( map_begin<map_end )
            my_array[map_begin] = id;
    }
    bool check_for_demand( task &t ) {
        if( !my_delay ) {
            if( map_mid<map_end ) {
                __TBB_ASSERT(my_max_depth>__TBB_Log2(map_end-map_mid), 0);
                return true;// do not do my_max_depth++ here, but be sure my_max_depth is big enough
            }
            if( flag_task::is_peer_stolen(t) ) {
                my_max_depth++;
                return true;
            }
        } else my_delay = false;
        return false;
    }
    bool divisions_left() { // part of old should_execute_range()
        return my_divisor > 1;
    }
    bool should_create_trap() {
        return true; // TODO: rethink for the stage after memorizing level
    }
    static const unsigned range_pool_size = __TBB_RANGE_POOL_CAPACITY;
};

class auto_partition_type: public auto_partition_type_base<auto_partition_type> {
public:
    auto_partition_type( const auto_partitioner& ) {}
    auto_partition_type( auto_partition_type& src, split)
      : auto_partition_type_base<auto_partition_type>(src, split()) {}
    static const unsigned range_pool_size = __TBB_RANGE_POOL_CAPACITY;
};

class simple_partition_type: public partition_type_base<simple_partition_type> {
public:
    simple_partition_type( const simple_partitioner& ) {}
    simple_partition_type( const simple_partition_type&, split ) {}
    //! simplified algorithm
    template<typename StartType, typename Range>
    void execute(StartType &start, Range &range) {
        while( range.is_divisible() )
            split_work( start );
        start.run_body( range );
    }
    //static const unsigned range_pool_size = 1; - not necessary because execute() is overridden
};

//! Backward-compatible partition for auto and affinity partition objects.
class old_auto_partition_type: public tbb::internal::partition_type_base {
    size_t num_chunks;
    static const size_t VICTIM_CHUNKS = 4;
public:
    bool should_execute_range(const task &t) {
        if( num_chunks<VICTIM_CHUNKS && t.is_stolen_task() )
            num_chunks = VICTIM_CHUNKS;
        return num_chunks==1;
    }
    old_auto_partition_type( const auto_partitioner& )
      : num_chunks(internal::get_initial_auto_partitioner_divisor()*__TBB_INITIAL_CHUNKS/4) {}
    old_auto_partition_type( const affinity_partitioner& )
      : num_chunks(internal::get_initial_auto_partitioner_divisor()*__TBB_INITIAL_CHUNKS/4) {}
    old_auto_partition_type( old_auto_partition_type& pt, split ) {
        num_chunks = pt.num_chunks = (pt.num_chunks+1u) / 2u;
    }
};

} // namespace interfaceX::internal
//! @endcond
} // namespace interfaceX

//! A simple partitioner
/** Divides the range until the range is not divisible.
    @ingroup algorithms */
class simple_partitioner {
public:
    simple_partitioner() {}
private:
    template<typename Range, typename Body, typename Partitioner> friend class serial::interface6::start_for;
    template<typename Range, typename Body, typename Partitioner> friend class interface6::internal::start_for;
    template<typename Range, typename Body, typename Partitioner> friend class interface6::internal::start_reduce;
    template<typename Range, typename Body, typename Partitioner> friend class internal::start_scan;
    // backward compatibility
    class partition_type: public internal::partition_type_base {
    public:
        bool should_execute_range(const task& ) {return false;}
        partition_type( const simple_partitioner& ) {}
        partition_type( const partition_type&, split ) {}
    };
    // new implementation just extends existing interface
    typedef interface6::internal::simple_partition_type task_partition_type;
};

//! An auto partitioner
/** The range is initial divided into several large chunks.
    Chunks are further subdivided into smaller pieces if demand detected and they are divisible.
    @ingroup algorithms */
class auto_partitioner {
public:
    auto_partitioner() {}

private:
    template<typename Range, typename Body, typename Partitioner> friend class serial::interface6::start_for;
    template<typename Range, typename Body, typename Partitioner> friend class interface6::internal::start_for;
    template<typename Range, typename Body, typename Partitioner> friend class interface6::internal::start_reduce;
    template<typename Range, typename Body, typename Partitioner> friend class internal::start_scan;
    // backward compatibility
    typedef interface6::internal::old_auto_partition_type partition_type;
    // new implementation just extends existing interface
    typedef interface6::internal::auto_partition_type task_partition_type;
};

//! An affinity partitioner
class affinity_partitioner: internal::affinity_partitioner_base_v3 {
public:
    affinity_partitioner() {}

private:
    template<typename Range, typename Body, typename Partitioner> friend class serial::interface6::start_for;
    template<typename Range, typename Body, typename Partitioner> friend class interface6::internal::start_for;
    template<typename Range, typename Body, typename Partitioner> friend class interface6::internal::start_reduce;
    template<typename Range, typename Body, typename Partitioner> friend class internal::start_scan;
    // backward compatibility - for parallel_scan only
    typedef interface6::internal::old_auto_partition_type partition_type;
    // new implementation just extends existing interface
    typedef interface6::internal::affinity_partition_type task_partition_type;
};

} // namespace tbb

#if defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    #pragma warning (pop)
#endif // warning 4244 is back
#undef __TBB_INITIAL_CHUNKS
#undef __TBB_RANGE_POOL_CAPACITY
#undef __TBB_INIT_DEPTH
#endif /* __TBB_partitioner_H */
