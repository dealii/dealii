/*
    Copyright 2005-2010 Intel Corporation.  All Rights Reserved.

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

#ifndef _TBB_tbb_statistics_H
#define _TBB_tbb_statistics_H

/**
    This file defines parameters of the internal statistics collected by the TBB
    library (currently by the task scheduler only).
    
    In __TBB_ARENA_PER_MASTER implementation statistics is accumulated in each 
    thread separately and is dumped when the scheduler instance in the given 
    thread is destroyed. For apps with multiple master threads or with the same 
    master repeatedly initializing and then deinitializing task scheduler this 
    results in TBB workers statistics getting unseparably mixed.
    
    Therefore in new __TBB_ARENA_PER_MASTER mode statistics is instead accumulated
    in arena slots, and should be dumped when arena gets destroyed. This separates
    statistics collected for each scheduler activity region in each master thread.

    With the current RML implementation (TBB 2.2, 3.0) to avoid complete loss of 
    statistics data during app shutdown (because of lazy workers deinitialization 
    logic) set __TBB_STATISTICS_EARLY_DUMP macro to write the statistics at the 
    moment a master thread deinitializes its scheduler. This may happen a little 
    earlier than the moment of arena destruction resulting in the following undesired
    (though usually tolerable) effects:
    - a few events related to unsuccessful stealing or thread pool activity may be lost,
    - statistics may be substantially incomplete in case of FIFO tasks used in 
      the FAF mode.

    Macro __TBB_STATISTICS_STDOUT and global variable __TBB_ActiveStatisticsGroups
    defined below can be used to configure the statistics output.

    To add new counter:
    1) Insert it into the appropriate group range in statistics_counters;
    2) Insert the corresponding field title into StatFieldTitles (preserving 
       relative order of the fields).

    To add new counters group:
    1) Insert new group bit flag into statistics_groups;
    2) Insert the new group title into StatGroupTitles (preserving 
       relative order of the groups).
    3) Add counter belonging to the new group as described above
**/

#include "tbb/tbb_stddef.h"

#ifndef __TBB_STATISTICS
#define __TBB_STATISTICS 0
#endif /* __TBB_STATISTICS */

#if __TBB_STATISTICS

#include <string.h>  // for memset

//! Dump counters into stdout as well.
/** By default statistics counters are written to the file "statistics.txt" only. **/
#define __TBB_STATISTICS_STDOUT 1

//! Dump statistics for an arena when its master completes
/** By default (when this macro is not set) the statistics is sent to output when
    arena object is destroyed. But with the current lazy workers termination
    logic default behavior may result in loosing all statistics output. **/
#define __TBB_STATISTICS_EARLY_DUMP 1

#define GATHER_STATISTIC(x) (x)

namespace tbb {
namespace internal {

//! Groups of statistics counters.
/** The order of enumerators must be the same as the order of the corresponding
    field groups in the statistics_counters structure. **/
enum statistics_groups {
    sg_task_allocation = 0x01,
    sg_task_execution = 0x02,
    sg_stealing = 0x04,
    sg_affinity = 0x08,
    sg_arena = 0x10,
    sg_market = 0x20,
    // List end marker. Insert new groups only before it.
    sg_end
};

//! Groups of counters to output
const uintptr_t __TBB_ActiveStatisticsGroups = sg_task_execution | sg_stealing | sg_affinity | sg_arena | sg_market;

//! A set of various statistics counters that are updated by the library on per thread basis.
/** All the fields must be of the same type (statistics_counters::counter_type).
    This is necessary to allow reinterpreting this structure as an array. **/
struct statistics_counters {
    typedef long counter_type;

    // Group: sg_task_allocation
    // Counters in this group can have negative values as the tasks migrate across 
    // threads while the associated counters are updated in the current thread only
    // to avoid data races
    
    //! Number of tasks allocated and not yet destroyed
    counter_type active_tasks;
    //! Number of task corpses stored for future reuse
    counter_type free_list_length;
    //! Number of big tasks allocated during the run
    /** To find total number of tasks malloc'd, compute (big_tasks+small_task_count) */
    counter_type big_tasks;
    
    // Group: sg_task_execution

    //! Number of tasks executed
    counter_type tasks_executed;
    //! Number of elided spawns
    counter_type spawns_bypassed;
    
    // Group: sg_stealing

    //! Number of tasks successfully stolen
    counter_type steals_committed;
    //! Number of failed stealing attempts
    counter_type steals_failed;
    //! Number of failed stealing attempts
    counter_type thieves_conflicts;
    //! Number of tasks received from mailbox

    // Group: sg_affinity

    counter_type mails_received;
    //! Number of affinitized tasks executed by the owner
    /** Goes as "revoked" in statistics printout. **/
    counter_type proxies_executed;
    //! Number of affinitized tasks intercepted by thieves 
    counter_type proxies_stolen;
    //! Number of proxy bypasses by thieves during stealing
    counter_type proxies_bypassed;
    //! Number of affinitized tasks executed by the owner via scheduler bypass mechanism
    counter_type affinity_ignored;

    // Group: sg_arena

    //! Number of times the state of arena switched between "full" and "empty"
    counter_type gate_switches;
    //! Number of times workers left an arena and returned into the market
    counter_type arena_roundtrips;
    //! Number of times workers left the market and returned into RML
    counter_type market_roundtrips;

    // Constructor and helpers

    statistics_counters() { reset(); }

    void reset () { memset( this, 0, sizeof(statistics_counters) ); }

    counter_type& field ( size_t index ) { return reinterpret_cast<counter_type*>(this)[index]; }

    const counter_type& field ( size_t index ) const { return reinterpret_cast<const counter_type*>(this)[index]; }

    static size_t size () { return sizeof(statistics_counters) / sizeof(counter_type); }

    const statistics_counters& operator += ( const statistics_counters& rhs ) {
        for ( size_t i = 0; i < size(); ++i )
            field(i) += rhs.field(i);
        return *this;
    }
}; // statistics_counters

static const size_t workers_counters_total = (size_t)-1;
static const size_t arena_counters_total = (size_t)-2;

void dump_statistics ( const statistics_counters& c, size_t id );

} // namespace internal
} // namespace tbb

#else /* !__TBB_STATISTICS */

#define GATHER_STATISTIC(x) ((void)0)

#endif /* !__TBB_STATISTICS */

#endif /* _TBB_tbb_statistics_H */
