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

#ifndef _TBB_mailbox_H
#define _TBB_mailbox_H

#include "tbb/tbb_stddef.h"
#include "tbb/cache_aligned_allocator.h"

#include "scheduler_common.h"

namespace tbb {
namespace internal {

class mail_outbox;

struct task_proxy : public task {
    static const intptr_t pool_bit = 1;
    static const intptr_t mailbox_bit = 2;
    /* All but two low-order bits represent a (task*).
       Two low-order bits mean:
       1 = proxy is/was/will be in task pool
       2 = proxy is/was/will be in mailbox */
    intptr_t task_and_tag;

    //! Pointer to next task_proxy in a mailbox
    task_proxy* next_in_mailbox;

    //! Mailbox to which this was mailed.
    mail_outbox* outbox;
};

//! Internal representation of mail_outbox, without padding.
class unpadded_mail_outbox {
protected:
    //! Pointer to first task_proxy in mailbox, or NULL if box is empty. 
    task_proxy* my_first;

    //! Pointer to pointer that will point to next item in the queue.  Never NULL.
    task_proxy** my_last;

    //! Owner of mailbox is not executing a task, and has drained its own task pool.
    bool my_is_idle;
};

//! Class representing where mail is put.
/** Padded to occupy a cache line. */
class mail_outbox: unpadded_mail_outbox {
    char pad[NFS_MaxLineSize-sizeof(unpadded_mail_outbox)];

    task_proxy* internal_pop() {
        //! No fence on load of my_first, because if it is NULL, there's nothing further to read from another thread.
        task_proxy* first = my_first;
        if( first ) {
            // There is a first item in the mailbox.  See if there is a second.
            if( task_proxy* second = __TBB_load_with_acquire(first->next_in_mailbox) ) {
                // There are at least two items, so first item can be popped easily.
                __TBB_store_with_release( my_first, second );
            } else {
                // There is only one item.  Some care is required to pop it.
                my_first = NULL;
                if( (task_proxy**)__TBB_CompareAndSwapW(&my_last, (intptr_t)&my_first,
                    (intptr_t)&first->next_in_mailbox)==&first->next_in_mailbox ) 
                {
                    // Successfully transitioned mailbox from having one item to having none.
                    __TBB_ASSERT(!first->next_in_mailbox,NULL);
                } else {
                    // Some other thread updated my_last but has not filled in result->next_in_mailbox
                    // Wait until first item points to second item.
                    atomic_backoff backoff;
                    while( !(second=const_cast<volatile task_proxy*>(first)->next_in_mailbox) ) 
                        backoff.pause();
                    my_first = second;
                } 
            }
        }
        return first;
    }
public:
    friend class mail_inbox;

    //! Push task_proxy onto the mailbox queue of another thread.
    /** Implementation is wait-free. */
    void push( task_proxy& t ) {
        __TBB_ASSERT(&t, NULL);
        t.next_in_mailbox = NULL; 
        task_proxy** link = (task_proxy**)__TBB_FetchAndStoreW(&my_last,(intptr_t)&t.next_in_mailbox);
        // No release fence required for the next store, because there are no memory operations 
        // between the previous fully fenced atomic operation and the store.
        *link = &t;
    }

    //! Construct *this as a mailbox from zeroed memory.
    /** Raise assertion if *this is not previously zeored, or sizeof(this) is wrong.
        This method is provided instead of a full constructor since we know the objecxt
        will be constructed in zeroed memory. */
    void construct() {
        __TBB_ASSERT( sizeof(*this)==NFS_MaxLineSize, NULL );
        __TBB_ASSERT( !my_first, NULL );
        __TBB_ASSERT( !my_last, NULL );
        __TBB_ASSERT( !my_is_idle, NULL );
        my_last=&my_first;
    }

    //! Drain the mailbox 
    intptr_t drain() {
        intptr_t k = 0;
        // No fences here because other threads have already quit.
        for( ; task_proxy* t = my_first; ++k ) {
            my_first = t->next_in_mailbox;
            NFS_Free((char*)t - task_prefix_reservation_size);
        }
        return k;  
    }

    //! True if thread that owns this mailbox is looking for work.
    bool recipient_is_idle() {
        return my_is_idle;
    }
}; // class mail_outbox

//! Class representing source of mail.
class mail_inbox {
    //! Corresponding sink where mail that we receive will be put.
    mail_outbox* my_putter;
public:
    //! Construct unattached inbox
    mail_inbox() : my_putter(NULL) {}

    //! Attach inbox to a corresponding outbox. 
    void attach( mail_outbox& putter ) {
        __TBB_ASSERT(!my_putter,"already attached");
        my_putter = &putter;
    }
    //! Detach inbox from its outbox
    void detach() {
        __TBB_ASSERT(my_putter,"not attached");
        my_putter = NULL;
    }
    //! Get next piece of mail, or NULL if mailbox is empty.
    task_proxy* pop() {
        return my_putter->internal_pop();
    }
    //! Indicate whether thread that reads this mailbox is idle.
    /** Raises assertion failure if mailbox is redundantly marked as not idle. */
    void set_is_idle( bool value ) {
        if( my_putter ) {
            __TBB_ASSERT( my_putter->my_is_idle || value, "attempt to redundantly mark mailbox as not idle" );
            my_putter->my_is_idle = value;
        }
    }
    //! Indicate whether thread that reads this mailbox is idle.
    bool is_idle_state ( bool value ) const {
        return !my_putter || my_putter->my_is_idle == value;
    }

#if DO_ITT_NOTIFY
    //! Get pointer to corresponding outbox used for ITT_NOTIFY calls.
    void* outbox() const {return my_putter;}
#endif /* DO_ITT_NOTIFY */ 
}; // class mail_inbox

} // namespace internal
} // namespace tbb

#endif /* _TBB_mailbox_H */
