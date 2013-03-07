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

#include "ittnotify_config.h"

#ifndef ITT_STUB
#  define ITT_STUB ITT_STUBV
#endif /* ITT_STUB */

#ifndef ITTAPI
#  define ITTAPI CDECL
#endif /* ITTAPI */

#ifndef LIBITTAPI
#  define LIBITTAPI /* nothing */
#endif /* LIBITTAPI */

#ifndef ITT_FORMAT_DEFINED
#  ifndef ITT_FORMAT
#    define ITT_FORMAT
#  endif /* ITT_FORMAT */
#  ifndef ITT_NO_PARAMS
#    define ITT_NO_PARAMS
#  endif /* ITT_NO_PARAMS */
#endif /* ITT_FORMAT_DEFINED */

/*
 * parameters for macro expected:
 * ITT_STUB(api, type, func_name, arguments, params, func_name_in_dll, group, printf_fmt)
 */
/* public */
ITT_STUBV(ITTAPI, void, pause,  (void), (ITT_NO_PARAMS), pause,  __itt_group_control | __itt_group_legacy, "no args")
ITT_STUBV(ITTAPI, void, resume, (void), (ITT_NO_PARAMS), resume, __itt_group_control | __itt_group_legacy, "no args")

#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUBV(ITTAPI, void, thread_set_nameA, (const char    *name), (ITT_FORMAT name), thread_set_nameA, __itt_group_thread, "\"%s\"")
ITT_STUBV(ITTAPI, void, thread_set_nameW, (const wchar_t *name), (ITT_FORMAT name), thread_set_nameW, __itt_group_thread, "\"%S\"")
#else  /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
ITT_STUBV(ITTAPI, void, thread_set_name,  (const char    *name), (ITT_FORMAT name), thread_set_name,  __itt_group_thread, "\"%s\"")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
ITT_STUBV(ITTAPI, void, thread_ignore, (void), (ITT_NO_PARAMS), thread_ignore, __itt_group_thread, "no args")

#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUBV(ITTAPI, void, sync_createA, (void *addr, const char    *objtype, const char    *objname, int attribute), (ITT_FORMAT addr, objtype, objname, attribute), sync_createA, __itt_group_sync | __itt_group_fsync, "%p, \"%s\", \"%s\", %x")
ITT_STUBV(ITTAPI, void, sync_createW, (void *addr, const wchar_t *objtype, const wchar_t *objname, int attribute), (ITT_FORMAT addr, objtype, objname, attribute), sync_createW, __itt_group_sync | __itt_group_fsync, "%p, \"%S\", \"%S\", %x")
ITT_STUBV(ITTAPI, void, sync_renameA, (void *addr, const char    *name), (ITT_FORMAT addr, name), sync_renameA, __itt_group_sync | __itt_group_fsync, "%p, \"%s\"")
ITT_STUBV(ITTAPI, void, sync_renameW, (void *addr, const wchar_t *name), (ITT_FORMAT addr, name), sync_renameW, __itt_group_sync | __itt_group_fsync, "%p, \"%S\"")
#else  /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
ITT_STUBV(ITTAPI, void, sync_create,  (void *addr, const char    *objtype, const char    *objname, int attribute), (ITT_FORMAT addr, objtype, objname, attribute), sync_create,  __itt_group_sync | __itt_group_fsync, "%p, \"%s\", \"%s\", %x")
ITT_STUBV(ITTAPI, void, sync_rename,  (void *addr, const char    *name), (ITT_FORMAT addr, name), sync_rename,  __itt_group_sync | __itt_group_fsync, "%p, \"%s\"")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
ITT_STUBV(ITTAPI, void, sync_destroy,    (void *addr), (ITT_FORMAT addr), sync_destroy,   __itt_group_sync | __itt_group_fsync, "%p")

ITT_STUBV(ITTAPI, void, sync_prepare,    (void* addr), (ITT_FORMAT addr), sync_prepare,   __itt_group_sync,  "%p")
ITT_STUBV(ITTAPI, void, sync_cancel,     (void *addr), (ITT_FORMAT addr), sync_cancel,    __itt_group_sync,  "%p")
ITT_STUBV(ITTAPI, void, sync_acquired,   (void *addr), (ITT_FORMAT addr), sync_acquired,  __itt_group_sync,  "%p")
ITT_STUBV(ITTAPI, void, sync_releasing,  (void* addr), (ITT_FORMAT addr), sync_releasing, __itt_group_sync,  "%p")

ITT_STUBV(ITTAPI, void, fsync_prepare,   (void* addr), (ITT_FORMAT addr), sync_prepare,   __itt_group_fsync, "%p")
ITT_STUBV(ITTAPI, void, fsync_cancel,    (void *addr), (ITT_FORMAT addr), sync_cancel,    __itt_group_fsync, "%p")
ITT_STUBV(ITTAPI, void, fsync_acquired,  (void *addr), (ITT_FORMAT addr), sync_acquired,  __itt_group_fsync, "%p")
ITT_STUBV(ITTAPI, void, fsync_releasing, (void* addr), (ITT_FORMAT addr), sync_releasing, __itt_group_fsync, "%p")

ITT_STUBV(ITTAPI, void, model_site_begin,          (__itt_model_site *site, __itt_model_site_instance *instance, const char *name), (ITT_FORMAT site, instance, name), model_site_begin, __itt_group_model, "%p, %p, \"%s\"")
ITT_STUBV(ITTAPI, void, model_site_end,            (__itt_model_site *site, __itt_model_site_instance *instance),                   (ITT_FORMAT site, instance),       model_site_end,   __itt_group_model, "%p, %p")
ITT_STUBV(ITTAPI, void, model_task_begin,          (__itt_model_task *task, __itt_model_task_instance *instance, const char *name), (ITT_FORMAT task, instance, name), model_task_begin, __itt_group_model, "%p, %p, \"%s\"")
ITT_STUBV(ITTAPI, void, model_task_end,            (__itt_model_task *task, __itt_model_task_instance *instance),                   (ITT_FORMAT task, instance),       model_task_end,   __itt_group_model, "%p, %p")
ITT_STUBV(ITTAPI, void, model_lock_acquire,        (void *lock), (ITT_FORMAT lock), model_lock_acquire, __itt_group_model, "%p")
ITT_STUBV(ITTAPI, void, model_lock_release,        (void *lock), (ITT_FORMAT lock), model_lock_release, __itt_group_model, "%p")
ITT_STUBV(ITTAPI, void, model_record_allocation,   (void *addr, size_t size), (ITT_FORMAT addr, size), model_record_allocation,   __itt_group_model, "%p, %d")
ITT_STUBV(ITTAPI, void, model_record_deallocation, (void *addr),              (ITT_FORMAT addr),       model_record_deallocation, __itt_group_model, "%p")
ITT_STUBV(ITTAPI, void, model_induction_uses,      (void* addr, size_t size), (ITT_FORMAT addr, size), model_induction_uses,      __itt_group_model, "%p, %d")
ITT_STUBV(ITTAPI, void, model_reduction_uses,      (void* addr, size_t size), (ITT_FORMAT addr, size), model_reduction_uses,      __itt_group_model, "%p, %d")
ITT_STUBV(ITTAPI, void, model_observe_uses,        (void* addr, size_t size), (ITT_FORMAT addr, size), model_observe_uses,        __itt_group_model, "%p, %d")
ITT_STUBV(ITTAPI, void, model_clear_uses,          (void* addr),              (ITT_FORMAT addr),       model_clear_uses,          __itt_group_model, "%p")
ITT_STUBV(ITTAPI, void, model_disable_push,        (__itt_model_disable x),   (ITT_FORMAT x),          model_disable_push,        __itt_group_model, "%p")
ITT_STUBV(ITTAPI, void, model_disable_pop,         (void),                    (ITT_NO_PARAMS),         model_disable_pop,         __itt_group_model, "no args")

#ifndef __ITT_INTERNAL_BODY
#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUB(ITTAPI, __itt_counter, counter_createA, (const char    *name, const char    *domain), (ITT_FORMAT name, domain), counter_createA, __itt_group_counter, "\"%s\", \"%s\"")
ITT_STUB(ITTAPI, __itt_counter, counter_createW, (const wchar_t *name, const wchar_t *domain), (ITT_FORMAT name, domain), counter_createW, __itt_group_counter, "\"%s\", \"%s\"")
#else  /* ITT_PLATFORM==ITT_PLATFORM_WIN */
ITT_STUB(ITTAPI, __itt_counter, counter_create,  (const char    *name, const char    *domain), (ITT_FORMAT name, domain), counter_create,  __itt_group_counter, "\"%s\", \"%s\"")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#endif /* __ITT_INTERNAL_BODY */
ITT_STUBV(ITTAPI, void, counter_destroy,   (__itt_counter id),                           (ITT_FORMAT id),        counter_destroy,   __itt_group_counter, "%p")
ITT_STUBV(ITTAPI, void, counter_inc,       (__itt_counter id),                           (ITT_FORMAT id),        counter_inc,       __itt_group_counter, "%p")
ITT_STUBV(ITTAPI, void, counter_inc_delta, (__itt_counter id, unsigned long long value), (ITT_FORMAT id, value), counter_inc_delta, __itt_group_counter, "%p, %lu")

#ifndef __ITT_INTERNAL_BODY
ITT_STUB(ITTAPI, __itt_caller, stack_caller_create,  (void),         (ITT_NO_PARAMS), stack_caller_create,  __itt_group_stitch, "no args")
#endif /* __ITT_INTERNAL_BODY */
ITT_STUBV(ITTAPI, void, stack_caller_destroy,     (__itt_caller id), (ITT_FORMAT id), stack_caller_destroy, __itt_group_stitch, "%p")
ITT_STUBV(ITTAPI, void, stack_callee_enter,       (__itt_caller id), (ITT_FORMAT id), stack_callee_enter,   __itt_group_stitch, "%p")
ITT_STUBV(ITTAPI, void, stack_callee_leave,       (__itt_caller id), (ITT_FORMAT id), stack_callee_leave,   __itt_group_stitch, "%p")

#ifndef __ITT_INTERNAL_BODY
#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUB(ITTAPI, __itt_frame, frame_createA, (const char    *domain), (ITT_FORMAT domain), frame_createA, __itt_group_frame, "\"%s\"")
ITT_STUB(ITTAPI, __itt_frame, frame_createW, (const wchar_t *domain), (ITT_FORMAT domain), frame_createW, __itt_group_frame, "\"%s\"")
#else  /* ITT_PLATFORM==ITT_PLATFORM_WIN */
ITT_STUB(ITTAPI, __itt_frame, frame_create,  (const char    *domain), (ITT_FORMAT domain), frame_create,  __itt_group_frame, "\"%s\"")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#endif /* __ITT_INTERNAL_BODY */
ITT_STUBV(ITTAPI, void, frame_begin,         (__itt_frame frame),     (ITT_FORMAT frame),  frame_begin,   __itt_group_frame, "%p")
ITT_STUBV(ITTAPI, void, frame_end,           (__itt_frame frame),     (ITT_FORMAT frame),  frame_end,     __itt_group_frame, "%p")

#ifndef __ITT_INTERNAL_BODY
#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUB(LIBITTAPI, __itt_event, event_createA, (const char    *name, int namelen), (ITT_FORMAT name, namelen), event_createA, __itt_group_mark | __itt_group_legacy, "\"%s\", %d")
ITT_STUB(LIBITTAPI, __itt_event, event_createW, (const wchar_t *name, int namelen), (ITT_FORMAT name, namelen), event_createW, __itt_group_mark | __itt_group_legacy, "\"%S\", %d")
#else  /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
ITT_STUB(LIBITTAPI, __itt_event, event_create,  (const char    *name, int namelen), (ITT_FORMAT name, namelen), event_create,  __itt_group_mark | __itt_group_legacy, "\"%s\", %d")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
ITT_STUB(LIBITTAPI, int,  event_start,          (__itt_event event),                (ITT_FORMAT event),         event_start,   __itt_group_mark | __itt_group_legacy, "%d")
ITT_STUB(LIBITTAPI, int,  event_end,            (__itt_event event),                (ITT_FORMAT event),         event_end,     __itt_group_mark | __itt_group_legacy, "%d")
#endif /* __ITT_INTERNAL_BODY */

#ifndef __ITT_INTERNAL_BODY
#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUB(ITTAPI, __itt_heap_function, heap_function_createA, (const char    *name, const char    *domain), (ITT_FORMAT name, domain), heap_function_createA, __itt_group_heap, "\"%s\", \"%s\"")
ITT_STUB(ITTAPI, __itt_heap_function, heap_function_createW, (const wchar_t *name, const wchar_t *domain), (ITT_FORMAT name, domain), heap_function_createW, __itt_group_heap, "\"%s\", \"%s\"")
#else  /* ITT_PLATFORM==ITT_PLATFORM_WIN */
ITT_STUB(ITTAPI, __itt_heap_function, heap_function_create,  (const char    *name, const char    *domain), (ITT_FORMAT name, domain), heap_function_create,  __itt_group_heap, "\"%s\", \"%s\"")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#endif /* __ITT_INTERNAL_BODY */
ITT_STUBV(ITTAPI, void, heap_allocate_begin,   (__itt_heap_function h, size_t size, int initialized),             (ITT_FORMAT h, size, initialized),       heap_allocate_begin, __itt_group_heap, "%p, %lu, %d")
ITT_STUBV(ITTAPI, void, heap_allocate_end,     (__itt_heap_function h, void* addr, size_t size, int initialized), (ITT_FORMAT h, addr, size, initialized), heap_allocate_end,   __itt_group_heap, "%p, %p, %lu, %d")
ITT_STUBV(ITTAPI, void, heap_free_begin,       (__itt_heap_function h, void* addr), (ITT_FORMAT h, addr), heap_free_begin, __itt_group_heap, "%p, %p")
ITT_STUBV(ITTAPI, void, heap_free_end,         (__itt_heap_function h, void* addr), (ITT_FORMAT h, addr), heap_free_end,   __itt_group_heap, "%p, %p")
ITT_STUBV(ITTAPI, void, heap_reallocate_begin, (__itt_heap_function h, void* addr, size_t new_size, int initialized),                 (ITT_FORMAT h, addr, new_size, initialized),           heap_reallocate_begin, __itt_group_heap, "%p, %p, %lu, %d")
ITT_STUBV(ITTAPI, void, heap_reallocate_end,   (__itt_heap_function h, void* addr, void* new_addr, size_t new_size, int initialized), (ITT_FORMAT h, addr, new_addr, new_size, initialized), heap_reallocate_end,   __itt_group_heap, "%p, %p, %p, %lu, %d")
ITT_STUBV(ITTAPI, void, heap_internal_access_begin, (void), (ITT_NO_PARAMS), heap_internal_access_begin, __itt_group_heap, "no args")
ITT_STUBV(ITTAPI, void, heap_internal_access_end,   (void), (ITT_NO_PARAMS), heap_internal_access_end,   __itt_group_heap, "no args")

/* legacy */
#ifndef __ITT_INTERNAL_BODY
#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUB(LIBITTAPI, int,  thr_name_setA, (const char    *name, int namelen), (ITT_FORMAT name, namelen), thr_name_setA, __itt_group_thread | __itt_group_legacy, "\"%s\", %d")
ITT_STUB(LIBITTAPI, int,  thr_name_setW, (const wchar_t *name, int namelen), (ITT_FORMAT name, namelen), thr_name_setW, __itt_group_thread | __itt_group_legacy, "\"%S\", %d")
#else  /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
ITT_STUB(LIBITTAPI, int,  thr_name_set,  (const char    *name, int namelen), (ITT_FORMAT name, namelen), thr_name_set,  __itt_group_thread | __itt_group_legacy, "\"%s\", %d")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
ITT_STUBV(LIBITTAPI, void, thr_ignore,   (void),                             (ITT_NO_PARAMS),            thr_ignore,    __itt_group_thread | __itt_group_legacy, "no args")

#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUBV(ITTAPI, void, sync_set_nameA, (void *addr, const char    *objtype, const char    *objname, int attribute), (ITT_FORMAT addr, objtype, objname, attribute), sync_set_nameA, __itt_group_sync | __itt_group_fsync | __itt_group_legacy, "%p, \"%s\", \"%s\", %x")
ITT_STUBV(ITTAPI, void, sync_set_nameW, (void *addr, const wchar_t *objtype, const wchar_t *objname, int attribute), (ITT_FORMAT addr, objtype, objname, attribute), sync_set_nameW, __itt_group_sync | __itt_group_fsync | __itt_group_legacy, "%p, \"%S\", \"%S\", %x")
#else  /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
ITT_STUBV(ITTAPI, void, sync_set_name,  (void *addr, const char    *objtype, const char    *objname, int attribute), (ITT_FORMAT addr, objtype, objname, attribute), sync_set_name,  __itt_group_sync | __itt_group_fsync | __itt_group_legacy, "p, \"%s\", \"%s\", %x")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */

#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUB(LIBITTAPI, int, notify_sync_nameA, (void *p, const char    *objtype, int typelen, const char    *objname, int namelen, int attribute), (ITT_FORMAT p, objtype, typelen, objname, namelen, attribute), notify_sync_nameA, __itt_group_sync | __itt_group_fsync | __itt_group_legacy, "%p, \"%s\", %d, \"%s\", %d, %x")
ITT_STUB(LIBITTAPI, int, notify_sync_nameW, (void *p, const wchar_t *objtype, int typelen, const wchar_t *objname, int namelen, int attribute), (ITT_FORMAT p, objtype, typelen, objname, namelen, attribute), notify_sync_nameW, __itt_group_sync | __itt_group_fsync | __itt_group_legacy, "%p, \"%S\", %d, \"%S\", %d, %x")
#else  /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
ITT_STUB(LIBITTAPI, int, notify_sync_name,  (void *p, const char    *objtype, int typelen, const char    *objname, int namelen, int attribute), (ITT_FORMAT p, objtype, typelen, objname, namelen, attribute), notify_sync_name,  __itt_group_sync | __itt_group_fsync | __itt_group_legacy, "%p, \"%s\", %d, \"%s\", %d, %x")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */

ITT_STUBV(LIBITTAPI, void, notify_sync_prepare,   (void *p), (ITT_FORMAT p), notify_sync_prepare,   __itt_group_sync | __itt_group_fsync | __itt_group_legacy, "%p")
ITT_STUBV(LIBITTAPI, void, notify_sync_cancel,    (void *p), (ITT_FORMAT p), notify_sync_cancel,    __itt_group_sync | __itt_group_fsync | __itt_group_legacy, "%p")
ITT_STUBV(LIBITTAPI, void, notify_sync_acquired,  (void *p), (ITT_FORMAT p), notify_sync_acquired,  __itt_group_sync | __itt_group_fsync | __itt_group_legacy, "%p")
ITT_STUBV(LIBITTAPI, void, notify_sync_releasing, (void *p), (ITT_FORMAT p), notify_sync_releasing, __itt_group_sync | __itt_group_fsync | __itt_group_legacy, "%p")
#endif /* __ITT_INTERNAL_BODY */

ITT_STUBV(LIBITTAPI, void, memory_read,   (void *addr, size_t size), (ITT_FORMAT addr, size), memory_read,   __itt_group_legacy, "%p, %lu")
ITT_STUBV(LIBITTAPI, void, memory_write,  (void *addr, size_t size), (ITT_FORMAT addr, size), memory_write,  __itt_group_legacy, "%p, %lu")
ITT_STUBV(LIBITTAPI, void, memory_update, (void *addr, size_t size), (ITT_FORMAT addr, size), memory_update, __itt_group_legacy, "%p, %lu")

ITT_STUB(LIBITTAPI, __itt_state_t,     state_get,    (void),                                    (ITT_NO_PARAMS),   state_get,    __itt_group_legacy, "no args")
ITT_STUB(LIBITTAPI, __itt_state_t,     state_set,    (__itt_state_t s),                         (ITT_FORMAT s),    state_set,    __itt_group_legacy, "%d")
ITT_STUB(LIBITTAPI, __itt_obj_state_t, obj_mode_set, (__itt_obj_prop_t p, __itt_obj_state_t s), (ITT_FORMAT p, s), obj_mode_set, __itt_group_legacy, "%d, %d")
ITT_STUB(LIBITTAPI, __itt_thr_state_t, thr_mode_set, (__itt_thr_prop_t p, __itt_thr_state_t s), (ITT_FORMAT p, s), thr_mode_set, __itt_group_legacy, "%d, %d")

/* internal */
#ifndef __ITT_INTERNAL_BODY
#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUB(ITTAPI, __itt_mark_type, mark_createA, (const char    *name), (ITT_FORMAT name), mark_createA, __itt_group_mark, "\"%s\"")
ITT_STUB(ITTAPI, __itt_mark_type, mark_createW, (const wchar_t *name), (ITT_FORMAT name), mark_createW, __itt_group_mark, "\"%S\"")
#else  /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
ITT_STUB(ITTAPI, __itt_mark_type, mark_create,  (const char    *name), (ITT_FORMAT name), mark_create,  __itt_group_mark, "\"%s\"")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#endif /* __ITT_INTERNAL_BODY */
#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUB(ITTAPI, int,  markA,        (__itt_mark_type mt, const char    *parameter), (ITT_FORMAT mt, parameter), markA, __itt_group_mark, "%d, \"%s\"")
ITT_STUB(ITTAPI, int,  markW,        (__itt_mark_type mt, const wchar_t *parameter), (ITT_FORMAT mt, parameter), markW, __itt_group_mark, "%d, \"%S\"")
#else  /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
ITT_STUB(ITTAPI, int,  mark,         (__itt_mark_type mt, const char    *parameter), (ITT_FORMAT mt, parameter), mark,  __itt_group_mark, "%d, \"%s\"")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
ITT_STUB(ITTAPI, int,  mark_off, (__itt_mark_type mt), (ITT_FORMAT mt), mark_off, __itt_group_mark, "%d")
#if ITT_PLATFORM==ITT_PLATFORM_WIN
ITT_STUB(ITTAPI, int,  mark_globalA, (__itt_mark_type mt, const char    *parameter), (ITT_FORMAT mt, parameter), mark_globalA, __itt_group_mark, "%d, \"%s\"")
ITT_STUB(ITTAPI, int,  mark_globalW, (__itt_mark_type mt, const wchar_t *parameter), (ITT_FORMAT mt, parameter), mark_globalW, __itt_group_mark, "%d, \"%S\"")
#else  /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
ITT_STUB(ITTAPI, int,  mark_global,  (__itt_mark_type mt, const char    *parameter), (ITT_FORMAT mt, parameter), mark_global,  __itt_group_mark, "%d, \"%S\"")
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
ITT_STUB(ITTAPI, int,  mark_global_off, (__itt_mark_type mt),                        (ITT_FORMAT mt),            mark_global_off, __itt_group_mark, "%d")

/* prototype */
/* empty so far */

/* hidden */
#ifndef __ITT_INTERNAL_BODY
ITT_STUB(ITTAPI, const char*, api_version, (void), (ITT_NO_PARAMS), api_version, __itt_group_all & ~__itt_group_legacy, "no args")
#endif /* __ITT_INTERNAL_BODY */
