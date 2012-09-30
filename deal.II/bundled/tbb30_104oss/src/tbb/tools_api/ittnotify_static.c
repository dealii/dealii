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

#if ITT_PLATFORM==ITT_PLATFORM_WIN
#include <windows.h>
#else /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
#include <pthread.h>
#include <dlfcn.h>
#include <errno.h>
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "disable_warnings.h"

#define INTEL_NO_MACRO_BODY 
#include "ittnotify.h"
#include "legacy/ittnotify.h"
#include "internal/ittnotify.h"
#include "prototype/ittnotify.h"

#include "ittnotify_types.h"

#ifndef INTEL_ITTNOTIFY_PREFIX
#define INTEL_ITTNOTIFY_PREFIX __itt_
#endif /* INTEL_ITTNOTIFY_PREFIX */
#ifndef INTEL_ITTNOTIFY_POSTFIX
#define INTEL_ITTNOTIFY_POSTFIX _ptr_
#endif /* INTEL_ITTNOTIFY_POSTFIX */

#define _N_(n) ITT_JOIN(INTEL_ITTNOTIFY_PREFIX,n)

#ifndef CDECL
#if ITT_PLATFORM==ITT_PLATFORM_WIN
#define CDECL __cdecl
#else /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#define CDECL
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#endif /* CDECL */

#ifndef STDCALL
#if ITT_PLATFORM==ITT_PLATFORM_WIN
#define STDCALL __stdcall
#else /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
#define STDCALL
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#endif /* STDCALL */

#if ITT_PLATFORM==ITT_PLATFORM_WIN
typedef FARPROC   FPTR;
typedef DWORD     TIDT;
#else /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
typedef void*     FPTR;
typedef pthread_t TIDT;
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */

/* OS communication functions */
#if ITT_PLATFORM==ITT_PLATFORM_WIN
typedef HMODULE lib_t;
typedef CRITICAL_SECTION mutex_t;
#else /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
typedef void* lib_t;
typedef pthread_mutex_t mutex_t;
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */

static volatile long ittnotify_init = 0;
static lib_t ittnotify_lib = NULL;
static __itt_error_notification_t* error_handler = NULL;

#if ITT_OS==ITT_OS_WIN
static const char* ittnotify_lib_name = "libittnotify.dll";
#elif ITT_OS==ITT_OS_LINUX
static const char* ittnotify_lib_name = "libittnotify.so";
#elif ITT_OS==ITT_OS_MAC
static const char* ittnotify_lib_name = "libittnotify.dylib";
#else
#error Unsupported or unknown OS.
#endif

#ifndef LIB_VAR_NAME
#if ITT_ARCH==ITT_ARCH_IA32
#define LIB_VAR_NAME INTEL_LIBITTNOTIFY32
#else
#define LIB_VAR_NAME INTEL_LIBITTNOTIFY64
#endif
#endif /* LIB_VAR_NAME */

#if ITT_PLATFORM==ITT_PLATFORM_WIN
#define __itt_get_proc(lib, name) GetProcAddress(lib, name)
#define __itt_mutex_init(mutex)   InitializeCriticalSection(mutex)
#define __itt_mutex_lock(mutex)   EnterCriticalSection(mutex)
#define __itt_mutex_unlock(mutex) LeaveCriticalSection(mutex)
#define __itt_load_lib(name)      LoadLibraryA(name)
#define __itt_unload_lib(handle)  FreeLibrary(handle)
#define __itt_system_error()      (int)GetLastError()
#define __itt_fstrcmp(s1, s2)     lstrcmpA(s1, s2)
#define __itt_fstrlen(s)          lstrlenA(s)
#define __itt_fstrcpyn(s1, s2, l) lstrcpynA(s1, s2, l)
#define __itt_thread_id()         GetCurrentThreadId()
#define __itt_thread_yield()      SwitchToThread()
#ifndef ITT_SIMPLE_INIT
static int __itt_interlocked_increment(volatile int* ptr)
{
    ITT_BUILD_ASSERT(sizeof(int) == sizeof(long));
    return InterlockedIncrement((volatile long *)ptr);
}
#endif /* ITT_SIMPLE_INIT */
#else /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
#define __itt_get_proc(lib, name) dlsym(lib, name)
#define __itt_mutex_init(mutex)   \
    {                                                                                        \
        pthread_mutexattr_t mutex_attr;                                                      \
        int error_code = pthread_mutexattr_init(&mutex_attr);                                \
        if (error_code)                                                                      \
            __itt_report_error(__itt_error_system, "pthread_mutexattr_init", error_code);    \
        error_code = pthread_mutexattr_settype(&mutex_attr, PTHREAD_MUTEX_RECURSIVE);        \
        if (error_code)                                                                      \
            __itt_report_error(__itt_error_system, "pthread_mutexattr_settype", error_code); \
        error_code = pthread_mutex_init(mutex, &mutex_attr);                                 \
        if (error_code)                                                                      \
            __itt_report_error(__itt_error_system, "pthread_mutex_init", error_code);        \
        error_code = pthread_mutexattr_destroy(&mutex_attr);                                 \
        if (error_code)                                                                      \
            __itt_report_error(__itt_error_system, "pthread_mutexattr_destroy", error_code); \
    }
#define __itt_mutex_lock(mutex)   pthread_mutex_lock(mutex)
#define __itt_mutex_unlock(mutex) pthread_mutex_unlock(mutex)
#define __itt_load_lib(name)      dlopen(name, RTLD_LAZY)
#define __itt_unload_lib(handle)  dlclose(handle)
#define __itt_system_error()      errno
#define __itt_fstrcmp(s1, s2)     strcmp(s1, s2)
#define __itt_fstrlen(s)          strlen(s)
#define __itt_fstrcpyn(s1, s2, l) strncpy(s1, s2, l)
#define __itt_thread_id()         pthread_self()
#define __itt_thread_yield()      sched_yield()
#if ITT_ARCH==ITT_ARCH_IA64
#ifdef __INTEL_COMPILER
#define __TBB_machine_fetchadd4(addr, val) __fetchadd4_acq((void *)addr, val)
#else  /* __INTEL_COMPILER */
// TODO: Add Support for not Intel compilers for IA64
#endif /* __INTEL_COMPILER */
#else /* ITT_ARCH!=ITT_ARCH_IA64 */
#ifndef ITT_SIMPLE_INIT
static int __TBB_machine_fetchadd4(volatile void* ptr, int addend)
{
    int result;
    __asm__ __volatile__("lock\nxaddl %0,%1"
                          : "=r"(result),"=m"(*(int *)ptr)
                          : "0"(addend), "m"(*(int *)ptr)
                          : "memory");
    return result;
}
#endif // ITT_SIMPLE_INIT
#endif /* ITT_ARCH==ITT_ARCH_IA64 */
#ifndef ITT_SIMPLE_INIT
static int __itt_interlocked_increment(volatile int* ptr)
{
    return __TBB_machine_fetchadd4(ptr, 1) + 1;
}
#endif /* ITT_SIMPLE_INIT */
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */

const int _N_(err) = 0;

typedef int (__itt_init_ittlib_t)(const char*, __itt_group_id);

/* this define used to control initialization function name. */
#ifndef __itt_init_ittlib_name
static int _N_(init_ittlib)(const char*, __itt_group_id);
static __itt_init_ittlib_t* __itt_init_ittlib_ptr = _N_(init_ittlib);
#define __itt_init_ittlib_name __itt_init_ittlib_ptr
#endif /* __itt_init_ittlib_name */

/* building pointers to imported funcs */
#undef ITT_STUBV
#undef ITT_STUB
#define ITT_STUB(api,type,name,args,params,ptr,group,format)      \
    static type api ITT_JOIN(_N_(name),_init) args;               \
    typedef type api name##_t args;                               \
    extern "C" name##_t* ITTNOTIFY_NAME(name);                    \
    name##_t* ITTNOTIFY_NAME(name) = ITT_JOIN(_N_(name),_init);   \
    static type api ITT_JOIN(_N_(name),_init) args                \
    {                                                             \
        if (__itt_init_ittlib_name(NULL, __itt_group_none)        \
            && ITTNOTIFY_NAME(name)                               \
            && ITTNOTIFY_NAME(name) != ITT_JOIN(_N_(name),_init)) \
            return ITTNOTIFY_NAME(name) params;                   \
        else                                                      \
            return (type)0;                                       \
    }

#define ITT_STUBV(api,type,name,args,params,ptr,group,format)     \
    static type api ITT_JOIN(_N_(name),_init) args;               \
    typedef type api name##_t args;                               \
    extern "C" name##_t* ITTNOTIFY_NAME(name);                    \
    name##_t* ITTNOTIFY_NAME(name) = ITT_JOIN(_N_(name),_init);   \
    static type api ITT_JOIN(_N_(name),_init) args                \
    {                                                             \
        if (__itt_init_ittlib_name(NULL, __itt_group_none)        \
            && ITTNOTIFY_NAME(name)                               \
            && ITTNOTIFY_NAME(name) != ITT_JOIN(_N_(name),_init)) \
            ITTNOTIFY_NAME(name) params;                          \
        else                                                      \
            return;                                               \
    }

/* Define types and *_init functions. */
#include "ittnotify_static.h"

ITT_GROUP_LIST(group_list);

typedef struct __itt_group_alias_
{
    const char*    env_var;
    __itt_group_id groups;
} __itt_group_alias;

static __itt_group_alias group_alias[] = {
    { "KMP_FOR_TPROFILE", (__itt_group_id)(__itt_group_control | __itt_group_thread | __itt_group_sync  | __itt_group_mark) },
    { "KMP_FOR_TCHECK",   (__itt_group_id)(__itt_group_control | __itt_group_thread | __itt_group_fsync | __itt_group_mark) },
    { NULL,               (__itt_group_none) }
};

typedef struct __itt_func_map_
{
    const char*    name;
    void**         func_ptr;
    __itt_group_id group;
} __itt_func_map;

#define __ptr_(pname,name,group) {ITT_TO_STR(ITT_JOIN(__itt_,pname)), (void**)(void*)&ITTNOTIFY_NAME(name), (__itt_group_id)(group)},
#undef ITT_STUB
#undef ITT_STUBV
#define ITT_STUB(api,type,name,args,params,nameindll,group,format) __ptr_(nameindll,name,group)
#define ITT_STUBV ITT_STUB

static __itt_func_map func_map[] = {
#include "ittnotify_static.h"
    {NULL, NULL, __itt_group_none}
};

#ifndef ITT_SIMPLE_INIT

#undef ITT_STUBV
#undef ITT_STUB
#define ITT_STUBV(api,type,name,args,params,ptr,group,format) \
ITT_EXTERN_C type api _N_(name) args                          \
{                                                             \
    if (ITTNOTIFY_NAME(name))                                 \
        ITTNOTIFY_NAME(name) params;                          \
    else                                                      \
        return;                                               \
}

#define ITT_STUB(api,type,name,args,params,ptr,group,format) \
ITT_EXTERN_C type api _N_(name) args                         \
{                                                            \
    if (ITTNOTIFY_NAME(name))                                \
        return ITTNOTIFY_NAME(name) params;                  \
    else                                                     \
        return (type)0;                                      \
}

/* Define ITT functions. */
#include "ittnotify_static.h"

#endif /* ITT_SIMPLE_INIT */

static const char* __itt_fsplit(const char* s, const char* sep, const char** out, int* len)
{
    int i;
    int j;

    if (!s || !sep || !out || !len)
        return 0;

    for (i = 0; s[i]; i++)
    {
        int b = 0;
        for (j = 0; sep[j]; j++)
            if (s[i] == sep[j])
            {
                b = 1;
                break;
            }
        if (!b)
            break;
    }

    if (!s[i])
        return 0;

    *len = 0;
    *out = s + i;

    for (; s[i]; i++, (*len)++)
    {
        int b = 0;
        for (j = 0; sep[j]; j++)
            if (s[i] == sep[j])
            {
                b = 1;
                break;
            }
        if (b)
            break;
    }

    for (; s[i]; i++)
    {
        int b = 0;
        for (j = 0; sep[j]; j++)
            if (s[i] == sep[j])
            {
                b = 1;
                break;
            }
        if (!b)
            break;
    }

    return s + i;
}

#ifdef ITT_NOTIFY_EXT_REPORT
ITT_EXTERN_C void _N_(error_handler)(__itt_error_code, va_list args);
#endif /* ITT_NOTIFY_EXT_REPORT */

static void __itt_report_error(__itt_error_code code, ...)
{
    va_list args;
    va_start( args, code );
    if (error_handler != NULL)
        error_handler(code, args);
#ifdef ITT_NOTIFY_EXT_REPORT
    _N_(error_handler)(code, args);
#endif /* ITT_NOTIFY_EXT_REPORT */
    va_end(args);
}

static const char* __itt_get_env_var(const char* name)
{
#define MAX_ENV_VALUE_SIZE 4086
    static char  env_buff[MAX_ENV_VALUE_SIZE];
    static char* env_value = (char*)&env_buff;

    if (name != NULL)
    {
#if ITT_PLATFORM==ITT_PLATFORM_WIN
        size_t max_len = MAX_ENV_VALUE_SIZE - ((size_t)env_value - (size_t)&env_buff);
        DWORD rc = GetEnvironmentVariableA(name, env_value, (DWORD)max_len);
        if (rc >= max_len)
        {
            __itt_report_error(__itt_error_env_too_long, name, (size_t)rc - 1, (size_t)(max_len - 1));
        }
        else if (rc > 0)
        {
            char* ret = env_value;
            env_value += rc + 1;
            return ret;
        }
        else
        {
            /* If environment variable is empty, GetEnvirornmentVariables() returns zero (number of   */
            /* characters (not including terminating null), and GetLastError() returns ERROR_SUCCESS. */
            DWORD err = GetLastError();
            if (err == ERROR_SUCCESS)
                return env_value;

            if (err != ERROR_ENVVAR_NOT_FOUND)
                __itt_report_error(__itt_error_cant_read_env, name, (int)err);
        }
#else  /* ITT_PLATFORM!=ITT_PLATFORM_WIN */
        char* env = getenv(name);
        if (env != NULL)
        {
            size_t len = strlen(env);
            size_t max_len = MAX_ENV_VALUE_SIZE - ((size_t)env_value - (size_t)&env_buff);
            if (len < max_len)
            {
                char* ret = env_value;
                strncpy(env_value, env, len + 1);
                env_value += len + 1;
                return ret;
            } else
                __itt_report_error(__itt_error_env_too_long, name, (size_t)len, (size_t)(max_len - 1));
        }
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
    }
    return NULL;
}

static const char* __itt_get_lib_name()
{
    const char* lib_name = __itt_get_env_var(ITT_TO_STR(LIB_VAR_NAME));
    return (lib_name == NULL) ? ittnotify_lib_name : lib_name;
}

#ifndef min
#define min(a,b) (a) < (b) ? (a) : (b)
#endif /* min */

static __itt_group_id __itt_get_groups()
{
    int i;
    __itt_group_id res = __itt_group_none;

    const char* var_name  = "INTEL_ITTNOTIFY_GROUPS";
    const char* group_str = __itt_get_env_var(var_name);
    if (group_str != NULL)
    {
        int len;
        char gr[255];
        const char* chunk;
        while ((group_str = __itt_fsplit(group_str, ",; ", &chunk, &len)) != NULL)
        {
            __itt_fstrcpyn(gr, chunk, sizeof(gr));

            gr[min((size_t)len, sizeof(gr) - 1)] = 0;

            for (i = 0; group_list[i].name != NULL; i++)
            {
                if (!__itt_fstrcmp(gr, group_list[i].name))
                {
                    res = (__itt_group_id)(res | group_list[i].id);
                    break;
                }
            }
        }
        /* TODO: !!! Workaround for bug with warning for unknown group !!!
         * Should be fixed in new initialization scheme.
         * Now the following groups should be set always.
         */
        for (i = 0; group_list[i].id != __itt_group_none; i++)
            if (group_list[i].id != __itt_group_all && group_list[i].id > __itt_group_splitter)
                res = (__itt_group_id)(res | group_list[i].id);
        return res;
    }
    else
    {
        for (i = 0; group_alias[i].env_var != NULL; i++)
            if (__itt_get_env_var(group_alias[i].env_var) != NULL)
                return group_alias[i].groups;
    }

    return res;
}

static int __itt_is_legacy_lib(lib_t lib)
{
    if (lib == NULL)
        return 0; // if unknown assume NO

    if (__itt_get_proc(lib, "__itt_api_version"))
        return 0; // New interface - NO
    return 1; // It's legacy otherwise
}

#if ITT_PLATFORM==ITT_PLATFORM_WIN
#pragma warning(push)
#pragma warning(disable: 4054)
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */

/* ITT_EXTERN_C - should be exported after agreament
static void _N_(fini_ittlib)(void)
{
    int i;

    if (ittnotify_init)
    {
        // Clear all pointers
        for (i = 0; func_map[i].name != NULL; i++)
            *func_map[i].func_ptr = NULL;

        if (ittnotify_lib != NULL)
            __itt_unload_lib(ittnotify_lib);

        ittnotify_lib  = NULL;
        ittnotify_init = 0;
    }
}
*/

static int _N_(init_ittlib)(const char* lib_name, __itt_group_id groups)
{
    int i, ret = 0;
    static volatile TIDT current_thread = 0;

    if (!ittnotify_init)
    {
#ifndef ITT_SIMPLE_INIT
        static mutex_t mutex;
        static volatile int inter_counter = 0;
        static volatile int mutex_initialized = 0;

        if (!mutex_initialized)
        {
            if (__itt_interlocked_increment(&inter_counter) == 1)
            {
                __itt_mutex_init(&mutex);
                mutex_initialized = 1;
            }
            else
                while (!mutex_initialized)
                    __itt_thread_yield();
        }

        __itt_mutex_lock(&mutex);
#endif /* ITT_SIMPLE_INIT */

        if (!ittnotify_init)
        {
            if (current_thread == 0)
            {
                current_thread = __itt_thread_id();
                if (groups == __itt_group_none)
                    groups = __itt_get_groups();
                if (groups == __itt_group_none)
                {
                    // Clear all pointers
                    for (i = 0; func_map[i].name != NULL; i++ )
                        *func_map[i].func_ptr = NULL;
                }
                else
                {
                    __itt_group_id zero_group = __itt_group_none;
                    if (lib_name == NULL)
                        lib_name = __itt_get_lib_name();
                    ittnotify_lib = __itt_load_lib(lib_name);
                    if (ittnotify_lib != NULL)
                    {
                        if (__itt_is_legacy_lib(ittnotify_lib))
                            groups = __itt_group_legacy;

                        for (i = 0; func_map[i].name != NULL; i++)
                        {
                            if (func_map[i].group & groups)
                            {
                                *func_map[i].func_ptr = (void*)__itt_get_proc(ittnotify_lib, func_map[i].name);
                                if (*func_map[i].func_ptr == NULL)
                                {
                                    __itt_report_error(__itt_error_no_symbol, lib_name, func_map[i].name );
                                    zero_group = (__itt_group_id)(zero_group | func_map[i].group);
                                }
                            }
                            else
                                *func_map[i].func_ptr = NULL;
                        }

                        if (groups == __itt_group_legacy)
                        {
                            // Compatibility with legacy tools
                            ITTNOTIFY_NAME(sync_prepare)   = ITTNOTIFY_NAME(notify_sync_prepare);
                            ITTNOTIFY_NAME(sync_cancel)    = ITTNOTIFY_NAME(notify_sync_cancel);
                            ITTNOTIFY_NAME(sync_acquired)  = ITTNOTIFY_NAME(notify_sync_acquired);
                            ITTNOTIFY_NAME(sync_releasing) = ITTNOTIFY_NAME(notify_sync_releasing);
                        }
                    }
                    else
                    {
                        // Clear all pointers
                        for (i = 0; func_map[i].name != NULL; i++)
                            *func_map[i].func_ptr = NULL;

                        __itt_report_error(__itt_error_no_module, lib_name,
#if ITT_PLATFORM==ITT_PLATFORM_WIN
                            __itt_system_error()
#else  /* ITT_PLATFORM==ITT_PLATFORM_WIN */
                            dlerror()
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
                        );
                    }
#ifdef ITT_COMPLETE_GROUP
                    for (i = 0; func_map[i].name != NULL; i++)
                        if (func_map[i].group & zero_group)
                            *func_map[i].func_ptr = NULL;
#endif /* ITT_COMPLETE_GROUP */

                    /* evaluating if any function ptr is non empty */
                    for (i = 0; func_map[i].name != NULL; i++)
                    {
                        if (*func_map[i].func_ptr != NULL)
                        {
                            ret = 1;
                            break;
                        }
                    }
                }

                ittnotify_init = 1;
                current_thread = 0;
            }
        }

#ifndef ITT_SIMPLE_INIT
        __itt_mutex_unlock(&mutex);
#endif /* ITT_SIMPLE_INIT */
    }

    return ret;
}

ITT_EXTERN_C __itt_error_notification_t* _N_(set_error_handler)(__itt_error_notification_t* handler)
{
    __itt_error_notification_t* prev = error_handler;
    error_handler = handler;
    return prev;
}

#if ITT_PLATFORM==ITT_PLATFORM_WIN
#pragma warning(pop)
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
