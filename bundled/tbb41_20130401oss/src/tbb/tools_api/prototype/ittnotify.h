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

#ifndef _PROTOTYPE_ITTNOTIFY_H_
#define _PROTOTYPE_ITTNOTIFY_H_
/**
 * @file
 * @brief Prototype User API functions and types
 */

/** @cond exclude_from_documentation */
#ifndef ITT_OS_WIN
#  define ITT_OS_WIN   1
#endif /* ITT_OS_WIN */

#ifndef ITT_OS_LINUX
#  define ITT_OS_LINUX 2
#endif /* ITT_OS_LINUX */

#ifndef ITT_OS_MAC
#  define ITT_OS_MAC   3
#endif /* ITT_OS_MAC */

#ifndef ITT_OS
#  if defined WIN32 || defined _WIN32
#    define ITT_OS ITT_OS_WIN
#  elif defined( __APPLE__ ) && defined( __MACH__ )
#    define ITT_OS ITT_OS_MAC
#  else
#    define ITT_OS ITT_OS_LINUX
#  endif
#endif /* ITT_OS */

#ifndef ITT_PLATFORM_WIN
#  define ITT_PLATFORM_WIN 1
#endif /* ITT_PLATFORM_WIN */

#ifndef ITT_PLATFORM_POSIX
#  define ITT_PLATFORM_POSIX 2
#endif /* ITT_PLATFORM_POSIX */

#ifndef ITT_PLATFORM
#  if ITT_OS==ITT_OS_WIN
#    define ITT_PLATFORM ITT_PLATFORM_WIN
#  else
#    define ITT_PLATFORM ITT_PLATFORM_POSIX
#  endif /* _WIN32 */
#endif /* ITT_PLATFORM */

#include <stddef.h>
#include <stdarg.h>
#if ITT_PLATFORM==ITT_PLATFORM_WIN
#include <tchar.h>
#endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */

#ifndef CDECL
#  if ITT_PLATFORM==ITT_PLATFORM_WIN
#    define CDECL __cdecl
#  else  /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#    define CDECL /* nothing */
#  endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#endif /* CDECL */

#ifndef STDCALL
#  if ITT_PLATFORM==ITT_PLATFORM_WIN
#    define STDCALL __stdcall
#  else  /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#    define STDCALL /* nothing */
#  endif /* ITT_PLATFORM==ITT_PLATFORM_WIN */
#endif /* STDCALL */

#define ITTAPI_CALL    CDECL
#define LIBITTAPI_CALL /* nothing */

#define ITT_JOIN_AUX(p,n) p##n
#define ITT_JOIN(p,n)     ITT_JOIN_AUX(p,n)

#ifndef INTEL_ITTNOTIFY_PREFIX
#  define INTEL_ITTNOTIFY_PREFIX __itt_
#endif /* INTEL_ITTNOTIFY_PREFIX */
#ifndef INTEL_ITTNOTIFY_POSTFIX
#  define INTEL_ITTNOTIFY_POSTFIX _ptr_
#endif /* INTEL_ITTNOTIFY_POSTFIX */

#define ITTNOTIFY_NAME_AUX(n) ITT_JOIN(INTEL_ITTNOTIFY_PREFIX,n)
#define ITTNOTIFY_NAME(n)     ITTNOTIFY_NAME_AUX(ITT_JOIN(n,INTEL_ITTNOTIFY_POSTFIX))

#define ITTNOTIFY_VOID(n) (!ITTNOTIFY_NAME(n)) ? (void)0 : ITTNOTIFY_NAME(n)
#define ITTNOTIFY_DATA(n) (!ITTNOTIFY_NAME(n)) ?       0 : ITTNOTIFY_NAME(n)

#ifdef ITT_STUB
#undef ITT_STUB
#endif
#ifdef ITT_STUBV
#undef ITT_STUBV
#endif
#define ITT_STUBV(api,type,name,args,params)                      \
    typedef type (api* ITT_JOIN(ITTNOTIFY_NAME(name),_t)) args;   \
    extern ITT_JOIN(ITTNOTIFY_NAME(name),_t) ITTNOTIFY_NAME(name);
#define ITT_STUB ITT_STUBV

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */
/** @endcond */

/**
 * @defgroup prototype Prototype API
 * @{
 * @}
 */

/****************************************************************************
 * ??? group
 ****************************************************************************/

/** @cond exclude_from_documentation */
#ifdef __cplusplus
}
#endif /* __cplusplus */
/** @endcond */

#endif /* _PROTOTYPE_ITTNOTIFY_H_ */
