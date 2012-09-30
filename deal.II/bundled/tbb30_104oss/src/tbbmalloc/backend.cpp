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

// intrin.h available since VS2005
#if defined(_MSC_VER) && _MSC_VER >= 1400
#define __TBB_HAS_INTRIN_H 1
#else
#define __TBB_HAS_INTRIN_H 0
#endif

#if __TBB_HAS_INTRIN_H
#include <intrin.h>   /* for __cpuid */
#endif

#include "tbbmalloc_internal.h"
//! Define the main synchronization method
/** It should be specified before including LifoList.h */
#define FINE_GRAIN_LOCKS
#include "LifoList.h"


namespace rml {
namespace internal {

// If USE_MALLOC_FOR_LARGE_OBJECT is nonzero, then large allocations are done via malloc.
// Otherwise large allocations are done using the scalable allocator's block allocator.
// As of 06.Jun.17, using malloc is about 10x faster on Linux.
#if !_WIN32
#define USE_MALLOC_FOR_LARGE_OBJECT 1
#endif

/*********** Code to acquire memory from the OS or other executive ****************/

#if USE_DEFAULT_MEMORY_MAPPING
#include "MapMemory.h"
#else
/* assume MapMemory and UnmapMemory are customized */
#endif

#if USE_MALLOC_FOR_LARGE_OBJECT

// (get|free)RawMemory only necessary for the USE_MALLOC_FOR_LARGE_OBJECT case
void* getRawMemory (size_t size, bool useMapMem = false)
{
    void *object;

    if (useMapMem) 
        object = MapMemory(size);
    else
#if MALLOC_CHECK_RECURSION
    if (RecursiveMallocCallProtector::noRecursion())
        object = malloc(size);
    else if ( rml::internal::original_malloc_found )
        object = (*rml::internal::original_malloc_ptr)(size);
    else
        object = MapMemory(size);
#else
    object = malloc(size);
#endif /* MALLOC_CHECK_RECURSION */
    return object;
}

void freeRawMemory (void *object, size_t size, bool useMapMem)
{
    if (useMapMem)
        UnmapMemory(object, size);
    else
#if MALLOC_CHECK_RECURSION
    if (RecursiveMallocCallProtector::noRecursion())
        free(object);
    else if ( rml::internal::original_malloc_found )
        (*rml::internal::original_free_ptr)(object);
    else
        UnmapMemory(object, size);
#else
    free(object);
#endif /* MALLOC_CHECK_RECURSION */
}

#else /* USE_MALLOC_FOR_LARGE_OBJECT */

void* getRawMemory (size_t size, bool = false) { return MapMemory(size); }

void freeRawMemory (void *object, size_t size, bool) {
    UnmapMemory(object, size);
}

#endif /* USE_MALLOC_FOR_LARGE_OBJECT */

/********* End memory acquisition code ********************************/

static unsigned int getCPUid()
{
    unsigned int id;

#if (__ARCH_x86_32||__ARCH_x86_64) && (__linux__||__APPLE__||__FreeBSD__||__sun||__MINGW32__)
    int res;
 #if __ARCH_x86_32
    /* EBX used for PIC support. Having EAX in output operands 
       prevents ICC from crash like in __TBB_ICC_ASM_VOLATILE_BROKEN. */
    int _eax, _ecx, _edx;
    __asm__ ("xchgl %%ebx, %1\n\t"
             "cpuid\n\t"
             "xchgl %%ebx, %1\n\t"
             : "=a" (_eax), "=r" (res)
             : "a" (1) : "ecx", "edx");
 #else
    __asm__ ("cpuid\n\t"
             : "=b" (res)
             : "a" (1) );
 #endif // __ARCH_x86_32
    id = (res >> 24) & 0xff;
#elif _WIN32 || _WIN64
 #if __TBB_HAS_INTRIN_H
    int CPUInfo[4];
    __cpuid(CPUInfo, 1);
    id = (CPUInfo[1] >> 24) & 0xff;
 #else
    int res;
    _asm {
        push ebx
        push ecx
        mov  eax,1
        cpuid
        mov  res,ebx
        pop  ecx
        pop  ebx
    }
    id = (res >> 24) & 0xff;
 #endif
# else
    id = getThreadId();
#endif
    return id;
}


/* 
 * To decrease contention for free blocks, free blocks are split, and access
 * to them is based on process number.
 */
const int numOfFreeBlockLists = 4;

/*
 * This is a LIFO linked list that one can init, push or pop from
 */
static LifoList freeBlockList[numOfFreeBlockLists];

FreeBlocks freeBlocks;

bool FreeBlocks::bootstrap(RawAlloc myAlloc, RawFree myFree, size_t /*myReqSize*/)
{
    if (!myAlloc && !myFree) {
        rawAlloc = getRawMemory;
        rawFree = freeRawMemory;
        // Get virtual memory in pieces of this size: 0x0100000 is 1 megabyte decimal
        memReqSize = 0x0100000;
    } else
        MALLOC_ASSERT(0, "Not implemented yet.");
    return mallocBigBlock();
}

BlockI *FreeBlocks::get(bool startup)
{
    BlockI *bigBlock;
    // must not call getCPUid during malloc initialization 
    // because getCPUid can call malloc
    const unsigned myFreeList = startup? 0 : getCPUid()%numOfFreeBlockLists;
    unsigned currListIdx = myFreeList;

    do {
        if (bigBlock = (BlockI *) freeBlockList[currListIdx].pop()) {
            MALLOC_ITT_SYNC_ACQUIRED(freeBlockList+currListIdx);
            break;
        }
        currListIdx = (currListIdx+1) % numOfFreeBlockLists;
    } while (currListIdx != myFreeList);

    while (!bigBlock) {
        /* We are out of blocks so go to the OS and get another one */
        if (!mallocBigBlock()) return NULL;

        bigBlock = (BlockI *) freeBlockList[myFreeList].pop();
        if (bigBlock)
            MALLOC_ITT_SYNC_ACQUIRED(freeBlockList+myFreeList);
    }

    return bigBlock;
}

void FreeBlocks::put(BlockI *ptr, bool startup)
{
    unsigned myFreeList = startup? 0 : getCPUid()%numOfFreeBlockLists;
    MALLOC_ITT_SYNC_RELEASING(freeBlockList+myFreeList);
    freeBlockList[myFreeList].push((void **)ptr);
}

void FreeBlocks::putList(BlockI *head, BlockI *tail)
{
    unsigned myFreeList = getCPUid()%numOfFreeBlockLists;
    MALLOC_ITT_SYNC_RELEASING(freeBlockList+myFreeList);
    freeBlockList[myFreeList].pushList((void**)head, (void**)tail);
}

/*
 * Big Blocks are the blocks we get from the OS or some similar place using getMemory above.
 * They are placed on the freeBlockList once they are acquired.
 */
bool FreeBlocks::mallocBigBlock()
{
/* Divide the big block into smaller bigBlocks that hold that many blocks.
 * This is done since we really need a lot of blocks on the freeBlockList 
 * or there will be contention problems.
 */
    const unsigned int blocksPerBigBlock = 16/numOfFreeBlockLists;

    void *unalignedBigBlock = (*rawAlloc)(memReqSize, /*useMapMem=*/true);

    if (!unalignedBigBlock) {
        TRACEF(( "[ScalableMalloc trace] in mallocBigBlock, getMemory returns 0\n" ));
        /* We can't get any more memory from the OS or executive */
        return false;
    }

    void *alignedBigBlock = alignUp(unalignedBigBlock, blockSize);
    void *bigBlockCeiling = (void*)((uintptr_t)unalignedBigBlock + memReqSize);

    size_t bigBlockSplitSize = blocksPerBigBlock * blockSize;

    BlockI *splitBlock = (BlockI*)alignedBigBlock;

    // distribute alignedBigBlock between all freeBlockList elements
    for (unsigned currListIdx = 0;
         ((uintptr_t)splitBlock + blockSize) <= (uintptr_t)bigBlockCeiling;
         currListIdx = (currListIdx+1) % numOfFreeBlockLists) {
        void *splitEdge = (void*)((uintptr_t)splitBlock + bigBlockSplitSize);
        if( splitEdge > bigBlockCeiling) {
            splitEdge = alignDown(bigBlockCeiling, blockSize);
        }
        ((BlockI*)splitBlock)->initialize(splitEdge);
        MALLOC_ITT_SYNC_RELEASING(freeBlockList+currListIdx);
        freeBlockList[currListIdx].push((void**) splitBlock);
        splitBlock = (BlockI*)splitEdge;
    }

    TRACEF(( "[ScalableMalloc trace] in mallocBigBlock returning 1\n" ));
    return true;
}

} } // namespaces
