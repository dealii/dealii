INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckIncludeFiles)

#
# Check for various CXX features.
#



#
# Check whether the std::vector::iterator is just a plain pointer
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <vector>
  template <typename T> void f(T) {}
  template void f(int *);
  template void f(std::vector<int>::iterator);
  int main(){return 0;}
  "
  DEAL_II_VECTOR_ITERATOR_IS_NOT_POINTER)

IF(NOT DEAL_II_VECTOR_ITERATOR_IS_NOT_POINTER)
  SET(DEAL_II_VECTOR_ITERATOR_IS_POINTER 1)
ENDIF()


#
# Checks for various header files:
#

CHECK_INCLUDE_FILES("ostream" HAVE_STD_OSTREAM_HEADER) # TODO: remove this...
CHECK_INCLUDE_FILES("iosfwd" HAVE_STD_IOSFWD_HEADER)
CHECK_INCLUDE_FILES("stdint.h" HAVE_STDINT_H)
CHECK_INCLUDE_FILES("stdlib.h" HAVE_STDLIB_H)
CHECK_INCLUDE_FILES("strings.h" HAVE_STRINGS_H)
CHECK_INCLUDE_FILES("string.h" HAVE_STRING_H)
CHECK_INCLUDE_FILES("sys/stat.h" HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES("sys/syscall.h" HAVE_SYS_SYSCALL_H)
CHECK_INCLUDE_FILES("sys/times.h" HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILES("sys/types.h" HAVE_SYS_TYPES_H)

