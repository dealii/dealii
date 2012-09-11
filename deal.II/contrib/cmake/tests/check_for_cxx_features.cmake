INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckIncludeFiles)

#
# Check for various CXX features.
#



#
# gcc2.95 doesn't have the std::iterator class, but the standard
# requires it, so check whether we have to work around it
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <iterator>
  class MyIterator : public std::iterator<std::bidirectional_iterator_tag,int>{};
  int main(){return 0;}
  "
  HAVE_STD_ITERATOR_CLASS)


#
# Up to early gcc2.95 releases, the i/ostringstream classes were not
# available. check their availability, or whether we have to fall back
# to the old strstream classes.
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <sstream>
  std::istringstream i;
  std::ostringstream o;
  int main(){return 0;}
  "
  HAVE_STD_STRINGSTREAM)


#
# Check whether the numeric_limits classes are available
#
CHECK_CXX_SOURCE_COMPILES(
  "
  #include <limits>
  unsigned int i = std::numeric_limits<unsigned int>::min();
  int main(){return 0;}
  "
  HAVE_STD_NUMERIC_LIMITS)


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

CHECK_INCLUDE_FILES("ostream" HAVE_STD_OSTREAM_HEADER)
CHECK_INCLUDE_FILES("iosfwd" HAVE_STD_IOSFWD_HEADER)
CHECK_INCLUDE_FILES("stdint.h" HAVE_STDINT_H)
CHECK_INCLUDE_FILES("stdlib.h" HAVE_STDLIB_H)
CHECK_INCLUDE_FILES("strings.h" HAVE_STRINGS_H)
CHECK_INCLUDE_FILES("string.h" HAVE_STRING_H)
CHECK_INCLUDE_FILES("sys/stat.h" HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILES("sys/syscall.h" HAVE_SYS_SYSCALL_H)
CHECK_INCLUDE_FILES("sys/times.h" HAVE_SYS_TIMES_H)
CHECK_INCLUDE_FILES("sys/types.h" HAVE_SYS_TYPES_H)

