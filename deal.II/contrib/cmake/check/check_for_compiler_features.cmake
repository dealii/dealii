INCLUDE(CheckCXXSourceCompiles)
INCLUDE(CheckIncludeFiles)


#
# Check for various compiler features.
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
  SET(DEAL_II_VECTOR_ITERATOR_IS_POINTER TRUE)
ENDIF()
