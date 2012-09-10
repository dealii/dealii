INCLUDE(CheckCXXSourceCompiles)

#
# Check for various CXX Features.
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
