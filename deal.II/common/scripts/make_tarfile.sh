######################################################################
# We assume that the most recent checked out version of the branch
# is in subdir checkout and configured.
#
# Will create both tar files (source and documentation)
#
######################################################################

MAJOR=`perl -n -e "m/(\\d+)\\.(\\d+)\\.(\\d+)/; print \\$1" checkout/Version`
MINOR=`perl -n -e "m/(\\d+)\\.(\\d+)\\.(\\d+)/; print \\$2" checkout/Version`
PATCH=`perl -n -e "m/(\\d+)\\.(\\d+)\\.(\\d+)/; print \\$3" checkout/Version`

cvs -d `cat checkout/CVS/Root` export -kv -r Version-$MAJOR-$MINOR-$PATCH deal.II
tar czf deal.II-$MAJOR.$MINOR.$PATCH.tar.gz deal.II
cd deal.II
./configure --with-umfpack
make online-doc
cd ..
tar czf deal.doc-$MAJOR.$MINOR.$PATCH.tar.gz deal.II/doc
