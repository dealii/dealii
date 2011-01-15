######################################################################
# creates both tar files (source and documentation). call this script
# in an otherwise empty directory as follows:
#
#   bash make_tarfile.sh 6 0 0
#
# to generate the tar files for release 6.0.0
#
######################################################################

MAJOR=$1
MINOR=$2
PATCH=$3

svn export http://www.dealii.org/svn/dealii/tags/Version-$MAJOR-$MINOR-$PATCH/deal.II

# Generate distribution tarfile

tar czf deal.II-$MAJOR.$MINOR.$PATCH.tar.gz deal.II

# Generate stripped tarfile

mv deal.II/doc deal.II/examples .
tar czf deal.nodoc-$MAJOR.$MINOR.$PATCH.tar.gz deal.II
mv doc examples deal.II

# Generate documentation

cd deal.II
./configure --with-umfpack
make online-doc
cd ..
tar czf deal.doc-$MAJOR.$MINOR.$PATCH.tar.gz deal.II/doc

