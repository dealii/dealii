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

svn export http://wolfgang.math.tamu.edu/svn/public/deal.II/tags/Version-$MAJOR-$MINOR-$PATCH/deal.II

tar czf deal.II-$MAJOR.$MINOR.$PATCH.tar.gz deal.II
cd deal.II
./configure --with-umfpack
make online-doc
cd ..
tar czf deal.doc-$MAJOR.$MINOR.$PATCH.tar.gz deal.II/doc
