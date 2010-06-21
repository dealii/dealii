######################################################################
# $Id$
######################################################################

$/ = undef;

# Read source for web page
open TOC, "<toc.html.in";
$toc = <TOC>;
close TOC;

# Read generated map file
open MAP, "<steps.cmapx";
$map = <MAP>;
close MAP;

# Insert contents of map file for @@MAP@@
$toc =~ s/\@\@MAP\@\@/$map/;

print $toc;
