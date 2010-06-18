######################################################################
# $Id$
######################################################################

$/ = undef;

open TOC, "<toc.html.in";
$toc = <TOC>;
close TOC;

open MAP, "<steps.cmapx";
$map = <MAP>;
close MAP;

$toc =~ s/\@\@MAP\@\@/$map/;

print $toc;
