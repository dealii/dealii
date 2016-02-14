## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE at
## the top level of the deal.II distribution.
##
## ---------------------------------------------------------------------

use strict;

my $gallery_file = shift;
open GALLERY, "<$gallery_file";

my $gallery_dir = shift;


# Print the first part of gallery.h.in up until the point where we
# find the line with '@@GALLERY_LIST@@'
while (my $line = <GALLERY>)
{
  last if($line =~ m/\@\@GALLERY_LIST\@\@/);
  print $line;
}

# print the list of code gallery programs as a descriptor/description
# list
print "<dl>\n";
foreach my $gallery (sort @ARGV)
{
    my $gallery_underscore = $gallery;

    open AUTHOR, "<$gallery_dir/$gallery/doc/author";
    my $authors;
    while (my $line = <AUTHOR>) {
        chop $line;
        $authors .= $line . ", ";
    }

    # remove trailing whitespaces, as well as the trailing comma
    chop $authors;
    $authors =~ s/,$//;

    $gallery_underscore    =~ s/-/_/;
    print "  <dt><b>\@ref code_gallery_${gallery_underscore} \"$gallery\"</b> (by $authors)</dt>\n";
    print "    <dd>\n";
    open TOOLTIP, "<$gallery_dir/$gallery/doc/tooltip";
    while (my $line = <TOOLTIP>) {
        print "      $line";
    }
    print "    </dd>\n";
    print "\n";
}
print "</dl>\n";


# Then print the rest of code-gallery.h.in
while (my $line = <GALLERY>)
{
  print $line;
}
close GALLERY;
