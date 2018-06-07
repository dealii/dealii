## ---------------------------------------------------------------------
##
## Copyright (C) 2016 - 2018 by the deal.II authors
##
## This file is part of the deal.II library.
##
## The deal.II library is free software; you can use it, redistribute
## it, and/or modify it under the terms of the GNU Lesser General
## Public License as published by the Free Software Foundation; either
## version 2.1 of the License, or (at your option) any later version.
## The full text of the license can be found in the file LICENSE.md at
## the top level directory of deal.II.
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

# create a list of code gallery program descriptions. we will later
# output this as a descriptor/description list, but for the moment
# we will simply create each entry as a string, and insert it in
# a map keyed by the entry's name so that we can later output them
# in a way that looks sorted on the screen, rather than is sorted
# by the directory name (which looks pretty random to the human
# reader)
my %descriptions;
foreach my $gallery (@ARGV)
{
    my $gallery_underscore = $gallery;

    # Read the proper name of the program
    open ENTRYNAME, "<$gallery_dir/$gallery/doc/entry-name";
    my $entryname;
    while (my $line = <ENTRYNAME>) {
        chop $line;
        $entryname .= $line . " ";
    }
    chop $entryname;

    # Read the names of the authors, collate with commas, and at the
    # end chop the last comma off
    open AUTHOR, "<$gallery_dir/$gallery/doc/author";
    my $authors;
    while (my $line = <AUTHOR>) {
        chop $line;
        $authors .= $line . ", ";
    }
    chop $authors;
    $authors =~ s/,$//;

    $gallery_underscore    =~ s/-/_/;

    my $description;
    $description = "  <dt><b>\@ref code_gallery_${gallery_underscore} \"$entryname\"</b> (by $authors)</dt>\n";
    $description = $description . "    <dd>\n";
    open TOOLTIP, "<$gallery_dir/$gallery/doc/tooltip";
    while (my $line = <TOOLTIP>) {
        $description = $description . "      $line";
    }
    $description = $description . "    </dd>\n";
    $description = $description . "\n";

    # now insert this description into the map mentioned above
    $descriptions{$entryname} = $description;
}

# now print the entries generated above sorted by their keys
print "<dl>\n";
foreach my $key (sort keys %descriptions)
{
    print $descriptions{$key};
}
print "</dl>\n";


# Then print the rest of code-gallery.h.in
while (my $line = <GALLERY>)
{
  print $line;
}
close GALLERY;
