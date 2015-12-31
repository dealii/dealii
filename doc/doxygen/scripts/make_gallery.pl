## ---------------------------------------------------------------------
##
## Copyright (C) 2013, 2015 by the deal.II authors
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

if ($#ARGV < 1) {
  print "\nUsage: make_gallery.pl cmake_source_dir gallery_name gallery_dir gallery_src_files...\n";
  exit;
}

my $cmake_source_dir = shift(@ARGV);

my $gallery = shift(@ARGV);
my $gallery_underscore = $gallery;
$gallery_underscore    =~ s/-/_/;

my $gallery_dir = shift(@ARGV);
my $author_file = "$gallery_dir/doc/author";

my @src_files = @ARGV;

# read the names of authors; escape '<' and '>' as they
# appear in the email address. also trim trailing space and
# newlines
open AUTHORS, "<$author_file";
my $authors = <AUTHORS>;
$authors    =~ s/</&lt;/g; 
$authors    =~ s/>/&gt;/g; 
$authors    =~ s/\s*$//g;

print
"/**
  * \@page code_gallery_$gallery_underscore The $gallery code gallery program
\@htmlonly
<p align=\"center\"> 
  This program was contributed by $authors.
  <br>
  It comes without any warranty or support by its authors or the authors of deal.II.
</p>

\@endhtmlonly

This program consists of the following files (click to inspect):
";

foreach my $file (@src_files)
{ 
  print "- <a href=\"../code-gallery/$gallery/$file\">$file</a>\n";
}
print "\n";

print "*/\n";
