## ---------------------------------------------------------------------
##
## Copyright (C) 2013 - 2016 by the deal.II authors
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

if ($#ARGV < 1) {
  print "\nUsage: make_gallery.pl cmake_source_dir gallery_name gallery_dir gallery_src_files...\n";
  exit;
}

my $cmake_source_dir = shift(@ARGV);

my $gallery = shift(@ARGV);
my $gallery_underscore = $gallery;
$gallery_underscore    =~ s/-/_/;

my $gallery_dir = shift(@ARGV);

# next get the source files. sort all markdown files
# first so that we get to show them first. this makes
# sense because markdown files typically provide
# overview information
my @src_files = grep { $_ =~ m/.*\.(md|markdown)/ } @ARGV;
push @src_files, sort(grep { !($_ =~ m/.*\.(md|markdown)/) }@ARGV);


# read the names of authors; escape '<' and '>' as they
# appear in the email address. also trim trailing space and
# newlines
open AUTHORS, "<$gallery_dir/doc/author";
my $authors = <AUTHORS>;
$authors    =~ s/</&lt;/g; 
$authors    =~ s/>/&gt;/g; 
$authors    =~ s/\s*$//g;


# Read the proper name of the program
open ENTRYNAME, "<$gallery_dir/doc/entry-name";
my $entryname;
while (my $line = <ENTRYNAME>) {
    chop $line;
    $entryname .= $line . " ";
}
chop $entryname;


print
"/**
  * \@page code_gallery_$gallery_underscore The '$entryname' code gallery program
\@htmlonly
<p align=\"center\"> 
  This program was contributed by $authors.
  <br>
  It comes without any warranty or support by its authors or the authors of deal.II.
</p>

\@endhtmlonly

This program is part of the \@ref CodeGallery \"deal.II code gallery\" and
consists of the following files (click to inspect):
";

foreach my $file (@src_files)
{ 
  print "- <a href=\"../code-gallery/$gallery/$file\">$file</a>\n";
  if ($file =~ /.*\.(md|markdown|cc|cpp|cxx|c\+\+|h|hh|hxx)/) 
  {
      print "  (<a href=\"#ann-$file\">annotated version</a>)\n";
  }
}
print "\n";


# Next go through the list of files and see whether any of these are
# pictures we could show here:
my @picture_files;
foreach my $file (@src_files)
{ 
    if ($file =~ /.*\.(png|jpg|gif|svg)/)
    {
        push @picture_files, $file;
    }
}

if (@picture_files)
{
    print "<h1>Pictures from this code gallery program</h1>\n";
    print "<p align=\"center\">\n";
    print "<table>\n";

    # print four pictures per row
    while (@picture_files)
    {
        print "     <tr>\n";
        for my $i (0 .. 3)
        {
            if (@picture_files) 
            {
                print "       <td>\n";
                my $pic = pop(@picture_files);
                print "         <img width=\"250\" src=\"../code-gallery/$gallery/$pic\">\n";
                print "       </td>\n";
            }
        }
        print "     </tr>\n";
    }

    print "</table>\n";
    print "</p>\n";
}


# Then go through the list of files again and see which ones we can
# annotate and copy into the current document
foreach my $file (@src_files)
{ 
    # just copy markdown files as-is, but make sure we update links
    # that may be inlined. doxygen doesn't seem to understand the
    # ```...``` form of offset commands, so keep track of that as
    # well
    if ($file =~ /.*\.(md|markdown)/)
    {
        print "<a name=\"ann-$file\"></a>\n";
        print "<h1>Annotated version of $file</h1>\n";

        open MD, "<$gallery_dir/$file";
        my $incode = 0;
        while ($line = <MD>) 
        {
            # replace ``` markdown commands by doxygen equivalents
            while ($line =~ m/```/)
            {
                if ($incode == 0) {
                    $line =~ s/```/\@code{.sh}/;
                    $incode = 1;
                } else {
                    $line =~ s/```/\@endcode/;
                    $incode = 0;
                }
            }

            # update markdown links of the form "[text](./filename)"
            $line =~ s/(\[.*\])\(.\//\1\(..\/code-gallery\/$gallery\//g;
            print "$line";
        }

        print "\n\n";
    }

    # annotate source files
    if ($file =~ /.*\.(cc|cpp|cxx|c\+\+|h|hh|hxx)/)
    {
        print "<a name=\"ann-$file\"></a>\n";
        print "<h1>Annotated version of $file</h1>\n";

        system $^X, "$cmake_source_dir/doc/doxygen/scripts/program2doxygen", "$gallery_dir/$file";

        print "\n\n";
    }
}


# end the doxygen input file
print "*/\n";
