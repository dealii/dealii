## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2019 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

# This script sets the canonical webpage for all the webpages of the
# doxygen documentation. For more information see
# https://en.wikipedia.org/wiki/Canonical_link_element
#
# This script is invoked by CMake.

use strict;
use warnings;

use File::Find qw(find);

my @html;
my $link_start = '<link rel="canonical" href="https://www.dealii.org/current/doxygen/';

find(sub {
	  if ($_ =~ /\.html$/)
	  {
	       push(@html, $File::Find::name);
	  }
     }, "./");
chomp(@html);

for my $filename (@html) {
     if ($filename eq "./header.html")
     {
	  # there is no relevant content in the header
	  next;
     }

     open(my $file_handle, "+<:encoding(UTF-8)", $filename)
	  or die "Could not read $filename: $!";

     my $match = 0;
     while (my $line = <$file_handle>)
     {
	  if (index($line, '<link rel="canonical"') != -1)
	  {
	       $match = 1;
	       last;
	  }
     }

     my $new_text_data = "";
     if ( ! $match )
     {
	  seek($file_handle, 0, 0)
	       or die "Could not seek $filename: $!";

	  my $canonical_link_added = 0;
	  while (my $line = <$file_handle>)
	  {
	       $new_text_data .= $line;

	       # Do not add canonical link twice
	       if ( ! $canonical_link_added and index($line, '<head>') != -1)
	       {
		    die unless substr($filename, 0, 2) eq "./";

		    $new_text_data .= $link_start . substr($filename, 2) . '" />' . "\n";
		    $canonical_link_added = 1;
	       }
	  }
     }

     if (length($new_text_data) > 0)
     {
     	  # Truncate the file and write the new text
     	  seek($file_handle, 0, 0)
     	       or die "Could not seek $filename: $!";
	  print $file_handle $new_text_data;

     	  truncate($file_handle, tell($file_handle))
     	       or die "Could not truncate $filename: $!";
     	  close($file_handle)
     	       or die "Could not close $filename: $!";
     }
}
