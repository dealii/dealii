## ---------------------------------------------------------------------
##
## Copyright (C) 2014 - 2018 by the deal.II authors
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

#
# A PERL script that modifies the default-generated 'footer.html' file
# that doxygen provides for us and customizes it for our needs.
#

use Sys::Hostname;
my $host = hostname;

my $hosting = << 'EOT'
&nbsp;&nbsp;Hosting provided by&nbsp;
<a href="http://www.iwr.uni-heidelberg.de/"><img src="https://www.dealii.org/pictures/IWRlogo4.png" alt="IWR"></a>
<a href="http://www.uni-heidelberg.de/"><img src="https://www.dealii.org/pictures/UniLogo4.png" alt="Universität Heidelberg"></a>
EOT
    ;

if ($host eq "simweb")
{
    s/\$doxygenversion/\$doxygenversion $hosting/;
}
