$/='this should not appear in any of these stupid files';

$fn = $ARGV;
$guard = '__deal2__' . $fn;

$guard =~ s/\./_/g;

$head = << 'EOT'
//----------------------------  XXX  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  XXX  ---------------------------
#ifndef GGG
#define GGG


EOT
    ;

$head =~ s/XXX/$fn/g;
$head =~ s/GGG/$guard/g;

s/^.+\#define __[^\n]*/$head/s;
s/\/\*-+\s+$fn\s+-+\*\/\n//g;
s/\/\*\s*end of.*H\s*\*\///g;
s/\n[ \t]*\n[ \t]*\n\s*\n*/\n\n\n/g;
