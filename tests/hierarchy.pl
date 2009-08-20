# Given the name of a configuration (i.e. system-name+compiler-name), this script
# generates a list of files that might contain the stored output of a test. The
# purpose is so that we don't have to store the output of each test for each
# configuration: if the output for a given configuration is equal to that of
# a related configuration, then we fall back to the latter, whereas we use
# the former if it is really different. This way we get a chain of
# possible files.
#
# For example, for a test named "check_me", if we run on a linux x86 system with
# compiler gcc3.4, we may have a file with output stored as
#   check_me/cmp/i686-pc-linux-gnu+gcc3.4
# However, if it is the same as what we would have gotten with gcc3.3, there would
# be no point in storing it and we should fall back to the output of gcc3.3. If
# that output would have been the same as that for a generic compiler (which we
# take as gcc3.2 on a x86 box), then the file
#   check_me/cmp/i686-pc-linux-gnu+gcc3.3
# may not exist either.
#
# When we want to compare generated output against stored output, we should therefore
# first compare against
#   check_me/cmp/i686-pc-linux-gnu+gcc3.4
# and if that doesn't exist against
#   check_me/cmp/i686-pc-linux-gnu+gcc3.3
# and if that doesn't exist again
#   check_me/cmp/generic
# (this latter file should always exist).
#
# This program generates this chain of files to check against, for any given input
# configuration.
#
# Author: Wolfgang Bangerth, 2005


$hierarchy{"i686-pc-linux-gnu+gcc3.3"}        = "generic";
$hierarchy{"i686-pc-linux-gnu+gcc3.4"}        = "i686-pc-linux-gnu+gcc3.3";
$hierarchy{"i686-pc-linux-gnu+gcc4.0"}        = "i686-pc-linux-gnu+gcc3.4";
$hierarchy{"i686-pc-linux-gnu+gcc4.1"}        = "i686-pc-linux-gnu+gcc4.0";
$hierarchy{"i686-pc-linux-gnu+gcc4.2"}        = "i686-pc-linux-gnu+gcc4.1";
$hierarchy{"i686-pc-linux-gnu+gcc4.3"}        = "i686-pc-linux-gnu+gcc4.2";
$hierarchy{"i686-pc-linux-gnu+gcc4.4"}        = "i686-pc-linux-gnu+gcc4.3";
$hierarchy{"i686-pc-linux-gnu+gcc4.5"}        = "i686-pc-linux-gnu+gcc4.4";

$hierarchy{"i686-pc-linux-gnu+icc7"}          = "generic";
$hierarchy{"i686-pc-linux-gnu+icc7.1"}        = "i686-pc-linux-gnu+icc7";
$hierarchy{"i686-pc-linux-gnu+icc8"}          = "i686-pc-linux-gnu+icc7.1";
$hierarchy{"i686-pc-linux-gnu+icc9"}          = "i686-pc-linux-gnu+icc8";
$hierarchy{"i686-pc-linux-gnu+icc10"}         = "i686-pc-linux-gnu+icc9";

$hierarchy{"mips-sgi-irix6.5+MIPSpro7.4"}     = "generic";

$hierarchy{"x86_64-unknown-linux-gnu+gcc3.3"} = "generic";
$hierarchy{"x86_64-unknown-linux-gnu+gcc3.4"} = "x86_64-unknown-linux-gnu+gcc3.3";
$hierarchy{"x86_64-unknown-linux-gnu+gcc4.0"} = "x86_64-unknown-linux-gnu+gcc3.4";
$hierarchy{"x86_64-unknown-linux-gnu+gcc4.1"} = "x86_64-unknown-linux-gnu+gcc4.0";
$hierarchy{"x86_64-unknown-linux-gnu+gcc4.2"} = "x86_64-unknown-linux-gnu+gcc4.1";
$hierarchy{"x86_64-unknown-linux-gnu+gcc4.3"} = "x86_64-unknown-linux-gnu+gcc4.2";
$hierarchy{"x86_64-unknown-linux-gnu+gcc4.4"} = "x86_64-unknown-linux-gnu+gcc4.3";

$hierarchy{"x86_64-unknown-linux-gnu+icc9"}   = "x86_64-unknown-linux-gnu+gcc3.3";
$hierarchy{"x86_64-unknown-linux-gnu+icc10"}  = "x86_64-unknown-linux-gnu+icc9";

$hierarchy{"powerpc-apple-darwin8.8.0+gcc4.0"}  = "generic";
$hierarchy{"powerpc-apple-darwin8.10.0+gcc4.0"} = "powerpc-apple-darwin8.8.0+gcc4.0";

$hierarchy{"ia64-unknown-linux-gnu+icc9"}     = "x86_64-unknown-linux-gnu+icc9";
$hierarchy{"ia64-unknown-linux-gnu+icc10"}    = "ia64-unknown-linux-gnu+icc9";

# derive mac os x from mac os x because they use the same random number
# generator
$hierarchy{"i386-apple-darwin8.10.1+gcc4.0"}  = "powerpc-apple-darwin8.10.0+gcc4.0";


$configuration = $ARGV[0];

# first check whether the given configuration is known at all. if not,
# we don't know anything about the hierarchy we have to search, but we should
# certainly search under the name of this configuration as well as the generic
# one as well
if (! defined $hierarchy{$configuration}) {
    print "$configuration generic\n";
    exit;
}

# so the configuration is known. now output a list of files that we
# will have to check for results:
$name_list = "$configuration";
while (defined $hierarchy{$configuration}) {
    $name_list = "$name_list $hierarchy{$configuration}";
    $configuration = $hierarchy{$configuration};
}
print "$name_list\n";

