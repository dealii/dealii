#----------------------------  timer.cc  ---------------------------
#    $Id$
#    Version: $Name$
#
#    Copyright (C) 2003 by the deal.II authors
#
#    This file is subject to QPL and may not be  distributed
#    without copyright and license information. Please refer
#    to the file deal.II/doc/license.html for the  text  and
#    further information on this license.
#
#----------------------------  timer.cc  ---------------------------

# This is the script with which threads.cc is generated
#
# The idea of the generated file is: have classes X<N> which
# are not convertible into each other, and have functions with
# 0..10 arguments of types X<1>..X<10>, with reference and/or
# constant arguments, being virtual, static, or const. This way,
# we check all possible spawn(...) combinations and all possible
# combinations of arguments. Since the argument types are not
# convertible into each other, this also makes sure that there are
# no typos like "Arg4, Arg4" instead of "Arg4, Arg5", that might
# otherwise sneak in.

$N=10;

print "#include <base/thread_management.h>\n";
print "#include <base/logstream.h>\n";
print "#include <fstream>\n";
print "#include <iostream>\n";
print "template <int> struct X {};\n";
print "struct U {\n";


sub print_args {
    my $i = $_[0];
    my $pre = $_[1];
    my $post = $_[2];
    my $end = $_[3];
    for (my $j=1; $j<=$i; ++$j) {
	print $pre, $j, $post, ($j!=$i ? "," : "");
    }
    print ")${end} { \n";
    print "    deallog << __PRETTY_FUNCTION__ << std::endl;\n";
    print "    static X<0> x; return x;\n";
    print "  };\n";
}

for ($i=0; $i<=$N; ++$i) {
    print "  X<0> foo_$i (";
    print_args ($i, "X<", ">", "");

    print "  static X<0> static_foo_$i (";
    print_args ($i, "X<", ">", "");

    print "  X<0> & ref_foo_$i (";
    print_args ($i, "X<", ">", "");

    print "  static X<0> & static_ref_foo_$i (";
    print_args ($i, "X<", ">", "");

    print "  static const X<0> & static_const_ref_foo_$i (";
    print_args ($i, "X<", ">", "");

    print "  X<0> foo_ref_$i (";
    print_args ($i, "X<", ">&", "");

    print "  static X<0> static_foo_ref_$i (";
    print_args ($i, "X<", ">&", "");

    print "  X<0> & ref_foo_ref_$i (";
    print_args ($i, "X<", ">&", "");

    print "  const X<0> & const_ref_foo_$i (";
    print_args ($i, "X<", ">", "");

    print "  const X<0> & const_ref_foo_ref_$i (";
    print_args ($i, "X<", ">&", "");

    print "  const X<0> & const_ref_foo_const_ref_$i (";
    print_args ($i, "X<", ">&", "");

    print "  const X<0> & const_ref_foo_${i}_const (";
    print_args ($i, "X<", ">", "const");

    print "  const X<0> & const_ref_foo_ref_${i}_const (";
    print_args ($i, "X<", ">&", "const");

    print "  const X<0> & const_ref_foo_const_ref_${i}_const (";
    print_args ($i, "X<", ">&", "const");

    print "  virtual const X<0> & virtual_const_ref_foo_$i (";
    print_args ($i, "X<", ">", "");

    print "  virtual const X<0> & virtual_const_ref_foo_ref_$i (";
    print_args ($i, "X<", ">&", "");

    print "  virtual const X<0> & virtual_const_ref_foo_const_ref_$i (";
    print_args ($i, "X<", ">&", "");

    print "  virtual const X<0> & virtual_const_ref_foo_${i}_const (";
    print_args ($i, "X<", ">", "const");

    print "  virtual const X<0> & virtual_const_ref_foo_ref_${i}_const (";
    print_args ($i, "X<", ">&", "const");

    print "  virtual const X<0> & virtual_const_ref_foo_const_ref_${i}_const (";
    print_args ($i, "X<", ">&", "const");

    print "  static X<0> & static_ref_foo_ref_$i (";
    print_args ($i, "X<", ">&", "");

    print "  static const X<0> & static_const_ref_foo_ref_$i (";
    print_args ($i, "X<", ">&", "");

    print "  X<0> foo_const_ref_$i (";
    print_args ($i, "const X<", ">&", "");

    print "  static X<0> static_foo_const_ref_$i (";
    print_args ($i, "const X<", ">&", "");

    print "  X<0> & ref_foo_const_ref_$i (";
    print_args ($i, "const X<", ">&", "");

    print "  static X<0> & static_ref_foo_const_ref_$i (";
    print_args ($i, "const X<", ">&", "");

    print "  static const X<0> & static_const_ref_foo_const_ref_$i (";
    print_args ($i, "const X<", ">&", "");

    print "  X<0> foo_${i}_const (";
    print_args ($i, "X<", ">", " const");

    print "  X<0> & ref_foo_${i}_const (";
    print_args ($i, "X<", ">", " const");

    print "  X<0> foo_ref_${i}_const (";
    print_args ($i, "X<", ">&", " const");

    print "  X<0> & ref_foo_ref_${i}_const (";
    print_args ($i, "X<", ">&", " const");

    print "  X<0> foo_const_ref_${i}_const (";
    print_args ($i, "const X<", ">&", " const");

    print "  X<0> & ref_foo_const_ref_${i}_const (";
    print_args ($i, "const X<", ">&", " const");

    print "  virtual X<0> virtual_foo_${i}_const (";
    print_args ($i, "X<", ">", " const");

    print "  virtual X<0> & virtual_ref_foo_${i}_const (";
    print_args ($i, "X<", ">", " const");

    print "  virtual X<0> virtual_foo_ref_${i}_const (";
    print_args ($i, "X<", ">&", " const");

    print "  virtual X<0> & virtual_ref_foo_ref_${i}_const (";
    print_args ($i, "X<", ">&", " const");

    print "  virtual X<0> virtual_foo_const_ref_${i}_const (";
    print_args ($i, "const X<", ">&", " const");

    print "  virtual X<0> & virtual_ref_foo_const_ref_${i}_const (";
    print_args ($i, "const X<", ">&", " const");

    print "  virtual X<0> virtual_foo_$i (";
    print_args ($i, "X<", ">", "");

    print "  virtual X<0> & virtual_ref_foo_$i (";
    print_args ($i, "X<", ">", "");

    print "  virtual X<0> virtual_foo_ref_$i (";
    print_args ($i, "X<", ">&", "");

    print "  virtual X<0> & virtual_ref_foo_ref_$i (";
    print_args ($i, "X<", ">&", "");

    print "  virtual X<0> virtual_foo_const_ref_$i (";
    print_args ($i, "const X<", ">&", "");

    print "  virtual X<0> & virtual_ref_foo_const_ref_$i (";
    print_args ($i, "const X<", ">&", "");
}

print "};\n";


print "int main () {\n";
print "  std::ofstream logfile(\"threads.output\");\n";
print "  deallog.attach(logfile);\n";
print "  deallog.depth_console(0);\n";

print "  using namespace Threads;\n";
print "  ThreadGroup<X<0> > tg;\n";
print "  ThreadGroup<X<0>&> tgr;\n";
print "  ThreadGroup<const X<0>&> tgcr;\n";
print "  U u;\n";
for ($i=1; $i<=$N; ++$i) {
    print "X<$i> x$i;\n";
}
for ($i=0; $i<=$N; ++$i) {
    $arglist = "(";
    for ($j=1; $j<=$i; ++$j) {
	$arglist = $arglist . "x$j" . ($j!=$i ? "," : "");
    }
    $arglist = $arglist . ")";

#### where are the const ref functions??
    print << "END"
    tgr += spawn (u, &U::ref_foo_${i}) $arglist;
    tgr += spawn (u, &U::ref_foo_${i}_const) $arglist;
    tgr += spawn (u, &U::ref_foo_const_ref_${i}) $arglist;
    tgr += spawn (u, &U::ref_foo_const_ref_${i}_const) $arglist;
    tgr += spawn (u, &U::ref_foo_ref_${i}) $arglist;
    tgr += spawn (u, &U::ref_foo_ref_${i}_const) $arglist;
    tgcr += spawn (u, &U::const_ref_foo_${i}) $arglist;
    tgcr += spawn (u, &U::const_ref_foo_${i}_const) $arglist;
    tgcr += spawn (u, &U::const_ref_foo_const_ref_${i}) $arglist;
    tgcr += spawn (u, &U::const_ref_foo_const_ref_${i}_const) $arglist;
    tgcr += spawn (u, &U::const_ref_foo_ref_${i}) $arglist;
    tgcr += spawn (u, &U::const_ref_foo_ref_${i}_const) $arglist;
    tgcr += spawn (u, &U::virtual_const_ref_foo_${i}) $arglist;
    tgcr += spawn (u, &U::virtual_const_ref_foo_${i}_const) $arglist;
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_${i}) $arglist;
    tgcr += spawn (u, &U::virtual_const_ref_foo_const_ref_${i}_const) $arglist;
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_${i}) $arglist;
    tgcr += spawn (u, &U::virtual_const_ref_foo_ref_${i}_const) $arglist;
    tg += spawn (u, &U::foo_${i}) $arglist;
    tg += spawn (u, &U::foo_${i}_const) $arglist;
    tg += spawn (u, &U::foo_const_ref_${i}) $arglist;
    tg += spawn (u, &U::foo_const_ref_${i}_const) $arglist;
    tg += spawn (u, &U::foo_ref_${i}) $arglist;
    tg += spawn (u, &U::foo_ref_${i}_const) $arglist;
    tgr += spawn (u, &U::virtual_ref_foo_${i}) $arglist;
    tgr += spawn (u, &U::virtual_ref_foo_${i}_const) $arglist;
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_${i}) $arglist;
    tgr += spawn (u, &U::virtual_ref_foo_const_ref_${i}_const) $arglist;
    tgr += spawn (u, &U::virtual_ref_foo_ref_${i}) $arglist;
    tgr += spawn (u, &U::virtual_ref_foo_ref_${i}_const) $arglist;
    tg += spawn (u, &U::virtual_foo_${i}) $arglist;
    tg += spawn (u, &U::virtual_foo_${i}_const) $arglist;
    tg += spawn (u, &U::virtual_foo_const_ref_${i}) $arglist;
    tg += spawn (u, &U::virtual_foo_const_ref_${i}_const) $arglist;
    tg += spawn (u, &U::virtual_foo_ref_${i}) $arglist;
    tg += spawn (u, &U::virtual_foo_ref_${i}_const) $arglist;

    tgr += spawn (&U::static_ref_foo_${i}) $arglist;
    tgr += spawn (&U::static_ref_foo_const_ref_${i}) $arglist;
    tgr += spawn (&U::static_ref_foo_ref_${i}) $arglist;
    tgcr += spawn (&U::static_const_ref_foo_${i}) $arglist;
    tgcr += spawn (&U::static_const_ref_foo_const_ref_${i}) $arglist;
    tgcr += spawn (&U::static_const_ref_foo_ref_${i}) $arglist;
    tg += spawn (&U::static_foo_${i}) $arglist;
    tg += spawn (&U::static_foo_const_ref_${i}) $arglist;
    tg += spawn (&U::static_foo_ref_${i}) $arglist;
END
    ;
}
print "  tg.join_all();\n";
print "  tgr.join_all();\n";
print "}\n";

#############
