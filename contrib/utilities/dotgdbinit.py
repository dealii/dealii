## ------------------------------------------------------------------------
##
## SPDX-License-Identifier: LGPL-2.1-or-later
## Copyright (C) 2015 - 2022 by the deal.II authors
##
## This file is part of the deal.II library.
##
## Part of the source code is dual licensed under Apache-2.0 WITH
## LLVM-exception OR LGPL-2.1-or-later. Detailed license information
## governing the source code and code contributions can be found in
## LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
##
## ------------------------------------------------------------------------

#
# Instructions: Place a copy of this file, renamed as ".gdbinit", in your home
# directory to enable pretty-printing of various deal.II objects. If you already
# have a ".gdbinit" file or would like to manage multiple sets of pretty
# printers, then see the directions included in the Documentation, in the
# "Configuration for debugging via GDB" section in the "Information for users"
# category.
#
# This has only been tested with GDB 7.7.1 and newer, but it should work with
# slightly older versions of GDB (the Python interface was added in 7.0,
# released in 2009).
#
# Authors: Wolfgang Bangerth, 2015, David Wells, 2015 - 2018
#
set print pretty 1

python

import gdb
import re


def build_output_string(keys, accessor):
    """Build an output string of the form
    {
      a = foo,
      b = bar
    }
    where a and b are elements of keys and foo and bar are values of
    accessor (e.g., accessor['a'] = foo).

    Note that accessor need not be a dictionary (i.e., gdb.Values
    redefines __getitem__)."""
    return ("{\n" +
            ",\n".join(["  {} = {}".format(key, accessor[key])
                        for key in keys]) +
            "\n}")


class AlignedVectorPrinter(object):
    """Print a deal.II AlignedVector instance."""
    def __init__(self, typename, val):
        self.typename = typename
        self.val = val
        self.end = self.val['used_elements_end']
        # evaluate the get() method of the unique_ptr
        eval_string = "(*("+str(self.val['elements'].type)+"*)("+str(self.val['elements'].address)+")).get()"
        self.begin = gdb.parse_and_eval(eval_string);
        self.length = int(self.end - self.begin )

    def children(self):
        # The first entry (see the "Pretty Printing API" documentation of GDB)
        # in the tuple should be a name for the child, which should be nothing
        # (the array elements don't have individual names) here.
        return (("", (self.begin + count).dereference())
                for count in range(self.length))

    def to_string(self):
        return "AlignedVector<{}>({})".format(self.val.type.template_argument(0),
                                              self.length)

    @staticmethod
    def display_hint():
        """Provide a hint to GDB that this is an array-like container
        (so print values in sequence)."""
        return "array"


class PointPrinter(object):
    """Print a deal.II Point instance."""
    def __init__(self, typename, val):
        self.typename = typename
        self.val = val

    def to_string(self):
        return self.val['values']


class TensorPrinter(object):
    """Print a deal.II Tensor instance."""
    def __init__(self, typename, val):
        self.typename = typename
        self.val = val

    def to_string(self):
        if int(self.val.type.template_argument(0)) == 0:
            return self.val['value']
        else:
            return self.val['values']


class TriaIteratorPrinter(object):
    """Print a deal.II TriaIterator instance."""
    def __init__(self, typename, val):
        self.typename = typename
        self.val = val

    def to_string(self):
        keys = ['tria', 'present_level', 'present_index']
        if 'DoFHandler' in str(self.val.type.template_argument(0)):
            keys.insert(1, 'dof_handler')

        return build_output_string(keys, self.val['accessor'])


class VectorPrinter(object):
    """Print a deal.II Vector instance."""
    def __init__(self, typename, val):
        self.typename = typename
        self.val = val
        a_vec = self.val['values']
        self.end = a_vec['used_elements_end']
        # evaluate the elements.get() method of the AlignedVector member
        eval_string = "(*("+str(a_vec['elements'].type)+"*)("+str(a_vec['elements'].address)+")).get()"
        self.begin = gdb.parse_and_eval(eval_string);
        self.length = int(self.end - self.begin )

    def to_string(self):
        return ("Vector<{}>({})".format(self.val.type.template_argument(0),
                                        self.length) +
                build_output_string(['values'], self.val))


class QuadraturePrinter(object):
    """Print a deal.II Quadrature instance."""
    def __init__(self, typename, val):
        self.typename = typename
        self.val = val

    def to_string(self):
        return build_output_string(['quadrature_points', 'weights'], self.val)


class RxPrinter(object):
    """A "regular expression" printer which conforms to the
    "SubPrettyPrinter" protocol from gdb.printing."""
    def __init__(self, name, function):
        self.name = name
        self.function = function
        self.enabled = True

    def __call__(self, value):
        if self.enabled:
            return self.function(self.name, value)
        else:
            return None


class Printer(object):
    """A pretty-printer that conforms to the "PrettyPrinter" protocol
    from gdb.printing. It can also be used directly as an old-style
    printer."""
    def __init__(self, name):
        self.name = name
        self.subprinters = list()
        self.lookup = dict()
        self.enabled = True
        self.compiled_rx = re.compile('^([a-zA-Z0-9_:]+)<.*>$')

    def add(self, name, function):
        printer = RxPrinter(name, function)
        self.subprinters.append(printer)
        self.lookup[name] = printer

    @staticmethod
    def get_basic_type(object_type):
        # If it points to a reference, then get the reference.
        if object_type.code == gdb.TYPE_CODE_REF:
            object_type = object_type.target()

        object_type = object_type.unqualified().strip_typedefs()

        return object_type.tag

    def __call__(self, val):
        typename = self.get_basic_type(val.type)
        if typename:
            # All the types we match are template types, so we can use a
            # dictionary.
            match = self.compiled_rx.match(typename)
            if match:
                basename = match.group(1)
                if basename in self.lookup:
                    return self.lookup[basename](val)

        return None


dealii_printer = Printer("deal.II")


def register_dealii_printers():
    """Register deal.II pretty-printers with gdb."""
    printers = {
        AlignedVectorPrinter: ['AlignedVector'],
        PointPrinter: ['Point'],
        TensorPrinter: ['Tensor'],
        VectorPrinter: ['Vector'],
        TriaIteratorPrinter:
        ['TriaRawIterator', 'TriaIterator', 'TriaActiveIterator'],
        QuadraturePrinter:
        ['Quadrature', 'QGauss', 'QGaussLobatto', 'QMidpoint', 'QSimpson',
         'QTrapezoid',
         # The following name has been deprecated in deal.II 9.3 and can
         # be removed at a later time.
         'QTrapez',
         'QMilne', 'QWeddle', 'QGaussLog', 'QGaussLogR',
         'QGaussOneOverR', 'QSorted', 'QTelles', 'QGaussChebyshev',
         'QGaussRadauChebyshev', 'QIterated', 'QAnisotropic']
    }
    for printer, class_names in printers.items():
        for class_name in class_names:
            dealii_printer.add('dealii::' + class_name, printer)
    try:
        from gdb import printing
        printing.register_pretty_printer(gdb, dealii_printer)
    except ImportError:
        gdb.pretty_printers.append(dealii_printer)


register_dealii_printers()

end
