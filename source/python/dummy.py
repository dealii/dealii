# ---------------------------------------------------------------------
#
# Copyright (C) 2016 by the deal.II authors
#
# This file is part of the deal.II library.
#
# The deal.II library is free software; you can use it, redistribute
# it, and/or modify it under the terms of the GNU Lesser General
# Public License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# The full text of the license can be found in the file LICENSE at
# the top level of the deal.II distribution.
#
# ---------------------------------------------------------------------

"""This a dummy module that contains a global function and a class."""

__all__ = ['Dummy', 'bar']

def bar():
    """Print bar."""
    print('bar')

class Dummy:
    """This a dummy class that has only one member function."""
    
    def foo(self):
        """Print foo."""
        print('foo')
