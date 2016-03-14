#!/usr/bin/env python

## ---------------------------------------------------------------------
##
## Copyright (C) 2016 by the deal.II authors
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

## This script runs install_name_tool on each library under basedir (first argument
## on command line) and makes sure that all libraries are installed with absolute
## id name. This makes Mac OS X happy about not using DYLD_LIBRARY_PATH, which is generally
## a bad thing. 

import os
import sys
from subprocess import check_output as run

basedir = '/Applications/deal.II.app/Contents/Resources/opt/'
if len(sys.argv) > 1:
    basedir = sys.argv[1]

## Lib object
# 
# Stores the name of a library, its install name, its location and all its dependencies, as given by `otool -L`

class Lib(object):
    location = None
    install_name = None
    name = None
    deps = None
    def __str__(self):
        mystr = ""
        mystr +=  "Name        :"+str(self.name)
        mystr +="\nInstall name:"+str(self.install_name)
        mystr +="\nLocation    :"+str(self.location)
        mystr +="\nDeps        : ... "
        return mystr

# Walk the current tree, and extract all libraries
def get_libs():
    libraries = []

    for root, dirs, files in os.walk(basedir):
        for f in files:
            filename = os.path.join(root, f)
            if os.path.isfile(filename) and f.endswith("dylib"):
                a = Lib()
                a.name = f
                a.install_name = run(["otool", "-D", filename]).split('\n')[1]
                a.location= os.path.join(root, f)
                long_deps = run(["otool", "-L", a.location]).split('\n')
                a.deps = [dep.partition(' ')[0][1::] for dep in long_deps[2:-1]]
                libraries += [a]
    return libraries


# # Fix all install names first
# 
# Some will fail, because they are either stub files, or system files...

# In[ ]:

libraries = get_libs()
c = 0
failed = []
for c in range(len(libraries)):
    i = libraries[c]
    if os.path.islink(i.location):
        continue
    if i.install_name != i.location:
        try:
            run(["install_name_tool",'-id',i.location, i.location])
        except:
            print("Failed: ",i.name)
            print( "(",libraries[c].name,")")
            failed += [c]
print (failed)
print ('Removing failed libs...')
for fail in reversed(failed):
    print ("Removing from list",libraries[fail].name)
    del libraries[fail]


# # Fix all dependencies with absolute paths

for c  in xrange(len(libraries)):
    this_lib = libraries[c]
    if this_lib.install_name != this_lib.location:
        print('Not valid:', i.name)
    else:
        lib = this_lib.name
        command = ['install_name_tool']
        for dep in this_lib.deps:
            for loc in libraries:
                if dep.find(loc.name) != -1 and dep != loc.install_name:
                    command += ['-change', dep, loc.install_name]
                    break
        try:
            if len(command) != 1:
                command += [this_lib.location]
                print('Processing', lib)
                print('======================\n\n')
                print(" ".join(command))
                print('======================\n\n')
                run(command)
        except:
            print('\n\n*********************************************')
            print('Last command failed!')
            print(" ".join(command))
            print('*********************************************\n\n')




