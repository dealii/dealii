This folder contains changelog entries
======================================

Changes between different releases in deal.II are documented in the file 
"changes.h" that is automatically created from the files in the directories

  * "major", 
  * "minor" and 
  * "incompatibilities".

Each of these files contain a short description of a change, the authors' names
and the date where the latter two a separated by "<br>" from the first one. 
A typical file could look like this:

    New: We can do fancy stuff now.
    <br>
    (John Doe, YYYY/MM/DD)

and is named "YYYYMMDDJohnDoe". File names for multiple contributions from the 
same author on one day have a number appended, e.g. "YYYMMDDJohnDoe\_1".

incompatibilities
-----------------

The "incompatibilities" section contains modifications to the library that 
are incompatible with previous versions of the library, but which are
necessary for the future maintainability of the library. 

major
-----

Entries in the "major" section are reserved for substantial changes to the
library. In particular, this should be sufficiently compact. A good example
for this category would be adding a new tutorial program.

minor
-----

Finally, the "minor" section contains all the changes that do not fit in
the former two.
