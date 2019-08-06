CoalRe
======

BEAST 2 package for inference under the coalescent with reassortment,
applicable to segmented viral genomes.

This repository contains the source code for CoalRe.  It is mostly
of interest to phylogenetic methods developers.  If you are interested
in _using_ CoalRe, please visit the [Taming the BEAST tutorial page](https://taming-the-beast.org/tutorials/Reassortment-Tutorial/).

Building CoalRe
---------------

In order to build CoalRe from the source, you will need the following:

1. [OpenJDK](https://adoptopenjdk.net) v8 or later,
2. The Apache Ant build tool.

Once these are installed, open a shell in the root directory of this repository
and use

    $ ant

to build the package.

License
-------

CoalRe is free (as in freedom) software and is distributed under the terms of
version 3 of the GNU General Public License.  A copy of this license is found
in the file named `COPYING`.
