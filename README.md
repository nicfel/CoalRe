CoalRe
======

[![Build status](https://github.com/nicfel/CoalRe/actions/workflows/ci-publish.yml/badge.svg)](https://github.com/nicfel/CoalRe/actions/workflows/ci-publish.yml)


BEAST 2 package for inference under the coalescent with reassortment,
applicable to segmented viral genomes.

This repository contains the source code for CoalRe.  It is mostly
of interest to phylogenetic methods developers.  If you are interested
in _using_ CoalRe, please visit the [Taming the BEAST tutorial page](https://taming-the-beast.org/tutorials/Reassortment-Tutorial/).

Building CoalRe
---------------

In order to build CoalRe from the source, you will need the following:

1. [OpenJDK](https://adoptium.net) 25 or later,
2. The Apache Maven build tool.

Once these are installed, open a shell in the root directory of this repository
and use

    $ mvn package

to run the tests and build the BEAST package archive, which is written to
`target/CoalRe.v<version>.zip`.

License
-------

CoalRe is free (as in freedom) software and is distributed under the terms of
version 3 of the GNU General Public License.  A copy of this license is found
in the file named `COPYING`.
