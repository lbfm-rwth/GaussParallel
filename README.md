[![Build Status](https://github.com/lbfm-rwth/GaussParallel/workflows/CI/badge.svg?branch=master)](https://github.com/lbfm-rwth/GaussParallel/actions?query=workflow%3ACI+branch%3Amaster)
## The GaussParallel package
This package provides an implementation of a parallel Gaussian elimination algorithm as described in
[A parallel algorithm for Gaussian elemination over finite fields](https://arxiv.org/abs/1806.04211).

Please note that for small matrices the logistical overhead might be
significant and it might be more efficient to use the standard methods for
Gaussian elimination from the [GAUSS package](https://www.gap-system.org/Packages/gauss.html).

## Compatibility with old HPC-GAP version
Older versions of hpcgap still have working guards. Unfortunately, in those
versions you can't load this package via `LoadPackage`.
Let's assume that you can run your old HPC-GAP version via $HPCGAP-OLD.
Open a shell in the root folder of the `GaussParallel` package and run:
`$HPCGAP-OLD compatibility-for-old-hpcgap/read.g`

## Contact
Please submit bug reports, suggestions for improvements and patches via
the [issue tracker](https://github.com/lbfm-rwth/GaussParallel/issues)
or via email to
[Sergio Siccha](mailto:siccha@mathematik.uni-kl.de)
or
[Jendrik Brachter](mailto:brachter@cs.uni-kl.de).

## License
`GaussParallel` is free software you can redistribute it and/or modify it
under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 2 of the License, or (at your option) any
later version.
See file gpl-2.0.txt for further information.
