# Welcome to cellular-automaton

It contains the C++ based library CAM implementing a simple [cellular automaton method (CAM)](
https://en.wikipedia.org/wiki/Cellular_automaton). In this CAM, building units containing several cells/pixels are allowed to
move within a [von Neumann neighborhood (VNN)](
https://en.wikipedia.org/wiki/Von_Neumann_neighborhood) of range `jump_parameter` and larger
particles (composites), which are edge-connected sets of pixels, are allowed to move in a VNN depending on their
size, i.e., the number of the edge-connected pixels. The range of the VNN of larger particles can be
calculated as `jump_parameter` over the square-root of the larger particle's size. All particles,
i.e., single pixels and larger particles, move in such a way that the amount of their particle
neighbors is maximized.
Furthermore the following optional properties of a building unit and CAM can be defined:
  - Different shapes of the building unit (pre-implemented are spheres, planes/cubes and random particles based on spheres.)
  - Boundary face weights of the building unit (attraction and repulsion effects between building units or composites can be modelled)
  - Rotation of building units and composites (around the center or user-defined e.g. edge points of basal axis)
  - Stability of composites (The strength of a connection between two parts of a composite (sub-aggregates) is currently defined as number of contact faces between them divided by the size of each sub-aggregate. If the connection is considered as weak by both (below a certain threshold), they can move independently of each other)

# How to use cellular-automaton / CAM

There is two ways of using CAM:

- In `Python`, one has to import the module as illustrated in the Python examples. One has to run
  `cmake` and obtain all submodules, or the `setup.sh` script to use the Python version of the code.
  In any case, the selected compiler needs to be compatible with C++20.
- In `C++`, one has to include the file `include/CAM/cellular_automaton.hxx`. Then, you can run the
  CAM as shown in `cpp_example.cxx` provided that you compile it using `-std=gnu++20`. A possible
  compilation command for the `cpp_example.cxx` is `clang++-12 -std=gnu++20 -Wall -Wextra -pedantic 
  -Iinclude -O3 examples_cpp/cpp_example.cxx -o test`. Here, `clang++-12` can be replaced by any
  suitable compiler implementing C++20.


# Copyright, License, and Contribution Policy

This directory contains the CAM library.

The CAM library is copyrighted by the authors of `cellular-automaton`. This term refers to the
people listed at the very top of this page.

The CAM library is free software; you can use it, redistribute it, and/or modify it under the terms
of the <b>GNU Lesser General Public License</b> as published by the Free Software Foundation; either
<b>version 2.1</b> of the License, or (at your option) any later version. The full text of the GNU
Lesser General Public version 2.1 is quoted in [License.txt](License.txt).


## Contributions

As a contributor to this project, you agree that all of your contributions be governed by the
<b>Developer Certificate of Origin version 1.1</b>. The CAM project does not require copyright
assignments for contributions. This means that the copyright for code contributions in the CAM
project is held by its respective contributors who have each agreed to release their contributed
code under a compatible open source license (LGPL v2.1 for library code). The full text of the 
Developer Certificate of Origin version 1.1 is quoted in [DeveloperCertificateOfOrigin.txt](
DeveloperCertificateOfOrigin.txt).


## Referencing the library

In addition to the terms imposed by the LGPL v2.1 or later, we ask for the following courtesy:

> Every publication presenting numerical results obtained with the help of CAM should state the name
> of the library and cite one or more of the following references  
> - A. Kazarnikov, N. Ray, H. Haario, J. Lappalainen, and A. Rupp  
>   ***Parameter estimation for cellular automata***  
>   arXiv preprint, doi: [10.48550/arXiv.2301.13320](https://doi.org/10.48550/arXiv.2301.13320)
> - A. Rupp, M. Gahn, and G. Kanschat  
> ***Partial differential equations on hypergraphs and networks of surfaces: Derivation and hybrid
  discretizations***  
> ESAIM: Mathematical Modelling and Numerical Analysis, doi: [10.1051/m2an/2022011](
  https://doi.org/10.1051/m2an/2022011)

This is the usual, fair way of giving credit to contributors to a scientific result. In addition, it
helps us justify our effort in developing CAM as an academic undertaking. The last publication
refers to [HyperHDG](https://github.com/HyperHDG/), from which we borrowed the just-in-time C++
compilation.


## Contact

For further questions regarding licensing and commercial use please contact Andreas Rupp directly
using [Email](mailto:info@rupp.ink).


## Links

- The license can be found in [License.txt](License.txt). It contains the [GNU Lesser General Public
License version 2.1](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html).
- The developer certificate of origin can be found in 
[DeveloperCertificateOfOrigin.txt](DeveloperCertificateOfOrigin.txt). It contains the [Developer 
Certificate of Origin version 1.1](https://developercertificate.org/).
