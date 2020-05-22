// --------------------------------------------------------------------------
// Copyright (C) 2020 by Sebastian Stark
//
// This file is part of the GalerkinTools library
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 2
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

The main purpose of the GalerkinTools library is to help with the
assembly of finite element systems for formulations obtained based on the Galerkin approach. Though the library must
of course make certain restrictions on the class of problems for which it can be applied, it is believed to
be useful for the solution of many practical problems especially in the field of solid mechanics. Special focus of
the library is on problems involving different domain portions (e.g. different materials) as well as on problems
involving unknown fields living on interfaces or boundaries.

The library currently requires deal.II 9.2, lapacke and numdiff to be installed on your system (and cmake must be able to find them).
To use the direct sequential solver interface of the library, deal.II must additionally be configured with UMFPACK.
To use the distributed parallel capabilities of the library, deal.II must be configured with MPI, PETSc (including MUMPS for the parallel solver) and p4est.

Installation of the library is through cmake:

(1) place library source files into some folder /path/to/folder/GalerkinTools (you can use git clone https://github.com/starki0815/GalerkinTools.git for this)
(2) cd /path/to/folder/
(3) mkdir build
(4) cd build
(5) cmake -DCMAKE_INSTALL_PREFIX=/path/to/install/dir ../GalerkinTools
(6) make install
(7) optionally set an environment variable to GALERKIN_TOOLS_DIR=/path/to/install/dir
(8) optionally run the tests (first cd /path/to/folder/build, then ctest)
