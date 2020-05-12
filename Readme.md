// --------------------------------------------------------------------------
// Copyright (C) 2019 by Sebastian Stark
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

Minimum requirements: deal.II, lapacke, numdiff
For direct sequential solver: deal.II configured with UMFPACK
For distributed parallel computations: deal.II configured with MPI, PETSc, p4est
