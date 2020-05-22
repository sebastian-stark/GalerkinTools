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

#ifndef GALERKINTOOLS_LDR_H_
#define GALERKINTOOLS_LDR_H_

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>
#include <vector>

#include <galerkin_tools/config.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

namespace Auxiliary
{

/**
 * A method decomposing a matrix according to \f$\boldsymbol{C}=\boldsymbol{L}\boldsymbol{D}\boldsymbol{R}^\top\f$ such that
 * \f$\boldsymbol{R}=\boldsymbol{L}\f$ for symmetric matrices \f$\boldsymbol{C}\f$, and \f$\boldsymbol{D}\f$ diagonal with entries either @p -1 or @p 1.
 *
 * Internally, the method is based on a singular value decomposition as provided by the LAPACK library
 * (to interface to LAPACK, the LAPACKE function @p LAPACKE_dgesvd is used).
 *
 * @param[in]		C		\f$\boldsymbol{C}\f$
 *
 * @param[out]		D		vector with diagonal elements of \f$\boldsymbol{D}\f$ (elements are either @p -1 or @p 1)
 *
 * @param[out]		L		@p std::vector with rows of \f$\boldsymbol{L}\f$
 *
 * @param[out]		R		@p std::vector with rows of \f$\boldsymbol{R}\f$
 */
int
compute_ldr(	FullMatrix<double>& 			C,
				Vector<double>&					D,
				std::vector<Vector<double>>&	L,
				std::vector<Vector<double>>&	R);
}

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_LDR_H_ */
