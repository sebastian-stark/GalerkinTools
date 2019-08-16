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

#ifndef GALERKINTOOLS_TOTALPOTENTIAL_H_
#define GALERKINTOOLS_TOTALPOTENTIAL_H_

#include <vector>

#include <galerkin_tools/config.h>
#include <galerkin_tools/total_potential_contribution.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * Class collecting all TotalPotentialContribution objects \f$\Pi_i\f$ in the total potential \f$\Pi\f$
 * according to
 * \f{equation*}
 * \Pi = \Pi(H^\Omega_\rho, H^\Sigma_\tau, C_\iota) = \sum_i \Pi_i(H^\Omega_\rho, H^\Sigma_\tau, C_\iota).
 * \f}
 *
 * @tparam	spacedim	Spatial dimension
 */
template<unsigned int spacedim>
class TotalPotential
{
	private:

		/**
		 * Collection of the TotalPotentialContribution objects corresponding to the \f$\Pi_i\f$
		 */
		std::vector< SmartPointer<const TotalPotentialContribution<spacedim>> >
		total_potential_contributions;

		/**
		 * Maximum number of dependent variables in a ScalarFunctional (or ScalarFunctional<spacedim, spacedim>, respectively).
		 * This information is used internally to reserve an appropriate amount of memory for
		 * certain vectors and matrices.
		 */
		unsigned int
		max_dependent_vars = 0;

		/**
		 * AssemblyHelper is a friend of TotalPotential to get direct access to
		 * TotalPotential::total_potential_contributions and
		 * TotalPotential::max_dependent_vars
		 */
		template<unsigned int> friend class AssemblyHelper;

	public:

		/**
		 * Add a contribution \f$\Pi_i\f$ to the total potential (i.e. to TotalPotential::total_potential_contributions)
		 *
		 * @param[in]	total_potential_contribution	A contribution \f$\Pi_i\f$
		 */
		void
		add_total_potential_contribution(const TotalPotentialContribution<spacedim>& total_potential_contribution);

};

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_TOTALPOTENTIAL_H_ */
