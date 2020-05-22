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

#ifndef GALERKINTOOLS_TOTALPOTENTIALCONTRIBUTION_H_
#define GALERKINTOOLS_TOTALPOTENTIALCONTRIBUTION_H_

#include <vector>
#include <tuple>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <galerkin_tools/config.h>
#include <galerkin_tools/independent_field.h>
#include <galerkin_tools/scalar_functional.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * Objects of this class are used to define contributions to the total potential.
 *
 * It is assumed that the total potential can be written as
 * \f{equation*}
 * \Pi = \Pi(H^\Omega_\rho, H^\Sigma_\tau, C_\iota) = \sum_i \Pi_i(H^\Omega_\rho, H^\Sigma_\tau, C_\iota),
 * \f}
 * where the task of this class is to define a single contribution \f$\Pi_i\f$.
 *
 * The function \f$\Pi_i\f$ may, besides the current values of \f$C_\iota\f$, also depend
 * on an arbitrary number of sets of "reference values" of \f$C_\iota\f$. These reference values
 * can e.g. be the values of the \f$C_\iota\f$ at previous instants of time. When derivatives
 * of \f$\Pi_i\f$ w.r.t. its arguments \f$H^\Omega_\rho\f$, \f$H^\Sigma_\tau\f$, \f$C_\iota\f$ are
 * calculated, these reference values are generally regarded as fixed. "Reference values" of the
 * domain and interface related scalar functionals are currently not allowed for.
 *
 * In the general case, the user must implement a class inheriting from TotalPotentialContribution, which
 * in particular overwrites the method TotalPotentialContribution::get_potential_contribution() defining the
 * the form of the function \f$\Pi_i\f$. However, for two special (but very common)
 * cases the class TotalPotentialContribution can be used straight away:
 *
 * (1) \f$\Pi_i=H^\Omega_\rho\f$ - i.e., the total potential contribution is equal to a particular domain
 *     related scalar functional.
 *
 * (2) \f$\Pi_i=H^\Sigma_\tau\f$ - i.e., the total potential contribution is equal to a particular interface
 *     related scalar functional.
 *
 * In principle, it is sufficient that the first and second derivatives of \f$\Pi_i\f$ are known,
 * as the values themselves do not factor into the finite element system. Also, the class can still be used
 * if the function \f$\Pi_i\f$ does not even exist (then, the problem is simply defined in the sense
 * of a weak form by the "first derivative" of \f$\Pi_i\f$).
 *
 * The TotalPotentialContribution class inherits from Subscriptor in order to be
 * able to check that TotalPotentialContribution objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * @tparam	spacedim	Spatial dimension
 */
template<unsigned int spacedim>
class TotalPotentialContribution : public Subscriptor
{
public:

	/**
	 * This bool is @p true if (and only if), the TotalPotentialContribution is either of the form \f$\Pi_i=H^\Omega_\rho\f$ or
	 * \f$\Pi_i=H^\Sigma_\tau\f$. If this is the case, this class can be used straight away. If this is not the
	 * case, only derived classes overwriting the function TotalPotentialContribution::get_potential_contribution()
	 * can be used.
	 */
	const bool
	is_primitive;

	/**
	 * vector with the domain related scalar functionals \f$H^\Omega_\rho\f$ entering into \f$\Pi_i\f$
	 */
	const std::vector<const ScalarFunctional<spacedim, spacedim>*>
	H_omega;

	/**
	 * vector with the interface related scalar functionals \f$H^\Sigma_\tau\f$ entering into \f$\Pi_i\f$
	 */
	const std::vector<const ScalarFunctional<spacedim-1, spacedim>*>
	H_sigma;

	/**
	 * vector with the independent scalars \f$C_\iota\f$ entering into \f$\Pi_i\f$
	 */
	const std::vector<const IndependentField<0, spacedim>*>
	C;

	/**
	 * This function defines the form of \f$\Pi_i(H^\Omega_\rho, H^\Sigma_\tau, C_\iota)\f$ as well as the derivatives w.r.t. its arguments.
	 * It must be overwritten by the user in classes inheriting from TotalPotentialContribution. The standard implementation of the function
	 * will just terminate the program because it should be called under no circumstances (for cases with TotalPotentialContribution::is_primitive == @p true
	 * there is no necessity to call the function; and for all other cases the fact that the function is called indicates that the user has forgotten to overwrite
	 * it in a derived class).
	 *
	 * @param[in]	H_omega_H_sigma_C		This vector contains the values of \f$H^\Omega_\rho\f$, \f$H^\Sigma_\tau\f$, \f$C_\iota\f$ in the ordering defined by
	 * 										TotalPotentialContribution::H_omega, TotalPotentialContribution::H_sigma, TotalPotentialContribution::C (domain
	 * 										related scalar functionals are the first elements, followed by the interface related scalar functionals, and finally the
	 * 										independent scalars).
	 *
	 * @param[in]	C_ref_sets				Sets of reference values of \f$C_\iota\f$ (as many as are provided)
	 *
	 * @param[out]	Pi						Value of \f$\Pi_i\f$. As the value of this potential does not factor into the finite element system, it can be set to zero
	 * 										if desired (for example in cases where \f$\Pi_i\f$ cannot be expressed explicitly or where it does not even exist).
	 *
	 * @param[out]	Pi_1					First derivatives of \f$\Pi_i\f$ w.r.t. \f$H^\Omega_\rho\f$, \f$H^\Sigma_\tau\f$, \f$C_\iota\f$ (in the same ordering as in
	 * 										@p H_omega_H_sigma_C)
	 *
	 * @param[out]	Pi_2					Second derivatives of \f$\Pi_i\f$ w.r.t. \f$H^\Omega_\rho\f$, \f$H^\Sigma_\tau\f$, \f$C_\iota\f$ (in the same ordering as in
	 * 										@p H_omega_H_sigma_C). If there is a potential \f$\Pi_i\f$, this matrix will generally be symmetric. However, in principle
	 * 										the routines do also work for cases where \f$\Pi_i\f$ does not exist and @p Pi_2 is not symmetric.
	 *
	 * @param[in]	requested_quantities	Tuple indicating which quantities are actually to be computed
	 * 										(e.g. (@p true, @p false, @p true) indicates that @p Pi and @p Pi_2 are to be computed)
	 *
	 * @return								@p false if the evaluation of \f$\Pi_i\f$ and its derivatives was successful, and @p true
	 * 										if an error prevented the proper calculation of these quantities
	 */
	virtual
	bool
	get_potential_contribution(	const Vector<double>&				H_omega_H_sigma_C,
								const std::vector<Vector<double>>&	C_ref_sets,
								double&								Pi,
								Vector<double>&						Pi_1,
								FullMatrix<double>&					Pi_2,
								const std::tuple<bool,bool,bool>&	requested_quantities)
	const;


	/**
	 * Constructor, which is to be used in case that the potential contribution is neither of the form
	 * \f$\Pi_i=H^\Omega_\rho\f$ nor \f$\Pi_i=H^\Sigma_\tau\f$. This will set TotalPotentialContribution::is_primitive
	 * to false, which means that you'll have to implement TotalPotentialContribution::get_potential_contribution() in
	 * a derived class.
	 *
	 * @param[in]	H_omega		TotalPotentialContribution::H_omega
	 *
	 * @param[in]	H_sigma		TotalPotentialContribution::H_sigma
	 *
	 * @param[in]	C			TotalPotentialContribution::C
	 */
	TotalPotentialContribution(	const std::vector<const ScalarFunctional<spacedim, spacedim>*>&		H_omega,
								const std::vector<const ScalarFunctional<spacedim-1, spacedim>*>&	H_sigma,
								const std::vector<const IndependentField<0, spacedim>*>&			C);

	/**
	 * Constructor for the special case \f$\Pi_i=H^\Omega_\rho\f$
	 *
	 * @param[in]	H_omega		\f$H^\Omega_\rho\f$
	 */
	TotalPotentialContribution(	const ScalarFunctional<spacedim, spacedim>& H_omega);

	/**
	 * Constructor for the special case \f$\Pi_i=H^\Sigma_\tau\f$
	 *
	 * @param[in]	H_sigma		\f$H^\Sigma_\tau\f$
	 */
	TotalPotentialContribution(	const ScalarFunctional<spacedim-1, spacedim>& H_sigma);

	/**
	 * The destructor of TotalPotentialContribution essentially checks before destruction that the
	 * TotalPotentialContribution object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	virtual
	~TotalPotentialContribution();

};

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_TOTALPOTENTIALCONTRIBUTION_H_ */
