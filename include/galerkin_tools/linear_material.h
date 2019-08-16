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

#ifndef GALERKIN_TOOLS_LINEAR_MATERIAL_H_
#define GALERKIN_TOOLS_LINEAR_MATERIAL_H_

#include <galerkin_tools/config.h>
#include <galerkin_tools/assembly_helper.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * Class defining the scalar functional \f$H^\Omega_\rho\f$ according to
 * \f$H^\Omega_\rho=\int_\Omega \left[ \dfrac{1}{2} {\boldsymbol{e}^\Omega}^\top \boldsymbol{C} {\boldsymbol{e}^\Omega} + \boldsymbol{y}^\top \boldsymbol{e}^\Omega \right] \mathrm{d}V\f$,
 * where  \f$\boldsymbol{e}^\Omega\f$ is a vector with a selection of dependent fields \f$e^\Omega_\lambda\f$, \f$\boldsymbol{C}\f$ is a matrix, and \f$\boldsymbol{y}\f$ is a vector.
 */
template<unsigned int spacedim>
class LinearMaterialDomain : public ScalarFunctional<spacedim, spacedim>
{

private:

	/**
	 * \f$\boldsymbol{C}\f$
	 */
	const FullMatrix<double>
	C;

	/**
	 * \f$\boldsymbol{y}\f$
	 */
	const Vector<double>
	y;

	/**
	 * See ScalarFunctional<spacedim, spacedim>::get_h_omega for detailed information on this method.
	 *
	 * @param[in]	 	e_omega					\f$\boldsymbol{e}^\Omega\f$
	 *
	 * @param[in]	 	e_omega_ref_sets		Sets of reference values of \f$\boldsymbol{e}^\Omega\f$. None are required here.
	 *
	 * @param[inout]	hidden_vars				Values of "hidden" variables associated with the integrand \f$h^\Omega_\rho\f$. Here, no hidden variables are present.
	 *
	 * @param[in]		x						Location of material point
	 *
	 * @param[out] 		h_omega					\f$\dfrac{1}{2} {\boldsymbol{e}^\Omega}^\top \boldsymbol{C} {\boldsymbol{e}^\Omega} + \boldsymbol{y}^\top \boldsymbol{e}^\Omega\f$
	 *
	 * @param[out] 		h_omega_1				\f$\boldsymbol{C} {\boldsymbol{e}^\Omega} + \boldsymbol{y}^\top\f$
	 *
	 * @param[out] 		h_omega_2				\f$\boldsymbol{C}\f$
	 *
	 * @param[in] 		requested_quantities	A tuple indicating which quantities are actually to be computed
	 *
	 * @return									@p false if the evaluation of @p h_omega, @p h_omega_1, and @p h_omega_2 was successful, and @p true
	 * 											if an error prevented the proper calculation of these quantities
	 */
	bool
	get_h_omega(const Vector<double>& 				e_omega,
				const std::vector<Vector<double>>&	e_omega_ref_sets,
				Vector<double>&						hidden_vars,
				const Point<spacedim>&				x,
				double&								h_omega,
				Vector<double>&						h_omega_1,
				FullMatrix<double>&					h_omega_2,
				const std::tuple<bool, bool, bool>	requested_quantities)
	const;

public:

	/**
	 * The constructor of the class.
	 *
	 * @param[in]	e_omega					ScalarFunctional<spacedim, spacedim>::e_omega
	 *
	 * @param[in]	domain_of_integration	ScalarFunctional<spacedim, spacedim>::domain_of_integration
	 *
	 * @param[in]	quadrature				ScalarFunctional<spacedim, spacedim>::quadrature
	 *
	 * @param[in]	C						LinearMaterialDomain::C
	 *
	 * @param[in]	y						LinearMaterialDomain::y
	 *
	 * @param[in]	name					ScalarFunctional<spacedim, spacedim>::name
	 */
	LinearMaterialDomain(	const std::vector<DependentField<spacedim,spacedim>>	e_omega,
							const std::set<types::material_id>						domain_of_integration,
							const Quadrature<spacedim>								quadrature,
							const FullMatrix<double>								C,
							const Vector<double>									y,
							const std::string										name = "LinearMaterialDomain");
};

/**
 * Class defining the scalar functional \f$H^\Sigma_\tau\f$ according to
 * \f$H^\Sigma_\tau=\int_\Sigma \left[ \dfrac{1}{2} {\boldsymbol{e}^\Sigma}^\top \boldsymbol{C} {\boldsymbol{e}^\Sigma} + \boldsymbol{y}^\top \boldsymbol{e}^\Sigma \right] \mathrm{d}S\f$,
 * where  \f$\boldsymbol{e}^\Sigma\f$ is a vector with a selection of dependent fields \f$e^\Sigma_\nu\f$, \f$\boldsymbol{C}\f$ is a matrix, and \f$\boldsymbol{y}\f$ is a vector.
 */
template<unsigned int spacedim>
class LinearMaterialInterface : public GalerkinTools::ScalarFunctional<spacedim-1, spacedim>
{

private:

	/**
	 * \f$\boldsymbol{C}\f$
	 */
	const FullMatrix<double>
	C;

	/**
	 * \f$\boldsymbol{y}\f$
	 */
	const Vector<double>
	y;

	/**
	 * See ScalarFunctional::get_h_sigma for detailed information on this method.
	 *
	 * @param[in]	 	e_sigma					\f$\boldsymbol{e}^\Sigma\f$
	 *
	 * @param[in]	 	e_sigma_ref_sets		Sets of reference values of \f$\boldsymbol{e}^\Sigma\f$. None are required here.
	 *
	 * @param[inout]	hidden_vars				Values of "hidden" variables associated with the integrand \f$h^\Sigma_\tau\f$. Here, no hidden variables are present.
	 *
	 * @param[in]		x						Location of material point
	 *
	 * @param[in]		n						Normal vector pointing from - to + side of the interface \f$\Sigma\f$ at @p x. Not required for this particular scalar functional.
	 *
	 * @param[out] 		h_sigma					\f$\dfrac{1}{2} {\boldsymbol{e}^\Sigma}^\top \boldsymbol{C} {\boldsymbol{e}^\Sigma} + \boldsymbol{y}^\top \boldsymbol{e}^\Sigma\f$
	 *
	 * @param[out] 		h_sigma_1				\f$\boldsymbol{C} {\boldsymbol{e}^\Sigma} + \boldsymbol{y}^\top\f$
	 *
	 * @param[out] 		h_sigma_2				\f$\boldsymbol{C}\f$
	 *
	 * @param[in] 		requested_quantities	A tuple indicating which quantities are actually to be computed
	 *
	 * @return									@p false if the evaluation of @p h_sigma, @p h_sigma_1, and @p h_sigma_2 was successful, and @p true
	 * 											if an error prevented the proper calculation of these quantities
	 */
	bool
	get_h_sigma(const Vector<double>& 				e_sigma,
				const std::vector<Vector<double>>&	e_sigma_ref_sets,
				Vector<double>& 					hidden_vars,
				const Point<spacedim>& 				x,
				const Tensor<1,spacedim>& 			n,
				double& 							h_sigma,
				Vector<double>& 					h_sigma_1,
				FullMatrix<double>& 				h_sigma_2,
				const std::tuple<bool, bool, bool>	requested_quantities)
	const;

public:

	/**
	 * The constructor of the class.
	 *
	 * @param[in]	e_sigma					ScalarFunctional::e_sigma
	 *
	 * @param[in]	domain_of_integration	ScalarFunctional::domain_of_integration
	 *
	 * @param[in]	quadrature				ScalarFunctional::quadrature
	 *
	 * @param[in]	C						LinearMaterialInterface::C
	 *
	 * @param[in]	y						LinearMaterialInterface::y
	 *
	 * @param[in]	name					ScalarFunctional::name
	 */
	LinearMaterialInterface(	const std::vector<DependentField<spacedim-1,spacedim>>	e_sigma,
								const std::set<types::material_id>						domain_of_integration,
								const Quadrature<spacedim-1>							quadrature,
								const FullMatrix<double>								C,
								const Vector<double>									y,
								const std::string										name = "LinearMaterialInterface");
};

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKIN_TOOLS_LINEAR_MATERIAL_H_ */
