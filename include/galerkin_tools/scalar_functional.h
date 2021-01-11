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

#ifndef GALERKINTOOLS_SCALARFUNCTIONAL_H_
#define GALERKINTOOLS_SCALARFUNCTIONAL_H_

#include <vector>
#include <tuple>
#include <string>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/tensor.h>

#include <galerkin_tools/config.h>
#include <galerkin_tools/dependent_field.h>
#include <galerkin_tools/dof_handler_system.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * Class defining an interface related scalar functional
 * 		\f$H^\Sigma_\tau = \int_\Sigma h^\Sigma_\tau(e^\Sigma_\nu, \boldsymbol{X}) \mathrm{d}S\f$,
 * where \f$\tau \in T=\left\{1 \hdots N^\mathrm{H,\Sigma}\right\}\f$.
 *
 * "Hidden" variables involved in the definition of the functions \f$h^\Sigma_\tau\f$ are allowed for
 * (in order to allow for incorporation of e.g. classical plasticity formulations).
 *
 * The integrand \f$h^\Sigma_\tau\f$ may, besides the current values of \f$e^\Sigma_\nu\f$, also depend
 * on an arbitrary number of sets of "reference values" of \f$e^\Sigma_\nu\f$. These reference values
 * can e.g. be the values of the \f$e^\Sigma_\nu\f$ at previous instants of time. When derivatives
 * of \f$h^\Sigma_\tau\f$ w.r.t. the dependent variables are calculated, these reference values are
 * generally regarded as fixed.
 *
 * In principle, it is often sufficient that the first and second derivatives of \f$h^\Sigma_\tau\f$ are known,
 * as the values themselves do not factor into the finite element system.  However, this depends on the exact problem.
 *
 * The ScalarFunctional class inherits from Subscriptor in order to be
 * able to check that ScalarFunctional objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * @tparam dim		The dimension of the object on which the ScalarFunctional is defined. Currently,
 * 					only @p dim = @p spacedim-1 is considered, although this class would in principle also work for
 * 					scalar functionals on lower dimensional objects.
 *
 * @tparam spacedim	The spatial dimension of the problem
 */
template<unsigned int dim, unsigned int spacedim>
class ScalarFunctional : public Subscriptor
{

public:

	/**
	 * Set of types::material_id%s determining the interfacial domain of integration
	 * (i.e. those regions of the total interface \f$\Sigma\f$, where the integrand \f$h^\Sigma_\tau\f$
	 * is non-zero).
	 */
	const std::set<types::material_id>
	domain_of_integration;

	/**
	 * A @p Function, which must be supplied by the user for determination
	 * of non-zero initial values of the hidden variables (applies only if ScalarFunctional::n_hidden > 0).
	 *
	 * The number of components of ScalarFunctional::initial_vals_hidden must be equal to ScalarFunctional::n_hidden
	 *
	 * As a minimum requirement, ScalarFunctional::initial_vals_hidden must implement the method Function::value().
	 */
	const SmartPointer<const Function<spacedim>>
	initial_vals_hidden;

	/**
	 * %Vector containing the dependent fields \f$e^\Sigma_\nu\f$ whereupon the
	 * integrand \f$h^\Sigma_\tau\f$ depends. The ordering in this vector defines the ordering
	 * in vectors and matrices related to derivatives of \f$h^\Sigma_\tau\f$ w.r.t. the
	 * dependent fields \f$e^\Sigma_\nu\f$.
	 */
	const std::vector<DependentField<dim, spacedim>>
	e_sigma;

	/**
	 * This member exists only to allow for using the same code for domain and interface related scalar functionals..
	 */
	const std::vector<DependentField<dim, spacedim>>
	e_omega= std::vector<DependentField<dim, spacedim>>();

	/**
	 * Quadrature rule when integrating over the domain determined by ScalarFunctional::domain_of_integration
	 */
	const Quadrature<dim>
	quadrature;

	/**
	 * Number of "hidden" variables associated with the integrand \f$h^\Sigma_\tau\f$ of the scalar functional
	 */
	const unsigned int
	n_hidden;

	/**
	 * A name for the scalar functional
	 */
	const std::string
	name;

	/**
	 * The number of sets of reference values of the dependent fields \f$e^\Sigma_\nu\f$ involved
	 * in the definition of the integrand \f$h^\Sigma_\tau\f$
	 */
	const unsigned int
	n_ref_sets;

	/**
	 * The constructor of the class
	 *
	 * @param[in]	e_sigma					ScalarFunctional::e_sigma
	 *
	 * @param[in] 	domain_of_integration	ScalarFunctional::domain_of_integration
	 *
	 * @param[in]	quadrature				ScalarFunctional::quadrature
	 *
	 * @param[in]	name					ScalarFunctional::name
	 *
	 * @param[in]	n_ref_sets				ScalarFunctional::n_ref_sets
	 *
	 * @param[in]	n_hidden				ScalarFunctional::n_hidden
	 *
	 * @param[in]	initial_vals_hidden		ScalarFunctional::initial_vals_hidden
	 * 										(if @p n_hidden>0 and this argument is omitted, the initial values will be set to zero)
	 */
	ScalarFunctional(	const std::vector<DependentField<dim, spacedim>> 	e_sigma,
						const std::set<types::material_id> 					domain_of_integration,
						const Quadrature<dim> 								quadrature,
						const std::string 									name,
						const unsigned int 									n_ref_sets = 0,
						const unsigned int 									n_hidden = 0,
						const Function<spacedim> *const 					initial_vals_hidden = nullptr);

	/**
	 * %Function for evaluation of integrand \f$h^\Sigma_\tau\f$ and computation of
	 * first and second derivatives w.r.t. the dependent fields \f$e^\Sigma_\nu\f$.
	 *
	 * Derived user defined scalar functionals implement this function or it's other version in case local dependent fields are involved.
	 *
	 * @param[inout] 	e_sigma					Values of \f$e^\Sigma_\nu\f$ (the ordering is defined by the ordering of
	 * 											the dependent fields in ScalarFunctional::e_sigma). Note that this vector is non-const, as in the case
	 * 											of local dependent fields updated values of these may be returned.
	 *
	 * @param[in]	 	e_sigma_ref_sets		Sets of reference values of \f$e^\Sigma_\nu\f$ (as required according to ScalarFunctional::n_ref_sets)
	 *
	 * @param[inout]	hidden_vars				Values of "hidden" variables associated with the integrand \f$h^\Sigma_\tau\f$. This vector has the size
	 * 											ScalarFunctional::n_hidden. The inputs are the previously known values of the hidden variables, which
	 * 											may be overwritten by the function (the user can decide during assembly of the finite element system
	 * 											whether the returned updated values of the hidden variables are really used to update the stored
	 * 											values of the hidden variables or whether they are discarded)
	 *
	 * @param[in]		x						Location of material point at which integrand \f$h^\Sigma_\tau\f$ and derivatives thereof are evaluated
	 *
	 * @param[in]		n						Normal vector pointing from - to + side of the interface \f$\Sigma\f$ at @p x
	 *
	 * @param[out] 		h_sigma					Current value of \f$h^\Sigma_\tau\f$.
	 *
	 * @param[out] 		h_sigma_1				First derivatives of \f$h^\Sigma_\tau\f$ w.r.t. \f$e^\Sigma_\nu\f$ (the ordering is defined by the ordering of
	 * 											the dependent fields in ScalarFunctional::e_sigma, and @p h_sigma_1 will already be initialized to the correct
	 * 											size if it is called by AssemblyHelper objects)
	 *
	 * @param[out] 		h_sigma_2				Second derivatives of \f$h^\Sigma_\tau\f$ w.r.t. \f$e^\Sigma_\nu\f$ (the ordering is defined by the ordering of
	 * 											the dependent fields in ScalarFunctional::e_sigma, and @p h_sigma_2 will already be initialized to the correct
	 * 											size if it is called by AssemblyHelper objects). If \f$h^\Sigma_\tau\f$ exists, this matrix will generally be
	 * 											symmetric. However, in principle the routines do also work for cases where \f$h^\Sigma_\tau\f$ does not exist and
	 * 											@p h_sigma_2 is not symmetric.
	 *
	 * @param[in] 		requested_quantities	A tuple indicating which quantities are actually to be computed
	 * 											(e.g. (@p true, @p false, @p true) indicates that @p h_sigma and @p h_sigma_2 are to be computed)
	 *
	 * @return									@p false if the evaluation of the integrand \f$h^\Sigma_\tau\f$ and its derivatives was successful, and @p true
	 * 											if an error prevented the proper calculation of these quantities (e.g. because a dependent field, which should
	 * 											be non-negative, was actually negative)
	 */
	virtual
	bool
	get_h_sigma(Vector<double>& 						e_sigma,
				const std::vector<Vector<double>>&		e_sigma_ref_sets,
				Vector<double>&							hidden_vars,
				const Point<spacedim>&					x,
				const Tensor<1, spacedim>&				n,
				double&									h_sigma,
				Vector<double>&							h_sigma_1,
				FullMatrix<double>& 					h_sigma_2,
				const std::tuple<bool, bool, bool> 		requested_quantities)
	const = 0;

	/**
	 * This function (@see ScalarFunctional<spacedim, spacedim>::get_h_omega) exists only to allow for using the same code for domain and interface related scalar functionals.
	 * It should never be called and therefore throws an exception.
	 */
	bool
	get_h_omega(Vector<double>& 						e_omega,
				const std::vector<Vector<double>>&		e_omega_ref_sets,
				Vector<double>&							hidden_vars,
				const Point<spacedim>&					x,
				double&									h_omega,
				Vector<double>&							h_omega_1,
				FullMatrix<double>& 					h_omega_2,
				const std::tuple<bool, bool, bool> 		requested_quantities)
	const;

	/**
	 * %Function for evaluation of the maximum permissible step length
	 *
	 * @param[in]	 	e_sigma					Values of \f$e^\Sigma_\nu\f$ (the ordering is defined by the ordering of
	 * 											the dependent fields in ScalarFunctional::e_sigma)
	 *
	 * @param[in]	 	e_sigma_ref_sets		Sets of reference values of \f$e^\Sigma_\nu\f$ (as required according to ScalarFunctional::n_ref_sets)
	 *
	 * @param[in]	 	delta_e_sigma			Values of \f$\Delta e^\Sigma_\nu\f$
	 *
	 * @param[in]		hidden_vars				Values of "hidden" variables associated with the integrand \f$h^\Sigma_\tau\f$. This vector has the size
	 * 											ScalarFunctional::n_hidden.
	 *
	 * @param[in]		x						Location of material point
	 *
	 * @param[in]		n						Normal vector pointing from - to + side of the interface \f$\Sigma\f$ at @p x
	 *
	 * @return									The maximum \f$\alpha\f$ such that \f$e^\Sigma_\nu + \alpha \Delta e^\Sigma_\nu\f$ is a permissible state;
	 * 											the standard implementation returns DBL_MAX
	 */
	virtual
	double
	get_maximum_step(	const Vector<double>& 					e_sigma,
						const std::vector<Vector<double>>&		e_sigma_ref_sets,
						const Vector<double>& 					delta_e_sigma,
						const Vector<double>&					hidden_vars,
						const Point<spacedim>&					x,
						const Tensor<1, spacedim>&				n)
	const;

	/**
	 * This function (@see ScalarFunctional<spacedim, spacedim>::get_maximum_step) exists only to allow for using the same code for domain and interface related scalar functionals.
	 * It should never be called and therefore throws an exception.
	 */
	double
	get_maximum_step(	const Vector<double>& 					e_omega,
						const std::vector<Vector<double>>&		e_omega_ref_sets,
						const Vector<double>& 					delta_e_omega,
						const Vector<double>&					hidden_vars,
						const Point<spacedim>&					x)
	const;

	/**
	 * %Function comparing the computed derivatives of the integrand \f$h^\Sigma_\tau\f$ provided by ScalarFunctional::get_h_sigma()
 	 * with corresponding numerically computed finite difference based derivatives.
 	 *
 	 * In the case of the first derivative, the numerical derivatives
 	 * are obtained based on the values for \f$h^\Sigma_\tau\f$ provided by ScalarFunctional::get_h_sigma(); and in the
 	 * case of the second derivatives, the numerical derivatives are obtained based on the values for the first derivatives of
 	 * \f$h^\Sigma_\tau\f$ provided by ScalarFunctional::get_h_sigma() ).
 	 * In both cases, a simple forward finite difference approach is used. Generally, the first numerical derivative can only
 	 * be "correct" if the computation of the value of \f$h^\Sigma_\tau\f$ is correctly implemented in the function ScalarFunctional::get_h_sigma();
 	 * and likewise the second numerical derivative can only be "correct" if the computation of the first derivative of \f$h^\Sigma_\tau\f$ is correctly
 	 * implemented in the function ScalarFunctional::get_h_sigma().
 	 *
 	 * This function is essentially meant for testing of user defined functions ScalarFunctional::get_h_sigma().
	 *
	 * @param[in]	 	e_sigma					Values of \f$e^\Sigma_\nu\f$ (the ordering is defined by the ordering of
	 * 											the dependent fields in ScalarFunctional::e_sigma)
	 *
	 * @param[in]	 	e_sigma_ref_sets		Sets of reference values of \f$e^\Sigma_\nu\f$ (as required according to ScalarFunctional::n_ref_sets)
	 *
	 * @param[inout]	hidden_vars				Values of "hidden" variables associated with the integrand \f$h^\Sigma_\tau\f$. This vector has the size
	 * 											ScalarFunctional::n_hidden.
	 *
	 * @param[in]		x						Location of material point at which integrand \f$h^\Sigma_\tau\f$ and derivatives thereof are evaluated
	 *
	 * @param[in]		n						Normal vector pointing from - to + side of the interface \f$\Sigma\f$ at @p x
	 *
	 * @param[in]		detailed_printout_file	A file to which detailed printout is written if requested
	 *
	 * @param[in]		epsilon					Step width for finite difference computation
	 */
	void
	compare_derivatives_with_numerical_derivatives(	Vector<double>& 					e_sigma,
													const std::vector<Vector<double>>& 	e_sigma_ref_sets,
													Vector<double>& 					hidden_vars,
													const Point<spacedim>& 				x,
													const Tensor<1,spacedim>& 			n,
													const std::string 					detailed_printout_file = "",
													const double 						epsilon = 1e-8)
	const;

	/**
	 * The destructor of ScalarFunctional essentially checks before destruction that the
	 * ScalarFunctional object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	virtual ~ScalarFunctional();
};

/**
 * Class defining a domain related scalar functional
 * 		\f$H^\Omega_\rho = \int_\Omega h^\Omega_\rho(e^\Omega_\lambda, \boldsymbol{X}) \mathrm{d}V\f$,
 * where \f$\rho \in P=\left\{1 \hdots N^\mathrm{H,\Omega}\right\}\f$.
 *
 * "Hidden" variables involved in the definition of the functions \f$h^\Omega_\rho\f$ are allowed for
 * (in order to allow for incorporation of e.g. classical plasticity formulations).
 *
 * The integrand \f$h^\Omega_\rho\f$ may, besides the current values of \f$e^\Omega_\lambda\f$, also depend
 * on an arbitrary number of sets of "reference values" of \f$e^\Omega_\lambda\f$. These reference values
 * can e.g. be the values of the \f$e^\Omega_\lambda\f$ at previous instants of time. When derivatives
 * of \f$h^\Omega_\rho\f$ w.r.t. the dependent variables are calculated, these reference values are
 * generally regarded as fixed.
 *
 * In principle, it is often sufficient that the first and second derivatives of \f$h^\Omega_\rho\f$ are known,
 * as the values themselves do not factor into the finite element system. However, this depends on the exact problem.
 *
 * The ScalarFunctional<spacedim, spacedim> class inherits from Subscriptor in order to be
 * able to check that ScalarFunctional<spacedim, spacedim> objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * @tparam spacedim	The spatial dimension of the problem
 */
template<unsigned int spacedim>
class ScalarFunctional<spacedim, spacedim> : public Subscriptor
{

public:

	/**
	 * Set of types::material_id%s determining the domain of integration
	 * (i.e. those regions of the domain \f$\Omega\f$, where the integrand \f$h^\Omega_\rho\f$
	 * is non-zero).
	 */
	const std::set<types::material_id>
	domain_of_integration;

	/**
	 * A Function, which must be supplied by the user for determination
	 * of non-zero initial values of the hidden variables (applies only if ScalarFunctional<spacedim, spacedim>::n_hidden > 0).
	 *
	 * The number of components of ScalarFunctional<spacedim, spacedim>::initial_vals_hidden must be equal to ScalarFunctional<spacedim, spacedim>::n_hidden
	 *
	 * As a minimum requirement, ScalarFunctional<spacedim, spacedim>::initial_vals_hidden must implement the method Function::value().
	 */
	const SmartPointer<const Function<spacedim>>
	initial_vals_hidden;

	/**
	 * %Vector containing the dependent fields \f$e^\Omega_\lambda\f$ whereupon the
	 * integrand \f$h^\Omega_\rho\f$ depends. The ordering in this vector defines the ordering
	 * in vectors and matrices related to derivatives of \f$h^\Omega_\rho\f$ w.r.t. the
	 * dependent fields \f$e^\Omega_\lambda\f$.
	 */
	const std::vector<DependentField<spacedim, spacedim>>
	e_omega;

	/**
	 * This member exists only to allow for using the same code for domain and interface related scalar functionals..
	 */
	const std::vector<DependentField<spacedim, spacedim>>
	e_sigma = std::vector<DependentField<spacedim, spacedim>>();

	/**
	 * Quadrature rule when integrating over the domain determined by ScalarFunctional<spacedim, spacedim>::domain_of_integration
	 */
	const Quadrature<spacedim>
	quadrature;

	/**
	 * Number of "hidden" variables associated with the integrand \f$h^\Omega_\rho\f$ of the scalar functional
	 */
	const unsigned int
	n_hidden;

	/**
	 * A name for the scalar functional
	 */
	const std::string
	name;

	/**
	 * The number of sets of reference values of the dependent fields \f$e^\Omega_\lambda\f$ involved
	 * in the definition of the integrand \f$h^\Omega_\rho\f$
	 */
	const unsigned int
	n_ref_sets;

	/**
	 * The constructor of the class
	 *
	 * @param[in]	e_omega					ScalarFunctional<spacedim, spacedim>::e_omega
	 *
	 * @param[in] 	domain_of_integration	ScalarFunctional<spacedim, spacedim>::domain_of_integration
	 *
	 * @param[in]	quadrature				ScalarFunctional<spacedim, spacedim>::quadrature
	 *
	 * @param[in]	name					ScalarFunctional<spacedim, spacedim>::name
	 *
	 * @param[in]	n_ref_sets				ScalarFunctional<spacedim, spacedim>::n_ref_sets
	 *
	 * @param[in]	n_hidden				ScalarFunctional<spacedim, spacedim>::n_hidden
	 *
	 * @param[in]	initial_vals_hidden		ScalarFunctional<spacedim, spacedim>::initial_vals_hidden
	 * 										(if @p n_hidden>0 and this argument is omitted, the initial values will be set to zero)
	 */
	ScalarFunctional(	const std::vector<DependentField<spacedim, spacedim>>	e_omega,
						const std::set<types::material_id>						domain_of_integration,
						const Quadrature<spacedim> 								quadrature,
						const std::string 										name,
						const unsigned int 										n_ref_sets = 0,
						const unsigned int 										n_hidden = 0,
						const Function<spacedim> *const 						initial_vals_hidden = nullptr);

	/**
	 * %Function for evaluation of integrand \f$h^\Omega_\rho\f$ and computation of
	 * first and second derivatives w.r.t. the dependent fields \f$e^\Omega_\lambda\f$.
	 *
	 * This function is pure virtual, and, therefore, it is not possible to instantiate objects of this class.
	 * Rather, derived user defined classes must implement this function.
	 *
	 * @param[in]	 	e_omega					Values of \f$e^\Omega_\lambda\f$ (the ordering is defined by the ordering of
	 * 											the dependent fields in ScalarFunctional<spacedim, spacedim>::e_omega). Note that this vector is non-const, as in the case
	 * 											of local dependent fields updated values of these may be returned.
	 *
	 * @param[in]	 	e_omega_ref_sets		Sets of reference values of \f$e^\Omega_\lambda\f$ (as required according to ScalarFunctional<spacedim, spacedim>::n_ref_sets)
	 *
	 * @param[inout]	hidden_vars				Values of "hidden" variables associated with the integrand \f$h^\Omega_\rho\f$. This vector has the size
	 * 											ScalarFunctional<spacedim, spacedim>::n_hidden. The inputs are the previously known values of the hidden variables, which
	 * 											may be overwritten by the function (the user can decide during assembly of the finite element system
	 * 											whether the returned updated values of the hidden variables are really used to update the stored
	 * 											values of the hidden variables or whether they are discarded)
	 *
	 * @param[in]		x						Location of material point at which integrand \f$h^\Omega_\rho\f$ and derivatives thereof are evaluated
	 *
	 * @param[out] 		h_omega					Current value of \f$h^\Omega_\rho\f$.
	 *
	 * @param[out] 		h_omega_1				First derivatives of \f$h^\Omega_\rho\f$ w.r.t. \f$e^\Omega_\lambda\f$ (the ordering is defined by the ordering of
	 * 											the dependent fields in ScalarFunctional<spacedim, spacedim>::e_omega, and @p h_omega_1 will already be initialized to the correct
	 * 											size if it is called by AssemblyHelper objects)
	 *
	 * @param[out] 		h_omega_2				Second derivatives of \f$h^\Omega_\rho\f$ w.r.t. \f$e^\Omega_\lambda\f$ (the ordering is defined by the ordering of
	 * 											the dependent fields in ScalarFunctional<spacedim, spacedim>::e_omega, and @p h_omega_2 will already be initialized to the correct
	 * 											size if it is called by AssemblyHelper objects). If \f$h^\Omega_\rho\f$ exists, this matrix will generally be
	 * 											symmetric. However, in principle the routines do also work for cases where \f$h^\Omega_\rho\f$ does not exist and
	 * 											@p h_omega_2 is not symmetric.
	 *
	 * @param[in] 		requested_quantities	A tuple indicating which quantities are actually to be computed
	 * 											(e.g. (@p true, @p false, @p true) indicates that @p h_omega and @p h_omega_2 are to be computed)
	 *
	 * @return									@p false if the evaluation of the integrand \f$h^\Omega_\rho\f$ and its derivatives was successful, and @p true
	 * 											if an error prevented the proper calculation of these quantities (e.g. because a dependent field, which should
	 * 											be non-negative, was actually negative)
	 */
	virtual
	bool
	get_h_omega(	Vector<double>& 					e_omega,
					const std::vector<Vector<double>>&	e_omega_ref_sets,
					Vector<double>& 					hidden_vars,
					const Point<spacedim>& 				x,
					double& 							h_omega,
					Vector<double>& 					h_omega_1,
					FullMatrix<double>& 				h_omega_2,
					const std::tuple<bool, bool, bool> 	requested_quantities
					)
	const = 0;

	/**
	 * This function (@see ScalarFunctional::get_h_sigma) exists only to allow for using the same code for domain and interface related scalar functionals.
	 * It should never be called and therefore throws an exception.
	 */
	bool
	get_h_sigma(Vector<double>& 						e_sigma,
				const std::vector<Vector<double>>&		e_sigma_ref_sets,
				Vector<double>&							hidden_vars,
				const Point<spacedim>&					x,
				const Tensor<1, spacedim>&				n,
				double&									h_sigma,
				Vector<double>&							h_sigma_1,
				FullMatrix<double>& 					h_sigma_2,
				const std::tuple<bool, bool, bool> 		requested_quantities)
	const;

	/**
	 * %Function for evaluation of the maximum permissible step length
	 *
	 * @param[in]	 	e_omega					Values of \f$e^\Omega_\lambda\f$ (the ordering is defined by the ordering of
	 * 											the dependent fields in ScalarFunctional::e_omega)
	 *
	 * @param[in]	 	e_omega_ref_sets		Sets of reference values of \f$e^\Omega_\lambda\f$ (as required according to ScalarFunctional::n_ref_sets)
	 *
	 * @param[in]	 	delta_e_omega			Values of \f$\Delta e^\Omega_\lambda\f$
	 *
	 * @param[in]		hidden_vars				Values of "hidden" variables associated with the integrand \f$h^\Omega_\rho\f$. This vector has the size
	 * 											ScalarFunctional<spacedim, spacedim>::n_hidden.
	 *
	 * @param[in]		x						Location of material point
	 *
	 * @return									The maximum \f$\alpha\f$ such that \f$e^\Omega_\lambda + \alpha \Delta e^\Omega_\lambda\f$ is a permissible state;
	 * 											the standard implementation returns DBL_MAX
	 */
	virtual
	double
	get_maximum_step(	const Vector<double>& 					e_omega,
						const std::vector<Vector<double>>&		e_omega_ref_sets,
						const Vector<double>& 					delta_e_omega,
						const Vector<double>& 					hidden_vars,
						const Point<spacedim>& 					x)
	const;

	/**
	 * This function (@see ScalarFunctional::get_maximum_step) exists only to allow for using the same code for domain and interface related scalar functionals.
	 * It should never be called and therefore throws an exception.
	 */
	double
	get_maximum_step(	const Vector<double>& 					e_sigma,
						const std::vector<Vector<double>>&		e_sigma_ref_sets,
						const Vector<double>& 					delta_e_sigma,
						const Vector<double>& 					hidden_vars,
						const Point<spacedim>& 					x,
						const Tensor<1,spacedim>&				n)
	const;

	/**
	 * This allows for manually manipulating the contribution (related to this scalar functional) of a certain domain cell to the global finite element matrix and rhs.
	 *
	 * @param[in]		domain_cell												Reference to the cell
	 *
	 * @param[inout]	K_cell													The cell matrix. This matrix is indexed by the scalar functional related indexing (indices related to cells come first followed by indices related to independent scalars).
	 *
	 * @param[inout]	f_cell													The cell rhs. This vector is indexed by the scalar functional related indexing (indices related to cells come first followed by indices related to independent scalars).
	 *
	 * @param[in]		scalar_functional_indices_to_cell_shapefuns				This relates the scalar functional related indices to the cell shape function indices. In particular, the i-th component of this vector is the cell shape function
	 * 																			corresponding to the i-th scalar functional related index.
	 *
	 * @param[in]		scalar_functional_indices_to_independent_scalar_indices	This relates the scalar functional related indices to the independent scalar indices. In particular, the i-th component of this vector is the independent scalar index
	 * 																			corresponding to the (n_dofs_cell + i)-th scalar functional related index, where n_dofs_cell is the number of cell related shape functions for this scalar functional.
	 */
	virtual
	void
	modify_K_cell_f_cell(	const DomainCellDoFIterator<spacedim>&	domain_cell,
							FullMatrix<double>&						K_cell,
							Vector<double>&							f_cell,
							const std::vector<unsigned int>&		scalar_functional_indices_to_cell_shapefuns,
							const std::vector<unsigned int>&		scalar_functional_indices_to_independent_scalar_indices)
	const;

	/**
	 * %Function comparing the computed derivatives of the integrand \f$h^\Omega_\rho\f$ provided by ScalarFunctional<spacedim, spacedim>::get_h_omega()
 	 * with corresponding numerically computed finite difference based derivatives.
 	 *
 	 * In the case of the first derivative, the numerical derivatives
 	 * are obtained based on the values for \f$h^\Omega_\rho\f$ provided by ScalarFunctional<spacedim, spacedim>::get_h_omega(); and in the
 	 * case of the second derivatives, the numerical derivatives are obtained based on the values for the first derivatives of
 	 * \f$h^\Omega_\rho\f$ provided by ScalarFunctional<spacedim, spacedim>::get_h_omega() ).
 	 * In both cases, a simple forward finite difference approach is used. Generally, the first numerical derivative can only
 	 * be "correct" if the computation of the value of \f$h^\Omega_\rho\f$ is correctly implemented in the function ScalarFunctional<spacedim, spacedim>::get_h_omega();
 	 * and likewise the second numerical derivative can only be "correct" if the computation of the first derivative of \f$h^\Omega_\rho\f$ is correctly
 	 * implemented in the function ScalarFunctional<spacedim, spacedim>::get_h_omega().
 	 *
 	 * This function is essentially meant for testing of user defined functions ScalarFunctional<spacedim, spacedim>::get_h_omega().
	 *
	 * @param[in]	 	e_omega					Values of \f$e^\Omega_\lambda\f$ (the ordering is defined by the ordering of
	 * 											the dependent fields in ScalarFunctional<spacedim, spacedim>::e_omega)
	 *
	 * @param[in]	 	e_omega_ref_sets		Sets of reference values of \f$e^\Omega_\lambda\f$ (as required according to ScalarFunctional<spacedim, spacedim>::n_ref_sets)
	 *
	 * @param[inout]	hidden_vars				Values of "hidden" variables associated with the integrand \f$h^\Omega_\rho\f$. This vector has the size
	 * 											ScalarFunctional<spacedim, spacedim>::n_hidden.
	 *
	 * @param[in]		x						Location of material point at which integrand \f$h^\Omega_\rho\f$ and derivatives thereof are evaluated
	 *
	 * @param[in]		detailed_printout_file	A file to which detailed printout is written if requested
	 *
	 * @param[in]		epsilon					Step width for finite difference computation
	 */
	void
	compare_derivatives_with_numerical_derivatives(	Vector<double>& 					e_omega,
													const std::vector<Vector<double>>& 	e_omega_ref_sets,
													Vector<double>& 					hidden_vars,
													const Point<spacedim>& 				x,
													const std::string 					detailed_printout_file = "",
													const double 						epsilon = 1e-8)
	const;

	/**
	 * The destructor of ScalarFunctional<spacedim, spacedim> essentially checks before destruction that the
	 * ScalarFunctional<spacedim, spacedim> object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	virtual ~ScalarFunctional();
};

/**
 * An interface related scalar functional with functionality to eliminate local dependent fields.
 * This scalar functional uses a Newton-Raphson algorithm in order to eliminate local dependent fields in that the gradient of the scalar functional
 * w.r.t. the local dependent fields is zero. To achieve this, the values of the local dependent fields are modified in a way that the corresponding gradients become zero.
 *
 * The scalar functional may be constructed from several other scalar functionals (as a sum of these).
 *
 *
 * @tparam dim		The dimension of the object on which the ScalarFunctional is defined.
 *
 * @tparam spacedim	The spatial dimension of the problem
 */
template<unsigned int dim, unsigned int spacedim>
class ScalarFunctionalLocalElimination : public ScalarFunctional<dim, spacedim>
{

private:

	/**
	 * The scalar functionals which form this scalar functional (those scalar functionals are simply summed up)
	 */
	const std::vector<ScalarFunctional<dim,spacedim>*>
	scalar_functionals;

	/**
	 * This maps, for each underlying scalar functional, its dependent field indices to the dependent field indices of the combined scalar functional
	 */
	std::vector<std::vector<unsigned int>>
	map_dependent_fields;

	/**
	 * This contains the dependent field indices of the nonlocal dependent field (in the indexing of the combined scalar functional)
	 */
	std::vector<unsigned int>
	indices_nonlocal_dependent_fields;

	/**
	 * This contains the dependent field indices of the local dependent field (in the indexing of the combined scalar functional)
	 */
	std::vector<unsigned int>
	indices_local_dependent_fields;

	/**
	 * Safety distance to an inadmissible state within a single Newton-Raphson iteration during determination of the values of the local dependent fields (0 < @p safety_distance < 1.0).
	 * The Newton step length will be decreased such that the "distance" between the solution and the
	 * boundary of the domain of admissibility is decreased by @p safety_distance at most during a single iteration (1.0 would correspond to no safety distance at all).
	 * This is used to avoid ill-conditioning problems resulting from a too quick approach of the
	 * boundary of the domain of admissibility.
	 *
	 * Additionally, the routines will return an error if the total increment to the local variables (between the initial state and the result of the Newton-Raphson iteration)
	 * is such that the "distance" between the local solution and the boundary of the domain of admissibility is decreased by more than safety_distance.
	 */
	double
	safety_distance = 0.9;

	/**
	 * The threshold for the residual to be used during the Newton-Raphson iteration.
	 *
	 * The Hessian of the first Newton-Raphson step is used to determine a scaling for the residual such that each element of the residual is divided by the maximum norm of the corresponding row of the Hessian.
	 * The Newton-Raphson is terminated if the 2-norm of the (scaled) residual is less then sqrt(N)*threshold_residual, with N being the number of local dependent fields.
	 */
	double
	threshold_residual = 1e-10;

	/**
	 * maximum number of Newton-Raphson iterations
	 */
	unsigned int
	max_iter = 10;

	/**
	 * maximum number of cutbacks (= step bisections) during line search
	 */
	unsigned int
	max_cutbacks = 10;

	/**
	 * whether to use the bisection line search, which tries to ensure a decreasing residual between successive iterates
	 */
	bool
	use_line_search = true;

public:


	/**
	 * The constructor of the class.
	 *
	 * @param[in]	scalar_functionals		ScalarFunctionalLocalElimination::scalar_functionals
	 *
	 * @param[in]	name					ScalarFunctional::name
	 */
	ScalarFunctionalLocalElimination(	const std::vector<ScalarFunctional<dim,spacedim>*>	scalar_functionals,
										const std::string 									name = "ScalarFunctionalLocalElimination");

	/**
	 * @see ScalarFunctional::get_h_sigma
	 *
	 * This uses all get_h_sigma of the scalar functionals in ScalarFunctionalLocalElimination::scalar_functionals to construct the final result.
	 * In this context, a Newton-Raphson iteration is used to eliminate the local dependent fields such that @p h_sigma_1 is zero for all local dependent fields.
	 * The modified values of the local dependent fields need to be returned in @p e_sigma_local (the size of this quantity is equal to the number of local dependent fields and
	 * the order is according to their occurence in @p e_sigma)
	 */
	virtual
	bool
	get_h_sigma(Vector<double>&		 				e_sigma,
				const std::vector<Vector<double>>&	e_sigma_ref_sets,
				Vector<double>& 					hidden_vars,
				const Point<spacedim>& 				x,
				const Tensor<1,spacedim>& 			n,
				double& 							h_sigma,
				Vector<double>& 					h_sigma_1,
				FullMatrix<double>& 				h_sigma_2,
				const std::tuple<bool, bool, bool>	requested_quantities)
	const;

	/**
	 * @see ScalarFunctional::get_maximum_step
	 *
	 * This function calls the maximum step functions of all scalar functionals in ScalarFunctionalLocalElimination::scalar_functionals and returns the minimal possible maximum step.
	 * In this context, all values in @p delta_e_sigma corresponding to local dependent fields should be zero when calling this function since an external modification of these values is not reasonable
	 */
	virtual
	double
	get_maximum_step(	const Vector<double>& 					e_sigma,
						const std::vector<Vector<double>>&		e_sigma_ref_sets,
						const Vector<double>& 					delta_e_sigma,
						const Vector<double>&					hidden_vars,
						const Point<spacedim>&					x,
						const Tensor<1, spacedim>&				n)
	const;

	/**
	 * Sets ScalarFunctionalLocalElimination::safety_distance
	 */
	void
	set_safety_distance(const double safety_distance);

	/**
	 * Sets ScalarFunctionalLocalElimination::threshold_residual
	 */
	void
	set_threshold_residual(const double threshold_residual);

	/**
	 * Sets ScalarFunctionalLocalElimination::max_iter
	 */
	void
	set_max_iter(const unsigned int max_iter);

	/**
	 * Sets ScalarFunctionalLocalElimination::max_cutbacks
	 */
	void
	set_max_cutbacks(const unsigned int max_cutbacks);

	/**
	 * Sets ScalarFunctionalLocalElimination::use_line_search
	 */
	void
	set_use_line_search(const bool use_line_search);
};

/**
 * An interface related scalar functional with functionality to eliminate local dependent fields.
 * This scalar functional uses a Newton-Raphson algorithm in order to eliminate local dependent fields in that the gradient of the scalar functional
 * w.r.t. the local dependent fields is zero. To achieve this, the values of the local dependent fields are modified in a way that the corresponding gradients become zero.
 *
 * The scalar functional may be constructed from several other scalar functionals (as a sum of these).
 *
 *
 * @tparam dim		The dimension of the object on which the ScalarFunctional is defined.
 *
 * @tparam spacedim	The spatial dimension of the problem
 */
template<unsigned int spacedim>
class ScalarFunctionalLocalElimination<spacedim, spacedim> : public ScalarFunctional<spacedim, spacedim>
{

private:

	/**
	 * The scalar functionals which form this scalar functional (those scalar functionals are simply summed up)
	 */
	const std::vector<ScalarFunctional<spacedim,spacedim>*>
	scalar_functionals;

	/**
	 * This maps, for each underlying scalar functional, its dependent field indices to the dependent field indices of the combined scalar functional
	 */
	std::vector<std::vector<unsigned int>>
	map_dependent_fields;

	/**
	 * This contains the dependent field indices of the nonlocal dependent field (in the indexing of the combined scalar functional)
	 */
	std::vector<unsigned int>
	indices_nonlocal_dependent_fields;

	/**
	 * This contains the dependent field indices of the local dependent field (in the indexing of the combined scalar functional)
	 */
	std::vector<unsigned int>
	indices_local_dependent_fields;

	/**
	 * Safety distance to an inadmissible state within a single Newton-Raphson iteration during determination of the values of the local dependent fields (0 < @p safety_distance < 1.0).
	 * The Newton step length will be decreased such that the "distance" between the solution and the
	 * boundary of the domain of admissibility is decreased by @p safety_distance at most during a single iteration (1.0 would correspond to no safety distance at all).
	 * This is used to avoid ill-conditioning problems resulting from a too quick approach of the
	 * boundary of the domain of admissibility.
	 *
	 * Additionally, the routines will return an error if the total increment to the local variables (between the initial state and the result of the Newton-Raphson iteration)
	 * is such that the "distance" between the local solution and the boundary of the domain of admissibility is decreased by more than safety_distance.
	 */
	double
	safety_distance = 0.9;

	/**
	 * The threshold for the residual to be used during the Newton-Raphson iteration.
	 *
	 * The Hessian of the first Newton-Raphson step is used to determine a scaling for the residual such that each element of the residual is divided by the maximum norm of the corresponding row of the Hessian.
	 * The Newton-Raphson is terminated if the 2-norm of the (scaled) residual is less then sqrt(N)*threshold_residual, with N being the number of local dependent fields.
	 */
	double
	threshold_residual = 1e-10;

	/**
	 * maximum number of Newton-Raphson iterations
	 */
	unsigned int
	max_iter = 10;

	/**
	 * maximum number of cutbacks (= step bisections) during line search
	 */
	unsigned int
	max_cutbacks = 10;

	/**
	 * whether to use the bisection line search, which tries to ensure a decreasing residual between successive iterates
	 */
	bool
	use_line_search = true;

public:


	/**
	 * The constructor of the class.
	 *
	 * @param[in]	scalar_functionals		ScalarFunctionalLocalElimination<spacedim, spacedim>::scalar_functionals
	 *
	 * @param[in]	name					ScalarFunctional<spacedim, spacedim>::name
	 */
	ScalarFunctionalLocalElimination(	const std::vector<ScalarFunctional<spacedim,spacedim>*>	scalar_functionals,
										const std::string 										name = "ScalarFunctionalLocalElimination");

	/**
	 * @see ScalarFunctional<spacedim, spacedim>::get_h_omega
	 *
	 * This uses all get_h_omega of the scalar functionals in ScalarFunctionalLocalElimination<spacedim, spacedim>::scalar_functionals to construct the final result.
	 * In this context, a Newton-Raphson iteration is used to eliminate the local dependent fields such that @p h_omega_1 is zero for all local dependent fields.
	 * The modified values of the local dependent fields need to be returned in @p e_omega_local (the size of this quantity is equal to the number of local dependent fields and
	 * the order is according to their occurence in @p e_omega)
	 *
	 */
	virtual
	bool
	get_h_omega(	Vector<double>& 					e_omega,
					const std::vector<Vector<double>>&	e_omega_ref_sets,
					Vector<double>& 					hidden_vars,
					const Point<spacedim>& 				x,
					double& 							h_omega,
					Vector<double>& 					h_omega_1,
					FullMatrix<double>& 				h_omega_2,
					const std::tuple<bool, bool, bool> 	requested_quantities)
	const;

	/**
	 * @see ScalarFunctional<spacedim, spacedim>::get_maximum_step
	 *
	 * This function calls the maximum step functions of all scalar functionals in ScalarFunctionalLocalElimination::scalar_functionals and returns the minimal possible maximum step.
	 * In this context, all values in @p delta_e_omega corresponding to local dependent fields should be zero when calling this function since an external modification of these values is not reasonable
	 */
	virtual
	double
	get_maximum_step(	const Vector<double>& 					e_omega,
						const std::vector<Vector<double>>&		e_omega_ref_sets,
						const Vector<double>& 					delta_e_omega,
						const Vector<double>& 					hidden_vars,
						const Point<spacedim>& 					x)
	const;

	/**
	 * Sets ScalarFunctionalLocalElimination<spacedim,spacedim>::safety_distance
	 */
	void
	set_safety_distance(const double safety_distance);

	/**
	 * Sets ScalarFunctionalLocalElimination<spacedim,spacedim>::threshold_residual
	 */
	void
	set_threshold_residual(const double threshold_residual);


	/**
	 * Sets ScalarFunctionalLocalElimination<spacedim,spacedim>::max_iter
	 */
	void
	set_max_iter(const unsigned int max_iter);

	/**
	 * Sets ScalarFunctionalLocalElimination<spacedim,spacedim>::max_cutbacks
	 */
	void
	set_max_cutbacks(const unsigned int max_cutbacks);

	/**
	 * Sets ScalarFunctionalLocalElimination<spacedim,spacedim>::use_line_search
	 */
	void
	set_use_line_search(const bool use_line_search);

};



GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_SCALARFUNCTIONAL_H_ */
