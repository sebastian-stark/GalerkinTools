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

#ifndef GALERKINTOOLS_DEPENDENTFIELD_H_
#define GALERKINTOOLS_DEPENDENTFIELD_H_

#include <string>
#include <vector>

#include <galerkin_tools/config.h>
#include <galerkin_tools/independent_field.h>
#include <galerkin_tools/triangulation_system.h>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * This class allows for the definition of a term of the form
 * \f$a \dfrac{\partial^N f}{\partial x_{i_0} \hdots \partial x_{i_{N-1}}} \f$,
 * where \f$a\f$ is a constant coefficient, and \f$f\f$ is a function of the spatial
 * location \f$\boldsymbol{x}\f$. The function \f$f\f$ is derived \f$N\f$
 * times with respect to particular spatial coordinates, which is indicated by the notation
 * \f$\dfrac{\partial^N f}{\partial x_{i_0} \hdots \partial x_{i_{N-1}}}\f$. It is noted that the \f$i_k\f$
 * in the latter expression refer to fixed numbers. I.e.,
 * \f$\dfrac{\partial^N f}{\partial x_{i_0} \hdots \partial x_{i_{N-1}}}\f$ is scalar valued, with examples being
 * \f$\dfrac{\partial f}{\partial x_2}\f$ and \f$\dfrac{\partial^2 f}{\partial x_2 \partial x_1}\f$. "Zeroth"
 * derivatives are also allowed by this class (i.e., terms of the form \f$a f\f$).
 *
 * In particular, the function \f$f\f$ may either be a domain related independent field \f$u^\Omega_\epsilon\f$
 * or an interface related independent field \f$u^\Sigma_\eta\f$.
 *
 * @tparam	dim			The dimension of the object on which the function \f$f\f$ lives
 *
 * @tparam	spacedim	The spatial dimension of the problem
 */
template<unsigned int dim, unsigned int spacedim>
class DependentFieldTerm
{

public:

	/**
	 * The coefficient \f$a\f$
	 */
	const double
	coefficient;

	/**
	 * An IndependentField object, which defines together with
	 * DependentFieldTerm::component, the function \f$f\f$
	 */
	const SmartPointer<const IndependentField<dim, spacedim>>
	independent_field;

	/**
	 * A component of the (possibly vector valued) IndependentField DependentFieldTerm::independent_field, which
	 * defines together with DependentFieldTerm::independent_field the independent field \f$f\f$. If DependentFieldTerm::independent_field
	 * is scalar valued, DependentFieldTerm::component can only be zero.
	 *
	 */
	const unsigned int
	component;

	 /**
	  * A list of the spatial derivatives. The first element in the vector is \f$i_0\f$, the second
	  * element the derivative \f$i_1\f$, etc. The length of the vector is thus equal to the number of
	  * derivatives \f$N\f$; and in the special case that there are no derivatives, this vector will
	  * have the size zero.
	  */
	const std::vector<unsigned int>
	derivatives;

	/**
	 * The constructor of the class.
	 *
	 * @param[in]	coefficient			DependentFieldTerm::coefficient
	 * @param[in]	independent_field	DependentFieldTerm::independent_field
	 * @param[in]	component			DependentFieldTerm::component
	 * @param[in]	derivatives			DependentFieldTerm::derivatives
	 */
	DependentFieldTerm(	const double 							coefficient,
						const IndependentField<dim, spacedim>&	independent_field,
						const unsigned int						component = 0,
						const std::vector<unsigned int>			derivatives = std::vector<unsigned int>());

	/**
	 * @return DependentFieldTerm::derivatives[0]
	 */
	double
	first_derivative()
	const;

	/**
	 * @return	The size of DependentFieldTerm::derivatives (i.e. \f$N\f$)
	 */
	unsigned int
	n_derivatives()
	const;

	/**
	 * A comparison operator, which allows to use DependentFieldTerm in @p std::map and @p std::set.
	 * The comparison operator is based on a lexicographic ordering according to the members
	 * DependentFieldTerm::independent_field, DependentFieldTerm::component, DependentFieldTerm::derivatives.
	 * This means that the member DependentFieldTerm::coefficient does not factor into the comparison, which effectively
	 * means that two DependentFieldTerm objects are considered equal if they are associated with the same
	 * DependentFieldTerm::independent_field, DependentFieldTerm::component, DependentFieldTerm::derivatives (the same
	 * DependentFieldTerm::independent_field here really means "exactly the same object").
	 *
	 * @param[in]	dependent_field_2	The DependentFieldTerm to compare with
	 *
	 * @return 							Boolean indicating result of comparison
	 */
	bool
	operator<(const DependentFieldTerm& dependent_field_2)
	const;
};

/**
 * This class is used to define an interface related dependent field \f$e^\Sigma_\nu\f$ (\f$\nu \in N=\left\{1 \hdots N^\mathrm{e,\Sigma}\right\}\f$) according to
 * \f{equation*}
 * \begin{split}
 * e^\Sigma_\nu &=  \sum_{\eta \in H} \left[ a^\Sigma_{\nu\eta} u^\Sigma_\eta + b^\Sigma_{\nu\eta i} \dfrac{\partial u^\Sigma_\eta}{\partial x_i} \right]\\
 *              & + \sum_{\epsilon \in E} \left[ a^+_{\nu\epsilon} (u^\Omega_\epsilon)^+ + b^+_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^+ \right]\\
 *              & + \sum_{\epsilon \in E} \left[ a^-_{\nu\epsilon} (u^\Omega_\epsilon)^- + b^-_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^- \right]\\
 *              & + \sum_{\iota \in I}c^\Sigma_{\nu\iota} C_\iota + d^\Sigma_\nu.
 * \end{split}
 * \f}
 * where \f$a^\Sigma_{\nu\eta}\f$, \f$b^\Sigma_{\nu\eta i}\f$,
 * 		 \f$a^+_{\nu\epsilon}\f$, \f$b^+_{\nu\epsilon i}\f$,
 * 		 \f$a^-_{\nu\epsilon}\f$, \f$b^-_{\nu\epsilon i}\f$,
 * 		 \f$c^\Sigma_{\nu\iota}\f$, and \f$d^\Sigma_\nu\f$ are constants,
 * summation is implied over \f$i\f$,
 * and the notation \f$(\,)^-\f$ and \f$(\,)^+\f$, respectively, indicates whether a quantity
 * is evaluated on the + or - side of the interface
 * (the remaining notation is introduced in IndependentField and IndependentField<0, spacedim>).
 *
 * Although the class in principle also allows for higher spatial derivatives in the definition of
 * \f$e^\Sigma_\nu\f$, the documentation is restricted to first derivatives because most other parts of the GalerkinTools library make this restriction.
 *
 * @tparam dim		The dimension of the object on which the dependent field lives. Currently,
 * 					only @p dim = @p spacedim-1 is considered, although this class would in principle also work for
 * 					dependent fields on lower dimensional objects.
 *
 * @tparam spacedim	The spatial dimension of the problem
 */
template<unsigned int dim, unsigned int spacedim>
class DependentField
{

private:

	/**
	 * This set contains all summands making up the term
	 * \f$\sum_{\eta \in H} \left[ a^\Sigma_{\nu\eta} u^\Sigma_\eta + b^\Sigma_{\nu\eta i} \dfrac{\partial u^\Sigma_\eta}{\partial x_i} \right]\f$
	 * of the dependent field.
	 */
	std::set< DependentFieldTerm<dim, spacedim> >
	terms_interface;

	/**
	 * This set contains all summands making up the term
	 * \f$\sum_{\epsilon \in E} \left[ a^+_{\nu\epsilon} (u^\Omega_\epsilon)^+ + b^+_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^+ \right]\f$
	 * of the dependent field.
	 */
	std::set< DependentFieldTerm<dim+1, spacedim> >
	terms_neighbor_plus;

	/**
	 * This set contains all summands making up the term
	 * \f$\sum_{\epsilon \in E} \left[ a^-_{\nu\epsilon} (u^\Omega_\epsilon)^- + b^-_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^- \right]\f$
	 * of the dependent field.
	 */
	std::set< DependentFieldTerm<dim+1, spacedim> >
	terms_neighbor_minus;

	/**
	 * This set contains all summands making up the term \f$\sum_{\iota \in I}c^\Sigma_{\nu\iota} C_\iota\f$ of the dependent field.
	 */
	std::set< DependentFieldTerm<0, spacedim> >
	terms_independent_scalars;

	/**
	 * This is the constant \f$d^\Sigma_\nu\f$ of the dependent field.
	 */
	double
	constant=0.0;

public:

	/**
	 * The name of the dependent field, which is used to identify it in output, etc.
	 */
	const std::string name;

	/**
	 * Constructor for a DependentField.
	 *
	 * @param[in]	name	DependentField::name
	 */
	DependentField(	const std::string name );

	/**
	 * Add a term to DependentField::terms_interface.
	 *
	 * Although this method allows to add terms involving spatial derivatives of arbitrary order,
	 * its main purpose in the GalerkinTools library are the following two cases:
	 *
	 * (1)	Add a term \f$a^\Sigma_{\nu\eta} u^\Sigma_\eta\f$
	 * 		(@p coefficient is \f$a^\Sigma_{\nu\eta}\f$,
	 * 		 @p independent_field and @p component describe \f$u^\Sigma_\eta\f$,
	 * 		 @p derivatives is empty)
	 *
	 * (2)	Add a term  \f$b^\Sigma_{\nu \eta i} \dfrac{\partial u^\Sigma_\eta}{\partial x_i}\f$
	 * 		(@p coefficient is \f$b^\Sigma_{\nu \eta i}\f$,
	 * 		 @p independent_field and @p component describe \f$u^\Sigma_\eta\f$,
	 * 		 @p derivatives contains a single element with \f$i\f$)
	 *
	 * It is currently not allowed to add the same term (i.e. with the same
	 * @p independent_field, @p component, @p derivatives) two times.
	 *
	 * @param[in]	coefficient			DependentFieldTerm::coefficient of the term
	 *
	 * @param[in]	independent_field	DependentFieldTerm::independent_field of the term
	 *
	 * @param[in]	component			DependentFieldTerm::component of the term
	 *
	 * @param[in]	derivatives			DependentFieldTerm::derivatives of the term
	 */
	void
	add_term(	double									coefficient,
				const IndependentField<dim, spacedim>&	independent_field,
				const unsigned int						component = 0,
				const std::vector<unsigned int> 		derivatives=std::vector<unsigned int>());

	/**
	 * Add a term to DependentField::terms_neighbor_plus or DependentField::terms_neighbor_minus.
	 *
	 * Although this method allows to add terms involving spatial derivatives of arbitrary order,
	 * its main purpose in the GalerkinTools library are the following two cases:
	 *
	 * (1)	Add a term \f$a^+_{\nu\epsilon} (u^\Omega_\epsilon)^+\f$ or \f$a^-_{\nu\epsilon} (u^\Omega_\epsilon)^-\f$
	 * 		(@p coefficient is \f$a^+_{\nu\epsilon}\f$ or \f$a^-_{\nu\epsilon}\f$ depending on @p side,
	 * 		 @p independent_field and @p component describe \f$u^\Omega_\epsilon\f$,
	 * 		 @p derivatives is empty,
	 * 		 @p side specifies the side of the interface on which the independent field is evaluated)
	 *
	 * (2)	Add a term  \f$b^+_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^+\f$ or
	 * 					\f$b^-_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^-\f$
	 * 		(@p coefficient is \f$b^+_{\nu\epsilon i}\f$ or \f$b^-_{\nu\epsilon i}\f$ depending on @p side,
	 * 		 @p independent_field and @p component describe \f$u^\Omega_\epsilon\f$,
	 * 		 @p derivatives contains a single element with \f$i\f$,
	 * 		 @p side specifies the side of the interface on which the independent field is evaluated).
	 *
	 * It is currently not allowed to add the same term (i.e. with the same
	 * @p independent_field, @p component, @p derivatives, @p side) two times.
	 *
	 * @param[in]	coefficient			DependentFieldTerm::coefficient of the term
	 *
	 * @param[in]	independent_field	DependentFieldTerm::independent_field of the term
	 *
	 * @param[in]	component			DependentFieldTerm::component of the term
	 *
	 * @param[in]	derivatives			DependentFieldTerm::derivatives of the term
	 *
	 * @param[in]	side				Determines the side on which the independent field @p independent_field is evaluated
	 * 									(either ::InterfaceSide::@p minus or ::InterfaceSide::@p plus)
	 */
	void
	add_term(	double 										coefficient,
				const IndependentField<dim+1, spacedim>& 	independent_field,
				const unsigned int 							component,
				const std::vector<unsigned int> 			derivatives,
				const InterfaceSide 						side);

	/**
	 * Add a term \f$a^+_{\nu\epsilon} (u^\Omega_\epsilon)^+\f$ or \f$a^-_{\nu\epsilon} (u^\Omega_\epsilon)^-\f$
	 * to DependentField::terms_neighbor_plus or DependentField::terms_neighbor_minus, respectively
	 * (@p coefficient is \f$a^+_{\nu\epsilon}\f$ or \f$a^-_{\nu\epsilon}\f$ depending on @p side,
	 * 	@p independent_field and @p component describe \f$u^\Omega_\epsilon\f$,
	 * 	@p side specifies the side of the interface on which the independent field is evaluated).
	 *
	 * It is currently not allowed to add the same term (i.e. with the same
	 * @p independent_field, @p component, @p side) two times.
	 *
	 * @param[in]	coefficient			DependentFieldTerm::coefficient of the term (\f$a^+_{\nu\epsilon}\f$ or \f$a^-_{\nu\epsilon}\f$ depending on @p side)
	 *
	 * @param[in]	independent_field	DependentFieldTerm::independent_field of the term (describes together with @p component \f$u^\Omega_\epsilon\f$)
	 *
	 * @param[in]	component			DependentFieldTerm::component of the term (describes together with @p independent_field \f$u^\Omega_\epsilon\f$)
	 *
	 * @param[in]	side				Determines the side on which the independent field @p independent_field is evaluated
	 * 									(either ::InterfaceSide::@p minus or ::InterfaceSide::@p plus)
	 */
	void
	add_term(	double 										coefficient,
				const IndependentField<dim+1, spacedim>& 	independent_field,
				const unsigned int 							component,
				const InterfaceSide 						side);

	/**
	 * Add a term \f$b^\Sigma_{\nu \eta i} \dfrac{\partial u^\Sigma_\eta}{\partial x_i}\f$ to DependentField::terms_interface
	 * (@p coefficient is \f$b^\Sigma_{\nu \eta i}\f$,
	 * 	@p independent_field and @p component describe \f$u^\Sigma_\eta\f$,
	 *  @p derivative is \f$i\f$).
	 *
	 * It is currently not allowed to add the same term (i.e. with the same
	 * @p independent_field, @p component, @p derivative) two times.
	 *
	 * @param[in]	coefficient			DependentFieldTerm::coefficient of the term (\f$b^\Sigma_{\nu \eta i}\f$)
	 *
	 * @param[in]	independent_field	DependentFieldTerm::independent_field of the term  (describes together with @p component \f$u^\Sigma_\eta\f$)
	 *
	 * @param[in]	component			DependentFieldTerm::component of the term (describes together with @p independent_field \f$u^\Sigma_\eta\f$)
	 *
	 * @param[in]	derivative			The only element in DependentFieldTerm::derivatives of the term (\f$i\f$)
	 */
	void
	add_term(	double 									coefficient,
				const IndependentField<dim, spacedim>& 	independent_field,
				const unsigned int 						component,
				const unsigned int 						derivative);

	/**
	 * Add a term \f$b^+_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^+\f$ or
	 * 			  \f$b^-_{\nu\epsilon i} \left(\dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\right)^-\f$
	 * to DependentField::terms_neighbor_plus or DependentField::terms_neighbor_minus, respectively
	 * (@p coefficient is \f$b^+_{\nu\epsilon i}\f$ or \f$b^-_{\nu\epsilon i}\f$ depending on @p side,
	 * 	@p independent_field and @p component describe \f$u^\Omega_\epsilon\f$,
	 * 	@p derivative is \f$i\f$,
	 * 	@p side specifies the side of the interface on which the independent field is evaluated).
	 *
	 * It is currently not allowed to add the same term (i.e. with the same
	 * @p independent_field, @p component, @p derivative, @p side) two times.
	 *
	 * @param[in]	coefficient			DependentFieldTerm::coefficient of the term (\f$b^+_{\nu\epsilon i}\f$ or \f$b^-_{\nu\epsilon i}\f$ depending on @p side)
	 *
	 * @param[in]	independent_field	DependentFieldTerm::independent_field of the term (describes together with @p component \f$u^\Omega_\epsilon\f$)
	 *
	 * @param[in]	component			DependentFieldTerm::component of the term (describes together with @p independent_field \f$u^\Omega_\epsilon\f$)
	 *
	 * @param[in]	derivative			The only element in DependentFieldTerm::derivatives of the term (\f$i\f$)
	 *
	 * @param[in]	side				Determines the side on which the independent field @p independent_field is evaluated
	 * 									(either ::InterfaceSide::@p minus or ::InterfaceSide::@p plus)
	 */
	void
	add_term(	double 										coefficient,
				const IndependentField<dim+1, spacedim>& 	independent_field,
				const unsigned int 							component,
				const unsigned int 							derivative,
				const InterfaceSide 						side);

	/**
	 * Add a term \f$c^\Sigma_{\nu\iota} C_\iota\f$
	 * to DependentField::terms_independent_scalars
	 * (@p coefficient is \f$c^\Sigma_{\nu\iota}\f$,
	 *  @p independent_field is \f$C_\iota\f$).
	 *
	 * It is currently not allowed to add the same term (i.e. with the same
	 * @p independent_field) two times.
	 *
	 * @param[in]	coefficient			DependentFieldTerm::coefficient of the term (\f$c^\Sigma_{\nu\iota}\f$)
	 *
	 * @param[in]	independent_field	DependentFieldTerm::independent_field of the term (\f$C_\iota\f$)
	 */
	void
	add_term(	double 									coefficient,
				const IndependentField<0, spacedim>& 	independent_field);

	/**
	 * Set the constant DependentField::constant (which corresponds to \f$d^\Sigma_\nu\f$)
	 *
	 * It is currently not allowed to call this method again if DependentField::constant
	 * is already set to a non-zero value.
	 *
	 * @param[in]	constant			DependentField::constant (\f$d^\Sigma_\nu\f$)
	 */
	void
	add_term(	double constant);

	/**
	 * @return A const reference to DependentField::terms_interface
	 */
	const std::set< DependentFieldTerm<dim, spacedim> >&
	get_terms_interface()
	const;

	/**
	 * @param[in]	side	Determines whether DependentField::terms_neighbor_minus or
	 * 						DependentField::terms_neighbor_plus is returned
	 *
	 * @return 				A const reference to DependentField::terms_neighbor_minus or
	 * 						DependentField::terms_neighbor_plus depending on @p side
	 */
	const std::set< DependentFieldTerm<dim+1, spacedim> >&
	get_terms_neighbor(const InterfaceSide side)
	const;

	/**
	 * @return A const reference to DependentField::terms_independent_scalars
	 */
	const std::set< DependentFieldTerm<0, spacedim> >&
	get_terms_independent_scalars()
	const;

	/**
	 * @return DependentField::constant (\f$d^\Sigma_\nu\f$)
	 */
	double
	get_constant()
	const;

	/**
	 * @return	A set containing all interface related IndependentField objects currently involved in the definition
	 * 			of the DependentField
	 */
	std::set<const IndependentField<dim, spacedim>*>
	get_independent_fields_interface()
	const;

	/**
	 * @return	A set containing all domain related IndependentField objects currently involved in the definition
	 * 			of the DependentField
	 */
	std::set<const IndependentField<dim+1, spacedim>*>
	get_independent_fields_neighbors()
	const;

	/**
	 * @return	A set containing all IndependentField<0, spacedim> objects currently involved in the definition
	 * 			of the DependentField
	 */
	std::set<const IndependentField<0, spacedim>*>
	get_independent_scalars()
	const;

	/**
	 * This methods prints the current definition of the DependentField<spacedim, spacedim> to @p stdout.
	 */
	void
	print()
	const;

};

/**
 * This class is used to define a domain related dependent field \f$e^\Omega_\lambda\f$ (\f$\lambda \in \Lambda=\left\{1 \hdots N^{\mathrm{e},\Omega}\right\}\f$) according to
 * \f{equation*}
 * 		e^\Omega_\lambda = \sum_{\epsilon \in E} \left[ a^\Omega_{\lambda\epsilon} u^\Omega_\epsilon + b^\Omega_{\lambda \epsilon i} \dfrac{\partial u^\Omega_\epsilon}{\partial x_i} \right] + \sum_{\iota \in I}c^\Omega_{\lambda\iota} C_\iota  + d^\Omega_\lambda\,
 * \f}
 * where \f$a^\Omega_{\lambda\epsilon}\f$, \f$b^\Omega_{\lambda \epsilon i}\f$, \f$c^\Omega_{\lambda\iota}\f$, and \f$d^\Omega_\lambda\f$ are constants and summation is implied over \f$i\f$
 * (the remaining notation is introduced in IndependentField and IndependentField<0, spacedim>). Although the class in principle also allows for higher spatial derivatives in the definition of
 * \f$e^\Omega_\lambda\f$, the documentation is restricted to first derivatives because most other parts of the GalerkinTools library make this restriction.
 *
 * @tparam spacedim	The spatial dimension of the problem
 */
template<unsigned int spacedim>
class DependentField<spacedim, spacedim>
{

private:

	/**
	 * This set contains all summands making up the term
	 * \f$\sum_{\epsilon \in E} \left[ a^\Omega_{\lambda\epsilon} u^\Omega_\epsilon + b^\Omega_{\lambda \epsilon i} \dfrac{\partial u^\Omega_\epsilon}{\partial x_i} \right]\f$
	 * of the dependent field.
	 */
	std::set<DependentFieldTerm<spacedim, spacedim>>
	terms_domain;

	/**
	 * This set contains all summands making up the term \f$\sum_{\iota \in I}c^\Omega_{\lambda\iota} C_\iota\f$ of the dependent field.
	 */
	std::set<DependentFieldTerm<0, spacedim>>
	terms_independent_scalars;

	/**
	 * This is the constant \f$d^\Omega_\lambda\f$ of the dependent field.
	 */
	double
	constant = 0.0;

public:

	/**
	 * The name of the dependent field, which is used to identify it in output, etc.
	 */
	const std::string name;

	/**
	 * Constructor for a DependentField<spacedim, spacedim>
	 *
	 * @param[in]	name	DependentField<spacedim, spacedim>::name
	 */
	DependentField(	const std::string name );

	/**
	 * Add a term to DependentField<spacedim, spacedim>::terms_domain.
	 *
	 * Although this method allows to add terms involving spatial derivatives of arbitrary order,
	 * its main purpose in the GalerkinTools library are the following two cases:
	 *
	 * (1)	Add a term \f$a^\Omega_{\lambda\epsilon} u^\Omega_\epsilon\f$
	 * 		(@p coefficient is \f$a^\Omega_{\lambda\epsilon}\f$,
	 * 		 @p independent_field and @p component describe \f$u^\Omega_\epsilon\f$,
	 * 		 @p derivatives is empty).
	 *
	 * (2)	Add a term  \f$b^\Omega_{\lambda \epsilon i} \dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\f$
	 * 		(@p coefficient is \f$b^\Omega_{\lambda \epsilon i}\f$,
	 * 		 @p independent_field and @p component describe \f$u^\Omega_\epsilon\f$,
	 * 		 @p derivatives contains a single element with \f$i\f$).
	 *
	 * It is currently not allowed to add the same term (i.e. with the same
	 * @p independent_field, @p component, @p derivatives) two times.
	 *
	 * @param[in]	coefficient			DependentFieldTerm::coefficient of the term
	 *
	 * @param[in]	independent_field	DependentFieldTerm::independent_field of the term
	 *
	 * @param[in]	component			DependentFieldTerm::component of the term
	 *
	 * @param[in]	derivatives			DependentFieldTerm::derivatives of the term
	 */
	void
	add_term(	double 										coefficient,
				const IndependentField<spacedim, spacedim>& independent_field,
				const unsigned int 							component = 0,
				const std::vector<unsigned int> 			derivatives=std::vector<unsigned int>());

	/**
	 * Add a term \f$b^\Omega_{\lambda \epsilon i} \dfrac{\partial u^\Omega_\epsilon}{\partial x_i}\f$
	 * to DependentField<spacedim, spacedim>::terms_domain
	 * (@p coefficient is \f$b^\Omega_{\lambda \epsilon i}\f$,
	 *  @p independent_field and @p component describe \f$u^\Omega_\epsilon\f$,
	 *  @p derivative is \f$i\f$).
	 *
	 * It is currently not allowed to add the same term (i.e. with the same
	 * @p independent_field, @p component, @p derivative) two times.
	 *
	 * @param[in]	coefficient			DependentFieldTerm::coefficient of the term (\f$b^\Omega_{\lambda \epsilon i}\f$)
	 *
	 * @param[in]	independent_field	DependentFieldTerm::independent_field of the term (describes together with @p component \f$u^\Omega_\epsilon\f$)
	 *
	 * @param[in]	component			DependentFieldTerm::component of the term (describes together with @p independent_field \f$u^\Omega_\epsilon\f$)
	 *
	 * @param[in]	derivative			The only element in DependentFieldTerm::derivatives of the term (\f$i\f$)
	 */
	void
	add_term(	double 										coefficient,
				const IndependentField<spacedim, spacedim>& independent_field,
				const unsigned int 							component,
				const unsigned int 							derivative);

	/**
	 * Add a term \f$c^\Omega_{\lambda\iota} C_\iota\f$
	 * to DependentField<spacedim, spacedim>::terms_independent_scalars
	 * (@p coefficient is \f$c^\Omega_{\lambda\iota}\f$,
	 *  @p independent_field is \f$C_\iota\f$).
	 *
	 * It is currently not allowed to add the same term (i.e. with the same
	 * @p independent_field) two times.
	 *
	 * @param[in]	coefficient			DependentFieldTerm::coefficient of the term (\f$c^\Omega_{\lambda\iota}\f$)
	 *
	 * @param[in]	independent_field	DependentFieldTerm::independent_field of the term (\f$C_\iota\f$)
	 */
	void
	add_term(	double 									coefficient,
				const IndependentField<0, spacedim>&	independent_field);

	/**
	 * Set the constant DependentField<spacedim, spacedim>::constant (which corresponds to \f$d^\Omega_\lambda\f$)
	 *
	 * It is currently not allowed to call this method again if DependentField<spacedim, spacedim>::constant
	 * is already set to a non-zero value.
	 *
	 * @param[in]	constant			DependentField<spacedim, spacedim>::constant (\f$d^\Omega_\lambda\f$)
	 */
	void
	add_term(double constant);

	/**
	 * @return A const reference to DependentField<spacedim, spacedim>::terms_domain
	 */
	const std::set< DependentFieldTerm<spacedim, spacedim> >&
	get_terms_domain()
	const;

	/**
	 * @return A const reference to DependentField<spacedim, spacedim>::terms_independent_scalars
	 */
	const std::set< DependentFieldTerm<0, spacedim> >&
	get_terms_independent_scalars()
	const;

	/**
	 * @return DependentField<spacedim, spacedim>::constant (\f$d^\Omega_\lambda\f$)
	 */
	double
	get_constant()
	const;

	/**
	 * @return	A set containing all IndependentField objects currently involved in the definition
	 * 			of the DependentField<spacedim, spacedim>
	 */
	std::set<const IndependentField<spacedim, spacedim>*>
	get_independent_fields_domain()
	const;

	/**
	 * @return	A set containing all IndependentField<0, spacedim> objects currently involved in the definition
	 * 			of the DependentField<spacedim, spacedim>
	 */
	std::set<const IndependentField<0, spacedim>*>
	get_independent_scalars()
	const;

	/**
	 * This methods prints the current definition of the DependentField to @p stdout.
	 */
	void
	print()
	const;

};

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_DEPENDENTFIELD_H_ */
