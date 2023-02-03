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

#include <galerkin_tools/dependent_field.h>

#include <iostream>
#include <stdio.h>
#include <sstream>

#include <deal.II/base/exceptions.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int dim, unsigned int spacedim>
DependentFieldTerm<dim, spacedim>::DependentFieldTerm(	const double 							coefficient,
														const IndependentField<dim, spacedim>&	independent_field,
														const unsigned int						component,
														const vector<unsigned int>				derivatives)
:
coefficient(coefficient),
independent_field(&independent_field),
component(component),
derivatives(derivatives)
{
	if(dim == 0)
	{
		Assert(component == 0, ExcMessage("An independent scalar does not have components and, therefore, component must be zero!"));
		Assert(derivatives.size() == 0, ExcMessage("An independent scalar does not have derivatives and, therefore, the size of derivatives must be zero!"));
	}
	Assert(component < independent_field.n_components, ExcMessage("The independent field does not have the component required by the definition of the dependent field!"));

	if(derivatives.size() > 0)
		for(const auto& derivative : derivatives)
		{
			(void)derivative;	//silence unused parameter compiler warnings in Release mode
			Assert(derivative < spacedim, ExcMessage("The requested spatial derivative does not exist."));
		}
}

template<unsigned int dim, unsigned int spacedim>
double
DependentFieldTerm<dim, spacedim>::first_derivative()
const
{
	Assert(derivatives.size() > 0, ExcMessage("The dependent field term does not involve a derivative!"));
	return derivatives[0];
}

template<unsigned int dim, unsigned int spacedim>
unsigned int
DependentFieldTerm<dim, spacedim>::n_derivatives()
const
{
	return derivatives.size();
}

template<unsigned int dim, unsigned int spacedim>
bool
DependentFieldTerm<dim, spacedim>::operator<(const DependentFieldTerm& dependent_field_2)
const
{
    return tie(coefficient, independent_field, component, derivatives) < tie(dependent_field_2.coefficient, dependent_field_2.independent_field, dependent_field_2.component, dependent_field_2.derivatives);
}

template<unsigned int dim, unsigned int spacedim>
bool
DependentFieldTerm<dim, spacedim>::operator==(const DependentFieldTerm& dependent_field_2)
const
{
    return tie(coefficient, independent_field, component, derivatives) == tie(dependent_field_2.coefficient, dependent_field_2.independent_field, dependent_field_2.component, dependent_field_2.derivatives);
}

template<unsigned int dim, unsigned int spacedim>
bool
DependentFieldComparatorWithoutCoefficient<dim, spacedim>::operator()(	const DependentFieldTerm<dim, spacedim>& dependent_field_1,
																		const DependentFieldTerm<dim, spacedim>& dependent_field_2)
const
{
    return tie(dependent_field_1.independent_field, dependent_field_1.component, dependent_field_1.derivatives) < tie(dependent_field_2.independent_field, dependent_field_2.component, dependent_field_2.derivatives);

}

template<unsigned int dim, unsigned int spacedim>
DependentField<dim, spacedim>::DependentField(const string name )
:
name(name)
{
}

template<unsigned int dim, unsigned int spacedim>
void
DependentField<dim, spacedim>::add_term(	double 									coefficient,
											const IndependentField<dim, spacedim>& 	independent_field,
											const unsigned int 						component,
											const vector<unsigned int> 				derivatives)
{
	Assert(is_local != true, ExcMessage("No further terms can be added to a local dependent field!"));
	if(independent_field.is_local)
	{
		Assert(coefficient == 1.0, ExcMessage("The coefficient must be 1.0 for a local dependent field!"));
		Assert(derivatives.size() == 0, ExcMessage("A local dependent field may not involve derivatives of a local independent field!"));
		Assert(terms_interface.size() + terms_neighbor_minus.size() + terms_neighbor_plus.size() + terms_independent_scalars.size() == 0, ExcMessage("A local independent field must be the only contribution to a local dependent field"));
		Assert(constant == 0.0, ExcMessage("A local independent field must be the only contribution to a local dependent field"));
		is_local = true;
		is_locally_eliminated = independent_field.is_locally_eliminated;
	}

	const DependentFieldTerm<dim, spacedim> term(coefficient, independent_field, component, derivatives);
	Assert(terms_interface.find(term) == terms_interface.end(), ExcMessage("You are trying to add a term to the dependent field, which already exists. This is not allowed!"));
	terms_interface.insert(term);
}

template<unsigned int dim, unsigned int spacedim>
void
DependentField<dim, spacedim>::add_term(	double 										coefficient,
											const IndependentField<dim+1, spacedim>& 	independent_field,
											const unsigned int 							component,
											const vector<unsigned int> 					derivatives,
											const InterfaceSide 						side)
{
	Assert(is_local != true, ExcMessage("No further terms can be added to a local dependent field!"));
	Assert(independent_field.is_local != true, ExcMessage("This function cannot be used in conjunction with a local independent field!"));
	const DependentFieldTerm<dim+1, spacedim> term(coefficient, independent_field, component, derivatives);
	if(side == InterfaceSide::minus)
	{
		Assert(terms_neighbor_minus.find(term) == terms_neighbor_minus.end(), ExcMessage("You are trying to add a term to the dependent field, which already exists. This is not allowed!"));
		terms_neighbor_minus.insert(term);
	}
	else
	{
		Assert(terms_neighbor_plus.find(term) == terms_neighbor_plus.end(), ExcMessage("You are trying to add a term to the dependent field, which already exists. This is not allowed!"));
		terms_neighbor_plus.insert(term);
	}
}

template<unsigned int dim, unsigned int spacedim>
void
DependentField<dim, spacedim>::add_term(	double 										coefficient,
											const IndependentField<dim+1, spacedim>& 	independent_field,
											const unsigned int 							component,
											const InterfaceSide 						side)
{
	Assert(is_local != true, ExcMessage("No further terms can be added to a local dependent field!"));
	Assert(independent_field.is_local != true, ExcMessage("This function cannot be used in conjunction with a local independent field!"));
	add_term(coefficient, independent_field, component, vector<unsigned int>(), side);
}

template<unsigned int dim, unsigned int spacedim>
void
DependentField<dim, spacedim>::add_term(	double 									coefficient,
											const IndependentField<dim, spacedim>& 	independent_field,
											const unsigned int 						component,
											const unsigned int 						derivative)
{
	Assert(is_local != true, ExcMessage("No further terms can be added to a local dependent field!"));
	Assert(independent_field.is_local != true, ExcMessage("This function cannot be used in conjunction with a local independent field!"));
	const vector<unsigned int> derivatives = {derivative};
	add_term(coefficient, independent_field, component, derivatives);
}

template<unsigned int dim, unsigned int spacedim>
void
DependentField<dim, spacedim>::add_term(	double 										coefficient,
											const IndependentField<dim+1, spacedim>& 	independent_field,
											const unsigned int 							component,
											const unsigned int 							derivative,
											const InterfaceSide 						side)
{
	Assert(is_local != true, ExcMessage("No further terms can be added to a local dependent field!"));
	Assert(independent_field.is_local != true, ExcMessage("This function cannot be used in conjunction with a local independent field!"));
	const vector<unsigned int> derivatives = {derivative};
	add_term(coefficient, independent_field, component, derivatives, side);
}

template<unsigned int dim, unsigned int spacedim>
void
DependentField<dim, spacedim>::add_term(double 									coefficient,
										const IndependentField<0, spacedim>& 	independent_field)
{
	Assert(is_local != true, ExcMessage("No further terms can be added to a local dependent field!"));
	const DependentFieldTerm<0, spacedim> term(coefficient, independent_field);
	Assert(terms_independent_scalars.find(term) == terms_independent_scalars.end(), ExcMessage("You are trying to add a term to the dependent field, which already exists. This is not allowed!"));
	terms_independent_scalars.insert(term);
}

template<unsigned int dim, unsigned int spacedim>
void
DependentField<dim, spacedim>::add_term(double constant)
{
	Assert(is_local != true, ExcMessage("No further terms can be added to a local dependent field!"));
	Assert(this->constant == 0.0, ExcMessage("You are trying to set the constant term a second time. This is not allowed except if the constant was zero before!"));
	this->constant = constant;
}

template<unsigned int dim, unsigned int spacedim>
const set< DependentFieldTerm<dim, spacedim>, DependentFieldComparatorWithoutCoefficient<dim, spacedim> >&
DependentField<dim, spacedim>::get_terms_interface()
const
{
	return terms_interface;
}

template<unsigned int dim, unsigned int spacedim>
const set< DependentFieldTerm<dim+1, spacedim>, DependentFieldComparatorWithoutCoefficient<dim+1, spacedim> >&
DependentField<dim, spacedim>::get_terms_neighbor(const InterfaceSide side)
const
{
	if(side == InterfaceSide::minus)
		return terms_neighbor_minus;
	else
		return terms_neighbor_plus;
}

template<unsigned int dim, unsigned int spacedim>
const set< DependentFieldTerm<0, spacedim>, DependentFieldComparatorWithoutCoefficient<0, spacedim> >&
DependentField<dim, spacedim>::get_terms_independent_scalars()
const
{
	return terms_independent_scalars;
}

template<unsigned int dim, unsigned int spacedim>
double
DependentField<dim, spacedim>::get_constant()
const
{
	return constant;
}

template<unsigned int dim, unsigned int spacedim>
set<const IndependentField<dim, spacedim>*>
DependentField<dim, spacedim>::get_independent_fields_interface()
const
{
	set<const IndependentField<dim, spacedim>*> independent_fields;
	for(const auto& term : terms_interface)
		independent_fields.insert(term.independent_field);
	return independent_fields;
}

template<unsigned int dim, unsigned int spacedim>
set<const IndependentField<dim+1, spacedim>*>
DependentField<dim, spacedim>::get_independent_fields_neighbors()
const
{
	set<const IndependentField<dim+1, spacedim>*> independent_fields;
	for(const auto& term : terms_neighbor_minus)
		independent_fields.insert(term.independent_field);
	for(const auto& term : terms_neighbor_plus)
		independent_fields.insert(term.independent_field);
	return independent_fields;
}

template<unsigned int dim, unsigned int spacedim>
set<const IndependentField<0, spacedim>*>
DependentField<dim, spacedim>::get_independent_scalars()
const
{
	set<const IndependentField<0, spacedim>*> independent_scalars;
	for(const auto& term : terms_independent_scalars)
		independent_scalars.insert(term.independent_field);
	return independent_scalars;
}

template<unsigned int dim, unsigned int spacedim>
void
DependentField<dim, spacedim>::print()
const
{
	cout <<this->name << "=";
	set<string> output_terms;
	for(const auto& term : terms_interface)
	{
		stringstream output_term;
		output_term  << term.coefficient << "*" << term.independent_field->name << "_" << term.component;
		if(term.n_derivatives() > 0)
		{
			output_term << ",";
			for(const auto derivative : term.derivatives)
				output_term << derivative;
		}
		output_terms.insert(output_term.str());
	}

	for(const auto& term : terms_independent_scalars)
	{
		stringstream output_term;
		output_term << term.coefficient << "*" << term.independent_field->name;
		output_terms.insert(output_term.str());
	}

	for(const auto& term : terms_neighbor_minus)
	{
		stringstream output_term;
		output_term << term.coefficient << "*" << term.independent_field->name << "_" << term.component;
		if(term.n_derivatives() > 0)
		{
			output_term << ",";
			for(const auto derivative : term.derivatives)
				output_term << derivative;
		}
		output_term << "(-)";
		output_terms.insert(output_term.str());
	}

	for(const auto& term : terms_neighbor_plus)
	{
		stringstream output_term;
		output_term << term.coefficient << "*" << term.independent_field->name << "_" << term.component;
		if(term.n_derivatives() > 0)
		{
			output_term << ",";
			for(const auto derivative : term.derivatives)
				output_term << derivative;
		}
		output_term << "(+)";
		output_terms.insert(output_term.str());
	}

	for(const auto& output_term : output_terms)
		cout << output_term << " + ";

	stringstream output_term;
	output_term << constant;
	cout << output_term.str() << endl;
}

template<unsigned int dim, unsigned int spacedim>
bool
DependentField<dim, spacedim>::get_is_local()
const
{
	return is_local;
}

template<unsigned int dim, unsigned int spacedim>
bool
DependentField<dim, spacedim>::get_is_locally_eliminated()
const
{
	return is_locally_eliminated;
}

template<unsigned int dim, unsigned int spacedim>
bool
DependentField<dim, spacedim>::operator<(const DependentField& dependent_field_2)
const
{
    return tie(terms_interface, terms_neighbor_plus, terms_neighbor_minus, terms_independent_scalars, constant) < tie(dependent_field_2.terms_interface, dependent_field_2.terms_neighbor_plus, dependent_field_2.terms_neighbor_minus, dependent_field_2.terms_independent_scalars, dependent_field_2.constant);
}

template<unsigned int dim, unsigned int spacedim>
bool
DependentField<dim, spacedim>::operator==(const DependentField& dependent_field_2)
const
{
    return tie(terms_interface, terms_neighbor_plus, terms_neighbor_minus, terms_independent_scalars, constant) == tie(dependent_field_2.terms_interface, dependent_field_2.terms_neighbor_plus, dependent_field_2.terms_neighbor_minus, dependent_field_2.terms_independent_scalars, dependent_field_2.constant);
}

template<unsigned int spacedim>
DependentField<spacedim, spacedim>::DependentField(const string name ):
	name(name)
{
}

template<unsigned int spacedim>
void
DependentField<spacedim, spacedim>::add_term(	double 										coefficient,
												const IndependentField<spacedim, spacedim>& independent_field,
												const unsigned int 							component,
												const vector<unsigned int> 					derivatives)
{
	Assert(is_local != true, ExcMessage("No further terms can be added to a local dependent field!"));
	if(independent_field.is_local)
	{
		Assert(coefficient == 1.0, ExcMessage("The coefficient must be 1.0 for a local independent field!"));
		Assert(derivatives.size() == 0, ExcMessage("A local dependent field may not involve derivatives of a local independent field!"));
		Assert(terms_domain.size() + terms_independent_scalars.size() == 0, ExcMessage("A local independent field must be the only contribution to a local dependent field"));
		Assert(constant == 0.0, ExcMessage("A local independent field must be the only contribution to a local dependent field"));
		is_local = true;
		is_locally_eliminated = independent_field.is_locally_eliminated;
	}
	const DependentFieldTerm<spacedim, spacedim> term(coefficient, independent_field, component, derivatives);
	Assert(terms_domain.find(term) == terms_domain.end(), ExcMessage("You are trying to add a term to the dependent field, which already exists. This is not allowed!"));
	terms_domain.insert(term);
}

template<unsigned int spacedim>
void
DependentField<spacedim, spacedim>::add_term(	double 										coefficient,
												const IndependentField<spacedim, spacedim>& independent_field,
												const unsigned int 							component,
												const unsigned int 							derivative)
{
	Assert(is_local != true, ExcMessage("No further terms can be added to a local dependent field!"));
	Assert(independent_field.is_local != true, ExcMessage("This function cannot be used in conjunction with a local independent field!"));
	const vector<unsigned int> derivatives = {derivative};
	add_term(coefficient, independent_field, component, derivatives);
}

template<unsigned int spacedim>
void
DependentField<spacedim, spacedim>::add_term(	double 									coefficient,
												const IndependentField<0, spacedim>& 	independent_field)
{
	Assert(is_local != true, ExcMessage("No further terms can be added to a local dependent field!"));
	const DependentFieldTerm<0, spacedim> term(coefficient, independent_field);
	Assert(terms_independent_scalars.find(term) == terms_independent_scalars.end(), ExcMessage("You are trying to add a term to the dependent field, which already exists. This is not allowed!"));
	terms_independent_scalars.insert(term);
}

template<unsigned int spacedim>
void
DependentField<spacedim, spacedim>::add_term(double constant)
{
	Assert(is_local != true, ExcMessage("No further terms can be added to a local dependent field!"));
	Assert(this->constant == 0.0, ExcMessage("You are trying to set the constant term a second time. This is not allowed except if the constant was zero before!"));
	this->constant = constant;
}

template<unsigned int spacedim>
const set< DependentFieldTerm<spacedim, spacedim>, DependentFieldComparatorWithoutCoefficient<spacedim, spacedim> >&
DependentField<spacedim, spacedim>::get_terms_domain()
const
{
	return terms_domain;
}

template<unsigned int spacedim>
const set< DependentFieldTerm<0, spacedim>, DependentFieldComparatorWithoutCoefficient<0, spacedim> >&
DependentField<spacedim, spacedim>::get_terms_independent_scalars()
const
{
	return terms_independent_scalars;
}

template<unsigned int spacedim>
double
DependentField<spacedim, spacedim>::get_constant()
const
{
	return constant;
}

template<unsigned int spacedim>
set<const IndependentField<spacedim, spacedim>*>
DependentField<spacedim, spacedim>::get_independent_fields_domain()
const
{
	set<const IndependentField<spacedim, spacedim>*> independent_fields;
	for(const auto& term : terms_domain)
		independent_fields.insert(term.independent_field);
	return independent_fields;
}

template<unsigned int spacedim>
set<const IndependentField<0, spacedim>*>
DependentField<spacedim, spacedim>::get_independent_scalars()
const
{
	set<const IndependentField<0, spacedim>*> independent_fields;
	for(const auto& term : terms_independent_scalars)
		independent_fields.insert(term.independent_field);
	return independent_fields;
}

template<unsigned int spacedim>
void
DependentField<spacedim, spacedim>::print()
const
{
	cout <<this->name << "=";
	set<string> output_terms;
	for(const auto& term : terms_domain)
	{
		stringstream output_term;
		output_term << term.coefficient << "*" << term.independent_field->name << "_" << term.component;
		if(term.n_derivatives() > 0)
		{
			output_term << ",";
			for(const auto derivative : term.derivatives)
				output_term << derivative;
		}
		output_terms.insert(output_term.str());
	}

	for(const auto& term : terms_independent_scalars)
	{
		stringstream output_term;
		output_term << term.coefficient << "*" << term.independent_field->name;
		output_terms.insert(output_term.str());
	}

	for(const auto& output_term : output_terms)
		cout << output_term << " + ";

	stringstream output_term;
	output_term << constant;
	cout << output_term.str() << endl;
}

template<unsigned int spacedim>
bool
DependentField<spacedim, spacedim>::get_is_local()
const
{
	return is_local;
}

template<unsigned int spacedim>
bool
DependentField<spacedim, spacedim>::get_is_locally_eliminated()
const
{
	return is_locally_eliminated;
}

template<unsigned int spacedim>
bool
DependentField<spacedim, spacedim>::operator<(const DependentField<spacedim,spacedim>& dependent_field_2)
const
{
    return tie(terms_domain, terms_independent_scalars, constant) < tie(dependent_field_2.terms_domain, dependent_field_2.terms_independent_scalars, dependent_field_2.constant);
}

template<unsigned int spacedim>
bool
DependentField<spacedim, spacedim>::operator==(const DependentField<spacedim,spacedim>& dependent_field_2)
const
{
    return tie(terms_domain, terms_independent_scalars, constant) == tie(dependent_field_2.terms_domain, dependent_field_2.terms_independent_scalars, dependent_field_2.constant);
}

template class DependentField<3,3>;
template class DependentField<2,2>;
template class DependentField<2,3>;
template class DependentField<1,2>;
template class DependentFieldTerm<3,3>;
template class DependentFieldTerm<2,2>;
template class DependentFieldTerm<2,3>;
template class DependentFieldTerm<1,2>;
template class DependentFieldTerm<0,3>;
template class DependentFieldTerm<0,2>;
template struct DependentFieldComparatorWithoutCoefficient<2,2>;
template struct DependentFieldComparatorWithoutCoefficient<3,3>;
template struct DependentFieldComparatorWithoutCoefficient<1,2>;
template struct DependentFieldComparatorWithoutCoefficient<2,3>;
template struct DependentFieldComparatorWithoutCoefficient<0,2>;
template struct DependentFieldComparatorWithoutCoefficient<0,3>;

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
