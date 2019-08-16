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

#include <galerkin_tools/independent_field.h>

#include <deal.II/fe/fe_nothing.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int dim, unsigned int spacedim>
IndependentField<dim, spacedim>::IndependentField(	const string						name,
													const FiniteElement<dim, spacedim>&	fe,
													const unsigned int					n_components,
													const set<types::material_id>		non_zero_regions,
													const Function<spacedim>* const		initial_vals)
:
name(name),
fe(fe.clone()),
n_components(n_components),
non_zero_regions(non_zero_regions),
initial_vals(initial_vals)
{
	Assert(	fe.n_components() == 1,
			ExcMessage("Vector valued finite elements are not allowed with this constructor, use the constructor without the argument n_components!"));
	if(initial_vals != nullptr)
		Assert(	initial_vals->n_components == this->n_components,
				ExcMessage("The dealii::Function object for the initial values must have exactly the same number of components as the independent field!"));
}

template<unsigned int dim, unsigned int spacedim>
IndependentField<dim, spacedim>::IndependentField(	const string						name,
													const FiniteElement<dim, spacedim>&	fe,
													const set<types::material_id>		non_zero_regions,
													const Function<spacedim>* const		initial_vals)
:
name(name),
fe(fe.clone()),
n_components(fe.n_components()),
non_zero_regions(non_zero_regions),
initial_vals(initial_vals)
{
	if(initial_vals != nullptr)
		Assert( initial_vals->n_components == this->n_components,
				ExcMessage("The dealii::Function object for the initial values must have exactly the same number of components as the independent field!"));
}

template<unsigned int dim, unsigned int spacedim>
IndependentField<dim, spacedim>::~IndependentField()
{
	Assert(	n_subscriptions() == 0,
			ExcMessage("You are about to destroy an IndependentField, which is currently in use! Make sure that all IndependentField objects live at least as long as the objects using them!"));
}

template<unsigned int spacedim>
IndependentField<0, spacedim>::IndependentField(const string	name,
												const double	initial_value)
:
name(name),
initial_value(initial_value)
{
}

template<unsigned int spacedim>
IndependentField<0, spacedim>::~IndependentField()
{
	Assert(	n_subscriptions() == 0,
			ExcMessage("You are about to destroy an IndependentField, which is currently in use! Make sure that all IndependentField objects live at least as long as the objects using them!"));
}

//Instantiations
template class IndependentField<3,3>;
template class IndependentField<2,2>;
template class IndependentField<2,3>;
template class IndependentField<1,2>;
template class IndependentField<0,3>;
template class IndependentField<0,2>;

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
