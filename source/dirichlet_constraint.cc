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

#include <galerkin_tools/dirichlet_constraint.h>

using namespace std;

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

template<unsigned int spacedim>
DirichletConstraint<spacedim>::DirichletConstraint(	const IndependentField<spacedim, spacedim>&	independent_field,
													const unsigned int							component,
													const InterfaceSide							side,
													const set<dealii::types::material_id>		domain_of_constraint,
													const Function<spacedim>* const 			constraint_inhomogeneity,
													const IndependentField<0, spacedim>*		independent_scalar,
													const Function<spacedim>* const 			coefficient_c)
:
domain_of_constraint(domain_of_constraint),
independent_field(&independent_field),
component(component),
side(side),
constraint_inhomogeneity(constraint_inhomogeneity),
independent_scalar(independent_scalar),
coefficient_c(coefficient_c)
{
	Assert(component < independent_field.n_components, ExcMessage("You are trying to constrain a degree of freedom which does not exist!"));
	if(constraint_inhomogeneity != nullptr)
		Assert(constraint_inhomogeneity->n_components==1, ExcMessage("The Function object for a constraint must have exactly one component!"));
}

template<unsigned int spacedim>
DirichletConstraint<spacedim>::~DirichletConstraint()
{
	Assert(n_subscriptions() == 0, ExcMessage("You are about to destroy a DirichletConstraint, which is currently in use! Make sure that all DirichletConstraint objects live at least as long as the objects using them!"));
}

template class DirichletConstraint<2>;
template class DirichletConstraint<3>;

GALERKIN_TOOLS_NAMESPACE_CLOSE
DEAL_II_NAMESPACE_CLOSE
