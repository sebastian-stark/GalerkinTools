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

#ifndef GALERKINTOOLS_TRIANGULATIONSYSTEM_H_
#define GALERKINTOOLS_TRIANGULATIONSYSTEM_H_

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/distributed/tria.h>

#include <galerkin_tools/config.h>

#include <boost/signals2.hpp>

DEAL_II_NAMESPACE_OPEN
GALERKIN_TOOLS_NAMESPACE_OPEN

/**
 * This enum describes the several refinement situations that are possible for a cell on an interface.
 * The following possibilities exist:
 *
 * (1)	The face of the domain cell on the minus side of the interface cell is more refined than
 * 		the face of the domain cell on the plus side of the interface. In this case, the interface
 * 		cell coincides with the domain cell on the minus side (just because the interface mesh is
 * 		assumed to follow the mesh of the more refined domain cells).
 *
 * (2)	The opposite case of (1). In this case, the interface cell coincides with the face of the
 * 		domain cell on the plus side of the interface cell
 *
 * (3)	The face of the domain cell on the minus side of the interface cell is equally refined as
 * 		the face of the domain cell on the plus side of the interface. In this case, the interface
 * 		cell coincides with both, the face of the domain cell on the minus side as well as the
 * 		face of the domain cell on the plus side.
 *
 * (4)	The interface is actually a boundary. I.e., there is only a domain cell on the minus side.
 * 		In this case, the face of the domain cell underlying the boundary is equally refined as the
 * 		corresponding interface (or, rather, boundary) cell.
 */
enum InterfaceRefinementCase
{
	minus_is_finer,
	plus_is_finer,
	equally_fine,
	at_boundary
};

/**
 * This enum describes the two sides of an interface.
 */
enum InterfaceSide
{
	minus,
	plus
};

/**
 * Objects of this class contain information about an interface cell and the neighboring domain cells.
 *
 * The InterfaceCellDomainCells class inherits from Subscriptor in order to be
 * able to check that InterfaceCellDomainCells objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * @tparam	spacedim	The spatial dimension into which the interface is embedded.
 */
template<unsigned int spacedim>
class InterfaceCellDomainCells : public Subscriptor
{

public:

	/**
	 * A convenience typedef for a deal.II TriaIterator referring to an interface cell
	 */
	typedef TriaIterator<CellAccessor<spacedim-1, spacedim>>
	InterfaceCell;

	/**
	 * A convenience typedef for a deal.II TriaIterator referring to a domain cell
	 */
	typedef TriaIterator<CellAccessor<spacedim, spacedim>>
	DomainCell;

	/**
	 * A deal.II TriaIterator to the interface cell
	 */
	const InterfaceCell interface_cell;

	/**
	 * A deal.II TriaIterator to the domain cell adjacent to the interface cell
	 * on the minus side.
	 *
	 * InterfaceCellDomainCells::domain_cell_minus is refined to the same degree as
	 * InterfaceCellDomainCells::interface_cell if
	 * InterfaceCellDomainCells::refinement_case == @p minus_is_finer or if
	 * InterfaceCellDomainCells::refinement_case == @p equally_fine or if
	 * InterfaceCellDomainCells::refinement_case == @p at_boundary.
	 * Otherwise, it is once less refined than InterfaceCellDomainCells::interface_cell.
	 */
	const DomainCell domain_cell_minus;

	/**
	 * The face number of the face of the domain cell on the minus side
	 * lying on the interface
	 */
	const unsigned int face_minus;

	/**
	 * A deal.II TriaIterator to the domain cell adjacent to the interface cell
	 * on the plus side.
	 *
	 * InterfaceCellDomainCells::domain_cell_plus is refined to the same degree as
	 * InterfaceCellDomainCells::interface_cell if
	 * InterfaceCellDomainCells::refinement_case == @p plus_is_finer or if
	 * InterfaceCellDomainCells::refinement_case == @p equally_fine;
	 * and it is once less refined than InterfaceCellDomainCells::interface_cell if
	 * InterfaceCellDomainCells::refinement_case == @p minus_is_finer.
	 * Otherwise, the interface is actually a boundary, in which case
	 * InterfaceCellDomainCells::domain_cell_plus is set to the same cell as
	 * InterfaceCellDomainCells::domain_cell_minus.
	 */
	const DomainCell domain_cell_plus;

	/**
	 * The face number of the face of the domain cell on the plus side
	 * lying on the interface
	 */
	const unsigned int face_plus;

	/**
	 * The refinement situation at the interface cell. See ::InterfaceRefinementCase
	 * for further documentation.
	 */
	const InterfaceRefinementCase refinement_case;

	/**
	 * This member applies only to the case, where the domain mesh on one side of the
	 * interface cell is refined differently than the domain mesh on the other side (i.e., where
	 * InterfaceCellDomainCells::refinement_case == @p minus_is_finer or
	 * InterfaceCellDomainCells::refinement_case == @p plus_is_finer).
	 * It represents the subface of the domain cell face on the coarser side of the interface
	 * coinciding with the interface cell (and, by assumption on the properties of the meshes,
	 * with the domain cell face on the finer side of the interface as well).
	 *
	 * See also the documentation of the CellAccessor<dim, spacedim>::neighbor_of_coarser_neighbor()
	 * method.
	 */
	const unsigned int subface;

private:

	/**
	 * This method is used internally to convert the inputs of the public constructor into the
	 * information required to set up the const members of this class with the private constructor.
	 *
	 *  @param[in]	domain_cell		An iterator to the domain cell underlying the interface
	 * 								cell @p interface_cell
	 *
	 * @param[in]	face			The face number of the domain cell @p domain_cell underlying
	 * 								the interface @p interface_cell
	 *
	 * @param[in]	interface_cell	An iterator to the interface cell, which is defined by the pair
	 * 								@p domain_cell, @p face. The constructor will assert on whether the domain cell
	 * 								face defined by the pair @p domain_cell, @p face coincides with @p interface_cell
	 * 								by checking coincidence of the face / cell centers.
	 *
	 * @param[in]	interface_side	The side of interface on which @p domain_cell lies
	 *
	 * @return						A tuple containing the members InterfaceCellDomainCells::interface_cell,
	 * 								InterfaceCellDomainCells::domain_cell_minus, InterfaceCellDomainCells::face_minus,
	 * 								InterfaceCellDomainCells::domain_cell_plus, InterfaceCellDomainCells::face_plus,
	 * 								InterfaceCellDomainCells::refinement_case, InterfaceCellDomainCells::subface
	 */
	static
	std::tuple<const InterfaceCell, const DomainCell, const unsigned int, const DomainCell, const unsigned int, const InterfaceRefinementCase, const unsigned int>
	convert_constructor_inputs(	const DomainCell&		domain_cell,
								const unsigned int		face,
								const InterfaceCell&	interface_cell,
								const InterfaceSide		interface_side);

	/**
	 * This private constructor is called by the public constructor in order to set up the member variables.
	 * The conversion between the input parameters of the public constructor and this private constructor
	 * is done by the static member function InterfaceCellDomainCells::convert_constructor_inputs. Although
	 * this indirect approach is a bit awkward, it has the advantage of maintaining constness of the
	 * members of this class.
	 *
	 * @param[in]	input	A tuple containing the members InterfaceCellDomainCells::interface_cell,
	 * 						InterfaceCellDomainCells::domain_cell_minus, InterfaceCellDomainCells::face_minus,
	 * 						InterfaceCellDomainCells::domain_cell_plus, InterfaceCellDomainCells::face_plus,
	 * 						InterfaceCellDomainCells::refinement_case, InterfaceCellDomainCells::subface
	 */
	InterfaceCellDomainCells(const std::tuple<const InterfaceCell, const DomainCell, const unsigned int, const DomainCell, const unsigned int, const InterfaceRefinementCase, const unsigned int> input);

public:

	/**
	 * The public constructor of this class.
	 *
	 * @param[in]	domain_cell		An iterator to the domain cell underlying the interface
	 * 								cell @p interface_cell
	 *
	 * @param[in]	face			The face number of the domain cell @p domain_cell underlying
	 * 								the interface @p interface_cell
	 *
	 * @param[in]	interface_cell	An iterator to the interface cell, which is defined by the pair
	 * 								@p domain_cell, @p face. The constructor will assert on whether the domain cell
	 * 								face defined by the pair @p domain_cell, @p face coincides with @p interface_cell
	 * 								by checking coincidence of the face / cell centers.
	 *
	 * @param[in]	interface_side	The side of the interface cell on which @p domain_cell lies
	 */
	InterfaceCellDomainCells(	const DomainCell&		domain_cell,
								const unsigned int		face,
								const InterfaceCell&	interface_cell,
								const InterfaceSide		interface_side);

	/**
	 * The destructor of InterfaceCellDomainCells essentially checks before destruction that the
	 * InterfaceCellDomainCells object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	~InterfaceCellDomainCells();

	/**
	 * @return	This method returns a tuple with the types::material_id%s associated with
	 * 			InterfaceCellDomainCells::interface_cell, InterfaceCellDomainCells::domain_cell_minus,
	 * 			InterfaceCellDomainCells::domain_cell_plus (in this order).
	 *
	 * 			If the interface cell is actually at the boundary, the third member of the tuple will
	 * 			be returned as numbers::invalid_material_id
	 */
	std::tuple<const types::material_id, const types::material_id, const types::material_id>
	get_material_ids()
	const;

};

/**
 * Objects of this class contain all information about the domain and interface mesh and how both of them
 * are related to each other (i.e., which domain cell faces underly the interface mesh). They are also used
 * for the actual definition of the interface mesh based on the domain mesh. Moreover, objects of this class
 * are used to ensure that the interface mesh is refined in a consistent way with the domain mesh.
 *
 * The following general notes apply to domain and interface meshes:
 *
 * (1)	The domain mesh may be split into several portions, which are identified by types::material_id%s.
 * 		In the same way, the interface mesh may be split into several portions. The domain / interface
 * 		portions must generally be defined by the coarse mesh cells. I.e., child cells must always have the same
 * 		types::material_id as the corresponding coarse cell.
 *
 * (2)	Neighboring domain cells with two different types::material_id%s will share faces / nodes where
 * 		they meet if the domain triangulation has been set up accordingly.
 *
 * (3)	Neighboring interface cells with two different types::material_id%s will not share edges/nodes where
 * 		they meet. If degrees of freedom are to be shared for these situations, this can currently only be achieved by
 * 		appropriate constraints.
 *
 * (4)	Boundaries are generally considered as special cases of interfaces. The distinguishing feature of boundaries
 * 		is just that no mesh exists on one side of the interface (and, therefore, no domain related independent fields
 * 		live on this side).
 *
 * (5)	An interface mesh must be set up wherever interface related independent fields live and/or
 * 		Dirichlet type constraints on domain related fields are applied on interfaces using the built in functionality of
 * 		AssemblyHelper.
 *
 * (6)	The definition of the interface mesh is done via the domain mesh. In particular, each interface cell
 * 		is defined by a face of a domain cell. The domain cell used in this process defines the minus side of
 * 		the interface cell. Generally, only domain cells of the coarse mesh are allowed for the definition of
 * 		the interfaces. This implies that the coarse domain mesh must be set up before the interface mesh can
 * 		be defined.
 *
 * (7)	Each interface cell is always refined in the same way as the face of the more refined
 * 		domain cell	adjacent to the interface. Also, it is assumed that refinement of the interface mesh
 * 		always happens indirectly through refinement of the domain mesh (by the requirement set out
 * 		before). I.e., it is not possible to mark interface cells for refinement. Rather, the adjacent domain
 * 		cells must be marked for refinement if the interface mesh is to be refined.
 *
 * (8)	At present only isotropic refinement is allowed for. This choice has been made because many important
 * 		functionalities of the deal.II library are based on the very same assumption currently.
  *
 * (9)	The types::manifold_id%s of interface cells will generally be taken from the types::manifold_id%s
 * 		of the coarse mesh domain cell faces underlying the interface cells (on the minus side of the interface).
 *
 * (10)	In order to avoid undefined behavior when the TriangulationSystem is used with an AssemblyHelper object,
 * 		the definition of types::material_id%s, types::manifold_id%s, etc. should always be done
 * 		on the coarse domain mesh before the TriangulationSystem is constructed. Also it should be done before any
 * 		kind of refinement in order to ensure that ids are properly inherited to children (and, therefore, the
 * 		ids of the refined mesh stay consistent with the coarse mesh).
 *
 * (11)	Due to restrictions in the deal.II library, it is presently not allowed for that more than two interface cells
 * 		share a common edge (and the corresponding nodes). Currently, these situations can only be
 * 		treated by appropriately "chopping" the interface mesh into different portions (i.e., portions with different
 * 		types::material_id%s). As described above, this approach will create duplicate edges and nodes where
 * 		interface cells with different types::material_id%s meet, meaning that it may be necessary to
 * 		constrain dof%s together on these duplicate edges.
 *
 * The TriangulationSystem class inherits from Subscriptor in order to be
 * able to check that TriangulationSystem objects are only destroyed when they are
 * not needed anymore by other objects.
 *
 * @tparam	spacedim	The spatial dimension of the problem.
 */
template<unsigned int spacedim>
class TriangulationSystem : public Subscriptor
{

public:

	/**
	 * A convenience typedef for a deal.II TriaIterator referring to an interface cell
	 */
	typedef TriaIterator<CellAccessor<spacedim-1, spacedim>>
	InterfaceCell;

	/**
	 * A convenience typedef for a deal.II TriaIterator referring to a domain cell
	 */
	typedef TriaIterator<CellAccessor<spacedim, spacedim>>
	DomainCell;

	/**
	 * A convenience typedef for a deal.II TriaActiveIterator referring to an interface cell
	 */
	typedef TriaActiveIterator<CellAccessor<spacedim-1, spacedim>>
	ActiveInterfaceCell;

	/**
	 * A convenience typedef for a deal.II TriaActiveIterator referring to a domain cell
	 */
	typedef TriaActiveIterator<CellAccessor<spacedim, spacedim>>
	ActiveDomainCell;

private:

	/**
	 * %Mapping between pairs of (coarse domain cell, coarse domain cell face) defining interface cells and
	 * corresponding types::material_id%s. This data structures determines which interface cells
	 * are associated with the same interface portion; and it is used only before
	 * TriangulationSystem::close() is called.
	 */
	std::map< std::pair<const DomainCell, const unsigned int>, types::material_id >
	coarse_domain_faces_material_ids;

	/**
	 * This stores the manifolds used for the interfae triangulation
	 */
	std::map<dealii::types::manifold_id, const dealii::Manifold<spacedim-1, spacedim>*>
	interface_manifolds;

	/**
	 * This vector contains the interface cells and corresponding domain cells on the coarse mesh level.
	 * A peculiarity of the InterfaceCellDomainCells contained in the vector is that the domain cell faces
	 * on both sides of an interface cell are always refined to the same level. I.e., the InterfaceCellDomainCells
	 * are always related with either InterfaceCellDomainCells::refinement_case == @p equally_fine
	 * or InterfaceCellDomainCells::refinement_case == @p at_boundary.
	 */
	std::vector<InterfaceCellDomainCells<spacedim>>
	coarse_interface_cell_domain_cells;

	/**
	 * This vector contains the interface cells and corresponding domain cells of the active mesh.
	 */
	std::vector<InterfaceCellDomainCells<spacedim>>
	active_interface_cell_domain_cells;

	/**
	 * A list of connections with which this object connects to the triangulations to get information about when the triangulations change.
	 */
	std::vector<boost::signals2::connection>
	tria_listeners;

	/**
	 * This method is used internally to generate / update TriangulationSystem::active_interface_cell_domain_cells after meshing /
	 * mesh refinement.
	 *
	 * @param[in]	no_assert	If @p no_assert is @p true, no checking is performed whether the interface and domain mesh
	 * 							refinements are consistent. This functionality is necessary when generating
	 * 							@p TriangulationSystem::active_interface_cell_domain_cells after initial meshing because
	 * 							the domain mesh supplied to TriangulationSystem may already be refined. In the
	 * 							latter case, the corresponding interface mesh is refined in several refinement steps
	 * 							until it is consistent with the domain mesh. However, intermediate meshes in this process
	 * 							may not be consistent with the domain mesh and, thus, it must not checked for consistency.
	 */
	void
	generate_active_interface_cells_domain_cells(const bool no_assert = false);

	/**
	 * This internal auxiliary method is called by TriangulationSystem::generate_active_interface_cells_domain_cells() and
	 * does the actual work of generating / updating TriangulationSystem::active_interface_cell_domain_cells.
	 * For each pair of (@p domain_cell, @p face) defining an interface cell on the coarse level, the method
	 * does a recursion through the refinement levels and adds elements to
	 * TriangulationSystem::active_interface_cell_domain_cells only at the deepest (i.e. active) level.
	 *
	 * @param[in]	domain_cell		Domain cell underlying the interface cell @p interface_cell
	 *
	 * @param[in]	face			Face of @p domain_cell underlying interface cell @p interface_cell
	 *
	 * @param[in]	interface_cell	The interface cell corresponding to the pair (@p domain_cell, @p face)
	 *
	 * @param[in]	no_assert		see the documentation of the same parameter in
	 * 								TriangulationSystem::generate_active_interface_cells_domain_cells()
	 */
	void
	generate_active_interface_cells_domain_cells_recursion(	const DomainCell& 		domain_cell,
															const unsigned int& 	face,
															const InterfaceCell& 	interface_cell,
															const bool 				no_assert);

	/**
	 * This internal auxiliary method is used for checking whether all children of a certain coarse domain cell have the same
	 * types::material_id as the coarse mother cell. This method is used to assert on the case that the user
	 * supplies a refined domain mesh to TriangulationSystem, where the boundaries between different domain portions
	 * are not aligned with coarse mesh faces.
	 *
	 * The method works its way through the different refinement levels by recursion.
	 *
	 * @param[in]	domain_cell		The mother cell to be checked
	 * @return						If @p true is returned, all children of the coarse domain cell have the same
	 * 								types::material_id as the coarse mother cell; if @p false is returned,
	 * 								this is not the case.
	 */
	bool
	check_material_ids_recursion(const DomainCell& domain_cell)
	const;

	/**
	 * This function takes care that the interface cells are appropriately marked for refinement/coarsening if the domain mesh
	 * is going to be refined. The underlying mechanism is that the method is attached to the domain triangulation during
	 * construction of the TriangulationSystem; and it will always be called immediately before the domain mesh is
	 * actually refined (i.e., after all refinement flags of the domain mesh have been checked for consistency
	 * and mesh smoothing has been performed).
	 *
	 * @todo This function will also have to be involved in the transfer of hidden variables between meshes.
	 * 		 However, this is not implemented yet, meaning that hidden variables can currently not be combined with
	 * 		 mesh refinement during the computation.
	 */
	virtual
	void
	pre_refinement_domain();

	/**
	 * This internal function takes care that the interface cells are refined/coarsened after any refinement of the domain cells.
	 * The underlying mechanism is that the method is attached to the domain triangulation during
	 * construction of the TriangulationSystem; and it will always be called immediately after the domain mesh has been
	 * refined.
	 *
	 * @todo  This function will also have to be involved in the transfer of hidden variables between meshes.
	 * 		  However, this is not implemented yet, meaning that hidden variables can currently not be combined with
	 * 		  mesh refinement during the computation.
	 */
	virtual
	void
	post_refinement_domain();

	/**
	 * This function disconnects all listeners in TriangulationSystem::tria_listeners and reconnects
	 * them. This provides a mechanism to make sure that the functions related to TriangulationSystem::tria_listeners
	 * are called as the very last step. This function is for internal use by AssemblyHelper.
	 */
	void
	push_tria_listeners_to_end();

	/**
	 * If this is @p true, it is tried to adjust vertices of interface and domain triangulation such that they align properly
	 * (misalignment ma happen if inapproriate/inconsistent manifolds are used for domain and interface triangulation)
	 */
	const bool
	fix_vertex_positions = false;

	/**
	 * Make AssemblyHelper a friend to simplify some procedures
	 */
	template <unsigned int> friend class AssemblyHelper;

protected:

	/**
	 * This is the triangulation of the domain. For further information of the underlying assumptions,
	 * see the general documentation of this class.
	 */
	const SmartPointer<Triangulation<spacedim, spacedim>>
	tria_domain;

	/**
	 * This is the triangulation of the interfaces. For further information of the underlying assumptions,
	 * see the general documentation of this class.
	 */
	std::unique_ptr<Triangulation<spacedim-1, spacedim>>
	tria_interface;

	/**
	 * This boolean indicates whether the definition of the interface mesh has been finished already.
	 * As soon as all interface cells are defined, the method TriangulationSystem::close() must be
	 * called by the user, which will set TriangulationSystem::closed to @p true. Afterwards, no further
	 * interface cells can be added. Vice versa, TriangulationSystem objects are not associated with
	 * an actual triangulation for the interface before TriangulationSystem::close() is called,
	 * and are, therefore, of little use.
	 */
	bool
	closed = false;

	/**
	 * %Function generating the interface triangulation from the domain triangulation and setting up
	 * TriangulationSystem::coarse_interface_cell_domain_cells and TriangulationSystem::active_interface_cell_domain_cells
	 */
	void
	generate_tria_interface_from_tria_domain();


public:

	/**
	 * All functions attached to this signal will be called before the TriangulationSystem
	 * is going to be refined (but after all domain and interface cells are marked
	 * for refinement and coarsening).
	 */
	mutable boost::signals2::signal<void ()>
	pre_refinement;

	/**
	 * All functions attached to this signal will be called after the TriangulationSystem
	 * has been refined.
	 */
	mutable boost::signals2::signal<void ()>
	post_refinement;

	/**
	 * Constructor of the TriangulationSystem. Note that the constructor does not copy the @p tria_domain
	 * object supplied, but rather stores a pointer to it. So do not change it after construction of the
	 * TriangulationSystem unless you're knowing exactly what you are doing (the only change to the
	 * domain triangulation, which is admissible, is mesh refinement). No checking on external changes
	 * of the domain triangulation is currently performed and, therefore, doing so may lead to errors
	 * which are hard to debug.
	 *
	 * @param[in]	tria_domain					TriangulationSystem::tria_domain
	 *
	 * @param[in]	fix_vertex_positions		TriangulationSystem::fix_vertex_positions
	 */
	TriangulationSystem(Triangulation<spacedim, spacedim>& 	tria_domain,
						const bool 							fix_vertex_positions = false);

	/**
	 * The destructor of TriangulationSystem essentially checks before destruction that the
	 * TriangulationSystem object is not used by other objects. If this is the case, the program
	 * will be aborted.
	 */
	virtual
	~TriangulationSystem();

	/**
	 * @return	This returns a const reference to the domain triangulation.
	 */
	virtual
	const Triangulation<spacedim, spacedim>&
	get_triangulation_domain()
	const;

	/**
	 * @return	This returns a const reference to the interface triangulation.
	 */
	virtual
	const Triangulation<spacedim-1, spacedim>&
	get_triangulation_interface()
	const;

	/**
	 * @return	This returns a reference to the domain triangulation.
	 */
	virtual
	Triangulation<spacedim, spacedim>&
	get_triangulation_domain();

	/**
	 * @return	This returns a reference to the interface triangulation.
	 */
	virtual
	Triangulation<spacedim-1, spacedim>&
	get_triangulation_interface();

	/**
	 * %Function closing the interface definition. After this function has been called, no modifications to
	 * the interface triangulation are possible. This function does the real work of creating the interface
	 * triangulation and associating it with the underlying domain triangulation.
	 */
	virtual
	void
	close();

	/**
	 * This method adds an interface cell to the interface mesh, with the interface cell
	 * being defined through the pair (@p domain_cell, @p face). The @p domain_cell is assumed to be on the
	 * minus side of the interface.
	 *
	 * Note that this function can only be called before TriangulationSystem::close() has been called.
	 *
	 * Note also that the method overwrites previously added cells (i.e., if the same combination of
	 * cell and face has been added previously, the corresponding @p material_id will be overwritten).
	 *
	 * @param[in]	domain_cell	The domain cell underlying the interface cell to be defined
	 *
	 * @param[in]	face		The face of the domain cell underlying the interface cell to be defined
	 *
	 * @param[in]	material_id	The types::material_id which is assigned to the interface cell
	 * 							(thereby deciding to which interface portion the interface cell
	 * 							is linked)
	 */
	void
	add_interface_cell(	const DomainCell&			domain_cell,
						const unsigned int			face,
						const types::material_id	material_id);

	/**
	 * %Function to add several cells to the interface triangulation at once. This function mainly exists
	 * for compatibility with a premature version of the library. Essentially,
	 * TriangulationSystem::add_interface_cell() is called repeatedly for all elements of the
	 * input vector @p interface_cells.
	 *
	 * @param[in]	interface_cells		A vector defining the interface cell
	 * 									(first in tuple: 	the domain cell underlying the interface cell
	 * 														to be defined;
	 * 									 second in tuple: 	the face of the domain cell underlying the
	 * 									 					interface cell to be defined;
	 * 									 third in tuple:	the types::material_id which is
	 * 									 					assigned to the interface cell
	 * 									);
	 * 									see also TriangulationSystem::add_interface_cell()
	 */
	void
	add_interface_cells(std::vector< std::tuple<const DomainCell, const unsigned int, const types::material_id> > interface_cells);

	/**
	 * Attach a manifold to the interface (note that the types::manifold_id%s of interface cells
	 * will generally be taken from the types::manifold_id of the coarse mesh domain cell faces
	 * underlying the interface cells on the minus side of the interface). By attaching appropriate
	 * manifolds to the interface, it can be ensured that the refinement of the interface triangulation
	 * stays geometrically consistent with the refinement of the underlying domain cells. To achieve
	 * the latter, the codim 1 equivalents of the manifolds used by the domain triangulation must be added
	 * with the same types::manifold_id%s.
	 *
	 * This method is necessary because there is currently no way to directly obtain the codim 1
	 * equivalent from a manifold used on the domain.
	 *
	 * @param[in]	manifold_id 	The types::manifold_id with which @p manifold will be associated
	 *
	 * @param[in]	manifold		Reference to the codim 1 manifold to be used
	 */
	void
	set_interface_manifold(	const types::manifold_id				manifold_id,
							const Manifold<spacedim-1, spacedim>&	manifold );

	/**
	 * %Function to add several manifolds to the interface triangulation at once. This function mainly exists
	 * for compatibility with a premature version of the library. Essentially,
	 * TriangulationSystem::set_interface_manifold() is called repeatedly for all elements of the
	 * input vector @p manifolds.
	 *
	 * @param[in]	manifolds	A map between types::manifold_id%s and corresponding manifold objects
	 */
	void
	set_interface_manifolds(const std::map<types::manifold_id, std::reference_wrapper<const Manifold<spacedim-1,spacedim>>> manifolds);

	/**
	 * @return	An iterator to the first active interface cell. Note that this is not
	 * 			a deal.II iterator. Rather, it is an iterator to the first element
	 * 			in TriangulationSystem::active_interface_cell_domain_cells.
	 * 			However, through the returned iterator, deal.II iterators to the
	 * 			interface cells and adjacent domain cells can be obtained.
	 */
	typename std::vector<InterfaceCellDomainCells<spacedim>>::iterator
	interface_begin_active();

	/**
	 * @return	An iterator to past the last active interface cell. Note that this is not
	 * 			a deal.II iterator. Rather, it is an iterator to the past the end element
	 * 			in TriangulationSystem::active_interface_cell_domain_cells.
	 */
	typename std::vector<InterfaceCellDomainCells<spacedim>>::iterator
	interface_end_active();

	/**
	 * @return	This returns a const reference to TriangulationSystem::active_interface_cell_domain_cells,
	 * 			which is essentially be meant for range based loops instead of using the functions
	 * 			TriangulationSystem::interface_begin_active() and
	 * 			TriangulationSystem::interface_end_active() (hence the name of the function).
	 */
	const typename std::vector<InterfaceCellDomainCells<spacedim>>&
	interface_active_iterators()
	const;

	/**
	 * @return	An iterator to the first coarse interface cell. Note that this is not
	 * 			a deal.II iterator. Rather, it is an iterator to the first element
	 * 			in TriangulationSystem::coarse_interface_cell_domain_cells.
	 * 			However, through the returned iterator, deal.II iterators to the
	 * 			interface cells and adjacent domain cells can be obtained.
	 */
	typename std::vector<InterfaceCellDomainCells<spacedim>>::iterator
	interface_begin_coarse();

	/**
	 * @return	An iterator to past the last coarse interface cell. Note that this is not
	 * 			a deal.II iterator. Rather, it is an iterator to the past the end element
	 * 			in TriangulationSystem::coarse_interface_cell_domain_cells.
	 */
	typename std::vector<InterfaceCellDomainCells<spacedim>>::iterator
	interface_end_coarse();

	/**
	 * @return	This returns a const reference to TriangulationSystem::coarse_interface_cell_domain_cells,
	 * 			which is essentially be meant for range based loops instead of using the functions
	 * 			TriangulationSystem::interface_begin_coarse() and
	 * 			TriangulationSystem::interface_end_coarse() (hence the name of the function).
	 */
	const typename std::vector<InterfaceCellDomainCells<spacedim>>&
	interface_coarse_iterators()
	const;

	/**
	 * A convenience method writing the domain and interface triangulations to vtk files.
	 * Note that the extension .vtk will not be added automatically. So, it has to be
	 * included in the file names.
	 *
	 * If a file name is left empty, the corresponding mesh will not be written.
	 *
	 * @param[in]	file_name_domain	file name for the domain triangulation
	 *
	 * @param[in]	file_name_interface	file name for the interface triangulation
	 */
	void
	write_triangulations_vtk(	const std::string file_name_domain,
								const std::string file_name_interface)
	const;

	/**
	 * %Function for global mesh refinement. This method does nothing but calling
	 * the deal.II function refine_global() on the domain mesh, which will in
	 * turn trigger refinement of the interface by TriangulationSystem::pre_refinement_domain()
	 * and TriangulationSystem::post_refinement_domain().
	 */
	void
	refine_global(const unsigned int times = 1);

	/**
	 * %Function to execute coarsening and refinement.
	 * This method does nothing but calling the deal.II function execute_coarsening_and_refinement()
	 * on the domain mesh, which will in turn trigger refinement of the interface by
	 * TriangulationSystem::pre_refinement_domain() and TriangulationSystem::post_refinement_domain().
	 *
	 * @todo  This function will also have to be involved in the transfer of hidden variables between meshes.
	 * 		  However, this is not implemented yet, meaning that hidden variables can currently not be combined with
	 * 		  mesh refinement during the computation.
	 */
	void
	execute_coarsening_and_refinement();

	/**
	 * This method checks that the active interface cells are properly aligned with the
	 * underlying domain cell faces by comparison of cell / face centers.
	 *
	 * @param[in]	tol		The distance between the center of an interface cell
	 * 						and the center of the underlying domain cell face below
	 * 						which the centers are considered consistent.
	 *
	 * @return				If @p true is returned, all is in order and the interface
	 * 						mesh is consistent with the domain mesh. If @p false is returned,
	 * 						inconsistencies have been detected.
	 */
	bool
	check_active_interface_cell_domain_cells_consistency(const double tol = 1e-12)
	const;

	/**
	 * @return	The number of this processor and the total number of participating processors
	 * 			(as this is a sequential triangulation, the result is always (0, 1)
	 */
	virtual
	std::pair<const unsigned int, const unsigned int>
	get_this_proc_n_procs()
	const;
};

#ifdef DEAL_II_WITH_P4EST
#ifdef DEAL_II_WITH_MPI
namespace parallel
{

/**
 * This class derives from parallel::distributed::Triangulation and provides an interface triangulation appropriate for
 * use within the parallel::distributed::TriangulationSystem.
 * The main additional functionality of this class over parallel::distributed::Triangulation is that the subdomain ids
 * can be assigned manually. The triangulation can only work properly, if the partitioning is done consistently.
 */
template<unsigned int dim, unsigned int spacedim = dim>
class Triangulation : public dealii::parallel::DistributedTriangulationBase<dim, spacedim>
{
protected:

	virtual
	void
	update_cell_relations()
	override;

	virtual
	void
	save(const std::string& filename)
	const
	override;

	virtual
	void
	load(const std::string& filename)
	override;

	virtual
	void
	load(	const std::string&	filename,
		const bool		autopartition)
	override;

public:
	Triangulation(MPI_Comm mpi_communicator);

	/**
	 * This function must be called after any change of the subdomain id%s in order
	 * to make sure that cached variables are updated
	 */
	void
	finalize_subdomain_assignment();

	/**
	 * Return if multilevel hierarchy is supported and has been constructed (here always false).
	 */
	virtual
	bool
	is_multilevel_hierarchy_constructed()
	const;
};

/**
 * This is the parallel computing equivalent of the class TriangulationSystem (for general information
 * refer to the documentation of this class).
 *
 * For the domain triangulation the parallel::distributed::Triangulation class is used. In contrast,
 * for the interface triangulation the parallel::distributed::Triangulation class is used, which allows for the
 * necessary manual management of the partitioning of the triangulation.
 *
 * Regarding ownership of interface cells, the convention used is that the ownership of interface cells is
 * according to the ownership of the underlying domain cells on the minus side of the interface.
 *
 * In general, the complete coarse domain mesh with all types::material_id%s, types::manifold_id%s must be identically
 * defined on each processor (this ensures that no ids are lost during the repartitioning process after mesh refinement).
 *
 * @tparam	spacedim	The spatial dimension of the problem.
 */
template<unsigned int spacedim>
class TriangulationSystem : public dealii::GalerkinTools::TriangulationSystem<spacedim>
{

private:

	/**
	 * MPI communicator
	 */
	MPI_Comm
	mpi_communicator;

	/**
	 * In parallel, this function does nothing except triggering the TriangulationSystem::pre_refinement() signal.
	 *
	 * @todo This function will have to be involved in the transfer of hidden variables between meshes.
	 * 		 However, this is not implemented yet, meaning that hidden variables can currently not be combined with
	 * 		 mesh refinement during the computation.
	 */
	virtual
	void
	pre_refinement_domain();

	/**
	 * This internal function takes care that the interface mesh is updated after any refinement/coarsening of the domain cells.
	 * The underlying mechanism is that the method is attached to the domain triangulation during
	 * construction of the TriangulationSystem; and it will always be called immediately after the domain mesh has been
	 * refined.
	 *
	 * @todo  This function will also have to be involved in the transfer of hidden variables between meshes.
	 * 		  However, this is not implemented yet, meaning that hidden variables can currently not be combined with
	 * 		  mesh refinement during the computation.
	 */
	virtual
	void
	post_refinement_domain();

	/**
	 * This function assigns the interface subdomain ids according to the subdomain ids of
	 * the underlying domain cells on the minus side.
	 */
	void
	update_interface_subdomain_ids()
	const;

public:
	/**
	 * Constructor of the parallel::distributed::TriangulationSystem. Note that the constructor does not copy the @p tria_domain
	 * object supplied, but rather stores a pointer to it. So do not change it after construction of the
	 * parallel::distributed::TriangulationSystem unless you're knowing exactly what you are doing (the only change to the
	 * domain triangulation, which is admissible, is mesh refinement). No checking on external changes
	 * of the domain triangulation is currently performed and, therefore, doing so may lead to errors
	 * which are hard to debug.
	 *
	 * @param[in]	tria_domain				TriangulationSystem::tria_domain
	 *
	 * @param[in]	fix_vertex_positions	TriangulationSystem::fix_vertex_positions
	 */
	TriangulationSystem(dealii::parallel::distributed::Triangulation<spacedim, spacedim>&	tria_domain,
						const bool															fix_vertex_positions = false);

		/**
	 * Destructor
	 */
	virtual
	~TriangulationSystem() = default;

	/**
	 * @return	This returns a const reference to the domain triangulation.
	 */
	virtual
	const dealii::parallel::distributed::Triangulation<spacedim, spacedim>&
	get_triangulation_domain()
	const;

	/**
	 * @return	This returns a const reference to the interface triangulation.
	 */
	virtual
	const dealii::GalerkinTools::parallel::Triangulation<spacedim-1, spacedim>&
	get_triangulation_interface()
	const;

	/**
	 * @return	This returns a reference to the domain triangulation.
	 */
	virtual
	dealii::parallel::distributed::Triangulation<spacedim, spacedim>&
	get_triangulation_domain();

	/**
	 * @return	This returns a reference to the interface triangulation.
	 */
	virtual
	dealii::GalerkinTools::parallel::Triangulation<spacedim-1, spacedim>&
	get_triangulation_interface();

	/**
	 * %Function closing the interface definition. After this function has been called, no modifications to
	 * the interface triangulation are possible. This function does the real work of creating the interface
	 * triangulation and associating it with the underlying domain triangulation.
	 */
	virtual
	void
	close();

	/**
	 * Write triangulations in *.vtu format for each processor, and add a .pvtu file for visualization
	 * in Visit or Paraview that describes the collection of *.vtu files.
	 *
	 * Note that the extension .vtu will be added automatically.
	 *
	 * If a file name is left empty, the corresponding meshes will not be written.
	 *
	 * @param[in]	file_name_domain	file name for the domain triangulation
	 *
	 * @param[in]	file_name_interface	file name for the interface triangulation
	 */
	void
	write_meshes_per_processor_as_vtu(	const std::string file_name_domain,
										const std::string file_name_interface)
	const;

	/**
	 * @return	The number of this processor and the total number of participating processors
	 */
	virtual
	std::pair<const unsigned int, const unsigned int>
	get_this_proc_n_procs()
	const;

};

}
#endif // DEAL_II_WITH_MPI
#endif // DEAL_II_WITH_P4EST

GALERKIN_TOOLS_NAMESPACE_CLOSE

DEAL_II_NAMESPACE_CLOSE

#endif /* GALERKINTOOLS_TRIANGULATIONSYSTEM_H_ */
