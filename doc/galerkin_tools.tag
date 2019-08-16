<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile>
  <compound kind="file">
    <name>mainpage.dox</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/doc/</path>
    <filename>mainpage_8dox</filename>
  </compound>
  <compound kind="file">
    <name>assembly_helper.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>assembly__helper_8h</filename>
    <class kind="class">AssemblyHelper</class>
  </compound>
  <compound kind="file">
    <name>config.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>config_8h</filename>
  </compound>
  <compound kind="file">
    <name>dependent_field.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>dependent__field_8h</filename>
    <class kind="class">DependentFieldTerm</class>
    <class kind="class">DependentField</class>
    <class kind="class">DependentField&lt; spacedim, spacedim &gt;</class>
  </compound>
  <compound kind="file">
    <name>dirichlet_constraint.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>dirichlet__constraint_8h</filename>
    <class kind="class">DirichletConstraint</class>
  </compound>
  <compound kind="file">
    <name>dof_handler_system.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>dof__handler__system_8h</filename>
    <class kind="class">DoFHandlerSystem</class>
    <class kind="class">InterfaceCellDoFIterator</class>
    <class kind="class">DomainCellDoFIterator</class>
    <class kind="class">InterfaceCellDomainCellsDoF</class>
    <class kind="class">DoFHandlerSystem</class>
  </compound>
  <compound kind="file">
    <name>dof_renumbering.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>dof__renumbering_8h</filename>
    <class kind="class">DoFRenumbering</class>
    <class kind="class">DoFRenumberingOffset</class>
  </compound>
  <compound kind="file">
    <name>fe_values_interface.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>fe__values__interface_8h</filename>
    <class kind="class">FEValuesInterface</class>
  </compound>
  <compound kind="file">
    <name>independent_field.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>independent__field_8h</filename>
    <class kind="class">IndependentField</class>
    <class kind="class">IndependentField&lt; 0, spacedim &gt;</class>
  </compound>
  <compound kind="file">
    <name>ldr.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>ldr_8h</filename>
    <namespace>Auxiliary</namespace>
    <member kind="function">
      <type>int</type>
      <name>compute_ldr</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>abe796c1529d11eed08fd05bca82f3002</anchor>
      <arglist>(FullMatrix&lt; double &gt; &amp;C, Vector&lt; double &gt; &amp;D, std::vector&lt; Vector&lt; double &gt;&gt; &amp;L, std::vector&lt; Vector&lt; double &gt;&gt; &amp;R)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>linear_material.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>linear__material_8h</filename>
    <class kind="class">LinearMaterialDomain</class>
    <class kind="class">LinearMaterialInterface</class>
  </compound>
  <compound kind="file">
    <name>scalar_functional.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>scalar__functional_8h</filename>
    <class kind="class">ScalarFunctional</class>
    <class kind="class">ScalarFunctional&lt; spacedim, spacedim &gt;</class>
  </compound>
  <compound kind="file">
    <name>solver_wrapper.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>solver__wrapper_8h</filename>
    <class kind="class">SolverWrapper</class>
    <class kind="class">SolverWrapperUMFPACK</class>
    <class kind="class">BlockSolverWrapperUMFPACK</class>
    <class kind="class">SolverWrapperPETSc</class>
  </compound>
  <compound kind="file">
    <name>tools.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>tools_8h</filename>
    <namespace>Auxiliary</namespace>
    <member kind="function">
      <type>void</type>
      <name>convert_local_indices_to_global_indices</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>ad3c42d209f0ba8e9d4ce305060634bf1</anchor>
      <arglist>(const std::vector&lt; unsigned int &gt; &amp;dof_indices_local, std::vector&lt; unsigned int &gt; &amp;dof_indices_global, const std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>combine_dof_indices</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a1d90ebc8738df3d8c70b540034137019</anchor>
      <arglist>(const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_interface, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_minus, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_plus, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_C, std::vector&lt; unsigned int &gt; &amp;dof_indices_interface_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_minus_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_plus_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_C_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_global_combined)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_to_index_set</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>aca5fe5966aef4ba2293fb44095bfd86c</anchor>
      <arglist>(const DoFRenumbering &amp;dof_renumbering, const IndexSet &amp;in, IndexSet &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>renumber_constraints</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>aa6148bcbaf5e3003717b4dd2a4da15b3</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix, const DoFRenumbering &amp;dof_renumbering=DoFRenumbering(), const bool close=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_dof_renumbering_contiguous</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a4261fc726ff965166f3f53242918acea</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system, DoFRenumberingOffset &amp;dof_renumbering_offset)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_map_dofs</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a072f85e6d745ae3c532bb0612f4bd3ce</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dhs_1, const DoFHandlerSystem&lt; spacedim &gt; &amp;dhs_2, std::vector&lt; unsigned int &gt; &amp;map_dofs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>split_vector</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a59c01a6511bffd7442693e86cd194ef1</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;in, Vector&lt; double &gt; &amp;out_0, Vector&lt; double &gt; &amp;out_1, const unsigned int size_1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>split_matrix</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>af746a1d08b1135a3684fc990f7b1384d</anchor>
      <arglist>(const FullMatrix&lt; double &gt; &amp;in, FullMatrix&lt; double &gt; &amp;out_00, FullMatrix&lt; double &gt; &amp;out_01, FullMatrix&lt; double &gt; &amp;out_10, FullMatrix&lt; double &gt; &amp;out_11, const unsigned int size_1)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>communicate_bool</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>ad0ff28386be7b54b3487ae36e5a074fa</anchor>
      <arglist>(const bool local_bool, const MPI_Comm &amp;mpi_communicator)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>total_potential.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>total__potential_8h</filename>
    <class kind="class">TotalPotential</class>
  </compound>
  <compound kind="file">
    <name>total_potential_contribution.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>total__potential__contribution_8h</filename>
    <class kind="class">TotalPotentialContribution</class>
  </compound>
  <compound kind="file">
    <name>triangulation_system.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>triangulation__system_8h</filename>
    <class kind="class">InterfaceCellDomainCells</class>
    <class kind="class">TriangulationSystem</class>
    <class kind="class">parallel::TriangulationSystem</class>
    <namespace>parallel</namespace>
    <member kind="enumeration">
      <type></type>
      <name>InterfaceRefinementCase</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a4cfb8c5e21535951e919b6a6b1023af7</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>minus_is_finer</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a4cfb8c5e21535951e919b6a6b1023af7a725ee72f9a8d77b879deec24e24a1d7b</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>plus_is_finer</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a4cfb8c5e21535951e919b6a6b1023af7aee243748f7ee1b44506873ec0b8ae191</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>equally_fine</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a4cfb8c5e21535951e919b6a6b1023af7a37477769ebb40217b8b2dcb406f8b564</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>at_boundary</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a4cfb8c5e21535951e919b6a6b1023af7ad453630f99a3f3cec0aaea2d520ac3cf</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumeration">
      <type></type>
      <name>InterfaceSide</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a44f3c00e36c1d6e3c389ae693c09b435</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>minus</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a44f3c00e36c1d6e3c389ae693c09b435aba34d25d1a0e7f60052d7418d1ccf58b</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>plus</name>
      <anchorfile>triangulation__system_8h.html</anchorfile>
      <anchor>a44f3c00e36c1d6e3c389ae693c09b435a772311abb4290d65c05c49d38c4da6fb</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>two_block_matrix.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>two__block__matrix_8h</filename>
    <class kind="class">TwoBlockMatrix</class>
    <class kind="class">parallel::TwoBlockMatrix</class>
    <namespace>parallel</namespace>
  </compound>
  <compound kind="file">
    <name>two_block_sparsity_pattern.h</name>
    <path>/home/sst/FE/code/incrementalFE/galerkin_tools/include/galerkin_tools/</path>
    <filename>two__block__sparsity__pattern_8h</filename>
    <class kind="class">TwoBlockSparsityPattern</class>
  </compound>
  <compound kind="class">
    <name>AssemblyHelper</name>
    <filename>class_assembly_helper.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>AssemblyHelper</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae164448dcc5e9e8d2849c354212d6df5</anchor>
      <arglist>(const TotalPotential&lt; spacedim &gt; &amp;total_potential, TriangulationSystem&lt; spacedim &gt; &amp;tria_system, const Mapping&lt; spacedim, spacedim &gt; &amp;mapping_domain, const Mapping&lt; spacedim-1, spacedim &gt; &amp;mapping_interface, const std::set&lt; const IndependentField&lt; 0, spacedim &gt; * &gt; &amp;independent_scalars=std::set&lt; const IndependentField&lt; 0, spacedim &gt; * &gt;())</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~AssemblyHelper</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac82eca7b04aedf772028c6ff77245e9b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_initial_fields_vector</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ade14ba73257f862f4e9cde8d0d6df12e</anchor>
      <arglist>(VectorType &amp;initial_fields, const AffineConstraints&lt; double &gt; *constraints=nullptr) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>make_dirichlet_constraints</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab9602b90beefa27eba50d9d837c1bf7b</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix, const std::vector&lt; const DirichletConstraint&lt; spacedim &gt; * &gt; &amp;dirichlet_constraints, const AffineConstraints&lt; double &gt; &amp;constraints_ignore=AffineConstraints&lt; double &gt;()) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>generate_sparsity_pattern_by_simulation</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af45c789076262ceb24613c04427b654e</anchor>
      <arglist>(SparsityPatternType &amp;dsp_K, const AffineConstraints&lt; double &gt; &amp;constraints) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>assemble_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9d2f3d1152046e639acfa6f6ac317b37</anchor>
      <arglist>(const SolutionVectorType &amp;solution, const std::vector&lt; const SolutionVectorType * &gt; solution_ref_sets, const AffineConstraints&lt; double &gt; &amp;constraints, double &amp;potential_value, RHSVectorType &amp;f, MatrixType &amp;K, const std::tuple&lt; bool, bool, bool &gt; requested_quantities=std::make_tuple(true, true, true)) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_nonprimitive_scalar_functional_values</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a43bac43f3aecf18e08bfd1741cf30af6</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, std::map&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; *, double &gt; &amp;nonprimitive_scalar_functional_values_domain, std::map&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; *, double &gt; &amp;nonprimitive_scalar_functional_values_interface) const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_maximum_step_length</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aae924a79482fed55899a59052bd9103d</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, const VectorType &amp;delta_solution) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compare_derivatives_with_numerical_derivatives</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac87490e0d3d84b11e2f7187f8706dab6</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, const std::string detailed_printout_file=&quot;&quot;, const double epsilon=1e-8) const </arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const std::string, const std::string &gt;</type>
      <name>write_output_independent_fields</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a74187aa98464043ea6572c3ac345640e</anchor>
      <arglist>(const VectorType &amp;solution, const std::string file_name_domain, const std::string file_name_interface, const unsigned int file_index=0, const std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt; &amp;dp_domain=std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt;(), const std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt; &amp;dp_interface=std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt;(), const unsigned int n_subdivisions=1) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_assembly_helper_definition</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>afd598b93397e6af53f0e4e274e6f880e</anchor>
      <arglist>(const bool detailed_printout_shapefuns=true) const </arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const double, const double &gt;</type>
      <name>compute_distance_to_other_solution</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac7860831588d35d05009474f2a695e14</anchor>
      <arglist>(const VectorType &amp;solution, const VectorType &amp;other_solution, const AssemblyHelper&lt; spacedim &gt; &amp;other_assembly_helper, const Quadrature&lt; spacedim &gt; quadrature_domain, const Quadrature&lt; spacedim-1 &gt; quadrature_interface, const VectorTools::NormType norm_type=VectorTools::NormType::L2_norm, const ComponentMask component_mask_domain=ComponentMask(), const ComponentMask component_mask_interface=ComponentMask(), const double exponent=2.0) const </arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const double, const double &gt;</type>
      <name>compute_distance_to_exact_solution</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aca82c19b1abbc316d0bc563e04db727c</anchor>
      <arglist>(const VectorType &amp;solution, const Function&lt; spacedim &gt; &amp;exact_solution_domain, const Function&lt; spacedim &gt; &amp;exact_solution_interface, const Quadrature&lt; spacedim &gt; quadrature_domain, const Quadrature&lt; spacedim-1 &gt; quadrature_interface, const VectorTools::NormType norm_type=VectorTools::NormType::L2_norm, const ComponentMask component_mask_domain=ComponentMask(), const ComponentMask component_mask_interface=ComponentMask(), const double exponent=2.0) const </arglist>
    </member>
    <member kind="function">
      <type>const TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>get_triangulation_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a261ecb9213338856aa88c8ae60a44c78</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>get_triangulation_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a16f9d21a79922d4879e37916b414f7d0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>get_dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab8e95f0469f1595ab00c3fb48bcaf4fd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>get_dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a04523eef6062ced8c88d4c093b65df3d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; const IndependentField&lt; spacedim, spacedim &gt; *, const unsigned int &gt;</type>
      <name>get_u_omega_global_component_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a396b89981e546af6b9bc0e35634290b3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; const IndependentField&lt; spacedim-1, spacedim &gt; *, const unsigned int &gt;</type>
      <name>get_u_sigma_global_component_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aca4d34e08f177e8d075c86bb34906f2f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_u_omega_global_component_index</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a41bdbb21e3f77cf717c9f7465363e415</anchor>
      <arglist>(const IndependentField&lt; spacedim, spacedim &gt; &amp;u_omega) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_u_sigma_global_component_index</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9a8e7a9f29b275dc20811cac001bd18f</anchor>
      <arglist>(const IndependentField&lt; spacedim-1, spacedim &gt; &amp;u_sigma) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>system_size</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a7de6972444e41dadb8eaac8024b261f6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_n_stretched_rows</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5cd3242e02bc5cb8b74cf808df257da0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_n_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad6590021b351fac59dcc655e3d0da9ee</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_global_dof_index_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a035fabd9601baf9efade5393164ea370</anchor>
      <arglist>(const IndependentField&lt; 0, spacedim &gt; *independent_scalar) const </arglist>
    </member>
    <member kind="function">
      <type>const IndexSet</type>
      <name>get_locally_relevant_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae6c72b5ed3b1cd419d58e081562e0ee7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const IndexSet</type>
      <name>get_locally_owned_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a51d99905072b1e6d1aadc43e62c5af92</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; IndexSet &gt;</type>
      <name>get_locally_relevant_indices_blocks</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae405978288614436b0851e1d9047f084</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; IndexSet &gt;</type>
      <name>get_locally_owned_indices_blocks</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1d0898b738b49d1ed38448d5686e19ba</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_dof_index_at_point_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a53c369b25d3a595229a9834950200da7</anchor>
      <arglist>(const IndependentField&lt; spacedim, spacedim &gt; *u_omega, const unsigned int component, const Point&lt; spacedim &gt; p) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_dof_index_at_point_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8c3efbac750aa236f8b116af994c07ee</anchor>
      <arglist>(const IndependentField&lt; spacedim-1, spacedim &gt; *u_sigma, const unsigned int component, const Point&lt; spacedim &gt; p) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_dof_information</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a532f565a8725675fcca12c1f8c669a44</anchor>
      <arglist>(const unsigned int dof_index) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>convert_dependent_fields_to_shapefunctions</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aad21ced11a2d90c804827854c18f7f89</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_hidden_variables</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1cd5d35b136347876aa13b89749338a0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>distribute_dofs</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a953859cc4cd745a1b51fdec5418be682</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_fe_values_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab4c51b059d53f0086ebce3b60d4d2400</anchor>
      <arglist>(const typename hp::DoFHandler&lt; spacedim, spacedim &gt;::active_cell_iterator &amp;cell, const unsigned int internal_index, const bool nonprimitive=false) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_fe_values_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a07b7b359ed67e0f949aa8817853ad0c2</anchor>
      <arglist>(const InterfaceCellDomainCellsDoF&lt; spacedim &gt; &amp;interface_cell_domain_cells, const unsigned int internal_index, const bool nonprimitive=false) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>compute_e_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a78e22b83e1effe9b40d26ed150bebc7f</anchor>
      <arglist>(const unsigned int internal_index, const unsigned int scalar_functional_index, const unsigned int q_point, const Vector&lt; double &gt; &amp;solution_u_omega, const Vector&lt; double &gt; &amp;solution_C, Vector&lt; double &gt; &amp;e_omega, FullMatrix&lt; double &gt; &amp;de_omega_dsol_T, const bool compute_derivative=true, const bool ignore_constants=false) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>compute_e_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a17cff5d9e32bd2ceac43d9218edd0b9b</anchor>
      <arglist>(const unsigned int internal_index, const unsigned int scalar_functional_index, const unsigned int q_point, const Vector&lt; double &gt; &amp;solution_u_sigma, const Vector&lt; double &gt; &amp;solution_u_omega_minus, const Vector&lt; double &gt; &amp;solution_u_omega_plus, const Vector&lt; double &gt; &amp;solution_C, const std::vector&lt; unsigned int &gt; &amp;dof_indices_interface_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_minus_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_plus_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_C_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_combined, Vector&lt; double &gt; &amp;e_sigma, FullMatrix&lt; double &gt; &amp;de_sigma_dsol_T, const bool compute_derivative=true, const bool ignore_constants=false) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>get_nonprimitive_scalar_functional_values</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5ab3659137ee74754cc27d9c18e1f9db</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, Vector&lt; double &gt; &amp;nonprimitive_scalar_functional_values) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>std::pair&lt; const int, const int &gt;</type>
      <name>get_scalar_functional_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aad795811df548677cb883341f6c52001</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim, spacedim &gt; *scalar_functional) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>std::pair&lt; const int, const int &gt;</type>
      <name>get_scalar_functional_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae5791ac73405d479df790815abc06380</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim-1, spacedim &gt; *scalar_functional) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>get_dof_indices_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4cc746ddb7917fa0e9f7cdda05345b94</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;global_dof_indices_C) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>make_dirichlet_constraints_recursion</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5c4d81be4c30b8a40761ac5710734576</anchor>
      <arglist>(const typename TriangulationSystem&lt; spacedim &gt;::DomainCell &amp;domain_cell, const unsigned int face, const std::vector&lt; unsigned int &gt; &amp;shapefuns, const DirichletConstraint&lt; spacedim &gt; &amp;constraint, AffineConstraints&lt; double &gt; &amp;constraint_matrix, const AffineConstraints&lt; double &gt; &amp;constraints_ignore) const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const TotalPotential&lt; spacedim &gt;</type>
      <name>total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a748eed9d73b7437a4bf2dcd73108790b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>tria_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>add08a8a7bb9c9325fcc7d92bfce525d4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const SmartPointer&lt; const Mapping&lt; spacedim, spacedim &gt; &gt;</type>
      <name>mapping_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3fbb49461000dea8f64266f830709fad</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const SmartPointer&lt; const Mapping&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>mapping_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a055fde6217c18e62cd80188d0130c201</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, const unsigned int &gt;</type>
      <name>material_id_to_internal_index_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a10b3acf64bccc169ee14dc2505ce4b46</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; std::tuple&lt; const types::material_id, const types::material_id, const types::material_id &gt;, const unsigned int &gt;</type>
      <name>material_ids_to_internal_index_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a43b82de0ede96d03b9f7fd8740d81668</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>hp::FECollection&lt; spacedim, spacedim &gt;</type>
      <name>fe_collection_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af3803b0aad9853e6bf018c70be41e791</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>hp::FECollection&lt; spacedim-1, spacedim &gt;</type>
      <name>fe_collection_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8b4d224a9ecd2e926a8860829874d2a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, unsigned int &gt;</type>
      <name>material_id_to_fe_system_id_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3045f80801fc31920efd161a268aae8e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, unsigned int &gt;</type>
      <name>material_id_to_fe_system_id_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fea54137e3c1c5a514e39c9b2ad7926</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DoFHandlerSystem&lt; spacedim &gt;</type>
      <name>dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a885e660c749e91a35e3279643ebcd87f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; spacedim, spacedim &gt; &gt; &gt;</type>
      <name>u_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1952a054a839a7a683ca108013e7d976</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; spacedim-1, spacedim &gt; &gt; &gt;</type>
      <name>u_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a696fe649b3503561235aa1ccbf2ddeef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt; &gt;</type>
      <name>C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa5234a46be82cfe7d92678169d38f326</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_component_indices_u_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6dae4b6ae7934eaec1ad7baff258ce6e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_component_indices_u_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a992a53a1fcac8a393ca53fb8d504bdfe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_indices_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9a8f0e8ea8c67ce9429c16a2017cafdc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa6fa619e4c2582e95950e878cd06628e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt; &gt; &gt;</type>
      <name>scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a29aa77e0e8e6b35c94966ea88840e462</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; unsigned int &gt; &gt;</type>
      <name>scalar_functionals_domain_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fe78a019aec03cbeeb336d1d2874729</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; unsigned int &gt; &gt;</type>
      <name>scalar_functionals_interface_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a833383aa6d157157545204143897ed9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_domain_nonprimitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>acf05fab2ddf57769a103d82a4f2d1cd3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_interface_nonprimitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0d15b3ab0c7bec9fc4f40e532f8776f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_domain_primitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4f08790a2235e48ce19f5d8d965a7874</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_interface_primitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad99c75f32cf3f18aa1d4067ad8b56ae8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt; &gt;</type>
      <name>contributions_scalar_functionals_domain_total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aceaf7ba62dfe0fa06ecb15ee8c14da34</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt; &gt;</type>
      <name>contributions_scalar_functionals_interface_total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9e76874224ab4946218fdce9bdba0e03</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>n_scalar_functionals_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af7bcfc1db651535a7aefc6071a81e124</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>n_scalar_functionals_primitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af5e03e8e47a85dbc96444ef61525c454</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0bdb6e2e2f9623f3f10dfa2ebe8e234c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>abcad51e64347fc3141d2840a2835b46c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_domain_facewise</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae92560183f1d2060265f0744a84f0349</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1a26b40224e3f04e5168accc91486493</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab346e146cf91fb7a0688076551b37355</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a12299d82365553a21fef8529c8fe8a17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af07bb528fdd350e9b467b08dc44a03e7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_C_indices_scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0edd25820c92a25ae87fc240f4916804</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_C_indices_scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a311e176038ee2b7ca0719abb384ca57b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a13fef9096e5fc7b1e922ea78d7aa2c28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fbb532e798c2427af5285c2df10c9f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>c_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a75e6f76c0b12b91c5feb230251f0137f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; double &gt; &gt; &gt;</type>
      <name>d_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad93b109608d4425d318434e01cb6246c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa266cc07e9670319481da52d633d2583</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af58c9a1c7093edc306070913aa1b9be2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4461be378c9be0364ca23153c367d24c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6f51f8b4dfdae385a2e2fe2dd9e66cdb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8eb513a6239d94cc0f1a6a01a037c572</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab09dd07d3ec596525be83f8f5859bef7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>c_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab3d4370661ba726010ee687ef0e98140</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; double &gt; &gt; &gt;</type>
      <name>d_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a48d7d677120eb1c84b4983f470246e02</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a904a24f53b66e1c1ef89f1bb7989eb32</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae1e2643696005415d78421882ca80e8e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain_reinit</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a28b58551c6afc68c7beaaa2604bc6e92</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface_reinit</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac9789c5a00867744dd906b85580e3091</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain_reinit_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3c592ef0a148753891cc3e03fd08324c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface_reinit_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>afaa20027ee539ca8d9c40c317127e471</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::pair&lt; std::string, unsigned int &gt; &gt;</type>
      <name>component_names_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af5e29448f133863a1859be8bfbb300c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::pair&lt; std::string, unsigned int &gt; &gt;</type>
      <name>component_names_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a7ae6ae2ec356cbb7b830d968315d280c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; boost::signals2::connection &gt;</type>
      <name>tria_listeners</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a228cec028ab5126d25c3ebf0e12a17a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const unsigned int</type>
      <name>this_proc</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a2aad83ae1bfe5338794cf9b50848469a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const unsigned int</type>
      <name>n_procs</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a87945d87baf37637673fd124b3803fd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const TotalPotential&lt; spacedim &gt;</type>
      <name>total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a748eed9d73b7437a4bf2dcd73108790b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>tria_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>add08a8a7bb9c9325fcc7d92bfce525d4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const SmartPointer&lt; const Mapping&lt; spacedim, spacedim &gt; &gt;</type>
      <name>mapping_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3fbb49461000dea8f64266f830709fad</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const SmartPointer&lt; const Mapping&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>mapping_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a055fde6217c18e62cd80188d0130c201</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, const unsigned int &gt;</type>
      <name>material_id_to_internal_index_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a10b3acf64bccc169ee14dc2505ce4b46</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; std::tuple&lt; const types::material_id, const types::material_id, const types::material_id &gt;, const unsigned int &gt;</type>
      <name>material_ids_to_internal_index_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a43b82de0ede96d03b9f7fd8740d81668</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>hp::FECollection&lt; spacedim, spacedim &gt;</type>
      <name>fe_collection_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af3803b0aad9853e6bf018c70be41e791</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>hp::FECollection&lt; spacedim-1, spacedim &gt;</type>
      <name>fe_collection_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8b4d224a9ecd2e926a8860829874d2a1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, unsigned int &gt;</type>
      <name>material_id_to_fe_system_id_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3045f80801fc31920efd161a268aae8e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; types::material_id, unsigned int &gt;</type>
      <name>material_id_to_fe_system_id_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fea54137e3c1c5a514e39c9b2ad7926</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DoFHandlerSystem&lt; spacedim &gt;</type>
      <name>dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a885e660c749e91a35e3279643ebcd87f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; spacedim, spacedim &gt; &gt; &gt;</type>
      <name>u_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1952a054a839a7a683ca108013e7d976</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; spacedim-1, spacedim &gt; &gt; &gt;</type>
      <name>u_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a696fe649b3503561235aa1ccbf2ddeef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt; &gt;</type>
      <name>C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa5234a46be82cfe7d92678169d38f326</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_component_indices_u_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6dae4b6ae7934eaec1ad7baff258ce6e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_component_indices_u_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a992a53a1fcac8a393ca53fb8d504bdfe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>global_indices_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9a8f0e8ea8c67ce9429c16a2017cafdc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa6fa619e4c2582e95950e878cd06628e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt; &gt; &gt;</type>
      <name>scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a29aa77e0e8e6b35c94966ea88840e462</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; unsigned int &gt; &gt;</type>
      <name>scalar_functionals_domain_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fe78a019aec03cbeeb336d1d2874729</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; unsigned int &gt; &gt;</type>
      <name>scalar_functionals_interface_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a833383aa6d157157545204143897ed9e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_domain_nonprimitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>acf05fab2ddf57769a103d82a4f2d1cd3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_interface_nonprimitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0d15b3ab0c7bec9fc4f40e532f8776f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_domain_primitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4f08790a2235e48ce19f5d8d965a7874</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, const unsigned int &gt;</type>
      <name>scalar_functionals_interface_primitive_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad99c75f32cf3f18aa1d4067ad8b56ae8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; &gt;, std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt; &gt;</type>
      <name>contributions_scalar_functionals_domain_total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aceaf7ba62dfe0fa06ecb15ee8c14da34</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; SmartPointer&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; &gt;, std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt; &gt;</type>
      <name>contributions_scalar_functionals_interface_total_potential</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9e76874224ab4946218fdce9bdba0e03</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>n_scalar_functionals_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af7bcfc1db651535a7aefc6071a81e124</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>n_scalar_functionals_primitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af5e03e8e47a85dbc96444ef61525c454</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0bdb6e2e2f9623f3f10dfa2ebe8e234c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>abcad51e64347fc3141d2840a2835b46c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt; &gt;</type>
      <name>components_to_shapefuns_domain_facewise</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae92560183f1d2060265f0744a84f0349</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1a26b40224e3f04e5168accc91486493</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab346e146cf91fb7a0688076551b37355</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a12299d82365553a21fef8529c8fe8a17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_dof_indices_scalar_functionals_interface_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af07bb528fdd350e9b467b08dc44a03e7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_C_indices_scalar_functionals_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a0edd25820c92a25ae87fc240f4916804</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; unsigned int &gt; &gt; &gt;</type>
      <name>coupled_C_indices_scalar_functionals_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a311e176038ee2b7ca0719abb384ca57b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a13fef9096e5fc7b1e922ea78d7aa2c28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5fbb532e798c2427af5285c2df10c9f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>c_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a75e6f76c0b12b91c5feb230251f0137f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; double &gt; &gt; &gt;</type>
      <name>d_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad93b109608d4425d318434e01cb6246c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aa266cc07e9670319481da52d633d2583</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af58c9a1c7093edc306070913aa1b9be2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4461be378c9be0364ca23153c367d24c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_minus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a6f51f8b4dfdae385a2e2fe2dd9e66cdb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>a_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8eb513a6239d94cc0f1a6a01a037c572</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, const unsigned int, const unsigned int, std::vector&lt; unsigned int &gt; &gt; &gt; &gt; &gt; &gt;</type>
      <name>b_plus</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab09dd07d3ec596525be83f8f5859bef7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; std::vector&lt; std::tuple&lt; const double, unsigned int &gt; &gt; &gt; &gt; &gt;</type>
      <name>c_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab3d4370661ba726010ee687ef0e98140</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::vector&lt; double &gt; &gt; &gt;</type>
      <name>d_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a48d7d677120eb1c84b4983f470246e02</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a904a24f53b66e1c1ef89f1bb7989eb32</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::vector&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae1e2643696005415d78421882ca80e8e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain_reinit</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a28b58551c6afc68c7beaaa2604bc6e92</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface_reinit</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac9789c5a00867744dd906b85580e3091</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValues&lt; spacedim, spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_domain_reinit_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a3c592ef0a148753891cc3e03fd08324c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::set&lt; std::shared_ptr&lt; FEValuesInterface&lt; spacedim &gt; &gt; &gt; &gt;</type>
      <name>fe_values_interface_reinit_nonprimitive</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>afaa20027ee539ca8d9c40c317127e471</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::pair&lt; std::string, unsigned int &gt; &gt;</type>
      <name>component_names_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af5e29448f133863a1859be8bfbb300c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::pair&lt; std::string, unsigned int &gt; &gt;</type>
      <name>component_names_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a7ae6ae2ec356cbb7b830d968315d280c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; boost::signals2::connection &gt;</type>
      <name>tria_listeners</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a228cec028ab5126d25c3ebf0e12a17a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const unsigned int</type>
      <name>this_proc</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a2aad83ae1bfe5338794cf9b50848469a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const unsigned int</type>
      <name>n_procs</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a87945d87baf37637673fd124b3803fd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>convert_dependent_fields_to_shapefunctions</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aad21ced11a2d90c804827854c18f7f89</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_hidden_variables</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1cd5d35b136347876aa13b89749338a0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>distribute_dofs</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a953859cc4cd745a1b51fdec5418be682</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_fe_values_domain</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab4c51b059d53f0086ebce3b60d4d2400</anchor>
      <arglist>(const typename hp::DoFHandler&lt; spacedim, spacedim &gt;::active_cell_iterator &amp;cell, const unsigned int internal_index, const bool nonprimitive=false) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>initialize_fe_values_interface</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a07b7b359ed67e0f949aa8817853ad0c2</anchor>
      <arglist>(const InterfaceCellDomainCellsDoF&lt; spacedim &gt; &amp;interface_cell_domain_cells, const unsigned int internal_index, const bool nonprimitive=false) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>compute_e_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a78e22b83e1effe9b40d26ed150bebc7f</anchor>
      <arglist>(const unsigned int internal_index, const unsigned int scalar_functional_index, const unsigned int q_point, const Vector&lt; double &gt; &amp;solution_u_omega, const Vector&lt; double &gt; &amp;solution_C, Vector&lt; double &gt; &amp;e_omega, FullMatrix&lt; double &gt; &amp;de_omega_dsol_T, const bool compute_derivative=true, const bool ignore_constants=false) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>compute_e_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a17cff5d9e32bd2ceac43d9218edd0b9b</anchor>
      <arglist>(const unsigned int internal_index, const unsigned int scalar_functional_index, const unsigned int q_point, const Vector&lt; double &gt; &amp;solution_u_sigma, const Vector&lt; double &gt; &amp;solution_u_omega_minus, const Vector&lt; double &gt; &amp;solution_u_omega_plus, const Vector&lt; double &gt; &amp;solution_C, const std::vector&lt; unsigned int &gt; &amp;dof_indices_interface_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_minus_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_plus_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_C_dof_indices_combined, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_combined, Vector&lt; double &gt; &amp;e_sigma, FullMatrix&lt; double &gt; &amp;de_sigma_dsol_T, const bool compute_derivative=true, const bool ignore_constants=false) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>get_nonprimitive_scalar_functional_values</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5ab3659137ee74754cc27d9c18e1f9db</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, Vector&lt; double &gt; &amp;nonprimitive_scalar_functional_values) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>std::pair&lt; const int, const int &gt;</type>
      <name>get_scalar_functional_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aad795811df548677cb883341f6c52001</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim, spacedim &gt; *scalar_functional) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>std::pair&lt; const int, const int &gt;</type>
      <name>get_scalar_functional_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae5791ac73405d479df790815abc06380</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim-1, spacedim &gt; *scalar_functional) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>get_dof_indices_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a4cc746ddb7917fa0e9f7cdda05345b94</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;global_dof_indices_C) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>make_dirichlet_constraints_recursion</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5c4d81be4c30b8a40761ac5710734576</anchor>
      <arglist>(const typename TriangulationSystem&lt; spacedim &gt;::DomainCell &amp;domain_cell, const unsigned int face, const std::vector&lt; unsigned int &gt; &amp;shapefuns, const DirichletConstraint&lt; spacedim &gt; &amp;constraint, AffineConstraints&lt; double &gt; &amp;constraint_matrix, const AffineConstraints&lt; double &gt; &amp;constraints_ignore) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_initial_fields_vector</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ade14ba73257f862f4e9cde8d0d6df12e</anchor>
      <arglist>(VectorType &amp;initial_fields, const AffineConstraints&lt; double &gt; *constraints=nullptr) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>make_dirichlet_constraints</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab9602b90beefa27eba50d9d837c1bf7b</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix, const std::vector&lt; const DirichletConstraint&lt; spacedim &gt; * &gt; &amp;dirichlet_constraints, const AffineConstraints&lt; double &gt; &amp;constraints_ignore=AffineConstraints&lt; double &gt;()) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>generate_sparsity_pattern_by_simulation</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>af45c789076262ceb24613c04427b654e</anchor>
      <arglist>(SparsityPatternType &amp;dsp_K, const AffineConstraints&lt; double &gt; &amp;constraints) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>assemble_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9d2f3d1152046e639acfa6f6ac317b37</anchor>
      <arglist>(const SolutionVectorType &amp;solution, const std::vector&lt; const SolutionVectorType * &gt; solution_ref_sets, const AffineConstraints&lt; double &gt; &amp;constraints, double &amp;potential_value, RHSVectorType &amp;f, MatrixType &amp;K, const std::tuple&lt; bool, bool, bool &gt; requested_quantities=std::make_tuple(true, true, true)) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_nonprimitive_scalar_functional_values</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a43bac43f3aecf18e08bfd1741cf30af6</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, std::map&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; *, double &gt; &amp;nonprimitive_scalar_functional_values_domain, std::map&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; *, double &gt; &amp;nonprimitive_scalar_functional_values_interface) const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_maximum_step_length</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aae924a79482fed55899a59052bd9103d</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, const VectorType &amp;delta_solution) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compare_derivatives_with_numerical_derivatives</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac87490e0d3d84b11e2f7187f8706dab6</anchor>
      <arglist>(const VectorType &amp;solution, const std::vector&lt; const VectorType * &gt; solution_ref_sets, const std::string detailed_printout_file=&quot;&quot;, const double epsilon=1e-8) const </arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const std::string, const std::string &gt;</type>
      <name>write_output_independent_fields</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a74187aa98464043ea6572c3ac345640e</anchor>
      <arglist>(const VectorType &amp;solution, const std::string file_name_domain, const std::string file_name_interface, const unsigned int file_index=0, const std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt; &amp;dp_domain=std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt;(), const std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt; &amp;dp_interface=std::vector&lt; dealii::SmartPointer&lt; const dealii::DataPostprocessor&lt; spacedim &gt;&gt;&gt;(), const unsigned int n_subdivisions=1) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_assembly_helper_definition</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>afd598b93397e6af53f0e4e274e6f880e</anchor>
      <arglist>(const bool detailed_printout_shapefuns=true) const </arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const double, const double &gt;</type>
      <name>compute_distance_to_other_solution</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ac7860831588d35d05009474f2a695e14</anchor>
      <arglist>(const VectorType &amp;solution, const VectorType &amp;other_solution, const AssemblyHelper&lt; spacedim &gt; &amp;other_assembly_helper, const Quadrature&lt; spacedim &gt; quadrature_domain, const Quadrature&lt; spacedim-1 &gt; quadrature_interface, const VectorTools::NormType norm_type=VectorTools::NormType::L2_norm, const ComponentMask component_mask_domain=ComponentMask(), const ComponentMask component_mask_interface=ComponentMask(), const double exponent=2.0) const </arglist>
    </member>
    <member kind="function">
      <type>std::pair&lt; const double, const double &gt;</type>
      <name>compute_distance_to_exact_solution</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aca82c19b1abbc316d0bc563e04db727c</anchor>
      <arglist>(const VectorType &amp;solution, const Function&lt; spacedim &gt; &amp;exact_solution_domain, const Function&lt; spacedim &gt; &amp;exact_solution_interface, const Quadrature&lt; spacedim &gt; quadrature_domain, const Quadrature&lt; spacedim-1 &gt; quadrature_interface, const VectorTools::NormType norm_type=VectorTools::NormType::L2_norm, const ComponentMask component_mask_domain=ComponentMask(), const ComponentMask component_mask_interface=ComponentMask(), const double exponent=2.0) const </arglist>
    </member>
    <member kind="function">
      <type>const TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>get_triangulation_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a261ecb9213338856aa88c8ae60a44c78</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>TriangulationSystem&lt; spacedim &gt; &amp;</type>
      <name>get_triangulation_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a16f9d21a79922d4879e37916b414f7d0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>get_dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ab8e95f0469f1595ab00c3fb48bcaf4fd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>get_dof_handler_system</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a04523eef6062ced8c88d4c093b65df3d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; const IndependentField&lt; spacedim, spacedim &gt; *, const unsigned int &gt;</type>
      <name>get_u_omega_global_component_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a396b89981e546af6b9bc0e35634290b3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::map&lt; const IndependentField&lt; spacedim-1, spacedim &gt; *, const unsigned int &gt;</type>
      <name>get_u_sigma_global_component_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>aca4d34e08f177e8d075c86bb34906f2f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_u_omega_global_component_index</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a41bdbb21e3f77cf717c9f7465363e415</anchor>
      <arglist>(const IndependentField&lt; spacedim, spacedim &gt; &amp;u_omega) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_u_sigma_global_component_index</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a9a8e7a9f29b275dc20811cac001bd18f</anchor>
      <arglist>(const IndependentField&lt; spacedim-1, spacedim &gt; &amp;u_sigma) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>system_size</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a7de6972444e41dadb8eaac8024b261f6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_n_stretched_rows</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a5cd3242e02bc5cb8b74cf808df257da0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_n_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ad6590021b351fac59dcc655e3d0da9ee</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_global_dof_index_C</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a035fabd9601baf9efade5393164ea370</anchor>
      <arglist>(const IndependentField&lt; 0, spacedim &gt; *independent_scalar) const </arglist>
    </member>
    <member kind="function">
      <type>const IndexSet</type>
      <name>get_locally_relevant_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae6c72b5ed3b1cd419d58e081562e0ee7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const IndexSet</type>
      <name>get_locally_owned_indices</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a51d99905072b1e6d1aadc43e62c5af92</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; IndexSet &gt;</type>
      <name>get_locally_relevant_indices_blocks</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>ae405978288614436b0851e1d9047f084</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; IndexSet &gt;</type>
      <name>get_locally_owned_indices_blocks</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a1d0898b738b49d1ed38448d5686e19ba</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_dof_index_at_point_omega</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a53c369b25d3a595229a9834950200da7</anchor>
      <arglist>(const IndependentField&lt; spacedim, spacedim &gt; *u_omega, const unsigned int component, const Point&lt; spacedim &gt; p) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_dof_index_at_point_sigma</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a8c3efbac750aa236f8b116af994c07ee</anchor>
      <arglist>(const IndependentField&lt; spacedim-1, spacedim &gt; *u_sigma, const unsigned int component, const Point&lt; spacedim &gt; p) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_dof_information</name>
      <anchorfile>class_assembly_helper.html</anchorfile>
      <anchor>a532f565a8725675fcca12c1f8c669a44</anchor>
      <arglist>(const unsigned int dof_index) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BlockSolverWrapperUMFPACK</name>
    <filename>class_block_solver_wrapper_u_m_f_p_a_c_k.html</filename>
    <base>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::BlockVector&lt; double &gt;, TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt; &gt;, TwoBlockSparsityPattern &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_block_solver_wrapper_u_m_f_p_a_c_k.html</anchorfile>
      <anchor>acb9f3b8c9a7996cad23a4d5a663432d7</anchor>
      <arglist>(const TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt;&gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::BlockVector&lt; double &gt; &amp;f_stretched, const bool symmetric=false) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DependentField</name>
    <filename>class_dependent_field.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>DependentField</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a5a47b8f5b0e51b29bf6f3ca6e18fda9c</anchor>
      <arglist>(const std::string name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a0726cab14197769ad07bd8f071e395ea</anchor>
      <arglist>(double coefficient, const IndependentField&lt; dim, spacedim &gt; &amp;independent_field, const unsigned int component=0, const std::vector&lt; unsigned int &gt; derivatives=std::vector&lt; unsigned int &gt;())</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a852ed477f6ee71ad85a6690d4a5e56fe</anchor>
      <arglist>(double coefficient, const IndependentField&lt; dim+1, spacedim &gt; &amp;independent_field, const unsigned int component, const std::vector&lt; unsigned int &gt; derivatives, const InterfaceSide side)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>acb8979c891a45ed7326eb47df41c0129</anchor>
      <arglist>(double coefficient, const IndependentField&lt; dim+1, spacedim &gt; &amp;independent_field, const unsigned int component, const InterfaceSide side)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a1aaaf072334d4e06640383ae81766a39</anchor>
      <arglist>(double coefficient, const IndependentField&lt; dim, spacedim &gt; &amp;independent_field, const unsigned int component, const unsigned int derivative)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a43de305294f94c38c88fa4768b3f0155</anchor>
      <arglist>(double coefficient, const IndependentField&lt; dim+1, spacedim &gt; &amp;independent_field, const unsigned int component, const unsigned int derivative, const InterfaceSide side)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>adf9e4863ca954a9fa70c57aef68849f8</anchor>
      <arglist>(double coefficient, const IndependentField&lt; 0, spacedim &gt; &amp;independent_field)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a72196bff18c2c28a1e811877e739f4ee</anchor>
      <arglist>(double constant)</arglist>
    </member>
    <member kind="function">
      <type>const std::set&lt; DependentFieldTerm&lt; dim, spacedim &gt; &gt; &amp;</type>
      <name>get_terms_interface</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a5fab03134453b1dd4999b90b50a74480</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::set&lt; DependentFieldTerm&lt; dim+1, spacedim &gt; &gt; &amp;</type>
      <name>get_terms_neighbor</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a3cd8757a280d4ef9232cb304a68ea6f4</anchor>
      <arglist>(const InterfaceSide side) const </arglist>
    </member>
    <member kind="function">
      <type>const std::set&lt; DependentFieldTerm&lt; 0, spacedim &gt; &gt; &amp;</type>
      <name>get_terms_independent_scalars</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a7e17fcef1db1822cac8c12212e3e8755</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_constant</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a9f255eb459c862a9cbe785084430445b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::set&lt; const IndependentField&lt; dim, spacedim &gt; * &gt;</type>
      <name>get_independent_fields_interface</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>acc7f03e2157ecf784540ec6d17250476</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::set&lt; const IndependentField&lt; dim+1, spacedim &gt; * &gt;</type>
      <name>get_independent_fields_neighbors</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a80a0bcfebadff71fdf10fec84097c821</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::set&lt; const IndependentField&lt; 0, spacedim &gt; * &gt;</type>
      <name>get_independent_scalars</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>ab0bff11e19a2d7035ace9f010b833f14</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a3f1fd18437bdbb4c9d83b9b841ed242d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a698b255f131279d1edb97b3c525153a0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; dim, spacedim &gt; &gt;</type>
      <name>terms_interface</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>aa6326d64fd7828c935946136ec5ec122</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; dim+1, spacedim &gt; &gt;</type>
      <name>terms_neighbor_plus</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a05ba04a9f0f50fd881055f3abf46c9f5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; dim+1, spacedim &gt; &gt;</type>
      <name>terms_neighbor_minus</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a9ac1d390137d52ef6c58b230d540ce38</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; 0, spacedim &gt; &gt;</type>
      <name>terms_independent_scalars</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>aa1ba265eae76489d22d3eecd595931c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>constant</name>
      <anchorfile>class_dependent_field.html</anchorfile>
      <anchor>a32b37c78e04a16b6b606442f156c8ca9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DependentField&lt; spacedim, spacedim &gt;</name>
    <filename>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>DependentField</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a3b2090f3b7c66d84dbb41ec5e677ecfe</anchor>
      <arglist>(const std::string name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a9e0b93f36ef8bc10ee9028576db57390</anchor>
      <arglist>(double coefficient, const IndependentField&lt; spacedim, spacedim &gt; &amp;independent_field, const unsigned int component=0, const std::vector&lt; unsigned int &gt; derivatives=std::vector&lt; unsigned int &gt;())</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a79f3b134f21759492662b2f7173fff8e</anchor>
      <arglist>(double coefficient, const IndependentField&lt; spacedim, spacedim &gt; &amp;independent_field, const unsigned int component, const unsigned int derivative)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a25b4485644ccc8e2aeafddceadb33950</anchor>
      <arglist>(double coefficient, const IndependentField&lt; 0, spacedim &gt; &amp;independent_field)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_term</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a24fe77c93937f9563b3dac9e4520fd9b</anchor>
      <arglist>(double constant)</arglist>
    </member>
    <member kind="function">
      <type>const std::set&lt; DependentFieldTerm&lt; spacedim, spacedim &gt; &gt; &amp;</type>
      <name>get_terms_domain</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>abe2d6beee95fae305add085afa7b6119</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const std::set&lt; DependentFieldTerm&lt; 0, spacedim &gt; &gt; &amp;</type>
      <name>get_terms_independent_scalars</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a24de9a650280f55800c5b22a1017f60d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_constant</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>acede82b94119d1a86f08030dbded477e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::set&lt; const IndependentField&lt; spacedim, spacedim &gt; * &gt;</type>
      <name>get_independent_fields_domain</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ac5527c93285b39f3adc70709243448f0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::set&lt; const IndependentField&lt; 0, spacedim &gt; * &gt;</type>
      <name>get_independent_scalars</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a793e428507763e02a79fb4c6773ada7b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>afa7738e59a3352eb4242d6b31fe3ff5f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a99a47f4c10f0e472dbdc2b0945712e23</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; spacedim, spacedim &gt; &gt;</type>
      <name>terms_domain</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ace23fb05e85b4c967793e8ffce58f332</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::set&lt; DependentFieldTerm&lt; 0, spacedim &gt; &gt;</type>
      <name>terms_independent_scalars</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a54061701af0eea90b02e56612f1491cc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>constant</name>
      <anchorfile>class_dependent_field_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a24efc9c0928896be871908f050c406fc</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DependentFieldTerm</name>
    <filename>class_dependent_field_term.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>DependentFieldTerm</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>a149a09f602b0baadba9033d6a8da448c</anchor>
      <arglist>(const double coefficient, const IndependentField&lt; dim, spacedim &gt; &amp;independent_field, const unsigned int component=0, const std::vector&lt; unsigned int &gt; derivatives=std::vector&lt; unsigned int &gt;())</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>first_derivative</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>a9c477ae87acdea43960f399ca2dd10b9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_derivatives</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>aa1c7aeb391135d4cafa1f208c15517ec</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator&lt;</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>ab9934d7ad41e3a52073cc11293e7178b</anchor>
      <arglist>(const DependentFieldTerm &amp;dependent_field_2) const </arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>coefficient</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>a59a9183a32ac55fb728f3797b68a9f8f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const IndependentField&lt; dim, spacedim &gt; &gt;</type>
      <name>independent_field</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>a89d1c3fea36e6fe105232097a321e095</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>component</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>ac6f3ac40d4ee2c8b9f9bbdfa34079b74</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; unsigned int &gt;</type>
      <name>derivatives</name>
      <anchorfile>class_dependent_field_term.html</anchorfile>
      <anchor>af09c5452c3e8e71e9ee99db304b90135</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DirichletConstraint</name>
    <filename>class_dirichlet_constraint.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>DirichletConstraint</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>a14e75a7f8bcc8a18daa04a5c1164fe87</anchor>
      <arglist>(const IndependentField&lt; spacedim, spacedim &gt; &amp;independent_field, const unsigned int component, const InterfaceSide side, const std::set&lt; types::material_id &gt; domain_of_constraint, const Function&lt; spacedim &gt; *const constraint_inhomogeneity=nullptr, const IndependentField&lt; 0, spacedim &gt; *independent_scalar=nullptr, const Function&lt; spacedim &gt; *const coefficient_c=nullptr)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~DirichletConstraint</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>ac6b6b7caf4694000486abff377fd1440</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const std::set&lt; types::material_id &gt;</type>
      <name>domain_of_constraint</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>a258b6ff11b206f966bb03943bb11f469</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const IndependentField&lt; spacedim, spacedim &gt; &gt;</type>
      <name>independent_field</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>abbd7945a973ed93d1d773307393ffde3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>component</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>a7e3c4d0e0906af1c81b88e05e41bdafc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const InterfaceSide</type>
      <name>side</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>ae049d107664d3bf23d287ec77545b6f3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>constraint_inhomogeneity</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>af22d2bca23999bdb4f2b67e8982d29c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const IndependentField&lt; 0, spacedim &gt; &gt;</type>
      <name>independent_scalar</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>a8793f0d41a9c6e88638d7cf5d52a6fdd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>coefficient_c</name>
      <anchorfile>class_dirichlet_constraint.html</anchorfile>
      <anchor>adea2dc6126a633b297eb7dfce832b9ef</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DoFHandlerSystem</name>
    <filename>class_do_f_handler_system.html</filename>
    <templarg>spacedim</templarg>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>InterfaceCell</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a64bfb9dab5f61c2583896a7d086cb1f2</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>hp::DoFHandler&lt; spacedim-1, spacedim &gt;::cell_iterator</type>
      <name>InterfaceCellDoF</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a82e65f8c260f076489bd92615221a1d0</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim, spacedim &gt; &gt;</type>
      <name>DomainCell</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a84df955a5d8da7e1439a10f2a64811c0</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>hp::DoFHandler&lt; spacedim, spacedim &gt;::cell_iterator</type>
      <name>DomainCellDoF</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a01c54960bba8faa935d32f2d3d2fe914</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DoFHandlerSystem</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>acb31af36eec0abac1d6e764e3d08c7fe</anchor>
      <arglist>(const TriangulationSystem&lt; spacedim &gt; &amp;tria_system)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~DoFHandlerSystem</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a5a2f0389f4f4216b9841d9c518cb107e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>distribute_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ac6e0950ae9d140b0a9185bcfffe167d6</anchor>
      <arglist>(const hp::FECollection&lt; spacedim, spacedim &gt; &amp;fe_collection_domain, const hp::FECollection&lt; spacedim-1, spacedim &gt; &amp;fe_collection_interface, const unsigned int n_additional_dofs=0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_fe</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a1cbcc00826aca94736d3aebef157da72</anchor>
      <arglist>(const hp::FECollection&lt; spacedim, spacedim &gt; &amp;fe_collection_domain, const hp::FECollection&lt; spacedim-1, spacedim &gt; &amp;fe_collection_interface)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCellsDoF&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_begin_active</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a88d93ac03defd1b296238de4eb122841</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCellsDoF&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_end_active</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aecbc22a86943c1bc921a13ec7e2333f2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; InterfaceCellDomainCellsDoF&lt; spacedim &gt; &gt; &amp;</type>
      <name>interface_active_iterators</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ae22105dc8db1090f12be9ecaed384228</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DomainCellDoFIterator&lt; spacedim &gt;</type>
      <name>domain_begin_active</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aeee77ce076a84afdae132910c72de3b4</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>DomainCellDoFIterator&lt; spacedim &gt;</type>
      <name>domain_end_active</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>afb624d9bed38a49be5ebbf16dd2abcd5</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>IteratorRange&lt; DomainCellDoFIterator&lt; spacedim &gt; &gt;</type>
      <name>domain_active_iterators</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a01759f9734ca7e762685d7c85fc30cad</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_dofs_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a69aca46be1a419d4ae3b66d57bed8b8e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_dofs_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ac9edbf2dd12e85e83f9a945bd43065ab</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_dofs_additional</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a516048c499fac349be1546daaca38d89</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a62460181681846997229c994c9df9a7a</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a0934c0d1d8ee44f1d5ad8771716db0e7</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;dof_indices) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_dof_index</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a3b13d23ab9e956c71529f415a9dbcd09</anchor>
      <arglist>(const unsigned int &amp;dof_index) const </arglist>
    </member>
    <member kind="function">
      <type>const hp::DoFHandler&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_dof_handler_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aa57297d063cb453c85d4c1135d0847b0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const hp::DoFHandler&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_dof_handler_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a445c8ac85aada74eb741dd5a07d9a231</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>hp::DoFHandler&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_dof_handler_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a47b35524f131eb242ade5a674096be90</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>hp::DoFHandler&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_dof_handler_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a7e006440056a4a246c2bdf14347c62d9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const IndexSet &amp;</type>
      <name>get_locally_owned_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a4f4bced9fc691a16e97648856b42577d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const IndexSet &amp;</type>
      <name>get_locally_relevant_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a15566db3dbb5d3ab2e4732354432f56e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>attach_dof_renumbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ab31e887efe0c3703e3a777287c250ed0</anchor>
      <arglist>(const DoFRenumbering &amp;dof_renumbering)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>make_hanging_node_constraints</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a8708acd1db19ae6965e885ebbe11d262</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>split_vector</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>abe00061acc73a319e705939ed965018d</anchor>
      <arglist>(const VectorType &amp;in_vect, VectorType &amp;out_vect_domain, VectorType &amp;out_vect_interface, VectorType &amp;out_vect_C) const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; unsigned int &gt; &amp;</type>
      <name>get_n_dofs_per_processor</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a741d168ac52d3591d687736a7a4f6dd6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_single_dof_index_component_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a0acf8d35c183ac73bb445e3bac3dd59a</anchor>
      <arglist>(const unsigned int component) const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_single_dof_index_component_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aa25553114a5cbe59607b8b9df2162c3e</anchor>
      <arglist>(const unsigned int component) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>update_interface_domain_relation</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a237b112398864f03ef97873f276960f5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>set_locally_owned_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ab24852414baf7b391b67d21644bd66f8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>set_locally_relevant_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>affa92a0bb132688295b748e0c5436d2e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>set_locally_owned_dofs_standard_numbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a56c9cd6d17075c1b7a6d24d5d28357b2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>set_locally_relevant_dofs_standard_numbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a51f57afbe35de6bdd39a979911eaa838</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>split_vector_implementation</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a5f235d3e4fa9dbd19cd541defdbaeacc</anchor>
      <arglist>(const LinearAlgebra::distributed::Vector&lt; double &gt; &amp;in_vect, LinearAlgebra::distributed::Vector&lt; double &gt; &amp;out_vect, const unsigned int window_begin, const unsigned int window_end) const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>split_vector_implementation</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a62a1c2f3a805f985203e5e131e9a547b</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;in_vect, Vector&lt; double &gt; &amp;out_vect, const unsigned int window_begin, const unsigned int window_end) const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const SmartPointer&lt; const TriangulationSystem&lt; spacedim &gt; &gt;</type>
      <name>tria_system</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a06d93193cb47591db138cd8f41953796</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::shared_ptr&lt; hp::DoFHandler&lt; spacedim, spacedim &gt; &gt;</type>
      <name>dof_handler_domain</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ac3c43d8113395b0011179231ff6c58aa</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::shared_ptr&lt; hp::DoFHandler&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>dof_handler_interface</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aa9480c1fcf0d9170c026ef6611074d06</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>dof_indices_C</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a15c27ca30c5905b44691779c755cd69c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>locally_owned_dof_indices_C</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a7b2e77b2c718b0b13d861da7c0530b28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; InterfaceCellDomainCellsDoF&lt; spacedim &gt; &gt;</type>
      <name>active_interface_cell_domain_cells</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>af0119b14377300f7f0457a18bb7dcd67</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>SmartPointer&lt; const DoFRenumbering &gt;</type>
      <name>dof_renumbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>aef6159c606a24ac7daadcd3fe082b3b6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>IndexSet</type>
      <name>locally_owned_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>ad72a701a3581187eec846c831ba384c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>IndexSet</type>
      <name>locally_relevant_dofs</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a6b2de0e80cf0e67e62391cc84f8a7be5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; unsigned int &gt;</type>
      <name>n_dofs_per_processor</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>acfe853db91a67d1f78e08d62d91bab1b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>IndexSet</type>
      <name>locally_owned_dofs_standard_numbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a257ddb680d9f8276f7b802884dd8e577</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>IndexSet</type>
      <name>locally_relevant_dofs_standard_numbering</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a0e6744c40aada76d695f2c564f4c5040</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; boost::signals2::connection &gt;</type>
      <name>tria_listeners</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a286776a935dacb3d1c79a90f4ca7e5d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>InterfaceCellDoFIterator</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a266cb42b2fa4af49d4df0d938f962749</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>DomainCellDoFIterator</name>
      <anchorfile>class_do_f_handler_system.html</anchorfile>
      <anchor>a22fa60ad60906aacbbe21d3b5704ebfc</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DoFRenumbering</name>
    <filename>class_do_f_renumbering.html</filename>
    <base>Subscriptor</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>convert_dof_indices</name>
      <anchorfile>class_do_f_renumbering.html</anchorfile>
      <anchor>ad37a052986c62df79a1ebc8abea13e89</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;dof_indices) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt;</type>
      <name>convert_range</name>
      <anchorfile>class_do_f_renumbering.html</anchorfile>
      <anchor>ae95ba2f7f04a172812f086fc75b61c15</anchor>
      <arglist>(const unsigned int range_begin, const unsigned int range_end) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~DoFRenumbering</name>
      <anchorfile>class_do_f_renumbering.html</anchorfile>
      <anchor>a0443d7cb614fa161e4201ff7a6bdf541</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DoFRenumberingOffset</name>
    <filename>class_do_f_renumbering_offset.html</filename>
    <base>DoFRenumbering</base>
    <member kind="function">
      <type>void</type>
      <name>add_range</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>a39333b4c3597c2b4a967f7af3cec0054</anchor>
      <arglist>(const unsigned int dof_start, const unsigned int dof_end, const int offset)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>convert_dof_indices</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>a7b48f6f59b90c015b7e176148edbb797</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;dof_indices) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::vector&lt; std::pair&lt; const unsigned int, const unsigned int &gt; &gt;</type>
      <name>convert_range</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>ac066233b202f7982d29206c1221c1f6a</anchor>
      <arglist>(const unsigned int range_begin, const unsigned int range_end) const </arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; std::tuple&lt; const unsigned int, const unsigned int, const int &gt; &gt; &amp;</type>
      <name>get_dof_offsets</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>af92a21719157c5fd5e8ae42cf8af87f6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clear</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>ac7156c58d0aa5dcb57fa0fe3ac437453</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>a51ae1b70be0b0a2444c4bba71a272f71</anchor>
      <arglist>(std::ostream &amp;out) const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; std::tuple&lt; const unsigned int, const unsigned int, const int &gt; &gt;</type>
      <name>dof_offsets</name>
      <anchorfile>class_do_f_renumbering_offset.html</anchorfile>
      <anchor>a5df6c4b70b1394c3670ced634146c9a9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>DomainCellDoFIterator</name>
    <filename>class_domain_cell_do_f_iterator.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>DomainCellDoFIterator</name>
      <anchorfile>class_domain_cell_do_f_iterator.html</anchorfile>
      <anchor>a36416883d4c91dd0849520d3e205dba1</anchor>
      <arglist>(const TriaIterator&lt; CellAccessor&lt; spacedim, spacedim &gt;&gt; &amp;domain_cell, const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DomainCellDoFIterator</name>
      <anchorfile>class_domain_cell_do_f_iterator.html</anchorfile>
      <anchor>a4e987bc73bb2c8ec6974cc6427505091</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices</name>
      <anchorfile>class_domain_cell_do_f_iterator.html</anchorfile>
      <anchor>ad93f1a26a4efd4ca94af45d13e250aba</anchor>
      <arglist>(std::vector&lt; types::global_dof_index &gt; &amp;dof_indices) const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>dof_handler_system</name>
      <anchorfile>class_domain_cell_do_f_iterator.html</anchorfile>
      <anchor>a7431a0505f1e5b0b796c7f58cf9a060c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>FEValuesInterface</name>
    <filename>class_f_e_values_interface.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>FEValuesInterface</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>ae48863eabc0754e1cc033143fcd2d443</anchor>
      <arglist>(const Mapping&lt; spacedim-1, spacedim &gt; &amp;mapping_interface, const Mapping&lt; spacedim, spacedim &gt; &amp;mapping_domain, const FiniteElement&lt; spacedim-1, spacedim &gt; &amp;fe_interface, const FiniteElement&lt; spacedim, spacedim &gt; &amp;fe_domain_minus, const FiniteElement&lt; spacedim, spacedim &gt; &amp;fe_domain_plus, const Quadrature&lt; spacedim-1 &gt; &amp;quadrature, const UpdateFlags update_flags_interface, const UpdateFlags update_flags_domain_minus, const UpdateFlags update_flags_domain_plus)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a334ed0ef179516854690e7e3fa273aa0</anchor>
      <arglist>(const InterfaceCellDomainCellsDoF&lt; spacedim &gt; &amp;interface_cell_domain_cells)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a875a17d471069ddb97cbd8e6f7aba2b8</anchor>
      <arglist>(const InterfaceCellDomainCells&lt; spacedim &gt; &amp;interface_cell_domain_cells)</arglist>
    </member>
    <member kind="function">
      <type>const Quadrature&lt; spacedim-1 &gt; &amp;</type>
      <name>get_quadrature</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>afb7be190894c93ccf7ffdd4ff9bb6447</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const FEValues&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_fe_values_interface</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a6d70108a06e5604bac341818da2d8223</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const FEFaceValuesBase&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_fe_values_domain</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a8e222051f35763d2aec65d1840d2581a</anchor>
      <arglist>(const InterfaceSide &amp;interface_side) const </arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_quadrature_points</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a100afc348b15432f9e241ff56ebffa8d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>InterfaceRefinementCase</type>
      <name>initialized_interface_refinement_case</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a4a14b4fa181df533a6c7e7971fad5c4e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const FEFaceValuesBase&lt; spacedim, spacedim &gt; *</type>
      <name>initialized_fe_values_domain_minus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a48f1a887c14615d665408f156d879980</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const FEFaceValuesBase&lt; spacedim, spacedim &gt; *</type>
      <name>initialized_fe_values_domain_plus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>ab1f7524d517a79cd5e4adbd43017f415</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>FEValues&lt; spacedim-1, spacedim &gt;</type>
      <name>fe_values_interface</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a75f4ee97ebaefa2d1ded52c445953473</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>FEFaceValues&lt; spacedim, spacedim &gt;</type>
      <name>fe_face_values_domain_minus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a1ce4354bfe2e852fdc5b7917614404a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>FEFaceValues&lt; spacedim, spacedim &gt;</type>
      <name>fe_face_values_domain_plus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a6d6820b66a2694d327238c6b3c7fb30c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>FESubfaceValues&lt; spacedim, spacedim &gt;</type>
      <name>fe_subface_values_domain_minus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>abe003c242c807dd162c396d812bbfe77</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>FESubfaceValues&lt; spacedim, spacedim &gt;</type>
      <name>fe_subface_values_domain_plus</name>
      <anchorfile>class_f_e_values_interface.html</anchorfile>
      <anchor>a184e705efd975db536299ad619585a65</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IndependentField</name>
    <filename>class_independent_field.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>IndependentField</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a03c2bec3886170c288c59381886afd90</anchor>
      <arglist>(const std::string name, const FiniteElement&lt; dim, spacedim &gt; &amp;fe, const unsigned int n_components, const std::set&lt; types::material_id &gt; non_zero_regions, const Function&lt; spacedim &gt; *const initial_vals=nullptr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IndependentField</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>af953c9ecc5965d5bcbd61fab9e0b6e2b</anchor>
      <arglist>(const std::string name, const FiniteElement&lt; dim, spacedim &gt; &amp;fe, const std::set&lt; types::material_id &gt; non_zero_regions, const Function&lt; spacedim &gt; *const initial_vals=nullptr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IndependentField</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>af83ee8c1600bb079a3075244a7b39481</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>ae05f8565e4ce1a70b5b833555dc084b5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::unique_ptr&lt; const FiniteElement&lt; dim, spacedim &gt; &gt;</type>
      <name>fe</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a1c583665b7710bd3b815b03ba026b6d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_components</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a7b19ea8c30d72cf27f05669de61f30a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::set&lt; types::material_id &gt;</type>
      <name>non_zero_regions</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a4e09e114870c0b3761bc2e32916e5850</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>initial_vals</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a274c902785d2937a6065f7e09f3976c3</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IndependentField&lt; 0, spacedim &gt;</name>
    <filename>class_independent_field_3_010_00_01spacedim_01_4.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>IndependentField</name>
      <anchorfile>class_independent_field_3_010_00_01spacedim_01_4.html</anchorfile>
      <anchor>a8a486862088446941ea97fc471c5d44c</anchor>
      <arglist>(const std::string name, const double initial_value=0.0)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IndependentField</name>
      <anchorfile>class_independent_field_3_010_00_01spacedim_01_4.html</anchorfile>
      <anchor>a653c01f2a6aeffba44542d9788433647</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_independent_field_3_010_00_01spacedim_01_4.html</anchorfile>
      <anchor>a05ecdcc8310253f055fbc59abaa2bc90</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const double</type>
      <name>initial_value</name>
      <anchorfile>class_independent_field_3_010_00_01spacedim_01_4.html</anchorfile>
      <anchor>a8c4c434806dadc5b885102322e27357c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_components</name>
      <anchorfile>class_independent_field_3_010_00_01spacedim_01_4.html</anchorfile>
      <anchor>a248c2226570c71914f5219e7e3052561</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>IndependentField&lt; spacedim, spacedim &gt;</name>
    <filename>class_independent_field.html</filename>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>IndependentField</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a03c2bec3886170c288c59381886afd90</anchor>
      <arglist>(const std::string name, const FiniteElement&lt; dim, spacedim &gt; &amp;fe, const unsigned int n_components, const std::set&lt; types::material_id &gt; non_zero_regions, const Function&lt; spacedim &gt; *const initial_vals=nullptr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>IndependentField</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>af953c9ecc5965d5bcbd61fab9e0b6e2b</anchor>
      <arglist>(const std::string name, const FiniteElement&lt; dim, spacedim &gt; &amp;fe, const std::set&lt; types::material_id &gt; non_zero_regions, const Function&lt; spacedim &gt; *const initial_vals=nullptr)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~IndependentField</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>af83ee8c1600bb079a3075244a7b39481</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>ae05f8565e4ce1a70b5b833555dc084b5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::unique_ptr&lt; const FiniteElement&lt; dim, spacedim &gt; &gt;</type>
      <name>fe</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a1c583665b7710bd3b815b03ba026b6d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_components</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a7b19ea8c30d72cf27f05669de61f30a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::set&lt; types::material_id &gt;</type>
      <name>non_zero_regions</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a4e09e114870c0b3761bc2e32916e5850</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>initial_vals</name>
      <anchorfile>class_independent_field.html</anchorfile>
      <anchor>a274c902785d2937a6065f7e09f3976c3</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>InterfaceCellDoFIterator</name>
    <filename>class_interface_cell_do_f_iterator.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>InterfaceCellDoFIterator</name>
      <anchorfile>class_interface_cell_do_f_iterator.html</anchorfile>
      <anchor>a41a74a50c801c67a5393d568ed2768b7</anchor>
      <arglist>(const TriaIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt;&gt; &amp;interface_cell, const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InterfaceCellDoFIterator</name>
      <anchorfile>class_interface_cell_do_f_iterator.html</anchorfile>
      <anchor>a0dc67a5bf4cbf83e786af1fb28b6dd7c</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices</name>
      <anchorfile>class_interface_cell_do_f_iterator.html</anchorfile>
      <anchor>a7cd69c855243c0d8fb5df97911411edd</anchor>
      <arglist>(std::vector&lt; types::global_dof_index &gt; &amp;dof_indices) const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const DoFHandlerSystem&lt; spacedim &gt; &amp;</type>
      <name>dof_handler_system</name>
      <anchorfile>class_interface_cell_do_f_iterator.html</anchorfile>
      <anchor>a62045b61faee901392a7326b624d2b02</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>InterfaceCellDomainCells</name>
    <filename>class_interface_cell_domain_cells.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>InterfaceCell</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>ac6c8ada1cde14364575ef0191a278bf8</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim, spacedim &gt; &gt;</type>
      <name>DomainCell</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>aea8f8f65d0f5021da9c5b8e147b6d587</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InterfaceCellDomainCells</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a016b875566033b0da13d774081eaeaf4</anchor>
      <arglist>(const DomainCell &amp;domain_cell, const unsigned int face, const InterfaceCell &amp;interface_cell, const InterfaceSide interface_side)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~InterfaceCellDomainCells</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a7dae568b2d3c80762d0e1af52c93f279</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::tuple&lt; const types::material_id, const types::material_id, const types::material_id &gt;</type>
      <name>get_material_ids</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a9b5388da2f3a61d4340ae4d1b979d6d1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable">
      <type>const InterfaceCell</type>
      <name>interface_cell</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>acabc5a62be3f7742c3c334ae1777fbd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const DomainCell</type>
      <name>domain_cell_minus</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>ae610a6a7b0ff4a421ee051cb873bbc22</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>face_minus</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a51a9eaf54de991f5fa432220b09205cd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const DomainCell</type>
      <name>domain_cell_plus</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a72c7faaed3a84c546d47960f1064f9df</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>face_plus</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a8780275f79c7137ef65df9b7d4039017</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const InterfaceRefinementCase</type>
      <name>refinement_case</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>ab1b5469ca5c40256942ea179abbba92c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>subface</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>aba7f3e048c6985006988d716cb31aa7d</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type></type>
      <name>InterfaceCellDomainCells</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a4c0fd73a697833a1257b9341982897de</anchor>
      <arglist>(const std::tuple&lt; const InterfaceCell, const DomainCell, const unsigned int, const DomainCell, const unsigned int, const InterfaceRefinementCase, const unsigned int &gt; input)</arglist>
    </member>
    <member kind="function" protection="private" static="yes">
      <type>static std::tuple&lt; const InterfaceCell, const DomainCell, const unsigned int, const DomainCell, const unsigned int, const InterfaceRefinementCase, const unsigned int &gt;</type>
      <name>convert_constructor_inputs</name>
      <anchorfile>class_interface_cell_domain_cells.html</anchorfile>
      <anchor>a3f892a73f1f38f7412e549b9e62a8928</anchor>
      <arglist>(const DomainCell &amp;domain_cell, const unsigned int face, const InterfaceCell &amp;interface_cell, const InterfaceSide interface_side)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>InterfaceCellDomainCellsDoF</name>
    <filename>class_interface_cell_domain_cells_do_f.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>InterfaceCell</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a55783857b93a31baac5a6e6ee3dad22f</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>hp::DoFHandler&lt; spacedim-1, spacedim &gt;::cell_iterator</type>
      <name>InterfaceCellDoF</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a65f6a62d58c8378667cab3df94251ab6</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim, spacedim &gt; &gt;</type>
      <name>DomainCell</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a0d0b12a444adb340e1ffc2dd23d77902</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>hp::DoFHandler&lt; spacedim, spacedim &gt;::cell_iterator</type>
      <name>DomainCellDoF</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a464c96faf99349555566e4e469fdd34e</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>InterfaceCellDomainCellsDoF</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a79a88327d88717a2c3381fe7131b435d</anchor>
      <arglist>(const InterfaceCellDomainCells&lt; spacedim &gt; &amp;interface_cell_domain_cells, const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~InterfaceCellDomainCellsDoF</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a090e1536777928b78dff7fa7b1021a1e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::tuple&lt; const types::material_id, const types::material_id, const types::material_id &gt;</type>
      <name>get_material_ids</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>af2edacffb796130c2653eea32e6dc71f</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices_local_global_interface</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>af71f68397bb1a7f1ce028746b250a0d4</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global, std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global_minus, std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global_plus) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_dof_indices_local_global_interface</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>af4d1bb9c64294a83121282863ef0d207</anchor>
      <arglist>(std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global) const </arglist>
    </member>
    <member kind="variable">
      <type>const InterfaceCellDoFIterator&lt; spacedim &gt;</type>
      <name>interface_cell</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a5fe7922af91545598d89c77ba42a59da</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const DomainCellDoFIterator&lt; spacedim &gt;</type>
      <name>domain_cell_minus</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>ac8c84df4298d2d5be54dc8bf81e48c4b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>face_minus</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a01fc12497d44b06c7a325ea8f0b2ad28</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const DomainCellDoFIterator&lt; spacedim &gt;</type>
      <name>domain_cell_plus</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>a8fd9dbd3a1e53023ecca7917b004bab1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>face_plus</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>afcb7abea75dfc6e310493a6ff9aadcf6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const InterfaceRefinementCase</type>
      <name>refinement_case</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>aebb7e5f13d079fc83f98f67bcfcc6de3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>subface</name>
      <anchorfile>class_interface_cell_domain_cells_do_f.html</anchorfile>
      <anchor>ac0f33cbe60afdc4f5856ea39a085ed1a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LinearMaterialDomain</name>
    <filename>class_linear_material_domain.html</filename>
    <templarg>spacedim</templarg>
    <base>ScalarFunctional&lt; spacedim, spacedim &gt;</base>
    <member kind="function">
      <type></type>
      <name>LinearMaterialDomain</name>
      <anchorfile>class_linear_material_domain.html</anchorfile>
      <anchor>a7f91a980e9475691513443cd6b2efbc1</anchor>
      <arglist>(const std::vector&lt; DependentField&lt; spacedim, spacedim &gt;&gt; e_omega, const std::set&lt; types::material_id &gt; domain_of_integration, const Quadrature&lt; spacedim &gt; quadrature, const FullMatrix&lt; double &gt; C, const Vector&lt; double &gt; y, const std::string name=&quot;LinearMaterialDomain&quot;)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>get_h_omega</name>
      <anchorfile>class_linear_material_domain.html</anchorfile>
      <anchor>ac071e6886b4e661442af7a003d5f2b2a</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, double &amp;h_omega, Vector&lt; double &gt; &amp;h_omega_1, FullMatrix&lt; double &gt; &amp;h_omega_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const FullMatrix&lt; double &gt;</type>
      <name>C</name>
      <anchorfile>class_linear_material_domain.html</anchorfile>
      <anchor>a8ac95fcf4f77790670b8520ba9593a7c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const Vector&lt; double &gt;</type>
      <name>y</name>
      <anchorfile>class_linear_material_domain.html</anchorfile>
      <anchor>a7ea9ab6930c0b0aa826e809ef245b0e2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>LinearMaterialInterface</name>
    <filename>class_linear_material_interface.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>LinearMaterialInterface</name>
      <anchorfile>class_linear_material_interface.html</anchorfile>
      <anchor>ac10a5c148dd81edab41d22dd207c78a2</anchor>
      <arglist>(const std::vector&lt; DependentField&lt; spacedim-1, spacedim &gt;&gt; e_sigma, const std::set&lt; types::material_id &gt; domain_of_integration, const Quadrature&lt; spacedim-1 &gt; quadrature, const FullMatrix&lt; double &gt; C, const Vector&lt; double &gt; y, const std::string name=&quot;LinearMaterialInterface&quot;)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>get_h_sigma</name>
      <anchorfile>class_linear_material_interface.html</anchorfile>
      <anchor>a5c7730aa5b8b175950ba6ce670419450</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n, double &amp;h_sigma, Vector&lt; double &gt; &amp;h_sigma_1, FullMatrix&lt; double &gt; &amp;h_sigma_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const FullMatrix&lt; double &gt;</type>
      <name>C</name>
      <anchorfile>class_linear_material_interface.html</anchorfile>
      <anchor>ae0b1ef211453e5fd3e1416ff37c315db</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>const Vector&lt; double &gt;</type>
      <name>y</name>
      <anchorfile>class_linear_material_interface.html</anchorfile>
      <anchor>a3864513d7662e4d1c91e606467befd59</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>ScalarFunctional</name>
    <filename>class_scalar_functional.html</filename>
    <templarg>dim</templarg>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>ScalarFunctional</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>ac0618ada3b80400784ee8553a66aade9</anchor>
      <arglist>(const std::vector&lt; DependentField&lt; dim, spacedim &gt;&gt; e_sigma, const std::set&lt; types::material_id &gt; domain_of_integration, const Quadrature&lt; dim &gt; quadrature, const std::string name, const unsigned int n_ref_sets=0, const unsigned int n_hidden=0, const Function&lt; spacedim &gt; *const initial_vals_hidden=nullptr)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>get_h_sigma</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a1b9874b2fd591c844ecfcd1db8212c54</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n, double &amp;h_sigma, Vector&lt; double &gt; &amp;h_sigma_1, FullMatrix&lt; double &gt; &amp;h_sigma_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>get_maximum_step</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a31e4e48bb968b12fde76530c343b433d</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, const Vector&lt; double &gt; &amp;delta_e_sigma, const Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compare_derivatives_with_numerical_derivatives</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a95174310efe3b02d78ef47979846ff51</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_sigma, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_sigma_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const Tensor&lt; 1, spacedim &gt; &amp;n, const std::string detailed_printout_file=&quot;&quot;, const double epsilon=1e-8) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ScalarFunctional</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>aed36f35e6d2c9f9a93ee2749d01f2a51</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const std::set&lt; types::material_id &gt;</type>
      <name>domain_of_integration</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>ae3b6dd6934e1cd55fcc55cf344179407</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>initial_vals_hidden</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a602d0bc2c945822c6b756fc63183ae2b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; DependentField&lt; dim, spacedim &gt; &gt;</type>
      <name>e_sigma</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a86662b03a63219227993a2c6c07aefc1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const Quadrature&lt; dim &gt;</type>
      <name>quadrature</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>adea9ff214aeb2a1d8c3712a9d2433883</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_hidden</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a8b1617930242870f22eef5e306cb717f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a4d184688053b3443d10e228e4a8eba60</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_ref_sets</name>
      <anchorfile>class_scalar_functional.html</anchorfile>
      <anchor>a7e12423f4b29e9e0aaa0f7f9c2d1c0eb</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>ScalarFunctional&lt; spacedim, spacedim &gt;</name>
    <filename>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function">
      <type></type>
      <name>ScalarFunctional</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a4c4903b402908f26b8c8b433a7cbeb76</anchor>
      <arglist>(const std::vector&lt; DependentField&lt; spacedim, spacedim &gt;&gt; e_omega, const std::set&lt; types::material_id &gt; domain_of_integration, const Quadrature&lt; spacedim &gt; quadrature, const std::string name, const unsigned int n_ref_sets=0, const unsigned int n_hidden=0, const Function&lt; spacedim &gt; *const initial_vals_hidden=nullptr)</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual bool</type>
      <name>get_h_omega</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a629bfeae4d8ea364fc3f72fea8016ac8</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, double &amp;h_omega, Vector&lt; double &gt; &amp;h_omega_1, FullMatrix&lt; double &gt; &amp;h_omega_2, const std::tuple&lt; bool, bool, bool &gt; requested_quantities) const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>get_maximum_step</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>aba0e5304e9786bf28a25483d467e5d70</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, const Vector&lt; double &gt; &amp;delta_e_omega, const Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compare_derivatives_with_numerical_derivatives</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ab7f0d81df4bb8604f0ebc270b64d7f68</anchor>
      <arglist>(Vector&lt; double &gt; &amp;e_omega, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;e_omega_ref_sets, Vector&lt; double &gt; &amp;hidden_vars, const Point&lt; spacedim &gt; &amp;x, const std::string detailed_printout_file=&quot;&quot;, const double epsilon=1e-8) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~ScalarFunctional</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a42b56519b3e4338d248c1e29e29831a6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const std::set&lt; types::material_id &gt;</type>
      <name>domain_of_integration</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>aa192395f822a64f60df43bf9d36c2f3a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const SmartPointer&lt; const Function&lt; spacedim &gt; &gt;</type>
      <name>initial_vals_hidden</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ae3282d5182360e0030e4cc5e02fbe2eb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; DependentField&lt; spacedim, spacedim &gt; &gt;</type>
      <name>e_omega</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>adfed9b70b743ba245a39c3e63b951f96</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const Quadrature&lt; spacedim &gt;</type>
      <name>quadrature</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>ab83ee3ae077b211137824b006098382e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_hidden</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a7df6711471715f907bc9911449c5c825</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::string</type>
      <name>name</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>a195248af3821548af3000872e9e6d00e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const unsigned int</type>
      <name>n_ref_sets</name>
      <anchorfile>class_scalar_functional_3_01spacedim_00_01spacedim_01_4.html</anchorfile>
      <anchor>acee2c3c289e5b2b680996facc2f79e78</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapper</name>
    <filename>class_solver_wrapper.html</filename>
    <templarg>SolutionVectorType</templarg>
    <templarg>RHSVectorType</templarg>
    <templarg>MatrixType</templarg>
    <templarg>SparsityPatternType</templarg>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a75a599630086e2edee14eb43ae071733</anchor>
      <arglist>(const MatrixType &amp;K_stretched, SolutionVectorType &amp;solution, const RHSVectorType &amp;f_stretched, const bool symmetric) const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~SolverWrapper</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a738e65e298fb0b13327e853fc9683b6d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapper&lt; dealii::LinearAlgebra::distributed::Vector&lt; double &gt;, dealii::PETScWrappers::MPI::BlockVector, parallel::TwoBlockMatrix&lt; dealii::PETScWrappers::MPI::SparseMatrix &gt;, TwoBlockSparsityPattern &gt;</name>
    <filename>class_solver_wrapper.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a75a599630086e2edee14eb43ae071733</anchor>
      <arglist>(const parallel::TwoBlockMatrix&lt; dealii::PETScWrappers::MPI::SparseMatrix &gt; &amp;K_stretched, dealii::LinearAlgebra::distributed::Vector&lt; double &gt; &amp;solution, const dealii::PETScWrappers::MPI::BlockVector &amp;f_stretched, const bool symmetric) const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~SolverWrapper</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a738e65e298fb0b13327e853fc9683b6d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::BlockVector&lt; double &gt;, TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt; &gt;, TwoBlockSparsityPattern &gt;</name>
    <filename>class_solver_wrapper.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a75a599630086e2edee14eb43ae071733</anchor>
      <arglist>(const TwoBlockMatrix&lt; dealii::SparseMatrix&lt; double &gt; &gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::BlockVector&lt; double &gt; &amp;f_stretched, const bool symmetric) const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~SolverWrapper</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a738e65e298fb0b13327e853fc9683b6d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::Vector&lt; double &gt;, dealii::SparseMatrix&lt; double &gt;, SparsityPattern &gt;</name>
    <filename>class_solver_wrapper.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a75a599630086e2edee14eb43ae071733</anchor>
      <arglist>(const dealii::SparseMatrix&lt; double &gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::Vector&lt; double &gt; &amp;f_stretched, const bool symmetric) const =0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~SolverWrapper</name>
      <anchorfile>class_solver_wrapper.html</anchorfile>
      <anchor>a738e65e298fb0b13327e853fc9683b6d</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapperPETSc</name>
    <filename>class_solver_wrapper_p_e_t_sc.html</filename>
    <base>SolverWrapper&lt; dealii::LinearAlgebra::distributed::Vector&lt; double &gt;, dealii::PETScWrappers::MPI::BlockVector, parallel::TwoBlockMatrix&lt; dealii::PETScWrappers::MPI::SparseMatrix &gt;, TwoBlockSparsityPattern &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper_p_e_t_sc.html</anchorfile>
      <anchor>ac1033bb655ccc7459c4549e5ae4988c5</anchor>
      <arglist>(const parallel::TwoBlockMatrix&lt; dealii::PETScWrappers::MPI::SparseMatrix &gt; &amp;K_stretched, dealii::LinearAlgebra::distributed::Vector&lt; double &gt; &amp;solution, const dealii::PETScWrappers::MPI::BlockVector &amp;f_stretched, const bool symmetric=false) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>SolverWrapperUMFPACK</name>
    <filename>class_solver_wrapper_u_m_f_p_a_c_k.html</filename>
    <base>SolverWrapper&lt; dealii::Vector&lt; double &gt;, dealii::Vector&lt; double &gt;, dealii::SparseMatrix&lt; double &gt;, SparsityPattern &gt;</base>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>solve</name>
      <anchorfile>class_solver_wrapper_u_m_f_p_a_c_k.html</anchorfile>
      <anchor>a80ed8c6186251e670960945e9ebc0c5b</anchor>
      <arglist>(const dealii::SparseMatrix&lt; double &gt; &amp;K_stretched, dealii::Vector&lt; double &gt; &amp;solution, const dealii::Vector&lt; double &gt; &amp;f_stretched, const bool symmetric=false) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TotalPotential</name>
    <filename>class_total_potential.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type>void</type>
      <name>add_total_potential_contribution</name>
      <anchorfile>class_total_potential.html</anchorfile>
      <anchor>ace8b22eef7789f5c7eb78606ba1f110c</anchor>
      <arglist>(const TotalPotentialContribution&lt; spacedim &gt; &amp;total_potential_contribution)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; SmartPointer&lt; const TotalPotentialContribution&lt; spacedim &gt; &gt; &gt;</type>
      <name>total_potential_contributions</name>
      <anchorfile>class_total_potential.html</anchorfile>
      <anchor>a5a14ce0e2fabf8116566aa67fb11db35</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>max_dependent_vars</name>
      <anchorfile>class_total_potential.html</anchorfile>
      <anchor>a800f9366116679fd0f7d3173a3bfc539</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>AssemblyHelper</name>
      <anchorfile>class_total_potential.html</anchorfile>
      <anchor>af4019c2e39cc934d646aaa35c3c52773</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TotalPotentialContribution</name>
    <filename>class_total_potential_contribution.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>get_potential_contribution</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>a0d281fceeb90ece5c4d2655df5eb9948</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;H_omega_H_sigma_C, const std::vector&lt; Vector&lt; double &gt;&gt; &amp;C_ref_sets, double &amp;Pi, Vector&lt; double &gt; &amp;Pi_1, FullMatrix&lt; double &gt; &amp;Pi_2, const std::tuple&lt; bool, bool, bool &gt; &amp;requested_quantities) const </arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TotalPotentialContribution</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>a1932d6d7c269344542fc8e9221de4448</anchor>
      <arglist>(const std::vector&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; * &gt; &amp;H_omega, const std::vector&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; * &gt; &amp;H_sigma, const std::vector&lt; const IndependentField&lt; 0, spacedim &gt; * &gt; &amp;C)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TotalPotentialContribution</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>ad64c449f4cae79fb8462f4679c50d606</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim, spacedim &gt; &amp;H_omega)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TotalPotentialContribution</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>a4625ad0c798469b9da9bcb62cc0d00f6</anchor>
      <arglist>(const ScalarFunctional&lt; spacedim-1, spacedim &gt; &amp;H_sigma)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TotalPotentialContribution</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>af9ecccfc930f6826be4a903194541d5f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>const bool</type>
      <name>is_primitive</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>a45bfb25a7693c949c26e223cf4a1a1e7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; const ScalarFunctional&lt; spacedim, spacedim &gt; * &gt;</type>
      <name>H_omega</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>a15191539345978a3d0c7293bd7ecaa91</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; const ScalarFunctional&lt; spacedim-1, spacedim &gt; * &gt;</type>
      <name>H_sigma</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>aac404e3a8493d9170541e34bd96673d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>const std::vector&lt; const IndependentField&lt; 0, spacedim &gt; * &gt;</type>
      <name>C</name>
      <anchorfile>class_total_potential_contribution.html</anchorfile>
      <anchor>adea8f8f88243adec43df300e8c8d4593</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TriangulationSystem</name>
    <filename>class_triangulation_system.html</filename>
    <templarg>spacedim</templarg>
    <base>Subscriptor</base>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>InterfaceCell</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a1d62a56e335cf19f4a4ec7932fbaef09</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaIterator&lt; CellAccessor&lt; spacedim, spacedim &gt; &gt;</type>
      <name>DomainCell</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>af53de5ec80a16d9cb167660b2832b240</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaActiveIterator&lt; CellAccessor&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>ActiveInterfaceCell</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a4c3f97884b5478ebd6771433c5273ac5</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>TriaActiveIterator&lt; CellAccessor&lt; spacedim, spacedim &gt; &gt;</type>
      <name>ActiveDomainCell</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a6aa358f6803facdc7e7ee72698a5f4cb</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TriangulationSystem</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a427810bad1fb5457979b8948bb776539</anchor>
      <arglist>(Triangulation&lt; spacedim, spacedim &gt; &amp;tria_domain)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TriangulationSystem</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a9b0dc550fe1566ab1ba6668a868d0a70</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const Triangulation&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_triangulation_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a38345fa4fb7725e66bf0762a1ef0b135</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const Triangulation&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_triangulation_interface</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>aa6615acfade0126dc55d168253c49c7c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Triangulation&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_triangulation_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a4c4e1d0589f962d92532586d808c7da2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual Triangulation&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_triangulation_interface</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a73ef25cdc88395ef6d1ebcf6792cb97f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>close</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>ad3605f3f59fbd55942288026107c4e6d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_interface_cell</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a004ced64c3160c7589b6e7bf12f58019</anchor>
      <arglist>(const DomainCell &amp;domain_cell, const unsigned int face, const types::material_id material_id)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_interface_cells</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a2bb42253dc27cab092c7c76c66b83585</anchor>
      <arglist>(std::vector&lt; std::tuple&lt; const DomainCell, const unsigned int, const types::material_id &gt; &gt; interface_cells)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_interface_manifold</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a54803f193d6f237b2039ba6f5bca09b9</anchor>
      <arglist>(const types::manifold_id manifold_id, const Manifold&lt; spacedim-1, spacedim &gt; &amp;manifold)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_interface_manifolds</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>afb5564d01d3914b15b552f95648d231b</anchor>
      <arglist>(const std::map&lt; types::manifold_id, std::reference_wrapper&lt; const Manifold&lt; spacedim-1, spacedim &gt;&gt;&gt; manifolds)</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_begin_active</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>abbbd02a7e813a0604d072480bda93b1a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_end_active</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>abcab5ca315d720879468cf92a820c29a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt; &amp;</type>
      <name>interface_active_iterators</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a6e3f8531f86b7404e5cb28f5e88e9cd2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_begin_coarse</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>aaf68847ea61f015fce319bb9810833f4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;::iterator</type>
      <name>interface_end_coarse</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a1595839f52f5f5b9d7af31836f9c1f0f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt; &amp;</type>
      <name>interface_coarse_iterators</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>aff36f319e3c837aa9b0522bc63f83bf3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_triangulations_vtk</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a2b90f1509c910930c1404548fa860d8c</anchor>
      <arglist>(const std::string file_name_domain, const std::string file_name_interface) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>refine_global</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a8aa29c77fd798ddd8a5e20806ea21558</anchor>
      <arglist>(const unsigned int times=1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>execute_coarsening_and_refinement</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a62be2563cc8d810a71941e15490f9840</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>check_active_interface_cell_domain_cells_consistency</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a53b005b496d27529bce82f1ed05fed0e</anchor>
      <arglist>(const double tol=1e-12) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::pair&lt; const unsigned int, const unsigned int &gt;</type>
      <name>get_this_proc_n_procs</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>abf6f1a3377be207410b4fd1a9568fb76</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable">
      <type>boost::signals2::signal&lt; void()&gt;</type>
      <name>pre_refinement</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>ad3b16b86e7cbe4800fb5a40114cbc5fb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>boost::signals2::signal&lt; void()&gt;</type>
      <name>post_refinement</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a8813b73cab76c8b77c47c00c39233f76</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>generate_tria_interface_from_tria_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>aee020d825bdf108d69f8e5b71b7ad665</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>const SmartPointer&lt; Triangulation&lt; spacedim, spacedim &gt; &gt;</type>
      <name>tria_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a68bafddc70652cb7c64701c74e86279f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::unique_ptr&lt; Triangulation&lt; spacedim-1, spacedim &gt; &gt;</type>
      <name>tria_interface</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>af145fb3f6a3b3ef9ed92111f9fca5b7d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>bool</type>
      <name>closed</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>aef5ae937cc00b7954357bb8ae88b1f73</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>generate_active_interface_cells_domain_cells</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a4214cf789cf39910305db9fe554249a0</anchor>
      <arglist>(const bool no_assert=false)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>generate_active_interface_cells_domain_cells_recursion</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a449457dd44d9c9d9aaa5754532e8537e</anchor>
      <arglist>(const DomainCell &amp;domain_cell, const unsigned int &amp;face, const InterfaceCell &amp;interface_cell, const bool no_assert)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>check_material_ids_recursion</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a677ea8e08a897b475375eaf4ef70be37</anchor>
      <arglist>(const DomainCell &amp;domain_cell) const </arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual void</type>
      <name>pre_refinement_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>ae1862e6da3157dc8d539fdc0439e9f48</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual void</type>
      <name>post_refinement_domain</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a8435489384095f687363d200ccfce628</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>push_tria_listeners_to_end</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>aca3dc94f1006a603f56ac12c730ae242</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; std::pair&lt; const DomainCell, const unsigned int &gt;, types::material_id &gt;</type>
      <name>coarse_domain_faces_material_ids</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a9a14a8a68bda679330a169413f0e2dfb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::map&lt; dealii::types::manifold_id, const dealii::Manifold&lt; spacedim-1, spacedim &gt; * &gt;</type>
      <name>interface_manifolds</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a130322b74b66f27d8611ce92ad241c6d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;</type>
      <name>coarse_interface_cell_domain_cells</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a5b909c472fe2508da44e89b017aa146c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; InterfaceCellDomainCells&lt; spacedim &gt; &gt;</type>
      <name>active_interface_cell_domain_cells</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a516c7a253cefbc5e208714538b21424d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>std::vector&lt; boost::signals2::connection &gt;</type>
      <name>tria_listeners</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>a713f97bee5570de1f5571d86a41e83f6</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>AssemblyHelper</name>
      <anchorfile>class_triangulation_system.html</anchorfile>
      <anchor>af4019c2e39cc934d646aaa35c3c52773</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>parallel::TriangulationSystem</name>
    <filename>classparallel_1_1_triangulation_system.html</filename>
    <templarg>spacedim</templarg>
    <member kind="function">
      <type></type>
      <name>TriangulationSystem</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a2d0fa99c5da897bedc4d9df93d6cdaf7</anchor>
      <arglist>(dealii::parallel::distributed::Triangulation&lt; spacedim, spacedim &gt; &amp;tria_domain)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TriangulationSystem</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>ad12682d17e85ee169fd53cdbe0536f16</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const dealii::parallel::Triangulation&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_triangulation_domain</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>afcd68b5e999d3c656d728e6313d37957</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const dealii::parallel::Triangulation&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_triangulation_interface</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a6c80e2e3391aa784100984c6aae834e1</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual dealii::parallel::Triangulation&lt; spacedim, spacedim &gt; &amp;</type>
      <name>get_triangulation_domain</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a4901a2f0e26104d5a085253620b98f52</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual dealii::parallel::Triangulation&lt; spacedim-1, spacedim &gt; &amp;</type>
      <name>get_triangulation_interface</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>ab799af10ed705dcf903bcd78325762b0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>close</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a5a7ef73bec1fbe698480e1bd7c6cd8a7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_meshes_per_processor_as_vtu</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a58151f4fc7d91474a24c5666a18aa354</anchor>
      <arglist>(const std::string file_name_domain, const std::string file_name_interface) const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual std::pair&lt; const unsigned int, const unsigned int &gt;</type>
      <name>get_this_proc_n_procs</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a7f9b450a3667344b0a291d8e8dfdee18</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual void</type>
      <name>pre_refinement_domain</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a79e9789e83e12900c85cf8de0644271f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private" virtualness="virtual">
      <type>virtual void</type>
      <name>post_refinement_domain</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a951181f2ad877283d458fa19db42efb2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>update_interface_subdomain_ids</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>a2dfa7c2bf13b929a4d1abf685d5b1448</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="private">
      <type>MPI_Comm</type>
      <name>mpi_communicator</name>
      <anchorfile>classparallel_1_1_triangulation_system.html</anchorfile>
      <anchor>aff7cdcf04d5a4fb633d714130da893b0</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TwoBlockMatrix</name>
    <filename>class_two_block_matrix.html</filename>
    <templarg>MatrixType</templarg>
    <member kind="typedef">
      <type>double</type>
      <name>value_type</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a9f50fb1f98df3ae0048c3c5220aa265f</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TwoBlockMatrix</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>aa74bb5722415e61241b5f38c9df29d39</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TwoBlockMatrix</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a861cdc13c0fc6a65ec7d99e4d2349a92</anchor>
      <arglist>(const TwoBlockSparsityPattern &amp;sp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TwoBlockMatrix</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a0674b188ebc7bde3ee8f5eb9a94d6265</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a5eaa76e9909c600aaa256aa77b577391</anchor>
      <arglist>(const TwoBlockSparsityPattern &amp;sp)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>m</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ab15732600f7e6917f961dce5e24b61a6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a7675abf4c2532f48eb3b6bf153755a26</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a9b5d922689d1f60cbc59d2e424df2dc4</anchor>
      <arglist>(const unsigned int i, const unsigned int j, const double value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a7373229281e1c3712681c1ac742dffb7</anchor>
      <arglist>(const unsigned int row, const unsigned int n_cols, const unsigned int *col_indices, const double *values, const bool elide_zero_values=false, const bool col_indices_are_sorted=false)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a5457fd3b6228bed7a731bab2569a17a0</anchor>
      <arglist>(const unsigned int i, const unsigned int j, const double value)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a618495b0ad4cef72d95db7a1f380c20f</anchor>
      <arglist>(const unsigned int i, const unsigned int j) const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>el</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a67a13d33b3110985bfc604969fc454c9</anchor>
      <arglist>(const unsigned int i, const unsigned int j) const </arglist>
    </member>
    <member kind="function">
      <type>TwoBlockMatrix&lt; MatrixType &gt; &amp;</type>
      <name>operator=</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ac008679a935fb6577931058ac0a6be36</anchor>
      <arglist>(const double value)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compress</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a8de26ebbb410811059d162fb164cedce</anchor>
      <arglist>(const VectorOperation::values operation)</arglist>
    </member>
    <member kind="function">
      <type>const MatrixType &amp;</type>
      <name>get_A</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a304eac57c72f23cbcf18a04e5bc64508</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const MatrixType &amp;</type>
      <name>get_B</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ab3807e975c3a5de9779a80c8d61fad96</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const MatrixType &amp;</type>
      <name>get_C</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ac0b0e0bc16412dfd73a8c6a3f782b8d9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const MatrixType &amp;</type>
      <name>get_D</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>af905a0e4b0504e45937a06df768a2ae2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>MatrixType &amp;</type>
      <name>get_A</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ab30f0e3b8f803456c3a2f1ce4f753081</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MatrixType &amp;</type>
      <name>get_B</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ae41ca5bbc0495c8ffb1604da488a33b7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MatrixType &amp;</type>
      <name>get_C</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>aeaa4102bb137700d7bbe0262ee16e8ee</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>MatrixType &amp;</type>
      <name>get_D</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a3725800849e26e2192c0ac6520777dab</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_block_0_size</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a688b1a69284ec154c493a99ec7cad704</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>get_block_1_size</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a9ecb16be5aa5ebbb916d256a9ded7c3e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>MatrixType</type>
      <name>A</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a91b02f00cf56b398e1624227e6efc8c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>MatrixType</type>
      <name>B</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ae8edfe1cd2298741586dcc46a4d09ebf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>MatrixType</type>
      <name>C</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a847505a77448124928540be0367d8747</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>MatrixType</type>
      <name>D</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a72c598f36d3e8cce12821a09920ded2d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SparsityPattern</type>
      <name>sp_A</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>af40216911c197c73687f7d3a1ff4e68e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SparsityPattern</type>
      <name>sp_B</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a854c7e77623a7c0012aa0ceeb60c147a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SparsityPattern</type>
      <name>sp_C</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a45d6e534ad775577c86e9ac9a116014e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>SparsityPattern</type>
      <name>sp_D</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>adb1833352096be79279134b218659665</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned int</type>
      <name>block_0_size</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>ae862db333ed43f1650d4926f237cfc43</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned int</type>
      <name>block_1_size</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a60ddea1d26313aac35cfb2ae9f192cb2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned int</type>
      <name>total_dimension</name>
      <anchorfile>class_two_block_matrix.html</anchorfile>
      <anchor>a7d63b4c4c1e499c9f9913444dc5fd6b6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>parallel::TwoBlockMatrix</name>
    <filename>classparallel_1_1_two_block_matrix.html</filename>
    <templarg>MatrixType</templarg>
    <member kind="function">
      <type></type>
      <name>TwoBlockMatrix</name>
      <anchorfile>classparallel_1_1_two_block_matrix.html</anchorfile>
      <anchor>a93990616301e671f9a967fbc49029674</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TwoBlockMatrix</name>
      <anchorfile>classparallel_1_1_two_block_matrix.html</anchorfile>
      <anchor>afa9b3f1c3d74a211c334c8f6e59b52f6</anchor>
      <arglist>(const TwoBlockSparsityPattern &amp;sp, const IndexSet &amp;locally_owned_indices, const MPI_Comm mpi_communicator=MPI_COMM_WORLD)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>classparallel_1_1_two_block_matrix.html</anchorfile>
      <anchor>a0e112695034c69aa8a430efd9e7c37ef</anchor>
      <arglist>(const TwoBlockSparsityPattern &amp;sp, const IndexSet &amp;locally_owned_indices, const MPI_Comm mpi_communicator=MPI_COMM_WORLD)</arglist>
    </member>
    <member kind="function">
      <type>dealii::GalerkinTools::parallel::TwoBlockMatrix&lt; MatrixType &gt; &amp;</type>
      <name>operator=</name>
      <anchorfile>classparallel_1_1_two_block_matrix.html</anchorfile>
      <anchor>a8e2ba56ebc9ee5ad1fe7b1d8604e9669</anchor>
      <arglist>(const double value)</arglist>
    </member>
    <member kind="function">
      <type>const MPI_Comm &amp;</type>
      <name>get_communicator</name>
      <anchorfile>classparallel_1_1_two_block_matrix.html</anchorfile>
      <anchor>a4894068aad11986f43c1221cd696dc2a</anchor>
      <arglist>() const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>TwoBlockSparsityPattern</name>
    <filename>class_two_block_sparsity_pattern.html</filename>
    <member kind="function">
      <type></type>
      <name>TwoBlockSparsityPattern</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a5bcc54c1c6b0d769807644f1c87fa136</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>TwoBlockSparsityPattern</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ab7240dc72bf5fdb32fc192f70e6c68b0</anchor>
      <arglist>(const IndexSet &amp;locally_relevant_indices, const unsigned int block_0_size)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~TwoBlockSparsityPattern</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>abe3cbf4d3ee10b66234c5f8b566a938f</anchor>
      <arglist>()=default</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ae18398b187a0dba1570abfec1a4e082d</anchor>
      <arglist>(const IndexSet &amp;locally_relevant_indices, const unsigned int block_0_size)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reinit</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a8c9a56a20123e8eaea65f9eb6efed250</anchor>
      <arglist>(const AssemblyHelper&lt; spacedim &gt; &amp;assembly_helper)</arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_rows</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a1f95d380bc078a6abf6ef67fcb076e0d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned int</type>
      <name>n_cols</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a33a0e1a01a788c73b474426128b18179</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ab3787657efc568999fbf24cb5208d46b</anchor>
      <arglist>(const unsigned int i, const unsigned int j)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_entries</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a97c69542e0e22e5653c8b4f7b4275635</anchor>
      <arglist>(const unsigned int row, ForwardIterator begin, ForwardIterator end, const bool indices_are_sorted=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>finalize</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a051b557d53ca8367228b7829c55d4af5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const SparsityPattern &amp;</type>
      <name>get_sp_A</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a595b3576cc8e02e24a795c9cd639af01</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const SparsityPattern &amp;</type>
      <name>get_sp_B</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a36c1dbb12ab6ac24fc30e793e59d2c3b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const SparsityPattern &amp;</type>
      <name>get_sp_C</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a156a0a521faae238543e76aa05a8bad7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>const SparsityPattern &amp;</type>
      <name>get_sp_D</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a0990cb1eabe18204dd46ac30e1cc0dda</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>exists</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a32bae81b9aec47fa34ddd40c963fb7e8</anchor>
      <arglist>(const unsigned int i, const unsigned int j) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>distribute</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a61ef6e3f3325cffcfc1d25efc41886b0</anchor>
      <arglist>(const IndexSet &amp;locally_owned_indices, const MPI_Comm mpi_communicator)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DynamicSparsityPattern</type>
      <name>dsp_A</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ac3bafc713fec0ec2232fc084f880c50c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DynamicSparsityPattern</type>
      <name>dsp_B</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ab1b2d6dc5d92f6353b80dd7a8d32ffa4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DynamicSparsityPattern</type>
      <name>dsp_C</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a444a6465b3ba1355fb675128c6e9800a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>DynamicSparsityPattern</type>
      <name>dsp_D</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a853d4000754f009a66853954d8c540bf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>SparsityPattern</type>
      <name>sp_A</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ad09712fc57a8cf5a2cd2cb9fc238a7cc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>SparsityPattern</type>
      <name>sp_B</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a95e817472f023d7237c752a18a24ac73</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>SparsityPattern</type>
      <name>sp_C</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>aee8da566a426045fb8101dba91295a27</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>SparsityPattern</type>
      <name>sp_D</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a7a820cba98001e45ba09551a306b14d7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>block_0_size</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a81789900ac1a9f5bfd7f85ad730eb03b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>block_1_size</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>ab215e3553e07b5f5872b1e01228968f2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned int</type>
      <name>total_dimension</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>aa560457e8b8bf7a06acf6b19987c5ee0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>finalized</name>
      <anchorfile>class_two_block_sparsity_pattern.html</anchorfile>
      <anchor>a7c930cd0dd5e1d6fb1550d6c3939c7aa</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>Auxiliary</name>
    <filename>namespace_auxiliary.html</filename>
    <member kind="function">
      <type>int</type>
      <name>compute_ldr</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>abe796c1529d11eed08fd05bca82f3002</anchor>
      <arglist>(FullMatrix&lt; double &gt; &amp;C, Vector&lt; double &gt; &amp;D, std::vector&lt; Vector&lt; double &gt;&gt; &amp;L, std::vector&lt; Vector&lt; double &gt;&gt; &amp;R)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>convert_local_indices_to_global_indices</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>ad3c42d209f0ba8e9d4ce305060634bf1</anchor>
      <arglist>(const std::vector&lt; unsigned int &gt; &amp;dof_indices_local, std::vector&lt; unsigned int &gt; &amp;dof_indices_global, const std::vector&lt; unsigned int &gt; &amp;dof_indices_local_global)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>combine_dof_indices</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a1d90ebc8738df3d8c70b540034137019</anchor>
      <arglist>(const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_interface, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_minus, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_plus, const std::vector&lt; unsigned int &gt; &amp;dof_indices_global_C, std::vector&lt; unsigned int &gt; &amp;dof_indices_interface_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_minus_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_plus_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_C_dof_indices_combined, std::vector&lt; unsigned int &gt; &amp;dof_indices_global_combined)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_to_index_set</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>aca5fe5966aef4ba2293fb44095bfd86c</anchor>
      <arglist>(const DoFRenumbering &amp;dof_renumbering, const IndexSet &amp;in, IndexSet &amp;out)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>renumber_constraints</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>aa6148bcbaf5e3003717b4dd2a4da15b3</anchor>
      <arglist>(AffineConstraints&lt; double &gt; &amp;constraint_matrix, const DoFRenumbering &amp;dof_renumbering=DoFRenumbering(), const bool close=true)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_dof_renumbering_contiguous</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a4261fc726ff965166f3f53242918acea</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dof_handler_system, DoFRenumberingOffset &amp;dof_renumbering_offset)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>compute_map_dofs</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a072f85e6d745ae3c532bb0612f4bd3ce</anchor>
      <arglist>(const DoFHandlerSystem&lt; spacedim &gt; &amp;dhs_1, const DoFHandlerSystem&lt; spacedim &gt; &amp;dhs_2, std::vector&lt; unsigned int &gt; &amp;map_dofs)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>split_vector</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>a59c01a6511bffd7442693e86cd194ef1</anchor>
      <arglist>(const Vector&lt; double &gt; &amp;in, Vector&lt; double &gt; &amp;out_0, Vector&lt; double &gt; &amp;out_1, const unsigned int size_1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>split_matrix</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>af746a1d08b1135a3684fc990f7b1384d</anchor>
      <arglist>(const FullMatrix&lt; double &gt; &amp;in, FullMatrix&lt; double &gt; &amp;out_00, FullMatrix&lt; double &gt; &amp;out_01, FullMatrix&lt; double &gt; &amp;out_10, FullMatrix&lt; double &gt; &amp;out_11, const unsigned int size_1)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>communicate_bool</name>
      <anchorfile>namespace_auxiliary.html</anchorfile>
      <anchor>ad0ff28386be7b54b3487ae36e5a074fa</anchor>
      <arglist>(const bool local_bool, const MPI_Comm &amp;mpi_communicator)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>parallel</name>
    <filename>namespaceparallel.html</filename>
    <class kind="class">parallel::TriangulationSystem</class>
    <class kind="class">parallel::TwoBlockMatrix</class>
    <member kind="function">
      <type>void</type>
      <name>transform</name>
      <anchorfile>namespaceparallel.html</anchorfile>
      <anchor>a6fc7c301550b8e89fddedc5f82d3fbe6</anchor>
      <arglist>(const InputIterator &amp;begin_in, const InputIterator &amp;end_in, OutputIterator out, Predicate &amp;predicate, const unsigned int grainsize)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>transform</name>
      <anchorfile>namespaceparallel.html</anchorfile>
      <anchor>a06ef7174350a0e70b082b7b6ad062dd6</anchor>
      <arglist>(const InputIterator1 &amp;begin_in1, const InputIterator1 &amp;end_in1, InputIterator2 in2, OutputIterator out, Predicate &amp;predicate, const unsigned int grainsize)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>transform</name>
      <anchorfile>namespaceparallel.html</anchorfile>
      <anchor>a9d2e003eb2c53199fbd7ed1721d97f3d</anchor>
      <arglist>(const InputIterator1 &amp;begin_in1, const InputIterator1 &amp;end_in1, InputIterator2 in2, InputIterator3 in3, OutputIterator out, Predicate &amp;predicate, const unsigned int grainsize)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>apply_to_subranges</name>
      <anchorfile>namespaceparallel.html</anchorfile>
      <anchor>aa860d510323b29ede22a7f69f5dc41ad</anchor>
      <arglist>(const RangeType &amp;begin, const typename identity&lt; RangeType &gt;::type &amp;end, const Function &amp;f, const unsigned int grainsize)</arglist>
    </member>
    <member kind="function">
      <type>ResultType</type>
      <name>accumulate_from_subranges</name>
      <anchorfile>namespaceparallel.html</anchorfile>
      <anchor>a7ffe536c55823d2dd1f60d12d2ca6afe</anchor>
      <arglist>(const Function &amp;f, const RangeType &amp;begin, const typename identity&lt; RangeType &gt;::type &amp;end, const unsigned int grainsize)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>GalerkinTools library</title>
    <filename>index</filename>
  </compound>
</tagfile>
