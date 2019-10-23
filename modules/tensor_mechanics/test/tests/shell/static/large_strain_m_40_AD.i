[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 1
  ny = 5
  xmin = 0.0
  xmax = 1.0
  ymin = 0.0
  ymax = 40.0
[]

[MeshModifiers]
  [./cnode]
    type = AddExtraNodeset
    coord = '0.0 0.0'
    new_boundary = 100
  [../]
[]

[Variables]
  [./disp_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_y]
    order = FIRST
    family = LAGRANGE
  [../]
  [./disp_z]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_x]
    order = FIRST
    family = LAGRANGE
  [../]
  [./rot_y]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[BCs]
  [./fixy1]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./fixz1]
    type = PresetBC
    variable = disp_z
    boundary = bottom
    value = 0.0
  [../]
  [./fixr1]
    type = PresetBC
    variable = rot_x
    boundary = bottom
    value = 0.0
  [../]
  [./fixr2]
    type = PresetBC
    variable = rot_y
    boundary = bottom
    value = 0.0
  [../]
  [./fixx1]
    type = PresetBC
    variable = disp_x
    boundary = bottom
    value = 0.0
  [../]
[]

[NodalKernels]
  [./force_y2]
    type = UserForcingFunctionNodalKernel
    variable = disp_z
    boundary = top
    function = force_y
  [../]
[]

[Functions]
  [./force_y]
    type = PiecewiseLinear
    x = '0.0 1.0 3.0'
    y = '0.0 1.0 1.0'
    scale_factor =0.147262156
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  line_search = 'none'
#nl_max_its = 2
  nl_rel_tol = 1e-10
  nl_abs_tol = 1e-8

  dt = 0.1
  dtmin = 0.1
  end_time = 3.0
[]

[Kernels]
  [./solid_disp_x]
    type = ADStressDivergenceShell
    block = '0'
    component = 0
    variable = disp_x
    order = SECOND
    large_strain = true
  [../]
  [./solid_disp_y]
    type = ADStressDivergenceShell
    block = '0'
    component = 1
    variable = disp_y
    order = SECOND
    large_strain = true
  [../]
  [./solid_disp_z]
    type = ADStressDivergenceShell
    block = '0'
    component = 2
    variable = disp_z
    order = SECOND
    large_strain = true
  [../]
  [./solid_rot_x]
    type = ADStressDivergenceShell
    block = '0'
    component = 3
    variable = rot_x
    order = SECOND
    large_strain = true
  [../]
  [./solid_rot_y]
    type = ADStressDivergenceShell
    block = '0'
    component = 4
    variable = rot_y
    order = SECOND
    large_strain = true
  [../]
[]

[Materials]
  [./elasticity]
    type = ADComputeIsotropicElasticityTensorShell
    youngs_modulus = 1800
    poissons_ratio = 0.0
    block = 0
    order = SECOND
  [../]
  [./strain]
    type = ADComputeFiniteShellStrain
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y'
    thickness = 1.0
    order = SECOND
  [../]
  [./stress]
    type = ADComputeShellStress
    block = 0
    order = SECOND
  [../]
[]

[Postprocessors]
  [./disp_z2]
    type = PointValue
    point = '1.0 40.0 0.0'
    variable = disp_z
  [../]
  [./disp_y2]
    type = PointValue
    point = '1.0 40.0 0.0'
    variable = disp_y
  [../]
  [./rot_x2]
    type = PointValue
    point = '1.0 40.0 0.0'
    variable = rot_x
  [../]
[]

[Outputs]
  exodus = true
[]
