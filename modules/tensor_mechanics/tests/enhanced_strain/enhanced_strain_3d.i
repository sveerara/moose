#This is a test of the piece-wise linear strain hardening model using the small strain formulation.
#The exact same problem was run in Abaqus with exactly the same result.

[Mesh]
  file = 3d_cuboid_square_10.e
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
[]


[AuxVariables]
  [./stress_yy]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]


[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
  [../]
[]

[AuxKernels]
  [./stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    index_i = 1
    index_j = 1
    variable = stress_yy
    execute_on = timestep_end
  [../]
 []

[NodalKernels]
  [./force1]
    type = ConstantRate
    rate = -666.6666666666666667
    boundary = 2
    variable = disp_x
  [../]
  [./force2]
    type = ConstantRate
    rate = 666.6666666666666667
    boundary = 3
    variable = disp_x
  [../]
[]

[BCs]
  [./x_left]
    type = DirichletBC
    variable = disp_x
    boundary = 1
    value = 0.0
  [../]

  [./y_left]
    type = DirichletBC
    variable = disp_y
    boundary = 1
    value = 0.0
  [../]

  [./z_left]
    type = DirichletBC
    variable = disp_z
    boundary = 1
    value = 0.0
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 1500.0
    poissons_ratio = 0.25
    block = 1
  [../]
  [./strain]
    type = ComputeSmallStrain
    displacements = 'disp_x disp_y disp_z'
    enhanced_strain = true
    block = 1
  [../]
  [./stress]
    type = ComputeLinearElasticStress
    block = 1
  [../]
[]

[Executioner]
  type = Steady

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'
[]

[Postprocessors]
  [./disp_1]
    type = PointValue
    point = '5.0 -1.0 -1.0'
    variable = disp_y
  [../]
  [./disp_2]
    type = PointValue
    point = '5.0 1.0 1.0'
    variable = disp_y
  [../]
[]

[Outputs]
  [./out]
    type = Exodus
    elemental_as_nodal = true
  [../]
[]
