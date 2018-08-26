# Test for small strain euler beam vibration in y direction

# An impulse load is applied at the end of a cantilever beam of length 4m.
# The beam is massless with a lumped mass at the end of the beam
# The properties of the cantilever beam are as follows:
# Young's modulus (E) = 1e4
# Shear modulus (G) = 4e7
# Shear coefficient (k) = 1.0
# Cross-section area (A) = 0.01
# Iy = 1e-4 = Iz
# Length (L)= 4 m
# mass (m) = 0.01899772

# For this beam, the dimensionless parameter alpha = kAGL^2/EI = 6.4e6
# Therefore, the beam behaves like a Euler-Bernoulli beam.

# The theoretical first frequency of this beam is:
# f1 = 1/(2 pi) * sqrt(3EI/(mL^3)) = 0.25

# This implies that the corresponding time period of this beam is 4s.

# The FEM solution for this beam with 10 element gives time periods of 4s with time step of 0.01s.
# A higher time step of 0.1 s is used in the test to reduce computational time.

# The time history from this analysis matches with that obtained from Abaqus.

# Values from the first few time steps are as follows:
# time   disp_y                vel_y                accel_y
# 0.0    0.0                   0.0                  0.0
# 0.1    0.0013076435060869    0.026152870121738    0.52305740243477
# 0.2    0.0051984378734383    0.051663017225289   -0.01285446036375
# 0.3    0.010269120909367     0.049750643493289   -0.02539301427625
# 0.4    0.015087433925158     0.046615616822532   -0.037307519138892
# 0.5    0.019534963888307     0.042334982440433   -0.048305168503101

[Mesh]
  type = GeneratedMesh
  xmin = 0.0
  xmax = 4.0
  nx = 10
  dim = 1
  displacements = 'disp_x disp_y disp_z'
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
  [./rot_z]
    order = FIRST
    family = LAGRANGE
  [../]
[]

[BCs]
  [./fixx1]
    type = DirichletBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./fixy1]
    type = DirichletBC
    variable = disp_y
    boundary = left
    value = 0.0
  [../]
  [./fixz1]
    type = DirichletBC
    variable = disp_z
    boundary = left
    value = 0.0
  [../]
  [./fixr1]
    type = DirichletBC
    variable = rot_x
    boundary = left
    value = 0.0
  [../]
  [./fixr2]
    type = DirichletBC
    variable = rot_y
    boundary = left
    value = 0.0
  [../]
  [./fixr3]
    type = DirichletBC
    variable = rot_z
    boundary = left
    value = 0.0
  [../]
  [./disp_x_eigen]
    type = EigenDirichletBC
    variable = disp_x
    boundary = 'left'
    use_displaced_mesh = true
  [../]
  [./disp_y_eigen]
    type = EigenDirichletBC
    variable = disp_y
    boundary = 'left'
    use_displaced_mesh = true
  [../]
  [./disp_z_eigen]
    type = EigenDirichletBC
    variable = disp_z
    boundary = 'left'
    use_displaced_mesh = true
  [../]
  [./rot_x_eigen]
    type = EigenDirichletBC
    variable = rot_x
    boundary = 'left'
    use_displaced_mesh = true
  [../]
  [./rot_y_eigen]
    type = EigenDirichletBC
    variable = rot_y
    boundary = 'left'
    use_displaced_mesh = true
  [../]
  [./rot_z_eigen]
    type = EigenDirichletBC
    variable = rot_z
    boundary = 'left'
    use_displaced_mesh = true
  [../]
[]

#[NodalKernels]
#  [./x_inertial]
#    type = NodalMassMatrix
#    variable = disp_x
#    #boundary = 'left right'
#    block = 0
#    mass = 0.01899772
#    extra_matrix_tags = 'eigen'
#    extra_vector_tags = 'eigen'
#  [../]
#  [./y_inertial]
#    type = NodalMassMatrix
#    variable = disp_y
#    #boundary = 'left right'
#    block = 0
#    mass = 0.01899772
#    extra_matrix_tags = 'eigen'
#    extra_vector_tags = 'eigen'
#  [../]
#  [./z_inertial]
#    type = NodalMassMatrix
#    variable = disp_z
#    #boundary = 'left right'
#    block = 0
#    mass = 0.01899772
#    extra_matrix_tags = 'eigen'
#    extra_vector_tags = 'eigen'
#  [../]
#  [./x_inertial_rot]
#    type = NodalMassMatrix
#    variable = rot_x
#    #boundary = 'left right'
#    block = 0
#    mass = 0.01899772
#    extra_matrix_tags = 'eigen'
#    extra_vector_tags = 'eigen'
#  [../]
#  [./y_inertial_rot]
#    type = NodalMassMatrix
#    variable = rot_y
#    #boundary = 'left right'
#    block = 0
#    mass = 0.01899772
#    extra_matrix_tags = 'eigen'
#    extra_vector_tags = 'eigen'
#  [../]
#  [./z_inertial_rot]
#    type = NodalMassMatrix
#    variable = rot_z
#    #boundary = 'left right'
#    block = 0
#    mass = 0.01899772
#    extra_matrix_tags = 'eigen'
#    extra_vector_tags = 'eigen'
#  [../]
#[]

[Functions]
  [./force]
    type = PiecewiseLinear
    x = '0.0 0.1 0.2 10.0'
    y = '0.0 1e-2  0.0  0.0'
  [../]
[]

[Preconditioning]
  [./smp]
    type = SMP
    full = true
  [../]
[]

[Executioner]
  type = Eigenvalue
  eigen_problem_type = gen_non_hermitian
  which_eigen_pairs = SMALLEST_MAGNITUDE
  n_eigen_pairs = 5
  n_basis_vectors = 18
  #solve_type = MF_Nonlinear_power
  #solve_type = krylovschur
  solve_type = jacobi_davidson
  #solve_type = MONOLITH_NEWTON
  petsc_options = '-eps_view'
[]

[Kernels]
  [./solid_disp_x]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 0
    variable = disp_x
  [../]
  [./solid_disp_y]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 1
    variable = disp_y
  [../]
  [./solid_disp_z]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 2
    variable = disp_z
  [../]
  [./solid_rot_x]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 3
    variable = rot_x
  [../]
  [./solid_rot_y]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 4
    variable = rot_y
  [../]
  [./solid_rot_z]
    type = StressDivergenceBeam
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    component = 5
    variable = rot_z
  [../]
  [./inertia_disp_x]
    type = NodalLumpedMassKernel
    block = '0'
    variable = disp_x
    extra_matrix_tags = 'eigen'
    extra_vector_tags = 'eigen'
    nodal_mass_file = 'nodal_mass.csv'
  [../]
  [./inertia_disp_y]
    type = NodalLumpedMassKernel
    block = '0'
    variable = disp_y
    extra_matrix_tags = 'eigen'
    extra_vector_tags = 'eigen'
    nodal_mass_file = 'nodal_mass.csv'
  [../]
  [./inertia_disp_z]
    type = NodalLumpedMassKernel
    block = '0'
    variable = disp_z
    extra_matrix_tags = 'eigen'
    extra_vector_tags = 'eigen'
    nodal_mass_file = 'nodal_mass.csv'
  [../]
[]

[Materials]
  [./elasticity]
    type = ComputeElasticityBeam
    youngs_modulus = 1.0e4
    poissons_ratio = -0.999875
    shear_coefficient = 1.0
    block = 0
  [../]
  [./strain]
    type = ComputeIncrementalBeamStrain
    block = '0'
    displacements = 'disp_x disp_y disp_z'
    rotations = 'rot_x rot_y rot_z'
    area = 0.01
    Ay = 0.0
    Az = 0.0
    Iy = 1.0e-4
    Iz = 1.0e-4
    y_orientation = '0.0 1.0 0.0'
  [../]
  [./stress]
    type = ComputeBeamResultants
    block = 0
  [../]
[]

[VectorPostprocessors]
  [./eigenvalues]
    type = Eigenvalues
    execute_on = 'timestep_end'
  [../]
[]

[Outputs]
  csv = true
  execute_on = 'timestep_end'
  [./console]
    type = Console
    outlier_variable_norms = false
  [../]
[]
