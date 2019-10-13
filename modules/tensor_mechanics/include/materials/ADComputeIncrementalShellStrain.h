//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADMaterial.h"
#include "libmesh/dense_matrix.h"

#define usingComputeIncrementalShellStrainMembers usingMaterialMembers;

// Forward Declarations
template <ComputeStage>
class ADComputeIncrementalShellStrain;

namespace libMesh
{
  class QGauss;
}

template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;

declareADValidParams(ADComputeIncrementalShellStrain);

template <ComputeStage compute_stage>
class ADComputeIncrementalShellStrain : public ADMaterial<compute_stage>
{
public:
  ADComputeIncrementalShellStrain(const InputParameters & parameters);
protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeProperties() override;
  virtual void computeLargeStrain();
  virtual void initialSetup() override;

  /// Number of coupled rotational variables
  unsigned int _nrot;

  /// Number of coupled displacement variables
  unsigned int _ndisp;

  /// Variable numbers corresponding to the rotational variables
  std::vector<unsigned int> _rot_num;

  /// Variable numbers corresponding to the displacement variables
  std::vector<unsigned int> _disp_num;

  /// Coupled variable for the shell thickness
  const VariableValue & _thickness;

  /// Flag to compute large strains
  const bool _large_strain;

  /// Strain increment in the covariant coordinate system
  std::vector<ADMaterialProperty(RankTwoTensor) *> _strain_increment;

  /// Total strain increment in the covariant coordinate system
  std::vector<ADMaterialProperty(RankTwoTensor) *> _total_strain;

  /// Old total strain increment in the covariant coordinate system
  std::vector<const MaterialProperty<RankTwoTensor> *> _total_strain_old;

  /// Reference to the nonlinear system object
  NonlinearSystemBase & _nonlinear_sys;

  /// Indices of solution vector corresponding to displacement DOFs in 3 directions at the 4 nodes
  std::vector<std::vector<unsigned int>> _soln_disp_index;

  /// Indices of solution vector corresponding to rotation DOFs in 2 directions at the 4 nodes
  std::vector<std::vector<unsigned int>> _soln_rot_index;

  /// Vector that stores the incremental solution at all the 20 DOFs in the 4 noded element.
  ADDenseVector _soln_vector;

  ADDenseVector _strain_vector;

  std::vector<const Node *> _nodes;

  /// Material property storing the normal to the element at the 4 nodes. Stored as a material property for convinience.
  ADMaterialProperty(RealVectorValue) & _node_normal;

  /// Material property storing the old normal to the element at the 4 nodes.
  const MaterialProperty<RealVectorValue> & _node_normal_old;

  /// Quadrature rule in the out of plane direction
  std::unique_ptr<QGauss> _t_qrule;

  /// Quadrature points in the out of plane direction in isoparametric coordinate system
  std::vector<Point> _t_points;

  /// Quadrature points in the in-plane direction in isoparametric coordinate system
  std::vector<Point> _2d_points;

  /// Derivatives of shape functions w.r.t isoparametric coordinates xi
  std::vector<std::vector<Real>> _dphidxi_map;

  /// Derivatives of shape functions w.r.t isoparametric coordinates eta
  std::vector<std::vector<Real>> _dphideta_map;

  /// Shape function value
  std::vector<std::vector<Real>> _phi_map;

  /// Derivative of global x, y and z w.r.t isoparametric coordinate xi
  std::vector<std::vector<ADRealVectorValue>> _dxyz_dxi;

  /// Derivative of global x, y and z w.r.t isoparametric coordinate eta
  std::vector<std::vector<ADRealVectorValue>> _dxyz_deta;

  /// Derivative of global x, y and z w.r.t isoparametric coordinate zeta
  std::vector<std::vector<ADRealVectorValue>> _dxyz_dzeta;

  /// First tangential vectors at nodes
  std::vector<ADRealVectorValue> _V1;

  /// First tangential vectors at nodes
  std::vector<ADRealVectorValue> _V2;

  /// B_matrix for small strain
  std::vector<ADMaterialProperty(DenseMatrix<Real>) *> _B;

  std::vector<ADMaterialProperty(DenseMatrix<Real>) *> _BNL;

  /// ge matrix for elasticity tensor conversion
  std::vector<ADMaterialProperty(RankTwoTensor) *> _ge;

  /// Material property containing jacobian of transformation
  std::vector<ADMaterialProperty(Real) *> _Jmap;

  usingMaterialMembers;
};
