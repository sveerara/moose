/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTESHELLSTRAIN_H
#define COMPUTESHELLSTRAIN_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "ColumnMajorMatrix.h"
#include "libmesh/quadrature.h"

/**
 * ComputeShellStrain defines a displacement and rotation strain increment and rotation
 * increment (=1), for small strains.
 */

// Forward Declarations
class ComputeShellStrain;

namespace libMesh
{
  class QGauss;
}

template <>
InputParameters validParams<ComputeShellStrain>();

class ComputeShellStrain : public Material
{
public:
  ComputeShellStrain(const InputParameters & parameters);

protected:
  virtual void computeProperties() override;
  virtual void initQpStatefulProperties() override;
  virtual void computeLargeStrain();

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
  MaterialProperty<std::vector<RankTwoTensor>> & _strain_increment;

  /// Total strain increment in the covariant coordinate system
  MaterialProperty<std::vector<RankTwoTensor>> & _total_strain;

  /// Old total strain increment in the covariant coordinate system
  const MaterialProperty<std::vector<RankTwoTensor>> & _total_strain_old;

  /// Reference to the nonlinear system object
  NonlinearSystemBase & _nonlinear_sys;

  /// Indices of solution vector corresponding to displacement DOFs in 3 directions at the 4 nodes
  std::vector<std::vector<unsigned int>> _soln_disp_index;

  /// Indices of solution vector corresponding to rotation DOFs in 2 directions at the 4 nodes
  std::vector<std::vector<unsigned int>> _soln_rot_index;

  /// Vector that stores the incremental solution at all the 20 DOFs in the 4 noded element.
  ColumnMajorMatrix _soln_vector;

  ColumnMajorMatrix _strain_vector;

  std::vector<const Node *> _nodes;

  /// Material property storing the normal to the element at the 4 nodes. Stored as a material property for convinience.
  MaterialProperty<RealVectorValue> & _node_normal;

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
  std::vector<std::vector<RealVectorValue>> _dxyz_dxi;

  /// Derivative of global x, y and z w.r.t isoparametric coordinate eta
  std::vector<std::vector<RealVectorValue>> _dxyz_deta;

  /// Derivative of global x, y and z w.r.t isoparametric coordinate zeta
  std::vector<std::vector<RealVectorValue>> _dxyz_dzeta;

  /// First tangential vectors at nodes
  std::vector<RealVectorValue> _V1;

  /// First tangential vectors at nodes
  std::vector<RealVectorValue> _V2;

  /// B_matrix for small strain
  MaterialProperty<std::vector<ColumnMajorMatrix>> & _B;
  MaterialProperty<std::vector<ColumnMajorMatrix>> & _BNL_new;

  MaterialProperty<std::vector<ColumnMajorMatrix>> * _BNL;
  const MaterialProperty<std::vector<ColumnMajorMatrix>> * _BNL_old;


  /// ge matrix for elasticity tensor conversion
  MaterialProperty<std::vector<RankTwoTensor>> & _ge;

  /// Material property containing jacobian of transformation
  MaterialProperty<std::vector<Real>> & _Jmap;

  MaterialProperty<std::vector<Real>> & _soln_vector_prop;

};

#endif // COMPUTESHELLSTRAIN_H
