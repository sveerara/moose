/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "StressDivergenceShell.h"

// MOOSE includes
#include "Assembly.h"
#include "Material.h"
#include "MooseVariable.h"
#include "SystemBase.h"
#include "RankTwoTensor.h"
#include "NonlinearSystem.h"
#include "MooseMesh.h"
#include "ArbitraryQuadrature.h"
#include "ColumnMajorMatrix.h"

#include "libmesh/quadrature.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/string_to_enum.h"

registerMooseObject("TensorMechanicsApp", StressDivergenceShell);

template <>
InputParameters
validParams<StressDivergenceShell>()
{
  InputParameters params = validParams<Kernel>();
  params.addClassDescription("Quasi-static stress divergence kernel for Shell element");
  params.addRequiredParam<unsigned int>(
      "component",
      "An integer corresponding to the direction "
      "the variable this kernel acts in. (0 for disp_x, "
      "1 for disp_y, 2 for disp_z, 3 for rot_x, 4 for rot_y)");
  params.addRequiredParam<std::string>("order", "Quadrature order in out of plane direction");
  params.set<bool>("use_displaced_mesh") = true;
  return params;
}

StressDivergenceShell::StressDivergenceShell(const InputParameters & parameters)
  : Kernel(parameters),
    _component(getParam<unsigned int>("component")),
    _stress(getMaterialProperty<std::vector<RankTwoTensor>>("stress")),
    _B_mat(getMaterialProperty<std::vector<ColumnMajorMatrix>>("B_matrix"))
{
  _t_qrule = libmesh_make_unique<QGauss>(1, Utility::string_to_enum<Order>(getParam<std::string>("order")));
  _t_weights = _t_qrule->get_weights();

  printf("stress div done \n");
}

void
StressDivergenceShell::computeResidual()
{
  DenseVector<Number> & re = _assembly.residualBlock(_var.number());
  mooseAssert(re.size() == 4, "StressDivergenceShell: Beam element must have two nodes only.");
  _local_re.resize(re.size());
  _local_re.zero();

  _q_weights = _qrule->get_weights();
  precalculateResidual();
  for (_i = 0; _i < _test.size(); ++_i)
    for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
      for (_qp_z = 0; _qp_z < _stress[_qp].size(); ++_qp_z)
        _local_re(_i) += _q_weights[_qp] * _t_weights[_qp_z] * computeQpResidual();

  re += _local_re;

  if (_has_save_in)
  {
    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (_i = 0; _i < _save_in.size(); ++_i)
      _save_in[_i]->sys().solution().add_vector(_local_re, _save_in[_i]->dofIndices());
  }
}

Real
StressDivergenceShell::computeQpResidual()
{
  Real residual = _stress[_qp][_qp_z](0,0) * _B_mat[_qp][_qp_z](0,_i + _component*4) + _stress[_qp][_qp_z](1,1) * _B_mat[_qp][_qp_z](1,_i + _component*4) + 2.0 * _stress[_qp][_qp_z](0,1) * _B_mat[_qp][_qp_z](2,_i + _component*4) + 2.0 * _stress[_qp][_qp_z](0,2) * _B_mat[_qp][_qp_z](3,_i + _component*4) + 2.0 * _stress[_qp][_qp_z](1,2) * _B_mat[_qp][_qp_z](4,_i + _component*4);

  return residual;
}

void
StressDivergenceShell::computeJacobian()
{
  DenseMatrix<Number> & ke = _assembly.jacobianBlock(_var.number(), _var.number());
  _local_ke.resize(ke.m(), ke.n());
  _local_ke.zero();

  for (_i = 0; _i < _test.size(); ++_i)
    _local_ke(_i,_i) = 1.0;

  ke += _local_ke;

  if (_has_diag_save_in)
  {
    unsigned int rows = ke.m();
    DenseVector<Number> diag(rows);
    for (unsigned int i = 0; i < rows; ++i)
      diag(i) = _local_ke(i, i);

    Threads::spin_mutex::scoped_lock lock(Threads::spin_mtx);
    for (unsigned int i = 0; i < _diag_save_in.size(); ++i)
      _diag_save_in[i]->sys().solution().add_vector(diag, _diag_save_in[i]->dofIndices());
  }
}
