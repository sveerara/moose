//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeShellStress.h"

registerMooseObject("TensorMechanicsApp", ComputeShellStress);

template <>
InputParameters
validParams<ComputeShellStress>()
{
  InputParameters params = validParams<ComputeStressBase>();
  params.addClassDescription("Compute stress using elasticity for finite strains");
  return params;
}

ComputeShellStress::ComputeShellStress(
    const InputParameters & parameters)
  : ComputeStressBase(parameters),
    GuaranteeConsumer(this),
    _strain_increment(getMaterialPropertyByName<RankTwoTensor>(_base_name + "strain_increment")),
    _rotation_increment(
        getMaterialPropertyByName<RankTwoTensor>(_base_name + "rotation_increment")),
    _stress_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "stress")),
    _elastic_strain_old(getMaterialPropertyOldByName<RankTwoTensor>(_base_name + "elastic_strain"))
{
}

void
ComputeShellStress::initialSetup()
{
  if (!hasGuaranteedMaterialProperty(_elasticity_tensor_name, Guarantee::ISOTROPIC))
    mooseError("ComputeShellStress can only be used with elasticity tensor materials "
               "that guarantee isotropic tensors.");
}

void
ComputeShellStress::computeQpStress()
{
  for (unsigned int t = 0; t < _t_points.size(); ++t)
    _stress[_qp][t] = _elasticity_tensor[_qp][t] * _total_strain[_qp][t];

  // Compute dstress_dstrain
  _Jacobian_mult[_qp] = (_elasticity_tensor[_qp][0] + _elasticity_tensor[_qp][1]); // This is NOT the exact jacobian
}
