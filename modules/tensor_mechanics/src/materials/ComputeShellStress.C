//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ComputeShellStress.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "ComputeIsotropicElasticityTensorShell.h"

registerMooseObject("TensorMechanicsApp", ComputeShellStress);

template <>
InputParameters
validParams<ComputeShellStress>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Compute stress using elasticity for finite strains");
  return params;
}

ComputeShellStress::ComputeShellStress(
    const InputParameters & parameters)
  : Material(parameters),
    _strain_increment(getMaterialPropertyByName<std::vector<RankTwoTensor>>( "strain_increment")),
    _stress(declareProperty<std::vector<RankTwoTensor>>("stress")),
    _stress_old(getMaterialPropertyOldByName<std::vector<RankTwoTensor>>("stress")),
    _elasticity_tensor(getMaterialProperty<std::vector<RankFourTensor>>("elasticity_tensor"))
{
  printf("stress done \n");
}

void
ComputeShellStress::initQpStatefulProperties()
{
  printf("in stress init \n");
  // initialize stress tensor to zero
  _stress[_qp].resize(_strain_increment[_qp].size());
  RankTwoTensor a;
  for (unsigned int i = 0; i < _strain_increment[_qp].size(); ++i)
    _stress[_qp][i] = a;

  printf("stress init done \n");
}

void
ComputeShellStress::computeQpProperties()
{
  for (unsigned int i = 0; i < _strain_increment[_qp].size(); ++i)
  {
    _stress[_qp][i] = _stress_old[_qp][i] + _elasticity_tensor[_qp][i] * _strain_increment[_qp][i];
  //  printf("stress , qp, t: %u, %u \n", _qp, i);
  //  _stress[_qp][i].print();
  }
}
