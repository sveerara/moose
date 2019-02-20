//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTESHELLSTRESS_H
#define COMPUTESHELLSTRESS_H

#include "Material.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "ComputeIsotropicElasticityTensorShell.h"

/**
 * ComputeShellStress is the base class for stress tensors
 */

class ComputeShellStress;

template <>
InputParameters validParams<ComputeShellStress>();

class ComputeShellStress : public Material
{
public:
  ComputeShellStress(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpProperties() override;

  /// Material property for strain increment
  const MaterialProperty<std::vector<RankTwoTensor>> & _strain_increment;

  /// Material property for current stress
  MaterialProperty<std::vector<RankTwoTensor>> & _stress;

  /// Material property for old stress
  const MaterialProperty<std::vector<RankTwoTensor>> & _stress_old;

  /// Material property for elasticity tensor
  const MaterialProperty<std::vector<RankFourTensor>> & _elasticity_tensor;
};

#endif // COMPUTESHELLSTRESS_H
