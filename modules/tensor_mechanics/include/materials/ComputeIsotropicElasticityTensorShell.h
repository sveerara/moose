//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTEISOTROPICELASTICITYTENSORSHELL_H
#define COMPUTEISOTROPICELASTICITYTENSORSHELL_H

#include "Material.h"

class ComputeIsotropicElasticityTensorShell;

template <>
InputParameters validParams<ComputeIsotropicElasticityTensorShell>();

/**
 * ComputeIsotropicElasticityTensor defines an elasticity tensor material for
 * isotropic materials.
 */
class ComputeIsotropicElasticityTensorShell : public Material
{
public:
  ComputeIsotropicElasticityTensorShell(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  Real _poissons_ratio;
  Real _shear_modulus;
  Real _youngs_modulus;

  /// Individual elasticity tensor
  RankFourTensor _Cijkl;

  /// Material property elasticity tensor
  MaterialProperty<std::vector<RankFourTensor>> & _elasticity_tensor;

  /// Material property for ge matrix
  const MaterialProperty<std::vector<RankTwoTensor>> & _ge;
};

#endif // COMPUTEISOTROPICELASTICITYTENSORSHELL_H
