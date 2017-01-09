/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTESMALLSTRAIN_H
#define COMPUTESMALLSTRAIN_H

#include "ComputeStrainBase.h"

/**
 * ComputeSmallStrain defines a strain tensor, assuming small strains.
 */
class ComputeSmallStrain : public ComputeStrainBase
{
public:
  ComputeSmallStrain(const InputParameters & parameters);

protected:
  virtual void computeProperties() override;

  virtual void computeEnhancedStrain();
  virtual void computeEnhancedStrainIncrement(unsigned int qp, RankTwoTensor & strain_increment);
  const MaterialProperty<RankFourTensor> & _Jacobian_mult;
  const bool _enhanced_strain;
  const unsigned int _disp_x;
  const unsigned int _disp_y;
  const unsigned int _disp_z;
  VariablePhiGradient _dphi;
  ColumnMajorMatrix _a;
  unsigned int _dim;
  unsigned int _mesh_dim;
};

#endif //COMPUTESMALLSTRAIN_H
