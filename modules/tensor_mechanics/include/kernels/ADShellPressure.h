//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ADKernelValue.h"

template <ComputeStage>
class ADShellPressure;

declareADValidParams(ADShellPressure);

template <ComputeStage compute_stage>
class ADShellPressure : public ADKernelValue<compute_stage>
{
public:
  ADShellPressure(const InputParameters & parameters);

protected:
  ADReal precomputeQpResidual() override;

private:
  const Real _value;

  const unsigned int _component;

  ADRealVectorValue _normal;

  usingKernelValueMembers;
};
