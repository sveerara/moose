//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef STRESSDIVERGENCESHELL_H
#define STRESSDIVERGENCESHELL_H

#include "Kernel.h"
#include "RankTwoTensor.h"
#include "RankFourTensor.h"
#include "libmesh/quadrature_gauss.h"

// Forward Declarations
class StressDivergenceShell;

namespace libMesh
{
  class QGauss;
}

template <>
InputParameters validParams<StressDivergenceShell>();

/**
 * StressDivergenceTensors mostly copies from StressDivergence.  There are small changes to use
 * RankFourTensor and RankTwoTensors instead of SymmElasticityTensors and SymmTensors.  This is done
 * to allow for more mathematical transparancy.
 */
class StressDivergenceShell : public Kernel
{
public:
  StressDivergenceShell(const InputParameters & parameters);

  virtual void computeJacobian() override;

protected:
  virtual void computeResidual() override;
  virtual Real computeQpResidual() override;

  const unsigned int _component;

  const MaterialProperty<std::vector<RankTwoTensor>> & _stress;
  const MaterialProperty<std::vector<ColumnMajorMatrix>> & _B_mat;

  /// Quadrature rule in the out of plane direction
  std::unique_ptr<QGauss> _t_qrule;

  /// Quadrature points in the out of plane direction in isoparametric coordinate system
  std::vector<Real> _t_weights;

  /// Qrule weights in isoparametric coordinate system
  std::vector<Real> _q_weights;

  /// qp index in out of plane direction
  unsigned int _qp_z;
};

#endif // STRESSDIVERGENCESHELL_H
