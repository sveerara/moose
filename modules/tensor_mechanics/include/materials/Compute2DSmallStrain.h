//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#ifndef COMPUTE2DSMALLSTRAIN_H
#define COMPUTE2DSMALLSTRAIN_H

#include "ComputeSmallStrain.h"

class Compute2DSmallStrain;

template <>
InputParameters validParams<Compute2DSmallStrain>();

/**
 * Compute2DSmallStrain defines a strain tensor, assuming small strains,
 * in 2D geometries / simulations.  ComputePlaneSmallStrain acts as a
 * base class for ComputePlaneSmallStrain and ComputeAxisymmetricRZSmallStrain
 * through the computeOutOfPlaneStrain method.
 */
class Compute2DSmallStrain : public ComputeSmallStrain
{
public:
  Compute2DSmallStrain(const InputParameters & parameters);

protected:
  void initialSetup() override;
  virtual void computeProperties() override;
  virtual void displacementIntegrityCheck() override;
  virtual Real computeOutOfPlaneStrain() = 0;

  /// Variable that specifies the out of plane direction
  const unsigned int _out_of_plane_direction;

  const bool _legacy_volumetric_locking_correction;

  /// Variable that specifies the in plane directions
  std::vector<unsigned int> _in_plane_direction;

  /// Whether to average the out of plane strain
  bool _ave_zz_strain;
};

#endif // COMPUTE2DSMALLSTRAIN_H
