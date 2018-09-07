/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "NodalMassMatrix.h"
#include "MooseVariable.h"
#include "AuxiliarySystem.h"
#include "MooseUtils.h"
#include "DelimitedFileReader.h"

registerMooseObject("TensorMechanicsApp", NodalMassMatrix);

template <>
InputParameters
validParams<NodalMassMatrix>()
{
  InputParameters params = validParams<NodalKernel>();
  params.addClassDescription("Computes the interial forces and mass proportional damping terms "
                             "corresponding to nodal mass.");
  params.addParam<Real>("mass", "Mass associated with the node");
  params.addParam<FileName>(
      "nodal_mass_file",
      "The file containing the nodal positions and the corresponding nodal masses.");
  return params;
}

NodalMassMatrix::NodalMassMatrix(const InputParameters & parameters)
  : NodalKernel(parameters),
    _mass(isParamValid("mass") ? getParam<Real>("mass") : 0.0)
{
  if (!isParamValid("nodal_mass_file") && !isParamValid("mass"))
    mooseError(
        "NodalMassMatrix: Please provide either mass or nodal_mass_file as input.");
  else if (isParamValid("nodal_mass_file") && isParamValid("mass"))
    mooseError("NodalMassMatrix: Please provide either mass or nodal_mass_file as input, "
               "not both.");

  if (isParamValid("nodal_mass_file"))
  {
    MooseUtils::DelimitedFileReader nodal_mass_file(getParam<FileName>("nodal_mass_file"));
    nodal_mass_file.setHeaderFlag(MooseUtils::DelimitedFileReader::HeaderFlag::OFF);
    nodal_mass_file.read();
    std::vector<std::vector<Real>> data = nodal_mass_file.getData();
    if (data.size() != 4)
      mooseError("NodalMassMatrix: The number of columns in ",
                 getParam<FileName>("nodal_mass_file"),
                 " should be 4.");

    unsigned int node_found = 0;
    const std::set<BoundaryID> bnd_ids = BoundaryRestrictable::boundaryIDs();
    for (auto & bnd_id : bnd_ids)
    {
      const std::vector<dof_id_type> & bnd_node_set = _mesh.getNodeList(bnd_id);
      for (auto & bnd_node : bnd_node_set)
      {
        const Node & node = _mesh.nodeRef(bnd_node);
        _node_id_to_mass[bnd_node] = 0.0;

        for (unsigned int i = 0; i < data[0].size(); ++i)
        {
          if (MooseUtils::absoluteFuzzyEqual(data[0][i], node(0), 1e-6) &&
              MooseUtils::absoluteFuzzyEqual(data[1][i], node(1), 1e-6) &&
              MooseUtils::absoluteFuzzyEqual(data[2][i], node(2), 1e-6))
          {
            _node_id_to_mass[bnd_node] = data[3][i];
            node_found += 1;
            break;
          }
        }
      }
    }
    if (node_found != data[0].size())
      mooseError("NodalMassMatrix: Out of ",
                 data[0].size(),
                 " nodal positions in ",
                 getParam<FileName>("nodal_mass_file"),
                 " only ",
                 node_found,
                 " nodes were found in the boundary.");
  }
}

Real
NodalMassMatrix::computeQpResidual()
{
  Real mass = 0.0;
  if (isParamValid("mass"))
    mass = _mass;
  else
    mass = _node_id_to_mass[_current_node->id()];

  return mass * _u[_qp];
}

Real
NodalMassMatrix::computeQpJacobian()
{
  Real mass = 0.0;
  if (isParamValid("mass"))
    mass = _mass;
  else
    mass = _node_id_to_mass[_current_node->id()];
  return mass;
}
