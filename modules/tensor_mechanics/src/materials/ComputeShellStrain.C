/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeShellStrain.h"
#include "MooseMesh.h"
#include "Assembly.h"
#include "NonlinearSystem.h"
#include "MooseVariable.h"
#include "ArbitraryQuadrature.h"
#include "ColumnMajorMatrix.h"

#include "libmesh/quadrature.h"
#include "libmesh/utility.h"
#include "libmesh/enum_quadrature_type.h"
#include "libmesh/fe_type.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/quadrature_gauss.h"

registerMooseObject("TensorMechanicsApp", ComputeShellStrain);

template <>
InputParameters
validParams<ComputeShellStrain>()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Compute a infinitesimal/large strain increment for the shell.");
  params.addRequiredCoupledVar(
      "rotations", "The rotations appropriate for the simulation geometry and coordinate system");
  params.addRequiredCoupledVar(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  params.addRequiredCoupledVar(
      "thickness",
      "Cross-section area of the beam. Can be supplied as either a number or a variable name.");
  params.addRequiredParam<std::string>("order", "Quadrature order in out of plane direction");
  return params;
}

ComputeShellStrain::ComputeShellStrain(const InputParameters & parameters)
  : Material(parameters),
    _nrot(coupledComponents("rotations")),
    _ndisp(coupledComponents("displacements")),
    _rot_num(_nrot),
    _disp_num(_ndisp),
    _thickness(coupledValue("thickness")),
    _strain_increment(declareProperty<std::vector<RankTwoTensor>>("strain_increment")),
    _total_strain(declareProperty<std::vector<RankTwoTensor>>("total_strain")),
    _total_strain_old(getMaterialProperty<std::vector<RankTwoTensor>>("total_strain")),
    _nonlinear_sys(_fe_problem.getNonlinearSystemBase()),
    _soln_disp_index(4),
    _soln_rot_index(4),
    _soln_vector(20, 1),
    _strain_vector(5,1),
    _node_normal(declareProperty<RealVectorValue>("node_normal")),
    _node_normal_old(getMaterialPropertyOldByName<RealVectorValue>("node_normal")),
    _V1(4),
    _V2(4),
    _B(declareProperty<std::vector<ColumnMajorMatrix>>("B_matrix")),
    _ge(declareProperty<std::vector<RankTwoTensor>>("ge_matrix"))
{
  // Checking for consistency between length of the provided displacements and rotations vector
  if (_ndisp != 3 || _nrot != 2)
    mooseError("ComputeShellStrain: The number of variables supplied in 'displacements' "
               "must be 3 and that in 'rotations' must be 2.");

  // fetch coupled variables and gradients (as stateful properties if necessary)
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    MooseVariable * disp_variable = getVar("displacements", i);
    _disp_num[i] = disp_variable->number();

    if (i < _nrot)
    {
      MooseVariable * rot_variable = getVar("rotations", i);
      _rot_num[i] = rot_variable->number();
    }
  }
  printf("strain done \n");
}

void
ComputeShellStrain::initQpStatefulProperties()
{
  printf("in strain initialization \n");
  _t_qrule = libmesh_make_unique<QGauss>(1, Utility::string_to_enum<Order>(getParam<std::string>("order")));
  _t_points = _t_qrule->get_points();

  // quadrature points in isoparametric space
  _2d_points = _qrule->get_points(); // would be in 2D

  unsigned int dim = _current_elem->dim();
  if ((dim != 2))
    mooseError("Shell element is implemented only for 2D Linear elements");

  // derivatives of shape functions (dphidxi, dphideta and dphidzeta) evaluated at quadrature points (in isoparametric space).
  FEType fe_type(Utility::string_to_enum<Order>("First"),
  Utility::string_to_enum<FEFamily>("LAGRANGE"));
  auto & fe = _subproblem.assembly(_tid).getFE(fe_type, dim);
  _dphidxi_map = fe->get_fe_map().get_dphidxi_map();
  _dphideta_map = fe->get_fe_map().get_dphideta_map();
  _phi_map = fe->get_fe_map().get_phi_map();

  printf("dphidxi map, dphideta_map, phi_map : %lu, %lu, %lu \n", _dphidxi_map.size(), _dphideta_map.size(), _phi_map.size());
  printf("dphidxi 0: %e, %e, %e, %e \n", _dphidxi_map[0][0], _dphidxi_map[0][1], _dphidxi_map[0][2], _dphidxi_map[0][3]);
  printf("dphidxi 1: %e, %e, %e, %e \n", _dphidxi_map[1][0], _dphidxi_map[1][1], _dphidxi_map[1][2], _dphidxi_map[1][3]);
  printf("dphidxi 2: %e, %e, %e, %e \n", _dphidxi_map[2][0], _dphidxi_map[2][1], _dphidxi_map[2][2], _dphidxi_map[2][3]);
  printf("dphidxi 3: %e, %e, %e, %e \n", _dphidxi_map[3][0], _dphidxi_map[3][1], _dphidxi_map[3][2], _dphidxi_map[3][3]);

  printf("dphideta 0: %e, %e, %e, %e \n", _dphideta_map[0][0], _dphideta_map[0][1], _dphideta_map[0][2], _dphideta_map[0][3]);
  printf("dphideta 1: %e, %e, %e, %e \n", _dphideta_map[1][0], _dphideta_map[1][1], _dphideta_map[1][2], _dphideta_map[1][3]);
  printf("dphideta 2: %e, %e, %e, %e \n", _dphideta_map[2][0], _dphideta_map[2][1], _dphideta_map[2][2], _dphideta_map[2][3]);
  printf("dphideta 3: %e, %e, %e, %e \n", _dphideta_map[3][0], _dphideta_map[3][1], _dphideta_map[3][2], _dphideta_map[3][3]);

  printf("phi 0: %e, %e, %e, %e \n", _phi_map[0][0], _phi_map[0][1], _phi_map[0][2], _phi_map[0][3]);
  printf("phi 1: %e, %e, %e, %e \n", _phi_map[1][0], _phi_map[1][1], _phi_map[1][2], _phi_map[1][3]);
  printf("phi 2: %e, %e, %e, %e \n", _phi_map[2][0], _phi_map[2][1], _phi_map[2][2], _phi_map[2][3]);
  printf("phi 3: %e, %e, %e, %e \n", _phi_map[3][0], _phi_map[3][1], _phi_map[3][2], _phi_map[3][3]);

  // Initialize node normals stored as material property but it is actually at the nodes.
  // So even if number of qp points is > 4, only the first 4 values will be accesssed.
  if (_2d_points.size() < 4)
    mooseError("ComputeShellStrain: Please use atleast 4 quadrature points in the planar direction.");

//  MooseVariable * disp_variable = getVar("displacements", 0);
//  const MooseArray<Point> & normal = disp_variable->normals();
// Todo: figure out how to get normals for now hard code normal
  RealVectorValue normal;
  normal(2) = 1.0;
  for (unsigned int k = 0; k < 4; ++k)
  {
    _node_normal[k] = normal;
  }

  _strain_increment[_qp].resize(_t_points.size());
  _total_strain[_qp].resize(_t_points.size());
  _B[_qp].resize(_t_points.size());
  _ge[_qp].resize(_t_points.size());

  _dxyz_dxi.resize(_2d_points.size());
  _dxyz_deta.resize(_2d_points.size());
  _dxyz_dzeta.resize(_2d_points.size());

  _dxyz_dxi[_qp].resize(_t_points.size());
  _dxyz_deta[_qp].resize(_t_points.size());
  _dxyz_dzeta[_qp].resize(_t_points.size());
  RankTwoTensor a;
  ColumnMajorMatrix b(5,20);
  RealVectorValue c;
  printf("done 3 \n");
  for (unsigned int t = 0; t < _t_points.size(); ++t)
  {
    printf("here \n");
    _strain_increment[_qp][t] = a;
    _total_strain[_qp][t] = a;
    _B[_qp][t] = b;
    _ge[_qp][t] = a;

    _dxyz_dxi[_qp][t] = c;
    _dxyz_deta[_qp][t] = c;
    _dxyz_dzeta[_qp][t] = c;
  }

  printf("strain init done \n");
}

void
ComputeShellStrain::computeProperties()
{
  // calculating derivatives of shape function is physical space (dphi/dx, dphi/dy, dphi/dz) at quadrature points
  // these are g_{i} in Dvorkin's paper
  std::vector<const Node *> nodes;
  for (unsigned int i = 0; i < 4; ++i)
    nodes.push_back(_current_elem->node_ptr(i));

  printf("nodes 0: %e %e \n", (*nodes[0])(0), (*nodes[0])(1));
  printf("nodes 1: %e %e \n", (*nodes[1])(0), (*nodes[1])(1));
  printf("nodes 2: %e %e \n", (*nodes[2])(0), (*nodes[2])(1));
  printf("nodes 3: %e %e \n", (*nodes[3])(0), (*nodes[3])(1));

  for (unsigned int i = 0; i < _2d_points.size(); ++i)
  {
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      for (unsigned int component = 0; component < 3; ++component)
      {
        _dxyz_dxi[i][j](component) = 0.0;
        _dxyz_deta[i][j](component) = 0.0;
        _dxyz_dzeta[i][j](component) = 0.0;
        for (unsigned int k = 0; k < nodes.size(); ++k)
        {
          _dxyz_dxi[i][j](component) += _dphidxi_map[k][i] * (*nodes[k])(component) + _t_points[j](0) / 2.0 * _thickness[i] * _dphidxi_map[k][i] * _node_normal_old[k](component);
          _dxyz_deta[i][j](component) += _dphideta_map[k][i] * (*nodes[k])(component) + _t_points[j](0) / 2.0 * _thickness[i] * _dphideta_map[k][i] * _node_normal_old[k](component);
          _dxyz_dzeta[i][j](component) += _thickness[i] * _phi_map[k][i] * _node_normal_old[k](component) / 2.0;
        }
      }
  /*    printf("dxyz_dxi: i, j, x, y, z: %u, %u, %e, %e, %e\n", i, j, _dxyz_dxi[i][j](0),_dxyz_dxi[i][j](1), _dxyz_dxi[i][j](2));
      printf("dxyz_deta: i, j, x, y, z: %u, %u, %e, %e, %e\n", i, j, _dxyz_deta[i][j](0),_dxyz_deta[i][j](1), _dxyz_deta[i][j](2));
      printf("dxyz_dzeta: i, j, x, y, z: %u, %u, %e, %e, %e\n", i, j, _dxyz_dzeta[i][j](0),_dxyz_dzeta[i][j](1), _dxyz_dzeta[i][j](2)); */
    }
  }
  /*    RankTwoTensor Jac;
      Jac(0,0) = _dxyz_dxi[i][j](0);
      Jac(0,1) = _dxyz_dxi[i][j](1);
      Jac(0,2) = _dxyz_dxi[i][j](2);
      Jac(1,0) = _dxyz_deta[i][j](0);
      Jac(1,1) = _dxyz_deta[i][j](1);
      Jac(1,2) = _dxyz_deta[i][j](2);
      Jac(2,0) = _dxyz_dzeta[i][j](0);
      Jac(2,1) = _dxyz_dzeta[i][j](1);
      Jac(2,2) = _dxyz_dzeta[i][j](2);

      _mapping_Jacobian[i][j] = Jac.det();
      RankTwoTensor Jacinv = Jac.inverse();
      // compute dphi_dx, dphi_dy and dphi_dz for each node at each qp
      for (unsigned int k = 0; k < node.size(); ++k)
      {
        _dphidx_map[k][i][j] = Jacinv(0,0) * _dphidxi_map[k][i] + Jacinv(0, 1) * _dphideta_map[k][i];
        _dphidy_map[k][i][j] = Jacinv(1,0) * _dphidxi_map[k][i] + Jacinv(1, 1) * _dphideta_map[k][i];
        _dphidz_map[k][i][j] = Jacinv(2,0) * _dphidxi_map[k][i] + Jacinv(2, 1) * _dphideta_map[k][i];

        _phi_dx_dzeta[k][i][j] = Jacinv(0,2) * _phi_map[k][i];
        _phi_dy_dzeta[k][i][j] = Jacinv(1,2) * _phi_map[k][i];
        _phi_dz_dzeta[k][i][j] = Jacinv(2,2) * _phi_map[k][i];
      }
    }
  }
 */

  // Fetch the incremental displacement (current - old) at the nodes
  const NumericVector<Number> & sol = *_nonlinear_sys.currentSolution();
  const NumericVector<Number> & sol_old = _nonlinear_sys.solutionOld();
  _soln_vector.zero();

  for (unsigned int j = 0; j < nodes.size(); ++j)
  {
    _soln_disp_index[j].resize(_ndisp);
    _soln_rot_index[j].resize(_nrot);

    for (unsigned int i = 0; i < _ndisp; ++i)
    {
      _soln_disp_index[j][i] = nodes[j]->dof_number(_nonlinear_sys.number(), _disp_num[i], 0);
      _soln_vector(j+i*nodes.size(), 0) = sol(_soln_disp_index[j][i]) - sol_old(_soln_disp_index[j][i]);
    }

    for (unsigned int i = 0; i < _nrot; ++i)
    {
      _soln_rot_index[j][i] = nodes[j]->dof_number(_nonlinear_sys.number(), _rot_num[i], 0);
      _soln_vector(j+12+i*nodes.size(),0) = sol(_soln_rot_index[j][i]) - sol_old(_soln_rot_index[j][i]);
    }
  }
  // compute nodal local axis
  RealGradient x2;
  x2(1) = 1;
  RealGradient x3;
  x3(2) = 1;

  for (unsigned int k = 0; k < nodes.size(); ++k)
  {
    _V1[k] = x2.cross(_node_normal_old[k]);
    _V1[k] /= x2.norm() * _node_normal_old[k].norm();

    // If x2 is parallel to node normal, set V1 to x3
    if (MooseUtils::absoluteFuzzyEqual(_V1[k].norm(), 0.0, 1e-6))
      _V1[k] = x3;

    _V2[k] = _node_normal_old[k].cross(_V1[k]);
  }

  // compute B matrix rows correspond to [ux1, ux2, ux3, ux4, uy1, uy2, uy3, uy4, uz1, uz2 uz3, uz4, a1, a2, a3, a4, b1, b2, b3, b4]
  for (unsigned int i = 0; i < _2d_points.size(); ++i)
  {
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
        for (unsigned int k = 0; k < nodes.size(); ++k)
        {
          // corresponding to strain(0,0)
          _B[i][j](0,k) = _dphidxi_map[k][i] * _dxyz_dxi[i][j](0);
          _B[i][j](0,4+k) = _dphidxi_map[k][i] * _dxyz_dxi[i][j](1);
          _B[i][j](0,8+k) = _dphidxi_map[k][i] * _dxyz_dxi[i][j](2);
          _B[i][j](0,12+k) = _dphidxi_map[k][i] * _t_points[j](0) / 2.0 * _thickness[i] * (-_V2[k] * _dxyz_dxi[i][j]);
          _B[i][j](0,16+k) = _dphidxi_map[k][i] * _t_points[j](0) / 2.0 * _thickness[i] * (_V1[k] * _dxyz_dxi[i][j]);

          // corresponding to strain(1,1)
          _B[i][j](1,k) = _dphideta_map[k][i] * _dxyz_deta[i][j](0);
          _B[i][j](1,4+k) = _dphideta_map[k][i] * _dxyz_deta[i][j](1);
          _B[i][j](1,8+k) = _dphideta_map[k][i] * _dxyz_deta[i][j](2);
          _B[i][j](1,12+k) = _dphideta_map[k][i] * _t_points[j](0) / 2.0 * _thickness[i] * (-_V2[k] * _dxyz_deta[i][j]);
          _B[i][j](1,16+k) = _dphideta_map[k][i] * _t_points[j](0) / 2.0 * _thickness[i] * (_V1[k] * _dxyz_deta[i][j]);

          // corresponding to strain(2,2) = 0

          //corresponding to strain(0,1)
          _B[i][j](2,k) = 0.5 * (_dphideta_map[k][i] * _dxyz_dxi[i][j](0) + _dphidxi_map[k][i] * _dxyz_deta[i][j](0));
          _B[i][j](2,4+k) = 0.5 * (_dphideta_map[k][i] * _dxyz_dxi[i][j](1) + _dphidxi_map[k][i] * _dxyz_deta[i][j](1));
          _B[i][j](2,8+k) = 0.5 * (_dphideta_map[k][i] * _dxyz_dxi[i][j](2) + _dphidxi_map[k][i] * _dxyz_deta[i][j](2));
          _B[i][j](2,12+k) = 0.25 * _t_points[j](0) * _thickness[i] * -_V2[k] * (_dphideta_map[k][i] * _dxyz_dxi[i][j] + _dphidxi_map[k][i] * _dxyz_deta[i][j]);
          _B[i][j](2,16+k) = 0.25 * _t_points[j](0) * _thickness[i] * _V1[k] * (_dxyz_deta[i][j] * _dphidxi_map[k][i] + _dxyz_dxi[i][j] * _dphideta_map[k][i]);
        }

        RealVectorValue g3_A = _thickness[i] / 4.0 * (_node_normal_old[2] + _node_normal_old[3]);
        RealVectorValue g3_C = _thickness[i] / 4.0 * (_node_normal_old[0] + _node_normal_old[1]);
        RealVectorValue g3_B = _thickness[i] / 4.0 * (_node_normal_old[0] + _node_normal_old[3]);
        RealVectorValue g3_D = _thickness[i] / 4.0 * (_node_normal_old[1] + _node_normal_old[2]);

        RealVectorValue g1_A = 0.5 * ((*nodes[2]) - (*nodes[3])) + _t_points[j](0)/4.0 * _thickness[i] * (_node_normal_old[2] - _node_normal_old[3]);
        RealVectorValue g1_C = 0.5 * ((*nodes[1]) - (*nodes[0])) + _t_points[j](0)/4.0 * _thickness[i] * (_node_normal_old[1] - _node_normal_old[0]);
        RealVectorValue g2_B = 0.5 * ((*nodes[3]) - (*nodes[0])) + _t_points[j](0)/4.0 * _thickness[i] * (_node_normal_old[3] - _node_normal_old[0]);
        RealVectorValue g2_D = 0.5 * ((*nodes[2]) - (*nodes[1])) + _t_points[j](0)/4.0 * _thickness[i] * (_node_normal_old[2] - _node_normal_old[1]);

      /*  printf("g1_A:%e, %e, %e \n", g1_A(0), g1_A(1), g1_A(2));
        printf("g1_C:%e, %e, %e \n", g1_C(0), g1_C(1), g1_C(2));
        printf("g2_B:%e, %e, %e \n", g2_B(0), g2_B(1), g2_B(2));
        printf("g2_D:%e, %e, %e \n", g2_D(0), g2_D(1), g2_D(2));
*/
        // corresponding to strain(0,2)
        for (unsigned int component = 0; component < 3; component ++)
        {
          _B[i][j](3,2+ component * 4) = 1.0/8.0 * (1.0 + _2d_points[i](1)) * g3_A(component);
          _B[i][j](3,3+ component * 4) = 1.0/8.0 * (1.0 + _2d_points[i](1)) * -g3_A(component);
          _B[i][j](3,1+ component * 4) = 1.0/8.0 * (1.0 - _2d_points[i](1)) * g3_C(component);
          _B[i][j](3,component * 4) = 1.0/8.0 * (1.0 - _2d_points[i](1)) * -g3_C(component);
        }
        _B[i][j](3, 14) = 1.0/8.0 * (1.0 + _2d_points[i](1)) * 0.5 * _thickness[i] * g1_A * -_V2[2];
        _B[i][j](3, 18) = 1.0/8.0 * (1.0 + _2d_points[i](1)) * 0.5 * _thickness[i] * g1_A * _V1[2];
        _B[i][j](3, 15) = 1.0/8.0 * (1.0 + _2d_points[i](1)) * 0.5 * _thickness[i] * g1_A * -_V2[3];
        _B[i][j](3, 19) = 1.0/8.0 * (1.0 + _2d_points[i](1)) * 0.5 * _thickness[i] * g1_A * _V1[3];

        _B[i][j](3, 13) = 1.0/8.0 * (1.0 - _2d_points[i](1)) * 0.5 * _thickness[i] * g1_C * -_V2[1];
        _B[i][j](3, 17) = 1.0/8.0 * (1.0 - _2d_points[i](1)) * 0.5 * _thickness[i] * g1_C * _V1[1];
        _B[i][j](3, 12) = 1.0/8.0 * (1.0 - _2d_points[i](1)) * 0.5 * _thickness[i] * g1_C * -_V2[0];
        _B[i][j](3, 16) = 1.0/8.0 * (1.0 - _2d_points[i](1)) * 0.5 * _thickness[i] * g1_C * _V1[0];

        // corresponding to strain(1,2)
        for (unsigned int component = 0; component < 3; component ++)
        {
          _B[i][j](4,2+ component * 4) = 1.0/8.0 * (1.0 + _2d_points[i](0)) * g3_D(component);
          _B[i][j](4,1+ component * 4) = 1.0/8.0 * (1.0 + _2d_points[i](0)) * -g3_D(component);
          _B[i][j](4,3+ component * 4) = 1.0/8.0 * (1.0 - _2d_points[i](0)) * g3_B(component);
          _B[i][j](4,component * 4) = 1.0/8.0 * (1.0 - _2d_points[i](0)) * -g3_B(component);
        }
        _B[i][j](4, 14) = 1.0/8.0 * (1.0 + _2d_points[i](0)) * 0.5 * _thickness[i] * g2_D * -_V2[2];
        _B[i][j](4, 18) = 1.0/8.0 * (1.0 + _2d_points[i](0)) * 0.5 * _thickness[i] * g2_D * _V1[2];
        _B[i][j](4, 13) = 1.0/8.0 * (1.0 + _2d_points[i](0)) * 0.5 * _thickness[i] * g2_D * -_V2[1];
        _B[i][j](4, 17) = 1.0/8.0 * (1.0 + _2d_points[i](0)) * 0.5 * _thickness[i] * g2_D * _V1[1];

        _B[i][j](4, 15) = 1.0/8.0 * (1.0 - _2d_points[i](0)) * 0.5 * _thickness[i] * g2_B * -_V2[3];
        _B[i][j](4, 19) = 1.0/8.0 * (1.0 - _2d_points[i](0)) * 0.5 * _thickness[i] * g2_B * _V1[3];
        _B[i][j](4, 12) = 1.0/8.0 * (1.0 - _2d_points[i](0)) * 0.5 * _thickness[i] * g2_B * -_V2[0];
        _B[i][j](4, 16) = 1.0/8.0 * (1.0 - _2d_points[i](0)) * 0.5 * _thickness[i] * g2_B * _V1[0];

        // compute strain increment in covariant coordinate system using B and _soln_vector
        _strain_vector = _B[i][j] * _soln_vector;
        _strain_increment[i][j](0,0) = _strain_vector(0,0);
        _strain_increment[i][j](1,1) = _strain_vector(1,0);
        _strain_increment[i][j](0,1) = _strain_vector(2,0);
        _strain_increment[i][j](0,2) = _strain_vector(3,0);
        _strain_increment[i][j](1,2) = _strain_vector(4,0);
        _strain_increment[i][j](1,0) = _strain_increment[i][j](0,1);
        _strain_increment[i][j](2,0) = _strain_increment[i][j](0,2);
        _strain_increment[i][j](2,1) = _strain_increment[i][j](1,2);
        _total_strain[i][j] = _total_strain_old[i][j] + _strain_increment[i][j];
        printf("strain_vector: %e, %e, %e, %e, %e \n", _strain_vector(0,0), _strain_vector(1,0), _strain_vector(2,0), _strain_vector(3,0), _strain_vector(4,0));
        // calculate gij for elasticity tensor
        RankTwoTensor gmn;
        for (unsigned int component = 0; component < 3; ++component)
        {
          gmn(0,0) += _dxyz_dxi[i][j](component) *  _dxyz_dxi[i][j](component);
          gmn(1,1) += _dxyz_deta[i][j](component) *  _dxyz_deta[i][j](component);
          gmn(2,2) += _dxyz_dzeta[i][j](component) *  _dxyz_dzeta[i][j](component);
          gmn(0,1) += _dxyz_dxi[i][j](component) *  _dxyz_deta[i][j](component);
          gmn(0,2) += _dxyz_dxi[i][j](component) *  _dxyz_dzeta[i][j](component);
          gmn(1,2) += _dxyz_deta[i][j](component) *  _dxyz_dzeta[i][j](component);
        }
        gmn(1,0) = gmn(0,1);
        gmn(2,0) = gmn(0,2);
        gmn(2,1) = gmn(1,2);

        RankTwoTensor gmninv = gmn.inverse();
        Real det = gmn.det();
        printf("det of g : %e \n", std::sqrt(det));
        // calculate ge
        RealVectorValue e3 = _dxyz_dzeta[i][j]/_dxyz_dzeta[i][j].norm();
        RealVectorValue e1 = _dxyz_deta[i][j].cross(e3)/_dxyz_deta[i][j].norm();
        RealVectorValue e2 = e3.cross(e1);

        _ge[i][j](0,0) = (gmninv * _dxyz_dxi[i][j]) * e1;
        _ge[i][j](0,1) = (gmninv * _dxyz_dxi[i][j]) * e2;
        _ge[i][j](0,2) = (gmninv * _dxyz_dxi[i][j]) * e3;
        _ge[i][j](1,0) = (gmninv * _dxyz_deta[i][j]) * e1;
        _ge[i][j](1,1) = (gmninv * _dxyz_deta[i][j]) * e2;
        _ge[i][j](1,2) = (gmninv * _dxyz_deta[i][j]) * e3;
        _ge[i][j](2,0) = (gmninv * _dxyz_dzeta[i][j]) * e1;
        _ge[i][j](2,1) = (gmninv * _dxyz_dzeta[i][j]) * e2;
        _ge[i][j](2,2) = (gmninv * _dxyz_dzeta[i][j]) * e3;

        printf("ge \n");
        _ge[i][j].print();
    }
  }
}


/* strain calc in covariant system without using B
  // compute strain increment in covariant coordinate system
  for (unsigned int i = 0; i < _3d_points.size(); ++i)
  {
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      _strain_increment[i][j](0, 0) = 0.0;
      _strain_increment[i][j](1, 1) = 0.0;
      _strain_increment[i][j](0, 1) = 0.0;
      RealGradient strain_11;
      RealGradient strain_22;
      for (unsigned int k = 0; k < _nodes.size(); ++k)
      {
          for (unsigned int component = 0; component < _ndisp.size(); ++component)
          {
            strain_00(component) += _dphidxi_map[k][i] * (_disp[k](component) + _t_points[j] / 2.0 * _thickness[i] * (-V2[k](component) * rot[k](0) + V1[k](component) * _rot[k](1));
            strain_11(component) += _dphideta_map[k][i] * (_disp[k](component) + _t_points[j] / 2.0 * _thickness[i] * (-V2[k](component) * rot[k](0) + V1[k](component) * _rot[k](1));
          }
      }
      for (unsigned int component = 0; component < _ndisp.size(); ++component)
      {
        _strain_increment[i][j](0,0) += _dxyz_dxi[i][j](component) * strain_00(component);
        _strain_increment[i][j](1,1) += _dxyz_deta[i][j](component) * strain_11(component);
        _strain_increment[i][j](0,1) += _dxyz_dxi[i][j](component) * strain_11(component) + _dxyz_deta[i][j](component) * strain_00(component);
      }
      _strain_increment[i][j](0,1)*= 0.5;
      _strain_increment[i][j](1,0) = _strain_increment[i][j](0,1);

      RealVectorValue g3_A = _thickness[i] / 4.0 * (_node_normal_old[2] + _node_normal_old[3]);
      RealVectorValue g3_C = _thickness[i] / 4.0 * (_node_normal_old[0] + _node_normal_old[1]);
      RealVectorValue g3_B = _thickness[i] / 4.0 * (_node_normal_old[0] + _node_normal_old[3]);
      RealVectorValue g3_D = _thickness[i] / 4.0 * (_node_normal_old[1] + _node_normal_old[2]);

      RealVectorValue g1_A = 0.5 * (node[2] - node[3]) + _t_points[j]/4.0 * (_node_normal_old[2] - _node_normal_old[3]);
      RealVectorValue g1_C = 0.5 * (node[1] - node[0]) + _t_points[j]/4.0 * (_node_normal_old[1] - _node_normal_old[0]);
      RealVectorValue g2_B = 0.5 * (node[3] - node[0]) + _t_points[j]/4.0 * (_node_normal_old[3] - _node_normal_old[0]);
      RealVectorValue g2_D = 0.5 * (node[2] - node[1]) + _t_points[j]/4.0 * (_node_normal_old[2] - _node_normal_old[1]);

      _strain_increment[i][j](0,2) = 0.0;
      _strain_increment[i][j](1,2) = 0.0;
      for (unsigned int component = 0; component < _ndisp.size(); ++component)
      {
        _strain_increment[i][j](0,2) += 1.0/8.0 * (1.0 + _3d_points[i][j](1)) * (g3_A(component) * (_disp[2](_component) - _disp[3](_component)) + 0.5 * thickness[i] * g1_A(component) * (-rot[2](0) * V2[2](component) + rot[2](1) * V1[2](component) -rot[3](0) * V2[3](component) + rot[3](1) * V1[3](component)))
        + 1.0/8.0 * (1.0 - _3d_points[i][j](1)) * (g3_C(component) * (_disp[1](_component) - _disp[0](_component)) + 0.5 * thickness[i] * g1_C(component) * (-rot[1](0) * V2[1](component) + rot[1](1) * V1[1](component) -rot[0](0) * V2[0](component) + rot[0](1) * V1[0](component)));

        _strain_increment[i][j](1,2) += 1.0/8.0 * (1.0 + _3d_points[i][j](0)) * (g3_D(component) * (_disp[2](_component) - _disp[1](_component)) + 0.5 * thickness[i] * g2_D(component) * (-rot[2](0) * V2[2](component) + rot[2](1) * V1[2](component) -rot[1](0) * V2[1](component) + rot[1](1) * V1[1](component)))
        + 1.0/8.0 * (1.0 - _3d_points[i][j](0)) * (g3_B(component) * (_disp[3](_component) - _disp[0](_component)) + 0.5 * thickness[i] * g2_B(component) * (-rot[3](0) * V2[3](component) + rot[3](1) * V1[3](component) -rot[0](0) * V2[0](component) + rot[0](1) * V1[0](component)));
      }
      _strain_increment[i][j](2,0) = _strain_increment[i][j](2,0);
      _strain_increment[i][j](2,1) = _strain_increment[i][j](1,2);

      _total_strain[i][j] = _total_strain_old[i][j] + _strain_increment[i][j];

      // calculate gij
      RankTwoTensor gmn;
      for (unsinged int component = 0; component < 3; ++component)
      {
        gmn(0,0) += _dxyz_dxi[i][j](component) *  _dxyz_dxi[i][j](component);
        gmn(1,1) += _dxyz_deta[i][j](component) *  _dxyz_deta[i][j](component);
        gmn(2,2) += _dxyz_dzeta[i][j](component) *  _dxyz_dzeta[i][j](component);
        gmn(0,1) += _dxyz_dxi[i][j](component) *  _dxyz_deta[i][j](component);
        gmn(0,2) += _dxyz_dxi[i][j](component) *  _dxyz_dzeta[i][j](component);
        gmn(1,2) += _dxyz_deta[i][j](component) *  _dxyz_dzeta[i][j](component);
      }
      gmn(1,0) = gmn(0,1);
      gmn(2,0) = gmn(0,2);
      gmn(2,1) = gmn(1,2);

      // calculate ge
      RealVectorValue e3 = _dxyz_dzeta[i][j]/_dxyz_deta[i][j].norm();
      RealVectorValue e1 = _dxyz_deta[i][j].cross(e3)/ (_dxyz_deta[i][j].cross(e3)).norm();
      RealVectorValue e2 = e3.cross(e1);

      _ge[i][j](0,0) = (gmn * _dxyz_dxi[i][j]) * e1;
      _ge[i][j](1,1) = (gmn * _dxyz_deta[i][j]) * e2;
      _ge[i][j](2,2) = (gmn * _dxyz_dzeta[i][j]) * e3;
      _ge[i][j](0,1) = (gmn * _dxyz_dxi[i][j]) * e2;
      _ge[i][j](0,2) = (gmn * _dxyz_dxi[i][j]) * e3;
      _ge[i][j](1,2) = (gmn * _dxyz_deta[i][j]) * e3;
      _ge[i][j](1,0) = (gmn * _dxyz_deta[i][j]) * e1;
      _ge[i][j](2,0) = (gmn * _dxyz_dzeta[i][j]) * e1;
      _ge[i][j](2,1) = (gmn * _dxyz_dzeta[i][j]) * e2;
    }
  }
*/

  /* Global coordinate strain
  // compute du/dx, du/dy and du/dz at qps for computing strain
  for (unsigned int i = 0; i < _3d_points.size(); ++i)
  {
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      for (unsigned int component = 0; component < _mesh.dimension(); ++component) //u_i
      {
        _du_dx[i][j][component] = 0.0;
        _du_dy[i][j][component] = 0.0;
        _du_dz[i][j][component] = 0.0;
        for (unsigned int k = 0; k < _nodes.size(); ++k)
        {
          _du_dx[i][j][component] +=  _dphidx_map[k][i][j] * _disp[k](component) + _thickness[i] / 2.0 * (_t_points[j](0) * _dphidx_map[k][i][j] + _phi_dx_dzeta[k][i][j]) * (-V2[k](component) * _rot[k](0) + V1[k](component) * _rot[k](1));
          _du_dy[i][j][component] +=  _dphidy_map[k][i][j] * _disp[k](component) + _thickness[i] / 2.0 * (_t_points[j](0) *  _dphidy_map[k][i][j] + _phi_dy_dzeta[k][i][j]) * (-V2[k](component) * _rot[k](0) + V1[k](component) * _rot[k](1));
          _du_dz[i][j][component] +=  _dphidz_map[k][i][j] * _disp[k](component) + _thickness[i] / 2.0 * (_t_points[j](0) *  _dphidz_map[k][i][j] + _phi_dz_dzeta[k][i][j]) * (-V2[k](component) * _rot[k](0) + V1[k](component) * _rot[k](1));
        }
      }

      _strain[i][j](0, 0) = _du_dx[i][j][0];
      _strain[i][j](1, 1) = _du_dy[i][j][1];
      _strain[i][j](2, 2) = _du_dz[i][j][2];
      _strain[i][j](1, 2) = 0.5 * (_du_dx[i][j][1] + _du_dy[i][j][0]);
      _strain[i][j](2, 1) = _strain[i][j](1, 2);
      _strain[i][j](0, 2) = 0.5 * (_du_dx[i][j][2] + _du_dz[i][j][0]);
      _strain[i][j](1, 2) = 0.5 * (_du_dy[i][j][2] + _du_dz[i][j][1]);
    }
  }

  Real avg_strain_0_2 = 0.0;
  Real avg_strain_1_2 = 0.0;
  Real vol = 0.0;
  for (unsigned int i = 0; i < _3d_points.size(); ++i)
  {
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      avg_strain_0_2 += _strain[i][j](0,2) * _wi[i] * _wj[j] * _mapping_Jacobian[i][j];
      avg_strain_1_2 += _strain[i][j](1,2) * _wi[i] * _wj[j] * _mapping_Jacobian[i][j];
      vol += _wi[i] * _wj[j] * _mapping_Jacobian[i][j];
    }
  }

  avg_strain_0_2 /= vol;
  avg_strain_1_2 /= vol;

  for (unsigned int i = 0; i < _3d_points.size(); ++i)
  {
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      _strain[i][j](0, 2) = avg_strain_0_2;
      _strain[i][j](1, 2) = avg_strain_1_2;
      _strain[i][j](2, 0) = avg_strain_0_2;
      _strain[i][j](2, 1) = avg_strain_1_2;
    }
  }
} */
