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
  params.addParam<bool>("large_strain", false, "Set to true to turn on finite strain calculations.");
  return params;
}

ComputeShellStrain::ComputeShellStrain(const InputParameters & parameters)
  : Material(parameters),
    _nrot(coupledComponents("rotations")),
    _ndisp(coupledComponents("displacements")),
    _rot_num(_nrot),
    _disp_num(_ndisp),
    _thickness(coupledValue("thickness")),
    _large_strain(getParam<bool>("large_strain")),
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
    _BNL_new(declareProperty<std::vector<ColumnMajorMatrix>>("BNLnew_matrix")),
    _BNL(_large_strain? &declareProperty<std::vector<ColumnMajorMatrix>>("BNL_matrix"):nullptr),
    _BNL_old(_large_strain? &getMaterialPropertyOld<std::vector<ColumnMajorMatrix>>("BNL_matrix"):nullptr),
    _ge(declareProperty<std::vector<RankTwoTensor>>("ge_matrix")),
    _Jmap(declareProperty<std::vector<Real>>("J_mapping")),
    _soln_vector_prop(declareProperty<std::vector<Real>>("soln_vector"))
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
}

void
ComputeShellStrain::initQpStatefulProperties()
{
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

  // Initialize node normals stored as material property but it is actually at the nodes.
  // So even if number of qp points is > 4, only the first 4 values will be accesssed.
  if (_2d_points.size() < 4)
    mooseError("ComputeShellStrain: Please use atleast 4 quadrature points in the planar direction.");

//  MooseVariable * disp_variable = getVar("displacements", 0);
//  const MooseArray<Point> & normal = disp_variable->normals();
// Todo: figure out how to get normals for now hard code normal
  for (unsigned int i = 0; i < 4; ++i)
     _nodes.push_back(_current_elem->node_ptr(i));

  RealVectorValue x = (*_nodes[1]-*_nodes[0]);
  RealVectorValue y = (*_nodes[3]-*_nodes[0]);
  RealVectorValue normal = x.cross(y);
  normal /= normal.norm();

  //normal(2) = 1.0;
  for (unsigned int k = 0; k < 4; ++k)
  {
    _node_normal[k] = normal;
  }

  _soln_vector_prop[_qp].resize(20);
  _strain_increment[_qp].resize(_t_points.size());
  _total_strain[_qp].resize(_t_points.size());
  _B[_qp].resize(_t_points.size());
  _BNL_new[_qp].resize(_t_points.size());
  _ge[_qp].resize(_t_points.size());
  _Jmap[_qp].resize(_t_points.size());

  _dxyz_dxi.resize(_2d_points.size());
  _dxyz_deta.resize(_2d_points.size());
  _dxyz_dzeta.resize(_2d_points.size());

  _dxyz_dxi[_qp].resize(_t_points.size());
  _dxyz_deta[_qp].resize(_t_points.size());
  _dxyz_dzeta[_qp].resize(_t_points.size());

  if (_large_strain)
    (*_BNL)[_qp].resize(_t_points.size());

  RankTwoTensor a;
  ColumnMajorMatrix b(5,20);
  ColumnMajorMatrix d(6,20);
  RealVectorValue c;

  for (unsigned int t = 0; t < _t_points.size(); ++t)
  {
    _strain_increment[_qp][t] = a;
    _total_strain[_qp][t] = a;
    _B[_qp][t] = b;
    _BNL_new[_qp][t] = d;
    _ge[_qp][t] = a;
    _Jmap[_qp][t] = 0;
    if (_large_strain)
      (*_BNL)[_qp][t] = b;

    _dxyz_dxi[_qp][t] = c;
    _dxyz_deta[_qp][t] = c;
    _dxyz_dzeta[_qp][t] = c;
  }
}

void
ComputeShellStrain::computeProperties()
{
  // calculating derivatives of shape function is physical space (dphi/dx, dphi/dy, dphi/dz) at quadrature points
  // these are g_{i} in Dvorkin's paper
  std::vector<const Node *> nodes;
  for (unsigned int i = 0; i < 4; ++i)
    nodes.push_back(_current_elem->node_ptr(i));

  RealVectorValue en;
  RealVectorValue ex;
  RealVectorValue ey;
  en(2) = 1.0;
  ex(0) = 1.0;
  ey(1) = 1.0;
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
          _dxyz_dxi[i][j](component) += _dphidxi_map[k][i] * ((*nodes[k])(component) /*+ sol_old(_soln_disp_index[k][component]) */) + _t_points[j](0) / 2.0 * _thickness[i] * _dphidxi_map[k][i] * en(component);
          _dxyz_deta[i][j](component) += _dphideta_map[k][i] * ((*nodes[k])(component) /*+ sol_old(_soln_disp_index[k][component]) */) + _t_points[j](0) / 2.0 * _thickness[i] * _dphideta_map[k][i] * en(component);
          _dxyz_dzeta[i][j](component) += _thickness[i] * _phi_map[k][i] * en(component) / 2.0;
        }
      }
    }
  }

  for (unsigned int i = 0; i < _2d_points.size(); ++i)
  {
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
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
      _Jmap[i][j] = std::sqrt(gmn.det());

      // calculate ge
      RealVectorValue e3 = _dxyz_dzeta[i][j]/_dxyz_dzeta[i][j].norm();
      RealVectorValue e1 = _dxyz_deta[i][j].cross(e3);
      e1 /= e1.norm();
      RealVectorValue e2 = e3.cross(e1);
      e2 /= e2.norm();

      _ge[i][j](0,0) = (gmninv * _dxyz_dxi[i][j]) * e1;
      _ge[i][j](0,1) = (gmninv * _dxyz_dxi[i][j]) * e2;
      _ge[i][j](0,2) = (gmninv * _dxyz_dxi[i][j]) * e3;
      _ge[i][j](1,0) = (gmninv * _dxyz_deta[i][j]) * e1;
      _ge[i][j](1,1) = (gmninv * _dxyz_deta[i][j]) * e2;
      _ge[i][j](1,2) = (gmninv * _dxyz_deta[i][j]) * e3;
      _ge[i][j](2,0) = (gmninv * _dxyz_dzeta[i][j]) * e1;
      _ge[i][j](2,1) = (gmninv * _dxyz_dzeta[i][j]) * e2;
      _ge[i][j](2,2) = (gmninv * _dxyz_dzeta[i][j]) * e3;

    //  _ge[i][j] = (_ge[i][j] + _ge[i][j].transpose())/2.0;
    }
  }

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

  for (unsigned int temp = 0; temp < 20; ++temp)
    _soln_vector_prop[0][temp] = _soln_vector(temp, 0);
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

    // update _node_normal
    if (_large_strain)
    {
      _node_normal[k] = -_V2[k] * _soln_vector(12+k) + _V1[k] * _soln_vector(16+k) + _node_normal_old[k];
      _node_normal[k] /= _node_normal[k].norm();
    }
//    printf("node_normal, k: %u, %e, %e, %e \n", k, _node_normal[k](0), _node_normal[k](1), _node_normal[k](2));
  }

  if (_large_strain)
  {
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
          _dxyz_dxi[i][j](component) += _dphidxi_map[k][i] * ((*nodes[k])(component) + sol_old(_soln_disp_index[k][component]) ) + _t_points[j](0) / 2.0 * _thickness[i] * _dphidxi_map[k][i] * _node_normal_old[k](component);
          _dxyz_deta[i][j](component) += _dphideta_map[k][i] * ((*nodes[k])(component) + sol_old(_soln_disp_index[k][component]) ) + _t_points[j](0) / 2.0 * _thickness[i] * _dphideta_map[k][i] * _node_normal_old[k](component);
          _dxyz_dzeta[i][j](component) += _thickness[i] * _phi_map[k][i] * _node_normal_old[k](component) / 2.0;
        }
      }
    }
  }
}
  // compute B matrix rows correspond to [ux1, ux2, ux3, ux4, uy1, uy2, uy3, uy4, uz1, uz2, uz3, uz4, a1, a2, a3, a4, b1, b2, b3, b4]
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

        if (_large_strain)
      {
        for (unsigned int component = 0; component < 3; ++component)
        {
          g1_A(component) += 0.5 * (sol_old(_soln_disp_index[2][component]) - sol_old(_soln_disp_index[3][component]));
          g1_C(component) += 0.5 * (sol_old(_soln_disp_index[1][component]) - sol_old(_soln_disp_index[0][component]));
          g2_B(component) += 0.5 * (sol_old(_soln_disp_index[3][component]) - sol_old(_soln_disp_index[0][component]));
          g2_D(component) += 0.5 * (sol_old(_soln_disp_index[2][component]) - sol_old(_soln_disp_index[1][component]));
        }
      }

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

        // compute BNL if large strain is required
        if (_large_strain)
          computeLargeStrain();

        // compute strain increment in covariant coordinate system using B and _soln_vector
        _strain_vector = _B[i][j] * _soln_vector;

      //  if (_large_strain)
      //    _strain_vector += (*_BNL)[i][j] * _soln_vector;

        _strain_increment[i][j](0,0) = _strain_vector(0,0);
        _strain_increment[i][j](1,1) = _strain_vector(1,0);
        _strain_increment[i][j](0,1) = _strain_vector(2,0);
        _strain_increment[i][j](0,2) = _strain_vector(3,0);
        _strain_increment[i][j](1,2) = _strain_vector(4,0);
        _strain_increment[i][j](1,0) = _strain_increment[i][j](0,1);
        _strain_increment[i][j](2,0) = _strain_increment[i][j](0,2);
        _strain_increment[i][j](2,1) = _strain_increment[i][j](1,2);
        _total_strain[i][j] = _total_strain_old[i][j] + _strain_increment[i][j];
    }
  }
}

void
ComputeShellStrain::computeLargeStrain()
{
  // compute BNL matrix - rows correspond to [ux1, ux2, ux3, ux4, uy1, uy2, uy3, uy4, uz1, uz2, uz3, uz4, a1, a2, a3, a4, b1, b2, b3, b4]
/*  _soln_vector.zero();
  for (unsigned int temp = 0; temp < 20; ++temp)
    _soln_vector(temp,0) = temp+1.0; */

  for (unsigned int i = 0; i < _2d_points.size(); ++i)
  {
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      (*_BNL)[i][j].zero();
      for (unsigned int k = 0; k < 4; ++k)
      {
        for (unsigned int p = 0; p < 4; ++p) // loop over nodes
        {
          // corresponding to strain(0,0)
          (*_BNL)[i][j](0,k) += _dphidxi_map[k][i] * _dphidxi_map[p][i] * (_soln_vector(p,0)
                             + _t_points[j](0)/2.0 * _thickness[i] *(-_soln_vector(p+12,0) * _V2[p](0)
                             + _soln_vector(p+16,0) * _V1[p](0)));
          (*_BNL)[i][j](0,4+k) += _dphidxi_map[k][i] * _dphidxi_map[p][i] * (_soln_vector(p+4,0)
                               + _t_points[j](0)/2.0 * _thickness[i] *(-_soln_vector(p+12,0) * _V2[p](1)
                               + _soln_vector(p+16,0) * _V1[p](1)));
          (*_BNL)[i][j](0,8+k) += _dphidxi_map[k][i] * _dphidxi_map[p][i] * (_soln_vector(p+8,0)
                               + _t_points[j](0)/2.0 * _thickness[i] *(-_soln_vector(p+12,0) * _V2[p](2)
                               + _soln_vector(p+16,0) * _V1[p](2)));
          (*_BNL)[i][j](0,12+k) += _t_points[j](0)/2.0 * _thickness[i] * _dphidxi_map[k][i] * _dphidxi_map[p][i] * (- (_V2[p](0) * _soln_vector(p,0) + _V2[p](1) * _soln_vector(p+4,0) + _V2[p](2) * _soln_vector(p+8,0)) +  _t_points[j](0) / 2.0 * _thickness[i] * _V2[k]
                                * (_V2[p] * _soln_vector(p+12,0) - _V1[p] * _soln_vector(p+16,0)));
          (*_BNL)[i][j](0,16+k) += _t_points[j](0)/2.0 * _thickness[i] * _dphidxi_map[k][i] * _dphidxi_map[p][i] * ((_V1[p](0) * _soln_vector(p,0) + _V1[p](1) * _soln_vector(p+4,0) + _V1[p](2) * _soln_vector(p+8,0)) +  _t_points[j](0) / 2.0 * _thickness[i] * _V1[k]
                                * (-_V2[p] * _soln_vector(p+12,0) + _V1[p] * _soln_vector(p+16,0)));

          // corresponding to strain(1,1)
          (*_BNL)[i][j](1,k) += _dphideta_map[k][i] * _dphideta_map[p][i] * (_soln_vector(p,0)
                            + _t_points[j](0)/2.0 * _thickness[i] *(-_soln_vector(p+12,0) * _V2[p](0)
                            + _soln_vector(p+16,0) * _V1[p](0)));
          (*_BNL)[i][j](1,4+k) += _dphideta_map[k][i] * _dphideta_map[p][i] * (_soln_vector(p+4,0)
                              + _t_points[j](0)/2.0 * _thickness[i] *(-_soln_vector(p+12,0) * _V2[p](1)
                              + _soln_vector(p+16,0) * _V1[p](1)));
          (*_BNL)[i][j](1,8+k) += _dphideta_map[k][i] * _dphideta_map[p][i] * (_soln_vector(p+8,0)
                               + _t_points[j](0)/2.0 * _thickness[i] *(-_soln_vector(p+12,0) * _V2[p](2)
                               + _soln_vector(p+16,0) * _V1[p](2)));
          (*_BNL)[i][j](1,12+k) += _t_points[j](0)/2.0 * _thickness[i] * _dphideta_map[k][i] * _dphideta_map[p][i] * (- (_V2[p](0) * _soln_vector(p,0) + _V2[p](1) * _soln_vector(p+4,0) + _V2[p](2) * _soln_vector(p+8,0)) +  _t_points[j](0) / 2.0 * _thickness[i] * _V2[k]
                                * (_V2[p] * _soln_vector(p+12,0) - _V1[p] * _soln_vector(p+16,0)));
          (*_BNL)[i][j](1,16+k) += _t_points[j](0)/2.0 * _thickness[i] * _dphideta_map[k][i] * _dphideta_map[p][i] * ((_V1[p](0) * _soln_vector(p,0) + _V1[p](1) * _soln_vector(p+4,0) + _V1[p](2) * _soln_vector(p+8,0)) +  _t_points[j](0) / 2.0 * _thickness[i] * _V1[k]
                                * (-_V2[p] * _soln_vector(p+12,0) + _V1[p] * _soln_vector(p+16,0)));

          // terms corresponding to strain(2,2) are 0.

          // corresponding to strain(0,1)
          (*_BNL)[i][j](2,k) += 0.5 *(_dphidxi_map[k][i] * _dphideta_map[p][i] + _dphideta_map[k][i] * _dphidxi_map[p][i]) * (_soln_vector(p,0)
                             + _t_points[j](0) / 2.0 * _thickness[i] *(-_soln_vector(p+12,0) * _V2[p](0)
                             + _soln_vector(p+16,0) * _V1[p](0)));
          (*_BNL)[i][j](2,4+k) += 0.5 *(_dphidxi_map[k][i] * _dphideta_map[p][i] + _dphideta_map[k][i] * _dphidxi_map[p][i]) * (_soln_vector(p+4,0)
                              + _t_points[j](0)/2.0 * _thickness[i] *(-_soln_vector(p+12,0) * _V2[p](1)
                              + _soln_vector(p+16,0) * _V1[p](1)));
          (*_BNL)[i][j](2,8+k) +=  0.5 *(_dphidxi_map[k][i] * _dphideta_map[p][i] + _dphideta_map[k][i] * _dphidxi_map[p][i]) * (_soln_vector(p+8,0)
                              + _t_points[j](0)/2.0 * _thickness[i] *(-_soln_vector(p+12,0) * _V2[p](2)
                              + _soln_vector(p+16,0) * _V1[p](2)));
          (*_BNL)[i][j](2,12+k) += _t_points[j](0) * 0.25 *(_dphidxi_map[k][i] * _dphideta_map[p][i] + _dphideta_map[k][i] * _dphidxi_map[p][i]) * _thickness[i]
                              * (-(_V2[k](0) * _soln_vector(p,0) + _V2[k](1) * _soln_vector(p+4,0) + _V2[k](2) * _soln_vector(p+8,0))
                              + _t_points[j](0) / 2.0 * _thickness[i] * _V2[k]
                              * (_V2[p] * _soln_vector(p+12,0) - _V1[p] * _soln_vector(p+16,0)));
          (*_BNL)[i][j](2,16+k) +=  _t_points[j](0) * 0.25 *(_dphidxi_map[k][i] * _dphideta_map[p][i] + _dphideta_map[k][i] * _dphidxi_map[p][i]) * _thickness[i]
                              * ((_V1[k](0) * _soln_vector(p,0) + _V1[k](1) * _soln_vector(p+4,0) + _V1[k](2) * _soln_vector(p+8,0))
                              + _t_points[j](0) / 2.0 * _thickness[i] * _V1[k]
                              * (-_V2[p] * _soln_vector(p+12,0) + _V1[p] * _soln_vector(p+16,0)));
        }
      }

      for (unsigned int component = 0; component < 3; ++component)
      {
        // corresponding to strain(0,2)
        (*_BNL)[i][j](3,2 + component * 4) += 1.0/32.0 * (1.0 + _2d_points[i](1)) * _thickness[i] * (-_soln_vector(12 + 2,0) * _V2[2](component) + _soln_vector(16 + 2,0) * _V1[2](component) -_soln_vector(12 + 3,0) * _V2[3](component) + _soln_vector(16 + 3,0) * _V1[3](component));
        (*_BNL)[i][j](3, 3 + component * 4) += - (*_BNL)[i][j](3, 2 + component * 4);

        (*_BNL)[i][j](3,1 + component * 4) += 1.0/32.0 * (1.0 - _2d_points[i](1)) * _thickness[i] * (-_soln_vector(12 + 1,0) * _V2[1](component) + _soln_vector(16 + 1,0) * _V1[1](component) -_soln_vector(12 + 0,0) * _V2[0](component) + _soln_vector(16 + 0,0) * _V1[0](component));
        (*_BNL)[i][j](3, component * 4) += - (*_BNL)[i][j](3, 1 + component * 4);

        // adding contributions corresponding to alpha 2 and 3 and beta 2 and 3
        (*_BNL)[i][j](3, 12+2) += - 1.0/32.0 * (1.0 + _2d_points[i](1)) * _thickness[i] * _V2[2](component) * (_soln_vector(2 + component*4,0) - _soln_vector(3 + component * 4,0));
        (*_BNL)[i][j](3, 16+2) += 1.0/32.0 * (1.0 + _2d_points[i](1)) * _thickness[i] * _V1[2](component) * (_soln_vector(2 + component*4,0) - _soln_vector(3 + component * 4,0));
        (*_BNL)[i][j](3, 12+3) += - 1.0/32.0 * (1.0 + _2d_points[i](1)) * _thickness[i] * _V2[3](component) * (_soln_vector(2 + component*4,0) - _soln_vector(3 + component * 4,0));
        (*_BNL)[i][j](3, 16+3) += 1.0/32.0 * (1.0 + _2d_points[i](1)) * _thickness[i] * _V1[3](component) * (_soln_vector(2 + component*4,0) - _soln_vector(3 + component * 4,0));

        // adding contributions corresponding to alpha 1 and 0 and beta 1 and 0
        (*_BNL)[i][j](3, 12+1) += - 1.0/32.0 * (1.0 - _2d_points[i](1)) * _thickness[i] * _V2[1](component) * (_soln_vector(1 + component*4,0) - _soln_vector(component * 4,0));
        (*_BNL)[i][j](3, 16+1) += 1.0/32.0 * (1.0 - _2d_points[i](1)) * _thickness[i] * _V1[1](component) * (_soln_vector(1 + component*4,0) - _soln_vector(component * 4,0));
        (*_BNL)[i][j](3, 12+0) += - 1.0/32.0 * (1.0 - _2d_points[i](1)) * _thickness[i] * _V2[0](component) * (_soln_vector(1 + component*4,0) - _soln_vector(component * 4,0));
        (*_BNL)[i][j](3, 16+0) += 1.0/32.0 * (1.0 - _2d_points[i](1)) * _thickness[i] * _V1[0](component) * (_soln_vector(1 + component*4,0) - _soln_vector(component * 4,0));

        // corresponding to strain(1,2)
        (*_BNL)[i][j](4,2 + component * 4) += 1.0/32.0 * (1.0 + _2d_points[i](0)) * _thickness[i] * (-_soln_vector(12 + 2,0) * _V2[2](component) + _soln_vector(16 + 2,0) * _V1[2](component) -_soln_vector(12 + 1,0) * _V2[1](component) + _soln_vector(16 + 1,0) * _V1[1](component));
        (*_BNL)[i][j](4, 1 + component * 4) += - (*_BNL)[i][j](3, 2 + component * 4);

        (*_BNL)[i][j](4,3 + component * 4) += 1.0/32.0 * (1.0 - _2d_points[i](0)) * _thickness[i] * (-_soln_vector(12 + 3,0) * _V2[3](component) + _soln_vector(16 + 3,0) * _V1[3](component) -_soln_vector(12 + 0,0) * _V2[0](component) + _soln_vector(16 + 0,0) * _V1[0](component));
        (*_BNL)[i][j](4, component * 4) += - (*_BNL)[i][j](3, 3 + component * 4);

        // adding contributions corresponding to alpha 2, 1 and beta 2 , 1
        (*_BNL)[i][j](4, 12+2) += - 1.0/32.0 * (1.0 + _2d_points[i](0)) * _thickness[i] * _V2[2](component) * (_soln_vector(2 + component*4,0) - _soln_vector(1 + component * 4,0));
        (*_BNL)[i][j](4, 16+2) += 1.0/32.0 * (1.0 + _2d_points[i](0)) * _thickness[i] * _V1[2](component) * (_soln_vector(2 + component*4,0) - _soln_vector(1 + component * 4,0));
        (*_BNL)[i][j](4, 12+1) += - 1.0/32.0 * (1.0 + _2d_points[i](0)) * _thickness[i] * _V2[1](component) * (_soln_vector(2 + component*4,0) - _soln_vector(1 + component * 4,0));
        (*_BNL)[i][j](4, 16+1) += 1.0/32.0 * (1.0 + _2d_points[i](0)) * _thickness[i] * _V1[1](component) * (_soln_vector(2 + component*4,0) - _soln_vector(1 + component * 4,0));

        // adding contributions corresponding to alpha 3, 0 and beta 3 , 0
        (*_BNL)[i][j](4, 12+3) += - 1.0/32.0 * (1.0 - _2d_points[i](0)) * _thickness[i] * _V2[3](component) * (_soln_vector(3 + component*4,0) - _soln_vector(component * 4,0));
        (*_BNL)[i][j](4, 16+3) += 1.0/32.0 * (1.0 - _2d_points[i](0)) * _thickness[i] * _V1[3](component) * (_soln_vector(3 + component*4,0) - _soln_vector(component * 4,0));
        (*_BNL)[i][j](4, 12+0) += - 1.0/32.0 * (1.0 - _2d_points[i](0)) * _thickness[i] * _V2[0](component) * (_soln_vector(3 + component*4,0) - _soln_vector(component * 4,0));
        (*_BNL)[i][j](4, 16+0) += 1.0/32.0 * (1.0 - _2d_points[i](0)) * _thickness[i] * _V1[0](component) * (_soln_vector(3 + component*4,0) - _soln_vector(component * 4,0));
      }

    /*  _BNL_new[i][j].zero();
      for (unsigned int k = 0; k < 4; ++k)
      {
        _BNL_new[i][j](0, k) = _dphidxi_map[k][i];
        _BNL_new[i][j](1, k + 4) = _dphidxi_map[k][i];
        _BNL_new[i][j](2, k + 8) = _dphidxi_map[k][i];
        _BNL_new[i][j](3, k) = _dphideta_map[k][i];
        _BNL_new[i][j](4, k + 4) = _dphideta_map[k][i];
        _BNL_new[i][j](5, k + 8) = _dphideta_map[k][i];

        for (unsigned int component = 0; component < 3; ++component)
        {
          _BNL_new[i][j](component, k+12) = - _dphidxi_map[k][i] * _thickness[i] * _t_points[j](0)/2.0 * _V2[k](component);
          _BNL_new[i][j](component+3, k+12) = - _dphideta_map[k][i] * _thickness[i] * _t_points[j](0)/2.0 * _V2[k](component);
          _BNL_new[i][j](component, k+16) = _dphidxi_map[k][i] * _thickness[i] * _t_points[j](0)/2.0 * _V1[k](component);
          _BNL_new[i][j](component+3, k+16) = _dphideta_map[k][i] * _thickness[i] * _t_points[j](0)/2.0 * _V1[k](component);
        }
      } */
    }
  }
}
