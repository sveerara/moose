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

#include "libmesh/quadrature.h"
#include "libmesh/utility.h"

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
    _total_strain(declareProperty<std::vector<RankTwoTensor>>)("total_strain"),
    _total_strain_old(declareMaterialPropertyOldByName<std::vector<RankTwoTensor>>)("total_strain"),
    _large_strain(getParam<bool>("large_strain")),
    _nonlinear_sys(_fe_problem.getNonlinearSystemBase()),
    _soln_disp_index_0(_ndisp),
    _soln_disp_index_1(_ndisp),
    _soln_rot_index_0(_ndisp),
    _soln_rot_index_1(_ndisp)
{
  // Checking for consistency between length of the provided displacements and rotations vector
  if (_ndisp != _nrot)
    mooseError("ComputeShellStrain: The number of variables supplied in 'displacements' "
               "and 'rotations' must match.");

  // fetch coupled variables and gradients (as stateful properties if necessary)
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    MooseVariable * disp_variable = getVar("displacements", i);
    _disp_num[i] = disp_variable->number();

    MooseVariable * rot_variable = getVar("rotations", i);
    _rot_num[i] = rot_variable->number();
  }

  for (unsigned int i = 0; i < _eigenstrain_names.size(); ++i)
  {
    _disp_eigenstrain[i] = &getMaterialProperty<RealVectorValue>("disp_" + _eigenstrain_names[i]);
    _rot_eigenstrain[i] = &getMaterialProperty<RealVectorValue>("rot_" + _eigenstrain_names[i]);
    _disp_eigenstrain_old[i] =
        &getMaterialPropertyOld<RealVectorValue>("disp_" + _eigenstrain_names[i]);
    _rot_eigenstrain_old[i] =
        &getMaterialPropertyOld<RealVectorValue>("rot_" + _eigenstrain_names[i]);
  }
}

void
ComputeShellStrain::initQpStatefulProperties()
{
  _t_qrule = new ArbitraryQuadrature(1, _order);
  _t_points = _t_qrule->get_points();

  // quadrature points in isoparametric space
  std::vector<Point> qp_points = _qrule->get_points(); // would be in 2D

   for (unsigned int i = 0; i < qp_points.size(); ++i)
     for (unsigned int j = 0; j < _t_points.size(); ++j)
       _3d_points[i][j](0) = qp_point[i](0);
       _3d_points[i][j](1) = qp_point[i](1);
       _3d_points[i][j](2) = _t_points[i](0);

  unsigned int dim = _current_elem->dim();
  unsigned int mesh_dim = _mesh.dimension();
  if ((dim != 2))
    mooseError("Shell element is implemented only for 2D Linear elements");

  // derivatives of shape functions (dphidxi, dphideta and dphidzeta) evaluated at quadrature points (in isoparametric space).
  FEBase * & fe = _subproblem.assembly(_tid).getFE(fe_type, dim);
  _dphidxi_map = fe->get_fe_map().get_dphidxi_map();
  _dphideta_map = fe->get_fe_map().get_dphideta_map();
  _phi_map = fe->get_fe_map().get_phi_map();

  // initialize node normals
}

void
ComputeShellStrain::computeProperties()
{
  // calculating derivatives of shape function is physical space (dphi/dx, dphi/dy, dphi/dz) at quadrature points
  std::vector<Node *> node;
  for (unsigned int i = 0; i < 4; ++i)
  node.push_back(_current_elem->get_node(i));

  for (unsigned int i = 0; i < _3d_points.size(); ++i)
  {
    for (unsigned int j = 0; j < _t_points.size(); ++j)
    {
      for (unsinged int component = 0; component < _mesh.dimension(); ++component)
      {
        _dxyz_dxi[i][j](component) = 0.0;
        for (unsigned int k = 0; k < node.size(); ++k)
        {
          _dxyz_dxi[i][j](component) += _dphidxi_map[k][i] * node[k](component) + _t_points[j](0) / 2.0 * _thickness[i] * _dphidxi_map[k][i] * _node_normal_old[k][component];
          _dxyz_deta[i][j](component) += _dphieta_map[k][i] * node[k](component) + _t_points[j](0) / 2.0 * _thickness[i] * _dphieta_map[k][i] * _node_normal_old[k][component];
          _dxyz_dzeta[i][j](component) += _thickness[i] * _phi_map[k][i] * _node_normal_old[k][component] / 2.0;
        }
      }
      RankTwoTensor Jac;
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

  // Fetch the incremental displacement (current - old) at the nodes
  const NumericVector<Number> & sol = *_nonlinear_sys.currentSolution();
  const NumericVector<Number> & sol_old = _nonlinear_sys.solutionOld();

  for (unsigned int j = 0; j < _nodes.size(); ++j)
  {
    for (unsigned int i = 0; i < _ndisp.size(); ++i)
    {
      _soln_disp_index[j][i] = node[j]->dof_number(_nonlinear_sys.number(), _disp_num[i], 0);
      _disp[j](i) = sol(_soln_disp_index[j][i]) - sol_old(_soln_disp_index[j][i]);
    }

    for (unsigned int i = 0; i < 2; ++i)
    {
      _soln_rot_index[j][i] = node[j]->dof_number(_nonlinear_sys.number(), _rot_num[i], 0);
      _rot[j](i) = sol(_soln_rot_index[j][i]) - sol_old(_soln_rot_index[j][i]);
    }
  }

  // compute nodal local axis
  RealGradient x2;
  x2(1) = 1;
  RealGradient x3;
  x3(2) = 1;

  for (unsigned int k = 0; k < nodes.size(); ++k)
  {
    _V1[k] = x2.cross(_node_normal_old[k]) / x2 * _node_normal_old[k];

    // If x2 is parallel to node normal, set V1 to x3
    if (MooseUtils::absoluteFuzzyEqual(_V1[k].norm(), 0.0, 1e-6))
      _V1[k] = x3;

    _V2[k] = _node_normal_old[k].cross(_V1[k]);
  }

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
  } */
}
