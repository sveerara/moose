/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/

#include "ComputeSmallStrain.h"
#include "Assembly.h"
#include "MooseMesh.h"
#include "ElasticityTensorTools.h"

// libmesh includes
#include "libmesh/quadrature.h"

template<>
InputParameters validParams<ComputeSmallStrain>()
{
  InputParameters params = validParams<ComputeStrainBase>();
  params.addClassDescription("Compute a small strain.");
  params.addParam<bool>("enhanced_strain", false, "Set to true to turn on incompatible mode elements.");
  return params;
}

ComputeSmallStrain::ComputeSmallStrain(const InputParameters & parameters) :
    ComputeStrainBase(parameters),
    _Jacobian_mult(getMaterialPropertyByName<RankFourTensor>(_base_name + "elasticity_tensor")),
    _enhanced_strain(getParam<bool>("enhanced_strain")),
    _disp_x(coupled("displacements", 0)),
    _disp_y(coupled("displacements", 1)),
    _disp_z(parameters.get<SubProblem *>("_subproblem")->mesh().dimension() == 3 ? coupled("displacements", 2) : 0),
    _dphi(_assembly.gradPhi()),
    _a(3, 1)
//    _dim(_current_elem->dim())
//    _mesh_dim(_mesh.dimension())
{
}

void
ComputeSmallStrain::computeProperties()
{
  if (_enhanced_strain)
    computeEnhancedStrain();
  Real volumetric_strain = 0.0;
  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    //strain = (grad_disp + grad_disp^T)/2
    RankTwoTensor grad_tensor((*_grad_disp[0])[_qp], (*_grad_disp[1])[_qp], (*_grad_disp[2])[_qp]);

    _total_strain[_qp] = ( grad_tensor + grad_tensor.transpose() )/2.0;

    if (_enhanced_strain)
      computeEnhancedStrainIncrement(_qp, _total_strain[_qp]);

    if (_volumetric_locking_correction)
      volumetric_strain +=  _total_strain[_qp].trace() * _JxW[_qp] * _coord[_qp];
  }

  if (_volumetric_locking_correction)
    volumetric_strain /= _current_elem_volume;

  for (_qp = 0; _qp < _qrule->n_points(); ++_qp)
  {
    if (_volumetric_locking_correction)
    {
      Real trace = _total_strain[_qp].trace();
      _total_strain[_qp](0,0) += (volumetric_strain - trace) / 3.0;
      _total_strain[_qp](1,1) += (volumetric_strain - trace) / 3.0;
      _total_strain[_qp](2,2) += (volumetric_strain - trace) / 3.0;
    }

    _mechanical_strain[_qp] = _total_strain[_qp];

    //Remove the Eigen strain
    for (auto es : _eigenstrains)
      _mechanical_strain[_qp] -= (*es)[_qp];
  }
}

void
ComputeSmallStrain::computeEnhancedStrain()
{
  // quadrature points in isoparametric space
    std::vector<Point> qp_points = _qrule->get_points();

    FEType fe_type(FIRST, LAGRANGE);
    unsigned int dim = _current_elem->dim();
    unsigned int mesh_dim = _mesh.dimension();
    if ((dim != 2) && (dim != 3))
      mooseError("Incompatible mode elements is implemented only for 2D and 3D Linear elements");

    // derivatives of shape functions (dphidxi, dphideta and dphidzeta) evaluated at quadrature points (in isoparametric space).
    std::vector<std::vector<Real> > dphidxi(dim, std::vector<Real>(_qrule->n_points(), 0.0));
    std::vector<std::vector<Real> > dphideta(dim, std::vector<Real>(_qrule->n_points(), 0.0));
    std::vector<std::vector<Real> > dphidzeta(dim, std::vector<Real>(_qrule->n_points(), 0.0));

    for (unsigned int qp_loop = 0; qp_loop < _qrule->n_points(); qp_loop++)
    {
      dphidxi[0][qp_loop] = -2.0 * qp_points[qp_loop](0); // phi_0 = 1-xi^2
      dphideta[1][qp_loop] = -2.0 * qp_points[qp_loop](1); // phi_1 = 1-eta^2
      if (dim == 3)
        dphidzeta[2][qp_loop] = -2.0 * qp_points[qp_loop](2); // phi_3 = 1-zeta^2
    }

    // calculating derivatives of shape function is physical space (dphi/dx, dphi/dy, dphi/dz) at quadrature points
  //  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
    FEBase * & fe = _subproblem.assembly(_tid).getFE(fe_type, dim);
    std::vector<Real> dxidx_map = fe->get_fe_map().get_dxidx();
    std::vector<Real> dxidy_map = fe->get_fe_map().get_dxidy();
    std::vector<Real> dxidz_map(dxidy_map.size(), 0.0);
    if (mesh_dim == 3)
      dxidz_map = fe->get_fe_map().get_dxidz();

    std::vector<Real> detadx_map = fe->get_fe_map().get_detadx();
    std::vector<Real> detady_map = fe->get_fe_map().get_detady();
    std::vector<Real> detadz_map(dxidy_map.size(), 0.0);
    if (mesh_dim == 3)
      detadz_map = fe->get_fe_map().get_detadz();

    std::vector<Real> dzetadx_map(dxidy_map.size(), 0.0);
    std::vector<Real> dzetady_map(dxidy_map.size(), 0.0);
    std::vector<Real> dzetadz_map(dxidy_map.size(), 0.0);
    if (dim ==3)
    {
      dzetadx_map = fe->get_fe_map().get_dzetadx();
      dzetady_map = fe->get_fe_map().get_dzetady();
      dzetadz_map = fe->get_fe_map().get_dzetadz();
    }

    std::vector<std::vector<Real> > dphidx(dim, std::vector<Real>(_qrule->n_points(), 0.0));
    std::vector<std::vector<Real> > dphidy(dim, std::vector<Real>(_qrule->n_points(), 0.0));
    std::vector<std::vector<Real> > dphidz(dim, std::vector<Real>(_qrule->n_points(), 0.0));

    for (unsigned i = 0; i < dim; i++)
    {
      for (unsigned qp_loop = 0; qp_loop < _qrule->n_points(); qp_loop++)
      {
        dphidx[i][qp_loop] = dphidxi[i][qp_loop] * dxidx_map[qp_loop] + dphideta[i][qp_loop] * detadx_map[qp_loop];
        dphidy[i][qp_loop] = dphidxi[i][qp_loop] * dxidy_map[qp_loop] + dphideta[i][qp_loop] * detady_map[qp_loop];

        if ((dim == 2) && (mesh_dim == 3))
          dphidz[i][qp_loop] = dphidxi[i][qp_loop] * dxidz_map[qp_loop] + dphideta[i][qp_loop] * detadz_map[qp_loop];

        if (dim == 3)
        {
          dphidx[i][qp_loop] += dphidzeta[i][qp_loop] * dzetadx_map[qp_loop];
          dphidy[i][qp_loop] += dphidzeta[i][qp_loop] * dzetady_map[qp_loop];
          dphidz[i][qp_loop] = dphidxi[i][qp_loop] * dxidz_map[qp_loop] + dphideta[i][qp_loop] * detadz_map[qp_loop] + dphidzeta[i][qp_loop] * dzetadz_map[qp_loop];
        }
      }
    }

    // obtaining average values of shape function derivatives over the element, i.e., volume integral of dphidx, dphidy dphidz
    std::vector<Real> average_dphidx(mesh_dim, 0.0);
    std::vector<Real> average_dphidy(mesh_dim, 0.0);
    std::vector<Real> average_dphidz(mesh_dim, 0.0);
//    VariablePhiGradient dphi;
    _dphi.resize(dim);
    Real volume = 0.0;

    for (unsigned int i = 0; i < dim; i++)
    {
      _dphi[i].resize(_qrule->n_points());
      volume = 0.0;
      for (unsigned int qp_loop = 0; qp_loop < _qrule->n_points(); qp_loop++)
      {
        average_dphidx[i] += dphidx[i][qp_loop] * _JxW[qp_loop] * _coord[qp_loop];
        average_dphidy[i] += dphidy[i][qp_loop] * _JxW[qp_loop] * _coord[qp_loop];
        if (mesh_dim == 3)
          average_dphidz[i] += dphidz[i][qp_loop] * _JxW[qp_loop] * _coord[qp_loop];

        volume += _JxW[qp_loop] * _coord[qp_loop];
      }
      average_dphidx[i] = average_dphidx[i]/volume;
      average_dphidy[i] = average_dphidy[i]/volume;
      if (mesh_dim == 3)
        average_dphidz[i] = average_dphidz[i]/volume;

      // Edit shape functions for incompatible mode dofs to ensure patch test passes

      for (unsigned int qp_loop = 0; qp_loop < _qrule->n_points(); qp_loop++)
      {
        dphidx[i][qp_loop] -= average_dphidx[i];
        dphidy[i][qp_loop] -= average_dphidy[i];
        if (mesh_dim == 3)
          dphidz[i][qp_loop] -= average_dphidz[i];

        _dphi[i][qp_loop](0) = dphidx[i][qp_loop];
        _dphi[i][qp_loop](1) = dphidy[i][qp_loop];
        if (mesh_dim == 3)
          _dphi[i][qp_loop](2) = dphidz[i][qp_loop];
        else if (mesh_dim == 2)
          _dphi[i][qp_loop](2) = 0.0;
      }
    }
    // Gradient of shape functions for the regular dofs
  //  const VariablePhiGradient test(_subproblem->assembly(_tid).gradPhi());
    std::vector<std::vector<RealGradient> > test(fe->get_dphi());

    // Form the matrices E and H
    // [K   E] [u] = [F]
    // [E^T H] [a] = [0]
    ColumnMajorMatrix E(mesh_dim * test.size(), mesh_dim * _dphi.size());
    ColumnMajorMatrix H(mesh_dim * _dphi.size(), mesh_dim * _dphi.size());
    ColumnMajorMatrix H_inv(mesh_dim * _dphi.size(), mesh_dim * _dphi.size());
    _a.reshape(mesh_dim * _dphi.size(), 1);

  //  ColumnMajorMatrix a(mesh_dim * _dphi.size(), 1);
    ColumnMajorMatrix K(mesh_dim * test.size(), mesh_dim * test.size());

    for (unsigned qp_loop = 0; qp_loop < _qrule->n_points(); qp_loop++)
    {
      for (unsigned int j = 0; j < _dphi.size(); j++)
      {
        // for matrix E
        for (unsigned int i = 0; i < test.size(); i++)
        {
          for (unsigned int ii = 0; ii < mesh_dim; ii++)
          {
            for (unsigned int jj = 0; jj < mesh_dim; jj++)
              E(i * mesh_dim + ii, j * mesh_dim + jj) += ElasticityTensorTools::elasticJacobian(_Jacobian_mult[qp_loop], ii, jj,
                                                                                test[i][qp_loop], _dphi[j][qp_loop]) * _JxW[qp_loop] * _coord[qp_loop];
          }
        }

        // for matrix H
        for (unsigned int i = 0; i < _dphi.size(); i++)
        {
          for (unsigned int ii = 0; ii < mesh_dim; ii++)
          {
            for (unsigned int jj = 0; jj < mesh_dim; jj++)
              H(i * mesh_dim + ii, j * mesh_dim + jj) += ElasticityTensorTools::elasticJacobian(_Jacobian_mult[qp_loop], ii, jj,
                                                                              _dphi[i][qp_loop], _dphi[j][qp_loop]) * _JxW[qp_loop] * _coord[qp_loop];
          }
        }
      }

      // for matrix K
      for (unsigned int j = 0; j < test.size(); j++)
      {
        for (unsigned int i = 0; i < test.size(); i++)
        {
          for (unsigned int ii = 0; ii < mesh_dim; ii++)
          {
            for (unsigned int jj = 0; jj < mesh_dim; jj++)
              K(i * mesh_dim + ii, j * mesh_dim + jj) += ElasticityTensorTools::elasticJacobian(_Jacobian_mult[qp_loop], ii, jj,
                                                                              test[i][qp_loop], test[j][qp_loop]) * _JxW[qp_loop] * _coord[qp_loop];
          }
        }
      }
    }

    H.inverse(H_inv);
    H_inv *= -1.0;

    // Obtain displacements at the nodes of the element
    ColumnMajorMatrix u(mesh_dim * test.size(), 1);
    NonlinearSystem & sys = _fe_problem.getNonlinearSystem();
    const NumericVector<Real> & current_solution = *sys.currentSolution();
    MooseVariable & var_x = sys.getVariable(_tid, _disp_x);
    MooseVariable & var_y = sys.getVariable(_tid, _disp_y);
    std::vector<dof_id_type> dof_indices_x = var_x.dofIndices();
    std::vector<dof_id_type> dof_indices_y = var_y.dofIndices();
    //std::vector<dof_id_type> dof_indices_x;
    //std::vector<dof_id_type> dof_indices_y;
    //var_x.getDofIndices(_current_elem, dof_indices_x);
    //var_y.getDofIndices(_current_elem, dof_indices_y);

    std::vector<dof_id_type> dof_indices_z;

    if (mesh_dim == 3)
    {
      MooseVariable & var_z = sys.getVariable(_tid, _disp_z);
      dof_indices_z = var_z.dofIndices();
    }

    std::vector<unsigned int> test_dof_map(4);
    for (unsigned int i = 0; i < dof_indices_x.size(); i++)
    {
      u(i * mesh_dim) = current_solution(dof_indices_x[i]);
      u(i * mesh_dim + 1) = current_solution(dof_indices_y[i]);
      if (mesh_dim == 3)
        u(i * mesh_dim + 2) = current_solution(dof_indices_z[i]);
    }

    // Calculate incompatible dofs
    ColumnMajorMatrix ET = E.transpose();
    _a = H_inv * ET * u;
}

void
ComputeSmallStrain::computeEnhancedStrainIncrement(unsigned int qp, RankTwoTensor & strain_increment)
{
     // Edit the strain increment tensor to include strains from incompatible modes
    unsigned int _dim = _current_elem->dim();
    unsigned int _mesh_dim = _mesh.dimension();
    for (unsigned int i = 0; i < _dim; i++)
    {
        strain_increment(0,0) += _dphi[i][qp](0) * _a(i * _mesh_dim);
        strain_increment(1,1) += _dphi[i][qp](1) * _a(i * _mesh_dim + 1);
        strain_increment(0,1) += 0.5 * (_dphi[i][qp](1) * _a(i * _mesh_dim) + _dphi[i][qp](0) * _a(i * _mesh_dim + 1));
        strain_increment(1,0) = strain_increment(0,1);
        if (_mesh_dim == 3)
        {
          strain_increment(2,2) += _dphi[i][qp](2) * _a(i * _mesh_dim + 2);
          strain_increment(0,2) += 0.5 * (_dphi[i][qp](2) * _a(i * _mesh_dim) + _dphi[i][qp](0) * _a(i * _mesh_dim + 2));
          strain_increment(1,2) += 0.5 * (_dphi[i][qp](2) * _a(i * _mesh_dim + 1) + _dphi[i][qp](1) * _a(i * _mesh_dim + 2));
          strain_increment(2,1) = strain_increment(1,2);
          strain_increment(2,0) = strain_increment(0,2);
        }
    }
  //K.print();
//    ColumnMajorMatrix G = K * u;

    // Calulate edited stiffness matrix
  //  _K = K - E * H_inv * ET;

}
