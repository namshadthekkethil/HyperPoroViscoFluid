
#include "HyperElasticModel.h"

DenseMatrix<double> HyperElasticModel::S, HyperElasticModel::P;
DenseVector<double> HyperElasticModel::Resid;
DenseMatrix<double> HyperElasticModel::Jacob;

DMAT HyperElasticModel::Dmat;
DHyper HyperElasticModel::D_hyper;

double HyperElasticModel::mu;

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;

HyperElasticModel::HyperElasticModel() {}
HyperElasticModel::~HyperElasticModel() {}

void HyperElasticModel::initialise_lde(EquationSystems &es,
                                       LargeDeformationElasticity &lde)
{
  lde.define_systems();
  es.parameters.set<Real>("nonlinear solver absolute residual tolerance") =
      InputParam::nonlinear_abs_tol;
  es.parameters.set<Real>("nonlinear solver relative residual tolerance") =
      InputParam::nonlinear_rel_tol;
  es.parameters.set<unsigned int>("nonlinear solver maximum iterations") =
      InputParam::nonlinear_max_its;

  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");
  system.nonlinear_solver->residual_object = &lde;
  system.nonlinear_solver->jacobian_object = &lde;

  LargeDeformationElasticity::pre_itr = &compute_pre_itr;
  LargeDeformationElasticity::param_qp_resid = &compute_param_resid_qp;
  LargeDeformationElasticity::param_qp_jacob = &compute_param_jacob_qp;
  LargeDeformationElasticity::param_dofi_resid = &compute_param_resid_dofi;
  LargeDeformationElasticity::param_dofi_jacob = &compute_param_jacob_dofi;
  LargeDeformationElasticity::param_dofi_dofj_jacob =
      &compute_param_jacob_dofi_dofj;
  LargeDeformationElasticity::boundary_force =
      &BoundaryCond::apply_pressure_force;
  LargeDeformationElasticity::torsion_resid =
      &BoundaryCond::apply_torsion_resid;
  LargeDeformationElasticity::torsion_jacob =
      &BoundaryCond::apply_torsion_jacob;

  Stabilisation::compute_mu = &compute_mu;
}

void HyperElasticModel::define_systems(EquationSystems &es)
{

  System &f0_system = es.add_system<System>("fiber direction");
  System &s0_system = es.add_system<System>("sheet direction");

  System &system_dis_tot = es.add_system<System>("displacement_tot");
  system_dis_tot.add_variable("u_tot", FIRST, LAGRANGE);
  system_dis_tot.add_variable("v_tot", FIRST, LAGRANGE);
#if (MESH_DIMENSION == 3)
  system_dis_tot.add_variable("w_tot", FIRST, LAGRANGE);
#endif
  system_dis_tot.add_variable("pressure_tot", FIRST, LAGRANGE);

  ExplicitSystem &J_system = es.add_system<ExplicitSystem>("JTotSystem");
  J_system.add_variable("JTotVar", CONSTANT, MONOMIAL);

  ExplicitSystem &pext_system = es.add_system<ExplicitSystem>("pExtSystem");
  pext_system.add_variable("pExtVar", CONSTANT, MONOMIAL);

  for (unsigned int d = 0; d < MESH_DIMENSION; ++d)
  {
    ostringstream os;
    os << "f0_" << d;
    f0_system.add_variable(os.str(), CONSTANT, MONOMIAL);
  }

  for (unsigned int d = 0; d < MESH_DIMENSION; ++d)
  {
    ostringstream os;
    os << "s0_" << d;
    s0_system.add_variable(os.str(), CONSTANT, MONOMIAL);
  }

  Incompress::define_systems(es);
  Stabilisation::define_systems(es);

  BoundaryCond::read_boundaries(es);
  BoundaryCond::define_systems(es);

  if (InputParam::viscoelastic == 1)
    ViscoElastic::define_systems(es);
}

void HyperElasticModel::init_hyperelastic_model(EquationSystems &es, int rank)
{
  if (InputParam::strain_model == 1)
  {
    HOModel::aHO = InputParam::aHO;
    HOModel::bHO = InputParam::bHO;
    HOModel::afHO = InputParam::afHO;
    HOModel::bfHO = InputParam::bfHO;
    HOModel::asHO = InputParam::asHO;
    HOModel::bsHO = InputParam::bsHO;
    HOModel::afsHO = InputParam::afsHO;
    HOModel::bfsHO = InputParam::bfsHO;
  }

  else if (InputParam::strain_model == 2)
  {
    NeoHook::G = InputParam::G;
  }

  else if (InputParam::strain_model == 3)
  {
    CellFibre::G = InputParam::G;
    CellFibre::alpha_fib = InputParam::alpha_fib;
  }

  if (InputParam::viscoelastic == 1)
  {
    ViscoElastic::tau_visela = InputParam::tau_visela;
    ViscoElastic::beta_visela = InputParam::beta_visela;
  }

  cout << "G=" << CellFibre::G << " alpha_fib=" << CellFibre::alpha_fib << endl;

  Incompress::beta_s = InputParam::beta_s;
  Incompress::dt = InputParam::dt;
  Incompress::time_itr = InputParam::time_itr;

  Stabilisation::rho_s = InputParam::rho_s;
  Stabilisation::alpha_stab = InputParam::alpha_stab;
  Stabilisation::dt = InputParam::dt;
  Stabilisation::time_itr = 0;

  initialize_material_axes(es, rank);

  GeomPar::init_geoPar();
  Stabilisation::init_stab(es);

  BoundaryCond::init_boundary_cond(es);
}

void HyperElasticModel::update_hyperelastic_model(EquationSystems &es)
{

  Incompress::solve_incomp(es);

  InputParam::time_itr++;

  Incompress::dt = InputParam::dt;
  Incompress::time_itr = InputParam::time_itr;

  Stabilisation::dt = InputParam::dt;
  Stabilisation::time_itr = InputParam::time_itr;

  Stabilisation::solve_prime(es);

  if (InputParam::viscoelastic == 1)
  {
    ViscoElastic::update_QS(es);
  }
}

void HyperElasticModel::initialize_material_axes(EquationSystems &es, int rank)
{
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  DenseMatrix<double> f0_data;
  f0_data.resize(mesh.n_elem(),dim);

  cout<<"mesh_n_elem="<<mesh.n_elem()<<endl;

  if (rank == 0)
  {
    ifstream f0_stream(InputParam::fibre_file_name);
    for (unsigned int e = 0; e < mesh.n_elem(); ++e)
    {
      int e_idx;
      f0_stream >> e_idx;
      for (unsigned int d = 0; d < dim; ++d)
        f0_stream >> f0_data(e,d);
    }
    f0_stream.close();
  }

  DenseMatrix<double> s0_data;
  s0_data.resize(mesh.n_elem(),dim);

  if (rank == 0)
  {
    ifstream s0_stream(InputParam::sheet_file_name);
    for (unsigned int e = 0; e < mesh.n_elem(); ++e)
    {
      int e_idx;
      s0_stream >> e_idx;

      for (unsigned int d = 0; d < dim; ++d)
        s0_stream >> s0_data(e,d);
    }
    s0_stream.close();
  }

  for(unsigned int e = 0; e < mesh.n_elem(); ++e)
  {
    for (unsigned int d = 0; d < dim; ++d)
    {
      MPI_Bcast(&f0_data(e,d), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&s0_data(e,d), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }
  }


  System &f0_system = es.get_system<System>("fiber direction");
  NumericVector<double> &f0_vec = *f0_system.solution;
  System &s0_system = es.get_system<System>("sheet direction");
  NumericVector<double> &s0_vec = *s0_system.solution;

  int f0_system_num = f0_system.number();
  int s0_system_num = s0_system.number();

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  MeshBase::const_element_iterator el_st = mesh.active_elements_begin();
  Elem *elem_st = *el_st;
  const int elem_id_st = elem_st->id();

  for (; el != end_el; ++el)
  {
    Elem *elem = *el;
    const int elem_id = elem->id();

    for (unsigned int d = 0; d < dim; ++d)
    {
      const int dof_index = elem->dof_number(f0_system_num, d, 0);
      f0_vec.set(dof_index, f0_data(elem_id - elem_id_st,d));
    }

    for (unsigned int d = 0; d < dim; ++d)
    {
      const int dof_index = elem->dof_number(s0_system_num, d, 0);
      s0_vec.set(dof_index, s0_data(elem_id - elem_id_st,d));
    }
  }
  f0_vec.close();
  s0_vec.close();
  f0_vec.localize(*f0_system.current_local_solution);
  s0_vec.localize(*s0_system.current_local_solution);

  return;
}

void HyperElasticModel::compute_PK2_hyper(EquationSystems &es,
                                          const Elem *elem)
{
  S.resize(MESH_DIMENSION, MESH_DIMENSION);
  if (InputParam::strain_model == 1)
  {
    HOModel::compute_HO_PK2(es, elem);
    S.add(1.0, HOModel::SIso);
    S.add(1.0, HOModel::Sf);
    S.add(1.0, HOModel::Ss);
    S.add(1.0, HOModel::Sfs);
  }

  else if (InputParam::strain_model == 2)
  {
    NeoHook::compute_NH_PK2(es, elem);
    S.add(1.0, NeoHook::SIso);
  }

  else if (InputParam::strain_model == 3)
  {
    CellFibre::compute_CF_PK2(es, elem);
    S.add(1.0, CellFibre::SIso);
    S.add(1.0, CellFibre::Sf);
    S.add(1.0, CellFibre::Ss);
  }
}

void HyperElasticModel::compute_PK2(EquationSystems &es, const Elem *elem)
{
  compute_PK2_hyper(es, elem);
  if (InputParam::viscoelastic == 1)
  {
    ViscoElastic::compute_Q(es, elem);
    S.resize(MESH_DIMENSION, MESH_DIMENSION);
    S.add(1.0, ViscoElastic::Q);
  }

  if (InputParam::strain_model == 1)
  {
    S.add(1.0, HOModel::Sp);
  }

  else if (InputParam::strain_model == 2)
  {
    S.add(1.0, NeoHook::Sp);
  }

  else if (InputParam::strain_model == 3)
  {
    S.add(1.0, CellFibre::Sp);
  }

  if (InputParam::active_tension == 1)
  {
    if (InputParam::ttime > InputParam::t_end_diastole)
    {
      ActiveContra::compute_PK2_active();
      S.add(1.0, ActiveContra::Sa);
    }
  }

  P.resize(MESH_DIMENSION, MESH_DIMENSION);
  P.add(1.0, GeomPar::F);
  P.right_multiply(S);
}

void HyperElasticModel::compute_D()
{
  if (InputParam::strain_model == 1)
  {
    HOModel::compute_HO_D_iso();
    HOModel::compute_HO_D_aniso();
    HOModel::compute_D_p();
  }

  else if (InputParam::strain_model == 2)
  {
    NeoHook::compute_NH_D_iso();
    NeoHook::compute_D_p();
  }

  else if (InputParam::strain_model == 3)
  {
    CellFibre::compute_CF_D_iso();
    CellFibre::compute_CF_D_aniso();
    CellFibre::compute_D_p();
  }

  if (InputParam::active_tension == 1)
    if (InputParam::ttime > InputParam::t_end_diastole)
      ActiveContra::compute_D_active();

  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++)
        {
          if (InputParam::strain_model == 1)
          {
            Dmat.D[i][j][k][l] =
                HOModel::Dmat.Diso[i][j][k][l] + HOModel::Dmat.Df[i][j][k][l] +
                HOModel::Dmat.Ds[i][j][k][l] + HOModel::Dmat.Dfs[i][j][k][l] +
                HOModel::Dmat.Dp[i][j][k][l];
          }
          else if (InputParam::strain_model == 2)
          {
            Dmat.D[i][j][k][l] =
                NeoHook::Dmat.Diso[i][j][k][l] + NeoHook::Dmat.Dp[i][j][k][l];
          }

          else if (InputParam::strain_model == 3)
          {
            Dmat.D[i][j][k][l] = CellFibre::Dmat.Diso[i][j][k][l] +
                                 CellFibre::Dmat.Df[i][j][k][l] +
                                 CellFibre::Dmat.Ds[i][j][k][l] +
                                 CellFibre::Dmat.Dp[i][j][k][l];
          }

          if (InputParam::viscoelastic == 1)
          {
            Dmat.D[i][j][k][l] = 0.0;
            if (InputParam::strain_model == 1)
            {
              Dmat.D[i][j][k][l] +=
                  exp(-InputParam::dt / (2.0 * ViscoElastic::tau_visela)) *
                      (HOModel::Dmat.Diso[i][j][k][l] +
                       HOModel::Dmat.Df[i][j][k][l] +
                       HOModel::Dmat.Ds[i][j][k][l] +
                       HOModel::Dmat.Dfs[i][j][k][l]) +
                  HOModel::Dmat.Dp[i][j][k][l];
            }
            else if (InputParam::strain_model == 2)
            {
              Dmat.D[i][j][k][l] +=
                  exp(-InputParam::dt / (2.0 * ViscoElastic::tau_visela)) *
                      (NeoHook::Dmat.Diso[i][j][k][l]) +
                  NeoHook::Dmat.Dp[i][j][k][l];
            }

            else if (InputParam::strain_model == 3)
            {
              Dmat.D[i][j][k][l] +=
                  (exp(-InputParam::dt / (2.0 * ViscoElastic::tau_visela))) *
                      (CellFibre::Dmat.Diso[i][j][k][l] +
                       CellFibre::Dmat.Df[i][j][k][l] +
                       CellFibre::Dmat.Ds[i][j][k][l]) +
                  CellFibre::Dmat.Dp[i][j][k][l];
            }
          }

          if (InputParam::active_tension == 1)
            if (InputParam::ttime > InputParam::t_end_diastole)
              Dmat.D[i][j][k][l] += ActiveContra::Dact.Da[i][j][k][l];

          D_hyper.FD[i][j][k][l] = 0.0;
        }

  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++)
          for (unsigned int m = 0; m < MESH_DIMENSION; m++)
          {
            D_hyper.FD[i][j][k][l] += GeomPar::F(i, m) * Dmat.D[m][j][k][l];
          }
}

void HyperElasticModel::compute_mu(EquationSystems &es, const Elem *elem)
{
  if (InputParam::strain_model == 1)
  {
    HOModel::compute_mu(es, elem);
    Stabilisation::mu = HOModel::mu;
  }
  else if (InputParam::strain_model == 2)
  {
    NeoHook::compute_mu(es, elem);
    Stabilisation::mu = NeoHook::mu;
  }
  else if (InputParam::strain_model == 3)
  {
    CellFibre::compute_mu(es, elem);
    Stabilisation::mu = CellFibre::mu;
  }

  if (InputParam::active_tension == 1)
    if (InputParam::ttime > InputParam::t_end_diastole)
    {
      ActiveContra::compute_mu();
      Stabilisation::mu += ActiveContra::mu;
    }
}

void HyperElasticModel::compute_pre_itr(EquationSystems &es)
{
  Unsteady::update_velocity(es);
  Unsteady::update_acceleration(es);
  Stabilisation::compute_tau(es);
}

void HyperElasticModel::compute_param_resid_qp(
    EquationSystems &es, const Elem *elem, unsigned int qp,
    const std::vector<std::vector<Real>> &phi,
    const std::vector<std::vector<RealGradient>> &dphi)
{
  GeomPar::compute_geoPar(es, elem, qp, phi, dphi);

  compute_PK2(es, elem);

  Incompress::compute_incomp_res(es, elem, qp, phi);

  Stabilisation::compute_dofi_cur(es, elem);
  Stabilisation::compute_prime();
  Stabilisation::compute_u_prime(es, elem, qp, phi);
}

void HyperElasticModel::compute_param_jacob_qp(
    EquationSystems &es, const Elem *elem, unsigned int qp,
    const std::vector<std::vector<Real>> &phi,
    const std::vector<std::vector<RealGradient>> &dphi)
{
  GeomPar::compute_geoPar(es, elem, qp, phi, dphi);
  TensorDer::compute_dJdF();
  TensorDer::compute_dJFinvTdF();

  compute_PK2(es, elem);
  compute_D();

  Stabilisation::compute_dofi_cur(es, elem);
  Stabilisation::compute_dofij_cur(es, elem);
  Stabilisation::compute_prime();
  Stabilisation::compute_u_prime(es, elem, 0, phi);
}

void HyperElasticModel::compute_param_resid_dofi(
    const std::vector<std::vector<Real>> &phi,
    const std::vector<std::vector<RealGradient>> &dphi, unsigned int qp,
    unsigned int dof_i)
{
  TensorDer::compute_gradNA(dphi, qp, dof_i);
  Stabilisation::compute_stab();

  Resid.resize(MESH_DIMENSION + 1);
  Resid(MESH_DIMENSION) = Incompress::incomp_res * phi[dof_i][qp];
  Resid(MESH_DIMENSION) +=
      Stabilisation::dofi_cur(dof_i) * Stabilisation::res_stab;

  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
  {

    Resid(i) += phi[dof_i][qp] * (InputParam::rho_s * GeomPar::acc_vec(i));
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
    {
      Resid(i) += (P(i, j) * dphi[dof_i][qp](j));
    }
  }
}

void HyperElasticModel::compute_param_jacob_dofi(
    const std::vector<std::vector<Real>> &phi,
    const std::vector<std::vector<RealGradient>> &dphi, unsigned int qp,
    unsigned int dof_i)
{
  TensorDer::compute_gradNA(dphi, qp, dof_i);
  TensorDer::compute_gradNAxyz(dphi, qp, dof_i);
  Stabilisation::compute_stab();
}

void HyperElasticModel::compute_param_jacob_dofi_dofj(
    const std::vector<std::vector<Real>> &phi,
    const std::vector<std::vector<RealGradient>> &dphi, unsigned int qp,
    unsigned int dof_i, unsigned int dof_j)
{
  TensorDer::compute_gradNB(dphi, qp, dof_j);
  TensorDer::compute_dJFinvTdu();
  Incompress::compute_incomp_res_der(phi, qp, dof_j);
  Stabilisation::compute_stab_der(phi, qp, dof_j);

  Jacob.resize(MESH_DIMENSION + 1, MESH_DIMENSION + 1);
  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
  {
    Jacob(MESH_DIMENSION, i) += Incompress::dincompdu(i) * phi[dof_i][qp];
    Jacob(MESH_DIMENSION, i) +=
        Stabilisation::dofij_cur(dof_i, dof_j) * Stabilisation::dstabdu(i);
  }

  Jacob(MESH_DIMENSION, MESH_DIMENSION) +=
      Incompress::dincompdp * phi[dof_i][qp];
  Jacob(MESH_DIMENSION, MESH_DIMENSION) +=
      Stabilisation::dofij_cur(dof_i, dof_j) * Stabilisation::dstabdp;

  Jacob(0, MESH_DIMENSION) +=
      phi[dof_j][qp] * GeomPar::detF *
      (MatVecOper::contractMat(GeomPar::FInvTra, TensorDer::gradNAx));
  Jacob(1, MESH_DIMENSION) +=
      phi[dof_j][qp] * GeomPar::detF *
      (MatVecOper::contractMat(GeomPar::FInvTra, TensorDer::gradNAy));
#if (MESH_DIMENSION == 3)
  Jacob(2, MESH_DIMENSION) +=
      phi[dof_j][qp] * GeomPar::detF *
      (MatVecOper::contractMat(GeomPar::FInvTra, TensorDer::gradNAz));
#endif

  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    Jacob(i, i) +=
        (InputParam::rho_s * (1.0 / (InputParam::dt * InputParam::dt))) *
        phi[dof_i][qp] * phi[dof_j][qp];

  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int m = 0; m < MESH_DIMENSION; m++)
        Jacob(i, i) += (dphi[dof_j][qp](m) * S(m, j) * dphi[dof_i][qp](j));

  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++)
        {
          Jacob(i, k) += (0.5 * D_hyper.FD[i][j][k][l] * dphi[dof_j][qp](l) *
                          dphi[dof_i][qp](j));

          Jacob(i, l) += (0.5 * D_hyper.FD[i][j][k][l] * dphi[dof_j][qp](k) *
                          dphi[dof_i][qp](j));

          for (unsigned int n = 0; n < MESH_DIMENSION; n++)
            Jacob(i, n) += (0.5 * D_hyper.FD[i][j][k][l] *
                            (dphi[dof_j][qp](k) * GeomPar::grad_u(n, l) +
                             dphi[dof_j][qp](l) * GeomPar::grad_u(n, k)) *
                            dphi[dof_i][qp](j));
        }
}

void HyperElasticModel::update_total_velocity_displacement(
    EquationSystems &es)
{
  NonlinearImplicitSystem &system_dis =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");
  System &system_dis_tot = es.get_system<System>("displacement_tot");
  System &system_dis_fine = es.get_system<System>("fineDisSystem");

  NumericVector<double> &dis_solution = *(system_dis.solution);
  dis_solution.close();
  dis_solution.localize(*system_dis.current_local_solution);

  NumericVector<double> &dis_tot_solution = *(system_dis_tot.solution);
  dis_tot_solution.close();
  dis_tot_solution.localize(*system_dis_tot.current_local_solution);

  NumericVector<double> &dis_fine_solution = *(system_dis_fine.solution);
  dis_fine_solution.close();
  dis_fine_solution.localize(*system_dis_fine.current_local_solution);

  libMesh::MeshBase &mesh = es.get_mesh();
  libMesh::MeshBase::const_node_iterator node_it = mesh.local_nodes_begin();
  libMesh::MeshBase::const_node_iterator node_end = mesh.local_nodes_end();

  for (; node_it != node_end; ++node_it)
  {
    const libMesh::Node *nd = *node_it;

    const unsigned int dof_num_dis_x =
        nd->dof_number(system_dis.number(), 0, 0);
    const unsigned int dof_num_dis_y =
        nd->dof_number(system_dis.number(), 1, 0);
#if (MESH_DIMENSION == 3)
    const unsigned int dof_num_dis_z =
        nd->dof_number(system_dis.number(), 2, 0);
#endif

    const unsigned int dof_num_dis_fine_x =
        nd->dof_number(system_dis_fine.number(), 0, 0);
    const unsigned int dof_num_dis_fine_y =
        nd->dof_number(system_dis_fine.number(), 1, 0);
#if (MESH_DIMENSION == 3)
    const unsigned int dof_num_dis_fine_z =
        nd->dof_number(system_dis_fine.number(), 2, 0);
#endif

    double ux_cur = dis_solution(dof_num_dis_x);
    double uxp_cur = dis_fine_solution(dof_num_dis_fine_x);
    double uy_cur = dis_solution(dof_num_dis_y);
    double uyp_cur = dis_fine_solution(dof_num_dis_fine_y);
#if (MESH_DIMENSION == 3)
    double uz_cur = dis_solution(dof_num_dis_z);
    double uzp_cur = dis_fine_solution(dof_num_dis_fine_z);
#endif

    const unsigned int dof_num_dis_tot_x =
        nd->dof_number(system_dis_tot.number(), 0, 0);
    const unsigned int dof_num_dis_tot_y =
        nd->dof_number(system_dis_tot.number(), 1, 0);
#if (MESH_DIMENSION == 3)
    const unsigned int dof_num_dis_tot_z =
        nd->dof_number(system_dis_tot.number(), 2, 0);
#endif

    system_dis_tot.solution->set(dof_num_dis_tot_x, ux_cur + uxp_cur);
    system_dis_tot.solution->set(dof_num_dis_tot_y, uy_cur + uyp_cur);
#if (MESH_DIMENSION == 3)
    system_dis_tot.solution->set(dof_num_dis_tot_z, uz_cur + uzp_cur);
#endif
  }

  system_dis_tot.solution->close();
  system_dis_tot.solution->localize(*system_dis_tot.current_local_solution);
}

void HyperElasticModel::compute_Jtot(EquationSystems &es)
{
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem &system =
      es.get_system<LinearImplicitSystem>("displacement_tot");

  unsigned int displacement_vars[dim];
  displacement_vars[0] = system.variable_number("u_tot");
  displacement_vars[1] = system.variable_number("v_tot");
  if (dim == 3)
    displacement_vars[2] = system.variable_number("w_tot");
  const unsigned int u_var = system.variable_number("u_tot");

  const DofMap &dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, CONSTANT);
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem &J_system = es.get_system<ExplicitSystem>("JTotSystem");
  const DofMap &J_dof_map = J_system.get_dof_map();
  unsigned int J_var;
  J_var = J_system.variable_number("JTotVar");

  // Storage for the stress dof indices on each element
  std::vector<std::vector<dof_id_type>> dof_indices_var(system.n_vars());
  std::vector<dof_id_type> J_dof_indices_var;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  // int count_Jn = 0;

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    for (unsigned int var = 0; var < dim; var++)
      dof_map.dof_indices(elem, dof_indices_var[var], displacement_vars[var]);

    const unsigned int n_var_dofs = dof_indices_var[0].size();

    fe->reinit(elem);

    double J_cur = 0.0;
    DenseMatrix<Number> F(dim, dim);
    DenseMatrix<Number> grad_u(dim, dim);
    for (unsigned int var_i = 0; var_i < dim; var_i++)
      for (unsigned int var_j = 0; var_j < dim; var_j++)
        for (unsigned int j = 0; j < n_var_dofs; j++)
          grad_u(var_i, var_j) +=
              dphi[j][0](var_j) *
              system.current_solution(dof_indices_var[var_i][j]);

    F = grad_u;
    for (unsigned int var = 0; var < dim; var++)
      F(var, var) += 1.;
    J_cur = MatVecOper::detMat(F);

    const int dof_index = elem->dof_number(J_system.number(), 0, 0);

    J_system.solution->set(dof_index, J_cur);
  }
  // cout<<count_Jn<<endl;

  // Should call close and update when we set vector entries directly
  J_system.solution->close();
  J_system.solution->localize(*J_system.current_local_solution);
}

void HyperElasticModel::compute_pext(EquationSystems &es)
{
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

  const unsigned int u_var = system.variable_number("u");

  const DofMap &dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, CONSTANT);
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  ExplicitSystem &pext_system = es.get_system<ExplicitSystem>("pExtSystem");
  const DofMap &pext_dof_map = pext_system.get_dof_map();

  std::vector<dof_id_type> pext_dof_indices_var;


  // To store the stress tensor on each element
  DenseMatrix<Number> elem_avg_stress_tensor(dim, dim);

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    fe->reinit(elem);

    // clear the stress tensor
    elem_avg_stress_tensor.resize(dim, dim);

    GeomPar::compute_geoPar(es, elem, 0, phi, dphi);
    compute_PK2(es,elem);

    elem_avg_stress_tensor.add(1.0, S);

    elem_avg_stress_tensor.scale(1. / GeomPar::detF);
    elem_avg_stress_tensor.left_multiply(GeomPar::F);
    elem_avg_stress_tensor.right_multiply_transpose(GeomPar::F);

    double pext_cur = 0.0;

    #if(MESH_DIMENSION == 2)
    pext_cur = (-1.0 / 2.0)*(elem_avg_stress_tensor(0,0)+elem_avg_stress_tensor(1,1));
    #elif(MESH_DIMENSION == 3)
    pext_cur = (-1.0 / 3.0)*(elem_avg_stress_tensor(0,0)+elem_avg_stress_tensor(1,1)+elem_avg_stress_tensor(2,2));
    #endif

    pext_dof_map.dof_indices(elem, pext_dof_indices_var, 0);
    pext_system.solution->set(pext_dof_indices_var[0], pext_cur);
  }

  // Should call close and update when we set vector entries directly
  pext_system.solution->close();
  pext_system.update();
  pext_system.solution->localize(*pext_system.current_local_solution);

}
