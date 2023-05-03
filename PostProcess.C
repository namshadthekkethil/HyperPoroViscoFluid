#include "PostProcess.h"

double PostProcess::force_p, PostProcess::force_ppore, PostProcess::force_visc;
double PostProcess::force_p_total, PostProcess::force_ppore_total, PostProcess::m_net_total, PostProcess::force_visc_total;

DenseMatrix<double> PostProcess::FP, PostProcess::FPPORE, PostProcess::NETM,
    PostProcess::FPRATE, PostProcess::FPPORERATE, PostProcess::FTOTALRATE;

vector<double> PostProcess::fp_time, PostProcess::fppore_time, PostProcess::ftotal_time;

string PostProcess::file_surface_force, PostProcess::file_FbyS, PostProcess::file_FbyS_log;

using namespace libMesh;
using namespace std;

PostProcess::PostProcess() {}

PostProcess::~PostProcess() {}

void PostProcess::compute_vA(EquationSystems &es, double &u_A, double &v_A,
                             double &p_A, double &I4f_A, double &ppore_A)
{
  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

  System &system_vel = es.get_system<System>("velocity");

  System &system_pre = es.get_system<System>("pMonoSystem");

  System &system_I4f = es.get_system<System>("TSystem");

#if (POROUS == 1)
  System &system_ppore = es.get_system<System>("pPoreSystem");
#endif

  Point p_qp(0.0, 5.0);

  u_A = 0.0;
  v_A = 0.0;

  u_A = system.point_value(1, p_qp);
  v_A = system_vel.point_value(1, p_qp);
  // p_A = system.point_value(MESH_DIMENSION,p_qp);
  p_A = system_pre.point_value(0, p_qp);

  I4f_A = system_I4f.point_value(2, p_qp);

#if (POROUS == 1)
  ppore_A = system_ppore.point_value(0, p_qp);
#endif

  cout << "vA=" << v_A << ", pA=" << p_A << endl;
}

void PostProcess::compute_mixture_volume(EquationSystems &es,
                                         double &mixture_volume)
{
  const MeshBase &mesh_cur = es.get_mesh();

  MeshBase::const_element_iterator el = mesh_cur.active_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh_cur.active_elements_end();

  mixture_volume = 0.0;

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    mixture_volume += elem->volume();
  }

  cout << "mixture_volume=" << mixture_volume << endl;
}

void PostProcess::compute_skeleton_volume(EquationSystems &es,
                                          EquationSystems &es_cur,
                                          double &J_total, double &m_total)
{
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  const MeshBase &mesh_cur = es_cur.get_mesh();

  System &m_system = es.get_system<System>("mMonoSystem");
  NumericVector<double> &m_vec = *m_system.solution;
  int m_system_num = m_system.number();

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  MeshBase::const_element_iterator el_cur =
      mesh_cur.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el_cur =
      mesh_cur.active_local_elements_end();

  J_total = 0.0, m_total = 0.0;

  double elem_cur_vol = 0.0, elem_vol = 0.0, m_sum = 0.0;

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    const Elem *elem_cur = *el_cur;

    elem_vol += elem->volume();
    elem_cur_vol += elem_cur->volume();

    const int dof_index = elem->dof_number(m_system.number(), 0, 0);

    m_sum += m_vec(dof_index) * elem->volume();

    ++el_cur;
  }

  J_total = elem_cur_vol / elem_vol;
  m_total = m_sum / elem_vol;
}

void PostProcess::compute_L2_norm(EquationSystems &es, double &l2_dis_tot,
                                  double &l2_vel_tot, double &l2_p_tot)
{
  MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  System &system_l2_dis = es.get_system<System>("l2_dis_system");
  const DofMap &dof_map_l2_dis = system_l2_dis.get_dof_map();

  System &system_l2_vel = es.get_system<System>("l2_vel_system");
  const DofMap &dof_map_l2_vel = system_l2_vel.get_dof_map();

  System &system_l2_p = es.get_system<System>("l2_p_system");
  const DofMap &dof_map_l2_p = system_l2_p.get_dof_map();

  System &system = es.get_system<System>("velocity");
  const unsigned int u_var = system.variable_number("vel_x");
  const DofMap &dof_map = system.get_dof_map();

  System &system_dis = es.get_system<System>("NonlinearElasticity");
  const unsigned int dx_var = system_dis.variable_number("u");
  const DofMap &dof_map_dis = system_dis.get_dof_map();

  std::vector<dof_id_type> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_var(dim + 1);
  std::vector<std::vector<dof_id_type>> dof_indices_dis(dim + 1);
  std::vector<std::vector<dof_id_type>> dof_indices_l2_dis(dim);
  std::vector<std::vector<dof_id_type>> dof_indices_l2_vel(dim);
  std::vector<dof_id_type> dof_indices_l2_p;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    dof_map.dof_indices(elem, dof_indices);
    dof_map_l2_p.dof_indices(elem, dof_indices_l2_p);
    for (unsigned int var = 0; var < dim + 1; var++)
    {
      dof_map.dof_indices(elem, dof_indices_var[var], var);
      dof_map_dis.dof_indices(elem, dof_indices_dis[var], var);
    }

    for (unsigned int var = 0; var < dim; var++)
    {
      dof_map_l2_dis.dof_indices(elem, dof_indices_l2_dis[var], var);
      dof_map_l2_vel.dof_indices(elem, dof_indices_l2_vel[var], var);
    }

    const unsigned int n_var_dofs = dof_indices_var[0].size();

    for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++)
    {

      double dx_sol = system_dis.current_solution(dof_indices_dis[0][dof_i]);
      double dy_sol = system_dis.current_solution(dof_indices_dis[1][dof_i]);
      double p_sol = system_dis.current_solution(dof_indices_dis[2][dof_i]);

      double u_sol = system.current_solution(dof_indices_var[0][dof_i]);
      double v_sol = system.current_solution(dof_indices_var[1][dof_i]);

      Point pj;
      pj = elem->point(dof_i);

      double x_cur = pj(0);
      double y_cur = pj(1);

      double omega_freq = 0.5 * M_PI * sqrt((2.0 * 5666666.67) / 1100.0);

      double dx_ana = 0.01 * sin(omega_freq * 0.01) *
                      (-sin(0.5 * M_PI * x_cur) * cos(0.5 * M_PI * y_cur));
      double dy_ana = 0.01 * sin(omega_freq * 0.01) *
                      (cos(0.5 * M_PI * x_cur) * sin(0.5 * M_PI * y_cur));

      double u_ana = 0.01 * omega_freq *
                     (-sin(0.5 * M_PI * x_cur) * cos(0.5 * M_PI * y_cur));
      double v_ana = 0.01 * omega_freq *
                     (cos(0.5 * M_PI * x_cur) * sin(0.5 * M_PI * y_cur));
      double p_ana = 0.0;

      double l1_dis_x = (dx_sol - dx_ana);
      double l1_dis_y = (dy_sol - dy_ana);

      double l1_vel_x = (u_sol - u_ana);
      double l1_vel_y = (v_sol - v_ana);

      double l1_p = (p_sol - p_ana);

      system_l2_dis.solution->set(dof_indices_l2_dis[0][dof_i], l1_dis_x);
      system_l2_dis.solution->set(dof_indices_l2_dis[1][dof_i], l1_dis_y);
      system_l2_vel.solution->set(dof_indices_l2_vel[0][dof_i], l1_vel_x);
      system_l2_vel.solution->set(dof_indices_l2_vel[1][dof_i], l1_vel_y);
      system_l2_p.solution->set(dof_indices_l2_p[dof_i], l1_p);
    }
  }
  system_l2_dis.solution->close();
  system_l2_vel.solution->close();
  system_l2_p.solution->close();
  system_l2_dis.solution->localize(*system_l2_dis.current_local_solution);
  system_l2_vel.solution->localize(*system_l2_vel.current_local_solution);
  system_l2_p.solution->localize(*system_l2_p.current_local_solution);

  l2_dis_tot = system_l2_dis.solution->l2_norm();
  l2_vel_tot = system_l2_vel.solution->l2_norm();
  l2_p_tot = system_l2_p.solution->l2_norm();
}

void PostProcess::compute_surface_ppore_m(EquationSystems &es, double &m_0,
                                          double &ppore_0, double &lambda_0,
                                          double &m_1, double &ppore_1,
                                          double &lambda_1)
{
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  System &m_system = es.get_system<System>("NonlinearElasticity");
  const DofMap &dof_map_m = m_system.get_dof_map();
  std::vector<dof_id_type> dof_indices_m;

  m_system.solution->localize(*m_system.current_local_solution);
  NumericVector<double> &m_data = *(m_system.current_local_solution);
  m_data.close();

  const unsigned int u_var = m_system.variable_number("m");

  FEType fe_type = dof_map_m.variable_type(u_var);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_type));
  QGauss qface(dim - 1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule(&qface);

  System &ppore_system = es.get_system<System>("pPoreSystem");
  const DofMap &dof_map_ppore = ppore_system.get_dof_map();
  std::vector<dof_id_type> dof_indices_ppore;

  ppore_system.solution->localize(*ppore_system.current_local_solution);
  NumericVector<double> &ppore_data = *(ppore_system.current_local_solution);
  ppore_data.close();

  System &lambda_system = es.get_system<System>("NonlinearElasticity");
  const DofMap &dof_map_lambda = lambda_system.get_dof_map();
  std::vector<dof_id_type> dof_indices_lambda;

  lambda_system.solution->localize(*lambda_system.current_local_solution);
  NumericVector<double> &lambda_data = *(lambda_system.current_local_solution);
  lambda_data.close();

  vector<double> m_vec;
  vector<double> ppore_vec;
  vector<double> lambda_vec;
  m_data.localize(m_vec);
  ppore_data.localize(ppore_vec);
  lambda_data.localize(lambda_vec);

  m_0 = 0.0, ppore_0 = 0.0, lambda_0 = 0.0;
  m_1 = 0.0, ppore_1 = 0.0, lambda_1 = 0.0;
  double area_0 = 0.0, area_1 = 0.0;

  MeshBase::const_element_iterator el = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    dof_map_m.dof_indices(elem, dof_indices_m, dim + 1);
    dof_map_ppore.dof_indices(elem, dof_indices_ppore, 0);
    dof_map_lambda.dof_indices(elem, dof_indices_lambda, dim);

    for (unsigned int side = 0; side < elem->n_sides(); side++)
      if (elem->neighbor_ptr(side) == libmesh_nullptr)
      {
        const std::vector<std::vector<Real>> &phi_face = fe_face->get_phi();
        const std::vector<Real> &JxW_face = fe_face->get_JxW();

        fe_face->reinit(elem, side);

        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
        {
          vector<boundary_id_type> bc_id_vec;
          mesh.boundary_info->boundary_ids(elem, side, bc_id_vec);
          short int bc_id =
              bc_id_vec[0]; // mesh.boundary_info->boundary_id(elem, side);

          if (bc_id == 0)
          {
            double m_cur = 0.0, ppore_cur = 0.0, lambda_cur = 0.0;
            for (std::size_t i = 0; i < phi_face.size(); i++)
            {
              m_cur += m_vec[dof_indices_m[i]] * phi_face[i][qp];
              ppore_cur += ppore_vec[dof_indices_ppore[i]] * phi_face[i][qp];
              lambda_cur += lambda_vec[dof_indices_lambda[i]] * phi_face[i][qp];
            }
            m_0 += m_cur * JxW_face[qp];
            ppore_0 += ppore_cur * JxW_face[qp];
            lambda_0 += lambda_cur * JxW_face[qp];
            area_0 += JxW_face[qp];
          }

          if (bc_id == 5)
          {
            double m_cur = 0.0, ppore_cur = 0.0, lambda_cur = 0.0;
            for (std::size_t i = 0; i < phi_face.size(); i++)
            {
              m_cur += m_vec[dof_indices_m[i]] * phi_face[i][qp];
              ppore_cur += ppore_vec[dof_indices_ppore[i]] * phi_face[i][qp];
              lambda_cur += lambda_vec[dof_indices_lambda[i]] * phi_face[i][qp];
            }
            m_1 += m_cur * JxW_face[qp];
            ppore_1 += ppore_cur * JxW_face[qp];
            lambda_1 += lambda_cur * JxW_face[qp];
            area_1 += JxW_face[qp];
          }
        }
      }
  }

  cout << "m_0=" << m_0 << " ppore_0=" << ppore_0 << " lambda_0=" << lambda_0
       << " m_1=" << m_1 << " ppore_1=" << ppore_1 << " lambda_1=" << lambda_1
       << endl;

  m_0 /= area_0;
  ppore_0 /= area_0;
  lambda_0 /= area_0;

  m_1 /= area_1;
  ppore_1 /= area_1;
  lambda_1 /= area_1;
}

void PostProcess::compute_surface_force(EquationSystems &es,
                                        EquationSystems &es_cur)
{
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  ExplicitSystem &stress_system = es.get_system<ExplicitSystem>("StressSystem");
  NumericVector<double> &stress_vec = *stress_system.solution;
  int stress_system_num = stress_system.number();

  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

  const unsigned int u_var = system.variable_number("u");

  const DofMap &dof_map = system.get_dof_map();

  ExplicitSystem &p_system = es.get_system<ExplicitSystem>("pMonoSystem");
  NumericVector<double> &p_vec = *p_system.solution;
  int p_system_num = p_system.number();

  // #if(POROUS == 1)
  System &system_ppore = es.get_system<System>("pPoreSystem");
  const DofMap &dof_map_ppore = system_ppore.get_dof_map();
  std::vector<dof_id_type> dof_indices_ppore;
  // #endif

  LinearImplicitSystem &darcy_system =
      es.get_system<LinearImplicitSystem>("darcy_velocity");
  const DofMap &dof_map_darcy = darcy_system.get_dof_map();
  std::vector<std::vector<dof_id_type>> dof_indices_darcy(dim);

  system.solution->localize(*system.current_local_solution);
  NumericVector<double> &p_data = *(system.current_local_solution);
  p_data.close();

  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_type));
  QGauss qface(dim - 1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule(&qface);

  /************************************ for current
   * ****************************/
  System &system_cur = es_cur.get_system<System>("NodeIdSystem");

  const unsigned int u_var_cur = system_cur.variable_number("NodeIdVar");

  const DofMap &dof_map_cur = system_cur.get_dof_map();

  FEType fe_type_cur = dof_map_cur.variable_type(u_var_cur);
  UniquePtr<FEBase> fe_cur(FEBase::build(dim, fe_type_cur));
  QGauss qrule_cur(dim, fe_type_cur.default_quadrature_order());
  fe_cur->attach_quadrature_rule(&qrule_cur);

  UniquePtr<FEBase> fe_face_cur(FEBase::build(dim, fe_type_cur));
  QGauss qface_cur(dim - 1, fe_type_cur.default_quadrature_order());
  fe_face_cur->attach_quadrature_rule(&qface_cur);
  /****************************************************************************/

  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  std::vector<dof_id_type> dof_indices_lambda;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  double f_x = 0.0, f_y = 0.0, f_z = 0.0;

  force_p = 0.0;
  force_ppore = 0.0;
  force_visc = 0.0;

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    dof_map.dof_indices(elem, dof_indices_lambda, MESH_DIMENSION);

    // #if(POROUS == 1)
    dof_map_ppore.dof_indices(elem, dof_indices_ppore, 0);
    // #endif

    for (unsigned int var = 0; var < dim; var++)
    {
      dof_map_darcy.dof_indices(elem, dof_indices_darcy[var], var);
    }

    double sigma_00 = 0.0, sigma_01 = 0.0, sigma_02 = 0.0;
    double sigma_11 = 0.0, sigma_12 = 0.0;
    double sigma_22 = 0.0;
    const int dof_index_00 = elem->dof_number(stress_system_num, 0, 0);
    const int dof_index_01 = elem->dof_number(stress_system_num, 1, 0);
#if (MESH_DIMENSION == 2)
    const int dof_index_11 = elem->dof_number(stress_system_num, 2, 0);
#elif (MESH_DIMENSION == 3)
    const int dof_index_02 = elem->dof_number(stress_system_num, 2, 0);
    const int dof_index_11 = elem->dof_number(stress_system_num, 3, 0);
    const int dof_index_12 = elem->dof_number(stress_system_num, 4, 0);
    const int dof_index_22 = elem->dof_number(stress_system_num, 5, 0);
#endif

    const int dof_index_p = elem->dof_number(p_system_num, 0, 0);
    double p_cur = -p_vec(dof_index_p);

    sigma_00 = stress_vec(dof_index_00); //-p_cur;
    sigma_01 = stress_vec(dof_index_01);
    sigma_11 = stress_vec(dof_index_11); //-p_cur;
#if (MESH_DIMENSION == 3)
    sigma_02 = stress_vec(dof_index_02);
    sigma_12 = stress_vec(dof_index_12);
    sigma_22 = stress_vec(dof_index_22); //-p_cur;
#endif

    for (unsigned int side = 0; side < elem->n_sides(); side++)
      if (elem->neighbor_ptr(side) == libmesh_nullptr)
      {

        const std::vector<std::vector<Real>> &phi_face = fe_face->get_phi();
        const std::vector<std::vector<RealGradient>> &dphi_face =
            fe_face->get_dphi();
        const std::vector<Real> &JxW_face = fe_face->get_JxW();

        const std::vector<Point> &normal_face = fe_face->get_normals();
        const std::vector<Point> &Xref = fe_face->get_xyz();
        const std::vector<vector<Point>> &tangent_face =
            fe_face->get_tangents();

        const std::vector<std::vector<Real>> &phi_face_cur =
            fe_face_cur->get_phi();
        const std::vector<Real> &JxW_face_cur = fe_face_cur->get_JxW();
        const std::vector<Point> &normal_face_cur = fe_face_cur->get_normals();

        fe_face->reinit(elem, side);
        fe_face_cur->reinit(elem, side);

        for (unsigned int qp = 0; qp < qface.n_points(); qp++)
        {
          vector<boundary_id_type> bc_id_vec;
          mesh.boundary_info->boundary_ids(elem, side, bc_id_vec);
          short int bc_id =
              bc_id_vec[0]; // mesh.boundary_info->boundary_id(elem, side);

          if (bc_id == 5)
          {

            double lambda_cur = 0.0, ppore_cur = 0.0;
            for (std::size_t i = 0; i < phi_face.size(); i++)
            {
              lambda_cur += phi_face[i][qp] *
                            system.current_solution(dof_indices_lambda[i]);

              // #if(POROUS == 1)
              ppore_cur += phi_face[i][qp] *
                           system_ppore.current_solution(dof_indices_ppore[i]);
              // #endif
            }
            double pressure_cur = 0.0;
            // #if(POROUS == 0)
            // pressure_cur = lambda_cur;
            // #elif(POROUS == 1)
            pressure_cur = -ppore_cur;
            // #endif

            f_x += normal_face_cur[qp].operator()(0) * JxW_face_cur[qp] *
                   (sigma_00 + pressure_cur);
            f_x +=
                normal_face_cur[qp].operator()(1) * JxW_face_cur[qp] * sigma_01;
            f_x +=
                normal_face_cur[qp].operator()(2) * JxW_face_cur[qp] * sigma_02;

            f_y +=
                normal_face_cur[qp].operator()(0) * JxW_face_cur[qp] * sigma_01;
            f_y += normal_face_cur[qp].operator()(1) * JxW_face_cur[qp] *
                   (sigma_11 + pressure_cur);
            f_y +=
                normal_face_cur[qp].operator()(2) * JxW_face_cur[qp] * sigma_12;

            f_z +=
                normal_face_cur[qp].operator()(0) * JxW_face_cur[qp] * sigma_02;
            f_z +=
                normal_face_cur[qp].operator()(1) * JxW_face_cur[qp] * sigma_12;
            f_z += normal_face_cur[qp].operator()(2) * JxW_face_cur[qp] *
                   (sigma_22 + pressure_cur);

            force_p += normal_face_cur[qp].operator()(0) * JxW_face_cur[qp] *
                       (lambda_cur);

            force_ppore += normal_face_cur[qp].operator()(0) *
                           JxW_face_cur[qp] * (-ppore_cur);

            if (InputParam::brinkman == 1)
            {

              DenseVector<double> gradWt(dim);
              for (std::size_t i = 0; i < phi_face.size(); i++)
              {
                for (unsigned int var_i = 0; var_i < dim; var_i++)
                {
                  for (unsigned int var_j = 0; var_j < dim; var_j++)
                  {
                    gradWt(var_i) += dphi_face[i][qp](var_i) *
                                     darcy_system.current_solution(
                                         dof_indices_darcy[var_j][i]) *
                                     tangent_face[qp][0].operator()(var_j);
                  }
                }
              }

              double dWdt = 0.0;
              for (unsigned int var_j = 0; var_j < dim; var_j++)
              {
                dWdt += gradWt(var_j) * normal_face[qp].operator()(var_j);
              }

              force_visc += InputParam::viscocity * dWdt *
                            tangent_face[qp][0].operator()(0) * JxW_face_cur[qp];
            }
          }
        }
      }
  }

  MPI_Allreduce(&force_p, &force_p_total, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&force_ppore, &force_ppore_total, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  MPI_Allreduce(&force_visc, &force_visc_total, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
}

void PostProcess::init_postprocess(int rank)
{
  FP.resize(InputParam::vtaubya_size, InputParam::vabyd_size);
  FPPORE.resize(InputParam::vtaubya_size, InputParam::vabyd_size);
  NETM.resize(InputParam::vtaubya_size, InputParam::vabyd_size);
  FPRATE.resize(InputParam::vtaubya_size, InputParam::vabyd_size);
  FPPORERATE.resize(InputParam::vtaubya_size, InputParam::vabyd_size);
  FTOTALRATE.resize(InputParam::vtaubya_size, InputParam::vabyd_size);

  string surface_force_name = "results/surface_force_";
  string FbyS_name = "results/FbyS_";
  string FbyS_log_name = "results/FbyS_log_";

  string fileend = ".dat";

  char numstr_vtaubya[21];
  sprintf(numstr_vtaubya, "%.0e", InputParam::Vtaubya(0));

  string file_underscore = "_";

  char numstr_vabyd[21];
  sprintf(numstr_vabyd, "%.0e", InputParam::VabyD(0));

  file_surface_force = surface_force_name + numstr_vtaubya + file_underscore + numstr_vabyd + fileend;
  file_FbyS = FbyS_name + numstr_vtaubya + file_underscore + numstr_vabyd + fileend;
  file_FbyS_log = FbyS_log_name + numstr_vtaubya + file_underscore + numstr_vabyd + fileend;

  if (rank == 0)
  {
    ofstream file1;
    file1.open(file_surface_force, ios::out);
    file1.close();

    ofstream file_F;
    file_F.open(file_FbyS, ios::out);
    file_F.close();

    ofstream file_F_log;
    file_F_log.open(file_FbyS_log, ios::out);
    file_F_log.close();
  }
}

void PostProcess::update_postprocess(EquationSystems &es,
                                     EquationSystems &es_cur, int rank)
{
  compute_surface_force(es, es_cur);
  compute_net_m(es);
  if (rank == 0)
  {
    ofstream file1;
    file1.open(file_surface_force, ios::app);
    file1 << InputParam::time_itr * InputParam::dt << " " << force_p_total << " "
          << force_ppore_total << " " << m_net_total << " " << InputParam::torsion_t <<" "<<force_visc_total<< endl;
    file1.close();
  }

  fp_time.push_back(force_p_total / InputParam::a_bead);
  fppore_time.push_back(force_ppore_total / InputParam::a_bead);
  ftotal_time.push_back((0.3 * force_p_total + 0.7 * force_ppore_total) / InputParam::a_bead);
}

void PostProcess::compute_net_m(EquationSystems &es)
{
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem &mmono_system = es.get_system<ExplicitSystem>("mMonoSystem");
  const DofMap &mmono_dof_map = mmono_system.get_dof_map();
  unsigned int mmono_var;
  mmono_var = mmono_system.variable_number("mMonoVar");

  ExplicitSystem &p_system = es.get_system<ExplicitSystem>("pMonoSystem");
  NumericVector<double> &p_vec = *p_system.solution;
  int p_system_num = p_system.number();

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  double m_net = 0.0, p_net = 0.0;

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    const int dof_index = elem->dof_number(mmono_system.number(), 0, 0);
    const int dof_index_p = elem->dof_number(p_system.number(), 0, 0);

    m_net += fabs(mmono_system.current_solution(dof_index));
    p_net += fabs(p_system.current_solution(dof_index_p));
  }

  m_net /= p_net;

  MPI_Allreduce(&m_net, &m_net_total, 1, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
}

void PostProcess::update_force_rate(int i, int j)
{
  double dt_nondim = InputParam::dt * (InputParam::V_bead / InputParam::a_bead);

  FPRATE(i, j) = 0.0;
  FPPORERATE(i, j) = 0.0;
  FTOTALRATE(i, j) = 0.0;
  for (int n = InputParam::n_solves + 1; n < fp_time.size(); n++)
  {
    FPRATE(i, j) += fabs(fp_time[n] - fp_time[n - 1]) / dt_nondim;
    FPPORERATE(i, j) += fabs(fppore_time[n] - fppore_time[n - 1]) / dt_nondim;
    FTOTALRATE(i, j) += fabs(ftotal_time[n] - ftotal_time[n - 1]) / dt_nondim;
  }

  FPRATE(i, j) /= fp_time[InputParam::n_solves];
  FPPORERATE(i, j) /= fppore_time[InputParam::n_solves];
  FTOTALRATE(i, j) /= ftotal_time[InputParam::n_solves];
}

void PostProcess::write_final(int rank)
{

  for (int jj = 0; jj < InputParam::vabyd_size; jj++)
  {
    for (int ii = 0; ii < InputParam::vtaubya_size; ii++)
    {
      if (rank == 0)
      {
        ofstream file_F;
        ofstream file_F_log;
        file_F.open(file_FbyS, ios::app);
        file_F_log.open(file_FbyS_log, ios::app);
        file_F << InputParam::Vtaubya(ii) << " " << InputParam::VabyD(jj) << " "
               << FP(ii, jj) << " " << FPPORE(ii, jj) << " " << NETM(ii, jj) << " "
               << FPRATE(ii, jj) << " " << FPPORERATE(ii, jj) << " " << FTOTALRATE(ii, jj) << endl;

        file_F_log << log(InputParam::Vtaubya(ii)) << " "
                   << log(InputParam::VabyD(jj)) << " " << FP(ii, jj) << " "
                   << FPPORE(ii, jj) << " " << NETM(ii, jj) << " "
                   << FPRATE(ii, jj) << " " << FPPORERATE(ii, jj) << " " << FTOTALRATE(ii, jj) << endl;
        file_F.close();
        file_F_log.close();
      }
    }
  }
}
