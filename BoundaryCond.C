
#include "BoundaryCond.h"

DenseVector<Number> BoundaryCond::Re_u, BoundaryCond::Re_v;
#if (MESH_DIMENSION == 3)
DenseVector<Number> BoundaryCond::Re_w;
#endif
unsigned int BoundaryCond::qf_points, BoundaryCond::num_nodes;
size_t BoundaryCond::phif_size;

DenseVector<double> BoundaryCond::force_x_t, BoundaryCond::force_y_t;

#if (MESH_DIMENSION == 3)
DenseVector<double> BoundaryCond::force_z_t;
#endif

DenseVector<double> BoundaryCond::pressure_t;

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;

BoundaryCond::BoundaryCond() {}
BoundaryCond::~BoundaryCond() {}

void BoundaryCond::read_boundaries(EquationSystems &es) {
  GetPot infile1("input.in");
  std::string input_file_name = infile1("input_file_name", "input_LV.in");
  // GetPot infile("input.in");
  GetPot infile(input_file_name);

  if (InputParam::force_x_size > 0) {
    force_x_t.resize(InputParam::force_x_size);
  }

  if (InputParam::force_y_size > 0) {
    force_y_t.resize(InputParam::force_y_size);
  }

#if (MESH_DIMENSION == 3)
  if (InputParam::force_z_size > 0) {
    force_z_t.resize(InputParam::force_z_size);
  }
#endif

  if (InputParam::pressure_size > 0) {
    pressure_t.resize(InputParam::pressure_size);
  }

  clamp_boundaries(es);
}

void BoundaryCond::clamp_boundaries(EquationSystems &es) {
  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

  const unsigned int u_var = system.variable_number("u");
  const unsigned int v_var = system.variable_number("v");
#if (MESH_DIMENSION == 3)
  const unsigned int w_var = system.variable_number("w");
#endif

  if (InputParam::clamp_x_size > 0) {
    std::set<boundary_id_type> clamped_boundaries;
    for (boundary_id_type nbc = 0; nbc < InputParam::clamp_x_size; nbc++) {
      clamped_boundaries.insert(InputParam::clamp_x_bcs(nbc));
    }
    std::vector<unsigned int> uvw;
    uvw.push_back(u_var);
    ZeroFunction<Number> zero;

    system.get_dof_map().add_dirichlet_boundary(
        DirichletBoundary(clamped_boundaries, uvw, zero, LOCAL_VARIABLE_ORDER));
  }

  if (InputParam::clamp_y_size > 0) {
    std::set<boundary_id_type> clamped_boundaries;
    for (boundary_id_type nbc = 0; nbc < InputParam::clamp_y_size; nbc++) {
      clamped_boundaries.insert(InputParam::clamp_y_bcs(nbc));
    }
    std::vector<unsigned int> uvw;
    uvw.push_back(v_var);
    ZeroFunction<Number> zero;

    system.get_dof_map().add_dirichlet_boundary(
        DirichletBoundary(clamped_boundaries, uvw, zero, LOCAL_VARIABLE_ORDER));
  }

#if (MESH_DIMENSION == 3)
  if (InputParam::clamp_z_size > 0) {
    std::set<boundary_id_type> clamped_boundaries;
    for (boundary_id_type nbc = 0; nbc < InputParam::clamp_z_size; nbc++) {
      clamped_boundaries.insert(InputParam::clamp_z_bcs(nbc));
    }
    std::vector<unsigned int> uvw;
    uvw.push_back(w_var);
    ZeroFunction<Number> zero;

    system.get_dof_map().add_dirichlet_boundary(
        DirichletBoundary(clamped_boundaries, uvw, zero, LOCAL_VARIABLE_ORDER));
  }
#endif
}

void BoundaryCond::define_systems(EquationSystems &es) {
  System &system_bcid = es.add_system<System>("bcidSystem");
  system_bcid.add_variable("bcidVar", FIRST, LAGRANGE);
}

void BoundaryCond::init_boundary_cond(EquationSystems &es) {

  const MeshBase &mesh = es.get_mesh();
  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

  const unsigned int u_var = system.variable_number("u");

  const DofMap &dof_map = system.get_dof_map();
  std::vector<std::vector<dof_id_type>> dof_indices_var(MESH_DIMENSION);
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();

  const Elem *elem = *el;
  for (unsigned int var = 0; var < MESH_DIMENSION; var++) {
    dof_map.dof_indices(elem, dof_indices_var[var], var);
  }

  num_nodes = dof_indices_var[0].size();

  Re_u.resize(num_nodes);
  Re_v.resize(num_nodes);
#if (MESH_DIMENSION == 3)
  Re_w.resize(num_nodes);
#endif

  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe(FEBase::build(MESH_DIMENSION, fe_type));
  QGauss qrule(MESH_DIMENSION, CONSTANT);
  fe->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(MESH_DIMENSION, fe_type));
  QGauss qface(MESH_DIMENSION - 1, CONSTANT);
  fe_face->attach_quadrature_rule(&qface);

  const std::vector<std::vector<Real>> &phi_face = fe_face->get_phi();

  fe->reinit(elem);
  unsigned int sidef = 0;
  fe_face->reinit(elem, sidef);

  qf_points = qface.n_points();
  phif_size = phi_face.size();

  compute_bcid(es);
}

void BoundaryCond::compute_pressure() {

  if (InputParam::ttime <= InputParam::t_load)
    for (boundary_id_type nbc = 0; nbc < InputParam::pressure_size; nbc++) {
      pressure_t(nbc) = InputParam::pressure_values(nbc) *
                        (InputParam::ttime / InputParam::t_load);
    }

  else if (InputParam::ttime > InputParam::t_load &&
           InputParam::ttime <= InputParam::t_end_diastole)
    for (boundary_id_type nbc = 0; nbc < InputParam::pressure_size; nbc++) {
      pressure_t(nbc) = InputParam::pressure_values(nbc);
    }

  else {
    double time_syst = InputParam::ttime - InputParam::t_end_diastole;

    if (time_syst < 0.15) {
      for (boundary_id_type nbc = 0; nbc < InputParam::pressure_size; nbc++) {
        pressure_t(nbc) = InputParam::pressure_values(nbc) +
                          InputParam::pressure_sys *
                              (1.0 - exp(-(time_syst * time_syst) / 0.004));
      }
    } else if (time_syst < 0.3) {
      for (boundary_id_type nbc = 0; nbc < InputParam::pressure_size; nbc++) {
        pressure_t(nbc) =
            (InputParam::pressure_values(nbc) + InputParam::pressure_sys) *
            (1.0 - exp(-((0.3 - time_syst) * (0.3 - time_syst)) / 0.004));
      }
    } else {
      for (boundary_id_type nbc = 0; nbc < InputParam::pressure_size; nbc++) {
        pressure_t(nbc) = 0.0;
      }
    }
  }
}

void BoundaryCond::compute_torsion() {
  if (InputParam::torsion_type == 4) {
    if(InputParam::V_bead * InputParam::ttime < InputParam::bead_disp*InputParam::a_bead)
      InputParam::torsion_t = InputParam::V_bead * InputParam::ttime;
    else
      InputParam::torsion_t = InputParam::bead_disp*InputParam::a_bead;
  }
}

void BoundaryCond::apply_pressure(
    EquationSystems &es, const Elem *elem, unsigned int side,
    const std::vector<std::vector<Real>> &phi_face,
    const std::vector<Real> &JxW_face, const std::vector<Point> &normal_face,
    const std::vector<std::vector<Real>> &phi_face_cur,
    const std::vector<Real> &JxW_face_cur,
    const std::vector<Point> &normal_face_cur,
    DenseSubVector<double> Re_var[]) {

  const MeshBase &mesh = es.get_mesh();

  for (unsigned int qp = 0; qp < qf_points; qp++) {
    vector<boundary_id_type> bc_id_vec;
    mesh.boundary_info->boundary_ids(elem, side, bc_id_vec);
    short int bc_id = bc_id_vec[0];

    for (boundary_id_type nbc = 0; nbc < InputParam::pressure_size; nbc++) {
      if (bc_id == InputParam::pressure_bcs(nbc)) {
        for (std::size_t i = 0; i < phi_face.size(); i++) {
          if (InputParam::pressure_current == 0) {
            Re_var[0](i) += -JxW_face[qp] * pressure_t(nbc) * phi_face[i][qp] *
                            normal_face[qp].operator()(0);
            Re_var[1](i) += -JxW_face[qp] * pressure_t(nbc) * phi_face[i][qp] *
                            normal_face[qp].operator()(1);
#if (MESH_DIMENSION == 3)
            Re_var[2](i) += -JxW_face[qp] * pressure_t(nbc) * phi_face[i][qp] *
                            normal_face[qp].operator()(2);
#endif
          } else if (InputParam::pressure_current == 1) {
            Re_var[0](i) += -JxW_face_cur[qp] * pressure_t(nbc) *
                            phi_face_cur[i][qp] *
                            normal_face_cur[qp].operator()(0);
            Re_var[1](i) += -JxW_face_cur[qp] * pressure_t(nbc) *
                            phi_face_cur[i][qp] *
                            normal_face_cur[qp].operator()(1);
#if (MESH_DIMENSION == 3)
            Re_var[2](i) += -JxW_face_cur[qp] * pressure_t(nbc) *
                            phi_face_cur[i][qp] *
                            normal_face_cur[qp].operator()(2);
#endif
          }
        }
      }
    }
  }
}

void BoundaryCond::apply_force_x(
    EquationSystems &es, const Elem *elem, unsigned int side,
    const std::vector<std::vector<Real>> &phi_face,
    const std::vector<Real> &JxW_face, const std::vector<Point> &normal_face,
    const std::vector<std::vector<Real>> &phi_face_cur,
    const std::vector<Real> &JxW_face_cur,
    const std::vector<Point> &normal_face_cur,
    DenseSubVector<double> Re_var[]) {

  const MeshBase &mesh = es.get_mesh();

  for (unsigned int qp = 0; qp < qf_points; qp++) {
    vector<boundary_id_type> bc_id_vec;
    mesh.boundary_info->boundary_ids(elem, side, bc_id_vec);
    short int bc_id = bc_id_vec[0];

    for (boundary_id_type nbc = 0; nbc < InputParam::force_x_size; nbc++) {
      if (bc_id == InputParam::force_x_bcs(nbc)) {
        for (std::size_t i = 0; i < phi_face.size(); i++) {
          Re_var[0](i) += -JxW_face[qp] * force_x_t(nbc) * phi_face[i][qp];
        }
      }
    }
  }
}

void BoundaryCond::apply_force_y(
    EquationSystems &es, const Elem *elem, unsigned int side,
    const std::vector<std::vector<Real>> &phi_face,
    const std::vector<Real> &JxW_face, const std::vector<Point> &normal_face,
    const std::vector<std::vector<Real>> &phi_face_cur,
    const std::vector<Real> &JxW_face_cur,
    const std::vector<Point> &normal_face_cur,
    DenseSubVector<double> Re_var[]) {

  const MeshBase &mesh = es.get_mesh();

  for (unsigned int qp = 0; qp < qf_points; qp++) {
    vector<boundary_id_type> bc_id_vec;
    mesh.boundary_info->boundary_ids(elem, side, bc_id_vec);
    short int bc_id = bc_id_vec[0];

    for (boundary_id_type nbc = 0; nbc < InputParam::force_y_size; nbc++) {
      if (bc_id == InputParam::force_y_bcs(nbc)) {
        for (std::size_t i = 0; i < phi_face.size(); i++) {
          Re_var[1](i) += -JxW_face[qp] * force_y_t(nbc) * phi_face[i][qp];
        }
      }
    }
  }
}

#if (MESH_DIMENSION == 3)
void BoundaryCond::apply_force_z(
    EquationSystems &es, const Elem *elem, unsigned int side,
    const std::vector<std::vector<Real>> &phi_face,
    const std::vector<Real> &JxW_face, const std::vector<Point> &normal_face,
    const std::vector<std::vector<Real>> &phi_face_cur,
    const std::vector<Real> &JxW_face_cur,
    const std::vector<Point> &normal_face_cur,
    DenseSubVector<double> Re_var[]) {

  const MeshBase &mesh = es.get_mesh();

  for (unsigned int qp = 0; qp < qf_points; qp++) {
    vector<boundary_id_type> bc_id_vec;
    mesh.boundary_info->boundary_ids(elem, side, bc_id_vec);
    short int bc_id = bc_id_vec[0];

    for (boundary_id_type nbc = 0; nbc < InputParam::force_z_size; nbc++) {
      if (bc_id == InputParam::force_z_bcs(nbc)) {
        for (std::size_t i = 0; i < phi_face.size(); i++) {
          Re_var[2](i) += -JxW_face[qp] * force_z_t(nbc) * phi_face[i][qp];
        }
      }
    }
  }
}
#endif

void BoundaryCond::apply_pressure_force(
    EquationSystems &es, const Elem *elem, unsigned int side,
    const std::vector<std::vector<Real>> &phi_face,
    const std::vector<Real> &JxW_face, const std::vector<Point> &normal_face,
    const std::vector<std::vector<Real>> &phi_face_cur,
    const std::vector<Real> &JxW_face_cur,
    const std::vector<Point> &normal_face_cur,
    DenseSubVector<double> Re_var[]) {

  if (InputParam::pressure_size > 0)
    apply_pressure(es, elem, side, phi_face, JxW_face, normal_face,
                   phi_face_cur, JxW_face_cur, normal_face_cur, Re_var);
  if (InputParam::force_x_size > 0)
    apply_force_x(es, elem, side, phi_face, JxW_face, normal_face, phi_face_cur,
                  JxW_face_cur, normal_face_cur, Re_var);
  if (InputParam::force_y_size > 0)
    apply_force_y(es, elem, side, phi_face, JxW_face, normal_face, phi_face_cur,
                  JxW_face_cur, normal_face_cur, Re_var);
#if (MESH_DIMENSION == 3)
  if (InputParam::force_z_size > 0)
    apply_force_z(es, elem, side, phi_face, JxW_face, normal_face, phi_face_cur,
                  JxW_face_cur, normal_face_cur, Re_var);
#endif
}

void BoundaryCond::compute_bcid(EquationSystems &es) {
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  System &system = es.get_system<System>("bcidSystem");
  unsigned int bcid_var = system.variable_number("bcidVar");
  const DofMap &dof_map = system.get_dof_map();
  std::vector<dof_id_type> dof_indices;

  FEType fe_type = dof_map.variable_type(bcid_var);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_type));
  QGauss qface(dim - 1, fe_type.default_quadrature_order());
  fe_face->attach_quadrature_rule(&qface);

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {
    const Elem *elem = *el;

    dof_map.dof_indices(elem, dof_indices, bcid_var);

    const unsigned int n_var_dofs = dof_indices.size();

    fe->reinit(elem);

    for (unsigned int side = 0; side < elem->n_sides(); side++)
      if (elem->neighbor_ptr(side) == libmesh_nullptr) {
        const std::vector<std::vector<Real>> &phi_face = fe_face->get_phi();
        fe_face->reinit(elem, side);
        vector<boundary_id_type> bc_id_vec;
        mesh.boundary_info->boundary_ids(elem, side, bc_id_vec);
        short int bc_id =
            bc_id_vec[0]; // mesh.boundary_info->boundary_id(elem, side);

        if (bc_id == InputParam::torsion_bc) {
          for (unsigned int j = 0; j < n_var_dofs; j++) {
            if (fabs(phi_face[j][0]) > 0.01)
              system.solution->set(dof_indices[j], 1000.0);
          }
        }
      }
  }

  system.solution->close();
  system.solution->localize(*system.current_local_solution);
}

void BoundaryCond::apply_torsion_resid(EquationSystems &es, const Elem *elem,
                                       DenseSubVector<double> Re_var[]) {

  const MeshBase &mesh = es.get_mesh();
  NonlinearImplicitSystem &system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");
  const unsigned int u_var = system.variable_number("u");
  const DofMap &dof_map = system.get_dof_map();
  std::vector<std::vector<dof_id_type>> dof_indices_var(MESH_DIMENSION);

  System &bcid_system = es.get_system<System>("bcidSystem");
  unsigned int bcid_var = bcid_system.variable_number("bcidVar");
  const DofMap &dof_map_bcid = bcid_system.get_dof_map();
  std::vector<dof_id_type> dof_indices_bcid;

  for (unsigned int var = 0; var < MESH_DIMENSION; var++) {
    dof_map.dof_indices(elem, dof_indices_var[var], var);
  }

  dof_map_bcid.dof_indices(elem, dof_indices_bcid, bcid_var);

  for (unsigned int dof_i = 0; dof_i < num_nodes; dof_i++) {
    double u_cur = system.current_solution(dof_indices_var[0][dof_i]);
    double v_cur = system.current_solution(dof_indices_var[1][dof_i]);
#if (MESH_DIMENSION == 3)
    double w_cur = system.current_solution(dof_indices_var[2][dof_i]);
#endif

    double boundary_cur = bcid_system.current_solution(dof_indices_bcid[dof_i]);

    if (boundary_cur > 900.0) {
      Point pj;
      pj = elem->point(dof_i);

      double Dis = sqrt(pow(pj(0), 2) + pow(pj(1), 2));
      double sinTheta = (pj(1)) / Dis;
      double cosTheta = (pj(0)) / Dis;

      double res0 = Re_var[0](dof_i);
      double res1 = Re_var[1](dof_i);

      if (InputParam::torsion_type == 0) {
        Re_var[0](dof_i) = res0 * cosTheta + res1 * sinTheta;
        Re_var[1](dof_i) = (-u_cur * sinTheta + v_cur * cosTheta) -
                           Dis * InputParam::torsion_t;
      } else if (InputParam::torsion_type == 1) {
        Re_var[0](dof_i) = u_cur + Dis * InputParam::torsion_t * sinTheta;
        Re_var[1](dof_i) = v_cur - Dis * InputParam::torsion_t * cosTheta;
      } else if (InputParam::torsion_type == 2) {
        Re_var[0](dof_i) = u_cur - InputParam::torsion_t;
      } else if (InputParam::torsion_type == 3) {
        Re_var[0](dof_i) = u_cur + InputParam::torsion_t *
                                       sin(0.5 * M_PI * pj(0)) *
                                       cos(0.5 * M_PI * pj(1));
        Re_var[1](dof_i) = v_cur - InputParam::torsion_t *
                                       cos(0.5 * M_PI * pj(0)) *
                                       sin(0.5 * M_PI * pj(1));
      } else if (InputParam::torsion_type == 4) {
        Re_var[0](dof_i) = 1.0 * (u_cur - InputParam::torsion_t);
        Re_var[1](dof_i) = 1.0 * v_cur;
#if (MESH_DIMENSION == 3)
        Re_var[2](dof_i) = 1.0 * w_cur;
#endif

      }

      else if (InputParam::torsion_type == 10) {
        Re_var[MESH_DIMENSION](dof_i) = 0.0;
      }
    }
  }
}

void BoundaryCond::apply_torsion_jacob(
    EquationSystems &es, const Elem *elem,
    DenseSubMatrix<Number> Ke_var[][MESH_DIMENSION + 1]) {
  for (unsigned int dof_i = 0; dof_i < num_nodes; dof_i++) {

    System &bcid_system = es.get_system<System>("bcidSystem");
    unsigned int bcid_var = bcid_system.variable_number("bcidVar");
    const DofMap &dof_map_bcid = bcid_system.get_dof_map();
    std::vector<dof_id_type> dof_indices_bcid;
    dof_map_bcid.dof_indices(elem, dof_indices_bcid, bcid_var);

    double boundary_cur = bcid_system.current_solution(dof_indices_bcid[dof_i]);

    if (boundary_cur > 900.0) {
      Point pj;
      pj = elem->point(dof_i);

      double Dis = sqrt(pow(pj(0), 2) + pow(pj(1), 2));
      double sinTheta = (pj(1)) / Dis;
      double cosTheta = (pj(0)) / Dis;

      if (InputParam::torsion_type == 0) {
        for (unsigned int dof_j = 0; dof_j < num_nodes; dof_j++) {
          double K0i[MESH_DIMENSION + 1], K1i[MESH_DIMENSION + 1];
          for (unsigned int i = 0; i < MESH_DIMENSION + 1; i++) {
            K0i[i] = Ke_var[0][i](dof_i, dof_j);
            K1i[i] = Ke_var[1][i](dof_i, dof_j);
          }
          for (unsigned int i = 0; i < MESH_DIMENSION + 1; i++) {
            Ke_var[0][i](dof_i, dof_j) = K0i[i] * cosTheta + K1i[i] * sinTheta;
            Ke_var[1][i](dof_i, dof_j) = 0.0;
          }
        }
        Ke_var[1][0](dof_i, dof_i) = -sinTheta;
        Ke_var[1][1](dof_i, dof_i) = cosTheta;
      } else if (InputParam::torsion_type == 1 ||
                 InputParam::torsion_type == 3) {
        for (unsigned int dof_j = 0; dof_j < num_nodes; dof_j++) {
          for (unsigned int i = 0; i < MESH_DIMENSION + 1; i++) {
            Ke_var[0][i](dof_i, dof_j) = 0.0;
            Ke_var[1][i](dof_i, dof_j) = 0.0;
          }
        }
        Ke_var[0][0](dof_i, dof_i) = 1.0;
        Ke_var[1][1](dof_i, dof_i) = 1.0;
      }

      else if (InputParam::torsion_type == 2) {
        for (unsigned int dof_j = 0; dof_j < num_nodes; dof_j++) {
          for (unsigned int i = 0; i < MESH_DIMENSION + 1; i++) {
            Ke_var[0][i](dof_i, dof_j) = 0.0;
          }
        }
        Ke_var[0][0](dof_i, dof_i) = 1.0;
      }

      else if (InputParam::torsion_type == 4) {
        for (unsigned int dof_j = 0; dof_j < num_nodes; dof_j++) {
          for (unsigned int i = 0; i < MESH_DIMENSION + 1; i++) {
            Ke_var[0][i](dof_i, dof_j) = 0.0;
            Ke_var[1][i](dof_i, dof_j) = 0.0;
#if (MESH_DIMENSION == 3)
            Ke_var[2][i](dof_i, dof_j) = 0.0;
#endif
          }
        }
        Ke_var[0][0](dof_i, dof_i) = 1.0;
        Ke_var[1][1](dof_i, dof_i) = 1.0;
#if (MESH_DIMENSION == 3)
        Ke_var[2][2](dof_i, dof_i) = 1.0;
#endif
      }

      else if (InputParam::torsion_type == 10) {
        for (unsigned int dof_j = 0; dof_j < num_nodes; dof_j++) {
          for (unsigned int i = 0; i < MESH_DIMENSION + 1; i++) {
            Ke_var[MESH_DIMENSION][i](dof_i, dof_j) = 0.0;
          }
        }
        Ke_var[MESH_DIMENSION][MESH_DIMENSION](dof_i, dof_i) = 1.0;
      }
    }
  }
}
