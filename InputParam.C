#include "InputParam.h"

using namespace libMesh;
using namespace std;

double InputParam::mesh_scale;
std::string InputParam::mesh_file_name, InputParam::fibre_file_name, InputParam::perm_file_name,
    InputParam::sheet_file_name;
unsigned int InputParam::mesh_centre;
unsigned int InputParam::n_solves,InputParam::n_total;
unsigned int InputParam::output_terminal;
double InputParam::ttime, InputParam::dt,InputParam::time_per,InputParam::omega;
int InputParam::time_itr;
int InputParam::inertia, InputParam::trans_soln, InputParam::brinkman, InputParam::aniso_perm;

double InputParam::nonlinear_abs_tol, InputParam::nonlinear_rel_tol;
unsigned int InputParam::nonlinear_max_its;

int InputParam::strain_model;

double InputParam::aHO, InputParam::bHO, InputParam::afHO, InputParam::bfHO,
    InputParam::asHO, InputParam::bsHO, InputParam::afsHO, InputParam::bfsHO;

double InputParam::G;

double InputParam::alpha_stab, InputParam::rho_s;

double InputParam::kappa_0;

double InputParam::viscocity;

    boundary_id_type InputParam::clamp_x_size,
    InputParam::clamp_y_size;
DenseVector<boundary_id_type> InputParam::clamp_x_bcs, InputParam::clamp_y_bcs;
#if (MESH_DIMENSION == 3)
boundary_id_type InputParam::clamp_z_size;
DenseVector<boundary_id_type> InputParam::clamp_z_bcs;
#endif

boundary_id_type InputParam::torsion_bc, InputParam::torsion_type;
double InputParam::torsion_t;

boundary_id_type InputParam::force_x_size, InputParam::force_y_size;
DenseVector<boundary_id_type> InputParam::force_x_bcs, InputParam::force_y_bcs;
DenseVector<double> InputParam::force_x_values, InputParam::force_y_values;
#if (MESH_DIMENSION == 3)
boundary_id_type InputParam::force_z_size;
DenseVector<boundary_id_type> InputParam::force_z_bcs;
DenseVector<double> InputParam::force_z_values;
#endif

boundary_id_type InputParam::pressure_size;
DenseVector<boundary_id_type> InputParam::pressure_bcs;
DenseVector<double> InputParam::pressure_values;

double InputParam::t_load, InputParam::t_end_diastole, InputParam::pressure_sys;
int InputParam::pressure_current;

double InputParam::beta_s;

int InputParam::active_tension;
double InputParam::Ta_max;

double InputParam::V_bead, InputParam::a_bead,InputParam::bead_disp;

double InputParam::alpha_fib;

int InputParam::viscoelastic;
double InputParam::tau_visela, InputParam::beta_visela;
double InputParam::permeability;

DenseVector<double> InputParam::Vtaubya, InputParam::VabyD;
int InputParam::vtaubya_size, InputParam::vabyd_size;

int InputParam::var_v_a;

int InputParam::porous;

InputParam::InputParam() {}

InputParam::~InputParam() {}

void InputParam::read_input() {
  GetPot infile1("input.in");
  std::string input_file_name = infile1("input_file_name", "input_LV.in");
  GetPot infile(input_file_name);

  nonlinear_abs_tol = infile("nonlinear_abs_tol", 1.e-8);
  nonlinear_rel_tol = infile("nonlinear_rel_tol", 1.e-8);
  nonlinear_max_its = infile("nonlinear_max_its", 50);

  mesh_file_name = infile("mesh_file_name", "LV_hex_0.625.e");
  mesh_scale = infile("mesh_scale", 0.1);
  mesh_centre = infile("mesh_centre", 1);
  n_solves = infile("n_solves", 10);
  dt = infile("dt", 0.01);

  aniso_perm = infile("aniso_perm", 0);
  perm_file_name = infile("perm_file_name", "frame_perm.e");

  time_per = infile("time_per", 0.8);

  output_terminal = infile("output_terminal", 0);

  inertia = infile("inertia", 0);
  trans_soln = infile("trans_soln", 1);

  strain_model = infile("strain_model", 1);

  fibre_file_name = infile("fibre_file_name", "fibreDir.txt");
  sheet_file_name = infile("sheet_file_name", "sheetDir.txt");

  aHO = infile("aHO", 1.0);
  bHO = infile("bHO", 1.0);
  afHO = infile("afHO", 1.0);
  bfHO = infile("bfHO", 1.0);
  asHO = infile("asHO", 1.0);
  bsHO = infile("bsHO", 1.0);
  afsHO = infile("afsHO", 1.0);
  bfsHO = infile("bfsHO", 1.0);

  n_total = infile("n_total", 10);

  G = infile("G", 1540.0);

  rho_s = infile("rho_s", 1.0);

  kappa_0 = infile("kappa_0", 0.0);

  alpha_stab = infile("alpha_stab", 0.028);

  bead_disp = infile("bead_disp", 0.1);

  clamp_x_size =
      cast_int<boundary_id_type>(infile.vector_variable_size("clamp_x_bcs"));
  clamp_y_size =
      cast_int<boundary_id_type>(infile.vector_variable_size("clamp_y_bcs"));

  if (clamp_x_size > 0) {
    clamp_x_bcs.resize(clamp_x_size);
    for (boundary_id_type nbc = 0; nbc < clamp_x_size; nbc++) {
      clamp_x_bcs(nbc) =
          cast_int<boundary_id_type>(infile("clamp_x_bcs", 1, nbc));
    }
  }

  if (clamp_y_size > 0) {
    clamp_y_bcs.resize(clamp_y_size);
    for (boundary_id_type nbc = 0; nbc < clamp_y_size; nbc++) {
      clamp_y_bcs(nbc) =
          cast_int<boundary_id_type>(infile("clamp_y_bcs", 1, nbc));
    }
  }

#if (MESH_DIMENSION == 3)
  clamp_z_size =
      cast_int<boundary_id_type>(infile.vector_variable_size("clamp_z_bcs"));
  if (clamp_z_size > 0) {
    clamp_z_bcs.resize(clamp_z_size);
    for (boundary_id_type nbc = 0; nbc < clamp_z_size; nbc++) {
      clamp_z_bcs(nbc) =
          cast_int<boundary_id_type>(infile("clamp_z_bcs", 1, nbc));
    }
  }
#endif

  torsion_bc = infile("torsion_bc", 10000);
  torsion_t = infile("torsion_value", 0.0);
  torsion_type = infile("torsion_type", 0);

  force_x_size =
      cast_int<boundary_id_type>(infile.vector_variable_size("force_x_bcs"));
  force_y_size =
      cast_int<boundary_id_type>(infile.vector_variable_size("force_y_bcs"));

  if (force_x_size > 0) {
    force_x_bcs.resize(force_x_size);
    force_x_values.resize(force_x_size);
    for (boundary_id_type nbc = 0; nbc < force_x_size; nbc++) {
      force_x_bcs(nbc) = infile("force_x_bcs", 1, nbc);
      force_x_values(nbc) = infile("force_x_values", 0.0, nbc);
    }
  }

  if (force_y_size > 0) {
    force_y_bcs.resize(force_y_size);
    force_y_values.resize(force_x_size);
    for (boundary_id_type nbc = 0; nbc < force_y_size; nbc++) {
      force_y_bcs(nbc) = infile("force_y_bcs", 1, nbc);
      force_y_values(nbc) = infile("force_y_values", 0.0, nbc);
    }
  }

#if (MESH_DIMENSION == 3)
  force_z_size =
      cast_int<boundary_id_type>(infile.vector_variable_size("force_z_bcs"));
  if (force_z_size > 0) {
    force_z_bcs.resize(force_z_size);
    force_z_values.resize(force_x_size);
    for (boundary_id_type nbc = 0; nbc < force_z_size; nbc++) {
      force_z_bcs(nbc) = infile("force_z_bcs", 1, nbc);
      force_z_values(nbc) = infile("force_z_values", 0.0, nbc);
    }
  }
#endif

  pressure_size =
      cast_int<boundary_id_type>(infile.vector_variable_size("pressure_bcs"));

  if (pressure_size > 0) {
    pressure_bcs.resize(pressure_size);
    pressure_values.resize(pressure_size);
    for (boundary_id_type nbc = 0; nbc < pressure_size; nbc++) {
      pressure_bcs(nbc) = infile("pressure_bcs", 1, nbc);
      pressure_values(nbc) = infile("pressure_values", 0.0, nbc);
    }
  }

  t_load = infile("t_load", 0.8);
  t_end_diastole = infile("t_end_diastole", 0.8);
  pressure_current = infile("pressure_current", 0);
  pressure_sys = infile("pressure_sys", -10666.0);

  beta_s = infile("beta_s", 1.0e6);

  active_tension = infile("active_tension", 0);
  Ta_max = infile("Ta_max", 0.0);

  V_bead = infile("V_bead", 1.0);
  a_bead = infile("a_bead", 1.0);

  alpha_fib = infile("alpha_fib", 1.0);

  viscoelastic = infile("viscoelastic", 0);
  tau_visela = infile("tau_visela", 1.0);
  beta_visela = infile("beta_visela", 1.0);

  permeability = infile("permeability", 0.1);

  vtaubya_size = cast_int<double>(infile.vector_variable_size("vtaubya"));

  vabyd_size = cast_int<double>(infile.vector_variable_size("vabyd"));

  if (vtaubya_size > 0) {
    Vtaubya.resize(vtaubya_size);
    for (int nbc = 0; nbc < vtaubya_size; nbc++) {
      Vtaubya(nbc) = infile("vtaubya", 1.0, nbc);
    }
  }

  if (vabyd_size > 0) {
    VabyD.resize(vabyd_size);
    for (int nbc = 0; nbc < vabyd_size; nbc++) {
      VabyD(nbc) = infile("vabyd", 1.0, nbc);
    }
  }

  var_v_a = infile("var_v_a", 0);

  porous = infile("porous", 0);
  brinkman = infile("brinkman", 0);
  viscocity = infile("viscocity", 0.001);

  infile.print();
}

void InputParam::read_mesh(Mesh &mesh) {
  mesh.allow_renumbering(false);
  mesh.read(mesh_file_name, NULL);
  MeshTools::Modification::scale(mesh, mesh_scale);

  if (mesh_centre == 1) {
    //MeshTools::BoundingBox bbox = MeshTools::bounding_box(mesh);

    BoundingBox bbox = MeshTools::create_bounding_box(mesh);

    libMesh::Point lower = bbox.min();
    libMesh::Point upper = bbox.max();
    libMesh::out << "         mesh bounding box lower = (" << lower(0) << " , "
                 << lower(1) << " , " << lower(2) << ") cm\n"
                 << "         mesh bounding box upper = (" << upper(0) << " , "
                 << upper(1) << " , " << upper(2) << ") cm\n";
    MeshTools::Modification::translate(mesh, -0.5 * (upper(0) + lower(0)),
                                       -0.5 * (upper(1) + lower(1)),
                                       -0.5 * (upper(2) + lower(2)));
  }
}

void InputParam::read_mesh_perm(ExodusII_IO & exo_io, Mesh &mesh)
{
  mesh.allow_renumbering(false);

  if (mesh.processor_id() == 0)
    exo_io.read(perm_file_name);
  MeshCommunication().broadcast(mesh);
  mesh.prepare_for_use();

  MeshTools::Modification::scale(mesh, mesh_scale);

  if (mesh_centre == 1)
  {
    // MeshTools::BoundingBox bbox = MeshTools::bounding_box(mesh);

    BoundingBox bbox = MeshTools::create_bounding_box(mesh);

    libMesh::Point lower = bbox.min();
    libMesh::Point upper = bbox.max();
    libMesh::out << "         mesh bounding box lower = (" << lower(0) << " , "
                 << lower(1) << " , " << lower(2) << ") cm\n"
                 << "         mesh bounding box upper = (" << upper(0) << " , "
                 << upper(1) << " , " << upper(2) << ") cm\n";
    MeshTools::Modification::translate(mesh, -0.5 * (upper(0) + lower(0)),
                                       -0.5 * (upper(1) + lower(1)),
                                       -0.5 * (upper(2) + lower(2)));
  }
}
