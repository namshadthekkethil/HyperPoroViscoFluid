#include "PoroElastic.h"

using namespace libMesh;
using namespace std;

string PoroElastic::approx_order, PoroElastic::fe_family;
double PoroElastic::phi_crit, PoroElastic::phi_0;
double PoroElastic::q1, PoroElastic::q2, PoroElastic::q3;
double PoroElastic::ttime, PoroElastic::dt;
double PoroElastic::source_const, PoroElastic::source_prssr,
    PoroElastic::sink_const, PoroElastic::sink_prssr;
double PoroElastic::mesh_volume;
double PoroElastic::permeability;
unsigned int PoroElastic::ppore_penalty_on;
unsigned int PoroElastic::inertia, PoroElastic::incomp_type,
    PoroElastic::pmat_type;
double PoroElastic::rho_s, PoroElastic::kappa_0;
DenseMatrix<double> PoroElastic::II;
int PoroElastic::time_itr;
int PoroElastic::viscous_on;
double PoroElastic::epsilon, PoroElastic::beta_visc;
DenseVector<double> PoroElastic::vasc_nodes_x, PoroElastic::vasc_nodes_y,
    PoroElastic::vasc_nodes_z, PoroElastic::edges_1, PoroElastic::edges_2;

int PoroElastic::brinkman;
double PoroElastic::viscocity;
int PoroElastic::read_source, PoroElastic::read_tree, PoroElastic::read_permeability;

double PoroElastic::coord_source[4][21990];
DenseMatrix<double> PoroElastic::porous_data;
DenseVector<vector<double>> PoroElastic::source_data;

vector<int> PoroElastic::near_vess;
vector<double> PoroElastic::source_vess;

PoroElastic::PoroElastic() {}

PoroElastic::~PoroElastic() {}

void PoroElastic::read_input(int rank)
{
  GetPot infile1("input.in");
  std::string input_file_name = infile1("input_file_name", "input_LV.in");
  GetPot infile(input_file_name);

  approx_order = infile("approx_order", "FIRST");
  fe_family = infile("fe_family", "LAGRANGE");

  phi_crit = infile("phi_crit", 0.001);
  phi_0 = infile("phi_0", 0.001);
  source_const = infile("source_const", 0.1);
  source_prssr = infile("source_prssr", 0.1);
  sink_const = infile("sink_const", 0.0);
  sink_prssr = infile("sink_prssr", 0.0);
  permeability = InputParam::permeability; //infile("permeability", 0.00000001);
  ppore_penalty_on = infile("ppore_penalty_on", 0);

  kappa_0 = infile("kappa_0", 0.0);

  q1 = infile("q1", 0.022e4);
  q2 = infile("q2", 1.009e4);
  q3 = infile("q3", 80.0);

  pmat_type = infile("pmat_type", 0);

  inertia = infile("inertia", 1);
  rho_s = infile("rho_s", 1.0);

  incomp_type = infile("incomp_type", 0);

  epsilon = infile("epsilon", 0.001);

  viscous_on = infile("viscous_on", 0);
  beta_visc = infile("beta_visc", 0.0);

  brinkman = infile("brinkman", 0);
  viscocity = infile("viscocity", 0.001);

  read_source = infile("read_source", 0);
  read_tree = infile("read_tree", 0);
  read_permeability = infile("read_permeability", 0);
}

void PoroElastic::define_systems(EquationSystems &es, int rank)
{

  read_input(rank);

  LinearImplicitSystem &system_norm =
      es.add_system<LinearImplicitSystem>("NormalSystem");
  system_norm.add_variable("NormalVar",
                           Utility::string_to_enum<Order>(approx_order),
                           Utility::string_to_enum<FEFamily>(fe_family));

  System &system_pmat = es.add_system<System>("pMatSystem");
  System &system_ppen = es.add_system<System>("pPenSystem");
  System &system_ppore = es.add_system<System>("pPoreSystem");
  System &system_source = es.add_system<System>("sourceSystem");

  System &system_mmono = es.add_system<System>("mMonoSystem");

  LinearImplicitSystem &system_darcy =
      es.add_system<LinearImplicitSystem>("darcy_velocity");
  LinearImplicitSystem &system_dmdt =
      es.add_system<LinearImplicitSystem>("dmdtSystem");

  LinearImplicitSystem &system_mexp =
      es.add_system<LinearImplicitSystem>("mExpSystem");
  system_mexp.add_variable("mExpVar", FIRST, LAGRANGE);
  system_mexp.attach_assemble_function(assemble_mexp);

  system_pmat.add_variable("pMatVar",
                           Utility::string_to_enum<Order>(approx_order),
                           Utility::string_to_enum<FEFamily>(fe_family));
  system_ppen.add_variable("pPenVar",
                           Utility::string_to_enum<Order>(approx_order),
                           Utility::string_to_enum<FEFamily>(fe_family));
  system_ppore.add_variable("pPoreVar",
                            Utility::string_to_enum<Order>(approx_order),
                            Utility::string_to_enum<FEFamily>(fe_family));
  system_darcy.add_variable("w_x", Utility::string_to_enum<Order>(approx_order),
                            Utility::string_to_enum<FEFamily>(fe_family));
  system_darcy.add_variable("w_y", Utility::string_to_enum<Order>(approx_order),
                            Utility::string_to_enum<FEFamily>(fe_family));
#if (MESH_DIMENSION == 3)
  system_darcy.add_variable("w_z", Utility::string_to_enum<Order>(approx_order),
                            Utility::string_to_enum<FEFamily>(fe_family));
#endif
  system_source.add_variable("sourceVar", CONSTANT, MONOMIAL);
  system_dmdt.add_variable("dmdtVar",
                           Utility::string_to_enum<Order>(approx_order),
                           Utility::string_to_enum<FEFamily>(fe_family));

  system_mmono.add_variable("mMonoVar", CONSTANT, MONOMIAL);

  LinearImplicitSystem &system_flow =
      es.add_system<LinearImplicitSystem>("flowSystem");
  system_flow.add_variable("mVar", FIRST, LAGRANGE);

  System &system_flow_old = es.add_system<System>("flowOldSystem");
  system_flow_old.add_variable("mOldVar", FIRST, LAGRANGE);

  LinearImplicitSystem &system_porous =
      es.add_system<LinearImplicitSystem>("porousSystem");
  system_porous.add_variable("wxVar", FIRST, LAGRANGE);
  system_porous.add_variable("wyVar", FIRST, LAGRANGE);
#if (MESH_DIMENSION == 3)
  system_porous.add_variable("wzVar", FIRST, LAGRANGE);
#endif
  system_porous.add_variable("mAddedVar", FIRST, LAGRANGE);

  LinearImplicitSystem &system_porous_old =
      es.add_system<LinearImplicitSystem>("porousOldSystem");
  system_porous_old.add_variable("wxOldVar", FIRST, LAGRANGE);
  system_porous_old.add_variable("wyOldVar", FIRST, LAGRANGE);
#if (MESH_DIMENSION == 3)
  system_porous_old.add_variable("wzOldVar", FIRST, LAGRANGE);
#endif
  system_porous_old.add_variable("mAddedOldVar", FIRST, LAGRANGE);

  System &phi_system = es.add_system<System>("phiSystem");
  phi_system.add_variable("phiVar", CONSTANT, MONOMIAL);

  System &system_K = es.add_system<System>("KSystem");
  system_K.add_variable("K00", CONSTANT, MONOMIAL);
  system_K.add_variable("K01", CONSTANT, MONOMIAL);
#if (MESH_DIMENSION == 3)
  system_K.add_variable("K02", CONSTANT, MONOMIAL);
#endif
  system_K.add_variable("K10", CONSTANT, MONOMIAL);
  system_K.add_variable("K11", CONSTANT, MONOMIAL);
#if (MESH_DIMENSION == 3)
  system_K.add_variable("K12", CONSTANT, MONOMIAL);
  system_K.add_variable("K20", CONSTANT, MONOMIAL);
  system_K.add_variable("K21", CONSTANT, MONOMIAL);
  system_K.add_variable("K22", CONSTANT, MONOMIAL);
#endif

System &system_flowlarge = es.add_system<System>("flowLargeSystem");
system_flowlarge.add_variable("flowLargeVar", CONSTANT, MONOMIAL);

  if (brinkman == 1)
  {
    LinearImplicitSystem &system_delw =
        es.add_system<LinearImplicitSystem>("systemDelW");
#if (MESH_DIMENSION == 2)
    system_delw.add_variable("delW00", FIRST, LAGRANGE);
    system_delw.add_variable("delW01", FIRST, LAGRANGE);
    system_delw.add_variable("delW10", FIRST, LAGRANGE);
    system_delw.add_variable("delW11", FIRST, LAGRANGE);
#elif (MESH_DIMENSION == 3)
    system_delw.add_variable("delW00", FIRST, LAGRANGE);
    system_delw.add_variable("delW01", FIRST, LAGRANGE);
    system_delw.add_variable("delW02", FIRST, LAGRANGE);
    system_delw.add_variable("delW10", FIRST, LAGRANGE);
    system_delw.add_variable("delW11", FIRST, LAGRANGE);
    system_delw.add_variable("delW12", FIRST, LAGRANGE);
    system_delw.add_variable("delW20", FIRST, LAGRANGE);
    system_delw.add_variable("delW21", FIRST, LAGRANGE);
    system_delw.add_variable("delW22", FIRST, LAGRANGE);
#endif

    system_delw.attach_assemble_function(assemble_delw);
  }

  System &system_aha = es.add_system<System>("ahaSystem");
  system_aha.add_variable("ahaVar", CONSTANT, MONOMIAL);

  II.resize(MESH_DIMENSION, MESH_DIMENSION);
  II(0, 0) = 1.0;
  II(1, 1) = 1.0;
#if (MESH_DIMENSION == 3)
  II(2, 2) = 1.0;
#endif

  system_norm.attach_assemble_function(assemble_normal);
  system_darcy.attach_assemble_function(assemble_darcy);
  system_dmdt.attach_assemble_function(assemble_dmdt);

  system_flow.attach_assemble_function(assemble_flow);
  system_porous.attach_assemble_function(assemble_porous);
}

double PoroElastic::compute_ppen(double jphi_cur, double pnorm_cur)
{
  if (ppore_penalty_on == 0)
    return 0.0;
  else
  {
    double dirac_delta = 0.0;
    dirac_delta = (1.0 / M_PI) *
                  (epsilon / ((phi_crit - jphi_cur) * (phi_crit - jphi_cur) +
                              epsilon * epsilon));
    return -dirac_delta * pnorm_cur;
  }
}

double PoroElastic::compute_pmat(double jphi_cur)
{
  if (pmat_type == 0)
    return 0.0;
  else if (pmat_type == 1)
    return (kappa_0 * (jphi_cur - phi_0));
  else if (pmat_type == 2)
    return q1 * exp(q3 * jphi_cur) + q2 * log(q3 * jphi_cur) -
           (q1 * exp(q3 * phi_0) + q2 * log(q3 * phi_0));
}

double PoroElastic::compute_dpmatdjphi(double jphi_cur)
{
  if (pmat_type == 0)
    return 0.0;
  else if (pmat_type == 1)
    return kappa_0;
  else if (pmat_type == 2)
    return q1 * q3 * exp(q3 * jphi_cur) + (q2 / jphi_cur);
}

double PoroElastic::compute_dppendjphi(double jphi_cur, double pnorm_cur)
{
  if (ppore_penalty_on == 0)
    return 0.0;
  else
    return pnorm_cur * (epsilon / M_PI) *
           ((2.0 * (jphi_cur - phi_crit)) /
            pow(pow(jphi_cur - phi_crit, 2) + epsilon * epsilon, 2));
}

double PoroElastic::compute_dppendlambda(double jphi_cur)
{
  if (ppore_penalty_on == 0)
    return 0.0;
  else
    return (1.0 / M_PI) *
           (epsilon / ((phi_crit - jphi_cur) * (phi_crit - jphi_cur) +
                       epsilon * epsilon));
}

double PoroElastic::compute_d2pmatdjphi2(double jphi_cur)
{
  if (pmat_type == 0)
    return 0.0;
  else if (pmat_type == 1)
    return 0.0;
  else if (pmat_type == 2)
    q1 *q3 *q3 *exp(q3 * jphi_cur) - (q2 / (jphi_cur * jphi_cur));
}

double PoroElastic::compute_ppore(double pmono_cur, double pmat_cur,
                                  double ppen_cur)
{
  return pmono_cur + pmat_cur + ppen_cur;
}

void PoroElastic::assemble_normal(
    EquationSystems &es, const std::string &libmesh_dbg_var(system_name))
{

  // Get a constant reference to the mesh object.
  const MeshBase &mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  LinearImplicitSystem &norm_system =
      es.get_system<LinearImplicitSystem>("NormalSystem");

  unsigned int norm_var;

  // Numeric ids corresponding to each variable in the system
  norm_var = norm_system.variable_number("NormalVar");

  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_vel_type = norm_system.variable_type(norm_var);

  // Build a Finite Element object of the specified type for
  // the velocity variables.
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule(dim, fe_vel_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_vel_type));
  QGauss qface(dim - 1,
               fe_vel_type.default_quadrature_order()); // Not sure what the
                                                        // best accuracy is here

  fe_face->attach_quadrature_rule(&qface);

  // Here we define some references to cell-specific data that
  // will be used to assemble the linear system.
  //
  // The element Jacobian * quadrature weight at each integration point.
  const std::vector<double> &JxW = fe_vel->get_JxW();

  // The element shape functions evaluated at the quadrature points.
  const std::vector<std::vector<double>> &phi = fe_vel->get_phi();

  // The element shape function gradients for the velocity
  // variables evaluated at the quadrature points.
  const std::vector<std::vector<RealGradient>> &dphi = fe_vel->get_dphi();

  const DofMap &dof_map = norm_system.get_dof_map();

  DenseMatrix<double> Ke;
  DenseVector<double> Fe;

  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_norm;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem *elem = *el;

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);
    dof_map.dof_indices(elem, dof_indices_norm, norm_var);

    unsigned int n_dofs;

    n_dofs = dof_indices.size();

    // Compute the element-specific data for the current
    // element.  This involves computing the location of the
    // quadrature points (q_point) and the shape functions
    // (phi, dphi) for the current element.
    fe_vel->reinit(elem);

    // Zero the element matrix and right-hand side before
    // summing them.  We use the resize member here because
    // the number of degrees of freedom might have changed from
    // the last element.  Note that this will be the case if the
    // element type is different (i.e. the last element was a
    // triangle, now we are on a quadrilateral).
    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    // Now we will build the element matrix and right-hand-side.
    // Constructing the RHS requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {

      for (unsigned int i = 0; i < n_dofs; i++)
      {

        for (unsigned int j = 0; j < n_dofs; j++)
        {
          Ke(i, j) += JxW[qp] * (dphi[i][qp] * dphi[j][qp]);
        }
      }

    } // end of the quadrature point qp-loop

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

          if (bc_id == 4096)
          {
            for (std::size_t i = 0; i < phi_face.size(); i++)
            {
              Fe(i) += JxW_face[qp] * phi_face[i][qp] * 1.0e12 * 1.0;
            }
          }
          if (bc_id == 4096 || bc_id == 4097)
          {
            for (std::size_t i = 0; i < phi_face.size(); i++)
            {
              for (std::size_t j = 0; j < phi_face.size(); j++)
              {
                Ke(i, j) +=
                    JxW_face[qp] * 1.0e12 * phi_face[i][qp] * phi_face[j][qp];
              }
            }
          }
        }
      }

    norm_system.matrix->add_matrix(Ke, dof_indices);
    norm_system.rhs->add_vector(Fe, dof_indices);
  } // end of element loop
}

double PoroElastic::compute_source(double ppore_cur, double norm_cur)
{
  // double source_cur = 0.0;
  // if(ttime<=0.2)
  //{
  // if(norm_cur<0.3333)
  // source_cur=-0.403*ttime+0.242;
  // else if(norm_cur <0.6667)
  // source_cur=-6.955*pow(ttime,3)+1.417*pow(ttime,2)-0.109*ttime+0.254;
  // else
  // source_cur=-12.28*pow(ttime,3)+3.072*pow(ttime,2)-0.056*ttime+0.281;
  // }
  // else if(ttime<=0.58)
  //{
  // if(norm_cur<0.3333)
  // source_cur=-517.448*pow(ttime,5)+820.825*pow(ttime,4)-506.238*pow(ttime,3)+151.751*pow(ttime,2)-22.262*ttime+1.416;
  // else if(norm_cur <0.6667)
  // source_cur=-5.812*pow(ttime,3)+6.092*pow(ttime,2)-2.168*ttime+0.46;
  // else
  // source_cur=70.187*pow(ttime,4)-98.525*pow(ttime,3)+50.313*pow(ttime,2)-11.215*ttime+1.19;
  // }
  // else
  //{
  // if(norm_cur<0.3333)
  // source_cur=2309.718*pow(ttime,4)-5968.808*pow(ttime,3)+5774.409*pow(ttime,2)-2479.339*ttime+398.439;
  // else if(norm_cur <0.6667)
  // source_cur=8.858*pow(ttime,3)-29.074*pow(ttime,2)+27.955*ttime-8.034;
  // else
  // source_cur=-107.978*pow(ttime,3)+188.636*pow(ttime,2)-105.431*ttime+19.122;
  // }

  // source_cur /= mesh_volume;

  // return source_cur;

  return source_const * (source_prssr - ppore_cur) -
         sink_const * (ppore_cur - sink_prssr);
}

double PoroElastic::compute_source_derivative(double m_cur, double pnorm_cur,
                                              double &dsourcedm,
                                              double &dsourcedlambda)
{
  dsourcedm = (-source_const - sink_const) *
              compute_dpmatdjphi((m_cur / rho_s) + phi_0);
  dsourcedlambda = (source_const + sink_const);

  dsourcedm += (-source_const - sink_const) *
               compute_dppendjphi((m_cur / rho_s) + phi_0, pnorm_cur);
  dsourcedlambda += (-source_const - sink_const) *
                    compute_dppendlambda((m_cur / rho_s) + phi_0);
}

void PoroElastic::compute_mesh_volume(EquationSystems &es)
{

  // Get a constant reference to the mesh object.
  const MeshBase &mesh = es.get_mesh();

  MeshBase::const_element_iterator el = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();

  mesh_volume = 0.0;
  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem *elem = *el;

    mesh_volume += elem->volume();
  } // end of element loop
}

void PoroElastic::compute_preitr_poro(EquationSystems &es)
{
  System &system_pmat = es.get_system<System>("pMatSystem");
  System &system_ppen = es.get_system<System>("pPenSystem");
  System &system_ppore = es.get_system<System>("pPoreSystem");
  System &system_pnorm = es.get_system<System>("NonlinearElasticity");
  LinearImplicitSystem &system_darcy =
      es.get_system<LinearImplicitSystem>("darcy_velocity");
  System &system_source = es.get_system<System>("sourceSystem");

  // System & system_flow = es.get_system<System> ("flowSystem");
  System &system_flow = es.get_system<System>("porousSystem");

  NumericVector<double> &pmat_solution = *(system_pmat.solution);
  pmat_solution.close();
  pmat_solution.localize(*system_pmat.current_local_solution);

  NumericVector<double> &ppen_solution = *(system_ppen.solution);
  ppen_solution.close();
  ppen_solution.localize(*system_ppen.current_local_solution);

  NumericVector<double> &ppore_solution = *(system_ppore.solution);
  ppore_solution.close();
  ppore_solution.localize(*system_ppore.current_local_solution);

  NumericVector<double> &pnorm_solution = *(system_pnorm.solution);
  pnorm_solution.close();

  NumericVector<double> &source_solution = *(system_source.solution);
  source_solution.close();
  source_solution.localize(*system_source.current_local_solution);

  NumericVector<double> &m_solution = *(system_flow.solution);
  m_solution.close();
  m_solution.localize(*system_flow.current_local_solution);

  libMesh::MeshBase &mesh = es.get_mesh();
  libMesh::MeshBase::const_node_iterator node_it = mesh.local_nodes_begin();
  libMesh::MeshBase::const_node_iterator node_end = mesh.local_nodes_end();

  for (; node_it != node_end; ++node_it)
  {
    const libMesh::Node *nd = *node_it;
    const unsigned int dof_num_pmat =
        nd->dof_number(system_pmat.number(), 0, 0);
    const unsigned int dof_num_ppen =
        nd->dof_number(system_ppen.number(), 0, 0);
    const unsigned int dof_num_ppore =
        nd->dof_number(system_ppore.number(), 0, 0);
    // const unsigned int m_dof_num = nd->dof_number(system_pnorm.number(),
    // MESH_DIMENSION+1, 0);
    const unsigned int pnorm_dof_num =
        nd->dof_number(system_pnorm.number(), MESH_DIMENSION, 0);
    const unsigned int source_dof_num =
        nd->dof_number(system_source.number(), 0, 0);

    // const unsigned int m_dof_num = nd->dof_number(system_flow.number(), 0,
    // 0);
    const unsigned int m_dof_num =
        nd->dof_number(system_flow.number(), MESH_DIMENSION, 0);

    // double m_cur = pnorm_solution(m_dof_num);
    double m_cur = m_solution(m_dof_num);
    double pnorm_cur = -pnorm_solution(pnorm_dof_num);

    double pmat_cur = compute_pmat((m_cur / rho_s) + phi_0);
    pmat_solution.set(dof_num_pmat, pmat_cur);

    double ppen_cur = 0.0;
    if (ppore_penalty_on == 1)
    {
      ppen_cur = compute_ppen((m_cur / rho_s) + phi_0, pnorm_cur);
      ppen_solution.set(dof_num_ppen, ppen_cur);
    }

    ppore_solution.set(dof_num_ppore, pnorm_cur + pmat_cur + ppen_cur);
    source_solution.set(source_dof_num,
                        compute_source(pnorm_cur + pmat_cur + ppen_cur, 0.0));
  }
  pmat_solution.close();
  pmat_solution.localize(*system_pmat.current_local_solution);
  ppen_solution.close();
  ppen_solution.localize(*system_ppen.current_local_solution);
  ppore_solution.close();
  ppore_solution.localize(*system_ppore.current_local_solution);
  source_solution.close();
  source_solution.localize(*system_source.current_local_solution);

  system_darcy.solve();
}

// void PoroElastic::assemble_darcy (EquationSystems & es,
// const std::string & libmesh_dbg_var(system_name))
//{
//// Get a constant reference to the mesh object.
// const MeshBase & mesh = es.get_mesh();

//// The dimension that we are running
// const unsigned int dim = mesh.mesh_dimension();

//// Get a reference to the Stokes system object.
// LinearImplicitSystem & darcy_system = es.get_system<LinearImplicitSystem>
// ("darcy_velocity"); System & ppore_system = es.get_system<System>
// ("pPoreSystem"); NonlinearImplicitSystem & displacement_system =
// es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

// unsigned int u_var,v_var,w_var,ppore_var;

//// Numeric ids corresponding to each variable in the system
// u_var = darcy_system.variable_number ("w_x");
// v_var = darcy_system.variable_number ("w_y");
// #if(MESH_DIMENSION == 3)
// w_var = darcy_system.variable_number ("w_z");
// #endif
// ppore_var = ppore_system.variable_number ("pPoreVar");

//// Get the Finite Element type for "u".  Note this will be
//// the same as the type for "v".
// FEType fe_vel_type = darcy_system.variable_type(u_var);

//// Build a Finite Element object of the specified type for
//// the velocity variables.
// UniquePtr<FEBase> fe_vel  (FEBase::build(dim, fe_vel_type));

//// A Gauss quadrature rule for numerical integration.
//// Let the FEType object decide what order rule is appropriate.
// QGauss qrule (dim, fe_vel_type.default_quadrature_order());

//// Tell the finite element objects to use our quadrature rule.
// fe_vel->attach_quadrature_rule (&qrule);

// UniquePtr<FEBase> fe_face (FEBase::build(dim, fe_vel_type));
// QGauss qface(dim-1, fe_vel_type.default_quadrature_order()); // Not sure what
// the best accuracy is here

// fe_face->attach_quadrature_rule (&qface);

// const std::vector<double> & JxW = fe_vel->get_JxW();
// const std::vector<std::vector<double> > & phi = fe_vel->get_phi();
// const std::vector<std::vector<RealGradient> > & dphi = fe_vel->get_dphi();

// const DofMap & dof_map = darcy_system.get_dof_map();
// const DofMap & dof_map_dis = displacement_system.get_dof_map();
// const DofMap & dof_map_ppore = ppore_system.get_dof_map();

//// Define data structures to contain the element matrix
//// and right-hand-side vector contribution.  Following
//// basic finite element terminology we will denote these
//// "Ke" and "Fe".
// DenseMatrix<double> Ke;
// DenseVector<double> Fe;

// #if(MESH_DIMENSION == 2)
// DenseSubMatrix<double>
// Kuu(Ke), Kuv(Ke),
// Kvu(Ke), Kvv(Ke);

// DenseSubVector<double>
// Fu(Fe),
// Fv(Fe);
// #elif(MESH_DIMENSION == 3)
// DenseSubMatrix<double>
// Kuu(Ke), Kuv(Ke), Kuw(Ke),
// Kvu(Ke), Kvv(Ke), Kvw(Ke),
// Kwu(Ke), Kwv(Ke), Kww(Ke);

// DenseSubVector<double>
// Fu(Fe),
// Fv(Fe),
// Fw(Fe);
// #endif

//// This vector will hold the degree of freedom indices for
//// the element.  These define where in the global system
//// the element degrees of freedom get mapped.
// std::vector<dof_id_type> dof_indices;
// std::vector<std::vector<dof_id_type> > dof_indices_uvw(dim);

// std::vector<std::vector<dof_id_type> > dof_indices_dis(dim);

// std::vector<dof_id_type> dof_indices_stress;
// std::vector<dof_id_type> dof_indices_ppore;

//// Now we will loop over all the elements in the mesh that
//// live on the local processor. We will compute the element
//// matrix and right-hand-side contribution.  Since the mesh
//// will be refined we want to only consider the ACTIVE elements,
//// hence we use a variant of the active_elem_iterator.
// MeshBase::const_element_iterator       el     =
// mesh.active_local_elements_begin(); const MeshBase::const_element_iterator
// end_el = mesh.active_local_elements_end();

// for ( ; el != end_el; ++el)
//{
//// Store a pointer to the element we are currently
//// working on.  This allows for nicer syntax later.
// const Elem * elem = *el;

//// Get the degree of freedom indices for the
//// current element.  These define where in the global
//// matrix and right-hand-side this element will
//// contribute to.
// dof_map.dof_indices (elem, dof_indices);

// for (unsigned int var=0; var<dim; var++)
//{
// dof_map.dof_indices (elem, dof_indices_uvw[var], var);
// dof_map_dis.dof_indices (elem, dof_indices_dis[var], var);
//}

// dof_map_ppore.dof_indices (elem, dof_indices_ppore, ppore_var);

// unsigned int n_dofs, n_u_dofs, n_v_dofs, n_w_dofs;

// n_dofs   = dof_indices.size();
// n_u_dofs = dof_indices_uvw[0].size();
// n_v_dofs = dof_indices_uvw[1].size();
// #if(MESH_DIMENSION == 3)
// n_w_dofs = dof_indices_uvw[2].size();
// #endif

// fe_vel->reinit  (elem);

// Ke.resize (n_dofs, n_dofs);
// Fe.resize (n_dofs);

// Kuu.reposition (u_var*n_u_dofs, u_var*n_u_dofs, n_u_dofs, n_u_dofs);
// Kuv.reposition (u_var*n_u_dofs, v_var*n_u_dofs, n_u_dofs, n_v_dofs);
// #if(MESH_DIMENSION == 3)
// Kuw.reposition (u_var*n_u_dofs, w_var*n_u_dofs, n_u_dofs, n_w_dofs);
// #endif

// Kvu.reposition (v_var*n_v_dofs, u_var*n_v_dofs, n_v_dofs, n_u_dofs);
// Kvv.reposition (v_var*n_v_dofs, v_var*n_v_dofs, n_v_dofs, n_v_dofs);
// #if(MESH_DIMENSION == 3)
// Kvw.reposition (v_var*n_v_dofs, w_var*n_v_dofs, n_v_dofs, n_w_dofs);
// #endif

// #if(MESH_DIMENSION == 3)
//{
// Kwu.reposition (w_var*n_u_dofs, u_var*n_u_dofs, n_w_dofs, n_u_dofs);
// Kwv.reposition (w_var*n_u_dofs, v_var*n_u_dofs, n_w_dofs, n_v_dofs);
// Kww.reposition (w_var*n_u_dofs, w_var*n_u_dofs, n_w_dofs, n_w_dofs);
//}
// #endif

// Fu.reposition (u_var*n_u_dofs, n_u_dofs);
// Fv.reposition (v_var*n_u_dofs, n_v_dofs);
// #if(MESH_DIMENSION == 3)
// Fw.reposition (w_var*n_u_dofs, n_w_dofs);
// #endif

//// Now we will build the element matrix and right-hand-side.
//// Constructing the RHS requires the solution and its
//// gradient from the previous timestep.  This must be
//// calculated at each quadrature point by summing the
//// solution degree-of-freedom values by the appropriate
//// weight functions.
// for (unsigned int qp=0; qp<qrule.n_points(); qp++)
//{
// DenseMatrix<Number> FF(dim, dim);
// DenseVector<double> delta_p(dim);
////VectorValue<double> delta_p(0.0,0.0,0.0);

// for (unsigned int var_i=0; var_i<dim; var_i++)
//{
// for (unsigned int j=0; j<n_u_dofs; j++)
// delta_p(var_i) +=
// dphi[j][qp](var_i)*(ppore_system.current_solution(dof_indices_ppore[j]));

//// Row is variable u1, u2, or u3, column is x, y, or z
// for (unsigned int var_j=0; var_j<dim; var_j++)
// for (unsigned int j=0; j<n_u_dofs; j++)
// FF(var_i,var_j) +=
// dphi[j][qp](var_j)*displacement_system.current_solution(dof_indices_dis[var_i][j]);

//}

// for (unsigned int var=0; var<dim; var++)
// FF(var, var) += 1.;

// double detF = MatVecOper::detMat(FF);

// DenseMatrix<double> Ft;
// Ft.resize(dim,dim);
// MatVecOper::transposeMat(FF,Ft);

// DenseMatrix<double> C,invC;
// C.resize(dim,dim);
// invC.resize(dim,dim);
// C.add(1.0,FF);
// C.left_multiply(Ft);

// MatVecOper::inverseMat(C,invC);

// DenseVector<double> darcy_vel_vec(dim);
// invC.vector_mult(darcy_vel_vec,delta_p);
// darcy_vel_vec.scale(-detF*permeability);

////cout<<"darcy_vel_vec="<<darcy_vel_vec<<endl;

////VectorValue<double> darcy_vel_vec(0.0,0.0,0.0);

////darcy_vel_vec = -detF*permeability*invC*delta_p;

// for (unsigned int dof_i=0; dof_i<n_u_dofs; dof_i++)
//{
// Fu(dof_i) += darcy_vel_vec(0)*phi[dof_i][qp]*JxW[qp];
// Fv(dof_i) += darcy_vel_vec(1)*phi[dof_i][qp]*JxW[qp];
// #if(MESH_DIMENSION == 3)
// Fw(dof_i) += darcy_vel_vec(2)*phi[dof_i][qp]*JxW[qp];
// #endif

//// Matrix contributions for the uu and vv couplings.
// for (unsigned int dof_j=0; dof_j<n_u_dofs; dof_j++)
//{
// Kuu(dof_i,dof_j)+=phi[dof_i][qp]*phi[dof_j][qp]*JxW[qp];
// Kvv(dof_i,dof_j)+=phi[dof_i][qp]*phi[dof_j][qp]*JxW[qp];
// #if(MESH_DIMENSION == 3)
//{
// Kww(dof_i,dof_j) += phi[dof_i][qp]*phi[dof_j][qp]*JxW[qp];
//}
// #endif
// }
//}
//} // end of the quadrature point qp-loop
// darcy_system.matrix->add_matrix (Ke, dof_indices);
// darcy_system.rhs->add_vector    (Fe, dof_indices);
//} // end of element loop
//}

void PoroElastic::assemble_dmdt(
    EquationSystems &es, const std::string &libmesh_dbg_var(system_name))
{
  // Get a constant reference to the mesh object.
  const MeshBase &mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  LinearImplicitSystem &dmdt_system =
      es.get_system<LinearImplicitSystem>("dmdtSystem");
  System &darcy_system = es.get_system<System>("darcy_velocity");
  System &source_system = es.get_system<System>("sourceSystem");

  System &m_system = es.get_system<System>("NonlinearElasticity");
  const DofMap &dof_map_m = m_system.get_dof_map();
  std::vector<dof_id_type> dof_indices_m;

  System &m_old_system = es.get_system<System>("displacement_old");
  const DofMap &dof_map_m_old = m_old_system.get_dof_map();
  std::vector<dof_id_type> dof_indices_m_old;

  System &m_n1_system = es.get_system<System>("displacement_n1");
  const DofMap &dof_map_m_n1 = m_n1_system.get_dof_map();
  std::vector<dof_id_type> dof_indices_m_n1;

  unsigned int u_var;

  // Numeric ids corresponding to each variable in the system
  u_var = dmdt_system.variable_number("dmdtVar");

  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_vel_type = dmdt_system.variable_type(u_var);

  // Build a Finite Element object of the specified type for
  // the velocity variables.
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule(dim, fe_vel_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_vel_type));
  QGauss qface(dim - 1,
               fe_vel_type.default_quadrature_order()); // Not sure what the
                                                        // best accuracy is here

  fe_face->attach_quadrature_rule(&qface);

  const std::vector<double> &JxW = fe_vel->get_JxW();
  const std::vector<std::vector<double>> &phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe_vel->get_dphi();

  const DofMap &dof_map = dmdt_system.get_dof_map();
  const DofMap &dof_map_darcy = darcy_system.get_dof_map();
  const DofMap &dof_map_source = source_system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<double> Ke;
  DenseVector<double> Fe;

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_uvw(dim);
  std::vector<dof_id_type> dof_indices_source;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem *elem = *el;

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);

    for (unsigned int var = 0; var < dim; var++)
    {
      dof_map_darcy.dof_indices(elem, dof_indices_uvw[var], var);
    }

    dof_map_source.dof_indices(elem, dof_indices_source, 0);

    dof_map_m.dof_indices(elem, dof_indices_m, MESH_DIMENSION + 1);
    dof_map_m_old.dof_indices(elem, dof_indices_m_old, MESH_DIMENSION + 1);
    dof_map_m_n1.dof_indices(elem, dof_indices_m_n1, MESH_DIMENSION + 1);

    unsigned int n_dofs;

    n_dofs = dof_indices.size();

    fe_vel->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    // Now we will build the element matrix and right-hand-side.
    // Constructing the RHS requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      double m_cur = 0.0, m_old = 0.0, m_n1 = 0.0, dmdt_cur = 0.0;
      for (unsigned int j = 0; j < n_dofs; j++)
      {
        m_cur += phi[j][qp] * (m_system.current_solution(dof_indices_m[j]));
        m_old +=
            phi[j][qp] * (m_old_system.current_solution(dof_indices_m_old[j]));
        m_n1 +=
            phi[j][qp] * (m_n1_system.current_solution(dof_indices_m_n1[j]));
        dmdt_cur += phi[j][qp] * (dmdt_system.current_solution(dof_indices[j]));
      }

      if (inertia == 1)
        dmdt_cur += (m_cur - m_old) / dt;
      else if (inertia == 2)
      {
        if (time_itr == 0)
          dmdt_cur += (m_cur - m_old) / dt;
        else
          dmdt_cur += (3.0 * m_cur - 4.0 * m_old + m_n1) / (2.0 * dt);
      }

      dmdt_cur *= (1.0 / rho_s);

      for (unsigned int dof_i = 0; dof_i < n_dofs; dof_i++)
      {
        Fe(dof_i) +=
            dmdt_cur * phi[dof_i][qp] *
            JxW[qp]; //(divWS_cur+divW_cur-S_cur)*phi[dof_i][qp]*JxW[qp];

        // Matrix contributions for the uu and vv couplings.
        for (unsigned int dof_j = 0; dof_j < n_dofs; dof_j++)
        {
          Ke(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
        }
      }
    } // end of the quadrature point qp-loop
    dmdt_system.matrix->add_matrix(Ke, dof_indices);
    dmdt_system.rhs->add_vector(Fe, dof_indices);
  } // end of element loop
}

double PoroElastic::incomp_cond_poro(double NA, double m_cur, double m_old,
                                     double m_n1, double dmdt_cur)
{
  double incomp_res = 0.0;
  if (incomp_type == 0 || inertia == 0)
  {
    incomp_res = -(m_cur / rho_s) * NA;
  }
  else if (incomp_type == 1)
  {
    if (inertia == 1)
      incomp_res = (m_cur - m_old) / dt;
    else if (inertia == 2)
    {
      if (time_itr == 0)
        incomp_res = (m_cur - m_old) / dt;
      else
        incomp_res = (3.0 * m_cur - 4.0 * m_old + m_n1) / (2.0 * dt);
    }

    incomp_res *= -(1.0 / rho_s) * NA;
    incomp_res += -dmdt_cur * NA;
  }

  // cout<<"incomp_res="<<incomp_res<<endl;

  return incomp_res;
}

double PoroElastic::compute_added_mass_resid(DenseVector<double> &gradNA,
                                             DenseVector<double> &w_cur,
                                             double NA, double detF,
                                             double m_cur, double m_old,
                                             double m_n1, double source_cur)
{
  double m_resid = 0.0;
  if (inertia == 1)
    m_resid = (m_cur - m_old) / dt;
  else if (inertia == 2)
  {
    if (time_itr == 0)
      m_resid = (m_cur - m_old) / dt;
    else
      m_resid = (3.0 * m_cur - 4.0 * m_old + m_n1) / (2.0 * dt);
  }

  m_resid *= (1.0 / rho_s) * NA;

  m_resid += -MatVecOper::contractVec(w_cur, gradNA) - source_cur * NA;
  return (m_resid);
}

void PoroElastic::mresid_derivative(
    DenseMatrix<double> &F, DenseMatrix<double> &FinvTra,
    DenseMatrix<double> &invC, DenseVector<double> &gradNA,
    DenseVector<double> &gradNB, DenseVector<double> &gradPpore,
    DenseVector<double> &gradm, double NA, double NB, double detF, double m_cur,
    DenseVector<double> &dmresiddu, double &dmresiddm, double &dmresiddlambda,
    double &dincompdm, double pnorm_cur)
{
  DenseVector<double> dWdux(MESH_DIMENSION);
  DenseVector<double> dWduy(MESH_DIMENSION);
  DenseVector<double> dWduz(MESH_DIMENSION);
  DenseVector<double> dWdm(MESH_DIMENSION);
  DenseVector<double> dWdlambda(MESH_DIMENSION);

  compute_w_derivative(F, FinvTra, invC, gradNB, gradPpore, gradm, detF, m_cur,
                       NB, dWdux, dWduy, dWduz, dWdm, dWdlambda, pnorm_cur);

  dmresiddu(0) = -MatVecOper::contractVec(dWdux, gradNA);
  dmresiddu(1) = -MatVecOper::contractVec(dWduy, gradNA);
#if (MESHDIMENSION == 3)
  dmresiddu(2) = -MatVecOper::contractVec(dWduz, gradNA);
#endif

  if (inertia == 1)
    dmresiddm = (1.0 / rho_s) * (1.0 / dt) * NA * NB;
  else if (inertia == 2)
  {
    if (time_itr == 0)
      dmresiddm = (1.0 / rho_s) * (1.0 / dt) * NA * NB;
    else
      dmresiddm = (1.0 / rho_s) * (1.0 / ((2.0 / 3.0) * dt)) * NA * NB;
  }

  if (incomp_type == 0)
    dincompdm = -(1.0 / rho_s) * NA * NB;
  else if (incomp_type == 1)
    dincompdm = -dmresiddm;

  dmresiddm += -MatVecOper::contractVec(dWdm, gradNA);
  dmresiddlambda = -MatVecOper::contractVec(dWdlambda, gradNA);

  double dsourcedm = 0.0, dsourcedlambda = 0.0;
  compute_source_derivative(m_cur, pnorm_cur, dsourcedm, dsourcedlambda);

  dmresiddm += -dsourcedm * NA * NB;
  dmresiddlambda += -dsourcedlambda * NA * NB;
}

void PoroElastic::compute_w_derivative(
    DenseMatrix<double> &F, DenseMatrix<double> &FinvTra,
    DenseMatrix<double> &invC, DenseVector<double> &gradNB,
    DenseVector<double> gradPpore, DenseVector<double> gradm, double detF,
    double m_cur, double NB, DenseVector<double> &dWdux,
    DenseVector<double> &dWduy, DenseVector<double> &dWduz,
    DenseVector<double> &dWdm, DenseVector<double> &dWdlambda,
    double pnorm_cur)
{
  DenseMatrix<double> gradNBx(MESH_DIMENSION, MESH_DIMENSION);
  DenseMatrix<double> gradNBy(MESH_DIMENSION, MESH_DIMENSION);
  DenseMatrix<double> gradNBz(MESH_DIMENSION, MESH_DIMENSION);
  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
  {
    gradNBx(0, i) = gradNB(i);
    gradNBy(1, i) = gradNB(i);
    if (MESH_DIMENSION == 3)
      gradNBz(2, i) = gradNB(i);
  }

  DenseVector<double> dJdu(MESH_DIMENSION);
  dJdu(0) = detF * MatVecOper::contractMat(FinvTra, gradNBx);
  dJdu(1) = detF * MatVecOper::contractMat(FinvTra, gradNBy);
#if (MESHDIMENSION == 3)
  dJdu(2) = detF * MatVecOper::contractMat(FinvTra, gradNBz);
#endif

  DenseMatrix<double> FtgradNBx(MESH_DIMENSION, MESH_DIMENSION);
  DenseMatrix<double> FtgradNBy(MESH_DIMENSION, MESH_DIMENSION);
  DenseMatrix<double> FtgradNBz(MESH_DIMENSION, MESH_DIMENSION);

  FtgradNBx.add(1.0, gradNBx);
  FtgradNBy.add(1.0, gradNBy);
#if (MESHDIMENSION == 3)
  FtgradNBz.add(1.0, gradNBz);
#endif

  FtgradNBx.left_multiply_transpose(F);
  FtgradNBy.left_multiply_transpose(F);
#if (MESHDIMENSION == 3)
  FtgradNBz.left_multiply_transpose(F);
#endif

  DenseMatrix<double> dCinvdux(MESH_DIMENSION, MESH_DIMENSION);
  DenseMatrix<double> dCinvduy(MESH_DIMENSION, MESH_DIMENSION);
  DenseMatrix<double> dCinvduz(MESH_DIMENSION, MESH_DIMENSION);

  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
  {
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
    {
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
      {
        for (unsigned int l = 0; l < MESH_DIMENSION; l++)
        {
          dCinvdux(i, j) +=
              -0.5 * (invC(i, k) * invC(j, l) + invC(i, l) * invC(j, k)) *
              FtgradNBx(k, l);
          dCinvduy(i, j) +=
              -0.5 * (invC(i, k) * invC(j, l) + invC(i, l) * invC(j, k)) *
              FtgradNBy(k, l);
#if (MESHDIMENSION == 3)
          dCinvduz(i, j) +=
              -0.5 * (invC(i, k) * invC(j, l) + invC(i, l) * invC(j, k)) *
              FtgradNBz(k, l);
#endif
        }
      }
    }
  }

  DenseMatrix<double> dWldux(MESH_DIMENSION, MESH_DIMENSION);
  DenseMatrix<double> dWlduy(MESH_DIMENSION, MESH_DIMENSION);
  DenseMatrix<double> dWlduz(MESH_DIMENSION, MESH_DIMENSION);

  dWldux.add(1.0, II);
  dWldux.left_multiply(invC);
  dWldux.scale(-dJdu(0) * permeability);
  dWldux.add(-detF * permeability, dCinvdux);

  dWlduy.add(1.0, II);
  dWlduy.left_multiply(invC);
  dWlduy.scale(-dJdu(1) * permeability);
  dWlduy.add(-detF * permeability, dCinvduy);

#if (MESHDIMENSION == 3)
  dWlduz.add(1.0, II);
  dWlduz.left_multiply(invC);
  dWlduz.scale(-dJdu(2) * permeability);
  dWlduz.add(-detF * permeability, dCinvduz);
#endif

  dWldux.vector_mult(dWdux, gradPpore);
  dWlduy.vector_mult(dWduy, gradPpore);
#if (MESHDIMENSION == 3)
  dWlduz.vector_mult(dWduz, gradPpore);
#endif

  DenseVector<double> dWdm1(MESH_DIMENSION);
  invC.vector_mult(dWdm1, gradNB);
  dWdm1.scale(-detF * permeability *
              (compute_dpmatdjphi((m_cur / rho_s) + phi_0) +
               compute_dppendjphi((m_cur / rho_s) + phi_0, pnorm_cur)));

  invC.vector_mult(dWdm, gradm);
  dWdm.scale(-detF * permeability *
             compute_d2pmatdjphi2((m_cur / rho_s) + phi_0) * NB);

  dWdm.add(1.0, dWdm1);

  invC.vector_mult(dWdlambda, gradNB);
  dWdlambda.scale(detF * permeability *
                  (1.0 - compute_dppendlambda((m_cur / rho_s) + phi_0)));
}

// void PoroElastic::compute_mmono(EquationSystems & es)
//{
// const MeshBase & mesh = es.get_mesh();
// const unsigned int dim = mesh.mesh_dimension();

// NonlinearImplicitSystem & system =
// es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

// const unsigned int u_var = system.variable_number ("u");

// const DofMap & dof_map = system.get_dof_map();
// FEType fe_type = dof_map.variable_type(u_var);
// UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
// QGauss qrule (dim, fe_type.default_quadrature_order());
// fe->attach_quadrature_rule (&qrule);

// const std::vector<Real> & JxW = fe->get_JxW();
// const std::vector<std::vector<Real> > & phi = fe->get_phi();
// const std::vector<std::vector<RealGradient> > & dphi = fe->get_dphi();

//// Also, get a reference to the ExplicitSystem
// ExplicitSystem & mmono_system = es.get_system<ExplicitSystem>("mMonoSystem");
// const DofMap & mmono_dof_map = mmono_system.get_dof_map();
// unsigned int mmono_var;
// mmono_var = mmono_system.variable_number ("mMonoVar");

//// Storage for the stress dof indices on each element
// std::vector<dof_id_type> dof_indices_m;
// std::vector<dof_id_type> dof_indices_mmono;

// MeshBase::const_element_iterator       el     =
// mesh.active_local_elements_begin(); const MeshBase::const_element_iterator
// end_el = mesh.active_local_elements_end();

// for ( ; el != end_el; ++el)
//{
// const Elem * elem = *el;

// dof_map.dof_indices (elem, dof_indices_m, dim+1);

// const unsigned int n_var_dofs = dof_indices_m.size();

// fe->reinit (elem);

// Number mmono_cur=0.0;

// for (unsigned int qp=0; qp<qrule.n_points(); qp++)
//{
// double mmono_ = 0.0;

// for (unsigned int j=0; j<n_var_dofs; j++)
// mmono_ += phi[j][qp] * system.current_solution(dof_indices_m[j]);

// mmono_cur += mmono_*JxW[qp];
////cout<<"J="<<F.det()<<endl;
//}

// mmono_cur *= (1./elem->volume());

// const int dof_index = elem->dof_number(mmono_system.number(),0,0);

// mmono_system.solution->set(dof_index, mmono_cur);
//}

//// Should call close and update when we set vector entries directly
// mmono_system.solution->close();
// mmono_system.update();
//}

void PoroElastic::assemble_flow(
    EquationSystems &es, const std::string &libmesh_dbg_var(system_name))
{
  // Get a constant reference to the mesh object.
  const MeshBase &mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  LinearImplicitSystem &flow_system =
      es.get_system<LinearImplicitSystem>("flowSystem");

  unsigned int u_var;

  System &K_system = es.get_system<System>("KSystem");
  const DofMap &K_dof_map = K_system.get_dof_map();

  System &system_source = es.get_system<System>("sourceSystem");
  unsigned int system_source_num = system_source.number();

  NumericVector<double> &K_solution = *(K_system.solution);
  K_solution.close();
  K_solution.localize(*K_system.current_local_solution);

  // Numeric ids corresponding to each variable in the system
  u_var = flow_system.variable_number("mVar");
  FEType fe_vel_type = flow_system.variable_type(u_var);
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));
  QGauss qrule(dim, CONSTANT);

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_vel_type));
  QGauss qface(dim - 1,
               CONSTANT); // Not sure what the
                          // best accuracy is here

  fe_face->attach_quadrature_rule(&qface);

  const std::vector<double> &JxW = fe_vel->get_JxW();
  const std::vector<std::vector<double>> &phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe_vel->get_dphi();

  const DofMap &dof_map = flow_system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<double> Ke;
  DenseVector<double> Fe;

  DenseSubMatrix<double> Kmm(Ke);

  DenseSubVector<double> Fm(Fe);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_m;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  double dmdx_left = 0.0, dmdx_right = 0.0;

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem *elem = *el;

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);
    dof_map.dof_indices(elem, dof_indices_m, 0);

    unsigned int n_dofs, n_u_dofs;

    n_dofs = dof_indices.size();
    n_u_dofs = dof_indices_m.size();

    fe_vel->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    Kmm.reposition(0, 0, n_u_dofs, n_u_dofs);

    Fm.reposition(0, n_u_dofs);

    const int dof_index_source = elem->dof_number(system_source_num, 0, 0);
    double source_cur = system_source.current_solution(dof_index_source);

#if (MESH_DIMENSION == 2)
    const int dof_index_00 = elem->dof_number(K_system.number(), 0, 0);
    double K_00 = K_solution(dof_index_00);
    const int dof_index_01 = elem->dof_number(K_system.number(), 1, 0);
    double K_01 = K_solution(dof_index_01);

    const int dof_index_10 = elem->dof_number(K_system.number(), 2, 0);
    double K_10 = K_solution(dof_index_10);
    const int dof_index_11 = elem->dof_number(K_system.number(), 3, 0);
    double K_11 = K_solution(dof_index_11);

#elif (MESH_DIMENSION == 3)
    const int dof_index_00 = elem->dof_number(K_system.number(), 0, 0);
    double K_00 = K_solution(dof_index_00);
    const int dof_index_01 = elem->dof_number(K_system.number(), 1, 0);
    double K_01 = K_solution(dof_index_01);
    const int dof_index_02 = elem->dof_number(K_system.number(), 2, 0);
    double K_02 = K_solution(dof_index_02);

    const int dof_index_10 = elem->dof_number(K_system.number(), 3, 0);
    double K_10 = K_solution(dof_index_10);
    const int dof_index_11 = elem->dof_number(K_system.number(), 4, 0);
    double K_11 = K_solution(dof_index_11);
    const int dof_index_12 = elem->dof_number(K_system.number(), 5, 0);
    double K_12 = K_solution(dof_index_12);

    const int dof_index_20 = elem->dof_number(K_system.number(), 6, 0);
    double K_20 = K_solution(dof_index_20);
    const int dof_index_21 = elem->dof_number(K_system.number(), 7, 0);
    double K_21 = K_solution(dof_index_21);
    const int dof_index_22 = elem->dof_number(K_system.number(), 8, 0);
    double K_22 = K_solution(dof_index_22);
#endif

    // Now we will build the element matrix and right-hand-side.
    // Constructing the RHS requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {

      double m_old = 0.0;
      for (unsigned int j = 0; j < n_u_dofs; j++)
      {
        m_old += phi[j][qp] * (flow_system.current_solution(dof_indices_m[j]));
      }

      GeomPar::compute_geoPar(es, elem, qp, phi, dphi);

      DenseVector<double> divJCinvdelW(dim);
      DenseVector<double> JFinvmMuKdivJCinvdelW(MESH_DIMENSION);
      if (brinkman == 1)
      {
        LinearImplicitSystem &system_delw =
            es.get_system<LinearImplicitSystem>("systemDelW");
        const DofMap &dof_map_delw = system_delw.get_dof_map();
        std::vector<std::vector<dof_id_type>> dof_indices_delw(dim * dim);

        for (unsigned int var = 0; var < dim * dim; var++)
        {
          dof_map_delw.dof_indices(elem, dof_indices_delw[var], var);
        }

        for (unsigned int var_i = 0; var_i < dim; var_i++)
        {

          // Row is variable u1, u2, or u3, column is x, y, or z
          for (unsigned int var_j = 0; var_j < dim; var_j++)
            for (unsigned int j = 0; j < n_u_dofs; j++)
            {
              divJCinvdelW(var_i) +=
                  dphi[j][qp](var_j) *
                  system_delw.current_solution(
                      dof_indices_delw[dim * var_i + var_j][j]);
            }
        }

        GeomPar::FInv.vector_mult(JFinvmMuKdivJCinvdelW, divJCinvdelW);
        JFinvmMuKdivJCinvdelW.scale(viscocity * permeability);
      }

      DenseMatrix<double> K_perm(MESH_DIMENSION, MESH_DIMENSION);
      /* K_perm(0, 0) = K_00;
      K_perm(0, 1) = K_01;
      K_perm(0, 2) = K_02;

      K_perm(1, 0) = K_10;
      K_perm(1, 1) = K_11;
      K_perm(1, 2) = K_12;

      K_perm(2, 0) = K_20;
      K_perm(2, 1) = K_21;
      K_perm(2, 2) = K_22; */

      K_perm(0, 0) = K_00;
      K_perm(0, 1) = K_01;
#if (MESH_DIMENSION == 3)
      K_perm(0, 2) = K_02;
#endif

      K_perm(1, 0) = K_10;
      K_perm(1, 1) = K_11;
#if (MESH_DIMENSION == 3)
      K_perm(1, 2) = K_12;

      K_perm(2, 0) = K_20;
      K_perm(2, 1) = K_21;
      K_perm(2, 2) = K_22;
#endif

      K_perm.left_multiply(GeomPar::FInv);
      K_perm.right_multiply(GeomPar::FInvTra);
      K_perm.scale(GeomPar::detF * permeability);

      DenseVector<double> kdivXP(MESH_DIMENSION);
      K_perm.vector_mult(kdivXP, GeomPar::grad_pre);

      for (unsigned int dof_i = 0; dof_i < n_u_dofs; dof_i++)
      {
        DenseVector<double> gradNA;
        gradNA.resize(MESH_DIMENSION);
        for (unsigned int i = 0; i < dim; i++)
        {
          gradNA(i) = dphi[dof_i][qp](i);
        }

        Fm(dof_i) += ((m_old / dt) * phi[dof_i][qp]) * JxW[qp] +
                     MatVecOper::contractVec(kdivXP, gradNA) * JxW[qp];

        Fm(dof_i) += source_cur * phi[dof_i][qp] * JxW[qp];

        if (brinkman == 1)
          Fm(dof_i) +=
              MatVecOper::contractVec(JFinvmMuKdivJCinvdelW, gradNA) * JxW[qp];

        // Matrix contributions for the uu and vv couplings.
        for (unsigned int dof_j = 0; dof_j < n_u_dofs; dof_j++)
        {
          DenseVector<double> gradNB;
          gradNB.resize(MESH_DIMENSION);
          for (unsigned int i = 0; i < dim; i++)
          {
            gradNB(i) = dphi[dof_j][qp](i);
          }

          DenseVector<double> kdivXNB(MESH_DIMENSION);
          K_perm.vector_mult(kdivXNB, gradNB);

          Kmm(dof_i, dof_j) +=
              (1.0 / dt) * phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp] +
              kappa_0 * MatVecOper::contractVec(kdivXNB, gradNA) * JxW[qp];
        }
      }
    } // end of the quadrature point qp-loop

    flow_system.matrix->add_matrix(Ke, dof_indices);
    flow_system.rhs->add_vector(Fe, dof_indices);
  } // end of element loop
}

void PoroElastic::assemble_porous(
    EquationSystems &es, const std::string &libmesh_dbg_var(system_name))
{
  // Get a constant reference to the mesh object.
  const MeshBase &mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  LinearImplicitSystem &darcy_system =
      es.get_system<LinearImplicitSystem>("porousSystem");
  NonlinearImplicitSystem &displacement_system =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");

  unsigned int u_var, v_var, w_var, m_var;

  // Numeric ids corresponding to each variable in the system
  u_var = darcy_system.variable_number("wxVar");
  v_var = darcy_system.variable_number("wyVar");
#if (MESH_DIMENSION == 3)
  w_var = darcy_system.variable_number("wzVar");
#endif
  m_var = darcy_system.variable_number("mAddedVar");

  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_vel_type = darcy_system.variable_type(u_var);

  // Build a Finite Element object of the specified type for
  // the velocity variables.
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));

  // A Gauss quadrature rule for numerical integration.
  // Let the FEType object decide what order rule is appropriate.
  QGauss qrule(dim, fe_vel_type.default_quadrature_order());

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_vel_type));
  QGauss qface(dim - 1,
               fe_vel_type.default_quadrature_order()); // Not sure what the
                                                        // best accuracy is here

  fe_face->attach_quadrature_rule(&qface);

  const std::vector<double> &JxW = fe_vel->get_JxW();
  const std::vector<std::vector<double>> &phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe_vel->get_dphi();

  const DofMap &dof_map = darcy_system.get_dof_map();
  const DofMap &dof_map_dis = displacement_system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<double> Ke;
  DenseVector<double> Fe;

#if (MESH_DIMENSION == 2)
  DenseSubMatrix<double> Kuu(Ke), Kuv(Ke), Kum(Ke), Kvu(Ke), Kvv(Ke), Kvm(Ke),
      Kmu(Ke), Kmv(Ke), Kmm(Ke);

  DenseSubVector<double> Fu(Fe), Fv(Fe), Fm(Fe);
#elif (MESH_DIMENSION == 3)
  DenseSubMatrix<double> Kuu(Ke), Kuv(Ke), Kuw(Ke), Kum(Ke), Kvu(Ke), Kvv(Ke),
      Kvw(Ke), Kvm(Ke), Kwu(Ke), Kwv(Ke), Kww(Ke), Kwm(Ke), Kmu(Ke), Kmv(Ke),
      Kmw(Ke), Kmm(Ke);

  DenseSubVector<double> Fu(Fe), Fv(Fe), Fw(Fe), Fm(Fe);
#endif

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_uvw(dim + 1);

  std::vector<std::vector<dof_id_type>> dof_indices_dis(dim + 1);

  std::vector<dof_id_type> dof_indices_stress;
  std::vector<dof_id_type> dof_indices_ppore;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem *elem = *el;

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);

    for (unsigned int var = 0; var < dim + 1; var++)
    {
      dof_map.dof_indices(elem, dof_indices_uvw[var], var);
      dof_map_dis.dof_indices(elem, dof_indices_dis[var], var);
    }

    unsigned int n_dofs, n_u_dofs, n_v_dofs, n_w_dofs;

    n_dofs = dof_indices.size();
    n_u_dofs = dof_indices_uvw[0].size();
    n_v_dofs = dof_indices_uvw[1].size();
#if (MESH_DIMENSION == 3)
    n_w_dofs = dof_indices_uvw[2].size();
#endif

    fe_vel->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    Kuu.reposition(u_var * n_u_dofs, u_var * n_u_dofs, n_u_dofs, n_u_dofs);
    Kuv.reposition(u_var * n_u_dofs, v_var * n_u_dofs, n_u_dofs, n_v_dofs);
#if (MESH_DIMENSION == 3)
    Kuw.reposition(u_var * n_u_dofs, w_var * n_u_dofs, n_u_dofs, n_w_dofs);
#endif
    Kum.reposition(u_var * n_u_dofs, m_var * n_u_dofs, n_u_dofs, n_u_dofs);

    Kvu.reposition(v_var * n_v_dofs, u_var * n_v_dofs, n_v_dofs, n_u_dofs);
    Kvv.reposition(v_var * n_v_dofs, v_var * n_v_dofs, n_v_dofs, n_v_dofs);
#if (MESH_DIMENSION == 3)
    Kvw.reposition(v_var * n_v_dofs, w_var * n_v_dofs, n_v_dofs, n_w_dofs);
#endif
    Kvm.reposition(v_var * n_v_dofs, m_var * n_v_dofs, n_v_dofs, n_v_dofs);

#if (MESH_DIMENSION == 3)
    {
      Kwu.reposition(w_var * n_u_dofs, u_var * n_u_dofs, n_w_dofs, n_u_dofs);
      Kwv.reposition(w_var * n_u_dofs, v_var * n_u_dofs, n_w_dofs, n_v_dofs);
      Kww.reposition(w_var * n_u_dofs, w_var * n_u_dofs, n_w_dofs, n_w_dofs);
      Kwm.reposition(w_var * n_u_dofs, m_var * n_u_dofs, n_w_dofs, n_w_dofs);
    }
#endif

    Kmu.reposition(m_var * n_u_dofs, u_var * n_u_dofs, n_u_dofs, n_u_dofs);
    Kmv.reposition(m_var * n_u_dofs, v_var * n_u_dofs, n_u_dofs, n_v_dofs);
#if (MESH_DIMENSION == 3)
    Kmw.reposition(m_var * n_u_dofs, w_var * n_u_dofs, n_u_dofs, n_w_dofs);
#endif
    Kmm.reposition(m_var * n_u_dofs, m_var * n_u_dofs, n_u_dofs, n_u_dofs);

    Fu.reposition(u_var * n_u_dofs, n_u_dofs);
    Fv.reposition(v_var * n_u_dofs, n_v_dofs);
#if (MESH_DIMENSION == 3)
    Fw.reposition(w_var * n_u_dofs, n_w_dofs);
#endif
    Fm.reposition(m_var * n_u_dofs, n_u_dofs);

    // Now we will build the element matrix and right-hand-side.
    // Constructing the RHS requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      DenseMatrix<Number> FF(dim, dim);
      DenseVector<double> delta_p(dim);
      // VectorValue<double> delta_p(0.0,0.0,0.0);

      for (unsigned int var_i = 0; var_i < dim; var_i++)
      {
        for (unsigned int j = 0; j < n_u_dofs; j++)
          delta_p(var_i) +=
              dphi[j][qp](var_i) *
              (displacement_system.current_solution(dof_indices_dis[dim][j]));

        // Row is variable u1, u2, or u3, column is x, y, or z
        for (unsigned int var_j = 0; var_j < dim; var_j++)
          for (unsigned int j = 0; j < n_u_dofs; j++)
            FF(var_i, var_j) +=
                dphi[j][qp](var_j) *
                displacement_system.current_solution(dof_indices_dis[var_i][j]);
      }

      double m_old = 0.0;
      for (unsigned int j = 0; j < n_u_dofs; j++)
        m_old += phi[j][qp] *
                 (darcy_system.current_solution(dof_indices_uvw[dim][j]));

      // cout<<"m_old="<<m_old<<endl;

      for (unsigned int var = 0; var < dim; var++)
        FF(var, var) += 1.;

      double detF = MatVecOper::detMat(FF);

      DenseMatrix<double> Ft;
      Ft.resize(dim, dim);
      MatVecOper::transposeMat(FF, Ft);

      DenseMatrix<double> invF;
      invF.resize(dim, dim);
      MatVecOper::inverseMat(FF, invF);

      DenseMatrix<double> invFt;
      invFt.resize(dim, dim);
      MatVecOper::transposeMat(invF, invFt);

      DenseMatrix<double> C, invC;
      C.resize(dim, dim);
      invC.resize(dim, dim);
      C.add(1.0, FF);
      C.left_multiply(Ft);

      MatVecOper::inverseMat(C, invC);

      DenseVector<double> darcy_vel_vec(dim);
      invC.vector_mult(darcy_vel_vec, delta_p);
      darcy_vel_vec.scale(-detF * permeability);

      DenseVector<double> KFinvTdelP(dim);
      invFt.vector_mult(KFinvTdelP, delta_p);
      KFinvTdelP.scale(-permeability);

      // cout<<"darcy_vel_vec="<<darcy_vel_vec<<endl;

      // VectorValue<double> darcy_vel_vec(0.0,0.0,0.0);

      // darcy_vel_vec = -detF*permeability*invC*delta_p;

      for (unsigned int dof_i = 0; dof_i < n_u_dofs; dof_i++)
      {
        DenseVector<double> gradNA;
        gradNA.resize(MESH_DIMENSION);
        for (unsigned int i = 0; i < dim; i++)
        {
          gradNA(i) = dphi[dof_i][qp](i);
        }

        Fu(dof_i) += -KFinvTdelP(0) * phi[dof_i][qp] * JxW[qp];
        Fv(dof_i) += -KFinvTdelP(1) * phi[dof_i][qp] * JxW[qp];
#if (MESH_DIMENSION == 3)
        Fw(dof_i) += -KFinvTdelP(2) * phi[dof_i][qp] * JxW[qp];
#endif

        // Fm(dof_i) += (1.0/rho_s)*(m_old/dt)*phi[dof_i][qp]*JxW[qp];

        // Matrix contributions for the uu and vv couplings.
        for (unsigned int dof_j = 0; dof_j < n_u_dofs; dof_j++)
        {
          DenseVector<double> NBx;
          NBx.resize(MESH_DIMENSION);
          NBx(0) = phi[dof_j][qp];

          DenseVector<double> NBy;
          NBy.resize(MESH_DIMENSION);
          NBy(1) = phi[dof_j][qp];

#if (MESH_DIMENSION == 3)
          DenseVector<double> NBz;
          NBz.resize(MESH_DIMENSION);
          NBz(2) = phi[dof_j][qp];
#endif

          DenseVector<double> JFinvNBx(dim);
          invF.vector_mult(JFinvNBx, NBx);
          JFinvNBx.scale(detF);

          DenseVector<double> JFinvNBy(dim);
          invF.vector_mult(JFinvNBy, NBy);
          JFinvNBy.scale(detF);

#if (MESH_DIMENSION == 3)
          DenseVector<double> JFinvNBz(dim);
          invF.vector_mult(JFinvNBz, NBz);
          JFinvNBz.scale(detF);
#endif

          Kmu(dof_i, dof_j) +=
              -MatVecOper::contractVec(JFinvNBx, gradNA) * JxW[qp];
          Kmv(dof_i, dof_j) +=
              -MatVecOper::contractVec(JFinvNBy, gradNA) * JxW[qp];
#if (MESH_DIMENSION == 3)
          Kmw(dof_i, dof_j) +=
              -MatVecOper::contractVec(JFinvNBz, gradNA) * JxW[qp];
#endif

          // Kmm(dof_i,dof_j) +=
          // (1.0/rho_s)*(1.0/dt)*phi[dof_i][qp]*phi[dof_j][qp]*JxW[qp];

          DenseVector<double> gradNB;
          gradNB.resize(MESH_DIMENSION);
          for (unsigned int i = 0; i < dim; i++)
          {
            gradNB(i) = dphi[dof_j][qp](i);
          }

          DenseVector<double> KFinvTgradNB(dim);
          invFt.vector_mult(KFinvTgradNB, gradNB);
          KFinvTgradNB.scale(permeability);

          Kuu(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * detF * JxW[qp];
          Kvv(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * detF * JxW[qp];
#if (MESH_DIMENSION == 3)
          {
            Kww(dof_i, dof_j) +=
                phi[dof_i][qp] * phi[dof_j][qp] * detF * JxW[qp];
          }
#endif

          Kum(dof_i, dof_j) +=
              kappa_0 * KFinvTgradNB(0) * phi[dof_i][qp] * detF * JxW[qp];
          Kvm(dof_i, dof_j) +=
              kappa_0 * KFinvTgradNB(1) * phi[dof_i][qp] * detF * JxW[qp];
#if (MESH_DIMENSION == 3)
          {
            Kwm(dof_i, dof_j) +=
                kappa_0 * KFinvTgradNB(2) * phi[dof_i][qp] * detF * JxW[qp];
          }
#endif

          DenseMatrix<double> JinvCgradNBx, JinvCgradNBy, JinvCgradNBz;
          JinvCgradNBx.resize(dim, dim);
          JinvCgradNBy.resize(dim, dim);
          JinvCgradNBz.resize(dim, dim);
          for (unsigned int i = 0; i < dim; i++)
          {
            JinvCgradNBx(0, i) = dphi[dof_j][qp](i);
            JinvCgradNBy(1, i) = dphi[dof_j][qp](i);
            if (dim == 3)
              JinvCgradNBz(2, i) = dphi[dof_j][qp](i);
          }

          JinvCgradNBx.left_multiply(invC);
          JinvCgradNBy.left_multiply(invC);
          if (dim == 3)
            JinvCgradNBz.left_multiply(invC);

          JinvCgradNBx.scale(detF);
          JinvCgradNBy.scale(detF);
          if (dim == 3)
            JinvCgradNBz.scale(detF);

          DenseVector<double> JinvCgradNBxVec(dim);
          DenseVector<double> JinvCgradNByVec(dim);
          DenseVector<double> JinvCgradNBzVec(dim);
          for (unsigned int i = 0; i < dim; i++)
          {
            JinvCgradNBxVec(i) = JinvCgradNBx(0, i);
            JinvCgradNByVec(i) = JinvCgradNBy(1, i);
            if (dim == 3)
              JinvCgradNBzVec(i) = JinvCgradNBz(2, i);
          }

          if (viscous_on == 1)
          {
            Kuu(dof_i, dof_j) +=
                beta_visc * MatVecOper::contractVec(JinvCgradNBxVec, gradNA) *
                JxW[qp];
            Kvv(dof_i, dof_j) +=
                beta_visc * MatVecOper::contractVec(JinvCgradNByVec, gradNA) *
                JxW[qp];
#if (MESH_DIMENSION == 3)
            Kww(dof_i, dof_j) +=
                beta_visc * MatVecOper::contractVec(JinvCgradNBzVec, gradNA) *
                JxW[qp];
#endif
          }
        }
      }

      // if(viscous_on == 1)
      {
        for (unsigned int side = 0; side < elem->n_sides(); side++)
          if (elem->neighbor_ptr(side) == libmesh_nullptr)
          {

            const std::vector<std::vector<Real>> &phi_face = fe_face->get_phi();
            const std::vector<Real> &JxW_face = fe_face->get_JxW();

            const std::vector<Point> &normal_face = fe_face->get_normals();
            const std::vector<Point> &Xref = fe_face->get_xyz();

            fe_face->reinit(elem, side);

            for (unsigned int qp = 0; qp < qface.n_points(); qp++)
            {
              vector<boundary_id_type> bc_id_vec;
              mesh.boundary_info->boundary_ids(elem, side, bc_id_vec);
              short int bc_id =
                  bc_id_vec[0]; // mesh.boundary_info->boundary_id(elem, side);

              // if( bc_id == 1 || bc_id == 2 || bc_id == 3 || bc_id == 4 ||
              // bc_id == 5) // This is the x=+ boundary
              //{
              // for (std::size_t i=0; i<phi_face.size(); i++)
              //{
              // for (std::size_t j=0; j<phi_face.size(); j++)
              //{
              // Kuu(i,j)
              // += 1.0e12*phi_face[i][qp]*phi_face[j][qp]*JxW_face[qp];
              // Kvv(i,j)
              // += 1.0e12*phi_face[i][qp]*phi_face[j][qp]*JxW_face[qp];
              // #if(MESH_DIMENSION == 3)
              // Kww(i,j)
              // += 1.0e12*phi_face[i][qp]*phi_face[j][qp]*JxW_face[qp]; #endif
              //}

              //}
              //}

              if (bc_id == 4096 || bc_id == 4097) // This is the x=+ boundary
              {
                for (std::size_t i = 0; i < phi_face.size(); i++)
                {
                  if (bc_id == 4096)
                    Fm(i) += 0.1 * 1.0e6 * phi_face[i][qp] * JxW_face[qp];
                  for (std::size_t j = 0; j < phi_face.size(); j++)
                  {
                    Kmm(i, j) += 1.0e6 * phi_face[i][qp] * phi_face[j][qp] *
                                 JxW_face[qp];
                  }
                }
              }
            }
          }
      }
    } // end of the quadrature point qp-loop
    darcy_system.matrix->add_matrix(Ke, dof_indices);
    darcy_system.rhs->add_vector(Fe, dof_indices);
  } // end of element loop
}

void PoroElastic::compute_mmono(EquationSystems &es)
{
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem &system =
      es.get_system<LinearImplicitSystem>("flowSystem");

  const unsigned int u_var = system.variable_number("mVar");

  const DofMap &dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, CONSTANT);
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<Real>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  // Also, get a reference to the ExplicitSystem
  ExplicitSystem &mmono_system = es.get_system<ExplicitSystem>("mMonoSystem");
  const DofMap &mmono_dof_map = mmono_system.get_dof_map();
  unsigned int mmono_var;
  mmono_var = mmono_system.variable_number("mMonoVar");

  // Storage for the stress dof indices on each element
  std::vector<dof_id_type> dof_indices_m;
  std::vector<dof_id_type> dof_indices_mmono;

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    dof_map.dof_indices(elem, dof_indices_m, 0);

    const unsigned int n_var_dofs = dof_indices_m.size();

    fe->reinit(elem);

    Number mmono_cur = 0.0;
    for (unsigned int j = 0; j < n_var_dofs; j++)
      mmono_cur += phi[j][0] * system.current_solution(dof_indices_m[j]);

    const int dof_index = elem->dof_number(mmono_system.number(), 0, 0);

    mmono_system.solution->set(dof_index, mmono_cur);
  }

  // Should call close and update when we set vector entries directly
  mmono_system.solution->close();
  mmono_system.update();
}

void PoroElastic::initialise_poroelastic(EquationSystems &es)
{
  compute_mesh_volume(es);
  if (read_permeability == 0)
    initialise_K(es);
  else
    read_porous_data(es);
}

void PoroElastic::update_poroelastic(EquationSystems &es)
{
  dt = InputParam::dt;
  ttime = InputParam::ttime;

  solve_flow_system(es);
  compute_mmono(es);
  update_ppore(es);
  solve_darcy(es);
  solve_mexp_system(es);

  time_itr++;
}

void PoroElastic::initialise_K(EquationSystems &es)
{

  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  System &K_system = es.get_system<System>("KSystem");
  const DofMap &K_dof_map = K_system.get_dof_map();

  FEType fe_type = K_dof_map.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, CONSTANT);
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Point> &Xref = fe->get_xyz();

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  // int count_Jn = 0;

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    fe->reinit(elem);

    // double x_elem = Xref[0].operator()(0);
    // double y_elem = Xref[0].operator()(1);
    // double z_elem = Xref[0].operator()(2);
    //
    // double distmin = 10.0, dum;
    // int vess_index = 0;
    // for (int im = 0; im < 17857; im++) {
    //   dum = fabs(FindDistanceToSegment(
    //       vasc_nodes_x(int(edges_1(im))), vasc_nodes_y(int(edges_1(im))),
    //       vasc_nodes_z(int(edges_1(im))), vasc_nodes_x(int(edges_2(im))),
    //       vasc_nodes_y(int(edges_2(im))), vasc_nodes_z(int(edges_2(im))),
    //       x_elem, y_elem, z_elem));
    //   if (dum < distmin) {
    //     distmin = dum;
    //     vess_index = im;
    //   }
    // }
    //
    double K00 = 0.0, K01 = 0.0, K02 = 0.0, K10 = 0.0, K11 = 0.0, K12 = 0.0,
           K20 = 0.0, K21 = 0.0, K22 = 0.0;
    // // compute_K(vess_index,K00,K01,K02,K10,K11,K12,K20,K21,K22);

    K00 = 1.0, K01 = 0.0, K02 = 0.0, K10 = 0.0, K11 = 1.0, K12 = 0.0, K20 = 0.0,
    K21 = 0.0, K22 = 1.0;

    // K00 = f0_e(0)*f0_e(0), K01 = f0_e(0)*f0_e(1), K02 = f0_e(0)*f0_e(2);
    // K10 = f0_e(1)*f0_e(0), K11 = f0_e(1)*f0_e(1), K12 = f0_e(1)*f0_e(2);
    // K20 = f0_e(2)*f0_e(0), K21 = f0_e(2)*f0_e(1), K22 = f0_e(2)*f0_e(2);

#if (MESH_DIMENSION == 2)
    const int dof_index_00 = elem->dof_number(K_system.number(), 0, 0);
    K_system.solution->set(dof_index_00, K00);
    const int dof_index_01 = elem->dof_number(K_system.number(), 1, 0);
    K_system.solution->set(dof_index_01, K01);

    const int dof_index_10 = elem->dof_number(K_system.number(), 2, 0);
    K_system.solution->set(dof_index_10, K10);
    const int dof_index_11 = elem->dof_number(K_system.number(), 3, 0);
    K_system.solution->set(dof_index_11, K11);

#elif (MESH_DIMENSION == 3)
    const int dof_index_00 = elem->dof_number(K_system.number(), 0, 0);
    K_system.solution->set(dof_index_00, K00);
    const int dof_index_01 = elem->dof_number(K_system.number(), 1, 0);
    K_system.solution->set(dof_index_01, K01);

    const int dof_index_02 = elem->dof_number(K_system.number(), 2, 0);
    K_system.solution->set(dof_index_02, K02);

    const int dof_index_10 = elem->dof_number(K_system.number(), 3, 0);
    K_system.solution->set(dof_index_10, K10);
    const int dof_index_11 = elem->dof_number(K_system.number(), 4, 0);
    K_system.solution->set(dof_index_11, K11);

    const int dof_index_12 = elem->dof_number(K_system.number(), 5, 0);
    K_system.solution->set(dof_index_12, K12);
    const int dof_index_20 = elem->dof_number(K_system.number(), 6, 0);
    K_system.solution->set(dof_index_20, K20);
    const int dof_index_21 = elem->dof_number(K_system.number(), 7, 0);
    K_system.solution->set(dof_index_21, K21);
    const int dof_index_22 = elem->dof_number(K_system.number(), 8, 0);
    K_system.solution->set(dof_index_22, K22);
#endif
  }

  K_system.solution->close();
  K_system.solution->localize(*K_system.current_local_solution);
}

void PoroElastic::compute_K(int iv, double &K00, double &K01, double &K02,
                            double &K10, double &K11, double &K12, double &K20,
                            double &K21, double &K22)
{
  double dist_vess = sqrt(
      pow(vasc_nodes_x(int(edges_1(iv))) - vasc_nodes_x(int(edges_2(iv))), 2) +
      pow(vasc_nodes_y(int(edges_1(iv))) - vasc_nodes_y(int(edges_2(iv))), 2) +
      pow(vasc_nodes_z(int(edges_1(iv))) - vasc_nodes_z(int(edges_2(iv))), 2));
  double dx_vess =
      vasc_nodes_x(int(edges_1(iv))) - vasc_nodes_x(int(edges_2(iv)));
  double dy_vess =
      vasc_nodes_y(int(edges_1(iv))) - vasc_nodes_y(int(edges_2(iv)));
  double dz_vess =
      vasc_nodes_z(int(edges_1(iv))) - vasc_nodes_z(int(edges_2(iv)));

  double cos_a = dx_vess / dist_vess;
  double cos_b = dy_vess / dist_vess;
  double cos_c = dz_vess / dist_vess;

  K00 = cos_a * cos_a;
  K01 = cos_a * cos_b;
  K02 = cos_a * cos_c;
  K10 = cos_b * cos_a;
  K11 = cos_b * cos_b;
  K12 = cos_b * cos_c;
  K20 = cos_c * cos_a;
  K21 = cos_c * cos_b;
  K22 = cos_c * cos_c;
}

double PoroElastic::FindDistanceToSegment(double x1, double y1, double z1,
                                          double x2, double y2, double z2,
                                          double pointX, double pointY,
                                          double pointZ)
{
  double diffX = x2 - x1;
  double diffY = y2 - y1;
  double diffZ = z2 - z1;
  if ((diffX == 0) && (diffY == 0) && (diffZ == 0))
  {
    diffX = pointX - x1;
    diffY = pointY - y1;
    diffZ = pointZ - z1;
    return sqrt(diffX * diffX + diffY * diffY + diffZ * diffZ);
  }

  double t =
      ((pointX - x1) * diffX + (pointY - y1) * diffY + (pointZ - z1) * diffZ) /
      (diffX * diffX + diffY * diffY + diffZ * diffZ);

  if (t < 0)
  {
    // point is nearest to the first point i.e x1 and y1
    diffX = pointX - x1;
    diffY = pointY - y1;
    diffZ = pointZ - z1;
  }
  else if (t > 1)
  {
    // point is nearest to the end point i.e x2 and y2
    diffX = pointX - x2;
    diffY = pointY - y2;
    diffZ = pointZ - z2;
  }
  else
  {
    // if perpendicular line intersect the line segment.
    diffX = pointX - (x1 + t * diffX);
    diffY = pointY - (y1 + t * diffY);
    diffZ = pointZ - (z1 + t * diffZ);
  }

  // returning shortest distance
  return sqrt(diffX * diffX + diffY * diffY + diffZ * diffZ);
}

void PoroElastic::assemble_darcy(
    EquationSystems &es, const std::string &libmesh_dbg_var(system_name))
{
  // Get a constant reference to the mesh object.
  const MeshBase &mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  LinearImplicitSystem &darcy_system =
      es.get_system<LinearImplicitSystem>("darcy_velocity");
  System &ppore_system = es.get_system<System>("pPoreSystem");

  unsigned int u_var, v_var, w_var, ppore_var;

  // Numeric ids corresponding to each variable in the system
  u_var = darcy_system.variable_number("w_x");
  v_var = darcy_system.variable_number("w_y");
#if (MESH_DIMENSION == 3)
  w_var = darcy_system.variable_number("w_z");
#endif
  ppore_var = ppore_system.variable_number("pPoreVar");

  System &K_system = es.get_system<System>("KSystem");
  const DofMap &K_dof_map = K_system.get_dof_map();

  NumericVector<double> &K_solution = *(K_system.solution);
  K_solution.close();
  K_solution.localize(*K_system.current_local_solution);

  // Get the Finite Element type for "u".  Note this will be
  // the same as the type for "v".
  FEType fe_vel_type = darcy_system.variable_type(u_var);
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));
  QGauss qrule(dim, CONSTANT);
  fe_vel->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_vel_type));
  QGauss qface(dim - 1, CONSTANT); // Not sure what the
                                   // best accuracy is here

  fe_face->attach_quadrature_rule(&qface);

  const std::vector<double> &JxW = fe_vel->get_JxW();
  const std::vector<std::vector<double>> &phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe_vel->get_dphi();

  const DofMap &dof_map = darcy_system.get_dof_map();
  const DofMap &dof_map_ppore = ppore_system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<double> Ke;
  DenseVector<double> Fe;

#if (MESH_DIMENSION == 2)
  DenseSubMatrix<double> Kuu(Ke), Kuv(Ke), Kvu(Ke), Kvv(Ke);

  DenseSubVector<double> Fu(Fe), Fv(Fe);
#elif (MESH_DIMENSION == 3)
  DenseSubMatrix<double> Kuu(Ke), Kuv(Ke), Kuw(Ke), Kvu(Ke), Kvv(Ke), Kvw(Ke),
      Kwu(Ke), Kwv(Ke), Kww(Ke);

  DenseSubVector<double> Fu(Fe), Fv(Fe), Fw(Fe);
#endif

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_uvw(dim);

  std::vector<dof_id_type> dof_indices_ppore;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem *elem = *el;

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);

    for (unsigned int var = 0; var < dim; var++)
    {
      dof_map.dof_indices(elem, dof_indices_uvw[var], var);
    }

    dof_map_ppore.dof_indices(elem, dof_indices_ppore, ppore_var);

    unsigned int n_dofs, n_u_dofs, n_v_dofs, n_w_dofs;

    n_dofs = dof_indices.size();
    n_u_dofs = dof_indices_uvw[0].size();
    n_v_dofs = dof_indices_uvw[1].size();
#if (MESH_DIMENSION == 3)
    n_w_dofs = dof_indices_uvw[2].size();
#endif

    fe_vel->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    Kuu.reposition(u_var * n_u_dofs, u_var * n_u_dofs, n_u_dofs, n_u_dofs);
    Kuv.reposition(u_var * n_u_dofs, v_var * n_u_dofs, n_u_dofs, n_v_dofs);
#if (MESH_DIMENSION == 3)
    Kuw.reposition(u_var * n_u_dofs, w_var * n_u_dofs, n_u_dofs, n_w_dofs);
#endif

    Kvu.reposition(v_var * n_v_dofs, u_var * n_v_dofs, n_v_dofs, n_u_dofs);
    Kvv.reposition(v_var * n_v_dofs, v_var * n_v_dofs, n_v_dofs, n_v_dofs);
#if (MESH_DIMENSION == 3)
    Kvw.reposition(v_var * n_v_dofs, w_var * n_v_dofs, n_v_dofs, n_w_dofs);
#endif

#if (MESH_DIMENSION == 3)
    {
      Kwu.reposition(w_var * n_u_dofs, u_var * n_u_dofs, n_w_dofs, n_u_dofs);
      Kwv.reposition(w_var * n_u_dofs, v_var * n_u_dofs, n_w_dofs, n_v_dofs);
      Kww.reposition(w_var * n_u_dofs, w_var * n_u_dofs, n_w_dofs, n_w_dofs);
    }
#endif

    Fu.reposition(u_var * n_u_dofs, n_u_dofs);
    Fv.reposition(v_var * n_u_dofs, n_v_dofs);
#if (MESH_DIMENSION == 3)
    Fw.reposition(w_var * n_u_dofs, n_w_dofs);
#endif

#if (MESH_DIMENSION == 2)
    const int dof_index_00 = elem->dof_number(K_system.number(), 0, 0);
    double K_00 = K_solution(dof_index_00);
    const int dof_index_01 = elem->dof_number(K_system.number(), 1, 0);
    double K_01 = K_solution(dof_index_01);

    const int dof_index_10 = elem->dof_number(K_system.number(), 2, 0);
    double K_10 = K_solution(dof_index_10);
    const int dof_index_11 = elem->dof_number(K_system.number(), 3, 0);
    double K_11 = K_solution(dof_index_11);

#elif (MESH_DIMENSION == 3)
    const int dof_index_00 = elem->dof_number(K_system.number(), 0, 0);
    double K_00 = K_solution(dof_index_00);
    const int dof_index_01 = elem->dof_number(K_system.number(), 1, 0);
    double K_01 = K_solution(dof_index_01);
    const int dof_index_02 = elem->dof_number(K_system.number(), 2, 0);
    double K_02 = K_solution(dof_index_02);

    const int dof_index_10 = elem->dof_number(K_system.number(), 3, 0);
    double K_10 = K_solution(dof_index_10);
    const int dof_index_11 = elem->dof_number(K_system.number(), 4, 0);
    double K_11 = K_solution(dof_index_11);
    const int dof_index_12 = elem->dof_number(K_system.number(), 5, 0);
    double K_12 = K_solution(dof_index_12);

    const int dof_index_20 = elem->dof_number(K_system.number(), 6, 0);
    double K_20 = K_solution(dof_index_20);
    const int dof_index_21 = elem->dof_number(K_system.number(), 7, 0);
    double K_21 = K_solution(dof_index_21);
    const int dof_index_22 = elem->dof_number(K_system.number(), 8, 0);
    double K_22 = K_solution(dof_index_22);
#endif

    DenseMatrix<double> K_perm(MESH_DIMENSION, MESH_DIMENSION);
    /*  K_perm(0, 0) = K_00;
     K_perm(2, 2) = K_22;
     K_perm(0, 1) = K_01;
     K_perm(0, 2) = K_02;

     K_perm(1, 0) = K_10;
     K_perm(1, 1) = K_11;
     K_perm(1, 2) = K_12;

     K_perm(2, 0) = K_20;
     K_perm(2, 1) = K_21; */

    K_perm(0, 0) = K_00;
    K_perm(0, 1) = K_01;
#if (MESH_DIMENSION == 3)
    K_perm(0, 2) = K_02;
#endif

    K_perm(1, 0) = K_10;
    K_perm(1, 1) = K_11;
#if (MESH_DIMENSION == 3)
    K_perm(1, 2) = K_12;

    K_perm(2, 0) = K_20;
    K_perm(2, 1) = K_21;
    K_perm(2, 2) = K_22;
#endif

    // K_perm.add(1.0,II);

    // Now we will build the element matrix and right-hand-side.
    // Constructing the RHS requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {

      GeomPar::compute_geoPar(es, elem, qp, phi, dphi);
      DenseVector<double> delta_ppore(dim);
      // VectorValue<double> delta_p(0.0,0.0,0.0);

      for (unsigned int var_i = 0; var_i < dim; var_i++)
      {
        for (unsigned int j = 0; j < n_u_dofs; j++)
          delta_ppore(var_i) +=
              dphi[j][qp](var_i) *
              (ppore_system.current_solution(dof_indices_ppore[j]));
      }

      DenseMatrix<double> KFInvTra(dim, dim);
      KFInvTra.add(1.0, K_perm);
      // KFInvTra.right_multiply(GeomPar::FInvTra);

      DenseVector<double> darcy_vel_vec(dim);
      KFInvTra.vector_mult(darcy_vel_vec, delta_ppore);
      darcy_vel_vec.scale(-permeability);

      for (unsigned int dof_i = 0; dof_i < n_u_dofs; dof_i++)
      {
        DenseVector<double> gradNA;
        gradNA.resize(MESH_DIMENSION);
        for (unsigned int i = 0; i < dim; i++)
        {
          gradNA(i) = dphi[dof_i][qp](i);
        }

        Fu(dof_i) += darcy_vel_vec(0) * phi[dof_i][qp] * JxW[qp];
        Fv(dof_i) += darcy_vel_vec(1) * phi[dof_i][qp] * JxW[qp];
#if (MESH_DIMENSION == 3)
        Fw(dof_i) += darcy_vel_vec(2) * phi[dof_i][qp] * JxW[qp];
#endif

        // Matrix contributions for the uu and vv couplings.
        for (unsigned int dof_j = 0; dof_j < n_u_dofs; dof_j++)
        {
          Kuu(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          Kvv(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
#if (MESH_DIMENSION == 3)
          {
            Kww(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          }
#endif
          if (brinkman == 1)
          {
            DenseMatrix<double> JinvCgradNBx, JinvCgradNBy, JinvCgradNBz;
            JinvCgradNBx.resize(dim, dim);
            JinvCgradNBy.resize(dim, dim);
            JinvCgradNBz.resize(dim, dim);
            for (unsigned int i = 0; i < dim; i++)
            {
              JinvCgradNBx(0, i) = dphi[dof_j][qp](i);
              JinvCgradNBy(1, i) = dphi[dof_j][qp](i);
#if (MESH_DIMENSION == 3)
              JinvCgradNBz(2, i) = dphi[dof_j][qp](i);
#endif
            }

            JinvCgradNBx.left_multiply(GeomPar::CInv);
            JinvCgradNBy.left_multiply(GeomPar::CInv);
#if (MESH_DIMENSION == 3)
            JinvCgradNBz.left_multiply(GeomPar::CInv);
#endif

            JinvCgradNBx.scale(GeomPar::detF);
            JinvCgradNBy.scale(GeomPar::detF);
#if (MESH_DIMENSION == 3)
            JinvCgradNBz.scale(GeomPar::detF);
#endif

            DenseVector<double> JinvCgradNBxVec(dim);
            DenseVector<double> JinvCgradNByVec(dim);
            DenseVector<double> JinvCgradNBzVec(dim);
            for (unsigned int i = 0; i < dim; i++)
            {
              JinvCgradNBxVec(i) = JinvCgradNBx(0, i);
              JinvCgradNByVec(i) = JinvCgradNBy(1, i);
#if (MESH_DIMENSION == 3)
              JinvCgradNBzVec(i) = JinvCgradNBz(2, i);
#endif
            }

            Kuu(dof_i, dof_j) +=
                viscocity * permeability *
                MatVecOper::contractVec(JinvCgradNBxVec, gradNA) * JxW[qp];
            Kvv(dof_i, dof_j) +=
                viscocity * permeability *
                MatVecOper::contractVec(JinvCgradNByVec, gradNA) * JxW[qp];
#if (MESH_DIMENSION == 3)
            Kww(dof_i, dof_j) +=
                viscocity * permeability *
                MatVecOper::contractVec(JinvCgradNBzVec, gradNA) * JxW[qp];
#endif
          }
        }
      }
    }

    if (brinkman == 1)
    {
      for (unsigned int side = 0; side < elem->n_sides(); side++)
        if (elem->neighbor_ptr(side) == libmesh_nullptr)
        {

          const std::vector<std::vector<Real>> &phi_face = fe_face->get_phi();
          const std::vector<Real> &JxW_face = fe_face->get_JxW();

          const std::vector<Point> &normal_face = fe_face->get_normals();
          const std::vector<Point> &Xref = fe_face->get_xyz();

          fe_face->reinit(elem, side);

          for (unsigned int qp = 0; qp < qface.n_points(); qp++)
          {
            // short int bc_id = mesh.boundary_info->boundary_id(elem, side);

            vector<boundary_id_type> bc_id_vec;
            mesh.boundary_info->boundary_ids(elem, side, bc_id_vec);
            short int bc_id =
                bc_id_vec[0];

            if (bc_id == 1 || bc_id == 2 || bc_id == 3 || bc_id == 4 ||
                bc_id == 5) // This is the x=+ boundary
            {
              for (std::size_t i = 0; i < phi_face.size(); i++)
              {
                for (std::size_t j = 0; j < phi_face.size(); j++)
                {
                  Kuu(i, j) +=
                      1.0e12 * phi_face[i][qp] * phi_face[j][qp] * JxW_face[qp];
                  Kvv(i, j) +=
                      1.0e12 * phi_face[i][qp] * phi_face[j][qp] * JxW_face[qp];
#if (MESH_DIMENSION == 3)
                  Kww(i, j) +=
                      1.0e12 * phi_face[i][qp] * phi_face[j][qp] * JxW_face[qp];
#endif
                }
              }
            }
          }
        }
    }
    darcy_system.matrix->add_matrix(Ke, dof_indices);
    darcy_system.rhs->add_vector(Fe, dof_indices);
  } // end of element loop
}

void PoroElastic::solve_flow_system(EquationSystems &es)
{
  LinearImplicitSystem &system_flow =
      es.get_system<LinearImplicitSystem>("flowSystem");

  system_flow.solve();
}

void PoroElastic::solve_mexp_system(EquationSystems &es)
{
  LinearImplicitSystem &system_mexp =
      es.get_system<LinearImplicitSystem>("mExpSystem");

  system_mexp.solve();
}

void PoroElastic::solve_darcy(EquationSystems &es)
{
  LinearImplicitSystem &system_darcy =
      es.get_system<LinearImplicitSystem>("darcy_velocity");

  system_darcy.solve();


  if (brinkman == 1)
  {
    LinearImplicitSystem &system_delw =
        es.get_system<LinearImplicitSystem>("systemDelW");

    system_delw.solve();
  }
}

void PoroElastic::update_ppore(EquationSystems &es)
{
  System &system_ppore = es.get_system<System>("pPoreSystem");
  System &system_pnorm = es.get_system<System>("NonlinearElasticity");

  System &system_flow = es.get_system<System>("flowSystem");

  NumericVector<double> &ppore_solution = *(system_ppore.solution);
  ppore_solution.close();
  ppore_solution.localize(*system_ppore.current_local_solution);

  NumericVector<double> &pnorm_solution = *(system_pnorm.solution);
  pnorm_solution.close();

  NumericVector<double> &m_solution = *(system_flow.solution);
  m_solution.close();
  m_solution.localize(*system_flow.current_local_solution);

  libMesh::MeshBase &mesh = es.get_mesh();
  libMesh::MeshBase::const_node_iterator node_it = mesh.local_nodes_begin();
  libMesh::MeshBase::const_node_iterator node_end = mesh.local_nodes_end();

  for (; node_it != node_end; ++node_it)
  {
    const libMesh::Node *nd = *node_it;

    const unsigned int dof_num_ppore =
        nd->dof_number(system_ppore.number(), 0, 0);
    const unsigned int pnorm_dof_num =
        nd->dof_number(system_pnorm.number(), MESH_DIMENSION, 0);
    const unsigned int m_dof_num = nd->dof_number(system_flow.number(), 0, 0);

    double m_cur = m_solution(m_dof_num);
    double pnorm_cur = -pnorm_solution(pnorm_dof_num);

    ppore_solution.set(dof_num_ppore, pnorm_cur + kappa_0 * m_cur);
  }
  ppore_solution.close();
  ppore_solution.localize(*system_ppore.current_local_solution);
}

void PoroElastic::assemble_delw(
    EquationSystems &es, const std::string &libmesh_dbg_var(system_name))
{
  const MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  LinearImplicitSystem &system =
      es.get_system<LinearImplicitSystem>("systemDelW");
  const unsigned int u_var = system.variable_number("delW00");
  const unsigned int v_var = system.variable_number("delW11");
#if (MESH_DIMENSION == 3)
  const unsigned int w_var = system.variable_number("delW22");
#endif

  const DofMap &dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(u_var);
  UniquePtr<FEBase> fe(FEBase::build(dim, fe_type));
  QGauss qrule(dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);

  const std::vector<Real> &JxW = fe->get_JxW();
  const std::vector<std::vector<double>> &phi = fe->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe->get_dphi();

  LinearImplicitSystem &darcy_system =
      es.get_system<LinearImplicitSystem>("darcy_velocity");
  const DofMap &dof_map_darcy = darcy_system.get_dof_map();

  // Storage for the stress dof indices on each element
  std::vector<dof_id_type> dof_indices;
  std::vector<std::vector<dof_id_type>> dof_indices_darcy(dim);

  DenseMatrix<double> Ke;
  DenseVector<double> Fe;

#if (MESH_DIMENSION == 2)
  DenseSubMatrix<double> K00(Ke), K01(Ke), K10(Ke), K11(Ke);

  DenseSubVector<double> F00(Fe), F01(Fe), F10(Fe), F11(Fe);
#elif (MESH_DIMENSION == 3)
  DenseSubMatrix<double> K00(Ke), K01(Ke), K02(Ke), K10(Ke), K11(Ke), K12(Ke),
      K20(Ke), K21(Ke), K22(Ke);

  DenseSubVector<double> F00(Fe), F01(Fe), F02(Fe), F10(Fe), F11(Fe), F12(Fe),
      F20(Fe), F21(Fe), F22(Fe);
#endif

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  // int count_Jn = 0;

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    dof_map.dof_indices(elem, dof_indices);

    for (unsigned int var = 0; var < dim; var++)
    {
      dof_map_darcy.dof_indices(elem, dof_indices_darcy[var], var);
    }

    const unsigned int n_dofs = dof_indices.size();
    const unsigned int n_var_dofs = dof_indices_darcy[0].size();

    fe->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

#if (MESH_DIMENSION == 2)
    K00.reposition(0 * n_var_dofs, 0 * n_var_dofs, n_var_dofs, n_var_dofs);
    K01.reposition(1 * n_var_dofs, 1 * n_var_dofs, n_var_dofs, n_var_dofs);

    K10.reposition(2 * n_var_dofs, 2 * n_var_dofs, n_var_dofs, n_var_dofs);
    K11.reposition(3 * n_var_dofs, 3 * n_var_dofs, n_var_dofs, n_var_dofs);

    F00.reposition(0 * n_var_dofs, n_var_dofs);
    F01.reposition(1 * n_var_dofs, n_var_dofs);
    F10.reposition(2 * n_var_dofs, n_var_dofs);
    F11.reposition(3 * n_var_dofs, n_var_dofs);
#endif

#if (MESH_DIMENSION == 3)
    K00.reposition(0 * n_var_dofs, 0 * n_var_dofs, n_var_dofs, n_var_dofs);
    K01.reposition(1 * n_var_dofs, 1 * n_var_dofs, n_var_dofs, n_var_dofs);
    K02.reposition(2 * n_var_dofs, 2 * n_var_dofs, n_var_dofs, n_var_dofs);

    K10.reposition(3 * n_var_dofs, 3 * n_var_dofs, n_var_dofs, n_var_dofs);
    K11.reposition(4 * n_var_dofs, 4 * n_var_dofs, n_var_dofs, n_var_dofs);
    K12.reposition(5 * n_var_dofs, 5 * n_var_dofs, n_var_dofs, n_var_dofs);

    K20.reposition(6 * n_var_dofs, 6 * n_var_dofs, n_var_dofs, n_var_dofs);
    K21.reposition(7 * n_var_dofs, 7 * n_var_dofs, n_var_dofs, n_var_dofs);
    K22.reposition(8 * n_var_dofs, 8 * n_var_dofs, n_var_dofs, n_var_dofs);

    F00.reposition(0 * n_var_dofs, n_var_dofs);
    F01.reposition(1 * n_var_dofs, n_var_dofs);
    F02.reposition(2 * n_var_dofs, n_var_dofs);

    F10.reposition(3 * n_var_dofs, n_var_dofs);
    F11.reposition(4 * n_var_dofs, n_var_dofs);
    F12.reposition(5 * n_var_dofs, n_var_dofs);

    F10.reposition(6 * n_var_dofs, n_var_dofs);
    F11.reposition(7 * n_var_dofs, n_var_dofs);
    F12.reposition(8 * n_var_dofs, n_var_dofs);
#endif

    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {
      DenseMatrix<Number> grad_w(dim, dim);
      // Row is variable u1, u2, or u3, column is x, y, or z
      for (unsigned int var_i = 0; var_i < dim; var_i++)
        for (unsigned int var_j = 0; var_j < dim; var_j++)
          for (unsigned int j = 0; j < n_var_dofs; j++)
          {
            grad_w(var_i, var_j) +=
                dphi[j][qp](var_j) *
                darcy_system.current_solution(dof_indices_darcy[var_i][j]);
          }

      GeomPar::compute_geoPar(es, elem, qp, phi, dphi);

      DenseMatrix<double> JinvCgradW(MESH_DIMENSION, MESH_DIMENSION);
      JinvCgradW.add(GeomPar::detF, grad_w);
      JinvCgradW.left_multiply(GeomPar::CInv);

      for (unsigned int dof_i = 0; dof_i < n_var_dofs; dof_i++)
      {
#if (MESH_DIMENSION == 2)
        F00(dof_i) += JinvCgradW(0, 0) * phi[dof_i][qp] * JxW[qp];
        F01(dof_i) += JinvCgradW(0, 1) * phi[dof_i][qp] * JxW[qp];
        F10(dof_i) += JinvCgradW(1, 0) * phi[dof_i][qp] * JxW[qp];
        F11(dof_i) += JinvCgradW(1, 1) * phi[dof_i][qp] * JxW[qp];

#elif (MESH_DIMENSION == 3)
        F00(dof_i) += JinvCgradW(0, 0) * phi[dof_i][qp] * JxW[qp];
        F01(dof_i) += JinvCgradW(0, 1) * phi[dof_i][qp] * JxW[qp];
        F02(dof_i) += JinvCgradW(0, 2) * phi[dof_i][qp] * JxW[qp];
        F10(dof_i) += JinvCgradW(1, 0) * phi[dof_i][qp] * JxW[qp];
        F11(dof_i) += JinvCgradW(1, 1) * phi[dof_i][qp] * JxW[qp];
        F12(dof_i) += JinvCgradW(1, 2) * phi[dof_i][qp] * JxW[qp];
        F20(dof_i) += JinvCgradW(2, 0) * phi[dof_i][qp] * JxW[qp];
        F21(dof_i) += JinvCgradW(2, 1) * phi[dof_i][qp] * JxW[qp];
        F22(dof_i) += JinvCgradW(2, 2) * phi[dof_i][qp] * JxW[qp];
#endif

        for (unsigned int dof_j = 0; dof_j < n_var_dofs; dof_j++)
        {
#if (MESH_DIMENSION == 2)
          K00(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          K01(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          K10(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          K11(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];

#elif (MESH_DIMENSION == 3)
          K00(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          K01(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          K02(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          K10(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          K11(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          K12(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          K20(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          K21(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
          K22(dof_i, dof_j) += phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
#endif
        }
      }
    }
    system.matrix->add_matrix(Ke, dof_indices);
    system.rhs->add_vector(Fe, dof_indices);
  }
}

void PoroElastic::update_source(EquationSystems &es, EquationSystems &es_fluid)
{
  const MeshBase &mesh_fluid = es_fluid.get_mesh();

  LinearImplicitSystem &system_fluid = es_fluid.get_system<LinearImplicitSystem>("flowSystem");
  NumericVector<double> &flow_data = *(system_fluid.current_local_solution);
  flow_data.close();

  vector<double> flow_vec;
  flow_data.localize(flow_vec);

  const DofMap &dof_map_fluid = system_fluid.get_dof_map();
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_p;

  libMesh::MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  System &system_source = es.get_system<System>("sourceSystem");
  unsigned int u_var = system_source.variable_number("sourceVar");
  unsigned int system_source_num = system_source.number();

  double sigma_g = 5.0;

  double a_const = 1.0 / (sigma_g * sqrt(2.0 * M_PI));
  double b_const = -1.0 / (2 * sigma_g * sigma_g);

  FEType fe_vel_type = system_source.variable_type(u_var);
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));
  QGauss qrule(dim, CONSTANT);
  fe_vel->attach_quadrature_rule(&qrule);

  const std::vector<Point> &Xref = fe_vel->get_xyz();

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    int i = elem->id();

    double x_elem = VesselFlow::mesh_data[i].x;
    double y_elem = VesselFlow::mesh_data[i].y;
#if (MESH_DIMENSION == 3)
    double z_elem = VesselFlow::mesh_data[i].z;
#endif
    double source_cur = 0.0;

    /* for (int j = 0; j < VesselFlow::pArt(0).size(); j++)
    {
      int n = VesselFlow::termNum[j];
      const Elem *elem_fluid = mesh_fluid.elem_ptr(n);

      dof_map_fluid.dof_indices(elem_fluid, dof_indices_u, 0);
      dof_map_fluid.dof_indices(elem_fluid, dof_indices_p, 1);

      double dist_2 = pow(x_elem - VesselFlow::vessels_in[n].x2, 2) +
                      pow(y_elem - VesselFlow::vessels_in[n].y2, 2) +
                      pow(z_elem - VesselFlow::vessels_in[n].z2, 2);
      source_cur += a_const * exp(b_const * (dist_2)) * flow_vec[dof_indices_u[1]] *
                    sqrt(VesselFlow::p_0 / VesselFlow::rho_v) * VesselFlow::L_v * VesselFlow::L_v * (1.0e-3 / mesh_volume);

      if (VesselFlow::venous_flow == 1)
      {
        const Elem *elem_fluid_v = mesh_fluid.elem_ptr(n + VesselFlow::vessels_in.size());

        dof_map_fluid.dof_indices(elem_fluid_v, dof_indices_u, 0);
        dof_map_fluid.dof_indices(elem_fluid_v, dof_indices_p, 1);

        source_cur += a_const * exp(b_const * (dist_2)) * flow_vec[dof_indices_u[1]] *
                      sqrt(VesselFlow::p_0 / VesselFlow::rho_v) * VesselFlow::L_v * VesselFlow::L_v * (1.0e-3 / mesh_volume);
      }

      // if(i==0)
      // cout<<"i="<<i<<" j="<<j<<" "<<a_const<<" "<<b_const<<" "<<a_const * exp(b_const * (dist_2))<<" "<<dist_2<<endl;
    } */

     source_cur = source_vess[i];
    // source_cur = near_vess[i];

    // cout<<i<<" "<<VesselFlow::mesh_data[i].elem_id<<endl;

    const int dof_index_source = elem->dof_number(system_source_num, 0, 0);
    system_source.solution->set(dof_index_source, source_cur);
  }

  system_source.solution->close();
  system_source.solution->localize(*system_source.current_local_solution);
}

void PoroElastic::update_flowlarge(EquationSystems &es, EquationSystems &es_fluid)
{
  const MeshBase &mesh_fluid = es_fluid.get_mesh();

  LinearImplicitSystem &system_fluid = es_fluid.get_system<LinearImplicitSystem>("flowSystem");
  NumericVector<double> &flow_data = *(system_fluid.current_local_solution);
  flow_data.close();

  vector<double> flow_vec;
  flow_data.localize(flow_vec);

  const DofMap &dof_map_fluid = system_fluid.get_dof_map();
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_p;

  libMesh::MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  System &system_source = es.get_system<System>("flowLargeSystem");
  unsigned int u_var = system_source.variable_number("flowLargeVar");
  unsigned int system_source_num = system_source.number();

  double sigma_g = 2.0;

  double a_const = 1.0 / (sigma_g * sqrt(2.0 * M_PI));
  double b_const = -1.0 / (2 * sigma_g * sigma_g);

  FEType fe_vel_type = system_source.variable_type(u_var);
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));
  QGauss qrule(dim, CONSTANT);
  fe_vel->attach_quadrature_rule(&qrule);

  const std::vector<Point> &Xref = fe_vel->get_xyz();

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    int i = elem->id();

    double x_elem = VesselFlow::mesh_data[i].x;
    double y_elem = VesselFlow::mesh_data[i].y;
#if (MESH_DIMENSION == 3)
    double z_elem = VesselFlow::mesh_data[i].z;
#endif
    double source_cur = 0.0;

    for (int j = 0; j < VesselFlow::vessels_in.size(); j++)
    {
      int n = j;
      const Elem *elem_fluid = mesh_fluid.elem_ptr(n);

      dof_map_fluid.dof_indices(elem_fluid, dof_indices_u, 0);
      dof_map_fluid.dof_indices(elem_fluid, dof_indices_p, 1);

      double dist_2 = pow(x_elem - VesselFlow::vessels_in[n].x2, 2) +
                      pow(y_elem - VesselFlow::vessels_in[n].y2, 2);
#if (MESH_DIMENSION == 3)
      dist_2 += pow(z_elem - VesselFlow::vessels_in[n].z2, 2);
#endif
      source_cur += a_const * exp(b_const * (dist_2)) * flow_vec[dof_indices_u[1]] *
                    sqrt(VesselFlow::p_0 / VesselFlow::rho_v) * VesselFlow::L_v * VesselFlow::L_v * (1.0e-3 / mesh_volume);

      // if(i==0)
      // cout<<"i="<<i<<" j="<<j<<" "<<a_const<<" "<<b_const<<" "<<a_const * exp(b_const * (dist_2))<<" "<<dist_2<<endl;
    }

    const int dof_index_source = elem->dof_number(system_source_num, 0, 0);
    system_source.solution->set(dof_index_source, source_cur);
  }

  system_source.solution->close();
  system_source.solution->localize(*system_source.current_local_solution);
}

void PoroElastic::read_porous_data(EquationSystems &es)
{
  const MeshBase &mesh = es.get_mesh();
  int num_elem;
  num_elem = mesh.n_elem();
  porous_data.resize(num_elem, 10);

  cout << "num_elem=" << num_elem << endl;

  ifstream file_porous;
  file_porous.open("porous_data.dat");
  for (int i = 0; i < num_elem; i++)
  {
    file_porous >> porous_data(i, 0) >> porous_data(i, 1) >>
        porous_data(i, 2) >> porous_data(i, 3) >> porous_data(i, 4) >>
        porous_data(i, 5) >> porous_data(i, 6) >> porous_data(i, 7) >>
        porous_data(i, 8) >> porous_data(i, 9);
  }
  file_porous.close();

  System &system = es.get_system<System>("phiSystem");
  unsigned int u_var = system.variable_number("phiVar");

  System &K_system = es.get_system<System>("KSystem");

  MeshBase::const_element_iterator el = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_elements_end();
  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    int current_elem = elem->id();

    const int dof_index = elem->dof_number(system.number(), 0, 0);
    system.solution->set(dof_index, porous_data(current_elem, 0) / 1000.0);

    const int dof_index_00 = elem->dof_number(K_system.number(), 0, 0);
    K_system.solution->set(dof_index_00, porous_data(current_elem, 1) * 4.3e-5);
    const int dof_index_01 = elem->dof_number(K_system.number(), 1, 0);
    K_system.solution->set(dof_index_01, porous_data(current_elem, 2) * 4.3e-5);
    const int dof_index_02 = elem->dof_number(K_system.number(), 2, 0);
    K_system.solution->set(dof_index_02, porous_data(current_elem, 3) * 4.3e-5);

    const int dof_index_10 = elem->dof_number(K_system.number(), 3, 0);
    K_system.solution->set(dof_index_10, porous_data(current_elem, 4) * 4.3e-5);
    const int dof_index_11 = elem->dof_number(K_system.number(), 4, 0);
    K_system.solution->set(dof_index_11, porous_data(current_elem, 5) * 4.3e-5);
    const int dof_index_12 = elem->dof_number(K_system.number(), 5, 0);
    K_system.solution->set(dof_index_12, porous_data(current_elem, 6) * 4.3e-5);

    const int dof_index_20 = elem->dof_number(K_system.number(), 6, 0);
    K_system.solution->set(dof_index_20, porous_data(current_elem, 7) * 4.3e-5);
    const int dof_index_21 = elem->dof_number(K_system.number(), 7, 0);
    K_system.solution->set(dof_index_21, porous_data(current_elem, 8) * 4.3e-5);
    const int dof_index_22 = elem->dof_number(K_system.number(), 8, 0);
    K_system.solution->set(dof_index_22, porous_data(current_elem, 9) * 4.3e-5);
  }

  system.solution->close();
  system.solution->localize(*system.current_local_solution);

  K_system.solution->close();
  K_system.solution->localize(*K_system.current_local_solution);
}

void PoroElastic::update_nearest_vessel()
{

  near_vess.resize(VesselFlow::mesh_data.size());

  for (int i = 0; i < VesselFlow::mesh_data.size(); i++)
  {
    double x_elem = VesselFlow::mesh_data[i].x;
    double y_elem = VesselFlow::mesh_data[i].y;
#if (MESH_DIMENSION == 3)
    double z_elem = VesselFlow::mesh_data[i].z;
#endif

    double dist_min = 1.0e10;
    int j_min = 0;
    for (int j = 0; j < VesselFlow::termNum.size(); j++)
    {
      int n = VesselFlow::termNum[j];

      double dist_j = pow(x_elem - VesselFlow::vessels[n].x2, 2) + pow(y_elem - VesselFlow::vessels[n].y2, 2);
#if (MESH_DIMENSION == 3)
      dist_j += pow(z_elem - VesselFlow::vessels[n].z2, 2);
#endif

      if (dist_j < dist_min)
      {
        dist_min = dist_j;
        j_min = n;
      }
    }

    near_vess[i] = j_min;
    
  }
}

void PoroElastic::update_source_vessel(EquationSystems &es)
{
  source_vess.resize(VesselFlow::mesh_data.size());

  const MeshBase &mesh = es.get_mesh();

  LinearImplicitSystem &system = es.get_system<LinearImplicitSystem>("flowSystem");
  NumericVector<double> &flow_data = *(system.current_local_solution);
  flow_data.close();

  const unsigned int u_var = system.variable_number("QVar");
  const unsigned int p_var = system.variable_number("pVar");

  vector<double> flow_vec;
  flow_data.localize(flow_vec);

  const DofMap &dof_map = system.get_dof_map();
  std::vector<dof_id_type> dof_indices_u;
  std::vector<dof_id_type> dof_indices_p;

  for (int i = 0; i < VesselFlow::mesh_data.size(); i++)
  {
    int n = near_vess[i];
    const Elem *elem = mesh.elem_ptr(n);

    dof_map.dof_indices(elem, dof_indices_u, u_var);
    dof_map.dof_indices(elem, dof_indices_p, p_var);

    source_vess[i] = flow_vec[dof_indices_u[1]] * sqrt(VesselFlow::p_0 / VesselFlow::rho_v) * VesselFlow::L_v * VesselFlow::L_v;

    if (VesselFlow::venous_flow == 1)
    {
      int n = near_vess[i] + VesselFlow::vessels_in.size();
      const Elem *elem = mesh.elem_ptr(n);

      dof_map.dof_indices(elem, dof_indices_u, 0);
      dof_map.dof_indices(elem, dof_indices_p, 1);

      source_vess[i] += flow_vec[dof_indices_u[1]] * sqrt(VesselFlow::p_0 / VesselFlow::rho_v) * VesselFlow::L_v * VesselFlow::L_v;
    }
  }
}

void PoroElastic::update_aha(EquationSystems &es)
{

  libMesh::MeshBase &mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();

  System &system_aha = es.get_system<System>("ahaSystem");


  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el)
  {
    const Elem *elem = *el;

    int i = elem->id();

    double x_elem = VesselFlow::mesh_data[i].x;
    double y_elem = VesselFlow::mesh_data[i].y;
#if (MESH_DIMENSION == 3)
    double z_elem = VesselFlow::mesh_data[i].z;
#endif

    double aha_cur = 0.0;

    double angle_e = atan2(y_elem, x_elem)*(180.0/3.14159)+180.0;

#if (MESH_DIMENSION == 3)
    if(z_elem > 40.0)
    {
      if(angle_e <=60.0)
        aha_cur = 1.0;
      else if (angle_e <= 120.0)
        aha_cur = 2.0;
      else if (angle_e <= 180.0)
        aha_cur = 3.0;
      else if (angle_e <= 240.0)
        aha_cur = 4.0;
      else if (angle_e <= 300.0)
        aha_cur = 5.0;
      else
        aha_cur = 6.0;
    }

    else if (z_elem > 20.0)
    {
      if (angle_e <= 60.0)
        aha_cur = 7.0;
      else if (angle_e <= 120.0)
        aha_cur = 8.0;
      else if (angle_e <= 180.0)
        aha_cur = 9.0;
      else if (angle_e <= 240.0)
        aha_cur = 10.0;
      else if (angle_e <= 300.0)
        aha_cur = 11.0;
      else
        aha_cur = 12.0;
    }

    else
    {
      if (angle_e <= 45.0)
        aha_cur = 13.0;
      else if (angle_e <= 135.0)
        aha_cur = 14.0;
      else if (angle_e <= 225.0)
        aha_cur = 15.0;
      else if (angle_e <= 315.0)
        aha_cur = 16.0;
      else
        aha_cur = 13.0;
    }
#endif


    const int dof_index_aha = elem->dof_number(system_aha.number(), 0, 0);
    system_aha.solution->set(dof_index_aha, aha_cur);
  }

  system_aha.solution->close();
  system_aha.solution->localize(*system_aha.current_local_solution);

  NumericVector<double> &aha_data = *(system_aha.current_local_solution);
  aha_data.close();

  vector<double> aha_vec;
  aha_data.localize(aha_vec);

  VesselFlow::ahaTerm.resize(VesselFlow::pArt(0).size());
  for (int n = 0; n < VesselFlow::pArt(0).size(); n++)
  {
    const Elem *elem = mesh.elem_ptr(VesselFlow::nearElemTer[n]);

    const int dof_index_aha = elem->dof_number(system_aha.number(), 0, 0);

    VesselFlow::ahaTerm[n] = aha_vec[dof_index_aha];

    //cout << "aha_elem=" << aha_vec[dof_index_aha] << " " << VesselFlow::ahaTerm[n] <<endl;
  }

  VesselFlow::ahaVolume.resize(16);

  MeshBase::const_element_iterator el_full = mesh.active_elements_begin();
  const MeshBase::const_element_iterator end_el_full =
      mesh.active_elements_end();
  for (; el_full != end_el_full; ++el_full)
  {
    const Elem *elem = *el_full;

    double vol_el = elem->volume();

    const int dof_index_aha = elem->dof_number(system_aha.number(), 0, 0);

    int aha_ind = aha_vec[dof_index_aha];

    VesselFlow::ahaVolume[aha_ind - 1] += vol_el;
  }
}

void PoroElastic::assemble_mexp(
    EquationSystems &es, const std::string &libmesh_dbg_var(system_name))
{
  // Get a constant reference to the mesh object.
  const MeshBase &mesh = es.get_mesh();

  // The dimension that we are running
  const unsigned int dim = mesh.mesh_dimension();

  // Get a reference to the Stokes system object.
  LinearImplicitSystem &flow_system =
      es.get_system<LinearImplicitSystem>("mExpSystem");

  unsigned int u_var;

  System &K_system = es.get_system<System>("KSystem");
  const DofMap &K_dof_map = K_system.get_dof_map();

  System &system_source = es.get_system<System>("sourceSystem");
  unsigned int system_source_num = system_source.number();

  LinearImplicitSystem &darcy_system =
      es.get_system<LinearImplicitSystem>("darcy_velocity");
  const DofMap &dof_map_darcy = darcy_system.get_dof_map();

  std::vector<std::vector<dof_id_type>> dof_indices_darcy(dim);

  // Numeric ids corresponding to each variable in the system
  u_var = flow_system.variable_number("mExpVar");
  FEType fe_vel_type = flow_system.variable_type(u_var);
  UniquePtr<FEBase> fe_vel(FEBase::build(dim, fe_vel_type));
  QGauss qrule(dim, CONSTANT);

  // Tell the finite element objects to use our quadrature rule.
  fe_vel->attach_quadrature_rule(&qrule);

  UniquePtr<FEBase> fe_face(FEBase::build(dim, fe_vel_type));
  QGauss qface(dim - 1,
               CONSTANT); // Not sure what the
                          // best accuracy is here

  fe_face->attach_quadrature_rule(&qface);

  const std::vector<double> &JxW = fe_vel->get_JxW();
  const std::vector<std::vector<double>> &phi = fe_vel->get_phi();
  const std::vector<std::vector<RealGradient>> &dphi = fe_vel->get_dphi();

  const DofMap &dof_map = flow_system.get_dof_map();

  // Define data structures to contain the element matrix
  // and right-hand-side vector contribution.  Following
  // basic finite element terminology we will denote these
  // "Ke" and "Fe".
  DenseMatrix<double> Ke;
  DenseVector<double> Fe;

  DenseSubMatrix<double> Kmm(Ke);

  DenseSubVector<double> Fm(Fe);

  // This vector will hold the degree of freedom indices for
  // the element.  These define where in the global system
  // the element degrees of freedom get mapped.
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_m;

  // Now we will loop over all the elements in the mesh that
  // live on the local processor. We will compute the element
  // matrix and right-hand-side contribution.  Since the mesh
  // will be refined we want to only consider the ACTIVE elements,
  // hence we use a variant of the active_elem_iterator.
  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  double dmdx_left = 0.0, dmdx_right = 0.0;

  for (; el != end_el; ++el)
  {
    // Store a pointer to the element we are currently
    // working on.  This allows for nicer syntax later.
    const Elem *elem = *el;

    // Get the degree of freedom indices for the
    // current element.  These define where in the global
    // matrix and right-hand-side this element will
    // contribute to.
    dof_map.dof_indices(elem, dof_indices);
    dof_map.dof_indices(elem, dof_indices_m, 0);

    for (unsigned int var = 0; var < dim; var++)
    {
      dof_map_darcy.dof_indices(elem, dof_indices_darcy[var], var);
    }

    unsigned int n_dofs, n_u_dofs;

    n_dofs = dof_indices.size();
    n_u_dofs = dof_indices_m.size();

    fe_vel->reinit(elem);

    Ke.resize(n_dofs, n_dofs);
    Fe.resize(n_dofs);

    Kmm.reposition(0, 0, n_u_dofs, n_u_dofs);

    Fm.reposition(0, n_u_dofs);

    const int dof_index_source = elem->dof_number(system_source_num, 0, 0);
    double source_cur = system_source.current_solution(dof_index_source);


    // Now we will build the element matrix and right-hand-side.
    // Constructing the RHS requires the solution and its
    // gradient from the previous timestep.  This must be
    // calculated at each quadrature point by summing the
    // solution degree-of-freedom values by the appropriate
    // weight functions.
    for (unsigned int qp = 0; qp < qrule.n_points(); qp++)
    {

      double m_old = 0.0;
      for (unsigned int j = 0; j < n_u_dofs; j++)
      {
        m_old += phi[j][qp] * (flow_system.current_solution(dof_indices_m[j]));
      }


      DenseVector<double> w_cur(dim),JFinvW(dim);
      for (unsigned int var_i = 0; var_i < dim; var_i++)
      {
        w_cur(var_i) = 0;
        for (unsigned int j = 0; j < n_u_dofs; j++)
        {
          w_cur(var_i) += phi[j][qp] * darcy_system.current_solution(dof_indices_darcy[var_i][j]);
        }
      }

      GeomPar::compute_geoPar(es, elem, qp, phi, dphi);


      GeomPar::FInv.vector_mult(JFinvW, w_cur);
      JFinvW.scale(GeomPar::detF);

          for (unsigned int dof_i = 0; dof_i < n_u_dofs; dof_i++)
      {
        DenseVector<double> gradNA;
        gradNA.resize(MESH_DIMENSION);
        for (unsigned int i = 0; i < dim; i++)
        {
          gradNA(i) = dphi[dof_i][qp](i);
        }

        Fm(dof_i) += ((m_old / dt) * phi[dof_i][qp]) * JxW[qp] +
                     MatVecOper::contractVec(JFinvW, gradNA) * JxW[qp];

        Fm(dof_i) += source_cur * phi[dof_i][qp] * JxW[qp];

        // Matrix contributions for the uu and vv couplings.
        for (unsigned int dof_j = 0; dof_j < n_u_dofs; dof_j++)
        {

          Kmm(dof_i, dof_j) +=
              (1.0 / dt) * phi[dof_i][qp] * phi[dof_j][qp] * JxW[qp];
        }
      }
    } // end of the quadrature point qp-loop

    flow_system.matrix->add_matrix(Ke, dof_indices);
    flow_system.rhs->add_vector(Fe, dof_indices);
  } // end of element loop
}