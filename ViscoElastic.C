
#include "ViscoElastic.h"

double ViscoElastic::tau_visela, ViscoElastic::beta_visela;
DenseMatrix<double> ViscoElastic::Q, ViscoElastic::QS;

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;

ViscoElastic::ViscoElastic() {}
ViscoElastic::~ViscoElastic() {}

void ViscoElastic::define_systems(EquationSystems &es) {

  ExplicitSystem &QS_system = es.add_system<ExplicitSystem>("QSSystem");
  QS_system.add_variable("QS_00", CONSTANT, MONOMIAL);
  QS_system.add_variable("QS_01", CONSTANT, MONOMIAL);
#if (MESH_DIMENSION == 3)
  QS_system.add_variable("QS_02", CONSTANT, MONOMIAL);
#endif
  QS_system.add_variable("QS_10", CONSTANT, MONOMIAL);
  QS_system.add_variable("QS_11", CONSTANT, MONOMIAL);
#if (MESH_DIMENSION == 3)
  {
    QS_system.add_variable("QS_12", CONSTANT, MONOMIAL);
    QS_system.add_variable("QS_20", CONSTANT, MONOMIAL);
    QS_system.add_variable("QS_21", CONSTANT, MONOMIAL);
    QS_system.add_variable("QS_22", CONSTANT, MONOMIAL);
  }
#endif
}

void ViscoElastic::update_QS(EquationSystems &es) {
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

  ExplicitSystem &QS_system = es.get_system<ExplicitSystem>("QSSystem");
  const DofMap &QS_dof_map = QS_system.get_dof_map();

  // To store the stress tensor on each element

  MeshBase::const_element_iterator el = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el =
      mesh.active_local_elements_end();

  for (; el != end_el; ++el) {
    const Elem *elem = *el;

    fe->reinit(elem);

    GeomPar::compute_geoPar(es, elem, 0, phi, dphi);
    HyperElasticModel::compute_PK2_hyper(es, elem);

    compute_Q(es, elem);

    QS.resize(MESH_DIMENSION, MESH_DIMENSION);
    QS.add(exp(-InputParam::dt / tau_visela), Q);
    QS.add(-exp(-InputParam::dt / (2.0 * tau_visela)) * beta_visela,
           HyperElasticModel::S);

#if (MESH_DIMENSION == 2)
    int dof_index = elem->dof_number(QS_system.number(), 0, 0);
    QS_system.solution->set(dof_index, QS(0, 0));
    dof_index = elem->dof_number(QS_system.number(), 1, 0);
    QS_system.solution->set(dof_index, QS(0, 1));
    dof_index = elem->dof_number(QS_system.number(), 2, 0);
    QS_system.solution->set(dof_index, QS(1, 0));
    dof_index = elem->dof_number(QS_system.number(), 3, 0);
    QS_system.solution->set(dof_index, QS(1, 1));
#endif

#if (MESH_DIMENSION == 3)
    int dof_index = elem->dof_number(QS_system.number(), 0, 0);
    QS_system.solution->set(dof_index, QS(0, 0));
    dof_index = elem->dof_number(QS_system.number(), 1, 0);
    QS_system.solution->set(dof_index, QS(0, 1));
    dof_index = elem->dof_number(QS_system.number(), 2, 0);
    QS_system.solution->set(dof_index, QS(0, 2));
    dof_index = elem->dof_number(QS_system.number(), 3, 0);
    QS_system.solution->set(dof_index, QS(1, 0));
    dof_index = elem->dof_number(QS_system.number(), 4, 0);
    QS_system.solution->set(dof_index, QS(1, 1));
    dof_index = elem->dof_number(QS_system.number(), 5, 0);
    QS_system.solution->set(dof_index, QS(1, 2));
    dof_index = elem->dof_number(QS_system.number(), 6, 0);
    QS_system.solution->set(dof_index, QS(2, 0));
    dof_index = elem->dof_number(QS_system.number(), 7, 0);
    QS_system.solution->set(dof_index, QS(2, 1));
    dof_index = elem->dof_number(QS_system.number(), 8, 0);
    QS_system.solution->set(dof_index, QS(2, 2));
#endif
  }

  // Should call close and update when we set vector entries directly
  QS_system.solution->close();
  QS_system.update();
}

void ViscoElastic::compute_Q(EquationSystems &es, const Elem *elem) {
  ExplicitSystem &QS_system = es.get_system<ExplicitSystem>("QSSystem");
  const DofMap &QS_dof_map = QS_system.get_dof_map();

  QS.resize(MESH_DIMENSION, MESH_DIMENSION);
#if (MESH_DIMENSION == 2)
  int dof_index_QS = elem->dof_number(QS_system.number(), 0, 0);
  QS(0, 0) = QS_system.current_solution(dof_index_QS);
  dof_index_QS = elem->dof_number(QS_system.number(), 1, 0);
  QS(0, 1) = QS_system.current_solution(dof_index_QS);
  dof_index_QS = elem->dof_number(QS_system.number(), 2, 0);
  QS(1, 0) = QS_system.current_solution(dof_index_QS);
  dof_index_QS = elem->dof_number(QS_system.number(), 3, 0);
  QS(1, 1) = QS_system.current_solution(dof_index_QS);
#endif

#if (MESH_DIMENSION == 3)
  int dof_index_QS = elem->dof_number(QS_system.number(), 0, 0);
  QS(0, 0) = QS_system.current_solution(dof_index_QS);
  dof_index_QS = elem->dof_number(QS_system.number(), 1, 0);
  QS(0, 1) = QS_system.current_solution(dof_index_QS);
  dof_index_QS = elem->dof_number(QS_system.number(), 2, 0);
  QS(0, 2) = QS_system.current_solution(dof_index_QS);
  dof_index_QS = elem->dof_number(QS_system.number(), 3, 0);
  QS(1, 0) = QS_system.current_solution(dof_index_QS);
  dof_index_QS = elem->dof_number(QS_system.number(), 4, 0);
  QS(1, 1) = QS_system.current_solution(dof_index_QS);
  dof_index_QS = elem->dof_number(QS_system.number(), 5, 0);
  QS(1, 2) = QS_system.current_solution(dof_index_QS);
  dof_index_QS = elem->dof_number(QS_system.number(), 6, 0);
  QS(2, 0) = QS_system.current_solution(dof_index_QS);
  dof_index_QS = elem->dof_number(QS_system.number(), 7, 0);
  QS(2, 1) = QS_system.current_solution(dof_index_QS);
  dof_index_QS = elem->dof_number(QS_system.number(), 8, 0);
  QS(2, 2) = QS_system.current_solution(dof_index_QS);
#endif

  Q.resize(MESH_DIMENSION, MESH_DIMENSION);
  Q.add(1.0, QS);
  Q.add(exp(-InputParam::dt / (2.0 * tau_visela)) * beta_visela,
        HyperElasticModel::S);
}
