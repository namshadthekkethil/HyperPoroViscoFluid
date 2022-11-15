// All rights reserved.
//

#ifndef included_HyperElasticModel
#define included_HyperElasticModel

// C++ include files that we need
#include <algorithm>
#include <iostream>
#include <math.h>
#include <sstream>

// Basic include file needed for the mesh functionality.
#include "libmesh/boundary_info.h"
#include "libmesh/const_function.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h"
#include "libmesh/libmesh.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/parsed_function.h"
#include "libmesh/perf_log.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/transient_system.h"
#include "libmesh/utility.h"
#include "libmesh/zero_function.h"

// For systems of equations the DenseSubMatrix
// and DenseSubVector provide convenient ways for
// assembling the element matrix and vector on a
// component-by-component basis.
#include "libmesh/dense_submatrix.h"
#include "libmesh/dense_subvector.h"

#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"

// The definition of a geometric element
#include "libmesh/elem.h"

#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

#include <map>

#include <fstream>
#include <iostream>

#include "libmesh/getpot.h"
#include "libmesh/string_to_enum.h"

#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/LargeDeformationElasticity.h"

#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/CellFibre.h"
#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/GeomPar.h"
#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/HOModel.h"
#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/Incompress.h"
#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/MatVecOper.h"
#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/NeoHook.h"
#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/Stabilisation.h"
#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/TensorDer.h"

#include "ActiveContra.h"
#include "BoundaryCond.h"
#include "Unsteady.h"
#include "ViscoElastic.h"

using namespace std;
using namespace libMesh;

typedef struct DHyper {
  double D[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
  double FD[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
} DHyper;

// DarcySolver is a static class that provides data and functions
// required to sovle the darcy flow.
class HyperElasticModel {

public:
  HyperElasticModel();
  ~HyperElasticModel();

  static DenseMatrix<double> S, P;
  static DenseVector<double> Resid;
  static DenseMatrix<double> Jacob;

  static DMAT Dmat;
  static DHyper D_hyper;

  static double mu;

  static void initialise_lde(EquationSystems &es,
                             LargeDeformationElasticity &lde);
  static void define_systems(EquationSystems &es);
  static void init_hyperelastic_model(EquationSystems &es, int rank);
  static void update_hyperelastic_model(EquationSystems &es);
  static void initialize_material_axes(EquationSystems &es, int rank);
  static void compute_PK2_hyper(EquationSystems &es, const Elem *elem);
  static void compute_PK2(EquationSystems &es, const Elem *elem);
  static void compute_D();
  static void compute_mu(EquationSystems &es, const Elem *elem);
  static void compute_pre_itr(EquationSystems &es);
  static void
  compute_param_resid_qp(EquationSystems &es, const Elem *elem, unsigned int qp,
                         const std::vector<std::vector<Real>> &phi,
                         const std::vector<std::vector<RealGradient>> &dphi);
  static void
  compute_param_jacob_qp(EquationSystems &es, const Elem *elem, unsigned int qp,
                         const std::vector<std::vector<Real>> &phi,
                         const std::vector<std::vector<RealGradient>> &dphi);
  static void
  compute_param_resid_dofi(const std::vector<std::vector<Real>> &phi,
                           const std::vector<std::vector<RealGradient>> &dphi,
                           unsigned int qp, unsigned int dof_i);
  static void
  compute_param_jacob_dofi(const std::vector<std::vector<Real>> &phi,
                           const std::vector<std::vector<RealGradient>> &dphi,
                           unsigned int qp, unsigned int dof_i);
  static void compute_param_jacob_dofi_dofj(
      const std::vector<std::vector<Real>> &phi,
      const std::vector<std::vector<RealGradient>> &dphi, unsigned int qp,
      unsigned int dof_i, unsigned int dof_j);
  static void update_total_velocity_displacement(EquationSystems &es);
  static void compute_Jtot(EquationSystems &es);
  static void compute_pext(EquationSystems & es);

  //---------------------------------------------------------
};

#endif
