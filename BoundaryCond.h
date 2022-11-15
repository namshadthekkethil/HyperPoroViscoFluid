// All rights reserved.
//

#ifndef included_BoundaryCond
#define included_BoundaryCond

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

#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/MatVecOper.h"
#include "InputParam.h"

using namespace std;
using namespace libMesh;

// DarcySolver is a static class that provides data and functions
// required to sovle the darcy flow.
class BoundaryCond {

public:
  BoundaryCond();
  ~BoundaryCond();

  static DenseVector<Number> Re_u, Re_v;
#if (MESH_DIMENSION == 3)
  static DenseVector<Number> Re_w;
#endif
  static unsigned int qf_points, num_nodes;
  static size_t phif_size;

  static DenseVector<double> force_x_t, force_y_t;

#if (MESH_DIMENSION == 3)
  static DenseVector<double> force_z_t;
#endif

  static DenseVector<double> pressure_t;

  static void read_boundaries(EquationSystems &es);
  static void init_boundary_cond(EquationSystems &es);
  static void compute_pressure();
  static void compute_torsion();
  static void clamp_boundaries(EquationSystems &es);
  static void apply_pressure(EquationSystems &es, const Elem *elem,
                             unsigned int side,
                             const std::vector<std::vector<Real>> &phi_face,
                             const std::vector<Real> &JxW_face,
                             const std::vector<Point> &normal_face,
                             const std::vector<std::vector<Real>> &phi_face_cur,
                             const std::vector<Real> &JxW_face_cur,
                             const std::vector<Point> &normal_face_cur,
                             DenseSubVector<double> Re_var[]);
  static void apply_force_x(EquationSystems &es, const Elem *elem,
                            unsigned int side,
                            const std::vector<std::vector<Real>> &phi_face,
                            const std::vector<Real> &JxW_face,
                            const std::vector<Point> &normal_face,
                            const std::vector<std::vector<Real>> &phi_face_cur,
                            const std::vector<Real> &JxW_face_cur,
                            const std::vector<Point> &normal_face_cur,
                            DenseSubVector<double> Re_var[]);
  static void apply_force_y(EquationSystems &es, const Elem *elem,
                            unsigned int side,
                            const std::vector<std::vector<Real>> &phi_face,
                            const std::vector<Real> &JxW_face,
                            const std::vector<Point> &normal_face,
                            const std::vector<std::vector<Real>> &phi_face_cur,
                            const std::vector<Real> &JxW_face_cur,
                            const std::vector<Point> &normal_face_cur,
                            DenseSubVector<double> Re_var[]);
#if (MESH_DIMENSION == 3)
  static void apply_force_z(EquationSystems &es, const Elem *elem,
                            unsigned int side,
                            const std::vector<std::vector<Real>> &phi_face,
                            const std::vector<Real> &JxW_face,
                            const std::vector<Point> &normal_face,
                            const std::vector<std::vector<Real>> &phi_face_cur,
                            const std::vector<Real> &JxW_face_cur,
                            const std::vector<Point> &normal_face_cur,
                            DenseSubVector<double> Re_var[]);
#endif

  static void
  apply_pressure_force(EquationSystems &es, const Elem *elem, unsigned int side,
                       const std::vector<std::vector<Real>> &phi_face,
                       const std::vector<Real> &JxW_face,
                       const std::vector<Point> &normal_face,
                       const std::vector<std::vector<Real>> &phi_face_cur,
                       const std::vector<Real> &JxW_face_cur,
                       const std::vector<Point> &normal_face_cur,
                       DenseSubVector<double> Re_var[]);

  static void define_systems(EquationSystems &es);
  static void compute_bcid(EquationSystems &es);
  static void apply_torsion_resid(EquationSystems &es, const Elem *elem,
                                  DenseSubVector<double> Re_var[]);
  static void
  apply_torsion_jacob(EquationSystems &es, const Elem *elem,
                      DenseSubMatrix<Number> Ke_var[][MESH_DIMENSION + 1]);
  //---------------------------------------------------------
};

#endif
