// The libMesh Finite Element Library.
// Copyright (C) 2002-2017 Benjamin S. Kirk, John W. Peterson, Roy H. Stogner

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

// <h1> Systems Example 7 - Large deformation elasticity (St. Venant-Kirchoff
// material) </h1> \author Lorenzo Zanon \author David Knezevic \date 2014
//
// In this example, we consider an elastic cantilever beam modeled as a St.
// Venant-Kirchoff material (which is an extension of the linear elastic
// material model to the nonlinear regime). The implementation presented here
// uses NonlinearImplicitSystem.
//
// We formulate the PDE on the reference geometry (\Omega) as opposed to the
// deformed geometry (\Omega^deformed). As a result (e.g. see Ciarlet's 3D
// elasticity book, Theorem 2.6-2) the PDE is given as follows:
//
//     \int_\Omega F_im Sigma_mj v_i,j = \int_\Omega f_i v_i + \int_\Gamma g_i
//     v_i ds
//
// where:
//  * F is the deformation gradient, F = I + du/dx (x here refers to reference
//  coordinates).
//  * Sigma is the second Piola-Kirchoff stress, which for the St. Venant
//  Kirchoff model is
//    given by Sigma_ij = C_ijkl E_kl, where E_kl is the strain,
//    E_kl = 0.5 * (u_k,l + u_l,k + u_m,k u_m,l).
//  * f is a body load.
//  * g is a surface traction on the surface \Gamma.
//
// In this example we only consider a body load (e.g. gravity), hence we set g =
// 0.

#ifndef included_PoroElastic
#define included_PoroElastic

// C++ include files that we need
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

// Various include files needed for the mesh & solver functionality.
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/dof_map.h"
#include "libmesh/elem.h"
#include "libmesh/equation_systems.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/fe.h"
#include "libmesh/getpot.h"
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/mesh_refinement.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/string_to_enum.h"
#include "libmesh/zero_function.h"

// The nonlinear solver and system we will be using
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"

#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

// #include "/home/staff1/nthekkethi/FEMSOURCEGENVELPOR/DarcySolver.h"
#include "../../source_codes/FEMLDE_7/GeomPar.h"
#include "../../source_codes/FEMLDE_7/MatVecOper.h"

#include "InputParam.h"

#include "VesselFlow.h"

// typedef struct Dcheck{
//	double
// D[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
// };
// #define MESH_DIMENSION 3

using namespace libMesh;
using namespace std;


class PoroElastic {
public:
  PoroElastic();
  ~PoroElastic();

  static string approx_order, fe_family;
  static double phi_crit, phi_0;
  static double q1, q2, q3;
  static double ttime, dt;
  static double source_const, source_prssr, sink_const, sink_prssr;
  static double mesh_volume;
  static unsigned int ppore_penalty_on;
  static double permeability;
  static unsigned int inertia;
  static unsigned int incomp_type, pmat_type;
  static double rho_s, kappa_0;
  static DenseMatrix<double> II;
  static int time_itr;
  static double epsilon;
  static double beta_visc;
  static int viscous_on;
  static DenseVector<double> vasc_nodes_x, vasc_nodes_y, vasc_nodes_z;
  static DenseVector<double> edges_1, edges_2;
  static int brinkman;
  static double viscocity;
  static double coord_source[4][21990];
  static DenseMatrix<double> porous_data;
  static DenseVector<vector<double>> source_data;
  static int read_source,read_tree,read_permeability;
  static vector<int> near_vess;
  static vector<double> source_vess;
  static double tau_vispe;


  static void read_input(int rank);
  static void define_systems(EquationSystems &es, int rank);
  static double compute_ppen(double jphi_cur, double pnorm_cur);
  static double compute_pmat(double jphi_cur);
  static double compute_dpmatdjphi(double jphi_cur);
  static double compute_dppendjphi(double jphi_cur, double pnorm_cur);
  static double compute_dppendlambda(double jphi_cur);
  static double compute_d2pmatdjphi2(double jphi_cur);
  static double compute_ppore(double pmono_cur, double pmat_cur,
                              double ppen_cur);
  static void assemble_normal(EquationSystems &es,
                              const std::string &libmesh_dbg_var(system_name));
  static double compute_source(double ppore_cur, double norm_cur);
  static void compute_mesh_volume(EquationSystems &es);
  static void compute_preitr_poro(EquationSystems &es);
  static void assemble_darcy(EquationSystems &es,
                             const std::string &libmesh_dbg_var(system_name));
  static void assemble_dmdt(EquationSystems &es,
                            const std::string &libmesh_dbg_var(system_name));
  static double incomp_cond_poro(double NA, double m_cur, double m_old,
                                 double m_n1, double dmdt_cur);
  static double compute_added_mass_resid(DenseVector<double> &gradNA,
                                         DenseVector<double> &w_cur, double NA,
                                         double detF, double m_cur,
                                         double m_old, double m_n1,
                                         double source_cur);
  static void
  mresid_derivative(DenseMatrix<double> &F, DenseMatrix<double> &FinvTra,
                    DenseMatrix<double> &invC, DenseVector<double> &gradNA,
                    DenseVector<double> &gradNB, DenseVector<double> &gradPpore,
                    DenseVector<double> &gradm, double NA, double NB,
                    double detF, double m_cur, DenseVector<double> &dmresiddu,
                    double &dmresiddm, double &dmresiddlambda,
                    double &dincompdm, double pnorm_cur);
  static void
  compute_w_derivative(DenseMatrix<double> &F, DenseMatrix<double> &FinvTra,
                       DenseMatrix<double> &invC, DenseVector<double> &gradNB,
                       DenseVector<double> gradPpore, DenseVector<double> gradm,
                       double detF, double m_cur, double NB,
                       DenseVector<double> &dWdux, DenseVector<double> &dWduy,
                       DenseVector<double> &dWduz, DenseVector<double> &dWdm,
                       DenseVector<double> &dWdlambda, double pnorm_cur);
  static double compute_source_derivative(double m_cur, double pnorm_cur,
                                          double &dsourcedm,
                                          double &dsourcedlambda);
  static void compute_mmono(EquationSystems &es);
  static void assemble_flow(EquationSystems &es,
                            const std::string &libmesh_dbg_var(system_name));
  static void assemble_porous(EquationSystems &es,
                              const std::string &libmesh_dbg_var(system_name));
  static void assemble_porous_p1p0(EquationSystems &es,
                                   const std::string &libmesh_dbg_var(system_name));
  static void initialise_K(EquationSystems &es);
  static void update_poroelastic(EquationSystems &es);
  static void initialise_poroelastic(EquationSystems &es);
  static void compute_K(int iv, double &K00, double &K01, double &K02,
                        double &K10, double &K11, double &K12, double &K20,
                        double &K21, double &K22);
  static double FindDistanceToSegment(double x1, double y1, double z1,
                                      double x2, double y2, double z2,
                                      double pointX, double pointY,
                                      double pointZ);
  static void solve_flow_system(EquationSystems &es);
  static void solve_darcy(EquationSystems &es);
  static void update_ppore(EquationSystems &es);
  static void assemble_delw(EquationSystems &es,
                            const std::string &libmesh_dbg_var(system_name));
  static void update_source(EquationSystems &es, EquationSystems &es_fluid);
  static void read_porous_data(EquationSystems &es);
  static void update_nearest_vessel();
  static void update_source_vessel(EquationSystems &es);

  static void update_flowlarge(EquationSystems &es, EquationSystems &es_fluid);

  static void update_aha(EquationSystems &es);

  static void assemble_mexp(EquationSystems &es,
                            const std::string &libmesh_dbg_var(system_name));
  static void solve_mexp_system(EquationSystems &es);
  static void define_heir_systems(EquationSystems &es);
  static void read_perm_data(EquationSystems &es, ExodusII_IO &exo_io);
  static void update_source_heir(EquationSystems &es);
  static void assemble_m_heir(
      EquationSystems &es, const std::string &libmesh_dbg_var(system_name));
};

#endif
