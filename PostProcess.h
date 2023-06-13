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

#ifndef included_PostProcess
#define included_PostProcess

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

#include "../../source_codes/FEMLDE_7/LargeDeformationElasticity.h"
#include "../../source_codes/FEMLDE_7/MatVecOper.h"

#include "InputParam.h"

// typedef struct Dcheck{
//	double
// D[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
// };
// #define MESH_DIMENSION 3

using namespace libMesh;
using namespace std;

class PostProcess {
public:
  PostProcess();
  ~PostProcess();

  static double force_p, force_ppore,force_visc;
  static double force_p_total, force_ppore_total, m_net_total, force_visc_total, force_pmono_total, force_ppore_mono_total;

  static DenseMatrix<double> FP, FPPORE,NETM,FPRATE,FPPORERATE,FTOTALRATE;
  static vector<double> fp_time,fppore_time,ftotal_time;
  static string file_surface_force, file_FbyS, file_FbyS_log;

  static void compute_vA(EquationSystems &es, double &u_A, double &v_A,
                         double &p_A, double &I4f_A, double &ppore_A);
  static void compute_mixture_volume(EquationSystems &es,
                                     double &mixture_volume);
  static void compute_skeleton_volume(EquationSystems &es,
                                      EquationSystems &es_cur, double &J_total,
                                      double &m_total);
  static void compute_L2_norm(EquationSystems &es, double &l2_dis_tot,
                              double &l2_vel_tot, double &l2_p_tot);
  static void compute_surface_ppore_m(EquationSystems &es, double &m_0,
                                      double &ppore_0, double &lambda_0,
                                      double &m_1, double &ppore_1,
                                      double &lambda_1);
  static void compute_surface_force(EquationSystems &es,
                                    EquationSystems &es_cur);
  static void init_postprocess(int rank);
  static void update_postprocess(EquationSystems &es, EquationSystems &es_cur,
                                 int rank);
  static void write_final(int rank);
  static void compute_net_m(EquationSystems &es);
  static void update_force_rate(int i, int j);
};

#endif
