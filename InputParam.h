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

#ifndef included_InputParam
#define included_InputParam

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
#include "libmesh/mesh_communication.h"

// The nonlinear solver and system we will be using
#include "libmesh/linear_implicit_system.h"
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"

#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_tools.h>

#include "../../source_codes/FEMLDE_7/GeomPar.h"
#include "../../source_codes/FEMLDE_7/LargeDeformationElasticity.h"
#include "../../source_codes/FEMLDE_7/MatVecOper.h"

// typedef struct Dcheck{
//	double
// D[MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION][MESH_DIMENSION];
// };
// #define MESH_DIMENSION 3

using namespace libMesh;
using namespace std;



class InputParam {
public:
  InputParam();
  ~InputParam();

  static double mesh_scale;
  static std::string mesh_file_name, fibre_file_name, sheet_file_name, perm_file_name,
      zone_file_name, zone_file_name_2, zone_file_name_3;
  static unsigned int mesh_centre;
  static unsigned int n_solves,n_total;
  static int write_data_skip, write_data_bound;
  static unsigned int output_terminal;
  static double ttime, dt, time_per,omega;
  static int time_itr;

  static int inertia, trans_soln, heirarchy, anis_perm, solve_hyper, second_order_elem;

  static int strain_model;

  static double nonlinear_abs_tol, nonlinear_rel_tol;
  static unsigned int nonlinear_max_its;

  static double aHO, bHO, afHO, bfHO, asHO, bsHO, afsHO, bfsHO;
  static double G;

  static double alpha_stab, rho_s;

  static double kappa_0;
  static vector<int> zone_parent, zone_parent_2;
  static vector<Point> zone_inlet, zone_inlet_2, zone_inlet_3;
  static vector<double> zone_volumes, zone_volumes_2, zone_flow, zone_flow_2,
      zone_zetadiff_2, zone_zeta_2, zone_beta_daught_2;
  static vector<int> zone_near_elem_2;

  static boundary_id_type clamp_x_size, clamp_y_size;
  static DenseVector<boundary_id_type> clamp_x_bcs, clamp_y_bcs;
#if (MESH_DIMENSION == 3)
  static boundary_id_type clamp_z_size;
  static DenseVector<boundary_id_type> clamp_z_bcs;
#endif

  static boundary_id_type torsion_bc, torsion_type;
  static double torsion_t;

  static boundary_id_type force_x_size, force_y_size;
  static DenseVector<boundary_id_type> force_x_bcs, force_y_bcs;
  static DenseVector<double> force_x_values, force_y_values;
#if (MESH_DIMENSION == 3)
  static boundary_id_type force_z_size;
  static DenseVector<boundary_id_type> force_z_bcs;
  static DenseVector<double> force_z_values;
#endif

  static boundary_id_type pressure_size;
  static DenseVector<boundary_id_type> pressure_bcs;
  static DenseVector<double> pressure_values;

  static double t_load, t_end_diastole, pressure_sys;
  static int pressure_current;

  static double beta_s;

  static int active_tension;
  static double Ta_max;

  static double V_bead, a_bead,bead_disp;

  static double alpha_fib;

  static int viscoelastic;
  static double tau_visela, beta_visela;
  static double permeability;

  static DenseVector<double> Vtaubya, VabyD;
  static int vtaubya_size, vabyd_size;

  static int var_v_a;

  static int porous,brinkman;
  static double viscocity;
  static int flow_solver;

  static void read_input();
  static void read_mesh(Mesh &mesh);
  static void read_mesh_perm(ExodusII_IO &exo_io, Mesh &mesh);
  static void read_zone_data();
};

#endif
