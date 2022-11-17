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

// C++ include files that we need
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
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
#include "libmesh/nonlinear_implicit_system.h"
#include "libmesh/nonlinear_solver.h"

#include <libmesh/mesh_modification.h>
#include <libmesh/mesh_refinement.h>
#include <libmesh/mesh_tools.h>

#include "libmesh/petsc_matrix.h"

#include "libmesh/mesh_tetgen_interface.h"

#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/LargeDeformationElasticity.h"

#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/GeomPar.h"
#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/HOModel.h"
#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/Incompress.h"
#include "/work/e642/e642/namshadth/source_codes/FEMLDE_7/MeshGen.h"

#include "BoundaryCond.h"
#include "HyperElasticModel.h"
#include "InputParam.h"
#include "PoroElastic.h"
#include "PostProcess.h"
#include "ViscoElastic.h"

#include "Admittance.h"
#include "VesselFlow.h"

#define FLUIDFLOW 1

using namespace libMesh;
using namespace std;
using namespace std::chrono;

void write_input()
{
  GetPot infile1("input.in");
  std::string input_file_name = infile1("input_file_name", "input_LV.in");
  // GetPot infile("input.in");
  GetPot infile(input_file_name);

  infile.print();
}

double mesh_jacobian(EquationSystems &es)
{
  System &tau_system = es.get_system<System>("JSystem");
  return (tau_system.solution->sum()) / (tau_system.solution->size());
}


void run_time_step_fluid(EquationSystems &es, Mesh &mesh, int rank,
                   LibMeshInit &init, int count_solid)
{
  int dt_ratio = (InputParam::dt/VesselFlow::dt);
  int count_per = 0;
  for (unsigned int count = (count_solid-1)*dt_ratio+1; count <= (count_solid)*dt_ratio; count++)
  {
    VesselFlow::time_itr = count;
    VesselFlow::ttime = count * VesselFlow::dt_v;

    count_per = fmod(count, VesselFlow::N_period);
    VesselFlow::time_itr_per = count_per;
    VesselFlow::ttime_dim = count_per * VesselFlow::dt;

    VesselFlow::n1_old(es);
    VesselFlow::old_new(es);

    auto start = high_resolution_clock::now();

    VesselFlow::solve_flow(es);

    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
    cout << "time for solving=" << duration.count() << endl;
    if (VesselFlow::venous_flow == 1)
    {
      auto start = high_resolution_clock::now();

      VesselFlow::update_partvein(es, rank);

      auto stop = high_resolution_clock::now();
      auto duration = duration_cast<microseconds>(stop - start);
      cout << "time taken for part=" << duration.count() << endl;

      if (VesselFlow::st_tree == 1)
      {
        auto start = high_resolution_clock::now();
        VesselFlow::update_qartvein(rank);
        auto stop = high_resolution_clock::now();

        auto duration = duration_cast<microseconds>(stop - start);
        cout << "time taken for qart=" << duration.count() << endl;
      }
    }

    cout << "count=" << count << " count_per=" << count_per << " t=" << VesselFlow::ttime << " t_per=" << count_per * VesselFlow::dt_v << " tdim=" << count * VesselFlow::dt << " tdim_per=" << VesselFlow::ttime_dim << endl;

    if (((count + 1) % 500 == 0))
    {
      VesselFlow::writeFlowDataTime(es, count, rank);
    }

    if (((count + 1) % 10 == 0))
    {

      VesselFlow::writeFlowDataBound(es, count, rank);
    }

    if ((count + 1) % VesselFlow::N_period == 0)
      VesselFlow::write_restart_data(es, VesselFlow::time_itr, rank);
  }

  

}

void run_time_step(EquationSystems &es, EquationSystems &es_cur, EquationSystems &es_fluid, Mesh &mesh,
                   Mesh &mesh_cur, Mesh &mesh_fluid, LargeDeformationElasticity &lde, int rank,
                   LibMeshInit &init)
{
  ExodusII_IO exo_io(mesh);

  GetPot infile1("input.in");
  std::string output_file_name = infile1("output_file_name", "output.e");

  string name = "output_";
  string fileend = ".e";
  string out_frame;
  char numstr[21];
  sprintf(numstr, "%d", 0);
  out_frame = name + numstr + fileend;

  string name_results = "results_";
  fileend = ".dat";
  string out_results = name_results + numstr + fileend;

  LinearImplicitSystem &system_incomp =
      es.get_system<LinearImplicitSystem>("incomp_system");

  exo_io.write_timestep(out_frame, es, 1, 0);

  InputParam::ttime = 0.0;

  int count_write = 1;

  InputParam::time_itr = 0;

  PoroElastic::time_itr = 0;


  #if (FLUIDFLOW == 1)
  VesselFlow::time_itr = 0;
  VesselFlow::writeFlowDataTime(es_fluid, 0, rank);
  VesselFlow::ttime = 0.0;
  VesselFlow::ttime_dim = 0.0;
  VesselFlow::update_qartvein(rank);
  #endif

  int count_per = 0;

  for (unsigned int count = 1; count < InputParam::n_total; count++)
  {
    #if(FLUIDFLOW == 1)
    HyperElasticModel::compute_pext(es);
    VesselFlow::update_pext(es);

    run_time_step_fluid(es_fluid, mesh_fluid, rank, init,count);
    #endif

    count_per = fmod(count, InputParam::n_solves);
    InputParam::ttime = count_per*InputParam::dt;

    BoundaryCond::compute_pressure();
    if (InputParam::torsion_type == 4)
      BoundaryCond::compute_torsion();

    ActiveContra::compute_Ta_t();

    if (InputParam::pressure_size > 0)
    {
      cout << "pressure=" << BoundaryCond::pressure_t(0)
           << " Ta=" << ActiveContra::Ta_t << " Ta_max=" << InputParam::Ta_max
           << endl;
    }

    lde.solve_lde();
    lde.compute_pmono();

    lde.move_mesh();

    if (InputParam::porous == 1)
      PoroElastic::update_poroelastic(es);

    if (InputParam::porous == 1)
      PostProcess::update_postprocess(es, es_cur, rank);

    HyperElasticModel::update_hyperelastic_model(es);

    if ((count ) % 10 == 0 || InputParam::trans_soln == 0)
    {
      lde.compute_stresses();
      lde.compute_J();
      lde.compute_pmono();

      //HyperElasticModel::compute_pext(es);

      HyperElasticModel::update_total_velocity_displacement(es);
      HyperElasticModel::compute_Jtot(es);

      count_write++;
      exo_io.write_timestep(out_frame, es, count_write, InputParam::ttime);
    }
  }

  #if(FLUIDFLOW == 1)
  VesselFlow::write_restart_data(es_fluid, VesselFlow::time_itr, rank);
  #endif
}

void solve_systems(LibMeshInit &init, int rank, int np)
{
  // read_input();

  if (InputParam::output_terminal == 0)
    freopen("output.txt", "w", stdout);

  write_input();

  Mesh mesh(init.comm());
  Mesh mesh_temp(init.comm());
  Mesh mesh_cur(init.comm());
  Mesh mesh_fluid(init.comm());

#if (FLUIDFLOW == 1)
  VesselFlow::initialise_1Dflow(mesh_fluid, rank, np, init);
  mesh_fluid.print_info();
  VesselFlow::update_mesh_data(mesh_temp);
  VesselFlow::update_nearest_elem();
#endif

  InputParam::read_mesh(mesh);
  InputParam::read_mesh(mesh_cur);
  mesh.print_info();

  EquationSystems equation_systems(mesh);
  EquationSystems equation_systems_cur(mesh_cur);
  EquationSystems equation_systems_fluid(mesh_fluid);

#if (FLUIDFLOW == 1)
  VesselFlow::define_systems(equation_systems_fluid);
#endif

  LargeDeformationElasticity lde(
      equation_systems, equation_systems_cur, InputParam::ttime, InputParam::dt,
      HyperElasticModel::Resid, HyperElasticModel::Jacob);

  cout << "EVERYTHING IS FINE EVERYTHING IS FINE EVERYTHING IS FINE" << endl;
  HyperElasticModel::initialise_lde(equation_systems, lde);

  HyperElasticModel::define_systems(equation_systems);

  if (InputParam::porous == 1)
    PoroElastic::define_systems(equation_systems, rank);

  equation_systems.init();
  equation_systems_cur.init();

#if (FLUIDFLOW == 1)
  equation_systems_fluid.init();
#endif

HyperElasticModel::init_hyperelastic_model(equation_systems,rank);
  if (InputParam::porous == 1)
    PoroElastic::initialise_poroelastic(equation_systems);
  lde.pre_solve();

#if (FLUIDFLOW == 1)
  equation_systems_fluid.parameters.set<unsigned int>(
      "nonlinear solver maximum iterations") = 100;

  VesselFlow::initialise_flow_data(equation_systems_fluid);

  LinearImplicitSystem &flow_system =
      equation_systems_fluid.get_system<LinearImplicitSystem>("flowSystem");
  MatSetOption((dynamic_cast<PetscMatrix<Number> *>(flow_system.matrix))->mat(),
               MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE);
#endif

  
  run_time_step(equation_systems, equation_systems_cur, equation_systems_fluid, 
              mesh, mesh_cur, mesh_fluid, lde, rank,init);

  if (InputParam::output_terminal == 0)
    fclose(stdout);
}

int main(int argc, char **argv)
{
  LibMeshInit init(argc, argv);

  // This example requires the PETSc nonlinear solvers
  // libmesh_example_requires(libMesh::default_solver_package() ==
  // PETSC_SOLVERS,
  //                          "--enable-petsc");

  // We use a 3D domain.
  libmesh_example_requires(LIBMESH_DIM > 2,
                           "--disable-1D-only --disable-2D-only");

  int rank, np;

  MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  InputParam::read_input();
  PostProcess::init_postprocess(rank);

  cout << "Initialise postprocess" << endl;

  auto start = high_resolution_clock::now();

  if (InputParam::var_v_a == 1)
  {
    for (int ii = 0; ii < InputParam::vtaubya_size; ii++)
    {
      for (int jj = 0; jj < InputParam::vabyd_size; jj++)
      {
        InputParam::V_bead =
            sqrt(InputParam::Vtaubya(ii) * InputParam::VabyD(jj) *
                 ((InputParam::permeability * InputParam::G) /
                  InputParam::tau_visela));
        InputParam::a_bead = (InputParam::V_bead * InputParam::tau_visela) /
                             InputParam::Vtaubya(ii);
        InputParam::mesh_scale = InputParam::a_bead / 0.5;
        InputParam::dt = ((InputParam::bead_disp * InputParam::a_bead) / InputParam::V_bead) /
                         InputParam::n_solves;

        cout << "V_bead=" << InputParam::V_bead
             << " InputParam::a_bead=" << InputParam::a_bead
             << " dt=" << InputParam::dt
             << " mesh_scale=" << InputParam::mesh_scale << endl;

        PostProcess::fp_time.resize(0);
        PostProcess::fppore_time.resize(0);
        PostProcess::ftotal_time.resize(0);

        solve_systems(init, rank, np);

        PostProcess::update_force_rate(ii, jj);

        PostProcess::FP(ii, jj) =
            PostProcess::force_p_total / InputParam::a_bead;
        PostProcess::FPPORE(ii, jj) =
            PostProcess::force_ppore_total / InputParam::a_bead;
        PostProcess::NETM(ii, jj) = PostProcess::m_net_total;
      }
    }
    PostProcess::write_final(rank);
  }

  else
  {
    solve_systems(init, rank, np);
  }
  auto stop = high_resolution_clock::now();
  auto duration = duration_cast<microseconds>(stop - start);

  cout << "Time taken by function: " << duration.count() << " microseconds"
       << endl;

  return 0;
}
