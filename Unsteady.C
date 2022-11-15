#include "Unsteady.h"

using namespace libMesh;
using namespace std;

Unsteady::Unsteady() {}

Unsteady::~Unsteady() {}

void Unsteady::update_velocity(EquationSystems &es) {
  NonlinearImplicitSystem &system_dis =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");
  System &system_dis_old = es.get_system<System>("displacement_old");
  System &system_dis_n1 = es.get_system<System>("displacement_n1");
  System &system_vel = es.get_system<System>("velocity");
  System &system_vel_old = es.get_system<System>("velocity_old");
  System &system_acc = es.get_system<System>("acceleration");

  if (InputParam::inertia == 1) {
    system_vel.solution->zero();
    system_vel.solution->add(1.0, *system_dis.current_local_solution);
    system_vel.solution->add(-1.0, *system_dis_old.current_local_solution);
    system_vel.solution->scale(1.0 / InputParam::dt);
  } else if (InputParam::inertia == 2) {
    if (InputParam::time_itr == 0) {
      system_vel.solution->zero();
      system_vel.solution->add(1.0, *system_dis.current_local_solution);
      system_vel.solution->add(-1.0, *system_dis_old.current_local_solution);
      system_vel.solution->scale(1.0 / InputParam::dt);
    } else {
      system_vel.solution->zero();
      system_vel.solution->add(3.0, *system_dis.current_local_solution);
      system_vel.solution->add(-4.0, *system_dis_old.current_local_solution);
      system_vel.solution->add(1.0, *system_dis_n1.current_local_solution);
      system_vel.solution->scale(1.0 / (2.0 * InputParam::dt));
    }
  } else if (InputParam::inertia == 3) {
    system_vel.solution->scale(-1.0);
    system_vel.solution->add(2.0 / InputParam::dt,
                             *system_dis.current_local_solution);
    system_vel.solution->add(-(2.0 / InputParam::dt),
                             *system_dis_old.current_local_solution);
  }

  else if (InputParam::inertia == 5) {
    system_vel.solution->zero();
    system_vel.solution->add(1.0, *system_dis.current_local_solution);
    system_vel.solution->add(-1.0, *system_dis_old.current_local_solution);
    system_vel.solution->scale(1.0 / (0.5 * InputParam::dt));
    system_vel.solution->add(-1.0, *system_vel_old.current_local_solution);
  }

  system_vel.solution->close();
  system_vel.solution->localize(*system_vel.current_local_solution);
}

void Unsteady::update_acceleration(EquationSystems &es) {
  NonlinearImplicitSystem &system_dis =
      es.get_system<NonlinearImplicitSystem>("NonlinearElasticity");
  System &system_dis_old = es.get_system<System>("displacement_old");
  System &system_dis_n1 = es.get_system<System>("displacement_n1");
  System &system_vel = es.get_system<System>("velocity");
  System &system_vel_old = es.get_system<System>("velocity_old");
  System &system_vel_n1 = es.get_system<System>("velocity_n1");
  System &system_acc = es.get_system<System>("acceleration");

  if (InputParam::inertia == 1) {
    system_acc.solution->zero();
    system_acc.solution->add(1.0, *system_vel.current_local_solution);
    system_acc.solution->add(-1.0, *system_vel_old.current_local_solution);
    system_acc.solution->scale(1.0 / InputParam::dt);
  }

  else if (InputParam::inertia == 2) {
    if (InputParam::time_itr == 0) {
      system_acc.solution->zero();
      system_acc.solution->add(1.0, *system_vel.current_local_solution);
      system_acc.solution->add(-1.0, *system_vel_old.current_local_solution);
      system_acc.solution->scale(1.0 / InputParam::dt);
    } else {
      system_acc.solution->zero();
      system_acc.solution->add(1.0, *system_vel.current_local_solution);
      system_acc.solution->add(-4.0 / 3.0,
                               *system_vel_old.current_local_solution);
      system_acc.solution->add(1.0 / 3.0,
                               *system_vel_n1.current_local_solution);
      system_acc.solution->scale(1.0 / ((2.0 / 3.0) * InputParam::dt));
    }
  }

  else if (InputParam::inertia == 3) {
    system_acc.solution->zero();
    system_acc.solution->add(1.0, *system_dis.current_local_solution);
    system_acc.solution->add(-1.0, *system_dis_old.current_local_solution);
    system_acc.solution->add(-InputParam::dt,
                             *system_vel.current_local_solution);
    system_acc.solution->scale(2.0 / (InputParam::dt * InputParam::dt));
  }

  else if (InputParam::inertia == 5) {
    system_acc.solution->zero();
    system_acc.solution->add(1.0, *system_vel.current_local_solution);
    system_acc.solution->add(-1.0, *system_vel_old.current_local_solution);
    system_acc.solution->scale(1.0 / InputParam::dt);
  }

  system_acc.solution->close();
  system_acc.solution->localize(*system_acc.current_local_solution);
}
