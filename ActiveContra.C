
#include "ActiveContra.h"

double ActiveContra::Ta_t, ActiveContra::T_active, ActiveContra::lambda,
    ActiveContra::stretch, ActiveContra::dTdI4f;
double ActiveContra::mu;

DenseMatrix<double> ActiveContra::Sa;
structAct ActiveContra::Dact;

using namespace std;

// Bring in everything from the libMesh namespace
using namespace libMesh;

ActiveContra::ActiveContra() {}
ActiveContra::~ActiveContra() {}

void ActiveContra::compute_Ta_t() {
  if (InputParam::ttime <= InputParam::t_end_diastole)
    Ta_t = 0.0;
  else if (InputParam::ttime < (InputParam::t_end_diastole + 0.15))
    Ta_t = InputParam::Ta_max *
           (1.0 - exp(-(((InputParam::ttime - InputParam::t_end_diastole) *
                         (InputParam::ttime - InputParam::t_end_diastole)) /
                        0.005)));
  else if (InputParam::ttime < InputParam::t_end_diastole + 0.3)
    Ta_t = InputParam::Ta_max *
           (1.0 -
            exp(-((((InputParam::t_end_diastole + 0.3) - InputParam::ttime) *
                   ((InputParam::t_end_diastole + 0.3) - InputParam::ttime)) /
                  0.005)));
  else
    Ta_t = 0.0;
}

void ActiveContra::compute_PK2_active() {
  lambda = sqrt(HOModel::I4f);
  stretch = max(min(lambda, 1.15), 0.8);

  T_active = Ta_t * (1.0 + 4.9 * (stretch - 1.0));

  Sa.resize(MESH_DIMENSION, MESH_DIMENSION);
  Sa.add(0.5 * T_active, HOModel::dI4fdE);
}

void ActiveContra::compute_D_active() {

  dTdI4f = Ta_t * (4.9 / (2.0 * stretch));

  for (unsigned int i = 0; i < MESH_DIMENSION; i++)
    for (unsigned int j = 0; j < MESH_DIMENSION; j++)
      for (unsigned int k = 0; k < MESH_DIMENSION; k++)
        for (unsigned int l = 0; l < MESH_DIMENSION; l++) {
          Dact.Da[i][j][k][l] =
              0.5 * dTdI4f * HOModel::dI4fdE(i, j) * HOModel::dI4fdE(k, l);
        }
}

void ActiveContra::compute_mu() {
  lambda = sqrt(HOModel::I4f);
  stretch = max(min(lambda, 1.15), 0.8);

  T_active = Ta_t * (1.0 + 4.9 * (stretch - 1.0));
  mu = 0.5 * T_active;
}
