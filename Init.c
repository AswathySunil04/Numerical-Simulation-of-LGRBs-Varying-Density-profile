#include "pluto.h"

/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
{
  #if EOS == IDEAL
   g_gamma = 4.0/3.0;
  #endif

  double rho0 = 1.2323829e18;
  double Rs = 1e11;
  double rs = sqrt(x1*x1 + x2*x2);

  if (rs <= Rs) {
    v[RHO] = rho0 * pow(rs,-1.5) * pow((Rs-rs)/Rs,3);
    v[VX1] = 0.0;
    v[VX2] = 0.0;
    v[VX3] = 0.0;
    v[PRS] = v[RHO]/1e6;
  } else {
    v[RHO] = 1.e-14;
    v[PRS] = 1.e-20;
    v[VX1] = 0.0;
    v[VX2] = 0.0;
  }

  v[TRC] = 0.0;
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
{
}

/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
{
  double Ljet, gj, rj, hj, pj, rhoj, V1, U1;
  int i, j, k;
  double *x1, *x2, *x3;

  double cur_time, engine_time;

  // ----------- Stellar + BH parameters -----------
  double rho_0 = 1.2323829e18;
  double Rs = 1e11;
  double alpha = 1.5;
  double M_BH = 5*1.989e33;

  // ----------- Spin evolution (clean initialization) -----------
  static double a_spin = 0.0;
  static double last_time = 0.0;
  static int initialized = 0;

  double a0 = 0.8;     // INITIAL SPIN
  double a_min = 0.2;  // saturation spin
  double k_loss = 50.0;

  cur_time = g_time;
  engine_time = 70*3e10;

  gj = 10.0;
  rj = 1.e8;
  hj = 100.0;

  x1 = grid[IDIR].x;
  x2 = grid[JDIR].x;
  x3 = grid[KDIR].x;

  V1 = sqrt(1.0 - 1.0/(gj*gj));
  U1 = 1.0/sqrt(1.0 - V1*V1);

  // ----------- Initialize spin ONLY once -----------
  if (!initialized) {
    a_spin = a0;
    last_time = cur_time;
    initialized = 1;
  }

  // ----------- Time in seconds -----------
  double t_sec;
  if (cur_time <= 0.0) t_sec = 0.0;
  else t_sec = cur_time / 3e10;

  // ----------- Accretion rate -----------
  double dm_dot = 0.0;

  if (t_sec > 1e-6) {

    double rt = pow(2.0*CONST_G*M_BH*t_sec*t_sec, 1.0/3.0);

    dm_dot = pow(2.0*CONST_G*M_BH,1.0/3.0) *
             pow(t_sec, -1.0/3.0) *
             ((2.66666666667*CONST_PI*rho_0)/(3.0*pow(Rs,3))) *
             (-pow(rt,5.0-alpha) +
              3.0*Rs*pow(rt,4.0-alpha) -
              3.0*Rs*Rs*pow(rt,3.0-alpha) +
              Rs*Rs*Rs*pow(rt,2.0-alpha));
  }

  // ----------- Spin-dependent efficiency -----------
  double eta_t = 1.063*pow(a_spin,4) + 0.395*pow(a_spin,2);

  // ----------- Jet luminosity -----------
  Ljet = dm_dot * CONST_c * CONST_c * eta_t;

  // ----------- FAST spin-down with saturation -----------
  if (cur_time > last_time) {

    double dt = (cur_time - last_time) / 3e10;
    if (dt < 0) dt = 0;

    double factor = a_spin - a_min;
    if (factor < 0.0) factor = 0.0;

    a_spin -= k_loss * factor * pow(a_spin,2.0) *
              (Ljet / (M_BH * CONST_c * CONST_c)) * dt;

    // enforce limits
    if (a_spin < a_min) a_spin = a_min;
    if (a_spin > 0.998) a_spin = 0.998;

    last_time = cur_time;
  }

  // ----------- Debug print -----------
  if (g_stepNumber % 500 == 0) {
    int rank = 0;
    #ifdef PARALLEL
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    #endif
    if (rank == 0) {
      printf("[Step %d] Time=%e | a=%f | Ljet=%e\n",
             g_stepNumber, cur_time, a_spin, Ljet);
    }
  }

  // ----------- Jet properties -----------
  rhoj = Ljet / (gj*gj*rj*rj*CONST_PI*CONST_c*CONST_c*hj*CONST_c);
  pj   = (hj - 1.0) * rhoj / 4.0;

  // ----------- Boundary condition -----------
  if (side == X2_BEG) {
    BOX_LOOP(box,k,j,i) {

      if (x1[i] <= rj && cur_time < engine_time) {

        d->Vc[TRC][k][j][i] = 1.0;
        d->Vc[RHO][k][j][i] = (rhoj > 1e-20) ? rhoj : 1e-20;
        d->Vc[PRS][k][j][i] = (pj   > 1e-20) ? pj   : 1e-20;

        d->Vc[VX1][k][j][i] = 0.0;
        d->Vc[VX2][k][j][i] = V1;

        #if USE_FOUR_VELOCITY == YES
          d->Vc[VX2][k][j][i] = U1 * V1;
        #endif

      } else {

        int jm = 2*JBEG - j - 1;

        d->Vc[TRC][k][j][i] = d->Vc[TRC][k][jm][i];
        d->Vc[RHO][k][j][i] = d->Vc[RHO][k][jm][i];
        d->Vc[PRS][k][j][i] = d->Vc[PRS][k][jm][i];
        d->Vc[VX1][k][j][i] = d->Vc[VX1][k][jm][i];
        d->Vc[VX2][k][j][i] = -d->Vc[VX2][k][jm][i];
      }
    }
  }
}
