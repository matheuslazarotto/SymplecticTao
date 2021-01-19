#ifndef _AUXF_H_
#define _AUXF_H_

#include <cassert>
#include <cmath>
#include <cstdlib>

#define PI M_PI

/** Dynamical system functions **/
double V(double x, double y, double U, double a);
double H(double x, double y, double px, double py, double U, double a);
double dqdt(double pi);
double dpdt(double qi, double qj, double U, double a);
void modulation(double *var=nullptr, const double *mod_inf=nullptr, const double *mod_top=nullptr, int system_dim=0);

/** System motion equations **/
int system_funct(double t, const double s[], double f[], void *args, bool set_pos = true, bool set_mom = true);

#endif