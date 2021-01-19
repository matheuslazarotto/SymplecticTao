#ifndef _ODE_SYS_H_
#define _ODE_SYS_H_

#include <cmath>

void ode_step_symplectic(double var[], double var_copy[], double &t, double dt, double w, 
                         double gamma, unsigned int system_dim, unsigned int order, bool lowest, 
                         double *args, int (*system_function)(double t, const double s[], double f[], 
                                                              void *args, bool set_pos, bool set_mom));
void ode_step_symplectic_lowest(double var[], double var_copy[], double &t, double dt, double w, unsigned int system_dim, 
                                double *args, int (*system_function)(double t, const double s[], double f[], 
                                                                     void *args, bool set_pos, bool set_mom));
void phi_a(double var[], double var_copy[], const double function[], const double function_copy[], double step, unsigned int system_dim);
void phi_b(double var[], double var_copy[], const double function[], const double function_copy[], double step, unsigned int system_dim);
void phi_c(double var[], double var_copy[], double step, double w, unsigned int system_dim);

#endif