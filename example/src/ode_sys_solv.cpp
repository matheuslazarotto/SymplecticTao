#include "ode_sys_solv.hpp"

#include <iostream>

void ode_step_symplectic(double var[], double var_copy[], double &t, double dt, double w, 
                         double gamma, unsigned int system_dim, unsigned int order, bool lowest, 
                         double *args, int (*system_function)(double t, const double s[], double f[], 
                                                              void *args, bool set_pos, bool set_mom))
{
    /* Recursively call 'triple jumps' of Phi maps down to the lowest order (2) */
    if (lowest)
    {
        ode_step_symplectic_lowest(var, var_copy, t, dt, w, system_dim, args, system_function);
    }
    else
    {
        order -= 2;
        if (order > 2)
        {
            double dt_s = gamma * dt;
            double dt_m = (1.0 - 2.0 * gamma) * dt;
            double local_gamma = 1.0 / (2.0 - pow(2.0, 1 / (order + 1)));
            ode_step_symplectic(var, var_copy, t, dt_s, w, local_gamma, system_dim, order, lowest, args, system_function);
            ode_step_symplectic(var, var_copy, t, dt_m, w, local_gamma, system_dim, order, lowest, args, system_function);
            ode_step_symplectic(var, var_copy, t, dt_s, w, local_gamma, system_dim, order, lowest, args, system_function);
        }
        else if (order == 2)
        {
            double dt_s = gamma * dt;
            double dt_m = (1.0 - 2.0 * gamma) * dt;
            ode_step_symplectic_lowest(var, var_copy, t, dt_s, w, system_dim, args, system_function);
            ode_step_symplectic_lowest(var, var_copy, t, dt_m, w, system_dim, args, system_function);
            ode_step_symplectic_lowest(var, var_copy, t, dt_s, w, system_dim, args, system_function);
        }
    }
}

void ode_step_symplectic_lowest(double var[], double var_copy[], double &t, 
                                double dt, double w, unsigned int system_dim, double *args, 
                                int (*system_function)(double t, const double s[], double f[], 
                                                       void *args, bool set_pos, bool set_mom))
{
    /** Performs one integration step of first order of 
     *  the Molei Tao symplectic integrator.
     *  Ref: M. Tao; Explicit symplectic approximation of 
     *              nonseparable Hamiltonians: algorithm 
     *              and long time performance. Arkiv (2016)
     *
     * The method is a symplectic explicit map for integration 
     * of hamiltonian systems. For a hamiltonian system Ha, with 
     * variables (q, p) in [var], it considers a copy of it, Hb, 
     * with variables (q', p') in [var_copy] and iterate the copied 
     * system, exchanging variables from Ha to Hb. The integrated 
     * hamiltonian is given by:
     *
     * H = Ha(q,p') + Hb(q',p) + w Hc(q,q',p,p') 
     *
     * where Hc is a coupling perturbation:
     * 
     * Hc = (|q - q'|^2 + |p - p'|^2) / 2
     * 
     * A integration step is made by the composition of the mappings 
     * Pa, Pb and Pc, as follows:
     * 
     * P = Pa(dt/2) * Pb(dt/2) * Pc(dt) * Pb(dt/2) * Pa(dt/2)
     * 
     * where
     * 
     * Pa: [q ]   [            q            ]
     *     [p ] = [ p - dt * dHa/dq(q,p')   ] 
     *     [q']   [ q' + dt * dHa/dp'(q,p') ]
     *     [p']   [            p'           ]
     * 
     * Pb: [q ]   [ q + dt * dHb/dp(q',p)   ]
     *     [p ] = [            p            ] 
     *     [q']   [            q'           ] 
     *     [p']   [ p' - dt * dHb/dq'(q',p) ]
     *
     * Pc: [q ]         [ (q+q') + R(dt)(q-q') ]
     *     [p ] = 1/2 * [ (p+p') + R(dt)(p-p') ]
     *     [q']         [ (q+q') - R(dt)(q-q') ]
     *     [p']         [ (p+p') - R(dt)(p-p') ]
     * where:
     * 
     * R(dt) := [ cos(2 w dt) I     sin(2 w dt) I ]
     *          [-sin(2 w dt) I     cos(2 w dt) I ]
     * 
     * In this function, the variables [q, p'] are stored in [var_a]
     * and [q', p] in [var_b].
     * 
     * For more details, check the reference.
     **/
    double step = 0.5 * dt;
    double f[system_dim];
    double f_copy[system_dim];

    /** Evaluate derivatives for Pa **/
    system_function(0.0, var, f, args, false, true);
    system_function(0.0, var_copy, f_copy, args, true, false);
    /** Apply Pa **/
    phi_a(var, var_copy, f, f_copy, step, system_dim);

    /** Evaluate derivatives for Pb **/
    system_function(0.0, var, f, args, true, false);
    system_function(0.0, var_copy, f_copy, args, false, true);
    /** Apply Pb **/
    phi_b(var, var_copy, f, f_copy, step, system_dim);
    
    /** Apply Pc **/
    phi_c(var, var_copy, dt, w, system_dim);
    
    /** Evaluate derivatives for Pb **/
    system_function(0.0, var, f, args, true, false);
    system_function(0.0, var_copy, f_copy, args, false, true);
    /** Apply Pb **/
    phi_b(var, var_copy, f, f_copy, step, system_dim);

    /** Evaluate derivates for Pa **/
    system_function(0.0, var, f, args, false, true);
    system_function(0.0, var_copy, f_copy, args, true, false);
    /** Apply Pa **/
    phi_a(var, var_copy, f, f_copy, step, system_dim);

    t += dt;
}

void phi_a(double var[], double var_copy[], const double function[], 
           const double function_copy[], double step, unsigned int system_dim)
{
    /** Apply the map Phi A to the variable vectors of 
     *  the Molei Tao symplectic integration method. **/
    unsigned int n = system_dim / 2;
    for (unsigned int k = 0; k < n; k++)
    {
        var[k+n] += step * function[k+n];
        var_copy[k] += step * function_copy[k];
    }
}

void phi_b(double var[], double var_copy[], const double function[], 
           const double function_copy[], double step, unsigned int system_dim)
{
    /** Apply the map Phi B to the variable vectors of 
     * the Molei Tao symplectic integration method. **/
    unsigned int n = system_dim / 2;
    for (unsigned int k = 0; k < n; k++)
    {
        var[k] += step * function[k];
        var_copy[k+n] += step * function_copy[k+n];
    }
}

void phi_c(double var[], double var_copy[], double step, double w, unsigned int system_dim)
{
    /** Apply the map Pc to the variable vectors of the 
     *  Molei Tao symplectic integration method. **/
    unsigned int n = system_dim / 2;
    double cwd = cos(2.0 * w * step);
    double swd = sin(2.0 * w * step);
    double tmp_var[system_dim];
    double tmp_var_copy[system_dim];

    for (unsigned int k = 0; k < n; k++)
    {
        tmp_var[k] = 0.5 * (var[k] + var_copy[k] + cwd * (var[k] - var_copy[k])
                                                 + swd * (var[k+n] - var_copy[k+n]));
        tmp_var[k+n] = 0.5 * (var[k+n] + var_copy[k+n] - swd * (var[k] - var_copy[k])
                                                       + cwd * (var[k+n] - var_copy[k+n]));
        tmp_var_copy[k] = 0.5 * (var[k] + var_copy[k] - cwd * (var[k] - var_copy[k])
                                                      - swd * (var[k+n] - var_copy[k+n]));
        tmp_var_copy[k+n] = 0.5 * (var[k+n] + var_copy[k+n] + swd * (var[k] - var_copy[k])
                                                            - cwd * (var[k+n] - var_copy[k+n]));
    }
    
    for (unsigned int k = 0; k < system_dim; k++)
    {
        var[k] = tmp_var[k];
        var_copy[k] = tmp_var_copy[k];
    }
}