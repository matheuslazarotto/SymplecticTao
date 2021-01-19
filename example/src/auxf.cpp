#include "auxf.hpp"

double V(double x, double y, double U, double a) 
{
    return U * (cos(x) * cos(x) + cos(y) * cos(y) + 2 * a * cos(x) * cos(y));
}

double H(double x, double y, double px, double py, double U, double a) 
{
    return px * px + py * py + V(x, y, U, a);
}

double dqdt(double pi) 
{ 
    return 2.0 * pi; 
}

double dpdt(double qi, double qj, double U, double a) 
{ 
    return 2.0 * U * sin(qi) * (cos(qi) + a * cos(qj)); 
}

void modulation(double *var, const double *mod_inf, const double *mod_top, int system_dim)
{
    /** Applies modulation for a variable vector 'var[system_dim]'
     * The bottom and up modulation values are given in vectors 
     * [mod_inf] and [mod_top], respectively. Each coordinate 
     * mod_inf[i]|top[i] gives the cut value for the i-th 
     * coordinate of [var]. **/
    for (int k = 0; k < system_dim; k += 1)
    {
        if (var[k] > mod_top[k])
        {
            var[k] -= fabs(mod_top[k] - mod_inf[k]);
        }
        else if (var[k] < mod_inf[k])
        {
            var[k] += fabs(mod_top[k] - mod_inf[k]);
        }
    }
}

int system_funct(double t, const double s[], double f[], void *args, bool set_pos, bool set_mom) 
{
    (void)(t); // Avoid unused parameter warning
    double *par = (double *)args;
    double a = par[0];
    double U = par[1];
    
    // Positions update
    if (set_pos)
    {
        f[0] = dqdt(s[2]);
        f[1] = dqdt(s[3]);
    }
    // Momenta update
    if (set_mom)
    {
        f[2] = dpdt(s[0], s[1], U, a);
        f[3] = dpdt(s[1], s[0], U, a);
    }
    
    return 0;
}