/*           Trajectory calculation by application of           */
/*           Molei Tao's symplectic integration method          */
/*                                                              */
/*                    Minimal usable example                    */
/*                                                              */  
/*  Method reference:                                           */
/*  M. Tao, "Explicit symplectic approximation of nonseparable  */   
/*           Hamiltonians: algorithm and long time performance",*/   
/*           ArXiv, 7 Sep 2016.                                 */
/*                                                              */
/*  Usage:                                                      */
/*  compile it with 'make' and run it with './symplectic_tao'   */
/*                                                              */
/*  See README.md file for further details.                     */

#include "stdio.h"
#include "auxf.hpp"
#include "ode_sys_solv.hpp"

using namespace std;

int main(int argc, char **argv) 
{
    /* System variables */
    const int dim = 4;                           /* System's dimension    */
    const double U = 20.0;                       /* System's energy scale */
    const double a = 0.25;                       /* System's coupling     */
    const double xo  = 0.1;                      /* Initial position x    */
    const double yo  = 1.5;                      /* Initial position y    */
    const double pxo = 4.0;                      /* Initial momentum px   */
    const double pyo = 0.1;                      /* Initial momentum py   */ 
    const double t_run = 100.0;                  /* Integration time      */

    /* Molei Tao method auxiliary variables */
    bool lowest = false;                         /* Lowest order (2) flag */
    const int order = 2;                         /* Integration order     */
    const double dt = 1e-4;                      /* Time step             */
    const double omega = 500.0;                  /* Biding factor         */
    const double gamma = 1.0 / (2.0 - pow(2.0, 1 / (order + 1)));
    double varn_c[] = {xo, yo, pxo, pyo};        /* Auxiliary copy vector */
    
    if (order % 2 != 0)
    {
        printf("Integration order must be EVEN!");
        exit(EXIT_FAILURE);
    }
    if (order == 2) { lowest = true; }

    /* Integration parameters */
    double mod_inf[] = {-PI, -PI, -1e15, -1e15}; /* Position x|y inferior modulation */
    double mod_top[] = { PI,  PI,  1e15,  1e15}; /* Position x|y superior modulation */
    double args[] = {a, U};                      /* Parameters array                 */
    double varn[] = {xo, yo, pxo, pyo};          /* Posterior (n+1) var. array       */
    double t = 0.0;                              /* Dynamic time                     */
    
    /** Out file header **/
    FILE *fout = fopen("traj.dat", "w");              /* Out file: trajectory points */
    fprintf(fout, "; Trajectory for \n");
    fprintf(fout, "; H(x,y) = px^2 + py^2 + U(cos^2(x) + cos^2(y) + 2 * a * cos(x) * cos(y))\n");
    fprintf(fout, "; U = %.0f\n", U);
    fprintf(fout, "; a = %f\n", a);
    fprintf(fout, "; E = %f\n", H(xo, yo, pxo, pyo, U, a));
    fprintf(fout, "; xo = %f\n", xo);
    fprintf(fout, "; yo = %f\n", yo);
    fprintf(fout, "; pxo = %f\n", pxo);
    fprintf(fout, "; pyo = %f\n", pyo);
    fprintf(fout, "; t = %f\n", t_run);
    fprintf(fout, "; method = sympl-tao\n");
    fprintf(fout, "# [t             x          y           px         py            H]\n");
    fprintf(fout, "%f     %f   %f   %f   %f     %.15f\n", 
                    t, xo, yo, pxo, pyo, H(xo,yo,pxo,pyo,U,a));

    /* Run */
    while (t <= t_run)
    {
        ode_step_symplectic(varn, varn_c, t, dt, omega, gamma, dim, order, lowest, args, system_funct);

        modulation(varn, mod_inf, mod_top, dim);
        modulation(varn_c, mod_inf, mod_top, dim);
 
        /* Print trajectory */
        fprintf(fout, "%f     %f   %f   %f   %f     %.15f\n", 
                        t, varn[0], varn[1], varn[2], varn[3], 
                         H(varn[0], varn[1], varn[2], varn[3], U, a));
    }

    fclose(fout);

    return 0;
}