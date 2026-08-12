#include "rebound.h"
#include <math.h>

double a = 1.0;
double ecc = 0.0;
double inc = 0.0;
double omega = 0.0;
double Omega = 0.0;
double M0 = 0.0;
double M, f;

double G; // gravitational constant
double n; // mean motion
double mp; // planet mass
double ms; // star mass
double t0 = 0.0;

double xp, yp, zp; // planet position
double vxp, vyp, vzp; // planet velocity
double xs, ys, zs; // star position

double rp, xps, yps, zps; // planet position relative to star
double dxp, dyp, dzp, rp2, rp3; // distance to planet
double dxs, dys, dzs, rs2, rs3; // distance to star

struct reb_simulation* r;

void heartbeat(struct reb_simulation* r){
    // other things to do in heartbeat
}

void setup_force(double planet_a, double planet_ecc, double planet_inc,
    double planet_omega, double planet_Omega, double planet_M0, 
    double sim_t0, double planet_mass, double star_mass, double Gconst){

    a = planet_a;
    ecc = planet_ecc;
    inc = planet_inc;
    omega = planet_omega;
    Omega = planet_Omega;
    M0 = planet_M0;
    t0 = sim_t0;
    mp = planet_mass;
    ms = star_mass;
    G = Gconst;
    n = sqrt(G * (mp + ms) / (a * a * a));

    // fake simulation to calculate planet position and velocity later
    r = reb_simulation_create();
    r->G = G;
    reb_simulation_add_fmt(r, "m", ms);
}

void get_planet_cartesian(double t, double* xp_out, double* yp_out, double* zp_out, double* vxp_out, double* vyp_out, double* vzp_out){
    M = M0 + n * (t - t0);
    
    reb_simulation_add_fmt(r, "m a e inc omega Omega M", mp, a, ecc, inc, omega, Omega, M);

    // move to COM
    reb_simulation_move_to_com(r);
    *xp_out = r->particles[1].x;
    *yp_out = r->particles[1].y;
    *zp_out = r->particles[1].z;
    *vxp_out = r->particles[1].vx;
    *vyp_out = r->particles[1].vy;
    *vzp_out = r->particles[1].vz;

    // remove planet
    reb_simulation_remove_particle(r, 1, 1);
    reb_simulation_move_to_hel(r);
}

void planet_star_force(struct reb_simulation* r){
    // FORCES DUE TO PLANET AND STAR ON MASSLESS PARTICLES
    // THE PARTICLES ARE ASSUMED TO BE IN THE CENTER OF MASS FRAME 
    // OF THE STAR-PLANET SYSTEM
    
    // mean anomaly at time t
    M = M0 + n * (r->t - t0);
    f = reb_M_to_f(ecc, M);

    // calculate position of planet at time t, relative to the star
    rp = a * (1. - ecc*ecc) / (1. + ecc*cos(f));
    xps = rp*(cos(Omega)*cos(omega+f)-sin(Omega)*sin(omega+f)*cos(inc));
    yps = rp*(sin(Omega)*cos(omega+f)+cos(Omega)*sin(omega+f)*cos(inc)); 
    zps = rp*(sin(omega+f)*sin(inc));

    // calculate positions using COM
    xp = xps * ms / (ms + mp);
    yp = yps * ms / (ms + mp);
    zp = zps * ms / (ms + mp);

    xs = -mp * xp / ms;
    ys = -mp * yp / ms;
    zs = -mp * zp / ms;

    // calculate acceleration on particles due to planet and star
    for(int i=0; i<r->N; i++){
        if(r->particles[i].m == 0.0){ // only apply to massless particles
            // distance to planet
            dxp = r->particles[i].x - xp;
            dyp = r->particles[i].y - yp;
            dzp = r->particles[i].z - zp;

            rp2 = sqrt(dxp*dxp + dyp*dyp + dzp*dzp);
            rp3 = rp2 * rp2 * rp2;

            // distance to star
            dxs = r->particles[i].x - xs;
            dys = r->particles[i].y - ys;
            dzs = r->particles[i].z - zs;

            rs2 = sqrt(dxs*dxs + dys*dys + dzs*dzs);
            rs3 = rs2 * rs2 * rs2;

            // add forces
            r->particles[i].ax -= G * (ms * dxs / rs3 + mp * dxp / rp3);
            r->particles[i].ay -= G * (ms * dys / rs3 + mp * dyp / rp3);
            r->particles[i].az -= G * (ms * dzs / rs3 + mp * dzp / rp3);
        }
    }
}