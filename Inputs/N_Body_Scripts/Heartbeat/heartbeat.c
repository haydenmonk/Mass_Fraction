#include "rebound.h"
#include <math.h>

const double a = 1.0;
const double ecc = 0.0;
const double inc = 0.0;
const double omega = 0.0;
const double Omega = 0.0;
const double M0 = 0.0;

const double mp = 1.e-3; // planet mass
const double ms = 1.0;

const double t0 = 0.0;

void heartbeat(struct reb_simulation* r){
    // other things to do in heartbeat
}

void planet_star_force(struct reb_simulation* r){
    // FORCES DUE TO PLANET AND STAR ON MASSLESS PARTICLES
    // THE PARTICLES ARE ASSUMED TO BE IN THE CENTER OF MASS FRAME 
    // OF THE STAR-PLANET SYSTEM

    double t = r->t;
    double G = r->G;
    
    // mean motion (planet)
    double n = sqrt(G * (mp + ms) / (a * a * a)); 
    // mean anomaly at time t
    double M = M0 + n * (t - t0);
    double f = reb_M_to_f(ecc, M);

    // calculate position of planet at time t, relative to the star
    double rp = a * (1. - ecc*ecc) / (1. + ecc*cos(f));
    double xps = rp*(cos(Omega)*cos(omega+f)-sin(Omega)*sin(omega+f)*cos(inc));
    double yps = rp*(sin(Omega)*cos(omega+f)+cos(Omega)*sin(omega+f)*cos(inc)); 
    double zps = rp*(sin(omega+f)*sin(inc));

    // calculate positions using COM
    double xp = xps * ms / (ms + mp);
    double yp = yps * ms / (ms + mp);
    double zp = zps * ms / (ms + mp);

    double xs = -mp * xp / ms;
    double ys = -mp * yp / ms;
    double zs = -mp * zp / ms;

    // calculate acceleration on particles due to planet and star
    for(int i=0; i<r->N; i++){
        if(r->particles[i].m == 0.0){ // only apply to massless particles
            // distance to planet
            double dxp = r->particles[i].x - xp;
            double dyp = r->particles[i].y - yp;
            double dzp = r->particles[i].z - zp;

            double rp2 = sqrt(dxp*dxp + dyp*dyp + dzp*dzp);
            double rp3 = rp2 * rp2 * rp2;

            // distance to star
            double dxs = r->particles[i].x - xs;
            double dys = r->particles[i].y - ys;
            double dzs = r->particles[i].z - zs;

            double rs2 = sqrt(dxs*dxs + dys*dys + dzs*dzs);
            double rs3 = rs2 * rs2 * rs2;

            // add forces
            r->particles[i].ax -= G * (ms * dxs / rs3 + mp * dxp / rp3);
            r->particles[i].ay -= G * (ms * dys / rs3 + mp * dyp / rp3);
            r->particles[i].az -= G * (ms * dzs / rs3 + mp * dzp / rp3);
        }
    }
}