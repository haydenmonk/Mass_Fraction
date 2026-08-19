import rebound
import numpy as np
import sys
import time
import fcntl
import os

from ctypes import cdll, c_double, POINTER, byref

clib = cdll.LoadLibrary("/scratch/group/p.phy260085.000/Mass_Fraction/Inputs/N_Body_Scripts/Heartbeat/heartbeat.so")

# function to get the planet object at a given time t
clib.get_planet_cartesian.argtypes = [
    c_double,
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
]
xp, yp, zp, vxp, vyp, vzp = c_double(), c_double(), c_double(), c_double(), c_double(), c_double()
clib.get_planet_cartesian.restype = None

def get_planet(t, m_planet):
    clib.get_planet_cartesian(
                t,
                byref(xp), byref(yp), byref(zp),
                byref(vxp), byref(vyp), byref(vzp),
            )
    px, py, pz = xp.value, yp.value, zp.value
    pvx, pvy, pvz = vxp.value, vyp.value, vzp.value

    planet = rebound.Particle(x=px, y=py, z=pz,
        vx=pvx, vy=pvy, vz=pvz,
        m=m_planet
    )

    return planet

def create_sim(m_planet, 
               N_particles, r_min, r_max, a_planet=1.0,
               m_star=1.0, i_difference=0.1,
               seed=10):
    
    rng = np.random.default_rng()
    
    sim = rebound.Simulation()
    sim.units=('yr', 'AU', 'Msun')
    sim.add(m=m_star)
    sim.add(m=m_planet, a=a_planet, primary=sim.particles[0])

    for i in range(N_particles):
        a = rng.uniform(r_min, r_max)
        e = rng.rayleigh(scale=0.01)
        inc = rng.uniform(-i_difference, i_difference)
        Omega = rng.uniform(0, 2*np.pi)
        omega = rng.uniform(0, 2*np.pi)
        M = rng.uniform(0, 2*np.pi)
        sim.add(m=0, a=a, e=e, inc=inc, Omega=Omega, omega=omega, M=M)
    sim.move_to_com()

    # remove star and planet
    sim.remove(0)
    sim.remove(0)

    # this force will act as the star and planet on the particles    
    # set the planet mass
    clib.setup_force(c_double(a_planet), # planet a
                     c_double(0.0),      # ecc
                     c_double(0.0),      # inc
                     c_double(0.0),      # omega
                     c_double(0.0),      # Omega
                     c_double(0.0),      # M0
                     c_double(sim.t),    # time
                     c_double(m_planet), # planet mass
                     c_double(m_star),   # star mass
                     c_double(sim.G)     # gravitational constant
    )
    sim.additional_forces = clib.planet_star_force


    sim.integrator = 'ias15'
    sim.boundary = 'open'
    sim.configure_box(10000.0)
    return sim

def hill_radius(m_planet, a_planet=1.0, m_star=1.0):
    return a_planet * (m_planet/(3*m_star))**(1/3)

def write_results(output_file, *values):
    with open(output_file, "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            line = ",".join(map(str, values))
            f.write(line + "\n")
            f.flush()
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        job_id = sys.argv[1]
        task_id = sys.argv[2]
        bash_id = sys.argv[3]
        param_file = sys.argv[4]
    else:
        job_id = "test"
        task_id = "0"
        bash_id = "0"
    task_id_int = int(task_id)

    file_prefix = f"{job_id}-{task_id}-{bash_id}"
    output_directory =  '/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/Ejection_Results'

    output_file=output_directory + f"/{job_id}_ejection_results.txt"
    archive_filename='/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/Sim_Archives' + f"{file_prefix}_sim.bin"

    tmax_mass_dict={
        0.00005: 5_000_000,
        0.0001: 5_000_000,
        0.00030027: 4_000_000,
        0.00071206: 1_000_000,
        0.00168856: 400_000,
        0.00400422: 180_000,
        0.0094955: 100_000,
        0.030027: 100_000,
        0.1: 100_000,
    }

    tmax_mass_list=list(tmax_mass_dict.items())
        
    # m_planet_arr=3.00274e-6*np.logspace(2, 3.5, 5)
    parameter_index=int(task_id) % len(tmax_mass_list)
    m_planet,tmax = tmax_mass_list[parameter_index]
    print(task_id_int,m_planet)
    # tmax_index=int(task_id) // (20*len(m_planet_arr))

    # m_planet=1e-3
    # a_planet=1.0
    N_particles=1
    HR=hill_radius(m_planet)
    r_min = 1.0 - 5*HR
    r_min=np.max([r_min, 0.0])
    r_max = 1.0 + 5*HR

    # tmax_arr=np.linspace(100, 200_000, 20)

    #tmax=tmax_arr[tmax_index]

    t_arr=np.arange(0,tmax,1e3)
    sim=create_sim(
        m_planet,
        N_particles, r_min, r_max
    )
    
    start_time=time.monotonic()
    archive_started=False

    for t in t_arr:

        try:
            sim.integrate(t)
        except rebound.NoParticles:
            write_results(output_file, m_planet, t, bash_id, 'Ejected')
            break

        if len(sim.particles) == 0:
            write_results(output_file, m_planet, t, bash_id, 'Ejected')
            break

        #sim.integrate(t)

        planet_centric_e=sim.particles[0].orbit(primary=get_planet(sim.t, m_planet)).e

        if planet_centric_e < 1.0:
            sim.integrate(t+1)
            
            planet_centric_e_2=sim.particles[0].orbit(primary=get_planet(sim.t, m_planet)).e
            if planet_centric_e_2 < 1.0:
                write_results(output_file, m_planet, t, bash_id, 'Captured')
                break
        
        end_time=time.monotonic()

        elapsed_time=end_time-start_time
        if elapsed_time > (60*60*6) and not archive_started:
            sim.save_to_file(archive_filename, walltime=(60*30))
            archive_started=True

