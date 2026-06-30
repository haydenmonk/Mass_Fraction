import rebound
import numpy as np
import sys
import time
import fcntl
import os

def create_sim(m_planet, 
               N_particles, r_min, r_max, a_planet=1.0,
               m_star=1.0, i_difference=0.1,
               seed=10):
    
    rng = np.random.default_rng()
    
    sim = rebound.Simulation()
    sim.units=('yr', 'AU', 'Msun')
    sim.add(m=m_star, hash='star')
    sim.add(m=m_planet, a=a_planet, hash='planet')

    sim.N_active = sim.N


    for i in range(N_particles):
        a = rng.uniform(r_min, r_max)
        e = rng.rayleigh(scale=0.01)
        inc = rng.uniform(-i_difference, i_difference)
        Omega = rng.uniform(0, 2*np.pi)
        omega = rng.uniform(0, 2*np.pi)
        M = rng.uniform(0, 2*np.pi)
        sim.add(m=0, a=a, e=e, inc=inc, Omega=Omega, omega=omega, M=M)
        print(M, Omega)
    sim.move_to_com()
    sim.integrator = 'ias15'
    sim.boundary = 'open'
    sim.configure_box(10000.0)
    return sim

def hill_radius(m_planet, a_planet=1.0, m_star=1.0):
    return a_planet * (m_planet/(3*m_star))**(1/3)

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
    output_directory = '/mnt/home/monkhayd/Simulations/Ejection_modeling/Mass_Fraction/Outputs/Ejection_Results'

    output_file=output_directory + f"/{job_id}_ejection_results.txt"


    m_planet_arr=3.00274e-6*np.logspace(2, 3.5, 5)
    mass_index=int(task_id) % len(m_planet_arr)
    m_planet = m_planet_arr[mass_index]

    tmax_index=int(task_id) // (20*len(m_planet_arr))

    # m_planet=1e-3
    # a_planet=1.0
    N_particles=1
    HR=hill_radius(m_planet)
    r_min = 1.0 - 5*HR
    r_max = 1.0 + 5*HR




    tmax_arr=np.linspace(100, 200_000, 20)
    #tmax=tmax_arr[tmax_index]
    
    sim=create_sim(
        m_planet,
        N_particles, r_min, r_max
    )

    archive_filename='/mnt/home/monkhayd/Simulations/Ejection_modeling/Mass_Fraction/Outputs/Sim_Archives/' + f"{file_prefix}_sim.bin"
    
    particle_hash = sim.particles[2].hash.value
    planet_hash = sim.particles[1].hash.value
    start_time=time.monotonic()
    
    previous_frac=0
    for tmax in tmax_arr:
        sim.integrate(tmax)
        
        end_time=time.monotonic()

        elapsed_time=end_time-start_time
        if elapsed_time > (60*60*6):
            sim.save_to_file(archive_filename, walltime=(60*30))



        print(len(sim.particles))
        particle=sim.particles[-1]
        print(particle.x, particle.y, particle.z)
        if sim.particles[-1].hash.value == planet_hash:
            result=0
        elif sim.particles[-1].hash.value == particle_hash:
            result=1
        else:
            result=2

        # with open(output_file, "w") as f:
        #     f.write(f"{m_planet},{result},{tmax}\n")

        with open(output_file, "a") as f:
            fcntl.flock(f, fcntl.LOCK_EX)
            try:
                f.write(f"{m_planet},{result},{tmax}\n")
                f.flush()
            finally:
                fcntl.flock(f, fcntl.LOCK_UN)
    


