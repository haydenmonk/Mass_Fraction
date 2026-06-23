import rebound
import numpy as np
import sys

def create_sim(m_planet, a_planet,
               N_particles, r_min, r_max,
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

def hill_radius(m_planet, a_planet, m_star=1.0):
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
    # task_id_int = int(task_id)


    file_prefix = f"{job_id}-{task_id}-{bash_id}"
    output_directory = '/mnt/home/monkhayd/Simulations/Ejection_modeling/Mass_Fraction/Outputs/Ejection_Results'

    output_file=output_directory + f"/{file_prefix}_ejection_results.txt"


    grid_dim=[1,1]
    m_planet_arr=3.00274e-6*np.logspace(2, 3.5, 3*grid_dim[0])
    a_planet_arr=np.logspace(0.5, 1.6, 3*grid_dim[1])

    A, M = np.meshgrid(a_planet_arr, m_planet_arr, indexing="ij")
    combos = np.column_stack([A.ravel(), M.ravel()])
    print(combos)
    
    param_index=int(task_id) % combos.shape[0]
    a_planet, m_planet = combos[param_index]

    # m_planet=1e-3
    # a_planet=1.0
    N_particles=1
    HR=hill_radius(m_planet, a_planet)
    r_min = a_planet - 5*HR
    r_max = a_planet + 5*HR

    period=a_planet**(3/2)
    tmax=2_000_000*period
    
    sim=create_sim(
        m_planet, a_planet,
        N_particles, r_min, r_max
    )
    
    particle_hash = sim.particles[2].hash.value
    planet_hash = sim.particles[1].hash.value

    sim.integrate(tmax)

    print(len(sim.particles))
    particle=sim.particles[-1]
    print(particle.x, particle.y, particle.z)
    if sim.particles[-1].hash.value == planet_hash:
        result=0
    elif sim.particles[-1].hash.value == particle_hash:
        result=1
    else:
        result=2

    with open(output_file, "w") as f:
        f.write(f"{m_planet},{a_planet},{result}\n")
    


