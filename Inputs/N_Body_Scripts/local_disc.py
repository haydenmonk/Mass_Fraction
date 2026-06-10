import rebound
import numpy as np

def create_sim(m_planet, a_planet,
               N_particles, r_min, r_max,
               m_star=1.0, i_difference=0.1,
               seed=10):
    
    rng = np.random.default_rng(seed)
    
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
    sim.move_to_com()
    sim.integrator = 'ias15'
    return sim

def hill_radius(m_planet, a_planet, m_star=1.0):
    return a_planet * (m_planet/(3*m_star))**(1/3)

if __name__ == "__main__":
    m_planet=1e-3
    a_planet=1.0
    N_particles=1
    HR=hill_radius(m_planet, a_planet)
    r_min = a_planet - 5*HR
    r_max = a_planet + 5*HR

    tmax=1.5e6
    
    sim=create_sim(
        m_planet, a_planet,
        N_particles, r_min, r_max
    )
    sim.boundary = 'open'
    sim.configure_box(10_000.0)
    particle_hash = sim.particles[2].hash.value
    planet_hash = sim.particles[1].hash.value

    sim.integrate(tmax)


    if sim.particles[-1].hash.value == planet_hash:
        print(0)
    elif sim.particles[-1].hash.value == particle_hash:
        print(1)
    else:
        print("Last hash is neither planet nor particle.")

    


