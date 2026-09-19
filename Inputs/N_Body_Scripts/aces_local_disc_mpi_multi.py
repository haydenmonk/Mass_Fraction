import rebound
import numpy as np
import sys
import time
import fcntl
import os
import pandas as pd

from mpi4py import MPI
from ctypes import cdll, c_double, POINTER, byref

import signal

STOP_FILE = "/tmp/stop_my_simulation"

clib = cdll.LoadLibrary("/scratch/group/p.phy260085.000//Mass_Fraction/Inputs/N_Body_Scripts/Heartbeat/heartbeat.so")


clib.get_planet_cartesian.argtypes = [
    c_double,
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
    POINTER(c_double),
]

xp, yp, zp = c_double(), c_double(), c_double()
vxp, vyp, vzp = c_double(), c_double(), c_double()

clib.get_planet_cartesian.restype = None


def get_planet(t, m_planet):
    clib.get_planet_cartesian(
        t,
        byref(xp),
        byref(yp),
        byref(zp),
        byref(vxp),
        byref(vyp),
        byref(vzp),
    )

    planet = rebound.Particle(
        x=xp.value,
        y=yp.value,
        z=zp.value,
        vx=vxp.value,
        vy=vyp.value,
        vz=vzp.value,
        m=m_planet,
    )

    return planet



# Simulation creation
# ============================================================

def create_sim(
    m_planet,
    r_min,
    r_max,
    a_planet=1.0,
    m_star=1.0,
    i_difference=0.1,
    seed=None,
):

    # IMPORTANT:
    # Each simulation has its own RNG.
    rng = np.random.default_rng()

    sim = rebound.Simulation()
    sim.units = ("yr", "AU", "Msun")

    sim.add(m=m_star)

    sim.add(
        m=m_planet,
        a=a_planet,
        primary=sim.particles[0],
    )

    # Exactly ONE massless particle per simulation
    a = rng.uniform(r_min, r_max)
    print(a)
    e = rng.rayleigh(scale=0.1)
    while e >= 1:
        e = rng.rayleigh(scale=0.1)
    inc = rng.uniform(-i_difference, i_difference)
    Omega = rng.uniform(0, 2 * np.pi)
    omega = rng.uniform(0, 2 * np.pi)
    M = rng.uniform(0, 2 * np.pi)

    sim.add(
        m=0,
        a=a,
        e=e,
        inc=inc,
        Omega=Omega,
        omega=omega,
        M=M,
    )

    sim.move_to_com()

    # Remove star and planet.
    # heartbeat.so supplies their gravity analytically.
    sim.remove(0)
    sim.remove(0)

    sim.additional_forces = clib.planet_star_force

    sim.integrator = "ias15"
    sim.boundary = "open"
    sim.configure_box(100_000.0)

    return sim


def hill_radius(m_planet, a_planet=1.0, m_star=1.0):
    return a_planet * (m_planet / (3 * m_star)) ** (1 / 3)

def chaotic_zone(mass_ratio, a_planet):
    return 1.3 * mass_ratio**(2/7) * a_planet
def write_results(output_file, *values):
    with open(output_file, "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)

        try:
            line = ",".join(map(str, values))
            f.write(line + "\n")
            f.flush()

        finally:
            fcntl.flock(f, fcntl.LOCK_UN)
# def calculate_semimajor_axis(t,particle,m_planet, m_star=1.0):
#     planet=get_planet(t, m_planet)

#     # Star state in the COM frame
#     xs = -(m_planet / m_star) * planet.x
#     ys = -(m_planet / m_star) * planet.y
#     zs = -(m_planet / m_star) * planet.z

#     vxs = -(m_planet / m_star) * planet.vx
#     vys = -(m_planet / m_star) * planet.vy
#     vzs = -(m_planet / m_star) * planet.vz

#     # Particle state relative to the star
#     dx = particle.x - xs
#     dy = particle.y - ys
#     dz = particle.z - zs

#     dvx = particle.vx - vxs
#     dvy = particle.vy - vys
#     dvz = particle.vz - vzs
#     r = np.sqrt(dx**2 + dy**2 + dz**2)

#     v2 = dvx**2 + dvy**2 + dvz**2

#     mu = G * m_star

#     specific_energy = 0.5 * v2 - mu / r

#     a = -mu / (2.0 * specific_energy)

#     return a
def calculate_a_e(t, particle, m_planet, m_star=1.0):
    planet = get_planet(t, m_planet)

    # Star position and velocity in the COM frame
    xs = -(m_planet / m_star) * planet.x
    ys = -(m_planet / m_star) * planet.y
    zs = -(m_planet / m_star) * planet.z

    vxs = -(m_planet / m_star) * planet.vx
    vys = -(m_planet / m_star) * planet.vy
    vzs = -(m_planet / m_star) * planet.vz

    # Particle state relative to the star
    r_vec = np.array([
        particle.x - xs,
        particle.y - ys,
        particle.z - zs
    ])

    v_vec = np.array([
        particle.vx - vxs,
        particle.vy - vys,
        particle.vz - vzs
    ])

    r = np.linalg.norm(r_vec)
    v2 = np.dot(v_vec, v_vec)

    mu = G * m_star

    # Specific orbital energy
    specific_energy = 0.5 * v2 - mu / r

    # Semimajor axis
    a = -mu / (2.0 * specific_energy)

    # Specific angular momentum
    h_vec = np.cross(r_vec, v_vec)
    h2 = np.dot(h_vec, h_vec)

    # Eccentricity
    e = np.sqrt(
        1.0 + 2.0 * specific_energy * h2 / mu**2
    )

    return a, e

def save_remaining_particles(states, comm, rank, output_file, m_planet, T_planet):
    """
    Save the current Cartesian state of every particle that is still active.

    Each MPI rank collects its own remaining particles, then all data are
    gathered to rank 0 and written to one CSV file.
    """

    local_rows = []

    for particle_id, state in states.items():

        # Skip particles that have already finished
        if state["done"]:
            continue

        sim = state["sim"]

        # Safety check
        if sim.N == 0:
            continue

        p = sim.particles[0]
        a,e=calculate_a_e(sim.t,p,m_planet)
        if a<0:
            sim.integrate(sim.t-T_planet)
            a_2, e_2=calculate_a_e(sim.t,p,m_planet)
            if a_2<0:
                a='Unbound'
            #didn't change e so that I can flag Unbound, e<1 as error in calcualtion

        local_rows.append({
 #           "particle_id": state["particle_id"],
            "m_planet": m_planet,
            "time": sim.t,
            "a": a,
            "e":e
        })

    # Send every rank's remaining particles to rank 0
    gathered_rows = comm.gather(local_rows, root=0)

    if rank == 0:

        rows = [
            row
            for rank_rows in gathered_rows
            for row in rank_rows
        ]

        df = pd.DataFrame(rows)
        df.to_csv(output_file, index=False)

        print(
            f"Saved {len(df)} remaining particles to {output_file}",
            flush=True,
        )



# Main
# ============================================================

if __name__ == "__main__":


    # MPI
    # --------------------------------------------------------

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()


    # For 8 MPI ranks
    # gives about 100 separate simulations per rank.
    N_total = 200

    
    fraction_record = []

    stop = False

    peak_slope = -np.inf
    t_peak_slope = None

    #Stopping criterion parameters
    slope_window = 5
    minimum_comparison_steps =20
    consecutive_required = 3
    consecutive_below = 0
    minimum_fraction = 0.1
    max_integration_time = 1e8

    if len(sys.argv) > 1:
        job_id = sys.argv[1]
        task_id = sys.argv[2]
        param_file = sys.argv[3]

    else:
        job_id = "test"
        task_id = "0"
        param_file = None

    task_id_int = int(task_id)

    file_prefix = f"{job_id}-{task_id}"

    output_directory = (
        "/scratch/group/p.phy260085.000//Mass_Fraction/Outputs/Ejection_Results"
    )

    output_file = (
        output_directory
        + f"/{job_id}_ejection_results.txt"
    )

    archive_filename = (
        "/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/Sim_Archives/"
        + f"{file_prefix}_sim.bin"
    )

    particle_filename = (
        "/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/Particle_Tracking/"
            + f"{file_prefix}_initial_particles.txt"
    )

    # --------------------------------------------------------
    # Planet mass / tmax
    # --------------------------------------------------------

    # tmax_mass_dict = {
    #     # 0.00005: 5_000_000,
    #     # 0.0001: 5_000_000,
    #     0.00030027: 4_000_000,
    #     0.00071206: 1_000_000,
    #     0.00168856: 400_000,
    #     0.00400422: 180_000,
    #     0.0094955: 100_000,
    #     0.030027: 100_000,
    #     # 0.1: 100_000,
    # }

    tmax_mass_dict= {
        # 0.1: 100_000,
        0.030027: 100_000,
        0.0094955: 100_000,
        0.00400422: 180_000,
        0.00168856: 400_000,
        0.00071206: 1_000_000,
        0.00030027: 4_000_000,
        # 0.0001: 5_000_000,
        # 0.00005: 5_000_000,
    }

    masses=np.logspace(-2,-3.5,6)
    t_maxes=(2.02731513e+01*masses**(-3/2) + 6.97903308e+04)*5

    # tmax_mass_list = list(tmax_mass_dict.items())

    # parameter_index = (
    #     int(task_id) % len(tmax_mass_list)
    # )
    parameter_index = (
            int(task_id) % len(masses)
        )

    # m_planet, tmax = tmax_mass_list[parameter_index]
    m_planet=masses[parameter_index]
    tmax=t_maxes[parameter_index]

    a_planet = 1.0
    m_star = 1.0
    planet_period = np.sqrt(
        a_planet**3 / (m_star + m_planet)
    )
    cadence = 1e3*planet_period
    HR = hill_radius(
        m_planet,
        a_planet=a_planet,
        m_star=m_star,
    )
    CZ=chaotic_zone(m_planet/m_star, a_planet)

    t_arr = np.arange(0, tmax*planet_period, cadence)


    # C force setup ONCE per MPI rank

    temp_sim = rebound.Simulation()
    temp_sim.units = ("yr", "AU", "Msun")

    G = temp_sim.G

    clib.setup_force(
        c_double(a_planet),
        c_double(0.0),       # e
        c_double(0.0),       # i
        c_double(0.0),       # omega
        c_double(0.0),       # Omega
        c_double(0.0),       # M0
        c_double(0.0),       # t0
        c_double(m_planet),
        c_double(m_star),
        c_double(G),
    )

    # simulations -> MPI ranks

    # rank 0: 0, 8, 16...
    # rank 1: 1, 9, 17...

    my_particle_ids = range(rank, N_total, size)
    states = {}

    e_scale=0.1
    p=0.99
    e_p=e_scale*np.sqrt(-2*np.log(1-p))
    e_eff=e_p
    

    for pid in my_particle_ids:

        rng = np.random.default_rng(pid + 10)
        if rng.integers(0, 2):
            r_min = a_planet
            #r_max = a_planet + 2*np.sqrt(3)* HR
            #r_max=a_planet+a_planet*1.7*m_planet**0.31
            #r_max= a_planet+1.8*a_planet*e_eff**(1/5)*m_planet**(1/5)
            #r_max=a_planet+CZ
            r_max=2.5
        else:
            #r_min = a_planet- 2*np.sqrt(3) * HR
            #r_min=a_planet-CZ
            #r_min=a_planet-1.8*a_planet*e_eff**(1/5)*m_planet**(1/5)
            r_max = a_planet
            r_min= 0.5
            # r_min = max(r_min, 0.0)


        sim = create_sim(
            m_planet=m_planet,
            r_min=r_min,
            r_max=r_max,
            a_planet=a_planet,
            m_star=m_star,
            seed=pid + 10,
        )

        a_i,e_i=calculate_a_e(sim.t, sim.particles[0], m_planet)
        write_results(particle_filename, a_i, e_i, pid)

        states[pid] = {
            "sim": sim,

            # done integrating
            "done": False,

            # details
            "ejected": False,
            "captured": False,
        }

    if rank == 0:
        print(f"Running {N_total} simulations with {size} MPI ranks")

    print(f"Rank {rank}: {len(states)} simulations")

    # Integration
    # ========================================================

    start_time = time.monotonic()

    # This preserves the old behavior
    archive_started = False

    for t in t_arr:

        # Each MPI rank advances ALL of its simulations
        # ----------------------------------------------------

        for pid, state in states.items():

            if state["done"]:
                continue

            sim = state["sim"]

            try:
                sim.integrate(t)
            except rebound.NoParticles:
                write_results(output_file,m_planet,t,"Ejected" )

                state["ejected"] = True
                state["done"] = True
                continue

            planet_centric_e = (
                sim.particles[0]
                .orbit(
                    primary=get_planet(
                        sim.t,
                        m_planet,
                    )
                )
                .e
            )

            if planet_centric_e < 1.0:
                sim.integrate(t + planet_period)
                planet_centric_e_2 = (
                    sim.particles[0]
                    .orbit(
                        primary=get_planet(
                            sim.t,
                            m_planet,
                        )
                    )
                    .e
                )

                if planet_centric_e_2 < 1.0:

                    write_results(
                        output_file,
                        m_planet,
                        t,
                        "Captured",
                    )

                    state["captured"] = True
                    state["done"] = True

            elapsed_time = (
                time.monotonic() - start_time
            )

            if (elapsed_time > 60 * 60 * 6 and not archive_started):
                #Remember that these can possibly be overwritten by other ranks
                sim.save_to_file(archive_filename, walltime=60 * 30)
                archive_started = True

        # All local simulations have now reached checkpoint t
        # ====================================================
        global_shutdown = local_shutdown = os.path.exists(STOP_FILE)

        if global_shutdown:
            if rank == 0:
                print(f"Saving remaining particles at t = {t}...")

            cancelled = True
            break

        local_ejected = sum(
            int(state["ejected"])
            for state in states.values()
        )

        # Global number ejected
        # ----------------------------------------------------

        total_ejected = comm.reduce(
            local_ejected,
            op=MPI.SUM,
            root=0,
        )

        # Global stopping criterion
        # ====================================================

        if rank == 0:

            fraction = total_ejected / N_total
            fraction_record.append(fraction)

            stop = False

            # Estimate peak slope
            # -----------------------------------------------

            if len(fraction_record) > slope_window:

                slope = (
                    fraction_record[-1]
                    - fraction_record[-(slope_window + 1)]
                ) / (slope_window * cadence)

                if slope > peak_slope:
                    peak_slope = slope
                    t_peak_slope = t

            # Comparison window
            # -----------------------------------------------

            if fraction >= minimum_fraction and t_peak_slope is not None:

                comparison_cadence = max(2 * t_peak_slope,minimum_comparison_steps * cadence)

                n_steps_back = int(round(comparison_cadence/ cadence))

                if len(fraction_record) > n_steps_back:

                    difference = (
                        fraction_record[-1]
                        - fraction_record[
                            -(n_steps_back + 1)
                        ]
                    )

                    if difference < 0.05:
                        consecutive_below += 1

                    else:
                        consecutive_below = 0

                    if (
                        consecutive_below
                        >= consecutive_required
                    ):
                        stop = True

                    print(
                        f"t={t:.0f}  "
                        f"fraction={fraction:.5f}  "
                        f"peak_slope={peak_slope:.3e}  "
                        f"t_peak={t_peak_slope:.0f}  "
                        f"comparison cadence={comparison_cadence:.0f}  "
                        f"difference={difference:.5f}  "
                        f"below={consecutive_below}/"
                        f"{consecutive_required}",
                        flush=True,
                    )

        else:
            stop = None
        if t >= max_integration_time:
            stop = True

        # ====================================================
        # Send rank 0's decision to everybody
        # ====================================================

        stop = comm.bcast(
            stop,
            root=0,
        )

        if stop:
            if rank == 0:
                print(
                    f"Stopping at t={t}",
                    flush=True,
                )



            break

    final_position_file = os.path.join(
        output_directory,
        f"remaining_particles_{job_id}_{task_id}.csv"
    )

    save_remaining_particles(
        states=states,
        comm=comm,
        rank=rank,
        output_file=final_position_file,
        m_planet=m_planet,
        T_planet=planet_period
    )
