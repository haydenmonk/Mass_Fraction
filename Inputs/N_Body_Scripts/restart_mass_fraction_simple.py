"""Restart the original Mass_Fraction code without changing its basic loop.

Run: mpirun -np 24 python restart_mass_fraction_simple.py 2215154 0
The second argument indexes MASS_DIRECTORIES below (the selected list only).

IMPORTANT: Original terminal results MUST include particle_id as a fourth CSV
field (mass,time,status,particle_id). The September 6 version of the original
code did not write this field. This script refuses ambiguous restarts.
"""
import os
import sys
import time
import fcntl
from pathlib import Path
from ctypes import cdll, c_double, POINTER, byref

import numpy as np
import pandas as pd
import rebound
from mpi4py import MPI

# ---------------- USER SETTINGS ----------------
ROOT = Path('/scratch/group/p.phy260085.000/Mass_Fraction')
MASS_DIRECTORIES = [
    # MANUALLY ENTER ORIGINAL ARCHIVE DIRECTORY NAMES HERE:
    '2215154_mass_1.99526231e-05',
    '2215154_mass_3.98107171e-05',
    '2215154_mass_7.94328235e-05',
    '2215154_mass_1.00000000e-05'
]
N_total = 800                 # Must match the original simulation for these masses
A_PLANET = 1.0              # Must match the original simulation
M_STAR = 1.0                 # Must match the original simulation
SAVE_INTERVAL = 2 * 60 * 60  # Seconds of wall time between restart checkpoints

# Original September 6 stopping parameters
slope_window = 5
minimum_comparison_steps = 20
consecutive_required = 3
minimum_fraction = 0.1
fraction_change_limit = 0.05
max_integration_time = 1e8
# -----------------------------------------------

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

clib = cdll.LoadLibrary(str(ROOT / 'Inputs/N_Body_Scripts/Heartbeat/heartbeat.so'))
clib.get_planet_cartesian.argtypes = [c_double] + [POINTER(c_double)] * 6
clib.get_planet_cartesian.restype = None
xp, yp, zp = c_double(), c_double(), c_double()
vxp, vyp, vzp = c_double(), c_double(), c_double()


def get_planet(t, m_planet):
    clib.get_planet_cartesian(float(t), byref(xp), byref(yp), byref(zp),
                              byref(vxp), byref(vyp), byref(vzp))
    return rebound.Particle(x=xp.value, y=yp.value, z=zp.value,
                            vx=vxp.value, vy=vyp.value, vz=vzp.value,
                            m=m_planet)


def load_sim(filename):
    sim = rebound.Simulation(str(filename))
    if sim.N != 1:
        raise ValueError(f'{filename}: expected one particle; found {sim.N}')
    # REBOUND does not restore Python/C callback pointers from an archive.
    sim.additional_forces = clib.planet_star_force
    return sim


def write_results(output_file, *values):
    with open(output_file, 'a') as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            f.write(','.join(map(str, values)) + '\n')
            f.flush()
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def read_old_results(filename, m_planet):
    """Return {pid: (status,time)}; NEVER infer identity from a 3-column log."""
    if not filename.exists():
        raise FileNotFoundError(filename)
    events = {}
    with open(filename) as f:
        for number, line in enumerate(f, 1):
            fields = line.strip().split(',')
            if not fields or len(fields) < 3:
                continue
            try:
                mass = float(fields[0]); event_t = float(fields[1])
            except ValueError:
                continue  # optional header
            if not np.isclose(mass, m_planet, rtol=1e-8):
                continue
            status = fields[2].strip()
            if status not in ('Ejected', 'Captured', 'Integrator Error'):
                continue
            if len(fields) < 4 or not fields[3].strip():
                raise ValueError(f'{filename}:{number}: missing particle_id. '
                                 'Original 3-column output cannot identify finished particles.')
            pid = int(fields[3]);
            if not 0 <= pid < N_total:
                raise ValueError(f'Invalid particle ID {pid} at {filename}:{number}')
            if pid in events:
                raise ValueError(f'Duplicate result for particle {pid} in {filename}')
            events[pid] = (status, event_t)
    return events


def calculate_a_e(t, particle, m_planet, G):
    planet = get_planet(t, m_planet)
    star_r = -(m_planet / M_STAR) * np.array([planet.x, planet.y, planet.z])
    star_v = -(m_planet / M_STAR) * np.array([planet.vx, planet.vy, planet.vz])
    r_vec = np.array([particle.x, particle.y, particle.z]) - star_r
    v_vec = np.array([particle.vx, particle.vy, particle.vz]) - star_v
    r = np.linalg.norm(r_vec)
    v2 = np.dot(v_vec, v_vec)
    mu = G * M_STAR
    energy = 0.5 * v2 - mu / r
    a = -mu / (2.0 * energy)
    h2 = np.dot(np.cross(r_vec, v_vec), np.cross(r_vec, v_vec))
    e = np.sqrt(max(0, 1 + 2 * energy * h2 / mu ** 2))
    return a, e


def save_remaining_particles(states, output_file, m_planet, G):
    local_rows = []
    for pid, state in states.items():
        if state['done']:
            continue
        sim = state['sim']
        a, e = calculate_a_e(sim.t, sim.particles[0], m_planet, G)
        local_rows.append(dict(particle_id=pid, m_planet=m_planet,
                               time=sim.t, a=a, e=e))
    gathered = comm.gather(local_rows, root=0)
    if rank == 0:
        rows = [row for group in gathered for row in group]
        pd.DataFrame(rows, columns=['particle_id', 'm_planet', 'time', 'a', 'e']).to_csv(output_file, index=False)


def update_stopping(fraction_record, t, fraction, cadence,
                    peak_slope, t_peak_slope, consecutive_below):
    """Exactly the slope/window comparison from the original script."""
    fraction_record.append(fraction)
    stop = False
    if len(fraction_record) > slope_window:
        slope = (fraction_record[-1] - fraction_record[-(slope_window + 1)]) / (slope_window * cadence)
        if slope > peak_slope:
            peak_slope, t_peak_slope = slope, t
    if fraction >= minimum_fraction and t_peak_slope is not None:
        comparison_cadence = max(2 * t_peak_slope, minimum_comparison_steps * cadence)
        n_steps_back = int(round(comparison_cadence / cadence))
        if len(fraction_record) > n_steps_back:
            difference = fraction_record[-1] - fraction_record[-(n_steps_back + 1)]
            if difference < fraction_change_limit:
                consecutive_below += 1
            else:
                consecutive_below = 0
            stop = consecutive_below >= consecutive_required
    return stop, peak_slope, t_peak_slope, consecutive_below

def find_checkpoint(pid):
    original = original_directory / f"particle_{pid}.bin"
    restarted = restart_directory / f"particle_{pid}.bin"

    if restarted.exists():
        return restarted

    if original.exists():
        return original

    raise FileNotFoundError(
        f"No checkpoint found for particle {pid}"
    )


if __name__ == '__main__':
    if len(sys.argv) != 3:
        raise SystemExit('Usage: python restart_mass_fraction_simple.py ORIGINAL_JOB_ID MASS_INDEX')
    job_id = sys.argv[1]
    mass_index = int(sys.argv[2])
    directory_name = MASS_DIRECTORIES[mass_index]
    if not directory_name.startswith(f'{job_id}_mass_'):
        raise ValueError(f'{directory_name} does not match original job ID {job_id}')
    m_planet = float(directory_name.split('_mass_', 1)[1])

    original_directory = ROOT / 'Outputs/Sim_Archives' / directory_name
    original_results = ROOT / 'Outputs/Ejection_Results' / f'{job_id}_ejection_results.txt'
    restart_directory = ROOT / 'Outputs/Sim_Archives_Restarted' / job_id / directory_name
    result_directory = ROOT / 'Outputs/Ejection_Results_Restarted' / job_id
    output_file = result_directory / f'{directory_name}_restarted_results.txt'
    final_position_file = result_directory / f'{directory_name}_remaining_particles.csv'

    if rank == 0:
        missing_directory = not original_directory.is_dir()
    else:
        missing_directory = None

    missing_directory = comm.bcast(missing_directory, root=0)

    if missing_directory:
        raise FileNotFoundError(original_directory)

    restart_directory.mkdir(parents=True, exist_ok=True)
    result_directory.mkdir(parents=True, exist_ok=True)

    old_events = read_old_results(original_results, m_planet)
    if output_file.exists():
        restarted_events = read_old_results(output_file, m_planet)

        duplicates = set(old_events) & set(restarted_events)

        if duplicates:
            raise ValueError(
                f"Duplicate terminal events: {sorted(duplicates)}"
            )

        old_events.update(restarted_events)

    planet_period = np.sqrt(A_PLANET ** 3 / (M_STAR + m_planet))
    cadence = 1e3 * planet_period
    temp_sim = rebound.Simulation()
    temp_sim.units = ('yr', 'AU', 'Msun')
    G = temp_sim.G
    clib.setup_force(c_double(A_PLANET), c_double(0), c_double(0),
                     c_double(0), c_double(0), c_double(0), c_double(0),
                     c_double(m_planet), c_double(M_STAR), c_double(G))

    # Assign particles exactly as in the original simulation.
    my_particle_ids = range(rank, N_total, size)
    states = {}
    local_greatest_t = -np.inf
    for pid in my_particle_ids:
        if pid in old_events:
            states[pid] = dict(sim=None, done=True, ejected=old_events[pid][0] == 'Ejected')
            continue

        filename = find_checkpoint(pid)
        sim = load_sim(filename)
        local_greatest_t = max(local_greatest_t, sim.t)
        states[pid] = dict(sim=sim, done=False, ejected=False,
                           save_file=restart_directory / f'particle_{pid}.bin',
                           last_save=time.monotonic())

    greatest_t = comm.allreduce(local_greatest_t, op=MPI.MAX)
    if not np.isfinite(greatest_t):
        if rank == 0:
            print('All particles already have terminal results; nothing to restart.', flush=True)
        sys.exit(0)
    if rank == 0:
        print(f'Restarting {directory_name} with {size} ranks; greatest_t={greatest_t:.8g}', flush=True)

    # Rebuild original stopping statistics on the original cadence grid,
    # counting ONLY ejections, using their original event timestamps.
    old_ejection_times = sorted(when for status, when in old_events.values() if status == 'Ejected')
    fraction_record = []
    peak_slope = -np.inf
    t_peak_slope = None
    consecutive_below = 0
    event_index = 0
    history_stop = False
    start_step = int(np.ceil((greatest_t - 1e-9 * cadence) / cadence))
    if rank == 0:
        for step in range(start_step):
            t = step * cadence
            while event_index < len(old_ejection_times) and old_ejection_times[event_index] <= t:
                event_index += 1
            history_stop, peak_slope, t_peak_slope, consecutive_below = update_stopping(
                fraction_record, t, event_index / N_total, cadence,
                peak_slope, t_peak_slope, consecutive_below)
        print(f'Historical ejections: {len(old_ejection_times)}, '
              f'fraction before catch-up: {fraction_record[-1] if fraction_record else 0:.5f}', flush=True)
        if history_stop:
            print('Note: reconstructed history already satisfies the original stopping condition.', flush=True)

    # Bring all surviving simulations to the latest checkpoint.
    # An ejection during catch-up is logged as an ejection (requested behavior).
    for pid, state in states.items():
        if state['done']:
            continue
        sim = state['sim']
        if sim.t < greatest_t:
            try:
                sim.integrate(greatest_t)
            except rebound.NoParticles:
                write_results(output_file, m_planet, greatest_t, 'Ejected', pid)
                state['ejected'] = True
                state['done'] = True

    # Checkpoint the synchronized states before subsequent integration.
    for state in states.values():
        if not state['done']:
            state['sim'].save_to_file(str(state['save_file']))
            state['last_save'] = time.monotonic()

    # Continue the original integration loop on the shared cadence grid.
    start_time = time.monotonic()
    for step in range(start_step, int(max_integration_time / cadence) + 1):
        t = step * cadence
        if t < greatest_t:
            continue
        for pid, state in states.items():
            if state['done']:
                continue
            sim = state['sim']
            try:
                if sim.t < t:
                    sim.integrate(t)
            except rebound.NoParticles:
                write_results(output_file, m_planet, t, 'Ejected', pid)
                state['ejected'] = True
                state['done'] = True
                continue

            planet_centric_e = sim.particles[0].orbit(primary=get_planet(sim.t, m_planet)).e
            if planet_centric_e < 1.0:
                try:
                    sim.integrate(t + planet_period)
                except rebound.NoParticles:
                    write_results(output_file, m_planet, t, 'Ejected', pid)
                    state['ejected'] = True
                    state['done'] = True
                    continue
                planet_centric_e_2 = sim.particles[0].orbit(primary=get_planet(sim.t, m_planet)).e
                if planet_centric_e_2 < 1.0:
                    write_results(output_file, m_planet, t, 'Captured', pid)
                    state['done'] = True

            if not state['done'] and time.monotonic() - state['last_save'] >= SAVE_INTERVAL:
                sim.save_to_file(str(state['save_file']))
                state['last_save'] = time.monotonic()

        # Cumulative counts: original historical ejections at or before t,
        # plus all newly ejected particles (all currently <= t).
        local_new = sum(state['ejected'] for state in states.values() if state['sim'] is not None)
        total_new = comm.reduce(local_new, op=MPI.SUM, root=0)
        active = comm.allreduce(sum(not state['done'] for state in states.values()), op=MPI.SUM)
        if rank == 0:
            old_count = sum(when <= t for when in old_ejection_times)
            fraction = (old_count + total_new) / N_total
            stop, peak_slope, t_peak_slope, consecutive_below = update_stopping(
                fraction_record, t, fraction, cadence,
                peak_slope, t_peak_slope, consecutive_below)
            stop = stop or active == 0 or t >= max_integration_time
            if step % 10 == 0 or stop:
                print(f't={t:.0f} fraction={fraction:.5f} peak={peak_slope:.3e} '
                      f't_peak={t_peak_slope} below={consecutive_below}/{consecutive_required} '
                      f'active={active} stop={stop}', flush=True)
        else:
            stop = None
        stop = comm.bcast(stop, root=0)
        if stop:
            break

    # Save survivors and a summary after the stopping criterion fires.
    for state in states.values():
        if not state['done']:
            state['sim'].save_to_file(str(state['save_file']))
    save_remaining_particles(states, final_position_file, m_planet, G)
    if rank == 0:
        print(f'Completed {directory_name}', flush=True)
