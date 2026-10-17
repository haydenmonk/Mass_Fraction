"""Resume a previous Mass_Fraction REBOUND run from its particle archives.

Launch with mpirun -np N python restart_mass_fraction.py JOB_ID MASS_INDEX.
Each MPI job processes exactly one sorted mass directory selected by MASS_INDEX.
Original archives/results are read-only. No convergence-based stopping is enabled.
"""

from pathlib import Path
from ctypes import cdll, c_double, POINTER, byref
import fcntl
import re
import sys
import time
import warnings

import numpy as np
import pandas as pd
import rebound
from mpi4py import MPI

# ===================== USER SETTINGS =====================
ROOT = Path("/scratch/group/p.phy260085.000/Mass_Fraction")
ORIGINAL_ARCHIVE_ROOT = ROOT / "Outputs/Sim_Archives"
HEARTBEAT_LIB = ROOT / "Inputs/N_Body_Scripts/Heartbeat/heartbeat.so"
N_TOTAL = 800
SAVE_INTERVAL = 2 * 60 * 60  # real elapsed seconds, per active particle
A_PLANET = 1.0
M_STAR = 1.0
SLOPE_WINDOW = 5
MINIMUM_COMPARISON_STEPS = 100
CONSECUTIVE_REQUIRED = 3
MINIMUM_FRACTION = 0.07
SLOPE_MINIMUM = 0.04
DIAGNOSTIC_EVERY_STEPS = 10
# =========================================================

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

clib = cdll.LoadLibrary(str(HEARTBEAT_LIB))
clib.get_planet_cartesian.argtypes = [c_double] + [POINTER(c_double)] * 6
clib.get_planet_cartesian.restype = None
xp, yp, zp = c_double(), c_double(), c_double()
vxp, vyp, vzp = c_double(), c_double(), c_double()


def get_planet(t, m_planet):
    clib.get_planet_cartesian(
        float(t), byref(xp), byref(yp), byref(zp),
        byref(vxp), byref(vyp), byref(vzp)
    )
    return rebound.Particle(
        x=xp.value, y=yp.value, z=zp.value,
        vx=vxp.value, vy=vyp.value, vz=vzp.value, m=m_planet
    )


def setup_force(m_planet):
    temp_sim = rebound.Simulation()
    temp_sim.units = ("yr", "AU", "Msun")
    G = temp_sim.G
    clib.setup_force(
        c_double(A_PLANET), c_double(0.0), c_double(0.0),
        c_double(0.0), c_double(0.0), c_double(0.0),
        c_double(0.0), c_double(m_planet), c_double(M_STAR),
        c_double(G),
    )
    return G


def discover_mass_directories(job_id):
    if not ORIGINAL_ARCHIVE_ROOT.is_dir():
        raise FileNotFoundError(ORIGINAL_ARCHIVE_ROOT)
    prefix = f"{job_id}_mass_"
    found = []
    for path in ORIGINAL_ARCHIVE_ROOT.iterdir():
        if path.is_dir() and path.name.startswith(prefix):
            suffix = path.name[len(prefix):]
            try:
                mass = float(suffix)
            except ValueError:
                continue
            found.append((mass, path))
    if not found:
        raise RuntimeError(f"No mass directories for JOB_ID={job_id!r}")
    return sorted(found)


def read_original_events(m_planet, original_results_file):
    """Return pid -> (status, event_time) for terminal events of this mass."""
    if not original_results_file.is_file():
        raise FileNotFoundError(original_results_file)
    df = pd.read_csv(
        original_results_file, header=None,
        names=["mass", "time", "status", "pid"],
        on_bad_lines="skip",
    )
    for field in ("mass", "time", "pid"):
        df[field] = pd.to_numeric(df[field], errors="coerce")
    df = df[
        np.isclose(df["mass"], m_planet, rtol=1e-8, atol=0)
        & df["status"].isin(["Ejected", "Captured", "Integrator Error"])
        & df["pid"].notna() & df["time"].notna()
    ].sort_values("time")
    events = {}
    for row in df.itertuples(index=False):
        pid = int(row.pid)
        if 0 <= pid < N_TOTAL:
            events.setdefault(pid, (row.status, float(row.time)))
    return events


def append_event(file, m_planet, t, status, pid):
    with open(file, "a") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        try:
            f.write(f"{m_planet},{t},{status},{pid}\n")
            f.flush()
        finally:
            fcntl.flock(f, fcntl.LOCK_UN)


def calculate_a_e(t, particle, m_planet, G):
    planet = get_planet(t, m_planet)
    star_r = -(m_planet / M_STAR) * np.array([planet.x, planet.y, planet.z])
    star_v = -(m_planet / M_STAR) * np.array([planet.vx, planet.vy, planet.vz])
    r_vec = np.array([particle.x, particle.y, particle.z]) - star_r
    v_vec = np.array([particle.vx, particle.vy, particle.vz]) - star_v
    r = np.linalg.norm(r_vec)
    v2 = np.dot(v_vec, v_vec)
    mu = G * M_STAR
    specific_energy = 0.5 * v2 - mu / r
    a = -mu / (2.0 * specific_energy)
    h2 = np.dot(np.cross(r_vec, v_vec), np.cross(r_vec, v_vec))
    e = np.sqrt(max(0.0, 1.0 + 2.0 * specific_energy * h2 / mu**2))
    return a, e


def new_diagnostics():
    return dict(fraction_record=[], peak_slope=-np.inf,
                t_peak_slope=None, consecutive_below=0)


def update_diagnostics(d, t, fraction, cadence):
    """Original stopping statistics; return original hypothetical stop decision."""
    record = d["fraction_record"]
    record.append(fraction)
    if len(record) > SLOPE_WINDOW:
        slope = (record[-1] - record[-(SLOPE_WINDOW + 1)]) / (SLOPE_WINDOW * cadence)
        if slope > d["peak_slope"]:
            d["peak_slope"] = slope
            d["t_peak_slope"] = t
    comparison_cadence, difference = None, None
    would_stop = False
    if fraction >= MINIMUM_FRACTION and d["t_peak_slope"] is not None:
        comparison_cadence = max(
            3 * d["t_peak_slope"], MINIMUM_COMPARISON_STEPS * cadence
        )
        n_steps_back = int(round(comparison_cadence / cadence))
        if len(record) > n_steps_back:
            difference = record[-1] - record[-(n_steps_back + 1)]
            if difference < SLOPE_MINIMUM:
                d["consecutive_below"] += 1
            else:
                d["consecutive_below"] = 0
            if d["consecutive_below"] >= CONSECUTIVE_REQUIRED:
                would_stop = True
    return would_stop, comparison_cadence, difference


def save_archive(state):
    state["sim"].save_to_file(str(state["archive_file"]))
    state["last_save"] = time.monotonic()


def run_mass(m_planet, original_dir, job_id):
    original_results_file = ROOT / "Outputs/Ejection_Results" / f"{job_id}_ejection_results.txt"
    restart_archive_root = ROOT / "Outputs/Sim_Archives_Restarted" / job_id
    restart_results_root = ROOT / "Outputs/Ejection_Results_Restarted" / job_id
    G = setup_force(m_planet)
    planet_period = np.sqrt(A_PLANET**3 / (M_STAR + m_planet))
    cadence = 1e3 * planet_period
    new_archive_dir = restart_archive_root / original_dir.name
    new_results_file = restart_results_root / f"{original_dir.name}_restarted_results.txt"
    new_remaining_file = restart_results_root / f"{original_dir.name}_remaining_particles.csv"
    new_archive_dir.mkdir(parents=True, exist_ok=True)
    restart_results_root.mkdir(parents=True, exist_ok=True)

    # Prevent accidentally appending a second run to the same restarted output.
    if rank == 0:
        any_outputs = new_results_file.exists() or new_remaining_file.exists() or any(new_archive_dir.glob("particle_*.bin"))
    else:
        any_outputs = None
    any_outputs = comm.bcast(any_outputs, root=0)
    if any_outputs:
        raise RuntimeError(f"Restart output already exists for {original_dir.name}; move it before rerunning.")

    completed = read_original_events(m_planet, original_results_file)
    states = {}
    local_start = np.inf
    for pid in range(rank, N_TOTAL, size):
        if pid in completed:
            status, when = completed[pid]
            states[pid] = dict(sim=None, done=True, ejected=(status == "Ejected"),
                               captured=(status == "Captured"), original_event=(status, when))
            continue
        source_file = original_dir / f"particle_{pid}.bin"
        if not source_file.exists():
            raise FileNotFoundError(f"No archive or terminal event for particle {pid}: {source_file}")
        sim = rebound.Simulation(str(source_file))
        sim.additional_forces = clib.planet_star_force
        if sim.N != 1:
            raise RuntimeError(f"Expected 1 particle, found {sim.N} in {source_file}")
        state = dict(sim=sim, done=False, ejected=False, captured=False,
                     original_event=None,
                     archive_file=new_archive_dir / f"particle_{pid}.bin",
                     last_save=time.monotonic())
        states[pid] = state
        local_start = min(local_start, sim.t)

    restart_time = comm.allreduce(local_start, op=MPI.MIN)
    n_active = comm.allreduce(sum(not s["done"] for s in states.values()), op=MPI.SUM)
    if rank == 0:
        print(f"Mass {m_planet:.8e}: {n_active} active; earliest checkpoint t={restart_time}", flush=True)
    if n_active == 0:
        return

    # Populate complete cumulative-ejection history from ORIGINAL event log,
    # so the slope diagnostics are not reset to zero on restart.
    all_original_ejection_times = sorted(t for status, t in completed.values() if status == "Ejected")
    d = new_diagnostics()
    start_index = max(0, int(np.floor(restart_time / cadence)))
    past_index = 0
    historical_count = 0
    for i in range(start_index):
        t_hist = i * cadence
        while (past_index < len(all_original_ejection_times)
               and all_original_ejection_times[past_index] <= t_hist):
            historical_count += 1
            past_index += 1
        update_diagnostics(d, t_hist, historical_count / N_TOTAL, cadence)

    # The first t may precede some particles' checkpoints. Never integrate backwards.
    i = start_index
    step = 0
    while True:  # no tmax or slope-triggered shutdown
        t = i * cadence
        for pid, state in states.items():
            if state["done"]:
                continue
            sim = state["sim"]
            if t <= sim.t:
                continue
            try:

                sim.integrate(t)
            except rebound.NoParticles:
                append_event(new_results_file, m_planet, t, "Ejected", pid)
                state["done"], state["ejected"] = True, True
                continue


            if time.monotonic() - state["last_save"] > SAVE_INTERVAL:
                save_archive(state)

            planet_centric_e = sim.particles[0].orbit(primary=get_planet(sim.t, m_planet)).e
            if planet_centric_e < 1.0:
                try:
                    sim.integrate(t + planet_period)
                except rebound.NoParticles:
                    append_event(new_results_file, m_planet, sim.t, "Ejected", pid)
                    state["done"], state["ejected"] = True, True
                    continue
                planet_centric_e_2 = sim.particles[0].orbit(primary=get_planet(sim.t, m_planet)).e
                if planet_centric_e_2 < 1.0:
                    append_event(new_results_file, m_planet, t, "Captured", pid)
                    state["done"], state["captured"] = True, True

        # Each rank reaches the MPI calls, even if it has no active particles.
        local_active = sum(not s["done"] for s in states.values())
        total_active = comm.allreduce(local_active, op=MPI.SUM)
        local_new_ejected = sum(s["ejected"] and s["original_event"] is None for s in states.values())
        new_ejected = comm.allreduce(local_new_ejected, op=MPI.SUM)
        if rank == 0:
            historical_count_now = sum(
                et <= t for et in all_original_ejection_times
            )
            fraction = (historical_count_now + new_ejected) / N_TOTAL

            would_stop, comparison_cadence, difference = update_diagnostics(
                d, t, fraction, cadence
            )

            if step % DIAGNOSTIC_EVERY_STEPS == 0:
                print(...)  # Keep your existing print statement

        else:
            would_stop = None

        would_stop = comm.bcast(would_stop, root=0)

        if would_stop or total_active == 0:
            for state in states.values():
                if not state["done"]:
                    save_archive(state)
            break

        i += 1
        step += 1

    # Save surviving particle positions (typically empty when all completed).
    rows = []
    for pid, state in states.items():
        if state["done"] or state["sim"] is None:
            continue
        sim = state["sim"]
        a, e = calculate_a_e(sim.t, sim.particles[0], m_planet, G)
        rows.append(dict(particle_id=pid, m_planet=m_planet, time=sim.t, a=a, e=e))
    gathered = comm.gather(rows, root=0)
    if rank == 0:
        pd.DataFrame([item for group in gathered for item in group],
                     columns=["particle_id", "m_planet", "time", "a", "e"]).to_csv(
                         new_remaining_file, index=False
                     )
        print(f"Completed mass={m_planet:.8e}", flush=True)
    comm.Barrier()


def main():
    if len(sys.argv) != 3:
        raise SystemExit("Usage: python restart_mass_fraction.py JOB_ID MASS_INDEX")
    job_id = sys.argv[1]
    if not job_id or "/" in job_id or job_id in (".", ".."):
        raise SystemExit("JOB_ID must be a simple directory-name component")
    try:
        mass_index = int(sys.argv[2])
    except ValueError as exc:
        raise SystemExit("MASS_INDEX must be an integer") from exc

    mass_dirs = discover_mass_directories(job_id)
    if not 0 <= mass_index < len(mass_dirs):
        raise SystemExit(
            f"MASS_INDEX={mass_index} outside 0..{len(mass_dirs) - 1}"
        )
    mass, directory = mass_dirs[mass_index]
    if rank == 0:
        print(f"Original JOB_ID={job_id}; array index={mass_index}; "
              f"selected mass={mass:.8e}; directory={directory}", flush=True)
    run_mass(mass, directory, job_id)


if __name__ == "__main__":
    main()
