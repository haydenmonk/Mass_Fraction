"""Continue the second-to-last REBOUND archive snapshot with heartbeat gravity.

Edit ARCHIVE_FILE and END_TIME, then run: python test_ias15_warning.py
The archive is read only. IAS15 RuntimeWarnings are logged, not raised.
"""
from ctypes import CDLL, c_double
from pathlib import Path
import re
import warnings

import rebound

# ------------------------ USER SETTINGS ------------------------
ARCHIVE_FILE = Path(
    "/scratch/group/p.phy260085.000/Mass_Fraction/Outputs/Sim_Archives/"
    "2199100_mass_1.00000000e-05/particle_410.bin"
)
END_TIME = 105000.0  # Absolute simulation time (years), NOT years to add

HEARTBEAT_LIB = Path(
    "scratch/group/p.phy260085.000/Mass_Fraction/"
    "Inputs/N_Body_Scripts/Heartbeat/heartbeat.so"
)
A_PLANET = 1.0
M_STAR = 1.0
# None means extract planet mass from ..._mass_1.00000000e-05/ directory.
M_PLANET = None
# ---------------------------------------------------------------


def main():
    if not ARCHIVE_FILE.is_file():
        raise FileNotFoundError(ARCHIVE_FILE)
    if not HEARTBEAT_LIB.is_file():
        raise FileNotFoundError(HEARTBEAT_LIB)

    mass = M_PLANET
    if mass is None:
        match = re.search(r"_mass_([0-9.eE+-]+)$", ARCHIVE_FILE.parent.name)
        if not match:
            raise ValueError("Cannot infer planetary mass from directory; set M_PLANET.")
        mass = float(match.group(1))

    archive = rebound.Simulationarchive(str(ARCHIVE_FILE))
    if len(archive) < 2:
        raise ValueError(f"Archive has only {len(archive)} snapshot(s); need at least two.")

    # Loading snapshots can issue REBOUND's expected function-pointer warning.
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r".*reset function pointers.*",
            category=RuntimeWarning,
        )
        sim = archive[-2]

    if sim.N != 1 or sim.particles[0].m != 0.0:
        raise ValueError("Expected a single massless particle in this archive.")
    if END_TIME <= sim.t:
        raise ValueError(f"END_TIME ({END_TIME}) must exceed snapshot time ({sim.t}).")

    clib = CDLL(str(HEARTBEAT_LIB))
    clib.setup_force.argtypes = [c_double] * 10
    clib.setup_force.restype = None
    clib.setup_force(
        A_PLANET, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0, mass, M_STAR, sim.G,
    )
    sim.additional_forces = clib.planet_star_force
    # Important: retain clib for the entire integration.

    print(f"Archive: {ARCHIVE_FILE}")
    print(f"Snapshots: {len(archive)} (loaded index -2)")
    print(f"Planet mass: {mass:.12g} Msun")
    print(f"Initial t={sim.t:.12g}, dt={sim.dt:.12g}, G={sim.G:.12g}")
    print(f"Target t={END_TIME:.12g}", flush=True)

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", RuntimeWarning)
        # Match production's default exact_finish_time=True.
        sim.integrate(END_TIME)

    print(f"Final t={sim.t:.12g}, dt={sim.dt:.12g}")
    p = sim.particles[0]
    print(f"Final xyz=({p.x:.12g}, {p.y:.12g}, {p.z:.12g})")
    print(f"Runtime warnings: {len(caught)}")
    for index, warning in enumerate(caught, 1):
        print(f"  {index}. {warning.category.__name__}: {warning.message}")
    print("Archive not modified.")


if __name__ == "__main__":
    main()
