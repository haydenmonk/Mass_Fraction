from mpi4py import MPI
import socket
import os

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print(
    f"hostname={socket.gethostname()} "
    f"SLURM_PROCID={os.environ.get('SLURM_PROCID')} "
    f"MPI rank={rank} of {size}",
    flush=True
)
