import os
import rebound
from pathlib import Path

heartbeat_file = "heartbeat.c"
reb_so_libpath = rebound.__libpath__
reb_h_libpath  = Path(rebound.__file__).parent / "rebound.h"

cp_so = f'cp {reb_so_libpath} librebound.so'
cp_h  = f'cp {reb_h_libpath} .'

compile_1 = f'gcc -c -O3 -fPIC {heartbeat_file} -o heartbeat.o'
# https://github.com/hannorein/rebound/issues/491
compile_2 = f"gcc -L. -shared heartbeat.o -o heartbeat.so -lrebound -Wl,-rpath='{os.path.dirname(os.path.abspath(__file__))}'"

all_commands = [cp_so, cp_h, compile_1, compile_2]

for command in all_commands: os.system(command)
