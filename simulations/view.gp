# usage: gnuplot -persist -c plot.gp snap_000100.bin
file = ARG1
nx = int(system("od -An -t d4 -j0 -N4 ".file))
ny = int(system("od -An -t d4 -j4 -N4 ".file))
t  = real(system("od -An -t f8 -j8 -N8 ".file))

set title sprintf("t = %.3f", t)
set size ratio -1
set palette defined (0 "blue", 0.5 "white", 1 "red")
set cbrange [-1:1]
plot file binary skip=16 array=(nx,ny) format='%float64' with image notitle

# usage: gnuplot -c live.gp
set size ratio -1
set palette defined (0 "blue", 0.5 "white", 1 "red")
set cbrange [-1:1]
last = ""
while (1) {
  file = system("ls -t snap_*.bin 2>/dev/null | head -1")
  if (strlen(file) > 0 && file ne last) {
    nx = int(system("od -An -t d4 -j0 -N4 ".file))
    ny = int(system("od -An -t d4 -j4 -N4 ".file))
    t  = real(system("od -An -t f8 -j8 -N8 ".file))
    set title sprintf("%s   t = %.3f", file, t)
    plot file binary skip=16 array=(nx,ny) format='%float64' with image notitle
    last = file
  }
  pause 0.3
}

