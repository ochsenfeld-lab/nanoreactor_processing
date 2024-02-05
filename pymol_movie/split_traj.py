from read_write_utils import *

# Split trajectory files in N (# time steps) trajectory files

traj = read_traj_file(sys.argv[1])
lines_per_file = len(traj[0]) + 2
smallfile = None
with open(sys.argv[1]) as bigfile:
    for lineno, line in enumerate(bigfile):
        if lineno % lines_per_file == 0:
            if smallfile:
                smallfile.close()
            small_filename = 'traj_files/small_traj_{}.xyz'.format(int(lineno/lines_per_file))
            smallfile = open(small_filename, "w")
        smallfile.write(line)
    if smallfile:
        smallfile.close()
