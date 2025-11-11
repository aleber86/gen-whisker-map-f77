import numpy as np

with open("gwm_2.5_eta_fixed.dat", "r") as file:
    array = np.loadtxt(file)
    array = array.astype(np.float64)

partition = 512/16

for i in range(16):
    arr_to_file = array[int(i*partition):int((i+1)*partition), :]
    with open(f"aux_{i}.dat", "w") as file_out:
        np.savetxt(file_out,arr_to_file )
