import subprocess
import numpy as np
import matplotlib.pyplot as plt



c_program_path = "./periodic_tst"

f = np.linspace(0.5,1.8,100)

e = 15

if e==4:

    A0 = 0.0177
    n = 1.4
    folder = "/data/topo/periodic_1024/larger_step/A_4/out_p"
    
    A = A0-(f-0.5)*0.007

elif e==15:

    A0 = 0.0105
    n = 1.4
    folder = "/data/topo/periodic_1024/larger_step/A_15/out_p"
    
    A = A0*(1/(f+0.5))**n


    
print("A = ",A0)
print("n = ",n)

for i in range(len(f)):

    arguments = [str(A[i]), str(f[i])]

    output_file = folder+str(i)

    try:

        with open(output_file, "w") as out_file:

            result = subprocess.run([c_program_path] + arguments, stdout=out_file, text=True, check=True)

        print("Program Output written to", output_file)

    except subprocess.CalledProcessError as e:

        print("Error:", e)
