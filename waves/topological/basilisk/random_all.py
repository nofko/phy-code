import subprocess
import numpy as np
import matplotlib.pyplot as plt



c_program_path = "./random"

f = np.linspace(0.5,1.8,50)

A0 = 0.0177
n = 1.4

#A= A0*(1/(f+0.5))**n

A = A0-(f-0.5)*0.007

print("A = ",A0)
print("n = ",n)


for i in range(len(f)):

    arguments = [str(A[i]), str(f[i])]

    output_file = "/data/topo/random_1024/A_4/out_p"+str(i)

    try:

        with open(output_file, "w") as out_file:

            result = subprocess.run([c_program_path] + arguments, stdout=out_file, text=True, check=True)

        print("Program Output written to", output_file)

    except subprocess.CalledProcessError as e:

        print("Error:", e)
