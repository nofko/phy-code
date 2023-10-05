import subprocess
import numpy as np
import matplotlib.pyplot as plt



c_program_path = "./periodic"

output_file = "data/out_p"


f = np.arange(5,18)/10
A = 0.0032*(1/(f+0.5))**1.7

for i in range(len(f)):

    arguments = [str(A[i]), str(f[i])]

    output_file = "data/out_p"+str(i)

    try:

        with open(output_file, "w") as out_file:

            result = subprocess.run([c_program_path] + arguments, stdout=out_file, text=True, check=True)

        print("Program Output written to", output_file)

    except subprocess.CalledProcessError as e:

        print("Error:", e)
