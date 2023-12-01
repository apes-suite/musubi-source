import numpy as np

Re = $!musubi_Re!$

data=np.loadtxt('tracking/channel2D_2D_cyl_coeff_p*.res', skiprows = 2, usecols = [0,1])
cD_max = np.nanmax(data[:,1])

data2 = np.loadtxt('tracking/channel2D_2D_cyl_coeff_p*.res', skiprows = 1, usecols = [0,2])
cL_max = np.nanmax(data2[:,1])
time = data2[np.argmax(data2[:,1]),0]
file = open("coeff_max.txt", "w")
if Re <= 20:
    file.write("t_0 = "+str(time)
        +'\n'+"cD_max = "+str(cD_max)
        +'\n'+"upper reference value = 2.97"
        +'\n'+"lower reference value = 2.93"
        +'\n'+"cL_max = "+str(cL_max)
        +'\n'+"upper reference value = 0.49"
        +'\n'+"lower reference value = 0.47"
    )

if Re >20 && Re <= 100    
    file.write("t_0 = "+str(time)
        +'\n'+"cD = "+str(cD_max)
        +'\n'+"upper reference value = 5.59"
        +'\n'+"lower reference value = 5.57"
        +'\n'+"cL = "+str(cL_max)
        +'\n'+"upper reference value = 0.0104"
        +'\n'+"lower reference value = 0.0110"
    )

file.close()
