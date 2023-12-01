import numpy as np

data=np.loadtxt('tracking/Re100/FlowAroundCyl_2D_cyl_coeff_p00000.res', skiprows = 2, usecols = [0,1])
cD_max = np.nanmax(data[:,1])

data2 = np.loadtxt('tracking/Re100/FlowAroundCyl_2D_cyl_coeff_p00000.res', skiprows = 1, usecols = [0,2])
cL_max = np.nanmax(data2[:,1])
time = data2[np.argmax(data2[:,1]),0]
file = open("coeff_max.txt", "w")
file.write("t_0 = "+str(time)+'\n'+"cD_max = "+str(cD_max)+'\n'+"cL_max = "+str(cL_max))
file.close()
