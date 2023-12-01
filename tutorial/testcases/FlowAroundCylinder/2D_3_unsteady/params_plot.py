show_plot = True
font_family = 'serif'
font_type = 'Times New Roman'
font_size = 14
write_output = True
plot = []

## velocityX over time 

plot.append(dict(kind = 'xy',data=['tracking/Re100/FlowAroundCyl_2D_probe_velocity_p00000.res'], col=[1,2], 
                       startplot = True, endplot = True, figname = 'velX_over_time',  ls = 'b', format = 'png',
                       title = 'probe velocityX over time',
                       xlabel = 'time (s)', ylabel = 'velocityX (${m}/{s}$)', label = 'velocityX' ))

plot.append(dict(kind = 'xy',data=['tracking/Re100/FlowAroundCyl_2D_probe_velocity_in_p00000.res'], col=[1,2], 
                       startplot = True, endplot = True, figname = 'velX_inlet_over_time',  ls = 'b', format = 'png',
                       title = 'probe velocityX at the inlet over time',
                       xlabel = 'time (s)', ylabel = 'velocityX (${m}/{s}$)', label = 'velocityX inlet' ))
## velocityY over time

plot.append(dict(kind = 'xy',data=['tracking/Re100/FlowAroundCyl_2D_probe_velocity_p00000.res'], col=[1,3], 
                       startplot = True, endplot = True, figname = 'velY_over_time',  ls = 'b', format = 'png',
                       title = 'probe velocityY over time',
                       xlabel = 'time (s)', ylabel = 'velocityY (${m}/{s}$)', label = 'velocityY' ))

plot.append(dict(kind = 'xy',data=['tracking/Re100/FlowAroundCyl_2D_probe_velocity_in_p00000.res'], col=[1,3], 
                       startplot = True, endplot = True, figname = 'velY_inlet_over_time',  ls = 'b', format = 'png',
                       title = 'probe velocityY at the inlet over time',
                       xlabel = 'time (s)', ylabel = 'velocityY (${m}/{s}$)', label = 'velocityY inlet' ))
## drag coefficient over time

plot.append(dict(kind = 'xy',data=['tracking/Re100/FlowAroundCyl_2D_cyl_coeff_p00000.res'], col=[1,2], 
                       startplot = True, endplot = False, ls = 'b', label = 'cD' ))
plot.append(dict(kind='axhline', endplot = False, val = 2.97, label = 'upper\_bound', ls = 'b--'))
plot.append(dict(kind='axhline', endplot = True, val = 2.93, label = 'lower\_bound', ls = 'b--',
                       figname = 'cD_over_time',title = 'drag coefficient over time',format = 'png',
                       xlabel = 'time (s)', ylabel = 'cD ( )',)) 
## lift coefficient over time

plot.append(dict(kind = 'xy',data=['tracking/Re100/FlowAroundCyl_2D_cyl_coeff_p00000.res'], col=[1,3], 
                       startplot = True, endplot = False, ls = 'b', label = 'cL' ))
plot.append(dict(kind='axhline', endplot = False, val = 0.49, label = 'upper\_bound', ls = 'b--'))
plot.append(dict(kind='axhline', endplot = True, val = 0.47, label = 'lower\_bound', 
                        ls = 'b--', xlabel = 'time (s)', ylabel = 'cL ( )', 
                        title = 'lift coefficient over time',figname = 'cL_over_time',  format = 'png',))

 
## rescirculation length

# plot.append(dict(kind='rLength', data=['tracking/Re100/FlowAroundCyl_2D_line_p00000_t*.res'],
#              col=[1,5], xlabel='x', ylabel='Velocity X',
#              startplot = True, endplot = True, ls = 'r', dx=0.0082, re_start =0.25, 
#              label = 'velocity x'))

## pressure difference over time
  
plot.append(dict(kind='deltaP', data=['tracking/Re100/FlowAroundCyl_2D_probe_pressure_p00000.res',
                      'tracking/Re100/FlowAroundCyl_2D_probe_pressure_p00000.res'], col=[1,2,3], 
                      startplot = True, endplot = False, label = '$\Delta p$', ls = 'r'))
plot.append(dict(kind='axhline', endplot = False, val = -0.105, label = 'upper\_bound', ls = 'b--'))
plot.append(dict(kind='axhline', endplot = True, val = -0.115, label = 'lower\_bound', ls = 'b--',
                      xlabel='time (s)', ylabel='$\Delta p ({N}/{m^{2}})$', figname = 'deltaP_over_time',format='png', #xmin=29.5,xmax=30,
                      convert = True, facs=[1.0,1.e-3], title = 'pressure difference over time', titlesize = 16))

# ## velocityY over length
# 
# plot.append(dict(kind='xy', data=['tracking/Re100/FlowAroundCyl_2D_line_p00000_t10.000E+00.res'],
#             col=[1,6], xlabel='x', ylabel='Velocity Y [${m}/{s}$]', figname = 'velY_over_length', format = 'png',
#             startplot = True, endplot = True, ls = '-b', title = 'velocity Y over length',
#             label = 'v'))
# 
# ## velocityX over length
# 
# plot.append(dict(kind='xy', data=['tracking/Re100/FlowAroundCyl_2D_line_p00000_t10.000E+00.res'],
#             col=[1,6], xlabel='x', ylabel='Velocity X [${m}/{s}$]', figname = 'velX_over_length', format = 'png',
#             startplot = True, endplot = True, ls = '-b', title = 'velocity X over length',
#             label = 'u'))

## velocity over y-axis
#  in the middle of the channel
# 
# plot.append(dict(kind = 'xy',data=['tracking/Re100/FlowAroundCyl_2D_line_velocity_middle_p00000_t3.000E+00.res'], col=[2,4], 
#                       startplot = True, endplot = True, figname = 'line_velX_middle',  ls = 'b', format = 'png',
#                       convert = True, facs = [1.0, 1.0], title = 'velocityX over vertical line in the middle of the channel',
#                       xlabel = 'y-coordinate [m]', ylabel = 'velocityX [${m}/{s}$]', label = 'velocityX' ))
# 
# # in the beginning
# plot.append(dict(kind = 'xy',data=['tracking/Re100/FlowAroundCyl_2D_line_velocity_in_p00000_t3.000E+00.res'], col=[2,4], 
#                       startplot = True, endplot = True, figname = 'line_velX_start',  ls = 'b', format = 'png',
#                       convert = True, facs = [1.0, 1.0], title = 'velocityX over vertical line in the entry of the channel',
#                       xlabel = 'y-coordinate [m]', ylabel = 'velocityX [${m}/{s}$]', label = 'velocityX' ))
# 
# 
# 
## frequency 
# 
# plot.append(dict(kind = 'FFT', data=['tracking/Re100/FlowAroundCyl_2D_probe_velocity_p00000.res'],
#                     col=[1,3],startplot = True, endplot = True, xlabel = 'frequency [$s^{-1}$]', ylabel = 'velocityY', ls ='r', 
#                     convert = True, facs = [1.0, 1.0],figname = 'frequency',format ='png', label = 'velocityY', title = 'velocityY over frequency'))
