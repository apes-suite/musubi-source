show_plot = True
font_family = 'serif'
font_type = 'Times New Roman'
font_size = 14
write_output = True
plot = []

colormap = 'jet'
colormap = 'gist_rainbow'

## Grid convergence

# ## drag coefficient for different nHeight
# 
# plot.append(dict(kind='axhline', startplot = True, endplot = False, val = 3.24, label = 'upper bound = 3.24'))
# plot.append(dict(kind='axhline', endplot = False, val = 3.22, ls = 'b--',label = 'lower bound = 3.22'))
# plot.append(dict(kind = 'xy',data=['sdr_nHeight64/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,2], 
#                         endplot = False, label = 'nHeight = 64'))
# plot.append(dict(kind = 'xy',data=['sdr_nHeight101/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,2], 
#                         endplot = False, label = 'nHeight = 101')) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight128/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,2], 
#                         endplot = False, label = 'nHeight = 128')) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight200/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,2], 
#                         endplot = False, label = 'nHeight = 200')) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight256/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,2], 
#                         endplot = False, label = 'nHeight = 256'))
# plot.append(dict(kind = 'xy',data=['sdr_nHeight512/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,2], 
#                         endplot = True, label = 'nHeight = 512',
#                         title = 'drag coefficient over time', xmin = 9, xmax = 10, ymin = 3.2, ymax = 3.4, 
#                         legend = dict(loc=1,ncol=3),
#                         xlabel = 'time (s)', ylabel = 'cD ()', format = 'png', figname = 'cwith_diff_nHeight'))

# ## velocityX over line
# 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight64/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,5], 
#                         startplot = True, endplot = False, label = 'nHeight = 64'))
# plot.append(dict(kind = 'xy',data=['sdr_nHeight101/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,5], 
#                         endplot = False, label = 'nHeight = 101')) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight128/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,5], 
#                         endplot = False, label = 'nHeight = 128')) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight200/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,5], 
#                         endplot = False, label = 'nHeight = 200')) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight256/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,5], 
#                         endplot = False, label = 'nHeight = 256'))
# plot.append(dict(kind = 'xy',data=['sdr_nHeight512/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,5], 
#                         endplot = True, label = 'nHeight = 512',
#                         title = 'velocityX over line',  
#                         legend = dict(loc=1,ncol=3),
#                         xlabel = 'x (m)', ylabel = 'u (m/s)', format = 'png', figname = 'velX_with_diff_nHeight'))

# ## velocityY over line
# 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight64/mus_Re100/tracking/channel2D_2D_line_p*.res'], col=[1,6], 
#                         startplot = True, endplot = False, label = 'nHeight = 64'))
# plot.append(dict(kind = 'xy',data=['sdr_nHeight101/mus_Re100/tracking/channel2D_2D_line_p*.res'], col=[1,6], 
#                         endplot = False, label = 'nHeight = 101')) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight128/mus_Re100/tracking/channel2D_2D_line_p*.res'], col=[1,6], 
#                         endplot = False, label = 'nHeight = 128')) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight200/mus_Re100/tracking/channel2D_2D_line_p*.res'], col=[1,6], 
#                         endplot = False, label = 'nHeight = 200')) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight256/mus_Re100/tracking/channel2D_2D_line_p*.res'], col=[1,6], 
#                         endplot = False, label = 'nHeight = 256'))
# plot.append(dict(kind = 'xy',data=['sdr_nHeight512/mus_Re100/tracking/channel2D_2D_line_p*.res'], col=[1,6], 
#                         endplot = True, label = 'nHeight = 512',
#                         xmin = 0, xmax = 2.2, ymin = -1.3, ymax = 1.3,
#                         title = 'velocityY over line', 
#                         legend = dict(loc=1,ncol=3),
#                         xlabel = 'x (m)', ylabel = 'v (m/s)', format = 'png', figname = 'velY_with_diff_nHeight'))

# ## lift coefficient for different nHeight
# 
# plot.append(dict(kind='axhline', startplot = True, endplot = False, val = 1.01, label = 'upper bound = 1.01'))
# plot.append(dict(kind='axhline', endplot = False, val = 0.99, ls = 'b--', label = 'lower bound = 0.99')) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight64/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,3], 
#                         endplot = False, label = 'nHeight = 64',)) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight101/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,3], 
#                         endplot = False, label = 'nHeight = 101',)) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight128/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,3], 
#                         endplot = False, label = 'nHeight = 128',)) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight200/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,3], 
#                         endplot = False, label = 'nHeight = 200',)) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight256/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,3], 
#                         endplot = False, label = 'nHeight = 256',)) 
# plot.append(dict(kind = 'xy',data=['sdr_nHeight512/mus_Re100/tracking/channel2D_2D_cyl_coeff_p*.res'], col=[1,3], 
#                         endplot = True, label = 'nHeight = 512', legend = dict(loc=1,ncol=3),
#                         title = 'lift coefficient over time',xmin = 9, xmax = 10, ymin = 0.8, ymax = 1.1, 
#                         xlabel = 'time (s)', ylabel = 'cL ()', format = 'png', figname = 'cL_with_diff_nHeight'))


## pressure difference for Re 100 and different nHeight

## Re = 100

plot.append(dict(kind='axhline', startplot = True, endplot = False, val = 2.50, label = 'upper bound = 2.50'))
plot.append(dict(kind='axhline', endplot = False, val = 2.46, ls = 'b--', label = 'lower bound = 2.46'))
plot.append(dict(kind='deltaP', data=['../prod/sdr_nHeight64/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res',
                      '../prod/sdr_nHeight64/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res'], col=[1,2,3], 
                      endplot = False, label = 'nHeight = 64'))
plot.append(dict(kind='deltaP', data=['sdr_nHeight101/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res',
                      'sdr_nHeight101/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res'], col=[1,2,3], 
                      endplot = False, label = 'nHeight = 101'))
plot.append(dict(kind='deltaP', data=['sdr_nHeight128/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res',
                      'sdr_nHeight128/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res'], col=[1,2,3], 
                      endplot = False, label = 'nHeight = 128'))
plot.append(dict(kind='deltaP', data=['sdr_nHeight200/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res',
                      'sdr_nHeight200/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res'], col=[1,2,3], 
                      endplot = False, label = 'nHeight = 200'))
plot.append(dict(kind='deltaP', data=['sdr_nHeight256/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res',
                      'sdr_nHeight256/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res'], col=[1,2,3], 
                      endplot = False, label = 'nHeight = 256'))
plot.append(dict(kind='deltaP', data=['sdr_nHeight512/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res',
                      'sdr_nHeight512/mus_Re100/tracking/channel2D_2D_probe_pressure_p00001.res'], col=[1,2,3], 
                      endplot = True, label = 'nHeight = 512',legend = dict(loc=1,ncol=3),
                      xlabel='time ()', ylabel='$\Delta p (Pa)$', xmin = 9, xmax = 10, ymin = 2.3, ymax = 2.7, 
                      figname = 'deltaP_over_time_Re100',format='png', title = 'pressure difference', titlesize = 16))

## probe pressure over line
## Re = 20

# plot.append(dict(kind='xy', data=['../prod/sdr_nHeight64/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,4], 
#                       startplot = True, endplot = False, label = 'nHeight = 64'))
# plot.append(dict(kind='xy', data=['sdr_nHeight101/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,4], 
#                       endplot = False, label = 'nHeight = 101'))
# plot.append(dict(kind='xy', data=['sdr_nHeight128/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,4], 
#                       endplot = False, label = 'nHeight = 128'))
# plot.append(dict(kind='xy', data=['sdr_nHeight200/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,4], 
#                       endplot = False, label = 'nHeight = 200'))
# plot.append(dict(kind='xy', data=['sdr_nHeight256/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,4], 
#                       endplot = False, label = 'nHeight = 256'))
# plot.append(dict(kind='xy', data=['sdr_nHeight512/mus_Re20/tracking/channel2D_2D_line_p*.res'], col=[1,4], 
#                       endplot = True, label = 'nHeight = 512', xlabel='x (m)', ylabel='p (Pa)',  
#                       xmin = 0, xmax = 2.2,# ymin = 0.10, ymax = 0.13,
#                       legend = dict(loc=1,ncol=3),
#                       figname = 'deltaP_over_line', format='png', title = 'probe pressure over channel', titlesize = 16))
# 
