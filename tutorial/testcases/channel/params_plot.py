plot  = []

## Examples:
##http://matplotlib.github.com/api/pyplot_api.html#matplotlib.pyplot.plot
## HINT: startplot == False, can be used to generate multiple plots
show_plot = True
font_family= 'serif'
font_type= 'Times New Roman'
font_size = 14
### deltaP
plot.append(dict(kind = 'xy',data=['tracking/channel_lbm_incompbgkacoustic_probePressure_p00000.res'],
           col=[1,2], startplot = True, endplot = True, 
           label = 'pressure',figname = 'pressureOverTime', xlabel = 'time', ylabel = 'pressure',
           ls = 'k-'))

plot.append(dict(kind = 'xy',data=['tracking/channel_lbm_incompbgkacoustic_dpdx_p00000_t76.000E+00.res'],
           col=[1,4], startplot = True, endplot = True, 
           label = 'pressure',figname = 'pressureOverLength', xlabel = 'x', ylabel = 'pressure',
           ls = 'r-'))

plot.append(dict(kind = 'xy',data=['tracking/channel_lbm_incompbgkacoustic_velocity_p00000_t76.000E+00.res'],
           col=[2,4], startplot = True, endplot = True, 
           label = 'velocityX',figname = 'velocityOverHeight', xlabel = 'x', ylabel = 'velocity',
           ls = 'r-'))

