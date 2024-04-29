require 'args'
require 'func'

variable = {
  {
    name = 'spcie_real',
    ncomponents = 1,
    vartype = 'st_fun',
    st_fun = gauss_pulse_real,
  }
}
restart   = {read = 'restart/simulation_lastHeader.lua'}
NOtracking  = {
  label     = 'main',
  folder    = 'harvest/',
  output    = { format = 'vtk' },
  variable  = {'spc1_density','spcie_real'},
  -- variable  = {'velocity_phy'},
  shape     = {kind='all'},
}
output = { 
  folder = 'harvest/', 
  {format = 'VTU'}
} 
