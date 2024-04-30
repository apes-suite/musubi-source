require 'musubi'
-- require 'params'

restart   = {read = 'restart/'..simulation_name..'_lastHeader.lua'}
variable = {
  {
     name = 'vel_x',
     ncomponents = 1,
     vartype = "operation",
     operation = {kind='extract',
                  input_varname={'velocity_phy'},
                  input_varindex = {1}
                 }
  },
  {
      name = 'uxval',
      ncomponents = 1,
      vartype = 'st_fun',
      st_fun = uxanaval,
  },
  {
      name = 'pval',
      ncomponents = 1,
      vartype = 'st_fun',
      st_fun = panaval,
  },
  {
    name = 'uxdiff',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind = 'difference',
      input_varname = {
        'uxval', 'vel_x'
      },
    },
  },
  {
    name = 'pdiff',
    ncomponents = 1,
    vartype = "operation",
    operation = {
      kind = 'difference',
      input_varname = {
        'pval', 'pressure_phy'
      },
    },
  },
}

tracking  = {
  label     = 'main',
  folder    = 'harvest/',
  output    = { format = 'ascii' },
  variable  = {'uxdiff', 'pdiff'},
  reduction = {'l2norm', 'l2norm'},
  -- variable  = {'velocity_phy'},
  shape     = {kind='all'},
}
