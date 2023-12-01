require 'common'
simulation_name = 'channel_cylinder_Re'..Re

input = {
  read = 'tracking/Re100/FlowAroundCyl_2D_global_FlowAroundCyl_lastHeader.lua',
}

output = {
  folder = 'harvest/',     -- Output location 
  { -- first output
    format = 'VTU',   -- Output format 
  },
}
