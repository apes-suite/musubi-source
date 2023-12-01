require "common"
require "musubi"

simulation_name = simulation_name..'Re'..Re

-- define the input
input = {
          read = tracking_fol..'FlowAroundCyl_2D_global_shape_FlowAroundCyl_lastHeader.lua',
        }

-- define the output
output = {  -- some general information
            folder = 'harvest/',     -- Output location 
--            folder = mesh,     -- Output location 

           { -- first output

            format = 'VTU',   -- Output format 

           },

         }


