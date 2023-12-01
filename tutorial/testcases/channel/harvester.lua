package.path = package.path .. ";../../seeder.lua"
require 'seeder'

name = label..'_l'..level
--filename = 'tracking/'..idLabel..'_'..label..'_hvsXY_l'..level..'_'..simName..'_lastHeader.lua'
filename = 'tracking/channel_lbm_incompbgkacoustic_hvsXY_channel_lastHeader.lua'
visualizeMesh = false

simulation_name = 'cylinder'

-- define the input
if visualizeMesh then
--! [input mesh visualisation]
input = {
  mesh = './mesh/'
  }
--! [input mesh visualisation]
else
--! [input mesh false]
input = {
  read = filename
  }
--! [input mesh false]  
end


-- define the output
--! [output]
output = {  -- some general information
    folder = 'harvest/',     -- Output location 
   { 
     format = 'VTU',   -- Output format 
   }    
}
--! [output]
