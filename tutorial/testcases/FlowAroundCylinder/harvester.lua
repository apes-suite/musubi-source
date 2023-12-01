--! [mesh]
require "common"
simulation_name = 'FlowAroundCylinder2D'

-- define the input
mesh = 'mesh/'

-- define the output
output_folder = 'harvest/'   -- Output location
output = {
    format = 'vtk'   -- Output format 
}
--! [mesh]
