-- Example script showing the harvesting configuration for Seeder mesh files --
--
-- Use this as input to sdr_harvesting to obtain visualization files.
-- -----------------------------------------------------------------------------

-- Mesh to visualize, this is the common mesh definition, as found in
-- treelm/config.lua.
-- Either a predefined mesh, or the prefix for a mesh on disk.
mesh = 'mesh/'

-- If restart table and mesh are present then priority is given to restart
restart = {
--   read = 'debug/final100_restart_ProtoData_header_100.lua'
}

-- Verbosity of logging:
logging = {level=2}

-- Define tracking objects to further restrict the created output data.
-- If no tracking is defined, all the information will be written 'as is' to
-- the output files as given by the output table below.
-- For example in the following, we only write out the treeIDs.
tracking = {
    --label = 'walls_periodicX',      -- a label for the tracking object,
    label = 'channel',     -- a label for the tracking object,
                           -- default: 'unnamed'
    folder = 'output/',    -- output folder name
                           -- if it ends in a '/', this is a directory,
                           -- which needs to exist!
    output = {
      format = 'vtk',      -- Format of the output data (set to harvest!)
      write_pvd = false    -- write pvd file containing all timesteps
                           -- Default is true
    },
--    variable = {'treeid','solidsphere'}, -- variables to track
    --variable = {'treeid'}, -- variables to track
    shape = {kind='all'},  -- shape to track, 'all' or canoND are typical
                           -- other options see treelm/config.lua for examples.
}

-- To dump complete mesh or restart file with out any tracking table
-- just default output table with basename, format as shown below
-- Configuration for the output to produce
output_folder = 'debug/'    -- Prefix to use for the produced files,
                            -- if it ends in a '/', this is a directory,
                            -- which needs to exist!

output = {
           format   = 'vtk',     -- Visualization file format, defaults to 'vtk'
                                 -- Currently only vtk is supported here.
           dataform = 'binary',   -- Form to write data in:
                                 --   - 'binary',  or
                                 --   - 'ascii'
                                 -- Default is 'binary'
           write_pvd = false     -- write pvd file containing all timesteps
                                 -- Default is true
}


