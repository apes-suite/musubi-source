-- Example script showing the harvesting configuration for Musubi restart files --
--
-- Use this as input to mus_harvesting to obtain visualization files.
-- -----------------------------------------------------------------------------

simulation_name = 'test'

-- Restart file to visualize
restart = {read = 'restart/test_header_0.000E+00.lua'}

-- Verbosity of logging:
logging = {level=1}

-- Define tracking objects to further restrict the created output data.
-- If no tracking is defined, all the information will be written 'as is' to
-- the output files as given by the output table below.
-- For example in the following, we only write out the treeIDs.
tracking = {
    label = 'allIDs',      -- a label for the tracking object,
                           -- default: 'unnamed'
    folder = 'output/',     -- output folder name
                           -- if it ends in a '/', this is a directory,
                           -- which needs to exist!
    output = { 
      format = 'vtk'       -- Format of the output data (set to harvest!)
    },
    variable = {'treeid','solidsphere','density','momentum'}, -- variables to track
    --variable = {'treeid'}, -- variables to track
    shape = {kind='all'},  -- shape to track, 'all' or canoND are typical
                           -- other options see treelm/config.lua for examples.
}

-- To dump complete mesh or restart file with out any tracking table
-- just default output table with basename, format as shown below
-- Configuration for the output to produce
NOoutput_folder = 'output'    -- Prefix to use for the produced files,
                            -- if it ends in a '/', this is a directory,
                            -- which needs to exist!

NOoutput = {   
           format   = 'vtk',     -- Visualization file format, defaults to 'vtk'
                                 -- Currently only vtk is supported here.
           dataform = 'binary'   -- Form to write data in:
                                 --   - 'binary',  or
                                 --   - 'ascii'
                                 -- Default is 'binary'
}


