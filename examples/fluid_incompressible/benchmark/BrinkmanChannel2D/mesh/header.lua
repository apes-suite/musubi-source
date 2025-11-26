 version = 1
 label = 'fluidTree'
 comment = 'BrinkmanChannel2D_seeder'
 boundingbox = {
    origin = {
         -31.250000000000000E-03,
         -31.250000000000000E-03,
         -78.125000000000000E-03 
    },
    length =   32.000000000000000E+00 
}
 nElems = 32768
 minLevel = 10
 maxLevel = 10
 nProperties = 2
 effBoundingbox = {
    origin = {
           0.000000000000000E+00,
           0.000000000000000E+00,
         -15.625000000000000E-03 
    },
    effLength = {
          16.000000000000000E+00,
           2.000000000000000E+00,
          31.250000000000000E-03 
    } 
}
 property = {
    {
        label = 'has boundaries',
        bitpos = 3,
        nElems = 32768 
    },
    {
        label = 'has qVal',
        bitpos = 8,
        nElems = 1148 
    } 
}
