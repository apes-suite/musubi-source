----- FUNCTIONS ------
-- Element-wise addition of two vectors of length len
function add(a,b,len)
  local c = {}
  for i = 1, len do
    c[i] = a[i] + b[i]
  end
  return c
end

----- PARAMETERS --------
simulation_name = 'barton_bfs_h64_Stk_1e-2'

-- NOTE:
-- Top -> z = H
-- Bottom -> z = 0
-- North -> y = W
-- South -> y = 0
-- East -> x = L
-- West -> x = 0

--Dimensions
h_lat       = 64
h           = 2.0
dx          = h/h_lat
x_inlet     = 0.0
x_outlet    = 38*h
dx_eps      = 0.01*dx
bnd_length  = (x_outlet - x_inlet) + 2*dx
nBnd_length = bnd_length/dx
level       = math.ceil(math.log(nBnd_length)/math.log(2))

--Flow parameters
Re          = 400          -- flow Reynolds = 2*u*h/nu
Stk         = 1.0e-2       -- particle Stokes number = rho_p*d_p^2*u/(18*mu*2h)
Rho         = 10           -- density ratio rho_particle / rho_fluid
Fr          = 0.67         -- Froude number = u / (2h*g)^0.5 

omega       = 1.85
nu_lat      = (1.0/omega-0.5)/3.0    
cs_lat      = math.sqrt(1.0/3.0)

um_lat      = Re*nu_lat/(2*h_lat)     -- mean inflow velocity [lat] 
g_lat       = (um_lat/Fr)^2/(2*h_lat) -- gravitational acceleration [lat]


-- Fluid properties (water at 20 degrees Celsius)
mu_phy      = 1.0e-3                 -- Pa*s
rho_phy     = 998                    -- kg/m^3
nu_phy      = mu_phy/rho_phy         -- m^2/s
dt          = nu_lat*dx^2/nu_phy     -- s

um          = um_lat*dx/dt          -- mean inflow velocity [m/s]
g           = g_lat*dx/dt^2         -- gravitational acceleration [m/s^2]
rho_p       = Rho*rho_phy
dp          = math.sqrt( (36*Stk*mu_phy*h)/(rho_p*um) ) -- particle diameter [m]
dp_lat      = dp/dx                        -- particle diameter [lat units]
Vp          = (math.pi*dp^3)/6       -- volume of particle [m^3]

mp          = rho_p*(math.pi*dp^3/6)
Fbuoyancy   = -g*(rho_p-rho_phy)*Vp

-- Initialization
t_ref_lat   = math.ceil( 38*h_lat / (um_lat) )
t_ramp_lat  = math.ceil(0.5*t_ref_lat) 
t_ramp_phy  = t_ramp_lat * dt

domain_origin = { 0.0, 0.0, -1.0 } 

-- Bounding cube
bnd_offset    = {-dx, -dx, -dx}
bnd_origin = add(domain_origin, bnd_offset, 3)
bnd_length  = (2^level)*dx
seed_orig = { 40.0, 2.0, 0.5*dx }

Re_phy = um*(2*h)/nu_phy
Re_lat = um_lat*(2*h_lat)/nu_lat

Stk_phy = (rho_p*dp^2*um)/(18*mu_phy*2*h)
Fr_phy  = um/math.sqrt((2*h*g))
