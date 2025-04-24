-- SEDI_SPHERE --
-- Validation case for hydrodynamic interactions (including moving particles across
-- the grid). A spherical particle is dropped in a tank (walls on all sides) at
-- several particle Reynolds numbers. The results can be compared to the experimental
-- work and simulations of:
-- [1] Ten Cate (2002). Particle Imaging Velocimetry experiments
-- and lattice-Boltzmann simulations on a single sphere
-- settling under gravity. Physics of Fluids. 14-11.

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
simulation_name = 'tencate_Re4_1'

-- NOTE:
-- Top -> z = H
-- Bottom -> z = 0
-- North -> y = W
-- South -> y = 0
-- East -> x = L
-- West -> x = 0

--Dimensions
Dia_p_phy   = 15e-3    -- particle diameter [m]
L           = 100e-3   -- Length of domain (x-extent)
W           = 100e-3   -- Width (y-extent)
H           = 160e-3   -- Height (z-extent)
Dia_p_lat   = 12
dx          = Dia_p_phy/Dia_p_lat
dx_eps      = 0.01*dx
bnd_length  = math.max(L,W,H) + 2*dx
nBnd_length = bnd_length/dx
level       = math.ceil(math.log(nBnd_length)/math.log(2))

domain_origin = { 0.0, 0.0, 0.0 } 

-- Bounding cube
bnd_offset    = {-dx, -dx, -dx}
bnd_origin = add(domain_origin, bnd_offset, 3)
bnd_length  = (2^level)*dx
seed_orig = { 0.5*L, 0.5*W, 0.5*H }

-- Flow parameters
omega       = 1.6
nu_lat      = (1.0/omega-0.5)/3.0    
cs_lat      = math.sqrt(1.0/3.0)


-- Fluid properties
mu_phy      = 212e-3                   -- Pa*s
rho_phy     = 965                   -- kg/m^3
nu_phy      = mu_phy/rho_phy         -- m^2/s
dt          = nu_lat*dx^2/nu_phy     -- s

--Particle properties 
rho_p_phy   = 1120                                        -- [kg/m3]
m_p_phy     = rho_p_phy*(4/3)*math.pi*(0.5*Dia_p_phy)^3   -- particle mass [kg]
us_phy      = 0.060                                       -- [m/s]
Re          = rho_phy*us_phy*Dia_p_phy/mu_phy   
St          = rho_p_phy*us_phy*Dia_p_phy/(9*rho_phy*nu_phy)

us_lat      = us_phy*dt/dx
Ma_lat      = us_lat/cs_lat

-- terminal velocity according to Stokes law
g_phy       = 9.81       -- gravitational acceleration [m/s^2]
F_gravity   = m_p_phy*g_phy                                    -- [N]
F_buoyancy  = g_phy*rho_phy*(4/3)*math.pi*(0.5*Dia_p_phy)^3    -- [N]

cs_phy      = cs_lat*dx/dt
p0          = rho_phy*cs_phy^2         -- ambient pressure
-- Initial particle position
x_p_phy     = 0.5*L
y_p_phy     = 0.5*W
z_p_phy     = 120e-3 + 0.5*Dia_p_phy
