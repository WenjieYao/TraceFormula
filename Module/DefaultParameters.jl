"""This part is used to define all parameters used"""
# Geometry parameters of the mesh
# Rectangular with center circle design domain
λ = 1.0           # Wavelength
L = 2.5           # Width of the rectangular domain
H = 2.5           # Height of the rectangular domain
rd = 0.5          # Radius of the design domain circle
rt = rd + 0.2     # Radius of the target circle
dpml = 0.5        # Thickness of PML

# Periodic cell
h1 = 1.0          # Height of the air region
h2 = λ / 2        # Height of the design region
h3 = 1.0          # Height of the substrate region
ht = 0.5 * h1     # Target position
hs = λ / 10       # Source position

# Rectangular with rectangular design domain
Hd = 1.5          # Height of the design region
Ld = 0.5          # Length of the design region
xd = -0.2         # Center off-set of the design region
xt = 0.2          # Distance of the target line

# Characteristic length (controls the resolution, smaller the finer)
resol = 25.0      # Number of points per wavelength
l1 = λ / resol    # Normal region
l2 = l1 / 2.0     # Design region
l3 = 2 * l1       # PML

# Physical parameters 
k = 2 * π / λ       # Wave number 
# Bloch wavevector
kb = VectorValue(0.0, 0.0)     
ω = k               # c=1
ϵ1 = 1.0            # Relative electric permittivity for material 1 (y > 0)
ϵ2 = 1.0            # Relative electric permittivity for material 2 (y < 0)
ϵ3 = 0.0            # Relative electric permittivity for potential waveguide
ϵd = 12.0           # Relative electric permittivity for design material
μ = 1.0             # Relative magnetic permeability for all materials
wg_center = [0, 0]  # Waveguide center if exist
wg_size = [0, 0]    # Waveguide size if exist
#LHp = [Inf, h1 + h2]
#LHn = [Inf, h3]

# PML parameters
R = 1e-10           # Tolerence for PML reflection 
σ1 = -3 / 4 * log(R) / dpml / √ϵ1
σ2 = -3 / 4 * log(R) / dpml / √ϵ2
σs = [σ1, σ2]
############  Optimization parameters #############
flag_f = true       # Turn on filter
flag_t = true       # Turn on threshold

# Filter and threshold paramters
r = [0.02 * λ, 0.02 * λ]  # Filter radius
β = 40.0                  # β∈[1,∞], threshold sharpness
η = 0.5                   # η∈[0,1], threshold center

α = 1.0 / (2 * 1000.0)    # Equivalent loss α = 1/2Q

# Number of subspace
K = 20

# Amplify g for NLopt
Amp = 1

# Sum over kx
nkx = 30
nparts = nkx / 2

Bρ = true          # Matrix B depend on parameters?
ρv = 0.5

# Foundary constraint parameters
c = 0#resol^4
lw = r[1]
ls = r[1]
ηe = fηe(lw / r[1])
ηd = fηd(lw / r[1])