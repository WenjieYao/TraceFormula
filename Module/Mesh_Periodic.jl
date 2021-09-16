"""
Gmsh function that creates a periodic rectangular domain

Example paramters
# Geometry parameters of the mesh
L = 1.0           # Length of the normal region
λ = 1.0           # Wavelength (aribitrary unit)
h1 = 1.0          # Height of the air region
h2 = λ / 2        # Height of the design region
h3 = 1.0          # Height of the substrate region
dpml = 0.5        # Thickness of the PML
ht = 0.5 * h1     # Target position
hs = λ / 10       # Source position

# Characteristic length (controls the resolution, smaller the finer)
resol = 25.0      # Number of points per wavelength
l1 = λ/resol      # Normal region
l2 = l1/2.0       # Design region
l3 = 2*l1         # PML

geo_param = PeriodicGeometry(L, h1, h2, h3, ht, hs, dpml, l1, l2, l3)

"""

struct PeriodicGeometry
    L::Float64           # Length of the normal region  
    h1::Float64          # Height of the air region
    h2::Float64          # Height of the design region
    h3::Float64          # Height of the substrate region
    ht::Float64          # Target position
    hs::Float64          # Source position
    dpml::Float64        # Thickness of the PML
    # Characteristic length (controls the resolution, smaller the finer)
    l1::Float64          # Normal region
    l2::Float64          # Design region
    l3::Float64          # PML 
end


function MeshGenerator(geo_param::PeriodicGeometry, meshfile_name::String)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.clear()
    gmsh.model.add("geometry")

    # Add points
    gmsh.model.geo.addPoint(-geo_param.L/2, -geo_param.h3, 0, geo_param.l1, 1)
    gmsh.model.geo.addPoint( geo_param.L/2, -geo_param.h3, 0, geo_param.l1, 2)
    gmsh.model.geo.addPoint( geo_param.L/2, -geo_param.hs, 0, geo_param.l1, 3)
    gmsh.model.geo.addPoint( geo_param.L/2, 0, 0, geo_param.l2, 4)
    gmsh.model.geo.addPoint( geo_param.L/2, geo_param.h2, 0, geo_param.l2, 5)
    gmsh.model.geo.addPoint( geo_param.L/2, geo_param.h2+geo_param.ht, 0, geo_param.l1, 6)
    gmsh.model.geo.addPoint( geo_param.L/2, geo_param.h2+geo_param.h1, 0, geo_param.l1, 7)
    gmsh.model.geo.addPoint(-geo_param.L/2, geo_param.h2+geo_param.h1, 0, geo_param.l1, 8)
    gmsh.model.geo.addPoint(-geo_param.L/2, geo_param.h2+geo_param.ht, 0, geo_param.l1, 9)
    gmsh.model.geo.addPoint(-geo_param.L/2, geo_param.h2, 0, geo_param.l2, 10)
    gmsh.model.geo.addPoint(-geo_param.L/2, 0, 0, geo_param.l2, 11)
    gmsh.model.geo.addPoint(-geo_param.L/2, -geo_param.hs, 0, geo_param.l1, 12)
    gmsh.model.geo.addPoint( geo_param.L/2, geo_param.h2+geo_param.h1+geo_param.dpml, 0, geo_param.l3, 13)
    gmsh.model.geo.addPoint(-geo_param.L/2, geo_param.h2+geo_param.h1+geo_param.dpml, 0, geo_param.l3, 14)
    gmsh.model.geo.addPoint(-geo_param.L/2, -geo_param.dpml-geo_param.h3, 0, geo_param.l3, 15)
    gmsh.model.geo.addPoint( geo_param.L/2, -geo_param.dpml-geo_param.h3, 0, geo_param.l3, 16)
    # Add lines
    gmsh.model.geo.addLine( 1,  2,  1)
    gmsh.model.geo.addLine( 2,  3,  2)
    gmsh.model.geo.addLine(12,  3,  3)
    gmsh.model.geo.addLine( 1, 12,  4)
    gmsh.model.geo.addLine( 3,  4,  5)
    gmsh.model.geo.addLine(11,  4,  6)
    gmsh.model.geo.addLine(12, 11,  7)
    gmsh.model.geo.addLine( 4,  5,  8)
    gmsh.model.geo.addLine(10,  5,  9)
    gmsh.model.geo.addLine(11, 10, 10)
    gmsh.model.geo.addLine( 5,  6, 11)
    gmsh.model.geo.addLine( 9,  6, 12)
    gmsh.model.geo.addLine(10,  9, 13)
    gmsh.model.geo.addLine( 6,  7, 14)
    gmsh.model.geo.addLine( 8,  7, 15)
    gmsh.model.geo.addLine( 9,  8, 16)
    gmsh.model.geo.addLine( 7, 13, 17)
    gmsh.model.geo.addLine(13, 14, 18)
    gmsh.model.geo.addLine( 8, 14, 19)
    gmsh.model.geo.addLine(15,  1, 20)
    gmsh.model.geo.addLine(15, 16, 21)
    gmsh.model.geo.addLine(16,  2, 22)
    # Construct curve loops and surfaces 
    gmsh.model.geo.addCurveLoop([1, 2, -3, -4], 1)
    gmsh.model.geo.addPlaneSurface([1], 1)
    gmsh.model.geo.addCurveLoop([3, 5, -6, -7], 2)
    gmsh.model.geo.addPlaneSurface([2], 2)
    gmsh.model.geo.addCurveLoop([6, 8, -9, -10], 3)
    gmsh.model.geo.addPlaneSurface([3], 3)
    gmsh.model.geo.addCurveLoop([9, 11, -12, -13], 4)
    gmsh.model.geo.addPlaneSurface([4], 4)
    gmsh.model.geo.addCurveLoop([12, 14, -15, -16], 5)
    gmsh.model.geo.addPlaneSurface([5], 5)
    gmsh.model.geo.addCurveLoop([15, 17, 18, -19], 6)
    gmsh.model.geo.addPlaneSurface([6], 6)
    gmsh.model.geo.addCurveLoop([-20, 21, 22, -1], 7)
    gmsh.model.geo.addPlaneSurface([7], 7)
    # Physical groups
    gmsh.model.addPhysicalGroup(0, [15,16,13,14], 1)
    gmsh.model.setPhysicalName(0, 1, "DirichletNodes")
    gmsh.model.addPhysicalGroup(1, [21,18], 2)
    gmsh.model.setPhysicalName(1, 2, "DirichletEdges")
    gmsh.model.addPhysicalGroup(0, [4,5,10,11], 3)
    gmsh.model.setPhysicalName(0, 3, "DesignNodes")
    gmsh.model.addPhysicalGroup(1, [6,9], 4)
    gmsh.model.setPhysicalName(1, 4, "DesignEdges")
    gmsh.model.addPhysicalGroup(2, [3], 5)
    gmsh.model.setPhysicalName(2, 5, "Design")
    gmsh.model.addPhysicalGroup(2, [1,2], 6)
    gmsh.model.setPhysicalName(2, 6, "Substrate")
    gmsh.model.addPhysicalGroup(2, [4,5], 7)
    gmsh.model.setPhysicalName(2, 7, "Air")
    gmsh.model.addPhysicalGroup(2, [7,6], 8)
    gmsh.model.setPhysicalName(2, 8, "PML")
    gmsh.model.addPhysicalGroup(1, [12], 9)
    gmsh.model.setPhysicalName(1, 9, "Target")
    gmsh.model.addPhysicalGroup(1, [3], 10)
    gmsh.model.setPhysicalName(1, 10, "Source")
    gmsh.model.geo.synchronize()
    
    # Set periodic mesh on the left and right side
    gmsh.model.mesh.setPeriodic(1, [2], [4],
                            [1, 0, 0, geo_param.L, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    gmsh.model.mesh.setPeriodic(1, [5], [7],
                            [1, 0, 0, geo_param.L, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    gmsh.model.mesh.setPeriodic(1, [8], [10],
                            [1, 0, 0, geo_param.L, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    gmsh.model.mesh.setPeriodic(1, [11], [13],
                            [1, 0, 0, geo_param.L, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    gmsh.model.mesh.setPeriodic(1, [14], [16],
                            [1, 0, 0, geo_param.L, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    gmsh.model.mesh.setPeriodic(1, [17], [19],
                            [1, 0, 0, geo_param.L, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    gmsh.model.mesh.setPeriodic(1, [22], [20],
                            [1, 0, 0, geo_param.L, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1])
    # We can then generate a 2D mesh...
    gmsh.model.mesh.generate(2)
    
    

    # ... and save it to disk
    gmsh.write(meshfile_name)
    gmsh.finalize()
end

