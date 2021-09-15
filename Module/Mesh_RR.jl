"""
Gmsh function that creates a rectangular domain with rectangular design domain

Example paramters
# Geometry parameters of the mesh
λ = 1.0           # Wavelength
L = 2.0           # Length of the normal region
H = 2.0           # Height of the normal region
Hd = 1.5          # Height of the design region
Ld = 0.5          # Length of the design region

xd = -0.2         # Center off-set of the design region
xt = 0.2          # Distance of the target line
dpml = 0.5        # Thickness of PML

# Characteristic length (controls the resolution, smaller the finer)
resol = 50.0      # Number of points per wavelength
l1 = λ/resol      # Normal region
l2 = l1/2.0       # Design region
l3 = 2*l1         # PML

geo_param = RecRecGeometry(L, H, Ld, Hd, xd, xt, dpml, l1, l2, l3)

"""

struct RecRecGeometry
    L::Float64           # Length of the normal region  
    H::Float64           # Height of the normal region
    Ld::Float64          # Length of the design region
    Hd::Float64          # Height of the design region
    xd::Float64          # Center (xd,0) of design region
    xt::Float64          # Distance of the target line
    dpml::Float64        # Thickness of the PML
    # Characteristic length (controls the resolution, smaller the finer)
    l1::Float64          # Normal region
    l2::Float64          # Design region
    l3::Float64          # PML 
end


function MeshGenerator(geo_param::RecRecGeometry, meshfile_name::String)
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.clear()
    gmsh.model.add("geometry") # name it whatever you want

    # Add points
    gmsh.model.geo.addPoint(-geo_param.L/2-geo_param.dpml, -geo_param.H/2-geo_param.dpml, 0, geo_param.l3,  1)
    gmsh.model.geo.addPoint(-geo_param.L/2     , -geo_param.H/2-geo_param.dpml, 0, geo_param.l3,  2)
    gmsh.model.geo.addPoint( geo_param.L/2     , -geo_param.H/2-geo_param.dpml, 0, geo_param.l3,  3)
    gmsh.model.geo.addPoint( geo_param.L/2+geo_param.dpml, -geo_param.H/2-geo_param.dpml, 0, geo_param.l3,  4)
    gmsh.model.geo.addPoint(-geo_param.L/2-geo_param.dpml, -geo_param.H/2     , 0, geo_param.l3,  5)
    gmsh.model.geo.addPoint(-geo_param.L/2     , -geo_param.H/2     , 0, geo_param.l1,  6)
    gmsh.model.geo.addPoint( geo_param.L/2     , -geo_param.H/2     , 0, geo_param.l1,  7)
    gmsh.model.geo.addPoint( geo_param.L/2+geo_param.dpml, -geo_param.H/2     , 0, geo_param.l3,  8)
    gmsh.model.geo.addPoint(-geo_param.L/2-geo_param.dpml, geo_param.H/2      , 0, geo_param.l3,  9)
    gmsh.model.geo.addPoint(-geo_param.L/2     , geo_param.H/2      , 0, geo_param.l1, 10)
    gmsh.model.geo.addPoint( geo_param.L/2     , geo_param.H/2      , 0, geo_param.l1, 11)
    gmsh.model.geo.addPoint( geo_param.L/2+geo_param.dpml, geo_param.H/2      , 0, geo_param.l3, 12)
    gmsh.model.geo.addPoint(-geo_param.L/2-geo_param.dpml, geo_param.H/2+geo_param.dpml , 0, geo_param.l3, 13)
    gmsh.model.geo.addPoint(-geo_param.L/2     , geo_param.H/2+geo_param.dpml , 0, geo_param.l3, 14)
    gmsh.model.geo.addPoint( geo_param.L/2     , geo_param.H/2+geo_param.dpml , 0, geo_param.l3, 15)
    gmsh.model.geo.addPoint( geo_param.L/2+geo_param.dpml, geo_param.H/2+geo_param.dpml , 0, geo_param.l3, 16)
    gmsh.model.geo.addPoint(-geo_param.Ld/2+geo_param.xd , -geo_param.Hd/2    , 0, geo_param.l2, 17)
    gmsh.model.geo.addPoint( geo_param.Ld/2+geo_param.xd , -geo_param.Hd/2    , 0, geo_param.l2, 18)
    gmsh.model.geo.addPoint( geo_param.Ld/2+geo_param.xd ,  geo_param.Hd/2    , 0, geo_param.l2, 19)
    gmsh.model.geo.addPoint(-geo_param.Ld/2+geo_param.xd ,  geo_param.Hd/2    , 0, geo_param.l2, 20)
    gmsh.model.geo.addPoint( geo_param.Ld/2+geo_param.xd+geo_param.xt , -geo_param.H/2.1   , 0, geo_param.l1, 21)
    gmsh.model.geo.addPoint( geo_param.Ld/2+geo_param.xd+geo_param.xt ,  geo_param.H/2.1   , 0, geo_param.l1, 22)
    gmsh.model.geo.addPoint(-geo_param.Ld/2+geo_param.xd-geo_param.xt, -geo_param.H/2.1   , 0, geo_param.l1, 23)
    gmsh.model.geo.addPoint(-geo_param.Ld/2+geo_param.xd-geo_param.xt,  geo_param.H/2.1   , 0, geo_param.l1, 24)
    # Add lines
    gmsh.model.geo.addLine(  1,  2,  1)
    gmsh.model.geo.addLine(  2,  6,  2)
    gmsh.model.geo.addLine(  6,  5,  3)
    gmsh.model.geo.addLine(  5,  1,  4)
    gmsh.model.geo.addLine(  2,  3,  5)
    gmsh.model.geo.addLine(  3,  7,  6)
    gmsh.model.geo.addLine(  7,  6,  7)
    gmsh.model.geo.addLine(  3,  4,  8)
    gmsh.model.geo.addLine(  4,  8,  9)
    gmsh.model.geo.addLine(  8,  7, 10)
    gmsh.model.geo.addLine(  6, 10, 11)
    gmsh.model.geo.addLine( 10,  9, 12)
    gmsh.model.geo.addLine(  9,  5, 13)
    gmsh.model.geo.addLine(  7, 11, 14)
    gmsh.model.geo.addLine( 11, 10, 15)
    gmsh.model.geo.addLine(  8, 12, 16)
    gmsh.model.geo.addLine( 12, 11, 17)
    gmsh.model.geo.addLine( 10, 14, 18)
    gmsh.model.geo.addLine( 14, 13, 19)
    gmsh.model.geo.addLine( 13,  9, 20)
    gmsh.model.geo.addLine( 11, 15, 21)
    gmsh.model.geo.addLine( 15, 14, 22)
    gmsh.model.geo.addLine( 12, 16, 23)
    gmsh.model.geo.addLine( 16, 15, 24)
    gmsh.model.geo.addLine( 17, 18, 25)
    gmsh.model.geo.addLine( 18, 19, 26)
    gmsh.model.geo.addLine( 19, 20, 27)
    gmsh.model.geo.addLine( 20, 17, 28)
    gmsh.model.geo.addLine( 21, 22, 29)
    gmsh.model.geo.addLine( 23, 24, 30)
    gmsh.model.geo.addLine( 21, 23, 31)
    gmsh.model.geo.addLine( 22, 24, 32)
    #gmsh.model.geo.addLine( 17, 18, 25)
    #gmsh.model.geo.addLine( 18, 19, 26)
    #gmsh.model.geo.addLine( 19, 17, 31)
    # Construct curve loops and surfaces 
    gmsh.model.geo.addCurveLoop([  1,  2,  3,  4], 1)
    gmsh.model.geo.addCurveLoop([  5,  6,  7, -2], 2)
    gmsh.model.geo.addCurveLoop([  8,  9, 10, -6], 3)
    gmsh.model.geo.addCurveLoop([ 11, 12, 13, -3], 4)
    gmsh.model.geo.addCurveLoop([ -7, 14, 15,-11], 5)
    gmsh.model.geo.addCurveLoop([-10, 16, 17,-14], 6)
    gmsh.model.geo.addCurveLoop([-12, 18, 19, 20], 7)
    gmsh.model.geo.addCurveLoop([-15, 21, 22,-18], 8)
    gmsh.model.geo.addCurveLoop([-17, 23, 24,-21], 9)
    gmsh.model.geo.addCurveLoop([ 25, 26, 27, 28], 10)
    #gmsh.model.geo.addCurveLoop([ 25, 26,31], 10)
    
    gmsh.model.geo.addPlaneSurface([ 1],  1)
    gmsh.model.geo.addPlaneSurface([ 2],  2)
    gmsh.model.geo.addPlaneSurface([ 3],  3)
    gmsh.model.geo.addPlaneSurface([ 4],  4)
    gmsh.model.geo.addPlaneSurface([ 5,10],  5)
    gmsh.model.geo.addPlaneSurface([ 6],  6)
    gmsh.model.geo.addPlaneSurface([ 7],  7)
    gmsh.model.geo.addPlaneSurface([ 8],  8)
    gmsh.model.geo.addPlaneSurface([ 9],  9)
    gmsh.model.geo.addPlaneSurface([10], 10)
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.embed(1, [29], 2, 5)
    gmsh.model.mesh.embed(1, [30], 2, 5)
    gmsh.model.mesh.embed(1, [31], 2, 5)
    gmsh.model.mesh.embed(1, [32], 2, 5)

    # Physical groups
    gmsh.model.addPhysicalGroup(0, [1,2,3,4,5,8,9,12,13,14,15,16], 1)
    gmsh.model.setPhysicalName(0, 1, "DirichletNodes")
    gmsh.model.addPhysicalGroup(0, [17,18,19,20], 2)
    gmsh.model.setPhysicalName(0, 2, "DesignNodes")
    gmsh.model.addPhysicalGroup(1, [1,4,5,8,9,13,16,20,19,22,24,23], 3)
    gmsh.model.setPhysicalName(1, 3, "DirichletEdges")
    gmsh.model.addPhysicalGroup(1, [25,26,27,28], 4)
    gmsh.model.setPhysicalName(1, 4, "DesignEdges")
    gmsh.model.addPhysicalGroup(2, [1,2,3,4,6,7,8,9], 5)
    gmsh.model.setPhysicalName(2, 5, "PML")
    gmsh.model.addPhysicalGroup(2, [10], 6)
    gmsh.model.setPhysicalName(2, 6, "Design")
    gmsh.model.addPhysicalGroup(2, [5], 7)
    gmsh.model.setPhysicalName(2, 7, "Air")
    gmsh.model.addPhysicalGroup(1, [29], 8)
    gmsh.model.setPhysicalName(1, 8, "Target")
    gmsh.model.addPhysicalGroup(1, [30], 9)
    gmsh.model.setPhysicalName(1, 9, "Source")
    gmsh.model.addPhysicalGroup(1, [31], 10)
    gmsh.model.setPhysicalName(1, 10, "Target1")
    gmsh.model.addPhysicalGroup(1, [32], 11)
    gmsh.model.setPhysicalName(1, 11, "Target2")
    gmsh.model.mesh.generate(2)
    
    # ... and save it to disk
    gmsh.write(meshfile_name)
    gmsh.finalize()
end