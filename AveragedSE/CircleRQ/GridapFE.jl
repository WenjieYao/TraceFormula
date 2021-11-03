"""Gridap finite element implementations"""

struct GridapParameters
    FE_V            # Finite element function space V for fields
    FE_U            # Finite element function space U for fields
    Ω               # Whole mesh domain
    dΩ              # Numerical integration for whole mesh domain
    dΓ_t            # Numerical integration for target boudary
end


function GridapFE(meshfile, order, degree, diritags, targettags)
    model = GmshDiscreteModel(meshfile)

    # Test and trial finite element function space
    # Scalar-valued shape functions,
    # but a complex vector of free DOFs to represent the solution.
    # (this automatically leads to complex sparse matrices when assembling)
    reffe = ReferenceFE(lagrangian, Float64, order)
    FE_V = TestFESpace(model, reffe, dirichlet_tags = diritags, vector_type = Vector{ComplexF64})
    FE_U = FE_V

    ############### Integration domain ################
    Ω = Triangulation(model)
    dΩ = Measure(Ω, degree)

    # Desgin/Target/Source line tags
    Γ_t = BoundaryTriangulation(model; tags = targettags)
    dΓ_t = Measure(Γ_t, degree)

    gridap = GridapParameters(FE_V, FE_U, Ω, dΩ, dΓ_t)
    return gridap
end